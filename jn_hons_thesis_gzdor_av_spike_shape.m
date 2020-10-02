
function [all_files_average_pulse_shape_normalized, mat_file_names] = jn_hons_thesis_gzdor_av_spike_shape ( pulse_threshold )

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % 
   % Name: jn_hons_thesis_gzdor_av_spike_shape.m 
   %
   % Usage: Has integer input pulse_threshold (use 1800) that is used to
   % detect all pulses for all the .mat files in the current directory. 
   %
   % INPUTS: pulse_threshold: (use 1800 as default) 
   %
   % OUTPUTS: 
   %         - plots of the average pulse shape for each input file overlaid,
   %         along with the average pulse shape for all input files
   %
   %         - all_files_average_pulse_shape_normalized: an array of the
   %         average pulse shape for input files 
   %
   % AUTHOR: Greg Zdor 
   %
   % DATE: September 2018 
   % 
   % Project: Assessing the Mean Neuronal Firing Rate Information Hypothesis via Mutual Information  
   %
   % Originally developed for Andrews University J.N. Honors Program Thesis Project 
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET STRUCTURE CONTAINING ALL .MAT DATA FILES OF INTEREST 
   
    current_folder = pwd; % get current folder 
    mat_files = dir(fullfile(current_folder,'*.mat')); % get all current .mat file info  
    mat_file_names = cell(length(mat_files), 1); % define cell array to hold names of mat files 

    for ii = 1 : length(mat_files) % loop through all mat files in currrent directory 
        mat_file_names{ii,1} = [mat_file_names{ii,1}; mat_files(ii,1).name]; % populate cell array with just file names 
    end 

    %% LOOP THROUGH ALL DATA FILES AND CALCULATE AND PLOT SPIKE SHAPE 

    average_pulse_shape_over_all_files_storage = []; 

    for ii = 1 : length(mat_file_names) % loop through all data files and calculate spike shape for each spike
       %{
        DEFINE CURRENT FILE, CONSTANTS, AND DETECTION THRESHOLDS
        %}

        current_file = mat_file_names{ii,1}; % get current file to load  
        load(current_file, 'ch1', 'ch2'); % load current file     
        clear pulse_ID_threshold; 
        clear new_pulse_ID_threshold;   
        pulse_ID_threshold = find (ch1 > pulse_threshold); % find ch1 values above spike detection threshold  
        delta = 35; 
        spike_peaks = 0;

        %{
        SECTION I: SLIDING DELTA WINDOW PERFORMS THRESHOLD PULSE DETECTION  
        %}

        for j=1 : length(pulse_ID_threshold)-1

            new_pulse_ID_threshold = pulse_ID_threshold(j)-delta+min(find(ch1(pulse_ID_threshold(j)-delta:pulse_ID_threshold(j)+delta)==max(ch1(pulse_ID_threshold(j)-delta:pulse_ID_threshold(j)+delta)))); %#ok<MXFND>

            if(sum(new_pulse_ID_threshold~=spike_peaks)==length(spike_peaks))

                spike_peaks = [spike_peaks new_pulse_ID_threshold]; % save new spike peaks data 
            end
        end

        spike_peaks = spike_peaks(2:length(spike_peaks));

        %{
        SECTION II: SELECT +/- delta to each of the above detected spike peaks 
        %} 

        for kk = 1 : length(spike_peaks)

            single_pulse_all_data_points(kk,:) = ch1(spike_peaks(kk)-delta : spike_peaks(kk) + delta);
        end
        %{
        SECTION III: PLOT ALL PULSES OVERLAID FOR EACH FILE 
        %}

        %{ 
        UNCOMMENT TO PLOT     

        figure(ii)

        for kk = 1 : length(spike_peaks)

        plot(single_pulse_all_data_points(kk,:))
        hold on
        end
        hold off
        title(['Pulses Overlaid for File: ',current_file],'Interpreter', 'none');
        xlabel('sample number'); 
        ylabel('signal amplitude'); 

        %} 

        %calculate average pulse shape for current data file 
        average_pulse_shape_current_file = mean(single_pulse_all_data_points);

        % save average_spike_shape to storage array 
        average_pulse_shape_over_all_files_storage = [average_pulse_shape_over_all_files_storage ; average_pulse_shape_current_file];  %#ok<AGROW>

    end 

    %% CALCULATE AVERAGE SPIKE SHAPE FOR ALL FILES 

    figure (ii + length(mat_file_names)*3 +1) 

    % plot average pulse shape for each data file, all on top of each other, hence hold on and hold off 


    for jj = 1 : size(average_pulse_shape_over_all_files_storage,1)
        hold on 
        subplot(2,1,1)
        plot(average_pulse_shape_over_all_files_storage(jj,:));    
    end 
    hold off 

    title(['Average Pulse Shapes for Each Of: ', num2str(ii), ' Files Considered ']);
    xlabel('sample number'); 
    ylabel('signal amplitude'); 

    % calculate the average pulse shape across all input files considered     
    all_files_average_pulse_shape = mean(average_pulse_shape_over_all_files_storage); 

    subplot(2,1,2)
    plot(all_files_average_pulse_shape);
    title(['Average Spike Shape for All: ', num2str(ii), ' Files Considered: Pulse Threshold: ', num2str(pulse_threshold)]);
    xlabel('sample number'); 
    ylabel('signal amplitude'); 
    legend('average pulse signal shape'); 

    % output to user 
    all_files_average_pulse_shape_normalized = all_files_average_pulse_shape / ... 
        max(all_files_average_pulse_shape); 

end % end function 
