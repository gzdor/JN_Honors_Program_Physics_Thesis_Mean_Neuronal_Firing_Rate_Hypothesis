
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % 
       % Name: jn_hons_thesis_gzdor_MI_pipeline_final_version.m 
       %
       % Usage: this script performs all the calculations for my honors
       % thesis
       %
       % INPUTS: none 
       %
       % OUTPUTS: plots  
       %   
       % AUTHOR: Greg Zdor 
       %
       % DATE: October 2018 
       %
       % LAST MODIFIED: 11.18.18 @ 8:33 PM  
       % 
       % Project: Assessing the Mean Neuronal Firing Rate Information Hypothesis via Mutual Information  
       %
       % Originally developed for Andrews University J.N. Honors Program Thesis Project 
       %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% LOAD DATA FILE 
        clear; % clear all existing variables
        tic; % start timer 
        %{
        - complete realizations are not including beginning or end points                                                                                   
        
        voc_m1ku299_t11_000850_9_realizations.mat           
        voc_m1ku299_t11_001452_2_realizations.mat           
        voc_m1ku299_t11_001554_2_realizations.mat           
        voc_m1ku299_t11_001708_1_realization.mat            
        voc_m1ku299_t11_003515_3_realizations.mat           
        voc_m1ku299_t11_013719_2_realizations.mat           
        voc_m1ku299_t11_013826_2_realizations.mat           
        voc_m1ku299_t11_013840_2_realizations.mat           
        voc_m1ku299_t11_013910_2_realizations.mat           
        voc_m1ku299_t11_014006_1_realization.mat            
        voc_m1ku299_t11_014023_3_realizations.mat           
        voc_m1ku299_t11_014042_3_realizations.mat           
        voc_m1ku299_t11_014126_2_realizations.mat           
        voc_m1ku299_t11_014155_3_realizations.mat           
        voc_m1ku299_t11_014243_3_realizations.mat           
        voc_m1ku299_t11_014704_3_realizations.mat           
        voc_m1ku299_t12_001003_3_realizations.mat           
        voc_m1ku342_t31_012622_2_realizations.mat           
        voc_m1ku342_t31_014126_3_realizations.mat           
        voc_m1ku342_t31_014759_1_realization.mat            
        voc_m1ku342_t31_015326_2_realizations.mat           
        voc_m1ku342_t31_020010_3_realizations.mat           
        voc_m1ku342_t32_003014_2_realizations.mat           
        voc_m1ku342_t32_003223_2_realizations.mat           
        voc_m1ku342_t32_003549_3_realizations.mat           
        voc_m1ku342_t32_010303_3_realizations.mat           
        voc_m1ku342_t32_010712_3_realizations.mat           
        voc_m1ku342_t32_015504_9_realizations.mat  
        %}

        %{ 
        9 realizations / file 
        voc_m1ku342_t32_015504_9_realizations.mat
        voc_m1ku299_t11_000850_9_realizations.mat
        %}
        
        
        
        %**************************************************************
        % Uncomment the following cell array and enter the names of the
        % files you wish to process; then uncomment the for loop in terms
        % of "uu" and run. 
        %**************************************************************
        
        %{
        files_loop_thru = {'voc_m1ku299_t11_000850_9_realizations.mat';...
            'voc_m1ku342_t32_015504_9_realizations.mat';...
            'voc_m1ku299_t11_003515_3_realizations.mat';...
            'voc_m1ku299_t11_014023_3_realizations.mat';...
            'voc_m1ku299_t11_014042_3_realizations.mat';...
            'voc_m1ku299_t11_014155_3_realizations.mat';...           
            'voc_m1ku299_t11_014243_3_realizations.mat';...
            'voc_m1ku299_t12_001003_3_realizations.mat';...
            'voc_m1ku342_t31_020010_3_realizations.mat';...
            'voc_m1ku342_t32_010712_3_realizations.mat';...
            'voc_m1ku342_t32_003549_3_realizations.mat';...
            'voc_m1ku342_t31_014126_3_realizations.mat'};
        
        store_MI_data_surrogate_over_tau_mult_data_files = cell(4,length(files_loop_thru));
        
        for uu = 1:length(files_loop_thru) 
        
        current_file = files_loop_thru{uu};
        
        %}
        % uncommment file to load 
        %current_file = 'voc_m1ku342_t32_015504_9_realizations.mat'; 
        current_file = 'voc_m1ku299_t11_000850_9_realizations.mat';
        load(current_file); 
       
        %for bb = 400:200:2400 % for incrementing through multiple bin sizes         
        %% DEFINE REALIZATION DETECTION PARAMETERS
        min_peak_height = 1500; 
        min_peak_distance = 80; 
        RS_thresh = -6000; % define realization_separation, hence RS  
        RS_min_separation = 25E3; %Fs = 50KHz, so set min separation width as 0.5 seconds 
        RS_thresh_data = find(ch1 < RS_thresh); % get ch1 data below threshold 
        realization_start_sample_num = []; 
        
        %% DEFINE BIN SIZE 
        bin_size = 1000; % change to bin_size = bb for incrementing through various bin size 
        close all; % close all open figures 
        % the number of surrogates to create 
        num_of_MI_values = 1E4;  % define how many MI outputs desired for surrogate data 
        %% CREATE NEW FOLDER TO SAVE OUTPUT OF PLOTS 
        new_folder_name = strcat('\plots_',current_file(1:end-4));
        currFolderName = pwd; % get current folder 
        new_folder_path = strcat(currFolderName,new_folder_name); 
        mkdir(new_folder_path);
        
        %% LOOP THROUGH ALL DATA AND FIND REALIZATION START TIMES 
        
        for jj = 1:length(RS_thresh_data) % find realization start points 
           if jj > 1 && jj <= length(RS_thresh_data) % non end points case 
               if RS_thresh_data(jj) > (RS_thresh_data(jj-1)+RS_min_separation)
                   realization_start_sample_num = [realization_start_sample_num RS_thresh_data(jj)]; 
               end 
           end % if non end points case 
           if jj == 1 
               realization_start_sample_num = [realization_start_sample_num RS_thresh_data(jj)];
           end 
        end % end for loop 
        
        %% PLOT REALIZATION START TIMES OVER TOP CH1, ALONG WITH SPECTROGRAM OF DATA 
        
        figure('units','normalized','outerposition',[0 0 1 1]) % plot realization points overlaid on ch1 
        subplot(2,1,1)
        hold on
        plot(ch1); 
        plot(realization_start_sample_num, ch1(realization_start_sample_num), 'o'); 
        xlabel('sample number');
        ylabel('ch1: signal amplitude'); 
        title(['Auditory Neuron Response for File: ', current_file, ' Realization Starts Circled',... 
           ], 'Interpreter','none'); 
        hold off 
        set(gca,'FontSize',20)
        set(gca, 'FontName', 'Times New Roman')

        subplot(2,1,2) % plot spectrogram of data 
        spectrogram(ch1,1024,5,1024*16,50000,'yaxis');
        xlabel('time (seconds)');
        ylabel('frequency (hertz)'); 
        title(['Spectrogram for file: ', current_file], 'Interpreter','none'); 
        set(gca,'FontSize',20)
        set(gca, 'FontName', 'Times New Roman')
        
        % save figure to data file specific folder 
        plot_name = 'whole_file_signal_mag_spect'; 
        figure_name = strcat(new_folder_path,'\','plot_',plot_name);
        print(figure_name,'-djpeg')
        
        
        %% PLOT DATA SAMPLE (FOR HONORS THESIS DOCUMENT) 
        figure('units','normalized','outerposition',[0 0 1 1]) 
        upper_lim = 3.75E5;
        lower_lim = 3.6E5; 
        sam_num = lower_lim:1:upper_lim;
        plot((sam_num./50000)-(lower_lim./50000),ch1(lower_lim:upper_lim)./max(ch1(lower_lim:upper_lim)),'b'); 
        xlabel('Time (seconds)'); 
        ylabel('Voltage (normalized)'); 
        title('Marmoset Auditory Neuron Voltage Response Over Time','fontsize',18); 
        legend('Voltage Response');  
        set(gca,'FontSize',20)
        set(gca, 'FontName', 'Times New Roman')
        
        % save figure to data file specific folder 
        plot_name = 'sample_signal'; 
        figure_name = strcat(new_folder_path,'\','plot_',plot_name);
        print(figure_name,'-djpeg')
        
        %% PLOT SPECTROGRAM AND SIGNAL MAGNITUDE DATA OF AUDITORY STIMULUS CH2 
        figure('units','normalized','outerposition',[0 0 1 1]) 
        Fs = 50E3; 
        subplot(2,1,1) 
        time = (1:length(ch2))/Fs;
        plot(time,ch2);
        xlabel('time (seconds)');
        ylabel('signal'); 
        title(['Auditory Neuron Stimulus Signal Over Time for File: ', current_file, ' Sampling Rate = ' num2str(Fs),' Hertz'],'Interpreter','none');
        set(gca,'FontSize',15)
        set(gca, 'FontName', 'Times New Roman')
        subplot(2,1,2) 
        spectrogram(ch2,1024,5,1024*16,Fs,'yaxis');
        colormap('jet');
        xlabel('time (seconds');
        ylabel('frequency (hertz)');
        title(['Auditory Stimulus Signal Spectrogram for File: ', current_file, ' Sampling Rate = ' num2str(Fs),' Hertz'],'Interpreter','none');
        set(gca,'FontSize',15)
        set(gca, 'FontName', 'Times New Roman')
        
        % save figure to data file specific folder 
        plot_name = 'auditory_stimulus_time_freq_signal'; 
        figure_name = strcat(new_folder_path,'\','plot_',plot_name);
        print(figure_name,'-djpeg')
        
        %% AUTO CORRELATION PLOT 
        % perform cross correlation 
        auto_corr = xcorr(ch1);
        % plot auto_corr
        figure('units','normalized','outerposition',[0 0 1 1]) 
        subplot(2,1,2)
        plot(auto_corr(1:round(length(auto_corr)/2)))
        xlabel('sample number'); 
        ylabel('auto correlation');
        title('Auto Correlation');
        set(gca,'FontSize',15)
        set(gca, 'FontName', 'Times New Roman')
        % axis([0 49999 0 0.1E10]);
        % plot actual data for comparison 
        subplot(2,1,1) 
        plot(ch1)
        xlabel('sample number'); 
        ylabel('signal');
        title('Signal Versus Sample Number');
        set(gca,'FontSize',15)
        set(gca, 'FontName', 'Times New Roman')
        
        % save figure to data file specific folder 
        plot_name = 'auto_correlation'; 
        figure_name = strcat(new_folder_path,'\','plot_',plot_name);
        print(figure_name,'-djpeg')
       
        figure('units','normalized','outerposition',[0 0 1 1]) 
        subplot(2,1,2)
        plot(auto_corr(1:round(length(auto_corr)/2)))
        xlabel('sample number'); 
        ylabel('auto correlation');
        title('Auto Correlation');
        axis([10E4 15E4 0 0.5E10]);
        set(gca,'FontSize',15)
        set(gca, 'FontName', 'Times New Roman')

        % plot actual data for comparison 
        subplot(2,1,1) 
        plot(ch1)
        xlabel('sample number'); 
        ylabel('signal');
        title('Signal Versus Sample Number');
        axis([10E4 15E4 min(ch1) max(ch1)]);
        set(gca,'FontSize',15)
        set(gca, 'FontName', 'Times New Roman')
        
        
        % save figure to data file specific folder 
        plot_name = 'auto_correlation_one_realization'; 
        figure_name = strcat(new_folder_path,'\','plot_',plot_name);
        print(figure_name,'-djpeg')
       

        %% CALCULATE REALIZATIONS WIDTHS 
         
        RS_widths = diff(realization_start_sample_num); % get all the realization widths  
        min_RS_width = min(RS_widths); % get the smallest realization width 
       %{ 
       NOTE: Partial realizations are ignored, that is, every last
       realization marker to the end of the file is cut off and
       ignored as it is incomplete 
       %} 
       all_realizations_per_file = []; 
       for jj = 1:length(RS_widths) % make all realizations the length of the min realization length of the realizations not including the last realization
            current_realization = ch1(realization_start_sample_num(jj):realization_start_sample_num(jj)+min_RS_width);  
            all_realizations_per_file = [all_realizations_per_file ;current_realization]; 
       end 
       
       %% FIND PULSE PEAKS VIA FINDPEAKS 
       spike_val_per_file = cell(1,length(RS_widths)); 
       spike_location_per_file = cell(1,length(RS_widths)); 
       
       for jj = 1:length(RS_widths) 
           current_realization = all_realizations_per_file(jj,:); 
           [spike_peak_val,spike_peak_loc] = findpeaks(current_realization,...
        'minPeakHeight',min_peak_height,'MinPeakDistance',min_peak_distance); % perform peak detection via peakfinder function 
           
            spike_val_per_file{1,jj} = spike_peak_val; % contains all spike peaks per file 
            spike_location_per_file{1,jj} = spike_peak_loc;
       end % end for loop through length(RS_widths)
       
       %% PLOT PULSE PEAKS, REALIZATION START MARKS OVERTOP CH1 
       for jj = 1:length(RS_widths) 
            % plot individual realization 
           if jj <= length(RS_widths) 
               figure('units','normalized','outerposition',[0 0 1 1])
               hold on 
               stem(spike_location_per_file{jj},...
                   spike_val_per_file{jj},'r'); 
               plot(ch1(realization_start_sample_num(jj):...
                   realization_start_sample_num(jj)+ min_RS_width),'b');
               num_sample_points = length(ch1(realization_start_sample_num(jj):...
                   realization_start_sample_num(jj)+min_RS_width)); 
               plot(min_peak_height*ones(num_sample_points,1),'--');
               xlabel('sample number'); 
               ylabel('signal'); 
               title(['Single Realization #: ', num2str(jj),' with Spike Peaks Circled (in red)'],'Interpreter','none');
               set(gca,'FontSize',15)
               set(gca, 'FontName', 'Times New Roman')
               hold off 
               
               % save figure to data file specific folder 
               plot_name = strcat('realization_',num2str(jj)); 
               figure_name = strcat(new_folder_path,'\','plot_',plot_name);
               print(figure_name,'-djpeg')
           end % end individual plots   
       end % end cycling through jj = ... max val of lengths(RS_widths) 
    
       %% CONVERT spike_peak_val cell array to binary matrix 
       %{ 
       1 = spike 
       0 = no spike 
       %} 
       truncate_factor = 1; % range of (0 1) decimal form of % of realization length to include 
       num_rows = size(all_realizations_per_file,1); 
       num_cols = size(all_realizations_per_file,2)*truncate_factor;
       binary_matrix = zeros(num_rows,num_cols);
       for jj = 1:num_rows 
           for kk = 1:num_cols
               current_realization_spike_loc = spike_location_per_file{1,jj};
               for mm = 1:length(current_realization_spike_loc) 
                   current_spike_loc = current_realization_spike_loc(mm);   
                   if current_spike_loc == kk 
                       binary_matrix(jj,kk) = 1; 
                   end 
               end 
           end 
       end
       
       
       %% FIND MEAN COLUMN WISE OF binary_matrix AND PLOT OVER TIME  
       prob_spike_over_time = mean(binary_matrix); 
       figure('units','normalized','outerposition',[0 0 1 1])
       plot(prob_spike_over_time)
       xlabel('sample number')
       ylabel('Probability of Pulse'); 
       title(['Probability of Pulse Over Time: ', current_file],'Interpreter','none');
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')
       % save figure to data file specific folder 
       plot_name = 'spike_spacing_all_realizations'; 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
       
       %% BIN DATA INTO LARGER TIME SEGMENTS TO SEE IF NON STATIONARITY OR NOT 
       length_realization = length(binary_matrix); 
       num_bins = ceil(length_realization/bin_size); 
       binned_binary_matrix = 0; 
       
       for jj = 1:num_bins-1 % bin data 
           if jj < num_bins && jj ~= 1 
               current_matrix_bin = binary_matrix(:,(bin_size*(jj-1))+1:bin_size*(jj)); 
               num_non_zero_elements = nnz(current_matrix_bin); 
               binned_binary_matrix = [binned_binary_matrix num_non_zero_elements]; 
           end 
           if jj == 1 
               current_matrix_bin = binary_matrix(:,1:bin_size*(jj)); 
               num_non_zero_elements = nnz(current_matrix_bin); 
               binned_binary_matrix = [binned_binary_matrix num_non_zero_elements]; 
           end 
           if jj >= length(binary_matrix) - bin_size && jj ~= 1
               current_matrix_bin = binary_matrix(:,length(binary_matrix)-bin_size:length(binary_matrix));
               num_non_zero_elements = nnz(current_matrix_bin); 
               binned_binary_matrix = [binned_binary_matrix num_non_zero_elements];
           end 
       end % end
       
       %% CALCULATE PROBABILITY OF SPIKE PER BIN AND PLOT 
       total_spikes = sum(binned_binary_matrix,2); % get total num of spikes across all realizations 
       probability_per_bin_matrix = binned_binary_matrix/total_spikes; 
       figure('units','normalized','outerposition',[0 0 1 1])
       subplot(2,1,1)
       plot(probability_per_bin_matrix)
       xlabel('Bin Number')
       ylabel('Probability of Pulse');
       title(['Probability of Pulse Versus Bin Number for bin size: ', num2str(bin_size), ' : ', current_file],'Interpreter','none');
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')
       subplot(2,1,2)
       plot(prob_spike_over_time)
       xlabel('sample number')
       ylabel('probability of Pulse'); 
       title(['Probability of Pulse Over Time: ', current_file],'Interpreter','none');
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')
       % save figure to data file specific folder 
       plot_name = strcat('probability_spike_per_binned',num2str(bin_size)); 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
       
       %% 
       %{
       ASSUMING STATIONARITY 
       = calculate MI over time 1st, then realization
       = (ind realization => sum over realizations => over tau)
       
       ASSUMING NON STATIONARITY 
       = calculate MI over realization, then time 
       = (bin all realizatons in time => calculate 1 MI value for given
       tau) 
       %} 
       
       % REASONING 
       % non stationarity approach is used, based off of figure 284
       % demonstrating firing rate is not stationary 
       
       %% CALCULATE MI OVER REALIZATION, THEN TIME 
       %{ 
       METHODOLOGY (steps): 
            1. first bin each realization individually 
            2. second, set a threshold that defines whether a bin contains
            a spike or not, and bin each realization indvidually 
            3. third, bin across all the rows for the length of the
            realizations
            4. calculate MI of this serieds as a function of tau 
            
            5. surrogate data generation 
                i. determine how many surrogate series you want 
                ii. using rand, generate 9 individual arrays each the
                length of the actual realizations 
                iii. third, determine the average firing rate (tot
                spikes)/(total time) for the actual data file and express
                this as a decimal in terms of (spikes / realization count) 
                iv. fourth, use this average firing rate as a threshold to
                bin each of the 9 realizations; 
                v. fifth, the binning is done on assigning a 0 to the bin
                if the average of the data in that bin is less than the
                average firing rate, and a 1 if it is more than the average
                firing rate 
                vi. for the binning, ensure to use the same binning size as
                you use on the actual data 
                vii. seventh, with the 9 binned surrogate arrays, bin these
                according to the threshold in 2. used for the actual data 
                viii. bin across the rows as in 3. for the actual data 
                ix. calculate MI on this series 
                x. do above steps i - ix for a large number of sets of 9
                surrogates so as to generate surrogates 
                xi. calculate the mean and standard deviation of the
                surrogates 
            6. plot MI(data) vs tau and MI(surrogate) vs tau 
            7. plot MI(data) vs tau and MI(surrogate) vs time (bin num) 
            8. 3D plot of MI(data), MI(surrogate) vs time and tau 
            9. calculate significance 
                S(significance) = (av_MI -
                av_MI(surrogate))/sigma(surrogate)  
            10. plot significance vs tau and time (bin num) 
       %} 
       
       %% DEFINE CONSTANTS
       % clearvars -except ch2 ch1 binary_matrix min_peak_distance min_peak_height...
       %     current_file min_peak_distance RS_thresh RS_min_separation RS_thresh_data...
       %    bin_size
       assumed_min_spike_spacing_in_data = 300; % units in sample number, or 6E-3 seconds
       
       %% CREATE SURROGATE DATA 
       % calculate the average number of spikes per sample number 
       percentage_realization_to_include = 1; % range of (0 1)
       realization_counts_to_include = percentage_realization_to_include*length(binary_matrix);
       tot_num_spikes = nnz(binary_matrix(:,1:realization_counts_to_include));
       average_firing_rate = tot_num_spikes/(percentage_realization_to_include.*length(ch1)); % units = spikes / sample number
       % define constants for surrogate data generation 
       num_realizations = size(binary_matrix,1);
       length_realization = length(binary_matrix);
       % define storage cell 
       realizations_over_MI_cell = cell(1,num_of_MI_values);
       % create num_of_MI_values number of random series for MI calculation 
       for jj = 1:num_of_MI_values 
           % create rand matrix 
           rand_matrix = rand([num_realizations,length_realization]);
           % treat each realization separately 
           all_realizations_one_MI = cell(1,num_realizations);
           for kk = 1:num_realizations 
               % bin according to average firing rate 
               %{
               below average firing rate = 1 
               above average firing rate = 0 
               %} 
               current_realization_surrogate = rand_matrix(kk,:);
               num_bins_exact = length(rand_matrix)/bin_size; 
               num_bins_floor = floor(num_bins_exact);
               realization_binary_matrix = 0;
               for mm = 1:num_bins_floor+1
                   % all bins between beginning and end 
                   if mm < num_bins_floor && mm ~= 1 
                       current_rand_matrix_bin = current_realization_surrogate(bin_size*(mm-1)+1:...
                           bin_size*(mm));
                       val = find(current_rand_matrix_bin < average_firing_rate);
                       num_val = numel(val);
                       if num_val > 1 
                           realization_binary_matrix = [realization_binary_matrix; 1];
                       else 
                           realization_binary_matrix = [realization_binary_matrix; 0];
                       end   
                   end % end case all bins between beginning and end   
                   % beginning bin 
                   if mm == 1 
                      current_rand_matrix_bin = current_realization_surrogate(1:bin_size);
                       val = find(current_rand_matrix_bin < average_firing_rate);
                       num_val = numel(val);
                       if num_val > 1
                           realization_binary_matrix = [realization_binary_matrix; 1];
                       else 
                           realization_binary_matrix = [realization_binary_matrix; 0];
                       end  
                   end
                   % ending bin
                   if mm >= num_bins_floor
                       current_rand_matrix_bin = current_realization_surrogate((length(current_realization_surrogate)-bin_size):...
                           (length(current_realization_surrogate)));
                       val = find(current_rand_matrix_bin < average_firing_rate);
                       num_val = numel(val);
                       if num_val > 1
                           realization_binary_matrix = [realization_binary_matrix; 1];
                       else 
                           realization_binary_matrix = [realization_binary_matrix; 0];
                       end
                   end 
               end % end binning current realization 
               realization_binary_matrix = realization_binary_matrix(3:end); % delete 0 initialization 
               all_realizations_one_MI{kk} = [all_realizations_one_MI{kk}, realization_binary_matrix];
           end % end looping through num_realizations
           realizations_over_MI_cell{jj} = [realizations_over_MI_cell{jj}; all_realizations_one_MI];
       end % end looping through num_MI_values

       % surrogate data contained in "realizations_over_MI_cell" 

       %% FROM SURROGATE DATA, GENERATE SERIES FOR MI CALCULATION 
       all_MI_series = zeros(length(realizations_over_MI_cell),length(realization_binary_matrix)).';
       for ii = 1:length(realizations_over_MI_cell)  
           most_outer_cell = realizations_over_MI_cell(1,ii);
           one_MI_matrix = zeros(num_realizations,length(realization_binary_matrix));
           for jj = 1:num_realizations 
               middle_cell = most_outer_cell{1,1}(1,jj);
               inner_cell = middle_cell{1,1};
               one_MI_matrix(jj,:) = inner_cell;
           end 
           sum_over_realization = sum(one_MI_matrix);
           all_MI_series(:,ii) = sum_over_realization;
       end 

       % surrogate data in all_MI_series variable 

       %% BIN ACTUAL DATA 
       % bin each realization individually using the same bin_size used
       % for the surrogates 
       num_bins_for_realization = ceil(length(binary_matrix)/bin_size);
       each_realiz_binned = cell(size(binary_matrix,1),num_bins_for_realization); 
       contains_spike_threshold = 1; % ceil(bin_size/assumed_min_spike_spacing_in_data);
       for ii = 1:size(binary_matrix,1) 
           current_data_realization = binary_matrix(ii,:);
           for jj = 1:num_bins_for_realization
               % case 1: bins between beginning and end 
               if jj < num_bins_for_realization && jj ~= 1
                   current_bin_segment = current_data_realization(bin_size*(jj-1)+1:...
                       bin_size*(jj)); 
                   num_non_zero_elements_in_bin = nnz(current_bin_segment);
                   % perform threshold filtering on nnz results
                   if num_non_zero_elements_in_bin > contains_spike_threshold
                       % then contains a spike, so assign a 1 
                       each_realiz_binned{ii,jj} = 1; 
                   end 
                   if num_non_zero_elements_in_bin <= contains_spike_threshold
                       % then contains a spike, so assign a 1 
                       each_realiz_binned{ii,jj} = 0; 
                   end 
               end % end case 1    
               % case 2: beginning bins 
               if jj == 1 
                   current_bin_segment = current_data_realization(1:bin_size);
                   num_non_zero_elements_in_bin = nnz(current_bin_segment);
                   % perform threshold filtering on nnz results
                   if num_non_zero_elements_in_bin > contains_spike_threshold
                       % then contains a spike, so assign a 1 
                       each_realiz_binned{ii,jj} = 1; 
                   end 
                   if num_non_zero_elements_in_bin <= contains_spike_threshold
                       % then contains a spike, so assign a 1 
                       each_realiz_binned{ii,jj} = 0; 
                   end 
               end % end case 2
               % case 3: ending bins 
               if jj == num_bins_for_realization 
                   current_bin_segment = current_data_realization(length(current_data_realization)-bin_size:length(current_data_realization));
                   num_non_zero_elements_in_bin = nnz(current_bin_segment);
                   % perform threshold filtering on nnz results
                   if num_non_zero_elements_in_bin > contains_spike_threshold
                       % then contains a spike, so assign a 1 
                       each_realiz_binned{ii,jj} = 1; 
                   end 
                   if num_non_zero_elements_in_bin <= contains_spike_threshold
                       % then contains a spike, so assign a 1 
                       each_realiz_binned{ii,jj} = 0; 
                   end 
               end % end case 3
           end % end cycling through current realization 
       end % end cycling through all realizations 

       % now sum each colummn along the rows of the each_realiz_binned
       actual_data = sum(cell2mat(each_realiz_binned)).'; 
       actual_data = actual_data(1:size(all_MI_series,1));
       % plot actual_data and average of surrogate_data 
       

       
       figure('units','normalized','outerposition',[0 0 1 1])
       hold on 
       actual_data_integral = trapz(actual_data); 
       mean_all_MI_series_integral = trapz(mean(all_MI_series,2));
       plot(actual_data,'b'); 
       plot(mean(all_MI_series,2),'r');
       xlabel('bins');
       ylabel('number of spikes / bin');
       title(['SPB Versus Bins: Data Integrated = ', num2str(actual_data_integral),...
       '; Surrogate Integrated = ' num2str(mean_all_MI_series_integral),' ',current_file,...
           ' Bin Size = ' num2str(bin_size)], 'Interpreter','none');
       legend('data','mean of surrogate');
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')
       hold off 
       
       % save figure 
       plot_name = strcat('spikes_per_bin_surrogate_data_over_bin_num',num2str(bin_size)); 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
       
       % create figure name for saving 
       descrip = 'binned_actual_average_surrogate_data_versus_bin_number';
       file_name = strcat(new_folder_path,'\','PLOT_',descrip,'_',num2str(bin_size),'.png');
       % create figure object and save figure to "file_name" 

       

       %% FUNCTION CALL TO CALCULATE MI OVER TAU FOR DATA, SURROGATES 
       % define tau bounds for calculating MI for data, surrogate 
       tau_lower_lim = 1; 
       tau_upper_lim = ceil(floor(length(binary_matrix)/bin_size) - 0.1.* floor(length(binary_matrix)/bin_size));
       surrogate_data = all_MI_series; 
       % check series for surrogate and actual data have same length 
       if size(surrogate_data,1) ~= size(actual_data,1) 
           error('surrogate_data and actual_data have different lengths: to fix, truncate or lengthen surrogate data series'); 
       end
       % function call to sub function to calculate MI over tau for
       % surrogate 
       MI_surrogate = cell2mat(MI_over_Tau(surrogate_data,tau_lower_lim,tau_upper_lim)); 
       % function call to sub function to calculate MI over tau for actual
       % data
       MI_data = cell2mat(MI_over_Tau(actual_data,tau_lower_lim,tau_upper_lim));

       %% PLOTS OF MI VERSUS TAU, MI VERSUS TIME, MI VERSUS TAU AND TIME
       % plot MI for data, surrogate vs tau 
       stan_dev = std(MI_surrogate);
       MS_mean = mean(MI_surrogate,1);
       figure('units','normalized','outerposition',[0 0 1 1]) 
       hold on 
       tau = tau_lower_lim:tau_upper_lim; 
       plot(tau,MI_data,'b'); 
       errorbar(tau,MS_mean,stan_dev,'r');
       plot(tau,MS_mean,'g');
       hold off 
       xlabel('tau (look ahead)'); 
       ylabel('Mutual Information'); 
       legend ('data','standard deviation of surrogate','mean of surrogate','Location','best'); 
       title(['Mutual Information Versus Tau for Surrogate and Data for file: ', current_file....
           ,' Bin Size = ', num2str(bin_size)],'Interpreter','none');
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')

       % store values for plotting over multiple data files, constant bin
       % size
       % store_MI_data_surrogate_over_tau_mult_data_files{1,uu} = MI_data;
       % store_MI_data_surrogate_over_tau_mult_data_files{2,uu} = MS_mean;
       % store_MI_data_surrogate_over_tau_mult_data_files{3,uu} = stan_dev;
       % store_MI_data_surrogate_over_tau_mult_data_files{4,uu} = tau;
       
       % save figure 
       plot_name = strcat('MI_data_surrgate_over_tau_',num2str(bin_size)); 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
       
       % plot MI surrogate distribution
       figure('units','normalized','outerposition',[0 0 1 1]) 
       hold on 
       surf(MI_surrogate)
       xlabel('tau')
       ylabel('surrogate number');
       zlabel('mutual information');
       colormap('jet');
       cb = colorbar; 
       title(cb,'mutual information');
       view([-45,75,25]);
       title(['Surrogate MI Distribution for: ', current_file...
           ,' Bin Size = ', num2str(bin_size)],'Interpreter','none');
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')

       % save figure 
       plot_name = strcat('surrogate_MI_distribution_surf_',num2str(bin_size)); 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
       
       % plot MI surrogate histogram over tau 
       surrogate_3D_hist_points = cell(3,size(MI_surrogate,2)); 
       surrogate_hist_data = cell(5,size(MI_surrogate,2)); 
       for ii = 1:size(MI_surrogate,2) 
           h = histogram(MI_surrogate(:,ii),512); 
           surrogate_hist_data{1,ii} = h.Data; 
           surrogate_hist_data{2,ii} = h.Values;
           surrogate_hist_data{3,ii} = h.BinWidth; 
           surrogate_hist_data{4,ii} = h.BinLimits;
           % x,y,z = tau,MI,number of occurances 
           y_val = cell(1,512); 
           for jj = 1:512 
               surrogate_3D_hist_points{1,ii} = ii; % x = tau 
               MI_y = h.BinLimits(1,1)+(h.BinWidth.*jj); % y = MI 
               surrogate_3D_hist_points{3,ii} = h.Values; % z = num of occurrances
               y_val{jj} = MI_y;
           end 
           surrogate_3D_hist_points{2,ii} = cell2mat(y_val(1:end));
       end 
       figure('units','normalized','outerposition',[0 0 1 1]) 
       hold on 
       for ii = 1:size(MI_surrogate,2) 
           x = ones(1,512)*ii; 
           y = (surrogate_3D_hist_points{2,ii}); 
           z = log2((surrogate_3D_hist_points{3,ii}));
           scatter3(x,y,z) 
       end 
       hold off 
       xlabel('tau');
       ylabel('MI');
       zlabel('log2(number of occurrances)'); 
       view([-45,15,85])
       title(['MI Versus Tau and Log(MI ocurrences) for: ', current_file....
           ,' Bin Size = ', num2str(bin_size),' Number of Surrogates = ' num2str(num_of_MI_values)],'Interpreter','none');
        set(gca,'FontSize',15)
        set(gca, 'FontName', 'Times New Roman')
       
       % save figure 
       plot_name = strcat('MI_surrogate_hist_distrib_',num2str(bin_size)); 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
       
       % plot significance over tau 
       figure('units','normalized','outerposition',[0 0 1 1])  
       hold on 
       % calculate sigma (SD) of surrogate) 
       significance = ((MI_data - mean(MI_surrogate)))./stan_dev; 
       xlabel('tau'); 
       ylabel('significance'); 
       title(['Significance Versus Tau for File: ', current_file....
           ,' Bin Size = ', num2str(bin_size)],'Interpreter','none');
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')
       hold on
       plot(abs(significance),'b');
       plot(significance,'r')
       plot(ones(1,length(significance)),'--');
       plot(2.*ones(1,length(significance)),'--');
       plot(3.*ones(1,length(significance)),'--');
       hold off 
       legend('ABS(significance)','significance','1 sigma: surrogate','2 sigma: surrogate','3 sigma: surrogate');

       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')
       % save figure 
       plot_name = strcat('significance_ov_tau_',num2str(bin_size)); 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
       
       %% CALCULATE MI OVER TIME (WITH VARIABLE TAU) FOR ACTUAL DATA  
       %{ 
       Steps: 
       1. starting with tau = 1, increment +1 the beginning tau for the MI
       calculation series, and then for each new beginning of the series,
       calculate MI over tau 
       2. For an individual tau_beginning = n: 
            i. bin binary matrix 
            ii. find number of num non zero elements in each bin across
            across all realizations 
            iii. set a threshold (more than 1 spike / bin) for each bin and
            filter and create a binary matrix
 
       %} 
       
       %% MI VERSUS TIME: BIN ACTUAL DATA 
       
       MI_over_time = cell(length(tau_lower_lim:tau_upper_lim-1),1);
       % increment the start of the time series for the MI calculation 
       for ii = tau_lower_lim:tau_upper_lim 
           % bin across all realizations 
           length_realization = length(binary_matrix); 
           num_bins = ceil(length_realization/bin_size); 
           num_nnz_elements = 0; 
           for jj = 1:num_bins-1 % bin data 
               if jj < num_bins && jj ~= 1 
                   current_matrix_bin = binary_matrix(:,(bin_size*(jj-1))+1:bin_size*(jj)); 
                   num_non_zero_elements = nnz(current_matrix_bin); 
                   num_nnz_elements = [num_nnz_elements num_non_zero_elements]; 
               end 
               if jj == 1 
                   current_matrix_bin = binary_matrix(:,1:bin_size*(jj)); 
                   num_non_zero_elements = nnz(current_matrix_bin); 
                   num_nnz_elements = [num_nnz_elements num_non_zero_elements]; 
               end 
               if jj >= length(binary_matrix) - bin_size && jj ~= 1
                   current_matrix_bin = binary_matrix(:,length(binary_matrix)-bin_size:length(binary_matrix));
                   num_non_zero_elements = nnz(current_matrix_bin); 
                   num_nnz_elements = [num_nnz_elements num_non_zero_elements];
               end 
           end % end
           num_nnz_elements = num_nnz_elements(2:end);
           thresh_containing_spike = 5; 
           num_nnz_elements_filtered = 0; 
           % assign 1 to bins with more than 1 spike 
           % assign 0 to bins with 1 spike or less 
           for kk = 1:length(num_nnz_elements) 
               if num_nnz_elements(kk) > thresh_containing_spike
                   num_nnz_elements_filtered = [num_nnz_elements_filtered 1]; 
               else 
                   num_nnz_elements_filtered = [num_nnz_elements_filtered 0];
               end 
           end % end thresholding each bin 
           num_nnz_elements_filtered = (num_nnz_elements_filtered(2:end)); 
           % function call to calculate MI over tau 
           tau_lower = ii; 
           % increment through tau for set beginning tau, calculating MI 
           MI_over_tau_const_tau_beginning = cell(size(num_nnz_elements_filtered,2),1); 
           for gg = tau_lower:tau_upper_lim
              bin_size_for_MI = 1; 
              MI_current = MI_JS(num_nnz_elements_filtered(ii:end),num_nnz_elements_filtered(ii:end),...
                  gg,bin_size_for_MI);
              MI_over_tau_const_tau_beginning{gg} = [MI_over_tau_const_tau_beginning{gg}; MI_current]; 
           end % end cycling through variable tau, constant beginning tau 
           MI_over_time{ii} = cell2mat(MI_over_tau_const_tau_beginning);     
       end % end cycling through tau beginning 

       %% PLOT OF MI OVER TIME (WITH VARIABLE TAU) 
       figure('units','normalized','outerposition',[0 0 1 1])
       hold on 
       for ii = 1:length(MI_over_time) 
           tau = 1:length(MI_over_time{ii}); 
           T_iith = ii.*ones(length(MI_over_time{ii}),1);
           MI_iith = MI_over_time{ii}; 
           plot3(T_iith,tau,MI_iith);
       end 
       hold off 
       view([-135 165 135]);
       xlabel('time (bin number)');
       ylabel('tau (look ahead)');
       zlabel('mutual information');
       title(['MI Versus Tau and Time: ', current_file....
           ,' Bin Size = ', num2str(bin_size), ' Contains Spike Threshold = ' ,num2str(thresh_containing_spike)],'Interpreter','none');
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')
       % save figure 
       plot_name = strcat('MI_ov_Time_Tau_',num2str(bin_size)); 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
       
       figure('units','normalized','outerposition',[0 0 1 1]) 
       hold on 
       for ii = 1:length(MI_over_time) 
           tau = 1:length(MI_over_time{ii}); 
           MI_iith = MI_over_time{ii}; 
           T_iith = ii.*ones(length(MI_over_time{ii}),1);
           plot(tau,MI_iith);
       end 
       hold off 
       xlabel('tau (look ahead)');
       ylabel('mutual information');
       title(['MI Versus Tau For All Time Spike Peaks Not Aligned: ', current_file....
           ,' Bin Size = ', num2str(bin_size), ' Contains Spike Threshold = ' ,num2str(thresh_containing_spike)],'Interpreter','none');
       
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')
       % save figure 
       plot_name = strcat('MI_Over_Tau_Time_Spikes_Not_Aligned',num2str(bin_size)); 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
       
       MI_over_time_aligned = MI_over_time.';
       figure('units','normalized','outerposition',[0 0 1 1]) 
       hold on 
       for ii = 1: 10% length(MI_over_time) 
           %tau = length(MI_over_time{ii}):-1:1;
           tau = 1:+1:length(MI_over_time{ii});
           MI_iith = MI_over_time{ii}; 
           MI_iith_reversed = flipud(MI_iith);
           plot(tau,MI_iith);
       end 
       hold off 
       % ax = gca; 
       % ax.XDir = 'reverse';
       xlabel('tau (look ahead)');
       ylabel('mutual information');
       title(['MI Versus Tau With Spike Peaks Aligned: ', current_file....
           ,' Bin Size = ', num2str(bin_size), ' Contains Spike Threshold = ' ,num2str(thresh_containing_spike)],'Interpreter','none');
       set(gca,'FontSize',15)
       set(gca, 'FontName', 'Times New Roman')
       % save figure 
       plot_name = strcat('MI_Over_Tau_Time_Spikes_Aligned',num2str(bin_size)); 
       figure_name = strcat(new_folder_path,'\','plot_',plot_name);
       print(figure_name,'-djpeg')
        
       % end % end cycling through bin_size
       % end % end cycling through multiple files 
       %{
        %% PLOT MI OVER VAU FOR MULTIPLE DATA FILES, BIN SIZE = CONSTANT 
        figure('units','normalized','outerposition',[0 0 1 1])
        for ii = 1:uu
            % get values to plot 
            tau_files = store_MI_data_surrogate_over_tau_mult_data_files{4,ii}; 
            MI_data_files = store_MI_data_surrogate_over_tau_mult_data_files{1,ii}; 
            stan_dev_files = store_MI_data_surrogate_over_tau_mult_data_files{3,ii}; 
            MS_mean_files = store_MI_data_surrogate_over_tau_mult_data_files{2,ii};
            % plot values 
            subplot(uu,1,ii);
            hold on
            plot(tau_files,MI_data_files,'r'); 
            errorbar(tau_files,MS_mean_files,stan_dev_files);
            plot(tau_files,MS_mean_files,'b');
            legend('data','sigma of surrogate','mean of surrogate','Location','bestoutside');
            title(files_loop_thru{ii},'Interpreter','none');
            axis([0 120 0 max(MI_data_files)]);
            hold off

        end 
        xlabel('tau (look ahead'); 
        ylabel('Mutual Information'); 
       
        % save figure 
        plot_name = 'MI_vs_Tau_Multiple_Files'; 
        print(plot_name,'-djpeg')
        %}
        %% END PROGRAM 
        
       disp('program run time:');        
       toc; % stop timer 
       

      
                     
        %% SUB FUNCTIONS 
        %{ 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %} 
        %% HELPER FUNCTION 1: CALCULATE MI OVER TAU FOR INPUT OF 2D OR 1D ARRAY 
        %{ 
        Description: 
        Inputs:
            1. data: where data is a 2D array containing the series for which
            to calculate MI for 
                dimensions: length of series (rows) x number of series
                (columns) 
            2. tau_lower_limit: the lower look ahead value 
            3. tau_upper_limt: the upper look ahead value 
        Outputs: 
        %} 
        function [MI] = MI_over_Tau(data,tau_lower_lim,tau_upper_lim) 
            % increment through all tau, calculating MI 
            %{
            DISCLAIMER: !!!MI is dependent on length of data!!!
            NOTE THAT MI CHANGES AS A FUNCTION OF THE LENGTH OF THE
            INPUT DATA SERIES; THUS MI WILL BE LARGER FOR LONGER LENGTH OF DATA
            SERIES AND SMALLER FOR SHORTER SERIES OF SERIES 
            %}
            tau_dimension = tau_upper_lim - tau_lower_lim; 
            number_of_series = size(data,2); 
            MI = cell(number_of_series+1,tau_dimension+1); % cell to store MI calculation results 
            for xx = tau_lower_lim:tau_upper_lim 
                if size(data,2) > 1 % for surrogate 
                    % calculate MI across all columns (where there are
                    % multiples series, as in the case of multiple
                    % surrogate files)
                    for yy = 1:size(data,2) 
                        tau_look_ahead = xx; % look ahead is simply the variable of incrementing 
                        bin_size_for_MI = 1; % set bin size = 1 so that number of bins = ceil((max(data) - min(data))/ 1 
                        MI_current = MI_JS(data(:,yy),data(:,yy),tau_look_ahead,bin_size_for_MI);
                        MI{yy,xx} = [MI{yy,xx}; MI_current]; 
                    end % end looping through all series for constant tau 
                end % end if statement 
                if size(data,2) == 1 % for data 
                    % data is 1 column; calculate MI for just 1 column 
                    tau_look_ahead = xx; % look ahead is simply the variable of incrementing 
                    bin_size_for_MI = 1; % set bin size = 1 so that number of bins = ceil((max(data) - min(data))/ 1 
                    % function call to calculate MI 
                    MI_current = MI_JS(data(:,1),data(:,1),tau_look_ahead,bin_size_for_MI);
                    MI{1,xx} = [MI{1,xx}; MI_current];
                end % end if statement 
            end % end for loop incrementing through tau 
       end % end function 1   
        
       %%
       % end of script 
 
       

        
        
        
        
           
           
           
       
       
       
       
       
       
       
       