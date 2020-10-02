function mutual = Mutual_Information_JS(a,b,k,binsize)
% INPUT
%   a: first series
%   b: second series (currently only set up for being identical to 'a')
%   k: look ahead
%   binsize: the size of the bins that you want, in same units as 'a'

% Script credit: to Jesse Snelling for originally creating; with Greg Zdor
% edits

%

nbins = ceil((max(a) - min(a))/binsize);
 
% disp(['nbins = ' num2str(nbins)]);

%Define lists
ak = a(1:end-k);
bk = b(1+k:end);
N = numel(ak); %number of events

is = ceil(ak/binsize);
js = ceil(bk/binsize);
ks = is + nbins*js;

uis = unique(is);
proba = histc(is(:),uis)/N;

ujs = unique(js);
probb = histc(js(:),ujs)/N;

uks = unique(ks);
probab = histc(ks(:),uks)/N;

ha = entropy_JS(proba);
hb = entropy_JS(probb);
hab = entropy_JS(probab);

mutual = ha + hb - hab;


