function [EYEmedian PROPagree] = AnalysisDISP_Fmedian(directions, results, raw_position, range) 

clear left right
samp=2;
raw = raw_position;
dir = directions;
time = -50:samp:range(2);
plot_each = 0; %whether to cycle through trials one by one
t1=0; t2=50; % wondow for zero shift
start_index = find(time==range(1));
stop_index = find(time==range(2));
windowP = 0;  % window for position (0 seems to be the best)
windowV = 50;  % window for speed (larger windows are better - 300 is the best) 
n_trials = size(raw_position,1);

% zero shift
for i=1:n_trials
    data(i,:) = raw(i,:)-mean(raw(i,t1+26:t2+26));
end

% filtering
filterWidth = 19;
data1=sgolayfilt(data',3,filterWidth)';  % double filtering
dataF=sgolayfilt(data1',3,filterWidth)';
vel = (diff(dataF')/(0.002))';
mx = max(max(data))*.8;mn = min(min(data))*.8;

% median analysis
nums = (dataF(:,stop_index)-dataF(:,start_index))';
criterion = nanmedian(nums);
lefts = nums(dir==1);
resultsL = results(dir==1);
resultsR = results(dir'==2);
rights = nums(dir'==2);
EYEmedian = (sum(lefts<criterion)+sum(rights>=criterion))/length(nums);
PROPagree = (sum((lefts<criterion)==resultsL)+sum((rights>=criterion)==resultsR))/length(nums);