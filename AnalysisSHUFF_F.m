function STDEV = AnalysisSHUFF_F(directions, raw_position, range, n_shuffles)

samp=2;
raw = raw_position;
dir = directions;
time = -50:samp:300;
t1=0; t2=50; % window for zero shift
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
EYE = zeros(1,n_shuffles);
for shuffle=1:n_shuffles
    dir = Shuffle(dir);
    % ROC analysis
    [~,~,~,EYE(shuffle)] = perfcurve(int2str(dir'),(dataF(:,stop_index)-dataF(:,start_index)),'2');
    if mod(shuffle,100)==0
        fprintf(1,'.');
    end
    if mod(shuffle,500)==0
        fprintf(1,' ');
    end
end
fprintf(1,'\n');
STDEV = std(EYE);