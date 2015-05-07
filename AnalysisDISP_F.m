function [PSY EYE n_trials time left_pos right_pos] = AnalysisDISP_F(directions, results, raw_position, range) 
% directions = doi_directions;
% results = doi_results;
% raw_position = doi_raw_position;
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

% plot individual trials
k1=1; k2=1;
for i=1:n_trials
    if dir(i) == 1
        if plot_each
            subplot(1,2,1)
            plot(time,data(i,:),'r'); hold on
            plot(time,dataF(i,:),'r'); hold off
            ylim([mn mx])
            subplot(1,2,2)
            plot(time(1:end-1),vel(i,:),'r');
        end
        left(k1,:) = data(i,:); leftF(k1,:) = dataF(i,:); LV(k1,:)= vel(i,:); k1 = k1+1;
    else
        if plot_each
            subplot(1,2,1)
            plot(time,data(i,:),'b'); hold on
            plot(time,dataF(i,:),'b'); hold off
            ylim([mn mx])
            subplot(1,2,2)
            plot(time(1:end-1),vel(i,:),'b');
        end
        right(k2,:) = data(i,:);  rightF(k2,:) = dataF(i,:);RV(k1,:)= vel(i,:);  k2 = k2+1;
    end
    if plot_each
        pause
    end
end

% % plot means
% figure
% subplot(1,2,1)
% plot(time,mean(left),'r'); hold on; 
% plot(time,mean(right),'b'); plot(time,mean(leftF),'r'); plot(time,mean(rightF),'b')
% subplot(1,2,2)
% plot(time(1:end-1),mean(LV),'r'); hold on
% plot(time(1:end-1),mean(RV),'b')
% figure;
% plot(time,mean(rightF)-mean(leftF));
% figure;
% plot(time(2:length(time)),diff(mean(rightF)-mean(leftF))/(1/500));

% ROC analysis
    [Xr,Yr,T,AUC(i)] = perfcurve(int2str(dir'),(dataF(:,stop_index)-dataF(:,start_index)),'2');
% figure
% plot(time(1:length(AUC)),AUC,'b'); hold on
% plot(time(1:length(AUCv)),AUCv,'r');
%[max(AUC) max(AUCv)]
PSY = mean(results); 
EYE = max(AUC);
left_pos = nanmean(leftF);
right_pos = nanmean(rightF);
left_sem = nanstd(leftF)/sqrt(50);
right_sem = nanstd(rightF)/sqrt(50);
% displacements = dataF(:,stop_index)-dataF(:,start_index);
% displacements = displacements';
% lefts = displacements(dir==1);
% rights = displacements(dir==2);
% plot(rand(1,length(lefts)),lefts,'bo')
% hold on
% plot(rand(1,length(rights)),rights,'ro')
% nanmean(lefts)
% nanmean(rights)
% nanstd(lefts)
% nanstd(rights)