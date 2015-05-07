function [pAgree agreeCI pIndep pCeil nSig nCriteria] = AnalysisDISP_Fmedian_ALL(directions, results, raw_position, range)
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

% median analysis
nums = (dataF(:,stop_index)-dataF(:,start_index))';
snip = isnan(nums);
nums(snip)=[];
dir(snip)=[];
results(snip)=[];
n_trials = n_trials-sum(snip);
[sortednums,ranknums] = sort(nums);
scatter_diff = zeros(length(nums),2);
trialno = 1:n_trials;
for criterion=1:n_trials
    eyeresp(nums<nums(ranknums(criterion)))=1;
    eyeresp(nums>=nums(ranknums(criterion)))=2;
    psyresp(results==1)=dir(results==1);
    psyresp(results==0)=3-dir(results==0);
    eyeresults = eyeresp==dir;
    EYEaccuracy(criterion) = mean(eyeresults);
    PROPagree(criterion) = mean(eyeresults==results);
    
%     ind_bound(criterion) = (mean(results).*EYEaccuracy(criterion)+(1-mean(results)).*(1-EYEaccuracy(criterion)));
%     corr_bound(criterion) = 1-abs(EYEaccuracy(criterion)-mean(results));
%     anticorr_bound(criterion) = 1-(1-EYEaccuracy(criterion))-(1-mean(results));
%     
%     r_ind_bound(criterion) = sum(eyeresp==2)/n_trials*sum(psyresp==2)/n_trials+sum(eyeresp==1)/n_trials*sum(psyresp==1)/n_trials;
%     r_corr_bound(criterion) = (min([sum(eyeresp==2) sum(psyresp==2)])+min([sum(eyeresp==1) sum(psyresp==1)]))/n_trials;
    
    %actual direction, modality, response, all of these are rates
    nL = sum(dir==1);
    numlefts(criterion)=sum(eyeresp==1);
    L_EYE_L = sum(eyeresp(dir==1)==1)/nL;
    L_PSY_L = sum(psyresp(dir==1)==1)/nL;
    L_EYE_R = sum(eyeresp(dir==1)==2)/nL;
    L_PSY_R = sum(psyresp(dir==1)==2)/nL;
    
    nR = sum(dir==2);
    R_EYE_L = sum(eyeresp(dir==2)==1)/nR;
    R_PSY_L = sum(psyresp(dir==2)==1)/nR;
    R_EYE_R = sum(eyeresp(dir==2)==2)/nR;
    R_PSY_R = sum(psyresp(dir==2)==2)/nR;
    
    sk_ind_bound(criterion) = ((L_EYE_L*L_PSY_L+L_EYE_R*L_PSY_R)*nL+(R_EYE_L*R_PSY_L+R_EYE_R*R_PSY_R)*nR)/n_trials;
    sk_corr_bound(criterion) = ((min([L_EYE_L L_PSY_L])+min([L_EYE_R L_PSY_R]))*nL+(min([R_EYE_L R_PSY_L])+min([R_EYE_R R_PSY_R]))*nR)/n_trials;
    
    sk_CI(criterion) = 1.96*sqrt(PROPagree(criterion)*(1-PROPagree(criterion))/n_trials);
end
critlist = 1:length(nums);
% figure
% subplot(2,1,1)
% h=plot(critlist,EYEaccuracy,'o',critlist,mean(results)*ones(size(critlist)),'o');
% set(h(1),'MarkerEdgeColor','none','MarkerFaceColor','b');
% set(h(2),'MarkerEdgeColor','none','MarkerFaceColor','r');
% xlabel('Number of LEFTWARD responses');
% ylabel('Proportion Correct')
% title('SB 110ms, Red=PSY, Blue=EYE')
% xlim([0 n_trials+1]);
% subplot(2,1,2)
% h=plot(critlist,sk_corr_bound,'o',critlist,sk_ind_bound,'o');
% set(h(1),'MarkerEdgeColor','none','MarkerFaceColor','r');
% set(h(2),'MarkerEdgeColor','none','MarkerFaceColor','b');
% hold on
% h = errorbar(critlist,PROPagree,sk_CI,'go');
% set(h(1),'MarkerEdgeColor','none','MarkerFaceColor','g');
% xlabel('Criterion / Number of LEFTWARD responses');
% ylabel('Proportion Agreed')
% %title('SB 110ms, Green=Proportion Agree, Blue=Independent, Red=Maximal Correlation')
% xlim([0 n_trials+1]);

% %Set output vars
pAgree = PROPagree(nL+1);
agreeCI = sk_CI(nL+1);
pIndep = sk_ind_bound(nL+1);
pCeil = sk_corr_bound(nL+1);
nSig = sum(((PROPagree+sk_CI)>=sk_ind_bound).*(sk_ind_bound>=(PROPagree-sk_CI)));
nCriteria = n_trials;