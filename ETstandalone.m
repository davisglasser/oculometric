clear; clc;
tic;
files = {'Data/SB/20121002/SB1529r'
    'Data/SB/20121002/SB1539r'
    'Data/SB/20121002/SB1546r'
    'Data/SB/20121002/SB1553r'
    'Data/SB/20121003/SB1512r'
    'Data/SB/20121003/SB1518r'
    'Data/SB/20121003/SB1544r'
    'Data/SB/20121003/SB1551r'
    'Data/SB/20121004/SB1346r'
    'Data/SB/20121004/SB1356r'
    'Data/SB/20121004/SB145r'};

dur_list = [50 60 70 80 90 110 125 150 175 200];
range = [0 300];
plot_stuff = 0;
n_shuffles = 00;

for i=1:length(files)
    [raw_position directions samp results durations RTs goofs(i)] = ET_readerF(files{i});
    % concatenate things
    if i==1
        directions1 = directions;
        raw_position1 = raw_position;
        durations1 = durations;
        results1 = results;
        RTs1 = RTs;
    else
        if size(raw_position1,2)>size(raw_position,2)
            raw_position1 = cat(1,raw_position1(:,1:size(raw_position,2)),raw_position);
        else
            raw_position1 = cat(1,raw_position1,raw_position(:,1:size(raw_position1,2)));
        end
        directions1 = cat(2,directions1,directions);
        results1 = cat(2,results1,results);
        durations1 = cat(2,durations1,durations);
        RTs1 = cat(2,RTs1,RTs);
    end
end
raw_position = raw_position1;
directions = directions1;
durations = durations1;
results = results1;
RTs = RTs1;
PSY = zeros(1,length(dur_list));
EYE = PSY;
n_trials = PSY;
rEYE = PSY;
resEYE = PSY;
if ~isempty(dur_list)
    for i=1:length(dur_list)
        [doi_directions doi_results doi_raw_position doi_RTs] = ET_doiF(dur_list(i), durations, raw_position, results, directions, RTs);
        dur_RTs(i,1)=mean(doi_RTs);
        dur_RTs(i,2)=std(double(doi_RTs));
        dur_RTs(i,3)=median(doi_RTs);
        [PSY(i) EYE(i) n_trials(i) time{i} left{i} right{i}] = AnalysisDISP_F(doi_directions, doi_results, doi_raw_position,range);
        [EYEmedian(i) AGREEmedian(i)]=AnalysisDISP_Fmedian(doi_directions, doi_results, doi_raw_position,range);
        [pAgree(i) agreeCI(i) pIndep(i) pCeil(i) nSig(i) nCriteria(i)] = AnalysisDISP_Fmedian_ALL(doi_directions, doi_results, doi_raw_position,range);
        resp_directions = zeros(1,length(doi_results));
        for j=1:length(doi_results)
            if doi_directions(j)==1
                if doi_results(j)==1
                    resp_directions(j)=1;
                else
                    resp_directions(j)=2;
                end
            else
                if doi_results(j)==1
                    resp_directions(j)=2;
                else
                    resp_directions(j)=1;
                end
            end
        end
        rEYE(i) = AnalysisDISP_F_EYEONLY(resp_directions, doi_results, doi_raw_position,range); %Analysis based on responses
        CI(i)=sqrt(PSY(i)*(1-PSY(i))/n_trials(i));
        if n_shuffles
            STDEV(i) = AnalysisSHUFF_F(doi_directions, doi_raw_position,range,n_shuffles);
            rSTDEV(i) = AnalysisSHUFF_F(resp_directions, doi_raw_position,range,n_shuffles);
            %[dur_list(i) PSY(i) EYE(i) STDEV(i) n_trials(i)]
        else
            %[dur_list(i) PSY(i) EYE(i) n_trials(i)]
        end
        
        start_index = find(time{1}==range(1));
        stop_index = find(time{1}==range(2));
        zeroed = doi_raw_position(:,start_index:stop_index)-repmat(doi_raw_position(:,start_index),1,(stop_index-start_index)+1);
        mean_difference{i}(1:(stop_index-start_index+1)) = nanmean(zeroed(doi_directions==2,:))-nanmean(zeroed(doi_directions==1,:));
        stdev_difference{i}(1:(stop_index-start_index+1)) = sqrt(nanstd(zeroed(doi_directions==1,:)).^2+nanstd(zeroed(doi_directions==2,:)).^2);
        clear lefties righties
        trialindex = 1:size(zeroed,1);
        lefties(:,1) = trialindex(doi_directions==1);
        lefties(:,2) = zeroed(doi_directions==1,size(zeroed,2));
        righties(:,1) = trialindex(doi_directions==2);
        righties(:,2) = zeroed(doi_directions==2,size(zeroed,2));
        
    end
    if plot_stuff
        plot(dur_list,PSY,'-ob',dur_list,EYE,'-or')
        figure
        plot(time{1},right{1}-left{1})
    end
    if n_shuffles
        output = [dur_list' PSY' CI' EYE' STDEV' pAgree' pIndep' pCeil']
    else
        %output = [dur_list' PSY' EYE' rEYE' n_trials']
        output = [dur_list' PSY' CI' EYE' pAgree' pIndep' pCeil']
    end
    datP = [dur_list' PSY' n_trials'];
    datE = [dur_list' EYE' n_trials'];
else
    [PSY EYE n_trials time left right] = AnalysisDISP_F(directions, results, raw_position,range);
    %[PSY EYE n_trials time left right] = AnalysisPEAK_dgF(directions, results, raw_position);
    [PSY EYE n_trials]
end
RT_output = [mean(double(RTs)) std(double(RTs)) median(double(RTs))];
%crit_output = [dur_list' pAgree' agreeCI' pIndep' pCeil' nCriteria'-nSig' nCriteria']
ElapsedTime = toc/60