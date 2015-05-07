function [raw_position directions samp results durations RTs goofs] = ET_readerF(filename)

% Read ASC conversions of EDF data for the OFR.
%
% Davis Glasser
% Last Updated: 7/18/2012

if nargin==0
filename = 'Data/PR/20120821/PR1532r';                                                     % File to be analyzed
end
samp = 500;                                                                 % Hz

goofs = 0;
trial_no = 0;                                                               % Set current trial to 0
fid = fopen(strcat(filename,'e.asc'));                                      % Open the events file
blink_flag=0;
fix_flag=0;
sac_flag=0;
samp = 1000/samp;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end                                           % If the file is over, stop
    
    % Things to check for:
    if strfind(tline,'START')                                               % START (Set #saccades/blinks to 0)
        C = textscan(tline,'%s %d %s %s');
        trial_no = trial_no+1;
        data(trial_no).start = C{2};
        data(trial_no).sac_count = 0;
        data(trial_no).blink_count = 0;
        data(trial_no).fix_count = 0;
    end
    
    if strfind(tline,'TRIAL')                                           % TRIAL (make sure it matches the expected # from START)
        if ~strfind(tline,strcat('TRIAL_',num2str(trial_no),'_'))
            fprintf('Check out Trial %3.0f!!! \n',trial_no);
            indices = strfind(tline,'_');
            trial_no = str2double(tline((indices(1))+1:(indices(2)-1)));
        end
    end
    
    if strfind(tline,'DURATION')                                        % DURATION
        indices = strfind(tline,'_');
        data(trial_no).duration = str2double(tline((indices(1))+1:(indices(2)-1)));
    end
    
    if strfind(tline,'FIX_')                                            % FIX
        indices = strfind(tline,'_');
        data(trial_no).fix = tline((indices(1))+1:(indices(2)-1));
    end
    
    if strfind(tline,'DIRECTION')                                       % DIRECTION
        if strfind(tline,'LEFT')
            data(trial_no).direction = 1;
        else
            data(trial_no).direction = 2;
        end
    end
    
    if strfind(tline,'FIXATION_ARRIVED')                                % FIX ARRIVED
        C = textscan(tline,'%s %d %s');
        data(trial_no).fix_arrived = C{2};
    end
    
    if strfind(tline,'SACCADE_TARGET_ONSET')                            % SACCADE TARGET ONSET
        C = textscan(tline,'%s %d %s');
        data(trial_no).sac_target_on = C{2};
    end
    
    if strfind(tline,'STIMULUS_TARGET_ARRIVED')                         % SACCADE TARGET ARRIVED
        C = textscan(tline,'%s %d %s');
        data(trial_no).stim_target_arrived = C{2};
    end
    
    if strfind(tline,'STIMULUS_ONSET')                                  % STIMULUS ONSET
        C = textscan(tline,'%s %d %s');
        data(trial_no).stim_onset = C{2};
    end
    
    if strfind(tline,'STIMULUS_OFFSET')                                 % STIMULUS OFFSET
        C = textscan(tline,'%s %d %s');
        data(trial_no).stim_offset = C{2};
        data(trial_no).actual_duration = double(data(trial_no).stim_offset-data(trial_no).stim_onset);
    end
    
    if strfind(tline,'RESPONSE')                                        % RESPONSE
        C = textscan(tline,'%s %d %s');
        data(trial_no).response_time = C{2};
        if strfind(tline,'_CORRECT')
            data(trial_no).result = 1;
        else
            data(trial_no).result = 0;
        end
        if fix_flag
            data(trial_no).fixations(data(trial_no).fix_count).offset = C{2};
        end
        if blink_flag
            data(trial_no).blinks(data(trial_no).blink_count).offset = C{2};
        end
        if sac_flag
            data(trial_no).saccades(data(trial_no).sac_count).offset = C{2};
        end
    end
    
    if strfind(tline,'END')  
    if strcmp(tline(1:3),'END')
        C = textscan(tline,'%s %d %s %s %d %d');
        data(trial_no).end = C{2};
    end
    end
    
    if strfind(tline,'SFIX')                                            % SFIX
        fix_flag = 1;
        data(trial_no).fix_count = data(trial_no).fix_count+1;
        C = textscan(tline,'%s %s %d');
        data(trial_no).fixations(data(trial_no).fix_count).onset = C{3};
    end
    
    if strfind(tline,'EFIX')                                            % EFIX
        if data(trial_no).fix_count == 0, data(trial_no).fix_count=1; end
        fix_flag = 0;
        C = textscan(tline,'%s %s %d %d %d %d %d %d');
        data(trial_no).fixations(data(trial_no).fix_count).onset = C{3};
        data(trial_no).fixations(data(trial_no).fix_count).offset = C{4};
        data(trial_no).fixations(data(trial_no).fix_count).duration = C{5};
    end
    
    if strfind(tline,'SBLINK')                                          % SBLINK
        blink_flag = 1;
        data(trial_no).blink_count = data(trial_no).blink_count+1;
        C = textscan(tline,'%s %s %d');
        data(trial_no).blinks(data(trial_no).blink_count).onset = C{3};
    end
    
    if strfind(tline,'EBLINK')                                          % EBLINK
        if data(trial_no).blink_count == 0, data(trial_no).blink_count=1; end 
        blink_flag = 0;
        C = textscan(tline,'%s %s %d %d %d');
        data(trial_no).blinks(data(trial_no).blink_count).onset = C{3};
        data(trial_no).blinks(data(trial_no).blink_count).offset = C{4};
        data(trial_no).blinks(data(trial_no).blink_count).duration = C{5};
    end
    
    if strfind(tline,'SSACC')                                           % SSACC
        C = textscan(tline,'%s %s %d');
        sac_flag = 1;
        data(trial_no).sac_count = data(trial_no).sac_count+1;
        data(trial_no).saccades(data(trial_no).sac_count).onset=C{3};
    end
    
    if strfind(tline,'ESACC')                                           % ESACC
        if data(trial_no).sac_count == 0, data(trial_no).sac_count=1; end                                 % If the start of the saccade was before the trial started
        sac_flag = 0;
        C = textscan(tline,'%s %s %d %d %d %f %f %f %f %f %d');
        data(trial_no).saccades(data(trial_no).sac_count).onset=C{3};
        data(trial_no).saccades(data(trial_no).sac_count).offset=C{4};
        data(trial_no).saccades(data(trial_no).sac_count).duration=C{5};
        %data(trial_no).saccades(data(trial_no).sac_count).start_pos(1:2)=[C{6} C{7}];
        %data(trial_no).saccades(data(trial_no).sac_count).stop_pos(1:2)=[C{8} C{9}];
        data(trial_no).saccades(data(trial_no).sac_count).amplitude = C{10};
        data(trial_no).saccades(data(trial_no).sac_count).peak_vel = C{11};
    end
end % End while loop

fclose(fid); % close events file

fid = fopen(strcat(filename,'s.asc'));
samples_cell  = textscan(fid,'%f %f %f %f','CollectOutput', 1,'TreatAsEmpty','.');
samples = samples_cell{1};
fclose(fid);% Analyze the samples

%Assign samples to trials, make fix, sacc, blink flags
for i=1:trial_no
    [y start_index] = min(abs(samples(:,1)-double(data(i).start)-.1));
    %start_index = find(samples(:,1)==data(i).start);
    %stop_index = find(samples(:,1)==data(i).end);
    [y stop_index] = min(abs(samples(:,1)-double(data(i).end)-.1));
%     while isempty(start_index)
%         start_index = find(samples(:,1)==data(i).start+1);
%     end
%     while isempty(stop_index)
%         stop_index = find(samples(:,1)==data(i).end-1);
%     end
    data(i).time = samples(start_index:stop_index,1)';
    data(i).raw_x = samples(start_index:stop_index,2)';
    data(i).raw_y = samples(start_index:stop_index,3)';
    data(i).missing_flag = isnan(data(i).raw_x);
    data(i).fix_flag = zeros(size(data(i).time));
    data(i).sac_flag = zeros(size(data(i).time));
    data(i).blink_flag = zeros(size(data(i).time));
    data(i).RT = data(i).response_time-data(i).stim_onset;
    for j=1:data(i).sac_count
        start_index = find(data(i).time==data(i).saccades(j).onset);
        stop_index = find(data(i).time==data(i).saccades(j).offset);
        data(i).sac_flag(start_index:stop_index) = 1;
    end
    for j=1:data(i).blink_count
        start_index = find(data(i).time==data(i).blinks(j).onset);
        if ~isempty(data(i).blinks(j).offset)
            stop_index = find(data(i).time==data(i).blinks(j).offset);
        else
            stop_index = length(data(i).time);
        end
        data(i).blink_flag(start_index:stop_index) = 1;
    end
    data(i).fix_flag = (1-data(i).blink_flag).*(1-data(i).sac_flag);
    data(i).raw_x_deg = (data(i).raw_x-640)/30;
    data(i).raw_y_deg = (data(i).raw_y-360)/30;
    [y data(i).stim_onset_index] = min(abs(data(i).time-double(data(i).stim_onset)-.1));
%     while isempty(find(data(i).time==data(i).stim_onset,1))
%         data(i).stim_onset=data(i).stim_onset+1;
%     end
%     data(i).stim_onset_index = find(data(i).time==data(i).stim_onset);
    if data(i).actual_duration<(ceil(data(i).duration/8.333333)*8.33333+16);
        if (length(data(i).fix_flag)-data(i).stim_onset_index)>200
            if mean(data(i).fix_flag(data(i).stim_onset_index:data(i).stim_onset_index+200/samp))==1
                data(i).OK = 1;
            else
                goofs = goofs+1;
                data(i).OK = 0;
            end
        else
            data(i).OK = 0;
        end
    else
        data(i).OK = 0;
    end
end
count=1;
for i=1:trial_no
    if data(i).OK
        if (length(data(i).raw_x_deg)-data(i).stim_onset_index)>175
            raw_position(count,1:(length(data(i).raw_x_deg)-data(i).stim_onset_index+26))=data(i).raw_x_deg(data(i).stim_onset_index-25:length(data(i).raw_x_deg));
            lengths(count) = (length(data(i).raw_x_deg)-data(i).stim_onset_index+26);
            directions(count)=data(i).direction;
            durations(count)=data(i).duration;
            results(count)=data(i).result;
            real_durations(count)=data(i).actual_duration;
            RTs(count)=data(i).RT;
            count = count+1;
        end
    end
end
raw_position = raw_position(:,1:min(lengths));