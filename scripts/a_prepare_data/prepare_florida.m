

%% description
% make sure that MFFMatlabIO plugin is installed!

% The beginning of the recording and the intervals between blocks 
% contain substantial noise and must be removed!! - check code below

% EGI system, 129 is reference
% no ET

% M1 (left mastoid) - Channel 57
% M2 (right mastoid) - Channel 100
% VEOGU - 25 (left above)
% VEOGL - 127 (left below)
% HEOGL - 32
% HEOGR - 125

% P3 - 52
% P7 - 58
% O1 - 70

% P4 - 92
% P8 - 96
% 02 - 83


% EEG data for subject 25 is missing
% subject 33: data for only 2 blocs and eye task. resting missing

% subject 10, block 1: 1 epoch is shorter ('2025-06-12T14:07:43.769000-04:00') - look for boundary event
% epoching [-0.7, 3] to recompute responens will remove this epoch, but 
% epoching [-0.7, 1] later during CDA computation will not, which will result 
% in missmatch between response length and number of epoch.
% here, i will add a row to beh file to accoount for this mismatch

%%
d = dir(fullfile(rawdatapath, 'data_florida/'));
d = d(~contains({d.name}, {'.', '..', '.DS'}));
result_folder = fullfile(formatted_data, 'florida');

% load and merge cda files
for sub = 1 : size(d, 1)
    
    d_sub_cda = dir(strcat(fullfile(d(sub).folder, d(sub).name), '/*cda*mff'));
    d_sub_cda_beh = dir(strcat(fullfile(d(sub).folder, d(sub).name), '/*cda*mat'));

    subjectID = d(sub).name;
    
    % dont continue, if EEG data is missing
    if isempty(d_sub_cda)
        continue;
    end

    % prepare to concat CDA data
    EEGconcat = {};
    BEHconcat = {};
    mkdir(fullfile(result_folder, subjectID))

    for f = 1 : size(d_sub_cda, 1)
        
        % load EEG
        mff_path = fullfile(d_sub_cda(f).folder, d_sub_cda(f).name);
        EEG = mff_import(mff_path);

        if strcmp(subjectID, '1') & f == 1 % remove becuase there is no cue information
            EEG.event(2) = [];
        end

        % trim whitespaces
        for i = 1:length(EEG.event)
            EEG.event(i).type = strtrim(EEG.event(i).type);
        end
           
        % load beh
        % subject 3 doesnt have beh file for 4th block. some eeg triggers are missing as well. recompute the entire beh file later.
        % subject 5 doesnt have beh file for 2th block. some eeg triggers are missing as well. recompute the entire beh file later.
        if not((strcmp(subjectID, '3') & f == 4) | (strcmp(subjectID, '5') & f == 2))
            load(fullfile(d_sub_cda_beh(f).folder, d_sub_cda_beh(f).name)); % behavioral
        end

        % check events length
        i_events = find(ismember({EEG.event.type}, {'21', '41', '61'}));
        
        % recompute beh data if some events are missing
        if length(i_events) < 144

            % first, we need longer segment than contains cues and responses
            tmpEEG = pop_epoch(EEG, {21, 41, 61}, [-0.7 3]);
            events = {tmpEEG.event.type};
            epochs = {tmpEEG.event.epoch};

            % init beh data
            allResponses = nan(tmpEEG.trials, 1); % 76, 65
            tmpResponses = nan(tmpEEG.trials, 1);
            allCorrect = nan(tmpEEG.trials, 1); % 0, 1
            trialSetSize = nan(1, tmpEEG.trials); % 2, 4, 6
            trialIfChange = nan(1, tmpEEG.trials); % 0, 1
            trialCuedSide = nan(1, tmpEEG.trials); % 0, 1

            % extract data 
            ev_onset = events(ismember(events, '21') | ismember(events, '41') |ismember(events, '61'));
            trialSetSize(ismember(ev_onset, '21')) = 2;
            trialSetSize(ismember(ev_onset, '41')) = 4;
            trialSetSize(ismember(ev_onset, '61')) = 6;

            ev_cue = events(ismember(events, '3') | ismember(events, '7'));
            trialCuedSide(ismember(ev_cue, '3')) = 0; % left
            trialCuedSide(ismember(ev_cue, '7')) = 1; % right

            ev_resp = events(ismember(events, '76') | ismember(events, '77') | ismember(events, '78') | ismember(events, '79'));
            epochs_of_resp = cell2mat(epochs(ismember(events, '76') | ismember(events, '77') | ismember(events, '78') | ismember(events, '79')));
            tmpResponses(epochs_of_resp) = str2double(ev_resp);
            allCorrect(tmpResponses == 76 | tmpResponses == 77) = 1;
            allCorrect(tmpResponses == 78 | tmpResponses == 79) = 0;
            
            trialIfChange(tmpResponses == 76 | tmpResponses == 79) = 0; % no change
            trialIfChange(tmpResponses == 77 | tmpResponses == 78) = 1; % change

            id2 = str2double(subjectID);
            if mod(id2,2) == 0 
                % L is change, A is no change.
                allResponses(tmpResponses == 77 | tmpResponses == 79) = 76; % double checked with existing data
                allResponses(tmpResponses == 76 | tmpResponses == 78) = 65;
            elseif mod(id2,2) == 1
                % L is no change, A is change.
                allResponses(tmpResponses == 77 | tmpResponses == 79) = 65;
                allResponses(tmpResponses == 76 | tmpResponses == 78) = 76;
            end

            if strcmp(subjectID, '10') & f == 1
                allResponses = [allResponses(1:30); NaN; allResponses(31:end)];
                allCorrect = [allCorrect(1:30); NaN; allCorrect(31:end)];
                trialSetSize = [trialSetSize(1:30), NaN, trialSetSize(31:end)];
                trialIfChange = [trialIfChange(1:30), NaN, trialIfChange(31:end)];
                trialCuedSide = [trialCuedSide(1:30), NaN, trialCuedSide(31:end)];
            end
            
            % replace
            beh.data.allResponses = allResponses;
            beh.data.allCorrect = allCorrect;
            beh.data.trialSetSize = trialSetSize;
            beh.data.trialIfChange = trialIfChange;
            beh.data.trialCuedSide = trialCuedSide;
        end

        % beginning of the data is long and contains a lot of noise. must be removed. 
        % consider that sometimes 11-15 and 91-95 events are missing.
        % find any possible event. 
        all_events = {'9', '10', '11', '12', '13', '14', '15', ...
                      '89', '90', '91', '92', '93', '94', '95', ...
                      '3', '7', '50', ...
                      '21', '41', '61', ...
                      '22', '42', '62', ...
                      '76','77','78','79', ...
                      };
        types = {EEG.event.type};
        i_all_events = find(ismember(types, all_events));
        latencies = cell2mat({EEG.event.latency});
        % find the seconds
        tmin = latencies(i_all_events(1)) - 100; % keep 200 ms before this stimulus srate = 500
        tmax = latencies(i_all_events(end)) + 100; % keep 200 ms after this stimulus
        % if you try to keep more than 200 ms, the huge artifacts remain in the data
        EEG = pop_select(EEG, 'point', [tmin tmax]);

        % maxFreqInterest = 80;
        % [eegpx,freq2] = pwelch(EEG.data(30,:),2000,1000,[],EEG.srate);
        % eegpx = 10*log10(eegpx);
        % [minval2 ind2]=min(abs(freq2-maxFreqInterest));
        % figure;plot(freq2(1:ind2),eegpx(1:ind2));

        % merge the blocks
        EEGconcat{f} = EEG;
        BEHconcat{f} = beh.data;  

    end

    % merge the blocks
    % if there are empty cells, remove them
    i_empty = ~cellfun('isempty',EEGconcat);
    EEGconcat = EEGconcat(i_empty);
    BEHconcat = BEHconcat(i_empty);

    E = EEGconcat{1};
    for ss = 2 : length(EEGconcat)
        E = pop_mergeset(E, EEGconcat{ss},  0);
    end
    % overwrite EEG
    EEG = E;

    % info for BIDS format
    EEG.tInfo.OldID = subjectID;
    EEG.tInfo.EEGReference = 'E129';
    EEG.tInfo.InstitutionName = 'University of South Florida';
    EEG.tInfo.PowerLineFrequency = 60;
    EEG.tInfo.EEGGround = 'PCz';
    EEG.tInfo.CapManufacturer = 'Magstim EGI';
    EEG.tInfo.SoftwareFilters = "n/a";

    % now concat BD data
    BEH = [BEHconcat{:}];

    newBEH = struct;
    field_names = fieldnames(BEH);
    for n = 1 : length(field_names)
        tmp = [];
        field_size = size(BEH(1).(field_names{n}));
        for ss = 1 : length(EEGconcat)
            if field_size(1) > field_size(2)
                tmp = [tmp; BEH(ss).(field_names{n})]; % concat verticaly
            elseif field_size(1) < field_size(2)
                tmp = [tmp, BEH(ss).(field_names{n})]; % concat horizonztally
            end
        end
        % overwrite the field with merged data
        newBEH.(field_names{n}) = tmp;
    end

    % append to merged EEG dataset
    EEG.beh = beh; % beh contains also experiment info
    EEG.beh.data = newBEH; % replace with merged data
    EEG.beh.triggers = trigger;
  
    % save CDA file
    task = [subjectID, '_CDA_EEG.mat'];
    save(fullfile(result_folder, subjectID, task), 'EEG', '-v7.3')

    % load Resting
    d_sub_resting = dir(strcat(fullfile(d(sub).folder, d(sub).name), '/*rest*mff'));

    % continue if empty
    if isempty(d_sub_resting)
        continue;
    end
    
    % load EEG
    mff_path = fullfile(d_sub_resting(1).folder, d_sub_resting(1).name);
    EEG = mff_import(mff_path);

    % info for BIDS format
    EEG.tInfo.OldID = subjectID;
    EEG.tInfo.EEGReference = 'E129';
    EEG.tInfo.InstitutionName = 'University of South Florida';
    EEG.tInfo.PowerLineFrequency = 60;
    EEG.tInfo.EEGGround = 'PCz';
    EEG.tInfo.CapManufacturer = 'Magstim EGI';
    EEG.tInfo.SoftwareFilters = "n/a";

    % trim whitespaces
    for i = 1:length(EEG.event)
        EEG.event(i).type = strtrim(EEG.event(i).type);
    end
    
    % remove noisy boundaries
    if not(isempty(EEG.event))
        types = {EEG.event.type};
        i_all_events = find(ismember(types, all_events));
        latencies = cell2mat({EEG.event.latency});
        % find the seconds
        tmin = latencies(i_all_events(1)) - 100; % keep 200 ms before this stimulus srate = 500
        tmax = latencies(i_all_events(end)) + 100; % keep 200 ms after this stimulus
        % if you try to keep more than 200 ms, the huge artifacts remain in the data
        EEG = pop_select(EEG, 'point', [tmin tmax]);
    end

    % save
    task = [subjectID, '_Resting_EEG.mat'];
    save(fullfile(result_folder, subjectID, task), 'EEG', '-v7.3')

    % load Eye movement task
    d_sub_eye = dir(strcat(fullfile(d(sub).folder, d(sub).name), '/*eye*mff'));
    
    % load EEG
    mff_path = fullfile(d_sub_eye(1).folder, d_sub_eye(1).name);
    EEG = mff_import(mff_path);

    % info for BIDS format
    EEG.tInfo.OldID = subjectID;
    EEG.tInfo.EEGReference = 'E129';
    EEG.tInfo.InstitutionName = 'University of South Florida';
    EEG.tInfo.PowerLineFrequency = 60;
    EEG.tInfo.EEGGround = 'PCz';
    EEG.tInfo.CapManufacturer = 'Magstim EGI';
    EEG.tInfo.SoftwareFilters = "n/a";

    % trim whitespaces
    for i = 1:length(EEG.event)
        EEG.event(i).type = strtrim(EEG.event(i).type);
    end

    % remove noisy boundaries
    types = {EEG.event.type};
    i_all_events = find(ismember(types, all_events));
    latencies = cell2mat({EEG.event.latency});
    % find the seconds
    tmin = latencies(i_all_events(1)) - 100; % keep 200 ms before this stimulus srate = 500
    tmax = latencies(i_all_events(end)) + 100; % keep 200 ms after this stimulus
    % if you try to keep more than 200 ms, the huge artifacts remain in the data
    EEG = pop_select(EEG, 'point', [tmin tmax]);

    task = [subjectID, '_Eye_EEG.mat'];
    save(fullfile(result_folder, subjectID, task), 'EEG', '-v7.3')

end
