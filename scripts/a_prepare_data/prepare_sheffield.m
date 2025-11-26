

%% description

% Reference; grounding: unreferenced/reference free; 
% CMS/DRL adjacent to POz. 

% there are always 2 files: 
% - original with 2048 Hz sampling rate 
% - downsampled to 512 Hz (deci)
% I will be using the 512 Hz from now on

% HEOG_chan = {'EXG1', 'EXG2'};  HEOGL HEOGR; Numbers: 65 66 
% VEOG_chan = {'EXG3', 'EXG4'}; VEOGU VEOGL; U is upper L is lower;  Numbers: 67 68
% Mastoids = {'EXG5', 'EXG6'}; L-M1 R-M2; Numbers: 69 70

% no resting state
% no eye tracker

% subject 3 has more EEG trials than responses. Block 1 was started 3 times. => exclude the first 2
% subject 13 has no trigger 9. 
% subject 16: behavioral data from block 1 was overwritten by block 4. => exclude subject
% subject 18: more triggers than responses in beh file. recompute responses

% visual inspection:
% 11 - large number of blinks and saccades
% 14 - large number of blinks and saccades
% 18 - large number of blinks and saccades
% 2 - large number of blinks and saccades
% 20 - eye calibration task rejects 170 trials
% 21 - drifts in the data (around 100-150 trials)
% 26 - large number of blinks and saccades


%% load data
data_path = fullfile(rawdatapath, 'data_sheffield');
d = dir([data_path, filesep, '*',filesep, '*Deci.bdf']);
result_folder = fullfile(formatted_data, 'sheffield');
mkdir(result_folder)

for f = 1 : size(d, 1)
   
    id = strsplit(d(f).folder, filesep);
    id = id{end};
    if startsWith(id, '0')
        id = id(2);
    end
    
    disp(id)
    
    % exclude subject 16
    if strcmp(id, '16')
        continue
    end
    
    EEG = pop_fileio(fullfile(d(f).folder, d(f).name));
    
    % rename channels
    EEG.chanlocs(65).labels = 'HEOGL';
    EEG.chanlocs(66).labels = 'HEOGR';
    
    EEG.chanlocs(67).labels = 'VEOGU'; % above left
    EEG.chanlocs(68).labels = 'VEOGL'; % below left
    
    EEG.chanlocs(69).labels = 'M1'; % left mastoid
    EEG.chanlocs(70).labels = 'M2';
    
    % remove remaining channels
    EEG = pop_select(EEG, 'nochannel', 71 : 80);
    
    % load chanloc file
    locspath = 'Biosemi64.locs';
    EEG = pop_chanedit(EEG, 'lookup', locspath);

    for chan = 65 : 68
        EEG.chanlocs(chan).type = 'eog';
    end
    for chan = 69 : 70
        EEG.chanlocs(chan).type = 'eeg';
    end
       
    % % sanity check
    % figure;
    % topoplot([], EEG.chanlocs, 'electrodes', 'labels')

    % info for BIDS format
    EEG.tInfo.OldID = id;
    EEG.tInfo.EEGReference = 'Reference free';
    EEG.tInfo.InstitutionName = 'University of Sheffield';
    EEG.tInfo.PowerLineFrequency = 50;
    EEG.tInfo.EEGGround = 'adjacent to POz';
    EEG.tInfo.CapManufacturer = 'Biosemi';
    EEG.tInfo.SoftwareFilters = "n/a";
    
    % Convert numbers to strings in EEG.event
    A = {EEG.event.type};
    for i = 1:numel(A)
        EEG.event(i).type = num2str(A{i});
    end

    % check for 'empty' event
    if any(strcmp({EEG.event.type}, 'empty'))
        i_empty = strcmp({EEG.event.type}, 'empty');
        [EEG.event(i_empty).type] = deal('999');
        [EEG.event(i_empty).timestamp] = deal([]);
    end
    if any(strcmp({EEG.urevent.type}, 'empty'))
        i_empty = strcmp({EEG.urevent.type}, 'empty');
        [EEG.urevent(i_empty).type] = deal('999');
        [EEG.urevent(i_empty).timestamp] = deal([]);
    end
    
    % extract CDA and resting
    EEGorig = EEG;
    
    mkdir(fullfile(result_folder, id))
    % Isolate eye calibration task
    i9 = find(ismember({EEGorig.event.type}, '9'));
    % special case of subject 13
    if strcmp(id, '13')
        i9 = 2;
    end
    i89 = find(ismember({EEGorig.event.type}, '89'));
    if not(isempty(i9)) & not(isempty(i89))
        EEG = pop_select(EEGorig, 'point', [EEGorig.event(i9).latency EEGorig.event(i89).latency]);
        save(fullfile(result_folder, id, [id '_Eye_EEG.mat']), 'EEG', '-v7.3')
    else
        disp(['Subject ' id ': Eye calibration task is missing...'])
    end
    
    % load and merge CDA behavioral data 
    BEHconcat = {};
    beh = [];
    trigger = [];
    d_beh = dir(fullfile(d(f).folder, '*ManyLabs*block*task.mat'));
    for b = 1 : size(d_beh, 1)
        load(fullfile(d_beh(b).folder, d_beh(b).name)); % behavioral
        BEHconcat{b} = beh.data;  
    end
    % now concat BD data
    BEH = [BEHconcat{:}];
    newBEH = struct;
    if not(isempty(BEH))
        field_names = fieldnames(BEH);
        for n = 1 : length(field_names)
            tmp = [];
            field_size = size(BEH(1).(field_names{n}));
            for ss = 1 : length(BEHconcat)
                if field_size(1) > field_size(2)
                    tmp = [tmp; BEH(ss).(field_names{n})]; % concat verticaly
                elseif field_size(1) < field_size(2)
                    tmp = [tmp, BEH(ss).(field_names{n})]; % concat horizonztally
                end
            end
            % overwrite the field with merged data
            newBEH.(field_names{n}) = tmp;
        end
    end

    %%%%%% some beh data is missing. recompute from EEG.events
        
    % double checked with existing data!!
    
    % RESP no change & correct = 76;
    % RESP change & correct = 77;
    % RESP no change & incorrect = 78;
    % RESP change & incorrect = 79;

    % mod(id,2) == 0 => % L is change, A is no change.
    % mod(id,2) == 1 => % L is no change, A is change.

    % KeyCodeA = 65;
    % KeyCodeL = 76; 

    if strcmp(id, '18')

        % first, we need longer segment than contains cues and responses
        tmpEEG = pop_epoch(EEGorig, {21, 41, 61}, [-0.7 3]);
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

        id2 = str2double(id);
        if mod(id2,2) == 0 
            % L is change, A is no change.
            allResponses(tmpResponses == 77 | tmpResponses == 79) = 76; % double checked with existing data
            allResponses(tmpResponses == 76 | tmpResponses == 78) = 65;
        elseif mod(id2,2) == 1
            % L is no change, A is change.
            allResponses(tmpResponses == 77 | tmpResponses == 79) = 65;
            allResponses(tmpResponses == 76 | tmpResponses == 78) = 76;
        end
        
        newBEH.allResponses = allResponses;
        newBEH.allCorrect = allCorrect;
        newBEH.trialSetSize = trialSetSize;
        newBEH.trialIfChange = trialIfChange;
        newBEH.trialCuedSide = trialCuedSide;
    end
    %%%%%% end

 
    % CDA
    i11 = find(ismember({EEGorig.event.type}, '11'));
    % special case subject 3
    if strcmp(id, '3')
        i11 = i11(end);
    end
    i95 = find(ismember({EEGorig.event.type}, '95'));
    if not(isempty(i11)) & not(isempty(i95))
        EEG = pop_select(EEGorig, 'point', [EEGorig.event(i11).latency EEGorig.event(i95).latency]); 
        % append  behavioral
        EEG.beh = beh; % beh contains also experiment info
        EEG.beh.data = newBEH; % replace with merged data
        EEG.beh.triggers = trigger;
        save(fullfile(result_folder, id, [id '_CDA_EEG.mat']), 'EEG', '-v7.3')
    else
        disp(['Subject ' id ': CDA task is missing...'])
    end
end


