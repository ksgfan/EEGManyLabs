

%% description
% BrainVision 32 chans
% 32 channels (31 EEG + diode)
% reference: Cz
% no resting
% no eye tracker
% As for the EOG - we used Fp1 Fp2 (VEOG) and FT9 FT10 (HEOG)
% triggers for CDA are start: 10 and end: 90 for each block

% 'P3', 'P7', 'O1', 'P4', 'P8', 'O2' available

% EEG data missing: 2, 30, 44
% beh data missing: 1, 23, 4, 3  => recreate beh data from triggers 
% subject 31 has more events (720) than responses in behavioral data.
% recreate the beh file from EEG.events

% subject 31: raw files .vhdr and .vmrk were manually opened and the
% content was changed everywhere from subjectID '30' to '31'. If you dont do
% it, the pop_loadbv will crash.

%% load data
data_path = fullfile(rawdatapath, 'data_reykjavik');
d = dir([data_path, filesep, '*',filesep, '*vhdr']);
result_folder = fullfile(formatted_data, 'reykjavik');
mkdir(result_folder)

% load chanloc file
locspath = 'BrainVision-AC-32.bvef';
chanlocs = loadbvef(locspath);
chanlocs(ismember({chanlocs.labels}, 'GND')).labels = 'Fpz';

for f = 23 : size(d, 1)
   
    id = strsplit(d(f).folder, filesep);
    id = id{end};
    if startsWith(id, '0')
        id = id(2);
    end
    
    disp(id)

    EEG = pop_loadbv(d(f).folder, d(f).name);
       
    % remove photiodiode
    el_diod = find(strcmp({EEG.chanlocs.labels}, 'diod')); 
    EEG = pop_select(EEG, 'nochannel', el_diod);
    
    % add reference
    EEG.data(end+1, :) = 0;
    EEG.nbchan = EEG.nbchan+1;
    EEG.chanlocs(EEG.nbchan).labels = 'Fz';
    for ch = 1 : EEG.nbchan
        EEG.chanlocs(ch).type = 'eeg';  
    end
    EEG.ref = 'Fz';
    
    % add locs for 'Fz'
    locFields = fieldnames(chanlocs);
    tf = ismember({chanlocs.labels}, 'Fz');
    tf_new = ismember({EEG.chanlocs.labels}, 'Fz');

    % corresponding entry in the new file
    src = chanlocs(tf);     
    for fld = 1:numel(locFields)
        fn = locFields{fld};
        if isfield(src, fn)
            EEG.chanlocs(tf_new).(fn) = src.(fn);
        end
    end
    EEG = eeg_checkset(EEG);

    % info for BIDS format
    EEG.tInfo.OldID = id;
    EEG.tInfo.EEGReference = 'Fz';
    EEG.tInfo.InstitutionName = 'University of Iceland';
    EEG.tInfo.PowerLineFrequency = 50;
    EEG.tInfo.EEGGround = 'Fpz';
    EEG.tInfo.CapManufacturer = 'BrainVision';
    EEG.tInfo.SoftwareFilters = "n/a";
    
    % % sanity check
    % figure;
    % topoplot([], EEG.chanlocs, 'electrodes', 'labels')
    
    % Remove 'S ' from EEG.event.type
    A = {EEG.event.type};
    for i = 1:numel(A)
        if startsWith(EEG.event(i).type, 'S')
            EEG.event(i).type = strip(num2str(A{i}(2:end)));
        end
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
    % Isolate eye calibration task - first trigger to the last 89
    i9 = 1;
    i89 = find(ismember({EEGorig.event.type}, '89'));
    if not(isempty(i9)) & not(isempty(i89))
        i89 = i89(end);
        EEG = pop_select(EEGorig, 'point', [EEGorig.event(i9).latency EEGorig.event(i89).latency]);
        save(fullfile(result_folder, id, [id '_Eye_EEG.mat']), 'EEG')
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
    
    
    % CDA - first trigger 10 to last trigger 90
    i10 = find(ismember({EEGorig.event.type}, '10'));
    i90 = find(ismember({EEGorig.event.type}, '90'));
    if not(isempty(i10)) & not(isempty(i90))
        i10 = i10(1);
        i90 = i90(end);
        EEG = pop_select(EEGorig, 'point', [EEGorig.event(i10).latency EEGorig.event(i90).latency]); 
        
        % subject 3, 12 have more trials than responses.
        if strcmp(id, '3') 
            i_events = find(ismember({EEG.event.type}, {'21', '41', '61'}));
            length(i_events)

            i_10 = find(ismember({EEG.event.type}, '10')); % subject 3 has 6 x 10. the fivth should be removed 

            % remove the first part
            EEG = pop_select(EEG, 'nopoint', [EEG.event(i_10(5)).latency, EEG.event(i_10(6)-1).latency]);

            % sanity check
            i_events = find(ismember({EEG.event.type}, {'21', '41', '61'}));
            length(i_events)
        end
        
        % now subject 12
        if strcmp(id, '12') 
            i_events = find(ismember({EEG.event.type}, {'21', '41', '61'}));
            length(i_events)

            i_10 = find(ismember({EEG.event.type}, '10')); % subject 12 has 6 x 10. the first should be removed 

            % remove the first part
            EEG = pop_select(EEG, 'nopoint', [EEG.event(i_10(1)).latency, EEG.event(i_10(2)-1).latency]);

            % sanity check
            i_events = find(ismember({EEG.event.type}, {'21', '41', '61'}));
            length(i_events)
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

        if strcmp(id, '1') | strcmp(id, '3') | strcmp(id, '4') | strcmp(id, '23') | strcmp(id, '31')

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
        
        % append  behavioral
        EEG.beh = beh; % beh contains also experiment info
        EEG.beh.data = newBEH; % replace with merged data
        EEG.beh.triggers = trigger;
        save(fullfile(result_folder, id, [id '_CDA_EEG.mat']), 'EEG', '-v7.3')
    else
        disp(['Subject ' id ': CDA task is missing...'])
    end   
end

