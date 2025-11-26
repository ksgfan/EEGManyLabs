

%% description

% 63 channels
% no eye tracker
% no M1, M2
% reference: FCz


%% load data
data_path = fullfile(rawdatapath, 'data_mainz/EEG_data');
beh_data_path = fullfile(rawdatapath, 'data_mainz/Behavioral_data');
d = dir([data_path,filesep, '*eeg']);
result_folder = fullfile(formatted_data, 'mainz');
mkdir(result_folder)

% load chanloc file
locspath = 'CACS-64_REF.bvef';
chanlocs = loadbvef(locspath);
chanlocs(ismember({chanlocs.labels}, 'REF')).labels = 'FCz';
chanlocs(ismember({chanlocs.labels}, 'GND')).labels = 'Fpz';

for f = 1 : size(d, 1)
   
    id = strsplit(d(f).name, '.');
    id = id{1};
    id = id(end-1:end);
    if startsWith(id, '0')
        id = id(2);
    end
    
    disp(id)
    
    EEG = pop_fileio(fullfile(d(f).folder, d(f).name));

    % add reference
    EEG.data(EEG.nbchan+1, :) = 0;
    EEG.nbchan = EEG.nbchan+1;
    EEG.chanlocs(EEG.nbchan).labels = 'FCz'; 
    EEG.chanlocs(EEG.nbchan).type = 'eeg'; 
    EEG.ref = 'FCz';

    % add chanlocs
    locFields = fieldnames(chanlocs);
    newChanlocs = EEG.chanlocs;
    labels_old = upper(string({EEG.chanlocs.labels}));
    labels_new = upper(string({chanlocs.labels}));

    % map old order to new locations
    [tf, idxNew] = ismember(labels_old, labels_new);

    for chans = 1:numel(newChanlocs)
        if tf(chans)
            % corresponding entry in the new file
            src = chanlocs(idxNew(chans));     
            for fld = 1:numel(locFields)
                fn = locFields{fld};
                if isfield(src, fn)
                    newChanlocs(chans).(fn) = src.(fn);
                end
            end
        else
            % keep as-is unmatched channels  
        end
    end
        
    % assign back
    EEG.chanlocs = newChanlocs;
    EEG = eeg_checkset(EEG);

    % % sanity check
    % figure;
    % topoplot([], chanlocs, 'electrodes', 'labels')
    % figure;
    % topoplot([], newChanlocs, 'electrodes', 'labels')
    % figure;
    % topoplot([], EEG.chanlocs, 'electrodes', 'labels')
    
    % info for BIDS format
    EEG.tInfo.OldID = id;
    EEG.tInfo.EEGReference = 'FCz';
    EEG.tInfo.InstitutionName = 'University of Mainz';
    EEG.tInfo.PowerLineFrequency = 50;
    EEG.tInfo.EEGGround = 'Fpz';
    EEG.tInfo.CapManufacturer = 'BrainProducts';
    EEG.tInfo.SoftwareFilters = "n/a";
    
     
%     % sanity check
    % figure;
    % topoplot([], chanlocs, 'electrodes', 'labels')

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
        save(fullfile(result_folder, id, [id '_Eye_EEG.mat']), 'EEG', '-v7.3')
    else
        disp(['Subject ' id ': Eye calibration task is missing...'])
    end
    
    % Resting
    i10 = find(ismember({EEGorig.event.type}, '10'));
    i90 = find(ismember({EEGorig.event.type}, '90'));
    if not(isempty(i10)) & not(isempty(i90))
        EEG = pop_select(EEGorig, 'point', [EEGorig.event(i10).latency EEGorig.event(i90).latency]); 
        save(fullfile(result_folder, id, [id '_Resting_EEG.mat']), 'EEG', '-v7.3')
    else
        disp(['Subject ' id ': Resting state is missing...'])
    end
    
    % load and merge CDA behavioral data 
    BEHconcat = {};
    beh = [];
    trigger = [];
    d_beh = dir(fullfile(beh_data_path,id, '*ManyLabs*block*.mat'));
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
    
    
    % CDA - first trigger 11 to last trigger 95
    i11 = find(ismember({EEGorig.event.type}, '11'));
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



