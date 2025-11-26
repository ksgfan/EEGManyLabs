

%% description

% no eye tracker
% TP9 and TP10 were glued on mastoids as M1, M2. It means that the actual
% TP9 and TP10 are missing. See code below.

% 'P3', 'P7', 'O1', 'P4', 'P8', 'O2' available

% subject 4, 5, 19, 24 - more responses than EEG data in CDA

% subject 4: started recording a few seconds too late in block 4
% - to not lose the whole block 4, I will reconstruct beh data from EEG triggers

% subject 5:  recording delayed in block 5, about 100 trials of EEG missing

% subject 19:  11 trials of EEG missing

% subject 24:  1 block of EEG missing

% subject 7: block 3 started twice, abbruch in block 4 - exclude

%% load data
data_path = fullfile(rawdatapath, 'data_zurich_sauseng');
d = dir([data_path, filesep, '*',filesep, '*eeg']);
result_folder = fullfile(formatted_data, 'zurich_sauseng');
mkdir(result_folder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % merge subject 18 and 80 - has to be done only once
% EEG1 = pop_fileio('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/18.eeg'));
% EEG2 = pop_fileio('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/80/80.eeg'));
% EEG = pop_mergeset(EEG1, EEG2);
% 
% mkdir(fullfile(/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/, 'z_old'))
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/18.eeg', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/18.vhdr', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/18.vmrk', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/80.eeg', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/80.vhdr', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/80.vmrk', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/z_old')
% 
% pop_writebva(EEG, /Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/18/18'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % merge subject 20a and 20b - has to be done only once
% EEG1 = pop_fileio('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/20_a.eeg');
% EEG2 = pop_fileio('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/20_b.eeg');
% EEG = pop_mergeset(EEG1, EEG2);
% 
% mkdir('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/20_a.eeg', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/20_a.vhdr', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/20_a.vmrk', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/20_b.eeg', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/20_b.vhdr', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/z_old')
% movefile('/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/20_b.vmrk', '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/z_old')
% 
% pop_writebva(EEG, '/Volumes/methlab_data/EEGManyLabs/raw_data/data_zurich_sauseng/20/20');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dir again, search for eeg and dat files
d = [dir(fullfile(data_path, '*', '*.eeg')); dir(fullfile(data_path, '*', '*.dat'))];

% load chanloc file
locspath = 'CACS-64_REF.bvef';
chanlocs = loadbvef(locspath);
chanlocs(ismember({chanlocs.labels}, 'REF')).labels = 'FCz';
chanlocs(ismember({chanlocs.labels}, 'GND')).labels = 'Fpz';


for f = 1 : size(d, 1)
   
    id = strsplit(d(f).folder, filesep);
    id = id{end};
    if startsWith(id, '0')
        id = id(2);
    end
    
    disp(id)
    
    if strcmp(id, '7') % exclude subject 7
        continue
    end
    
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
    EEG.tInfo.InstitutionName = 'University of Zurich (Prof. Sauseng)';
    EEG.tInfo.PowerLineFrequency = 50;
    EEG.tInfo.EEGGround = 'Fpz';
    EEG.tInfo.CapManufacturer = 'BrainProducts';
    EEG.tInfo.SoftwareFilters = "n/a";
    
    % Rename TP9/TP10 to M1/M2
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'TP9'))).labels = 'M1'; 
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'TP10'))).labels = 'M2';  

    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M1'))).theta = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M1'))).sph_radius = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M1'))).sph_phi = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M1'))).sph_theta = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M1'))).Z = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M1'))).Y = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M1'))).X = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M1'))).radius = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M2'))).theta = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M2'))).sph_radius = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M2'))).sph_phi = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M2'))).sph_theta = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M2'))).Z = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M2'))).Y = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M2'))).X = [];
    EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'M2'))).radius = [];

    for chan = 64 : 65
        EEG.chanlocs(chan).type = 'eog';
    end
    
    % Remove 'S ' from EEG.event.type
    A = {EEG.event.type};
    for i = 1:numel(A)
        if startsWith(EEG.event(i).type, 'S')
            EEG.event(i).type = strip(num2str(A{i}(2:end)));
        end
    end
    
    % make a copy
    EEGorig = EEG;
    
    % cut CDA 
    TASK =  '_CDA_EEG.mat';
    iSTART = find(ismember({EEGorig.event.type}, '11'));
    iEND = find(ismember({EEGorig.event.type}, '95'));
    
    if strcmp(id, '24') % subject 24 completed 4 blocks
        iEND = find(ismember({EEGorig.event.type}, '94'));
    end
    
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(iSTART(end)).latency, EEGorig.event(iEND(end)).latency]);
    
    % concat and append beh data to the CDA file
    if strcmp(TASK, '_CDA_EEG.mat')
    
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

        % append  behavioral
        EEG.beh = beh; % beh contains also experiment info
        EEG.beh.data = newBEH; % replace with merged data
        EEG.beh.triggers = trigger;
    end
    
    
    %%%%%% some EEG data is missing. match beh data with EEG triggers
        
    % double checked with existing data!!

    % RESP no change & correct = 76;
    % RESP change & correct = 77;
    % RESP no change & incorrect = 78;
    % RESP change & incorrect = 79;

    % mod(id,2) == 0 => % L is change, A is no change.
    % mod(id,2) == 1 => % L is no change, A is change.

    % KeyCodeA = 65;
    % KeyCodeL = 76; 

    if strcmp(id, '4') | strcmp(id, '5') | strcmp(id, '19') | strcmp(id, '24')

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
        
        % append new behavioral
        EEG.beh.data = newBEH; % replace with merged data
    end
    %%%%%% end
    
    % save
    mkdir(fullfile(result_folder, id))
    save(fullfile(result_folder, id, [id TASK]), 'EEG', '-v7.3')
    
    % cut eyecalib task
    TASK =  '_Eye_EEG.mat';
    iSTART = find(ismember({EEGorig.event.type}, '9'));
    iEND = find(ismember({EEGorig.event.type}, '89'));
    
    % deal with missing triggers 
    if strcmp(id, '20') % subject 20 does not have triggers 9 and 89
        [EEGorig.event(6).type] = '9';
        [EEGorig.event(127).type] = '89';
        iSTART = find(ismember({EEGorig.event.type}, '9'));
        iEND = find(ismember({EEGorig.event.type}, '89'));
    end
    
    % subjects 1, 15, 16, 19 dont not have trigger 9 
    if strcmp(id, '1') | strcmp(id, '15') | ...
            strcmp(id, '16') | strcmp(id, '19') | ...
            strcmp(id, '22') | strcmp(id, '23') | ...
            strcmp(id, '24') | strcmp(id, '25') | ...
            strcmp(id, '6')
        % replace 'empty' with '9'
        [EEGorig.event(1).type] = '9';
        [EEGorig.event(1).timestamp] = [];
        iSTART = find(ismember({EEGorig.event.type}, '9'));
    end

    if strcmp(id, '8')
        % replace 'empty' with '9'
        [EEGorig.event(2).type] = '9';
        [EEGorig.event(2).timestamp] = [];
        iSTART = find(ismember({EEGorig.event.type}, '9'));
    end
    
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(iSTART(end)).latency, EEGorig.event(iEND(end)).latency]);
    save(fullfile(result_folder, id, [id TASK]), 'EEG', '-v7.3')
    
    % cut resting
    TASK =  '_Resting_EEG.mat';
    iSTART = find(ismember({EEGorig.event.type}, '10'));
    
    if strcmp(id, '17') % subject 17 does not have trigger 10
        iSTART = find(ismember({EEGorig.event.type}, '89'));
    end
    
    iEND = find(ismember({EEGorig.event.type}, '90'));
    
    if strcmp(id, '18') % subject 18 does not have trigger 90
        iEND = find(ismember({EEGorig.event.type}, 'boundary'));
    end
        
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(iSTART(end)).latency, EEGorig.event(iEND(end)).latency]);
    save(fullfile(result_folder, id, [id TASK]), 'EEG', '-v7.3')
 
end


