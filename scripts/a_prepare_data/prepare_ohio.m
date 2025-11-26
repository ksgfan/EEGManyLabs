

%% description
% BrainVision 32 chans
% no resting
% Fp1 Fp2 (VEOG) and FT9 FT10 (HEOG)
% reference 'Cz'

% subject "sub02_cda-001.eeg" has to be manually renamed to "sub02_cda.eeg"

% subject 15 has more trials than responses. double check with EEG.event.
% later, double check with ET data, when provided

%% convert ET files (asci to mat)

%%%%%%% IMPORTANT !!!!!!!
% convert asci to mat - use the parseeyelink that is changed by Dawid (line282) should be: 
% test  = regexp(et.messages,'MSG\s+(\d+)\s+(.*)','tokens')'; => MSG not INPUT

% d = dir(strcat(fullfile(rawdatapath, 'data_ohio/*/', '*edf')));
% 
% for i = 1 : size(d, 1) 
% 
%     id = strsplit(d(i).folder, filesep);
%     id = id{end};
%     if startsWith(id, '0')
%         id = id(2);
%     end
% 
%     disp(id)
% 
%     % convert to text files
%     disp("CONVERTING EDF to ASCII...")
%     system([pathEdf2Asc ' "' fullfile(d(i).folder, d(i).name) '" -y']);
% 
%     % convert to mat files
%     filename = d(i).name(1:end-4); % remove .edf   
%     ascFile = strcat(filename, '.asc');
%     inputFile = fullfile(d(i).folder, ascFile);
%     x = strsplit(filename, '_');
% 
%     if strcmp(x{2}, '1') | strcmp(x{2}, '2') | strcmp(x{2}, '3') | strcmp(x{2}, '4') | strcmp(x{2}, '5')
%         newfilename = [x{1}, '_ManyLabsCDA_block' x{2} '_task_ET.mat'];
%     elseif strcmp(x{2}, 'Eye')
%         newfilename = [x{1}, '_Eye_ET.mat'];
%     elseif strcmp(x{2}, 'Res')
%         newfilename = [x{1}, '_Resting_ET.mat'];
%     end
% 
%     % convert
%     outputFile = fullfile(d(i).folder, newfilename);
%     ET = parseeyelink(inputFile, outputFile);
% 
% end

%% merge EEG and ET data, and save as mat file
d_eeg = dir(strcat(fullfile(rawdatapath, 'data_ohio/*/', '*cda*eeg')));
d_eye = dir(strcat(fullfile(rawdatapath, 'data_ohio/*/', '*cal*eeg')));

result_folder = fullfile(formatted_data, 'ohio');

% load chanloc file
locspath = 'BrainVision-AC-32.bvef';
chanlocs = loadbvef(locspath);
chanlocs(ismember({chanlocs.labels}, 'GND')).labels = 'Fpz';

for f = 1 : size(d_eeg, 1)

    id = strsplit(d_eeg(f).folder, filesep);
    id = id{end};
    if startsWith(id, '0')
        id = id(2);
    end
    
    disp(id)
    
    % load EEG data
    EEG = pop_fileio(fullfile(d_eeg(f).folder, d_eeg(f).name));
    
    % load chanloc file
    % first, add reference channel Cz
    EEG.data(32, :) = 0;
    EEG.nbchan = 32;
    EEG.chanlocs(32).labels = 'Cz';   
    EEG.chanlocs(32).type = 'eeg'; 
    EEG.ref = 'Cz';

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


    % info for BIDS format
    EEG.tInfo.OldID = id;
    EEG.tInfo.EEGReference = 'Cz';
    EEG.tInfo.InstitutionName = 'The Ohio State University';
    EEG.tInfo.PowerLineFrequency = 60;
    EEG.tInfo.EEGGround = 'Fpz';
    EEG.tInfo.CapManufacturer = 'BrainVision';
    EEG.tInfo.SoftwareFilters = "n/a";
    
    % % rename channels
    % EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'Fp1'))).labels = 'VEOGU'; % above left
    % EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'FT9'))).labels = 'HEOGL';
    % EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'FT10'))).labels = 'HEOGR';

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

    % file of subject 2 is huge. recording was probably not stopped after
    % the measurement. cut the data.
    if strcmp(id, '2')
        EEG = pop_select(EEG, 'point', [EEG.event(1).latency, EEG.event(3607).latency]);
    end
    
    % merge with ET. ET comes in blocks
    STARTtriggers = 11:15;
    ENDtriggers = 91:95;
    for et_file = 1 : 5
        etfile = fullfile(d_eeg(f).folder, strcat(id, '_ManyLabsCDA_block', num2str(et_file), '_task_ET.mat'));
        if exist(etfile, 'file') == 2
            EEG = pop_importeyetracker(EEG, etfile, [STARTtriggers(et_file) ENDtriggers(et_file)],[2 3],{'L_GAZE_X', 'L_GAZE_Y'},1,1,1,0, 10);
        end
    end
    
    % merging multiple times adds duplicate EEG chanlocs and data. we can
    % simply sumcolumn the X and Y gaze to obtain 1 time series
    i_x = ismember({EEG.chanlocs.labels}, 'L-GAZE-X');
    i_y = ismember({EEG.chanlocs.labels}, 'L-GAZE-Y');
    if not(isempty(find(i_x)))
        ETdataX = sum(EEG.data(i_x, :), 1);
        ETdataY = sum(EEG.data(i_y, :), 1);

        % remove reduntant channels and paste ET data 
        EEG.chanlocs(35:EEG.nbchan) = [];
        EEG.data(35:EEG.nbchan, :) = [];
        EEG.data(33:34, :) = [ETdataX; ETdataY];
        EEG.nbchan = 34;
    end
    
    % load and merge CDA behavioral data 
    BEHconcat = {};
    beh = [];
    trigger = [];
    d_beh = dir(fullfile(d_eeg(f).folder, '*ManyLabs*block*task.mat'));
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
    

    % subject 15 has more trials than responses.
    if strcmp(id, '15')
        i_events = find(ismember({EEG.event.type}, {'21', '41', '61'}));
        length(i_events)
        
        i_15 = find(ismember({EEG.event.type}, '15')); % 2 x 15. block5 was started twice
        
        % remove the first part
        EEG = pop_select(EEG, 'nopoint', [EEG.event(i_15(1)).latency, EEG.event(i_15(2)-1).latency]);
        
        % sanity check
        i_events = find(ismember({EEG.event.type}, {'21', '41', '61'}));
        length(i_events)
    end
        
    % append  behavioral and save
    EEG.beh = beh; % beh contains also experiment info
    EEG.beh.data = newBEH; % replace with merged data
    EEG.beh.triggers = trigger;
    mkdir(fullfile(result_folder, id))
    save(fullfile(result_folder, id, [id '_CDA_EEG.mat']), 'EEG', '-v7.3')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load eye calibration data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    EEG = pop_fileio(fullfile(d_eye(f).folder, d_eye(f).name));
    
    % load chanloc file
    % first, add reference channel Cz
    EEG.data(32, :) = 0;
    EEG.nbchan = 32;
    EEG.chanlocs(32).labels = 'Cz';   
    locspath = 'standard_1005.elc';
    EEG = pop_chanedit(EEG, 'lookup', locspath);

    % info for BIDS format
    EEG.tInfo.OldID = id;
    EEG.tInfo.EEGReference = 'Cz';
    EEG.tInfo.InstitutionName = 'The Ohio State University';
    EEG.tInfo.PowerLineFrequency = 60;
    EEG.tInfo.EEGGround = 'Fpz';
    EEG.tInfo.CapManufacturer = 'BrainVision';
    EEG.tInfo.SoftwareFilters = struct();
    
    % % rename channels
    % EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'Fp1'))).labels = 'VEOGU'; % above left
    % EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'FT9'))).labels = 'HEOGL';
    % EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'FT10'))).labels = 'HEOGR';
    
    % Remove 'S ' from EEG.event.type
    A = {EEG.event.type};
    for i = 1:numel(A)
        if startsWith(EEG.event(i).type, 'S')
            EEG.event(i).type = strip(num2str(A{i}(2:end)));
        end
    end

    % check for 'empty' event
    if any(strcmp({EEG.event.type}, 'empty'))
        i_empty = strcmp(A, 'empty');
        EEG.event(i_empty).type = '999';
        EEG.event(i_empty).timestamp = [];
    end
    if any(strcmp({EEG.urevent.type}, 'empty'))
        i_empty = strcmp(A, 'empty');
        EEG.urevent(i_empty).type = '999';
        EEG.urevent(i_empty).timestamp = [];
    end
    
    % merge with ET
    etfile = fullfile(d_eye(f).folder, strcat(id, '_Eye_ET.mat'));
    
    % merge
    if exist(etfile, 'file') == 2
        % subject 2: start trigger '9' is missing in EEG. merge using
        % different trigger
        if strcmp(id, '2')
            EEG = pop_importeyetracker(EEG, etfile, [104 89],[2 3],{'L_GAZE_X', 'L_GAZE_Y'},1,1,1,0);
        else
            EEG = pop_importeyetracker(EEG, etfile, [9 89],[2 3],{'L_GAZE_X', 'L_GAZE_Y'},1,1,1,0);
        end
    end
    
    % save
    save(fullfile(result_folder, id, [id '_Eye_EEG.mat']), 'EEG', '-v7.3')
    
end



