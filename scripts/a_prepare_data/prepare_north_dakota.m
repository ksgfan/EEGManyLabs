

%% description

% 65 =  VEOG (above left eye)
% 66 =  VEOG (above right eye)
% 67 = Left HEOG
% 68 = Right HEOG
% 69 = VEOG (below right eye)
% 70 = VEOG (below left eye)
% 71 = left mastoid
% 72 = right mastoid
%  
% The data I sent you is referenced just to Ch71 (left mastoid). 
% We typically re-reference to algebraic average of left and right mastoids during pre-processing.

% 'P3', 'P7', 'O1', 'P4', 'P8', 'O2' available

% ET has sampling rate of 250 Hz, EEG 500 Hz

data_path = fullfile(rawdatapath, 'data_north_dakota');
result_folder = fullfile(formatted_data, 'north_dakota');
mkdir(result_folder)
%% load data

% IMPORTANT !!!!!
% convert asci to mat - use the parseeyelink that is changed by Dawid (line282) should be: 
% test  = regexp(et.messages,'MSG\s+(\d+)\s+(.*)','tokens')'; => MSG not INPUT

d = dir([data_path, filesep, '*',filesep, '*bdf']);

for f = 1 : size(d, 1)
   
    id = strsplit(d(f).folder, filesep);
    id = id{end};
    if startsWith(id, '0')
        id = id(2);
    end
    
    disp(id)
    
    % convert to text files
    disp("CONVERTING EDF to ASCII...")
    edfFile = strcat(id, '.edf');
    system([pathEdf2Asc ' "' fullfile(d(f).folder, edfFile) '" -y']);
    ascFile = strcat(id, '.asc');
    inputFile = fullfile(d(f).folder, ascFile);
    % convert to mat file and save
    outputFile = fullfile(d(f).folder, strcat(id, '_ET.mat'));
    ET = parseeyelink(inputFile, outputFile);
    % no need for cutting of ET files, but check if something is missing
    i9 = find(ismember(ET.event(:, 2), 9));
    i89 = find(ismember(ET.event(:, 2), 89));
    if isempty(i9) | isempty(i89)
        disp(['Subject ' id ': Eye calibration task - ET - is missing...'])
    end
    i10 = find(ismember(ET.event(:, 2), 10));
    i90 = find(ismember(ET.event(:, 2), 90));
    if isempty(i10) | isempty(i90)
        disp(['Subject ' id ': Resting state - ET - is missing...'])
    end
    i11 = find(ismember(ET.event(:, 2), 11));
    i95 = find(ismember(ET.event(:, 2), 95));
    if isempty(i11) | isempty(i95)
        disp(['Subject ' id ': CDA task - ET - is missing...'])
    end
    
    % load EEG data
    EEG = pop_fileio(fullfile(d(f).folder, d(f).name));
    
    % remove channel 73
    EEG = pop_select(EEG, 'nochannel', 73);
    
    % rename EX channels
    EEG.chanlocs(65).labels = 'VEOGU'; % above left
    EEG.chanlocs(70).labels = 'VEOGL'; % below left

    EEG.chanlocs(67).labels = 'HEOGL';
    EEG.chanlocs(68).labels = 'HEOGR';

    EEG.chanlocs(71).labels = 'M1'; % left mastoid (old reference)
    EEG.chanlocs(72).labels = 'M2';

    for i = 1 : 72
        EEG.chanlocs(i).ref = 'M1';
    end

    % info for BIDS format
    EEG.tInfo.OldID = id;
    EEG.tInfo.EEGReference = 'Reference free';
    EEG.tInfo.InstitutionName = 'North Dakota State University';
    EEG.tInfo.PowerLineFrequency = 60;
    EEG.tInfo.EEGGround = 'adjacent to POz';
    EEG.tInfo.CapManufacturer = 'Biosemi';
    EEG.tInfo.SoftwareFilters = "n/a";

    % remove remianing EX channels
    i_ex = find(contains({EEG.chanlocs.labels}, 'EX'));
    EEG = pop_select(EEG, 'nochannel', i_ex);
    
    % load chanloc file
    locspath = 'BioSemi_64_ND.elp';
    EEG = pop_chanedit(EEG, 'lookup', locspath);

    for chan = 65 : 68
        EEG.chanlocs(chan).type = 'eog';
    end
    for chan = 69 : 70
        EEG.chanlocs(chan).type = 'eeg';
    end
    EEG.chanlocs(37).type = 'eeg';
    

    % % sanity check
    % figure;
    % topoplot([], EEG.chanlocs, 'electrodes', 'labels')

    % downsample already now to match ET frequency
    EEG = pop_resample(EEG, 250);
    
    % Convert numbers to strings in EEG.event
    A = {EEG.event.type};
    for i = 1:numel(A)
        EEG.event(i).type = num2str(A{i});
    end
    
    % extract CDA and resting
    EEGorig = EEG;
    
    % et file is the same for all tasks
    etfile = fullfile(d(f).folder, strcat(id, '_ET.mat'));
    % sometimes block-ending EEG triggers are missing, while ET triggers are not. load the ET
    % file to fix it.
    ET = load(etfile);
    STARTtriggers = 9:15;
    ENDtriggers = 89:95;
    for trig = 1 : length(ENDtriggers)
        iTrigEEGend = find(ismember({EEGorig.event.type}, num2str(ENDtriggers(trig))));
        iTrigEEGstart = find(ismember({EEGorig.event.type}, num2str(STARTtriggers(trig))));
        if isempty(iTrigEEGend) & not(isempty(iTrigEEGstart))
            % EEG.event is missing. recompute using ET.event and add it to EEG.event
            iTrigETend = find(ismember(ET.event(:, 2), ENDtriggers(trig)));
            iTrigETstart = find(ismember(ET.event(:, 2), STARTtriggers(trig)));
            % find the latency difference between start and end trigger
            lat_diff = ET.event(iTrigETend, 1) - ET.event(iTrigETstart, 1);
            % be carefull, because ET is in seconds and EEG in sampling
            % points
            lat_diff = lat_diff / 4;
            
            EEGorig.event(end+1).type = num2str(ENDtriggers(trig));
            EEGorig.event(end).latency = EEGorig.event(iTrigEEGstart).latency + lat_diff;
            EEGorig.event(end).duration = 0;
            EEGorig.event(end).urevent = EEGorig.event(end-1).urevent + 1;
        end
    end
    % sort the triggers
    [~, sortIdx] = sort([EEGorig.event.latency]);
    EEGorig.event = EEGorig.event(sortIdx);
    
    for ev = 1 : length(EEGorig.event)
        EEGorig.event(ev).urevent = ev;
    end
    EEGorig.urevent = EEGorig.event;
    
    % if trigger 95 was added, extend times and data
    i95 = find(ismember({EEGorig.event.type}, '95'));
    if EEGorig.event(i95).latency > length(EEGorig.times)
        tp_diff = ceil(EEGorig.event(i95).latency - length(EEGorig.times));
        EEGorig.times = [EEGorig.times, EEGorig.times(end) + 4 : 4 : EEGorig.times(end) + tp_diff*4];
        EEGorig.data = [EEGorig.data, zeros(70, tp_diff)];
        EEGorig.pnts = size(EEGorig.data, 2);
        EEGorig.xmax = EEG.times(end);
    end

    
    % merge ET and EEG
    % need to do it block by block, as the ET file contains pauses. 
    % if we do not account for the pauses, we get following warning:
    % Warning: For 3696 of 3713 EEG events (99.5%), no ET event of same name was found within plusminus 4 samples
    % This can occur, for example if some events were not transmitted to the ET (or EEG)
    for trig = 1 : length(ENDtriggers)
        try
            EEGorig = pop_importeyetracker(EEGorig, etfile, [STARTtriggers(trig) ENDtriggers(trig)],[2 3],{'L_GAZE_X', 'L_GAZE_Y'},1,1,1,0);
        catch
        end
    end
    
    % merging multiple times adds duplicate EEG chanlocs and data. we can
    % simply sumcolumn the X and Y gaze to obtain 1 time series
    i_x = ismember({EEGorig.chanlocs.labels}, 'L-GAZE-X');
    i_y = ismember({EEGorig.chanlocs.labels}, 'L-GAZE-Y');
    ETdataX = sum(EEGorig.data(i_x, :), 1);
    ETdataY = sum(EEGorig.data(i_y, :), 1);
    
    % remove reduntant channels and paste ET data 
    EEGorig.chanlocs(73:EEGorig.nbchan) = [];
    EEGorig.data(73:EEGorig.nbchan, :) = [];
    EEGorig.data(71:72, :) = [ETdataX; ETdataY];
    EEGorig.nbchan = 72;

    % cut into blocks and save
    mkdir(fullfile(result_folder, id))
    % Isolate eye calibration task
    i9 = find(ismember({EEGorig.event.type}, '9'));
    i89 = find(ismember({EEGorig.event.type}, '89')); % trigger 89 is missing
    if not(isempty(i9)) & not(isempty(i89))
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
 
    % CDA
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



%%

