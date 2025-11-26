

%% notes

% reference: Cpz

% visual inspection. Subjects:
% 13 - huge amplitudes at electrodes of interest
% 17 - saccade after almost each stimulus onset
% 29 - eye calibration task rejects 370 trials (ET did not record any saccades)
% 39 - saccade after almost each stimulus onset
% 41 - large number of blinks and saccades
% 42 - large number of blinks and saccades
% 45 - saccade after almost each stimulus onset
% 47 - +100 blinks and saccades after stimulus onset
% 49 - saccade after almost each stimulus onset


% 22 - large amplitudes in eye calib task

%% convert ET files (asci to mat)

% IMPORTANT !!!!!
% convert asci to mat - use the parseeyelink that is changed by Dawid (line282) should be: 
% test  = regexp(et.messages,'MSG\s+(\d+)\s+(.*)','tokens')'; => MSG not INPUT

d = dir(strcat(fullfile(rawdatapath, 'data_zurich/*/', '*asc')));

for i = 1 : size(d, 1) 

    inputFile = fullfile(d(i).folder, d(i).name);

    filename = d(i).name(1:end-4); % remove .asc   
    x = strsplit(filename, '_');

    if strcmp(x{2}, '1') | strcmp(x{2}, '2') | strcmp(x{2}, '3') | strcmp(x{2}, '4') | strcmp(x{2}, '5')
        newfilename = [x{1}, '_ManyLabsCDA_block' x{2} '_task_ET.mat'];
    elseif strcmp(x{2}, 'Eye')
        newfilename = [x{1}, '_Eye_ET.mat'];
    elseif strcmp(x{2}, 'Res')
        newfilename = [x{1}, '_Resting_ET.mat'];
    end

    % convert
    outputFile = fullfile(d(i).folder, newfilename);
    ET = parseeyelink(inputFile, outputFile);

end

%% merge EEG and ET data, and save as mat file
d = dir(strcat(fullfile(rawdatapath, 'data_zurich/*/', '*cnt')));
result_folder = fullfile(formatted_data, 'zurich_langer');

ids = {};
for f = 1 : size(d, 1)

    filePath = fullfile(d(f).folder, d(f).name);

    % load the data
    [EEGorig, command] = pop_loadeep_v4(filePath);

    % extract path and filename
    p = strsplit(filePath, filesep);
    filePath = fullfile(filesep, p{1:end-1});
    fileName = p{end};
    p = strsplit(fileName, '_');
    subjectID = p{1};
    
    disp(['Processing subject ' subjectID])

    % exclude photodiode data
    if EEGorig.nbchan > 128
        EEGorig = pop_select(EEGorig, 'nochannel', [129:EEGorig.nbchan]);
    end

    % add Ref channel and data, and load channel location file
    EEGorig.data(129, :) = 0;
    EEGorig.nbchan = 129;
    EEGorig.chanlocs(129).labels = 'CPz'; 
    EEG.ref = 'CPz';
    locspath = 'standard_1005.elc';
    EEGorig = pop_chanedit(EEGorig, 'lookup', locspath);

    % info for BIDS format
    EEGorig.tInfo.OldID = subjectID;
    EEGorig.tInfo.EEGReference = 'Cpz';
    EEGorig.tInfo.InstitutionName = 'University of Zurich (Prof. Langer)';
    EEGorig.tInfo.PowerLineFrequency = 50;
    EEGorig.tInfo.EEGGround = 'adjacent to M1';
    EEGorig.tInfo.CapManufacturer = 'ANT Neuro';
    EEGorig.tInfo.SoftwareFilters = "n/a";
   
    
%     % sanity check
%     figure;
%     topoplot([], EEGorig.chanlocs, 'electrodes', 'labels')    

    % find start and end triggers of each task
    % start - 10
    % end - 90
    
    % subject 10, 12: missing block 1 entirely
    if strcmp(subjectID, '10') | strcmp(subjectID, '12')
        START = [9, 10, 12, 13, 14, 15];
        END = [89, 90, 92, 93, 94, 95];
        BLOCKS = [NaN, NaN, 2, 3, 4, 5];
        numTasks = length(START);
    else
        START = [9, 10, 11, 12, 13, 14, 15];
        END = [89, 90, 91, 92, 93, 94, 95];
        BLOCKS = [NaN, NaN, 1, 2, 3, 4, 5];
        numTasks = length(START);
    end
    
    % prepare to concat EEG data
    EEGconcat = {};
    BEHconcat = {};
    mkdir(fullfile(result_folder, subjectID))
    for t = 1 : numTasks
        
        disp(t)

        iSTART = find(ismember({EEGorig.event.type}, num2str(START(t))));
        iEND = find(ismember({EEGorig.event.type}, num2str(END(t))));

        % cut the data
        EEG = pop_select(EEGorig, 'point', [EEGorig.event(iSTART(end)).latency, EEGorig.event(iEND(end)).latency]);

        % check if it's resting (triggers 20 and 30 are only in Resting EEG)
        if sum(ismember({EEG.event.type}, '20')) > 3
            task = [subjectID, '_Resting_EEG.mat'];
            etfile = fullfile(d(f).folder, strcat(subjectID, '_Resting_ET.mat'));
        elseif sum(ismember({EEG.event.type}, '55')) > 3
            task = [subjectID, '_Eye_EEG.mat'];
            etfile = fullfile(d(f).folder, strcat(subjectID, '_Eye_ET.mat'));
        else
            task = [subjectID, '_ManyLabsCDA', '_block', num2str(BLOCKS(t)), '_task_EEG.mat'];
            etfile = fullfile(d(f).folder, strcat(subjectID, '_ManyLabsCDA_block', num2str(BLOCKS(t)), '_task_ET.mat'));
            load(fullfile(d(f).folder, strcat(subjectID, '_ManyLabsCDA_block', num2str(BLOCKS(t)), '_task.mat'))); % behavioral
        end

        % merge ET with EEG
        EEG = pop_importeyetracker(EEG, etfile, [START(t) END(t)],[2 3],{'L_GAZE_X', 'L_GAZE_Y'},1,1,1,0);

        % save Eye and Resting to a file or prepare to concat CDA blocks
        if not(contains(task, 'ManyLabs'))
            save(fullfile(result_folder, subjectID, task), 'EEG', '-v7.3')
        else
            % merge the blocks
            EEGconcat{t} = EEG;
            BEHconcat{t} = beh.data;  
        end


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

end