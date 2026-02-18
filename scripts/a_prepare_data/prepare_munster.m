

%% description
% 73 channels
% reference free (biosemi)
% only CDA; no resting, no eye calibration task
% ET already merged

% Wir haben darin folgende Elektroden:

% 64 EEG-Elektroden
    % Channel 1 = Fp1
    % Channel 2-64 alle weiteren EEG-Elektroden
% 3 EOG-Elektroden:
    % Channel 65 = Afp9
    % Channel 66 = Afp10
    % Channel 67 = IO1
% 2 Mastoid-Elektroden 
    % Channel 68 = M1
    % Channel 69 = M2
% 2 bipolare EOGs, berechnet als:
    % Channel 70 = VEOG, d.h.: Channel 1 (Fp1) minus channel 67 (IO1)
    % Channel 71 = HEOG, d.h. Channel 65 (Afp9) minus channel 66 (Afp10)
% 2 Kanäle mit den Eye-Tracking Daten
    % Channel 72 = Eyegaze-X
    % Channel 73 = Eyegaze-Y


% subject 7 has only 608 trials, but all 720 responses - remove responses
% subject 12 has 580 trials and 720 responses (all above 576 are NaNs - remove)
% subject 27 has 438 trials and 720 responses (all above 432 are NaNs - remove)

% subject 10: O1 noisy
% subject 13 is very noisy (amps > 1000 mV) => from notes: M1 and M2 were very noisy
% subject 14: P3 and O1 are noisy
% subject 15: P7 noisy
% subject 20: P3 noisy
% subject 23: P3, P4, P7, P8, O1,O2 very noisy => from notes: M1 and M2 were very noisy
% we have to exclude them from the direct replication

%% load data
data_path = fullfile(rawdatapath, 'data_munster');
d = dir([data_path, filesep, '*',filesep, '*EEG*.mat']);
result_folder = fullfile(formatted_data, 'munster');
mkdir(result_folder)

for f = 1 : size(d, 1)
   
    id = strsplit(d(f).folder, filesep);
    id = id{end};
    if startsWith(id, '0')
        id = id(2);
    end
    
    disp(id)
    
    load(fullfile(d(f).folder, d(f).name));
    
    % rename channels
    % EEG.chanlocs(65).labels = 'HEOGL';
    % EEG.chanlocs(66).labels = 'HEOGR';
    % 
    % EEG.chanlocs(1).labels = 'VEOGU'; % above left
    % EEG.chanlocs(67).labels = 'VEOGL'; % below left
    
    % info for BIDS format
    EEG.tInfo.OldID = id;
    EEG.tInfo.EEGReference = 'Reference free';
    EEG.tInfo.InstitutionName = 'University of Münster';
    EEG.tInfo.PowerLineFrequency = 50;
    EEG.tInfo.EEGGround = 'adjacent to POz';
    
    EEG.tInfo.Manufacturer = "Biosemi";
    EEG.tInfo.ManufacturersModelName =;
    EEG.tInfo.CapManufacturer = "Biosemi";
    EEG.tInfo.CapManufacturersModelName =;

    EEG.tInfo.HardwareFilters.Highpass=struct('CutoffFrequency','n/a','Description','No hardware high-pass filter');
    EEG.tInfo.HardwareFilters.Lowpass=struct('CutoffFrequency',200,'Description','Built-in ADC sinc anti-aliasing filter (~0.4 x Fs at 500 Hz, fixed, not user configurable)');
    EEG.tInfo.HardwareFilters.Notch=struct('CutoffFrequency','n/a','Description','No hardware notch filter');
    EEG.tInfo.SoftwareFilters.Highpass=struct('CutoffFrequency','n/a','Description','No software high-pass filter');
    EEG.tInfo.SoftwareFilters.Lowpass=struct('CutoffFrequency','n/a','Description','No software low-pass filter');
    EEG.tInfo.SoftwareFilters.Notch=struct('CutoffFrequency','n/a','Description','No software notch filter');
    
    EEG.chanlocs(72).labels = 'L-GAZE-X'; 
    EEG.chanlocs(73).labels = 'L-GAZE-Y'; 
    
    % rename ET events to blink, saccade, fixation
    for ev = 1 : length(EEG.event)
        if strcmp(EEG.event(ev).type, '85') 
            EEG.event(ev).type = 'L_saccade';
        elseif strcmp(EEG.event(ev).type, '81')
            EEG.event(ev).type = 'L_fixation';
        elseif strcmp(EEG.event(ev).type, '89')
            EEG.event(ev).type = 'L_blink';
        end
    end
    
    % remove precomputed HEOG and VEOG
    EEG = pop_select(EEG, 'nochannel', [70, 71]);
    
    % % sanity check
    % figure;
    % topoplot([], EEG.chanlocs, 'electrodes', 'labels')
    
    % load and merge CDA behavioral data 
    d_beh = dir(fullfile(d(f).folder, '*.csv'));

    beh_csv = readtable(fullfile(d_beh(1).folder, d_beh(1).name)); % behavioral
    newBEH = struct;
    newBEH.allResponses = beh_csv.allResponses;
    newBEH.allCorrect = beh_csv.allCorrect;
    newBEH.trialSetSize = beh_csv.trialSetSize';
    newBEH.trialIfChange = beh_csv.trialIfChange';
    newBEH.trialCuedSide = beh_csv.trialCuedSide';
    newBEH.trialSOA = beh_csv.trialSOA';
    newBEH.TrialID = beh_csv.TrialID;
    
    % subject 7 has only 608 trials
    if strcmp(id, '7') | strcmp(id, '12') | strcmp(id, '27')
        i_events = find(ismember({EEG.event.type}, {'21', '41', '61'}));
        num_trials = length(i_events);
        
        % adjust beh data (remove all the nans)
        newBEH.allResponses = beh_csv.allResponses(1:num_trials);
        newBEH.allCorrect = beh_csv.allCorrect(1:num_trials);
        newBEH.trialSetSize = beh_csv.trialSetSize(1:num_trials)';
        newBEH.trialIfChange = beh_csv.trialIfChange(1:num_trials)';
        newBEH.trialCuedSide = beh_csv.trialCuedSide(1:num_trials)';
        newBEH.trialSOA = beh_csv.trialSOA(1:num_trials)';
        newBEH.TrialID = beh_csv.TrialID(1:num_trials);
    end

    % the data was recorded using the old version of the paradigm. We have
    % to change trigger values so that they are the same as for other labs.
    %TASK_START: 11, 12, .. 15
    %TASK_END: 91, ..95
    %CUE_LEFT: 3
    %CUE_RIGHT: 7
    %SETSIZE2: 21
    %SETSIZE4: 41
    %SETSIZE6: 61
    %RETENTION: 50
    %TEST2: 22
    %TEST4: 42
    %TEST6: 62
    %RESP_SAME_CORR: 76
    %RESP_DIFF_CORR: 77
    %RESP_SAME_INCORR: 78
    %RESP_DIFF_INCORR: 79

    % In Münster, we have:
    unique({EEG.event.type})
    % where:
    %TASK_START: always 10 
    %TASK_END: always 90
    %CUE_LEFT: 3
    %CUE_RIGHT: 7
    %RETENTION: 50
    %TEST: 27, 28, 47, 48, 67, 68

    %RESP_SAME: 77
    %RESP_DIFF: 78
    
    % 250, 80, - block start - can be removed
    % 91 - trial end - can be removed
    % remove them
    badTypes = {'250', '80', '91'};
    rm_idx = ismember({EEG.event.type}, badTypes);
    EEG.event(rm_idx) = [];

    % create new vars
    [EEG.event.type_old] = EEG.event.type;
    [EEG.event.type_new] = EEG.event.type;

    % blocs should start with triggers 11-15 and end with 91-95
    idx10 = find(ismember({EEG.event.type}, '10'));
    idx90 = find(ismember({EEG.event.type}, '90'));
    for k = 1:numel(idx10)
        EEG.event(idx10(k)).type_new = num2str(10 + k);
    end

    for k = 1:numel(idx90)
        EEG.event(idx90(k)).type_new = num2str(90 + k);
    end

    % find 'test' indices
    test_idx = find(ismember({EEG.event.type}, {'27', '28', '47', '48', '67', '68'}));
    test_idx2 = find(ismember({EEG.event.type}, {'27', '28'}));
    for i = test_idx2
        EEG.event(i).type_new = '22';
    end
    test_idx4 = find(ismember({EEG.event.type}, {'47', '48'}));
    for i = test_idx4
        EEG.event(i).type_new = '42';
    end
    test_idx6 = find(ismember({EEG.event.type}, {'67', '68'}));
    for i = test_idx6
        EEG.event(i).type_new = '62';
    end


    % responses
    old_resp_idx = find(ismember({EEG.event.type}, {'77', '78'}));
    % newBEH.allResponses % NaN means that response is missing
    missing_resp = isnan(newBEH.allResponses);
    % remove missing resp from beh file
    allResp = newBEH.allResponses(not(missing_resp));
    allCorr = newBEH.allCorrect(not(missing_resp));
    % now the three vectors have the same length
    for resp = 1 : length(old_resp_idx)

        old_resp = EEG.event(old_resp_idx(resp)).type_old;
        
        % if RESP_SAME
        if strcmp(old_resp, '77')
            if allCorr(resp) == 1
                EEG.event(old_resp_idx(resp)).type_new = '76'; % same, correct
            elseif allCorr(resp) == 0
                EEG.event(old_resp_idx(resp)).type_new = '78'; % same, incorrect
            end
        % if RESP_DIFF
        elseif strcmp(old_resp, '78')
            if allCorr(resp) == 1
                EEG.event(old_resp_idx(resp)).type_new = '77'; % diff, correct
            elseif allCorr(resp) == 0
                EEG.event(old_resp_idx(resp)).type_new = '79'; % diff, incorrect
            end
        end
    end

    % sanity check
    unique({EEG.event.type_new})

    % replace type
    [EEG.event.type] = EEG.event.type_new;

    % append  behavioral
    EEG.beh.data = newBEH; % replace with merged data
    mkdir(fullfile(result_folder, id))
    save(fullfile(result_folder, id, [id '_CDA_EEG.mat']), 'EEG', '-v7.3')

        
end

    
%%

