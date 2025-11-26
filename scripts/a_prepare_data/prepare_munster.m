

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
    EEG.tInfo.CapManufacturer = 'Biosemi';
    EEG.tInfo.SoftwareFilters = "n/a";
    
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
    
    % append  behavioral
    EEG.beh.data = newBEH; % replace with merged data
    mkdir(fullfile(result_folder, id))
    save(fullfile(result_folder, id, [id '_CDA_EEG.mat']), 'EEG', '-v7.3')

        
end

    
%%

