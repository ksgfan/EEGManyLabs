%% preprocess and identify bad epochs

p = pwd;
funpath = strsplit(p, filesep);
addpath(fullfile(strjoin(funpath(1:end-1), filesep), 'funs'))
initPaths;

%% list files
d = dir(fullfile(bidspath, 'sub-*'));

% load participants info
participants = readtable(fullfile(bidspath, 'participants.tsv'), 'FileType', 'text', 'Delimiter', '\t');

%% Params
params = struct;
params.heog.winpnts = 50; % window in time points. each window width is half of this value; 200 ms window at 250 Hz
params.heog.stepnts = 2; % step in time points
params.veog = 90;
params.blocking.winpnts = 50; 
params.blocking.stepnts = 12; % 50 ms
params.blocking.flat_samples = 15; % 60 ms
params.blocking.range_thresh = 0.1;    
params.highamp = 200;
params.crd.FlatlineCriterion = 5; %seconds
params.crd.LineNoiseCriterion = 4; %standard deviations
params.crd.ChannelCriterion = 0.85; %correlation
params.crd.BurstCriterion = 'off';
params.crd.WindowCriterion = 'off';
params.filter.high = 0.01;
params.filter.low = 80;
params.ica.thresh = 0.8;
params.high_var.thresh = 100;
params.high_var.prop = 0.3;
params.high_var.reject_ratio = 0.5;

%% preprocessing

for sub = 1 : size(d, 1)
    
    % clear variables
    clear ET EEG EYE d_beh beh
    
    % loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(d(sub).name)
    current_lab = participants.lab(sub);
    d_beh = dir(fullfile(d(sub).folder, d(sub).name, 'beh', '*json*'));
    % load data
    [STUDY, ALLEEG] = pop_importbids_dawid(bidspath, 'subjects', {d(sub).name});
    if length(ALLEEG) == 2
        EEG = ALLEEG(ismember({ALLEEG.task}, 'CDA'));
        EYE = ALLEEG(ismember({ALLEEG.task}, 'Eye'));
        % downsample
        if EEG.srate > 250
            EEG = pop_resample(EEG, 250);
            EYE = pop_resample(EYE, 250);
        end
    else
        EEG = ALLEEG(ismember({ALLEEG.task}, 'CDA'));
        % downsample
        if EEG.srate > 250
            EEG = pop_resample(EEG, 250);
        end
    end

    % load behavioral
    fid = fopen(fullfile(d_beh(1).folder, d_beh(1).name));
    raw = fread(fid, inf);
    fclose(fid);
    EEG.beh = jsondecode(char(raw'));

    % If the eye-calibration task is missing, create an empty struct so the 
    % downstream code can still run without errors
    if not(exist('EYE'))
        EYE.chanlocs = EEG.chanlocs;
    end
    
    % rename channels
    [EEG, EYE] = rename_channels(EEG, EYE, current_lab);
   
    % compute saccade thresholds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % additionally, the eye calibration task failed at some labs. set fixed
    % thresholds for now.
    if length(ALLEEG) == 2 & ...
         (not(strcmp(current_lab, 'Ohio')) & ...
          not(strcmp(current_lab, 'Zürich (Prof. Sauseng)')) & ...
          not(strcmp(current_lab, 'Florida')) ...
         )
        
        % compute thresholds
        eye_cal = eyeMovement(EYE, current_lab);
    else
        % if eye calibration is missing...
        eye_cal = struct;
        eye_cal.avgR1_plus_sd = 30;
        eye_cal.avgL1_minus_sd = -30;
        eye_cal.avgR1_95CI = 50;
        eye_cal.avgL1_95CI = -50;
    end
    
    % save info calibration task thresholds
    EEG.eye_cal = eye_cal;
    
    % Extract ET data
    ET_chan = {'L-GAZE-X', 'L-GAZE-Y'};
    if any(ismember({EEG.chanlocs.labels}, ET_chan))
        ET = pop_select(EEG, 'channel', ET_chan);
        % remove ET from EEG data
        EEG = pop_select(EEG, 'nochannel', [ET_chan]);
    end

    % bad channel detection: clean raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EEG = perform_crd(EEG, current_lab, params);

    % filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [EEG, fig1] = filter_data(EEG, current_lab, params);

    % ICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create a duplicate EEG structure. 
    % first will be used for "advanced" preprocessing, as in stage 1
    % second will undergo ICA as well (for supplement)
    EEGfiltered = EEG;
    
    % perform ICA
    EEGica = perform_ica(EEG, params);
    
    % % sanity check
    % EEGica.icasphere = EEGica.etc.icasphere_beforerms;
    % EEGica.icaweights = EEGica.etc.icaweights_beforerms;
    % EEGica.icawinv = pinv(EEGica.icaweights * EEGica.icasphere); 
    % EEGica.icaact = icaact(EEGica.data(EEGica.icachansind, :, :),EEGica.icaweights,0);
    % EEGica.etc.ic_classification.ICLabel.classifications = EEGica.etc.ic_classification.ICLabel.all_classifications;
    % addpath(fullfile(eeglabpath, 'plugins/ICLabel1.7/viewprops'))
    % pop_viewprops(EEGica, 0)

    % % plot before and after
    % figure;
    % plot(EEGfiltered.times, EEGfiltered.data(2, :))
    % hold on
    % plot(EEGica.times, EEGica.data(2, :))
    % legend('Fitlered', 'ICA')

    % struct with datasets
    EEGdatasets = {EEGfiltered, EEGica};
    
    % paths with result folder
    paths = {preprocessedAdvanced, preprocessedICA};
    
    % loop over them
    for dset = 1 : length(EEGdatasets)
        
        % select dataset
        EEG = EEGdatasets{dset};

        % remove channels with high variance and interpolate all bad chans
        EEG = perform_high_variance_and_interpolate(EEG, params);

        % rereference to algebraic average of left and right mastoid %%%%%%%%%%
        el_m1 = find(strcmp({EEG.chanlocs.labels}, 'M1')); % M1
        el_m2 = find(strcmp({EEG.chanlocs.labels}, 'M2')); % M2
        % some labs dont have M1/M2 - use T7/T8
        if isempty(el_m1) | isempty(el_m2)
            el_m1 = find(strcmp({EEG.chanlocs.labels}, 'T7')); % M1
            el_m2 = find(strcmp({EEG.chanlocs.labels}, 'T8')); % M2
        end
        EEG = pop_reref(EEG, [el_m1 el_m2], 'keepref', 'on');
        clear el_m1 el_m2

        % extract HEOG and VEOG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(current_lab, 'Dartmouth') | strcmp(current_lab, 'Zürich (Prof. Sauseng)')
            HEOG_chan = {'HEOG'};
        else
            HEOG_chan = {'HEOGR', 'HEOGL'};
        end
        HEOG = pop_select(EEG, 'channel', HEOG_chan);
        % filter - saccadic thresholds were computed on filtered data as well
        HEOG = pop_eegfiltnew(HEOG, 0.1, 40);

        VEOG_chan = {'VEOGU'};
        VEOG = pop_select(EEG, 'channel', VEOG_chan);
        % filter
        VEOG = pop_eegfiltnew(VEOG, 0.1, 40);

        % imagesc of preprocessed data
        fig2 = figure('Visible', 'off');
        imagesc(EEG.data, [-100, 100])
        mycolormap = customcolormap_preset('red-white-blue');
        colormap(mycolormap)
        colorbar; 

        % epoch data (also VEOG and HEOG for trial exclusion)
        % before: prepare HEOG and EOG channels:
        if  not(strcmp(current_lab, 'Dartmouth') | strcmp(current_lab, 'Zürich (Prof. Sauseng)'))
            HEOG.datasave = HEOG.data;
            HEOG.data = HEOG.data(1, :) - HEOG.data(2, :);
            HEOG.nbchan = 1;
            HEOG.chanlocs(2) = [];
        end

        % epoch data -200 to 1200 ms 
        HEOG = pop_epoch(HEOG, {21, 41, 61}, [-0.2 1.2]);
        VEOG = pop_epoch(VEOG, {21, 41, 61}, [-0.2 1.2]);
        EEG = pop_epoch(EEG, {21, 41, 61}, [-0.2 1.2]);
        if exist('ET')
            ET = pop_epoch(ET, {21, 41, 61}, [-0.2 1.2]);
            % save ET data in EEG structure
            EEG.ET = ET;
        end

        % baseline correction
        EEG = pop_rmbase(EEG, [-200, 0]);
        HEOG = pop_rmbase(HEOG, [-200, 0]);
        VEOG = pop_rmbase(VEOG, [-200, 0]);
        
        % mark bad epochs
        EEG = mark_bad_epochs(EEG, HEOG, VEOG, eye_cal, params);

        % save
        newPath = fullfile(paths{dset}, d(sub).name);
        mkdir(newPath)
        save(fullfile(newPath, 'EEG.mat'), 'EEG', '-v7.3')
        saveas(fig1, fullfile(newPath, 'zapline.png'));
        saveas(fig2, fullfile(newPath, 'imagesc_plot.png'));
    end
    close all
end
