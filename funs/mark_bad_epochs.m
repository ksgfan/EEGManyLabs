function [EEG] = mark_bad_epochs(EEG, HEOG, VEOG, eye_cal, params)


    % exclude trials based on calibration task: 2 windows approach   
    winpnts = params.heog.winpnts; % window in time points. each window width is half of this value
    stepnts = params.heog.stepnts; % step in time points
  
    excl_HEOG = [];
    for tr=1:size(HEOG.data,3)
        for j=1:stepnts:HEOG.pnts-(winpnts-1)
            w1  = HEOG.data(1, j:j+round(winpnts/2)-1,tr);
            w2  = HEOG.data(1, j+round(winpnts/2):j+winpnts-1 ,tr);
            if mean(w2) - mean(w1) < mean(eye_cal.avgL1_minus_sd) | mean(w2) - mean(w1) > mean(eye_cal.avgR1_plus_sd)
                excl_HEOG(tr) = 1;     
                break
            end
        end
    end
    EEG.excl_HEOG_windows = find(excl_HEOG);
    clear excl_HEOG
    
    % blink rejection: 90 mV threshold on the UNIPOLAR VEOG channel
    i_VEOGU = ismember({VEOG.chanlocs.labels}, 'VEOGU');
    excl_VEOG = max(VEOG.data(i_VEOGU, :, :), [], 2) > params.veog | min(VEOG.data(i_VEOGU, :, :), [], 2) < -params.veog;
    EEG.excl_VEOG = find(excl_VEOG);
    clear excl_VEOG

    % Adam 2018: blocking artifacts
    winpnts = params.blocking.winpnts; % 200 ms window at 250 Hz
    stepnts = params.blocking.stepnts; % 50 ms
    flat_samples = params.blocking.flat_samples; % 60 ms
    range_thresh = params.blocking.range_thresh;       
    chans = ismember({EEG.chanlocs.labels}, {'P3', 'P7', 'O1', 'P4', 'P8', 'O2'});
    excl_blocking = zeros(1, size(EEG.data, 3));
    % Loop through epochs
    for e = 1:size(EEG.data, 3)
        epoch_data = EEG.data(chans, :, e);
        for w = 1:stepnts:(size(epoch_data, 2) - winpnts + 1)
            window_data = epoch_data(:, w:(w + winpnts - 1));
            for ch = 1:size(window_data, 1)
                chan_data = window_data(ch, :);
                for fs = 1:(length(chan_data) - flat_samples + 1)
                    segment = chan_data(fs:(fs + flat_samples - 1));
                    % Check if the segment is flat
                    if range(segment) < range_thresh
                        excl_blocking(e) = 1;
                        break;
                    end
                end
                if excl_blocking(e), break; end
            end
            if excl_blocking(e), break; end
        end
    end
    EEG.excl_blocking = find(excl_blocking);
    clear excl_blocking
       
    % exclude tials with saccades larger than 1 degree and blinks, if
    % there is ET data aviaible
    if isfield(EEG, 'ET')
        types = {EEG.event.type};
        epochs = {EEG.event.epoch};
        sacc_amp = {EEG.event.sac_amplitude};
        epochsWithSaccades = unique(cell2mat(epochs(find(contains(types, '_saccade') & cell2mat(sacc_amp) > 1))));
        epochsWithAllSaccades = unique(cell2mat(epochs(find(contains(types, '_saccade')))));
        epochsWithBlinks = unique(cell2mat(epochs(find(contains(types, '_blink')))));
        EEG.excl_ET_Sacc = epochsWithSaccades;
        EEG.all_ET_Sacc = epochsWithAllSaccades;
        EEG.excl_ET_Blink = epochsWithBlinks;
        clear epochsWithSaccades epochsWithAllSaccades epochsWithBlinks
    end
    

    % exclude segments if any electrode of interest (HEOG, P3, T5/P7, O1; P4,
    % T6/P8, O2) had a peak to peak amplitude > 200 mV
    chans = ismember({EEG.chanlocs.labels}, {'HEOG', 'P3', 'P7', 'O1', 'P4', 'P8', 'O2'});
    threshold = params.highamp;
    diffs = squeeze( max(EEG.data(chans, :, :),[],2) - min(EEG.data(chans, :, :),[],2) );
    [chanidx, epochidx] = find(diffs > threshold);
    EEG.excl_AmpThresh = unique(epochidx);
    clear epochidx

end