function EEG = perform_high_variance_and_interpolate(EEG, params)

    % high variance criterion
    high_var = [];
    for ch = 1 : EEG.nbchan
        % Proportion of time points exceeding 100 mV threshold
        high_amp_prop = sum(EEG.data(ch, :) > params.high_var.thresh | EEG.data(ch, :) < -params.high_var.thresh) / EEG.pnts;
        % if 30% of the data has large amplitudes, exclude
        if high_amp_prop > params.high_var.prop 
            high_var(end+1) = ch;
        end
    end

    % keep channels that dont have chanlocs, especially we need the HEOG and VEOG
    no_chanlocs = find(cellfun(@isempty, {EEG.chanlocs.X}));

    % store the bad channels info
    EEG.high_var = high_var;
    badChannels = unique([EEG.high_var(:)', EEG.badChannelsCRD(:)']);
    % remove channels without chanlocs from badChannels (they cant be interpolated)
    badChannels = setdiff(badChannels, no_chanlocs);
    EEG.badChannelsFinal = badChannels;

    % interpolate bad channels, but only if at least 50% are good.
    % we will exclude subjects with too many bad channels anyway.
    if length(badChannels) < params.high_var.reject_ratio * size(EEG.data, 1)
        EEG.data(badChannels, :) = NaN;
        EEG = eeg_interp(EEG, badChannels, 'spherical');
    end


end