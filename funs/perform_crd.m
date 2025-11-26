function [EEG] = perform_crd(EEG, current_lab, params)
    
    % CRD
    params = params.crd;
    [EEG_CRD] = clean_artifacts(EEG, params);
    if isfield(EEG_CRD.etc, 'clean_channel_mask')
        CRD_goodChannels = EEG_CRD.etc.clean_channel_mask;
    else
        CRD_goodChannels = ones(1, EEG.nbchan);
    end

    % check if there is a reference channel in the data. It shoudnt be marked
    % as bad.
    i_ref = find(ismember({EEG.chanlocs.labels}, EEG.BIDS.tInfo.EEGReference));
    if not(isempty(i_ref))
        CRD_goodChannels(i_ref) = 1;
    end
    
    % dont remove HEOG channel in Dartmouth & Zürich (its already left - right and does not have chanlocs)
    if strcmp(current_lab, 'Dartmouth') | strcmp(current_lab, 'Zürich (Prof. Sauseng)')
        CRD_goodChannels(ismember({EEG.chanlocs.labels}, 'HEOG')) = 1;
    end
   
    % find the bad channels
    badChannelsCRD = find(not(CRD_goodChannels));
    EEG.badChannelsCRD = badChannelsCRD;

    % change to 0. ZapLine skips channels that are 0, but fails when
    % they are NaN
    if not(isempty(badChannelsCRD))
        EEG.data(badChannelsCRD, :) = 0;
    end

end