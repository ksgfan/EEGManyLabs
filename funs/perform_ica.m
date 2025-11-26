function [EEG] = perform_ica(EEG, params)

    % keep orig data (without 2 Hz HP-filter)
    EEGorig = EEG;

    % high pass for ICA
    EEG = pop_eegfiltnew(EEG, 2, 0);

    % exclude reference and bad channels from ICA
    i_ref = find(ismember({EEG.chanlocs.labels}, EEG.BIDS.tInfo.EEGReference));
    EEG.icachansind = 1 : EEG.nbchan;
    exclude_chans = [EEG.badChannelsCRD, i_ref];
    EEG.icachansind(exclude_chans) = [];

    % run ica
    EEG = pop_runica(EEG, 'icatype', 'runica', 'chanind', EEG.icachansind);

    % run iclabel
    EEG = iclabel(EEG);

    % save all classes
    EEG.etc.ic_classification.ICLabel.all_classifications = EEG.etc.ic_classification.ICLabel.classifications;

    % replace with non-filtered data
    EEG.data = EEGorig.data;
    
    % find comps:
    muscleComponents = find(EEG.etc.ic_classification.ICLabel.classifications(:, 2) > params.ica.thresh);
    eyeComponents = find(EEG.etc.ic_classification.ICLabel.classifications(:, 3) > params.ica.thresh);
    heartComponents = find(EEG.etc.ic_classification.ICLabel.classifications(:, 4) > params.ica.thresh);
    lineNoiseComponents = find(EEG.etc.ic_classification.ICLabel.classifications(:, 5) > params.ica.thresh);
    channelNoiseComponents = find(EEG.etc.ic_classification.ICLabel.classifications(:, 6) > params.ica.thresh);


    uni_comps = {muscleComponents, eyeComponents, ...
        heartComponents, lineNoiseComponents, channelNoiseComponents};
    components = unique(cat(1, uni_comps{:}));
    
    % remove comps
    EEG.icaact = []; % let the eeglab recompute the icaact using not filtered data
    
    % if no bad comps were found
    if not(isempty(eyeComponents))
        EEG = pop_subcomp(EEG, components);
    else
        EEG = EEG;
    end
end