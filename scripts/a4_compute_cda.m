%% compute CDA

p = pwd;
funpath = strsplit(p, filesep);
addpath(fullfile(strjoin(funpath(1:end-1), filesep), 'funs'))
initPaths;



%% load participants info
participants = readtable(fullfile(bidspath, 'participants.tsv'), 'FileType', 'text', 'Delimiter', '\t');

%% Compute CDA for all preprocessing pipelines

% list pipelines
pipelines = {preprocessedDirect, preprocessedAdvanced, preprocessedICA, preprocessedICA};
tabs = {'cda_table_direct', 'cda_table_advanced', 'cda_table_ica', 'cda_table_ica_keepall'};
% we need a table in long format with single trial CDA estimates for 
% data quality measurements 
tabs_quality = {'tab_data_quality_direct', ...
                'tab_data_quality_advanced', ...
                'tab_data_quality_ica', ...
                'tab_data_quality_ica_keepall'};

for pip = 1 : length(pipelines)

    % list files
    d = dir(fullfile(pipelines{pip}, 'sub-*'));
      
    % init empty tables
    varNames = {'ID', 'Lab', 'Age', 'Sex', 'Hand'};
    varTypes = {'cell', 'cell','double','cell', 'cell'};
    tab = table('Size',[0 length(varNames)], ...
                'VariableNames',varNames, ...
                'VariableTypes',varTypes);
    tab_data_quality =  table();
    
    % loop
    for sub = 1 : size(d, 1)
        
        % load data
        disp(d(sub).name)
        oldId = participants.old_id(sub);
        current_lab = participants.lab(sub);
        d_sub = dir(fullfile(d(sub).folder, d(sub).name, '*EEG*'));
        load(fullfile(d_sub(1).folder, d_sub(1).name))
    
        % save demog to table 
        tab.ID{sub} = d(sub).name;
        tab.Lab{sub} = participants.lab{sub};
        tab.Age(sub) = participants.age(sub);
        tab.Sex{sub} = participants.sex{sub};
        tab.Hand{sub} = participants.handedness{sub};
    
        % save eye calib thresholds to table
        tab.avgR1_plus_sd(sub) = EEG.eye_cal.avgR1_plus_sd;
        tab.avgL1_plus_sd(sub) = EEG.eye_cal.avgR1_plus_sd;
    
        % sanity check
        if not(size(EEG.data, 3) == length(EEG.beh.data.allResponses))
            error('Mismatch between number of trials and number of responses')
        end
    
        % compute accuracy for each set size
        count = 1;
        for sz = [2, 4 ,6]
            acc = nansum(EEG.beh.data.allCorrect(EEG.beh.data.trialSetSize == sz)) / length(EEG.beh.data.allCorrect(EEG.beh.data.trialSetSize == sz));
            tab.(strcat('Accuracy', num2str(sz)))(sub) = acc;
            count = count + 1;
        end
    
        % compute WM capacity
        [K, D_prime] = compute_wm_capacity(EEG, oldId, current_lab);
        tab.K2(sub) = K(1);
        tab.K4(sub) = K(2);
        tab.K6(sub) = K(3);
        tab.d_prime2(sub) = D_prime(1);
        tab.d_prime4(sub) = D_prime(2);
        tab.d_prime6(sub) = D_prime(3);
    
        %  ET data, if exist
        if isfield(EEG, 'ET')
            tab.excl_ET_Sacc(sub) = numel(EEG.excl_ET_Sacc);
            tab.excl_ET_Blink(sub) = numel(EEG.excl_ET_Blink);
        else
            tab.excl_ET_Sacc(sub) = NaN;
            tab.excl_ET_Blink(sub) = NaN;
        end
        
        % bad trials
        tab.excl_AmpThresh(sub) = numel(EEG.excl_AmpThresh);
        tab.excl_HEOG(sub) = numel(EEG.excl_HEOG_windows);
        tab.excl_VEOG(sub) = numel(EEG.excl_VEOG);
        tab.excl_blocking(sub) = numel(EEG.excl_blocking);
    
        % Bad trial rejection for direct replication (eye calibration task)
        if strcmp(tabs{pip}, 'cda_table_direct')
            excl_trial = unique([EEG.excl_AmpThresh; EEG.excl_HEOG_windows'; EEG.excl_VEOG; EEG.excl_blocking']);
            
        % bad trial rejection for advanced pipeline (ET data, if exist)
        elseif strcmp(tabs{pip}, 'cda_table_advanced') | strcmp(tabs{pip}, 'cda_table_ica')

            if isfield(EEG, 'ET')
                excl_trial = unique([EEG.excl_AmpThresh; EEG.excl_ET_Sacc'; EEG.excl_ET_Blink'; EEG.excl_blocking']);
            else
                excl_trial = unique([EEG.excl_AmpThresh; EEG.excl_HEOG_windows'; EEG.excl_VEOG; EEG.excl_blocking']);
            end
        % dont exclude any trials
        elseif strcmp(tabs{pip}, 'cda_table_ica_keepall')
            excl_trial = [];

        % else, error and double check
        else
            error('Bad trial rejection error')
        end

        % final count of removed trials
        tab.excl_total(sub) = numel(excl_trial);
        tab.total_trials(sub) = size(EEG.data, 3);
    

        % remove bad data and according task behavioural results (!!)
        EEG.data(:,:,excl_trial) = [];
        EEG.trials = size(EEG.data,3);
        stim.setsize = EEG.beh.data.trialSetSize;
        stim.SideRight = EEG.beh.data.trialCuedSide;
        stim.allCorrect = EEG.beh.data.allCorrect;
        stim.setsize(excl_trial)=[];
        stim.SideRight(excl_trial)=[];
        stim.allCorrect(excl_trial)=[];
    
        % Get electrode clusters for CDA
        Cluster_L = [];
        L1 = find(strcmp({EEG.chanlocs.labels}, 'P3'));
        L2 = find(strcmp({EEG.chanlocs.labels}, 'P7'));
        L6 = find(strcmp({EEG.chanlocs.labels}, 'O1'));
        % L7 = find(strcmp({EEG.chanlocs.labels}, 'PO3'));
        % L8 = find(strcmp({EEG.chanlocs.labels}, 'PO7'));
        Cluster_L = [L1,L2,L6];
        % Cluster_L = [L1,L2,L6,L7,L8];
        
        Cluster_R = [];
        R1 = find(strcmp({EEG.chanlocs.labels}, 'P4'));
        R2 = find(strcmp({EEG.chanlocs.labels}, 'P8'));
        R6 = find(strcmp({EEG.chanlocs.labels}, 'O2'));
        % R7 = find(strcmp({EEG.chanlocs.labels}, 'PO4'));
        % R8 = find(strcmp({EEG.chanlocs.labels}, 'PO8'));
        Cluster_R = [R1,R2,R6];
        % Cluster_R = [R1,R2,R6,R7,R8];
    
        % compute CDA for incorrect trials
        incorr = stim.allCorrect == 0;
        EEG_L_INCORR = EEG.data(:,:, find(stim.SideRight == 0 & incorr));
        EEG_R_INCORR = EEG.data(:,:,find(stim.SideRight == 1 & incorr));
    
        % left cues
        CDAL2_trialwise = nanmean(EEG_L_INCORR(Cluster_R,:,:),1) - nanmean(EEG_L_INCORR(Cluster_L,:,:),1);
        % right cues
        CDAR2_trialwise = nanmean(EEG_R_INCORR(Cluster_L,:,:),1) - nanmean(EEG_R_INCORR(Cluster_R,:,:),1);
    
        % mean
        CDAL2 = squeeze(nanmean(CDAL2_trialwise ,3));
        CDAR2 = squeeze(nanmean(CDAR2_trialwise ,3));
        CDA_INCORR = mean([CDAR2;CDAL2],1);
        tab.CDA_incorrect{sub} = CDA_INCORR;
    
        % calculate CDA: 
        % first divide EEG data by set-size and cue presentation:
        EEG_L_2 = EEG.data(:,:,find(stim.SideRight == 0&stim.setsize==2));
        EEG_L_4 = EEG.data(:,:,find(stim.SideRight == 0&stim.setsize==4));
        EEG_L_6 = EEG.data(:,:,find(stim.SideRight == 0&stim.setsize==6));
        EEG_R_2 = EEG.data(:,:,find(stim.SideRight == 1&stim.setsize==2));
        EEG_R_4 = EEG.data(:,:,find(stim.SideRight == 1&stim.setsize==4));
        EEG_R_6 = EEG.data(:,:,find(stim.SideRight == 1&stim.setsize==6));
    
        % calculate the CDA trial-wise, separate for left and right presented cues!
        % left cues
        CDAL2_trialwise = nanmean(EEG_L_2(Cluster_R,:,:) - EEG_L_2(Cluster_L,:,:), 1);
        CDAL4_trialwise = nanmean(EEG_L_4(Cluster_R,:,:) - EEG_L_4(Cluster_L,:,:),1);
        CDAL6_trialwise = nanmean(EEG_L_6(Cluster_R,:,:) - EEG_L_6(Cluster_L,:,:),1);
        CDAL2 = squeeze(nanmean( CDAL2_trialwise ,3));
        CDAL4 = squeeze(nanmean( CDAL4_trialwise ,3));
        CDAL6 = squeeze(nanmean( CDAL6_trialwise ,3));
        % right cues:
        CDAR2_trialwise = nanmean(EEG_R_2(Cluster_L,:,:) - EEG_R_2(Cluster_R,:,:),1);
        CDAR4_trialwise = nanmean(EEG_R_4(Cluster_L,:,:) - EEG_R_4(Cluster_R,:,:),1);
        CDAR6_trialwise = nanmean(EEG_R_6(Cluster_L,:,:) - EEG_R_6(Cluster_R,:,:),1);
        CDAR2 = squeeze(nanmean(CDAR2_trialwise ,3));
        CDAR4 = squeeze(nanmean(CDAR4_trialwise ,3));
        CDAR6 = squeeze(nanmean(CDAR6_trialwise ,3));
    
        CDA2 = mean([CDAR2;CDAL2],1);
        CDA4 = mean([CDAR4;CDAL4],1);
        CDA6 = mean([CDAR6;CDAL6],1);
    
        % save to table
        tab.CDAL2{sub} = CDAL2;
        tab.CDAR2{sub} = CDAR2;
        tab.CDAL4{sub} = CDAL4;
        tab.CDAR4{sub} = CDAR4;
        tab.CDAL6{sub} = CDAL6;
        tab.CDAR6{sub} = CDAR6;
        tab.CDA2{sub} = CDA2;
        tab.CDA4{sub} = CDA4;
        tab.CDA6{sub} = CDA6;
    
        % GA to reproduce Fig. 1
        tab.GAcontra1{sub} = nanmean(nanmean(EEG.data(Cluster_R,:,find(stim.SideRight == 0)), 1), 3); % ele right, stimuli left
        tab.GAcontra2{sub} = nanmean(nanmean(EEG.data(Cluster_L,:,find(stim.SideRight == 1)), 1), 3); % ele left, stimuli right
        tab.GAipsi1{sub} = nanmean(nanmean(EEG.data(Cluster_R,:,find(stim.SideRight == 1)), 1), 3); % ele right, stimuli right
        tab.GAipsi2{sub} = nanmean(nanmean(EEG.data(Cluster_L,:,find(stim.SideRight == 0)), 1), 3); % ele left, stimuli left
    
    
    
        % create a table with single trial CDA for data quality checks
        % time of interest for correlation with VWM
        TOI = find(EEG.times >= 300 & EEG.times <= 900);
    
        avg_CDAR2_contra = squeeze(mean(mean(EEG_R_2(Cluster_L,TOI,:), 2), 1));
        avg_CDAR2_ipsi = squeeze(mean(mean(EEG_R_2(Cluster_R,TOI,:), 2), 1));
        avg_CDAR4_contra = squeeze(mean(mean(EEG_R_4(Cluster_L,TOI,:), 2), 1));
        avg_CDAR4_ipsi = squeeze(mean(mean(EEG_R_4(Cluster_R,TOI,:), 2), 1));
        avg_CDAR6_contra = squeeze(mean(mean(EEG_R_6(Cluster_L,TOI,:), 2), 1));
        avg_CDAR6_ipsi = squeeze(mean(mean(EEG_R_6(Cluster_R,TOI,:), 2), 1));
    
        avg_CDAL2_contra = squeeze(mean(mean(EEG_L_2(Cluster_R,TOI,:), 2), 1));
        avg_CDAL2_ipsi = squeeze(mean(mean(EEG_L_2(Cluster_L,TOI,:), 2), 1));
        avg_CDAL4_contra = squeeze(mean(mean(EEG_L_4(Cluster_R,TOI,:), 2), 1));
        avg_CDAL4_ipsi = squeeze(mean(mean(EEG_L_4(Cluster_L,TOI,:), 2), 1));
        avg_CDAL6_contra = squeeze(mean(mean(EEG_L_6(Cluster_R,TOI,:), 2), 1));
        avg_CDAL6_ipsi = squeeze(mean(mean(EEG_L_6(Cluster_L,TOI,:), 2), 1));
    
        Contra = [ ...
            avg_CDAR2_contra; ...
            avg_CDAR4_contra; ...
            avg_CDAR6_contra; ...
            avg_CDAL2_contra; ...
            avg_CDAL4_contra; ...
            avg_CDAL6_contra ];
        
        Ipsi = [ ...
            avg_CDAR2_ipsi; ...
            avg_CDAR4_ipsi; ...
            avg_CDAR6_ipsi; ...
            avg_CDAL2_ipsi; ...
            avg_CDAL4_ipsi; ...
            avg_CDAL6_ipsi ];
    
        n_CDAR2 = numel(avg_CDAR2_contra);
        n_CDAR4 = numel(avg_CDAR4_contra);
        n_CDAR6 = numel(avg_CDAR6_contra);
        
        n_CDAL2 = numel(avg_CDAL2_contra);
        n_CDAL4 = numel(avg_CDAL4_contra);
        n_CDAL6 = numel(avg_CDAL6_contra);
    
        Condition = [ ...
            repmat({'CDAR2'}, n_CDAR2, 1); ...
            repmat({'CDAR4'}, n_CDAR4, 1); ...
            repmat({'CDAR6'}, n_CDAR6, 1); ...
            repmat({'CDAL2'}, n_CDAL2, 1); ...
            repmat({'CDAL4'}, n_CDAL4, 1); ...
            repmat({'CDAL6'}, n_CDAL6, 1) ...
        ];
    
        N = n_CDAR2 + n_CDAR4 + n_CDAR6 + n_CDAL2 + n_CDAL4 + n_CDAL6;
    
        % ID and Lab repeated for all rows
        ID_vec  = repmat({tab.ID{sub}},  N, 1);  
        Lab_vec = repmat({tab.Lab{sub}}, N, 1);   
        
        % Final table
        tmp_tab = table(ID_vec, Lab_vec, Condition, Contra, Ipsi, ...
            'VariableNames', {'ID', 'Lab', 'Condition', 'Contra', 'Ipsi'});
    
        % append to big table
        tab_data_quality = [tab_data_quality; tmp_tab];
    
        
    end
    
    % save the results
    save(fullfile(datapath, 'mat_files', strcat(tabs{pip}, '.mat')), 'tab')
    save(fullfile(datapath, 'mat_files', strcat(tabs_quality{pip}, '.mat')), 'tab_data_quality')
    writetable(tab_data_quality, fullfile(datapath, 'csv_files', strcat(tabs_quality{pip}, '.csv')))

end







