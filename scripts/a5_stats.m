%% stats, ANOVA, t-tests, correlations

p = pwd;
funpath = strsplit(p, filesep);
addpath(fullfile(strjoin(funpath(1:end-1), filesep), 'funs'))
initPaths;


%% load participants info
participants = readtable(fullfile(bidspath, 'participants.tsv'), 'FileType', 'text', 'Delimiter', '\t');


%% loop over all preprocessing pipelines

% list pipelines
pipelines = {preprocessedDirect, preprocessedAdvanced, preprocessedICA, preprocessedICA};
tabs = {'cda_table_direct', 'cda_table_advanced', 'cda_table_ica', 'cda_table_ica_keepall'};
pipes_labels = {'direct', 'advanced', 'ica', 'keep_all'};

for pip = 1 : length(pipelines)

    % load data
    load(fullfile(datapath, 'mat_files', strcat(tabs{pip}, '.mat')))
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % number of bad trials, exclude subjects
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    labs = unique(tab.Lab);
    exclude_vars = contains(tab.Properties.VariableNames,'excl_','IgnoreCase',true);
    exclude_vars = tab.Properties.VariableNames(exclude_vars);
    % loop over labs and compute:
    % avg number of rejected trials per excl. criterion
    % number of rejected subjects per lab
    % avg K score for set size 2, 4, 6
    % avg d prime for set size 2, 4, 6
    
    excl_tab = table;
    perf_tab = table;
    for lab = 1 : length(labs)
        excl_tab.Lab{lab} = labs{lab};
        perf_tab.Lab{lab} = labs{lab};

        idx = strcmp(tab.Lab, labs{lab});
        lab_tab = tab(idx, :);
        for e = 1 : length(exclude_vars)
            excl_tab.(strcat('sub_rejected_', exclude_vars{e}))(lab) = sum( (lab_tab.(exclude_vars{e}) ./ lab_tab.total_trials ) > 0.3); 
            excl_tab.(strcat('avg_', exclude_vars{e}))(lab) = round(mean(lab_tab.(exclude_vars{e})), 2); 
            excl_tab.(strcat('sd_', exclude_vars{e}))(lab) = round(std(lab_tab.(exclude_vars{e})), 2); 
        end
    
        excl_tab.N{lab} = size(lab_tab, 1);
        perf_tab.N{lab} = size(lab_tab, 1);

        % find number of subjects with low performance (K2 < 0)
        excl_tab.sub_rejected_perf(lab) = sum(lab_tab.K2 < 0);
    
        % Accuracy  
        perf_tab.avgAcc2(lab) = round(mean(lab_tab.Accuracy2), 2);
        perf_tab.sdAcc2(lab) = round(std(lab_tab.Accuracy2), 2);
        perf_tab.avgAcc4(lab) = round(mean(lab_tab.Accuracy4), 2);
        perf_tab.sdAcc4(lab) = round(std(lab_tab.Accuracy4), 2);
        perf_tab.avgAcc6(lab) = round(mean(lab_tab.Accuracy6), 2);
        perf_tab.sdAcc6(lab) = round(std(lab_tab.Accuracy6), 2);
        % K
        perf_tab.avgK2(lab) = round(mean(lab_tab.K2), 2);
        perf_tab.sdK2(lab) = round(std(lab_tab.K2), 2);
        perf_tab.avgK4(lab) = round(mean(lab_tab.K4), 2);
        perf_tab.sdK4(lab) = round(std(lab_tab.K4), 2);
        perf_tab.avgK6(lab) = round(mean(lab_tab.K6), 2);
        perf_tab.sdK6(lab) = round(std(lab_tab.K6), 2);
        % D'
        perf_tab.avgD2(lab) = round(mean(lab_tab.d_prime2), 2);
        perf_tab.sdD2(lab) = round(std(lab_tab.d_prime2), 2);
        perf_tab.avgD4(lab) = round(mean(lab_tab.d_prime4), 2);
        perf_tab.sdD4(lab) = round(std(lab_tab.d_prime4), 2);
        perf_tab.avgD6(lab) = round(mean(lab_tab.d_prime6), 2);
        perf_tab.sdD6(lab) = round(std(lab_tab.d_prime6), 2);
    end
    
    save(fullfile(datapath, 'mat_files', strcat('excl_tab_', pipes_labels{pip}, '.mat')), 'excl_tab')
    save(fullfile(datapath, 'mat_files', strcat('perf_tab_', pipes_labels{pip} ,'.mat')), 'perf_tab')
    
    % find subjects that should be excluded and exclude them
    excl = tab.K2 < 0 | ...                             % performance
           (tab.excl_total ./ tab.total_trials) > 0.3;  % bad trials
    
    % sanity check
    sum(excl)
    
    % table with clean data
    if strcmp(pipes_labels{pip}, 'keep_all')
        clean_tab = tab;
    else
        clean_tab = tab(not(excl), :);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute subject averages for models in R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % toi
    TOI = find(EEG.times >= 300 & EEG.times <= 900);
    % init table
    subject_averages =  table();
    % exxtract CDA
    CDA2 = cell2mat(clean_tab.CDA2);
    CDA4 = cell2mat(clean_tab.CDA4);
    CDA6 = cell2mat(clean_tab.CDA6);
    % add values to table
    subject_averages.ID = clean_tab.ID;
    subject_averages.Lab = clean_tab.Lab;
    subject_averages.Hand = clean_tab.Hand;
    subject_averages.Sex = clean_tab.Sex;
    subject_averages.CDA2 = mean(CDA2(:, TOI), 2);
    subject_averages.CDA4 = mean(CDA4(:, TOI), 2);
    subject_averages.CDA6 = mean(CDA6(:, TOI), 2);
    subject_averages.K = (clean_tab.K2 + clean_tab.K4 + clean_tab.K6) / 3;
    subject_averages.D_prime = (clean_tab.d_prime2 + clean_tab.d_prime4 + clean_tab.d_prime6) / 3;
    % save
    writetable(subject_averages, fullfile(datapath, 'csv_files', strcat('subject_averages_', pipes_labels{pip}, '.csv')))
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stats 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % correltion between K and D'
    [rho, pval] = corr(subject_averages.K, subject_averages.D_prime)
    
    % init table
    all_stats = table;
    
    % load 1 file for times
    load(fullfile(pipelines{pip}, 'sub-Vog001/EEG.mat'))
    labs = unique(participants.lab);
    
    % time of interest for correlation with VWM
    TOI = find(EEG.times >= 300 & EEG.times <= 900);
    
    for lab = 1 : length(labs)
        
        disp(labs{lab})
        all_stats.Lab{lab} = labs{lab};
    
        % extract lab data
        lab_tab = clean_tab(ismember(clean_tab.Lab, labs{lab}), :);
    
        % extract vars of interest
        GAcontra1 = cell2mat(lab_tab.GAcontra1);
        GAcontra2 = cell2mat(lab_tab.GAcontra2);
        GAipsi1 = cell2mat(lab_tab.GAipsi1);
        GAipsi2 = cell2mat(lab_tab.GAipsi2);
        CDA2 = cell2mat(lab_tab.CDA2);
        CDA4 = cell2mat(lab_tab.CDA4);
        CDA6 = cell2mat(lab_tab.CDA6);
        CDA_INCORR = cell2mat(lab_tab.CDA_incorrect);
        
        % compute means
        meanGAcontra1 = nanmean(GAcontra1, 1);
        meanGAcontra2 = nanmean(GAcontra2, 1);
        meanGAipsi1 = nanmean(GAipsi1, 1);
        meanGAipsi2 = nanmean(GAipsi2, 1);
        GAcontra = [GAcontra1 + GAcontra2] / 2;
        GAipsi = [GAipsi1 + GAipsi2] / 2;
        contra = nanmean(GAcontra(:, TOI), 2);
        ipsi = nanmean(GAipsi(:, TOI), 2);
        meanCDA2 = nanmean(CDA2, 1);
        meanCDA4 = nanmean(CDA4, 1);
        meanCDA6 = nanmean(CDA6, 1);
        % average over ROI
        meanROI2 = nanmean(CDA2(:, TOI), 2);
        meanROI4 = nanmean(CDA4(:, TOI), 2);
        meanROI6 = nanmean(CDA6(:, TOI), 2);
        % average over Subjects
        meanSUB2 = nanmean(meanROI2);
        meanSUB4 = nanmean(meanROI4);
        meanSUB6 = nanmean(meanROI6);
        % compute amplitude increase from 2 to 4 items and memory capacity
        ampIncrease_2_to_4 = meanROI4 - meanROI2;
        memoryCapacity = (lab_tab.K2 + lab_tab.K4 + lab_tab.K6) / 3; % average over all setsizes
        % 95% confidence interval
        ciCDA2 = 1.96 * (nanstd(meanROI2) / sqrt(size(meanROI2, 1))); % 95% CI
        ciCDA4 = 1.96 * (nanstd(meanROI4) / sqrt(size(meanROI4, 1))); % 95% CI
        ciCDA6 = 1.96 * (nanstd(meanROI6) / sqrt(size(meanROI6, 1))); % 95% CI
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % filter used only for visualization
        filt = 35;
        meanGAcontra1filt=eegfilt(meanGAcontra1,250,0,filt);
        meanGAcontra2filt=eegfilt(meanGAcontra2,250,0,filt);
        meanGAipsi1filt=eegfilt(meanGAipsi1,250,0,filt);
        meanGAipsi2filt=eegfilt(meanGAipsi2,250,0,filt);
        meanCDA2filt=eegfilt(meanCDA2,250,0,filt);
        meanCDA4filt=eegfilt(meanCDA4,250,0,filt);
        meanCDA6filt=eegfilt(meanCDA6,250,0,filt);
    
        % % contra vs ipsi
        % f1 = figure;
        % plot(EEG.times, meanGAcontra1filt);set(gca,'Ydir','reverse')
        % hold on
        % plot(EEG.times, meanGAcontra2filt);set(gca,'Ydir','reverse')
        % plot(EEG.times, meanGAipsi1filt);set(gca,'Ydir','reverse')
        % plot(EEG.times, meanGAipsi2filt);set(gca,'Ydir','reverse')
        % legend('Contra El - R','Contra El - L', 'Ipsi El - R', 'Ipsi El - L')
        % xlabel('Time (ms)')
        % axis square
        % set(gca, 'FontSize', 16)
        % f1.Position = [10 10 420 420];
    
        % ERP
        f1 = figure;
        subplot(2, 2, 1)
        plot(EEG.times, nanmean([meanGAcontra1filt; meanGAcontra2filt]), 'Color', [199,21,133]/255, 'LineWidth', 2);set(gca,'Ydir','reverse')
        hold on
        plot(EEG.times, nanmean([meanGAipsi1filt; meanGAipsi2filt]), 'Color', 'black', 'LineWidth', 2);set(gca,'Ydir','reverse')
        ylim([-6, 2])
        xlim([-200, 1000])
        xlabel('Time (ms)')
        legend('Contralateral','Ipsilateral')
        xlabel('Time (ms)')
        ylabel('Amplitude (μV)')
        set(gca, 'FontSize', 16)
        set(gcf, 'Color', 'white')
        axis square
        text(-0.25, 1.15, 'A' ,'FontSize', 20, 'Units', 'normalized')
        subplot(2, 2, 2)
        plot(EEG.times, meanCDA2filt,'b', 'LineWidth', 2);set(gca,'Ydir','reverse')
        hold on;
        plot(EEG.times, meanCDA4filt,'g', 'LineWidth', 2);set(gca,'Ydir','reverse')
        plot(EEG.times, meanCDA6filt,'r', 'LineWidth', 2);set(gca,'Ydir','reverse')
        legend({'Set size 2','Set size 4', 'Set size 6'}, 'NumColumns', 3, 'Position', ...
            [0.522916666666667,0.935732647814908,0.371428571428571,0.027634961439589])
        % ylim([-1.2, 1])
        xlim([-200, 1000])
        xlabel('Time (ms)')
        ylabel('Amplitude (μV)')
        set(gca, 'FontSize', 16)
        set(gcf, 'Color', 'white')
        axis square
        text(-0.25, 1.15, 'B' ,'FontSize', 20, 'Units', 'normalized')
    
        % Figure 3 - mean amplitude and visual memory capacity
        subplot(2, 2, 3)
        p = errorbar([2, 4, 6], [meanSUB2, meanSUB4, meanSUB6], [ciCDA2 ciCDA4 ciCDA6], '-dblack'); set(gca,'Ydir','reverse')
        xlabel('Memory array size'); ylabel('Mean amplitude (μV)')
        xlim([1, 7])
        p.MarkerSize = 8;
        p.MarkerFaceColor = [0 0 0];
        set(gca, 'FontSize', 16)
        axis square
        text(-0.25, 1.15, 'C' ,'FontSize', 20, 'Units', 'normalized')
        % Figure 3b
        subplot(2, 2, 4)
        p = scatter(memoryCapacity, ampIncrease_2_to_4, 50, 'dblack');set(gca,'Ydir','reverse')
        lsline
        xlabel('Memory capacity'); ylabel('Amplitude increase from two to four items')
        xlim([0.5, 4])
        p.MarkerFaceColor = [0 0 0];
        set(gca, 'FontSize', 16)
        set(gcf, 'Color', 'white')
        axis square
        box on
        text(-0.25, 1.15, 'D' ,'FontSize', 20, 'Units', 'normalized')
        f1.Position = [10 10 420*2 420*2];
    
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Statistics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1) outcome neutral test, contra vs ipsilateral sites
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % compute t-test and effect sizes
        % script adapted from: https://doi.org/10.1016/j.cortex.2025.05.014
        [mean_amps, ~, ~, stats] = custom_paired_t_test(contra, ipsi, 0.05);
        
        % save stats
        all_stats.num_subs(lab) = length(contra);
        all_stats.on_m_contra(lab) = nanmean(contra);
        all_stats.on_sd_contra(lab) = nanstd(contra);
        all_stats.on_m_ipsi(lab) = nanmean(ipsi);
        all_stats.on_sd_ipsi(lab) = nanstd(ipsi);
        all_stats.on_t_stat(lab) = stats.t;
        all_stats.on_df(lab) = stats.df;
        all_stats.on_sd(lab) = stats.diff_std;
        all_stats.on_p(lab) = stats.p;
        all_stats.on_ci1(lab) = stats.mean_diff - stats.diff_ci;
        all_stats.on_ci2(lab) = stats.mean_diff + stats.diff_ci;
        all_stats.on_cohensd(lab) = stats.dz.eff;
        all_stats.on_cohensd_se(lab) = stats.dz.se;
        all_stats.on_gz(lab) = stats.gz.eff;
        all_stats.on_gz_se(lab) = stats.gz.se;
    
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2) compare amplitude of correct and incorrect trials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        allSz = [CDA2 + CDA4 + CDA6];
        allSz = allSz / 3;
        m_corr = mean(allSz(:, TOI), 2);
        m_incorr = mean(CDA_INCORR(:, TOI), 2);
    
        % compute t-test and effect sizes
        [mean_amps, ~, ~, stats] = custom_paired_t_test(m_corr, m_incorr, 0.05);
        
        % save stats    
        all_stats.corrincorr_m_corr(lab) = nanmean(m_corr);
        all_stats.corrincorr_sd_corr(lab) = nanstd(m_corr);
        all_stats.corrincorr_m_incorr(lab) = nanmean(m_incorr);
        all_stats.corrincorr_sd_incorr(lab) = nanstd(m_incorr);
        all_stats.corrincorr_t_stat(lab) = stats.t;
        all_stats.corrincorr_df(lab) = stats.df;
        all_stats.corrincorr_sd(lab) = stats.diff_std;
        all_stats.corrincorr_p(lab) = p;
        all_stats.corrincorr_ci1(lab) = stats.mean_diff - stats.diff_ci;
        all_stats.corrincorr_ci2(lab) = stats.mean_diff + stats.diff_ci;
        all_stats.corrincorr_cohensd(lab) = stats.dz.eff;
        all_stats.corrincorr_cohensd_se(lab) = stats.dz.se;
        all_stats.corrincorr_gz(lab) = stats.gz.eff;
        all_stats.corrincorr_gz_se(lab) = stats.gz.se;
    
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3) pairwaise comparisions of setsizes: 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % repeated measures anova
        subjects = categorical(1:length(meanROI2))';
        setsize = table([2, 4, 6]', 'VariableNames', {'Load'});
        tbl = table(subjects, meanROI2, meanROI4, meanROI6, 'VariableNames', {'Subjects', 'SZ2', 'SZ4', 'SZ6'});
        rm = fitrm(tbl, 'SZ2-SZ6 ~ 1', 'WithinDesign', setsize);
        ranovaResults = ranova(rm);
        disp(ranovaResults);
        
        % save stats    
        all_stats.setsize_anova{lab} = ranovaResults;
        all_stats.setsize2_m(lab) = nanmean(meanROI2);
        all_stats.setsize2_sd(lab) = nanstd(meanROI2);
        all_stats.setsize4_m(lab) = nanmean(meanROI4);
        all_stats.setsize4_sd(lab) = nanstd(meanROI4);
        all_stats.setsize6_m(lab) = nanmean(meanROI6);
        all_stats.setsize6_sd(lab) = nanstd(meanROI6);
    
        % correlations
        all_stats.corr_2_4(lab) = corr(meanROI2, meanROI4);
        all_stats.corr_2_6(lab) = corr(meanROI2, meanROI6);
        all_stats.corr_4_6(lab) = corr(meanROI4, meanROI6);
    
        % post-hoc t-tests: 2 vs 4
        % compute t-test and effect sizes
        [mean_amps, ~, ~, stats] = custom_paired_t_test(meanROI2, meanROI4, 0.05);
        % save stats
        all_stats.ph_2_4_t_stat(lab, 1) = stats.t;
        all_stats.ph_2_4_df(lab) = stats.df;
        all_stats.ph_2_4_sd(lab) = stats.diff_std;
        all_stats.ph_2_4_p(lab) = p;
        all_stats.ph_2_4_ci1(lab) = stats.mean_diff - stats.diff_ci;
        all_stats.ph_2_4_ci2(lab) = stats.mean_diff + stats.diff_ci;
        all_stats.ph_2_4_cohensd(lab) = stats.dz.eff;
        all_stats.ph_2_4_cohensd_se(lab) = stats.dz.se;
        all_stats.ph_2_4_gz_se(lab) = stats.gz.eff;
        all_stats.ph_2_4_gz_se(lab) = stats.gz.se;

        
        % post-hoc t-tests: 2 vs 6
        % compute t-test and effect sizes
        [mean_amps, ~, ~, stats] = custom_paired_t_test(meanROI2, meanROI6, 0.05);
        % save stats
        all_stats.ph_2_6_t_stat(lab, 1) = stats.t;
        all_stats.ph_2_6_df(lab) = stats.df;
        all_stats.ph_2_6_sd(lab) = stats.diff_std;
        all_stats.ph_2_6_p(lab) = p;
        all_stats.ph_2_6_ci1(lab) = stats.mean_diff - stats.diff_ci;
        all_stats.ph_2_6_ci2(lab) = stats.mean_diff + stats.diff_ci;
        all_stats.ph_2_6_cohensd(lab) = stats.dz.eff;
        all_stats.ph_2_6_cohensd_se(lab) = stats.dz.se;
        all_stats.ph_2_6_gz_se(lab) = stats.gz.eff;
        all_stats.ph_2_6_gz_se(lab) = stats.gz.se;
        
        % post-hoc t-tests
        % compute t-test and effect sizes
        [mean_amps, ~, ~, stats] = custom_paired_t_test(meanROI4, meanROI6, 0.05);
        % save stats
        all_stats.ph_4_6_t_stat(lab, 1) = stats.t;
        all_stats.ph_4_6_df(lab) = stats.df;
        all_stats.ph_4_6_sd(lab) = stats.diff_std;
        all_stats.ph_4_6_p(lab) = p;
        all_stats.ph_4_6_ci1(lab) = stats.mean_diff - stats.diff_ci;
        all_stats.ph_4_6_ci2(lab) = stats.mean_diff + stats.diff_ci;
        all_stats.ph_4_6_cohensd(lab) = stats.dz.eff;
        all_stats.ph_4_6_cohensd_se(lab) = stats.dz.se;
        all_stats.ph_4_6_gz_se(lab) = stats.gz.eff;
        all_stats.ph_4_6_gz_se(lab) = stats.gz.se;
        
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % correlation between amplitude increase (set size 2 to 4) and memory capacity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [rho, pval] = corr(ampIncrease_2_to_4, memoryCapacity);
        % compute SE for correlations
        % https://online.ucpress.edu/collabra/article/9/1/87615/197169/A-Brief-Note-on-the-Standard-Error-of-the-Pearson
        rho_se = (1-rho^2) / sqrt(length(ampIncrease_2_to_4) - 3);
        % save stats
        all_stats.wm_corr_2_4_amp_m(lab) = nanmean(ampIncrease_2_to_4);
        all_stats.wm_corr_2_4_amp_sd(lab) = nanstd(ampIncrease_2_to_4);
        all_stats.wm_corr_2_4_wm_m(lab) = nanmean(memoryCapacity);
        all_stats.wm_corr_2_4_wm_sd(lab) = nanstd(memoryCapacity);
        all_stats.wm_corr_2_4_r(lab) = rho;
        all_stats.wm_corr_2_4_p(lab) = pval;
        all_stats.wm_corr_2_4_r_se(lab) = rho_se;
        
        % correlation, where K = 2 instead of average
        [rho, pval] = corr(ampIncrease_2_to_4, lab_tab.K2);
        % compute SE for correlations
        rho_se = (1-rho^2) / sqrt(length(ampIncrease_2_to_4) - 3);
        % save stats
        all_stats.wmk2_corr_2_4_amp_m(lab) = nanmean(ampIncrease_2_to_4);
        all_stats.wmk2_corr_2_4_amp_sd(lab) = nanstd(ampIncrease_2_to_4);
        all_stats.wmk2_corr_2_4_wm_m(lab) = nanmean(lab_tab.K2);
        all_stats.wmk2_corr_2_4_wm_sd(lab) = nanstd(lab_tab.K2);
        all_stats.wmk2_corr_2_4_r(lab) = rho;
        all_stats.wmk2_corr_2_4_p(lab) = pval;
        all_stats.wmk2_corr_2_4_r_se(lab) = rho_se;
        
        % correlation, where K = 4
        [rho, pval] = corr(ampIncrease_2_to_4, lab_tab.K4);
        rho_se = (1-rho^2) / sqrt(length(ampIncrease_2_to_4) - 3);
        % save stats
        all_stats.wmk4_corr_2_4_amp_m(lab) = nanmean(ampIncrease_2_to_4);
        all_stats.wmk4_corr_2_4_amp_sd(lab) = nanstd(ampIncrease_2_to_4);
        all_stats.wmk4_corr_2_4_wm_m(lab) = nanmean(lab_tab.K4);
        all_stats.wmk4_corr_2_4_wm_sd(lab) = nanstd(lab_tab.K4);
        all_stats.wmk4_corr_2_4_r(lab) = rho;
        all_stats.wmk4_corr_2_4_p(lab) = pval;
        all_stats.wmk4_corr_2_4_r_se(lab) = rho_se;
        
        % correlation, where K = 6
        [rho, pval] = corr(ampIncrease_2_to_4, lab_tab.K6);
        rho_se = (1-rho^2) / sqrt(length(ampIncrease_2_to_4) - 3);   
        % save stats
        all_stats.wmk6_corr_2_4_amp_m(lab) = nanmean(ampIncrease_2_to_4);
        all_stats.wmk6_corr_2_4_amp_sd(lab) = nanstd(ampIncrease_2_to_4);
        all_stats.wmk6_corr_2_4_wm_m(lab) = nanmean(lab_tab.K6);
        all_stats.wmk6_corr_2_4_wm_sd(lab) = nanstd(lab_tab.K6);
        all_stats.wmk6_corr_2_4_r(lab) = rho;
        all_stats.wmk6_corr_2_4_p(lab) = pval;
        all_stats.wmk6_corr_2_4_r_se(lab) = rho_se;
        
        % correlation between amplitude increase (set size 4 to 6) and memory capacity
        ampIncrease_4_to_6 = meanROI6 - meanROI4;
        [rho, pval] = corr(ampIncrease_4_to_6, memoryCapacity);
        rho_se = (1-rho^2) / sqrt(length(ampIncrease_4_to_6) - 3);
        
        % save stats
        all_stats.wm_corr_4_6_amp_m(lab) = nanmean(ampIncrease_4_to_6);
        all_stats.wm_corr_4_6_amp_sd(lab) = nanstd(ampIncrease_4_to_6);
        all_stats.wm_corr_4_6_wm_m(lab) = nanmean(memoryCapacity);
        all_stats.wm_corr_4_6_wm_sd(lab) = nanstd(memoryCapacity);
        all_stats.wm_corr_4_6_r(lab) = rho;
        all_stats.wm_corr_4_6_p(lab) = pval;
        all_stats.wm_corr_4_6_r_se(lab) = rho_se;
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % saving
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    save(fullfile(datapath, 'mat_files', strcat('all_stats_', pipes_labels{pip}, '.mat')), 'all_stats')
    
    % remove the ANOVA tab before saving as csv
    all_stats_csv = removevars(all_stats, 'setsize_anova');
    writetable(all_stats_csv, fullfile(datapath, 'csv_files', strcat('all_stats_', pipes_labels{pip}, '.csv')))

end

