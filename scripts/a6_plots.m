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

for pip = 3 : 3% length(pipelines)

    % load data
    load(fullfile(datapath, 'mat_files', strcat(tabs{pip}, '.mat')))

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

    % extract vars of interest
    GAcontra1 = cell2mat(clean_tab.GAcontra1);
    GAcontra2 = cell2mat(clean_tab.GAcontra2);
    GAipsi1 = cell2mat(clean_tab.GAipsi1);
    GAipsi2 = cell2mat(clean_tab.GAipsi2);
    CDA2 = cell2mat(clean_tab.CDA2);
    CDA4 = cell2mat(clean_tab.CDA4);
    CDA6 = cell2mat(clean_tab.CDA6);
    
    % compute means
    meanGAcontra1 = nanmean(GAcontra1, 1);
    meanGAcontra2 = nanmean(GAcontra2, 1);
    meanGAipsi1 = nanmean(GAipsi1, 1);
    meanGAipsi2 = nanmean(GAipsi2, 1);

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
    memoryCapacity = (clean_tab.K2 + clean_tab.K4 + clean_tab.K6) / 3; % average over all setsizes
    % 95% confidence interval
    ciCDA2 = 1.96 * (nanstd(meanROI2) / sqrt(size(meanROI2, 1))); % 95% CI
    ciCDA4 = 1.96 * (nanstd(meanROI4) / sqrt(size(meanROI4, 1))); % 95% CI
    ciCDA6 = 1.96 * (nanstd(meanROI6) / sqrt(size(meanROI6, 1))); % 95% CI


    % filter used only for visualization
    filt = 35;
    meanGAcontra1filt=eegfilt(meanGAcontra1,250,0,filt);
    meanGAcontra2filt=eegfilt(meanGAcontra2,250,0,filt);
    meanGAipsi1filt=eegfilt(meanGAipsi1,250,0,filt);
    meanGAipsi2filt=eegfilt(meanGAipsi2,250,0,filt);
    meanCDA2filt=eegfilt(meanCDA2,250,0,filt);
    meanCDA4filt=eegfilt(meanCDA4,250,0,filt);
    meanCDA6filt=eegfilt(meanCDA6,250,0,filt);


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

end


    