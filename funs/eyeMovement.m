%% Compute thresholds for saccade of 1deg


% 101 - left 6 degree saccade
% 102 - left 3 degree saccade
% 103 - right 3 degree saccade
% 104 - right 6 degree saccade

function eye_cal = eyeMovement(EYE, current_lab)

    % Extract ET data
    ET_chan = {'L-GAZE-X', 'L-GAZE-Y'};
    if any(ismember({EYE.chanlocs.labels}, ET_chan))
        % remove ET from EEG data
        EYE = pop_select(EYE, 'nochannel', [ET_chan]);
    end
    
    % bandpass filter 0.1 - 40 Hz
    EYE = pop_eegfiltnew(EYE, 0.1, 40);
    
    % rereference to mastoid
    el_m1 = find(strcmp({EYE.chanlocs.labels}, 'M1')); % M1
    el_m2 = find(strcmp({EYE.chanlocs.labels}, 'M2')); % M2
    % some labs dont have M1/M2 - use T7/T8
    if isempty(el_m1) | isempty(el_m2)
        el_m1 = find(strcmp({EYE.chanlocs.labels}, 'T7')); % M1
        el_m2 = find(strcmp({EYE.chanlocs.labels}, 'T8')); % M2
    end
    EYE = pop_reref(EYE, [el_m1 el_m2], 'keepref', 'on');
    
    % find HEOG and compute right - left
    if strcmp(current_lab, 'Dartmouth') | strcmp(current_lab, 'Zürich (Prof. Sauseng)')
        i = find(ismember({EYE.chanlocs.labels}, {'HEOG'}));
        EYE = pop_select(EYE, 'channel', [i]);
    else
        i = find(ismember({EYE.chanlocs.labels}, {'HEOGR', 'HEOGL'}));
        EYE = pop_select(EYE, 'channel', [i]);
        i_r = find(ismember({EYE.chanlocs.labels}, {'HEOGR'}));
        i_l = find(ismember({EYE.chanlocs.labels}, {'HEOGL'}));
        EYE.data = EYE.data(i_r, :) - EYE.data(i_l, :);
    end
    
    % %%%%%%%%%%% sanity checks: plot HEOG data
    % ttt = 1:78000;
    % lats = cell2mat({EYE.event.latency});
    % evs = {EYE.event.type};
    % f1 = figure;
    % hold on
    % plot(EYE.times(ttt), EYE.data(1, ttt))
    % hold on
    % f1.Position = [10, 10, 1600, 410];
    % legend('HEOG')
    % %%%%%%%%%%% end sanity checks
            
    
    EYE = pop_epoch(EYE, {'101', '102', '103', '104'}, [-0.2, 0.6]);
    EYE = pop_rmbase(EYE, [EYE.times(1), 0]);
    
    %% bad trial detection
    % remove trials that have large pre-stimulus activity
    i_bad = squeeze(sum(EYE.data(1, 1 : 60, :) > 50 | EYE.data(1, 1 : 60, :) < -50)) > 1 | squeeze(sum(EYE.data(1, :, :) < -300 | EYE.data(1, :, :) > 300)) > 1;
    
    EYE = pop_select(EYE, 'notrial', [find(i_bad)]);
    
    %% plotting
    CHAN = 1;
    
    epochs = {EYE.event.epoch};
    types = {EYE.event.type};
    % some events are numbers, other strings
    if isnumeric(cell2mat(types))
        types = cellfun(@num2str, types, 'UniformOutput', false);
    end
    i = find(ismember(types, '101'));
    iL6 = cell2mat(epochs(i));
    L6 = squeeze(mean(EYE.data(CHAN, :, iL6), 3));
    
    i = find(ismember(types, '102'));
    iL3 = cell2mat(epochs(i));
    L3 = squeeze(mean(EYE.data(CHAN, :, iL3), 3));
    
    i = find(ismember(types, '103'));
    iR3 = cell2mat(epochs(i));
    R3 = squeeze(mean(EYE.data(CHAN, :, iR3), 3));
    
    i = find(ismember(types, '104'));
    iR6 = cell2mat(epochs(i));
    R6 = squeeze(mean(EYE.data(CHAN, :, iR6), 3));
    
    
    % figure;
    % plot(EYE.times, squeeze(mean(EYE.data(CHAN, :, :), 1)))
    % xlabel('Time (ms)' )
    % ylabel('HEOG Voltage' )
    % 
    % figure;
    % plot(EYE.times, L6)
    % hold on
    % plot(EYE.times, L3)
    % plot(EYE.times, R3)
    % plot(EYE.times, R6)
    % xlabel('Time (ms)' )
    % ylabel('HEOG Voltage' )
    % legend('Left, 6 deg', 'Left, 3 deg', 'Right, 3 deg', 'Right, 6 deg')
    
    
    
    %% compute average amplitude between 300 and 400 ms
    TOI = EYE.times >= 300 & EYE.times <= 400;
    
    avgL6 = squeeze(mean(EYE.data(CHAN, TOI, iL6), 2));
    avgL3 = squeeze(mean(EYE.data(CHAN, TOI, iL3), 2));
    avgR3 = squeeze(mean(EYE.data(CHAN, TOI, iR3), 2));
    avgR6 = squeeze(mean(EYE.data(CHAN, TOI, iR6), 2));
    
    
    % avgL6(avgL6<0) = [];
    % avgL3(avgL3<0) = [];
    % 
    % avgR6(avgR6>0) = [];
    % avgR3(avgR3>0) = [];
    
    %% Approach 1: estimate the amplitude for 1 degree
    % avgL1 = (avgL6 - avgL3 ) / 3;
    % avgR1 = (avgR6 - avgR3 ) / 3;
    % 
    % figure;
    % errorbar([1, 3, 6], [mean(avgL1), mean(avgL3), mean(avgL6)], [std(avgL1), std(avgL3), std(avgL6)], 'Color', 'red')
    % hold on
    % errorbar([1, 3, 6], [mean(avgR1), mean(avgR3), mean(avgR6)], [std(avgR3), std(avgR3), std(avgR6)], 'Color', 'blue')
    % xlim([0, 7])
    % legend('Left saccades', 'Right saccades')
    % xlabel('Degrees')
    % ylabel('HEOG amplitude')
    % set(gca, 'FontSize', 16)
    
    
    %% Approach 2: use fit function and force the regression line to go trough 0
    
    % left
    ft1 = fittype({'x'}); %This creates a linear 'fittype' variable that is of the form f(a,x)=ax.
    x1 = [3*ones(1, length(avgL3)), 6*ones(1,length(avgL6))];
    y1 = [avgL3; avgL6];
    p1 = fit(x1',y1,ft1); %This creates a 'cfit' variable p that is your fitted function
    x_fit = linspace(0,10,11); %x-values to evaluate
    ci1 = predint(p1,x_fit,0.95,'obs');
    y1_fitted = feval(p1, x_fit); %y-values for the evaluated x-values
    avgL1 = y1_fitted(x_fit == 1); % find the value for 1 deg saccade
    
    % right
    ft2 = fittype({'x'}); %This creates a linear 'fittype' variable that is of the form f(a,x)=ax.
    x2 = [3*ones(1, length(avgR3)), 6*ones(1,length(avgR6))];
    y2 = [avgR3; avgR6];
    p2 = fit(x2',y2,ft2); %This creates a 'cfit' variable p that is your fitted function
    x_fit = linspace(0,10,11); %x-values to evaluate
    ci2 = predint(p2,x_fit,0.95,'obs');
    y2_fitted = feval(p2, x_fit); %y-values for the evaluated x-values
    avgR1 = y2_fitted(x_fit == 1); % find the value for 1 deg saccade
    
    % % 
    % figure;
    % hold on;
    % p1 = plot(x1,y1,'go'); hold on;
    % plot(x2,y2,'go')
    % p2 = plot(1, avgL1, 'ro');
    % p3 = plot(1, avgR1, 'bo');
    % p4 = plot(x_fit,y1_fitted,'r--'); plot(x_fit, ci1,'r.-'); % with 95% CI
    % p5 = plot(x_fit,y2_fitted,'b--'); plot(x_fit, ci2,'b.-'); 
    % xlabel('Degree' )
    % ylabel('HEOG Voltage' )
    % legend([p1, p2, p3, p4, p5], {'Raw data', 'Estimated', 'Estimated', 'Left saccades', 'Right saccades'});
    % 
    %% quantify uncertainty and add to the estimate
    eye_cal = struct;
    eye_cal.avgR1 = avgR1;
    eye_cal.avgL1 = avgL1;
    eye_cal.avgR3 = avgR3;
    eye_cal.avgL3 = avgL3;
    eye_cal.avgR6 = avgR6;
    eye_cal.avgL6 = avgL6;
    
    
    % EEG system from Dartmouth provides only 1 HEOG channel, which is already
    % pre-computed as HEOGL - HEOGR
    if strcmp(current_lab, 'Dartmouth') | strcmp(current_lab, 'Zürich (Prof. Sauseng)')
        
        eye_cal.avgR1_95CI = ci2(2, 1)*(-1);
        eye_cal.avgL1_95CI = ci1(2, 2)*(-1);
        
        eye_cal.avgR1_plus_sd = (eye_cal.avgR1 - std(abs([eye_cal.avgL3; eye_cal.avgR3])))*(-1);
        eye_cal.avgL1_minus_sd = (eye_cal.avgL1 + std(abs([eye_cal.avgL3; eye_cal.avgR3])))*(-1);
   
    else
        eye_cal.avgR1_95CI = ci2(2, 2);
        eye_cal.avgL1_95CI = ci1(2, 1);
        
        eye_cal.avgR1_plus_sd = eye_cal.avgR1 + std(abs([eye_cal.avgL3; eye_cal.avgR3])); % make it symmetric
        eye_cal.avgL1_minus_sd = eye_cal.avgL1 - std(abs([eye_cal.avgL3; eye_cal.avgR3])); % make it symmetric
    end
    
    %% combine all 2 figures for the manuscript
    
    % f1 = figure;
    % subplot(1, 2, 1)
    % plot(EYE.times, squeeze(mean(EYE.data(CHAN, :, :), 1)), 'Color', 'black')
    % ylim([-250, 250])
    % xlabel('Time (ms)' )
    % ylabel('HEOG Voltage (μV)' )
    % set(gca, 'FontSize', 16)
    % set(gcf, 'Color', 'white')
    % text(min(xlim) - 120, max(ylim) +30, 'A' ,'FontSize', 16)
    % 
    % hold on
    % p1 = plot(EYE.times, L6, 'LineWidth', 2, 'Color',[0.8500 0.3250 0.0980]);
    % p2 = plot(EYE.times, L3, 'LineWidth', 2, 'Color',[0.9290 0.540 0.1250]);
    % p3 = plot(EYE.times, R3, 'LineWidth', 2, 'Color', [0.3 0.3250 0.9]);
    % p4 = plot(EYE.times, R6, 'LineWidth', 2, 'Color', [0.3010 0.7450 0.9330]);
    % ylim([-250, 250])
    % xlabel('Time (ms)' )
    % ylabel('HEOG Voltage (μV)' )
    % legend([p1, p2, p3, p4], 'Left saccades, 6 deg', 'Left saccades, 3 deg', 'Right saccades, 3 deg', 'Right saccades, 6 deg')
    % set(gca, 'FontSize', 16)
    % set(gcf, 'Color', 'white')
    % 
    % subplot(1, 2, 2)
    % p1 = plot(x1,y1,'go'); hold on;
    % plot(x2,y2,'go')
    % p2 = plot(1, avgL1, 'ro');
    % p3 = plot(1, avgR1, 'bo');
    % p4 = plot(x_fit,y1_fitted,'r--'); plot(x_fit, ci1,'r.-'); % with 95% CI
    % p5 = plot(x_fit,y2_fitted,'b--'); plot(x_fit, ci2,'b.-'); 
    % xlabel('Degree of visual angle' )
    % ylabel('HEOG Voltage (μV)' )
    % ylim([-250, 250])
    % set(gca, 'FontSize', 16)
    % set(gcf, 'Color', 'white')
    % text(min(xlim) - 1.5, max(ylim) +30, 'B' ,'FontSize', 16)
    % legend([p1, p2, p3, p4, p5], {'Raw data', 'Estimated', 'Estimated', 'Left saccades', 'Right saccades'});
    % 
    % f1.Position = [10 10 560*2 420];

end