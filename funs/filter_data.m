function [EEG, fig1] = filter_data(EEG, current_lab, params)
    
    % filter
    high = params.filter.high;              
    low = params.filter.low;

    EEG = pop_eegfiltnew(EEG, high, low);

    % ZapLine Plus
    [clean, zaplineConfig, analyticsResults] = clean_data_with_zapline_plus(EEG.data,EEG.srate, 'plotResults', 0, 'saveSpectra', 1);
    if not(isempty(zaplineConfig.noisefreqs))
        EEG.data = clean;  
        EEG.zaplineConfig = zaplineConfig;
        EEG.analyticsResults = analyticsResults;
        % plot figure with results
        fig1 = figure('Visible', 'off');
        f = analyticsResults.frequencies;
        plot(f,mean(analyticsResults.rawSpectrumLog,2),'color','black','linewidth',1.5);
        hold on
        plot(f,mean(analyticsResults.cleanSpectrumLog,2),'color','red','linewidth',1.5);
        set(gca,'ygrid','on','xgrid','on');
        set(gca,'yminorgrid','on')
        set(gca,'fontsize',20)
        xlabel('frequency [Hz]');
        ylabel('Power [10*log10 \muV^2/Hz]');
        title('ZapLinePlus')
        legend('raw', 'clean')
        clear f
    else
        % If ZapLine Plus fails, compute ZapLine that removes 7 components
        if strcmp(current_lab, 'NorthDakota') | ...
            strcmp(current_lab, 'Ohio') | ...
            strcmp(current_lab, 'Florida') | ...
            strcmp(current_lab, 'Dartmouth')
    
            % USA
            zapline.freq = 60;
        else
            % Europe
            zapline.freq = 50;
        end
        zapline.ncomps = 7;
        fline = zapline.freq / EEG.srate; % line frequency normalised to srate
        clean = nt_zapline(EEG.data', fline , zapline.ncomps);
    
        % plot power spectrum of raw data
        segment_length = 8 * EEG.srate;     % 2-second windows
        overlap = segment_length / 2;   % 50% overlap
        nfft = max(256, 2^nextpow2(segment_length)); 
        all_pxx_raw = [];
        for ch = 1:EEG.nbchan
            [pxx, f] = pwelch(double(EEG.data(ch, :)), segment_length, overlap, nfft, EEG.srate);
            all_pxx_raw(:, ch) = pxx;  
        end
    
        % replace with clean data
        EEG.data = clean'; 
    
        all_pxx_clean = [];
        for ch = 1:EEG.nbchan
            [pxx, f] = pwelch(double(EEG.data(ch, :)), segment_length, overlap, nfft, EEG.srate);
            all_pxx_clean(:, ch) = pxx;  
        end
    
        % plot
        fig1 = figure('Visible', 'off');
        plot(f,10*log10(mean(all_pxx_raw,2)),'color','black','linewidth',1.5);
        hold on
        plot(f,10*log10(mean(all_pxx_clean,2)),'color','red','linewidth',1.5);
        set(gca,'ygrid','on','xgrid','on');
        set(gca,'yminorgrid','on')
        set(gca,'fontsize',20)
        xlabel('frequency [Hz]');
        ylabel('Power [10*log10 \muV^2/Hz]');
        title('ZapLine')
        legend('raw', 'clean')
        clear f pxx all_pxx_raw all_pxx_clean
    end
    
end