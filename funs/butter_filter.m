function EEG = butter_filter(EEG, low_cutoff, high_cutoff)

%
%   Applies a 2nd-order zero-phase Butterworth bandpass filter to EEG.data.
%   Defaults to 0.01–80 Hz bandpass if no cutoff values are provided.
%
%   Input:
%       EEG         - EEGLAB EEG structure
%       low_cutoff  - low cutoff in Hz 
%       high_cutoff - high cutoff in Hz
%
%   Output:
%       EEG - Filtered EEG structure

    order = 2;  % Filter order 
    nyquist = EEG.srate / 2;

    % Normalize frequencies
    low = low_cutoff / nyquist;
    high = high_cutoff / nyquist;

    % Design filter
    [b, a] = butter(order, [low high]);

    % Convert data to double for filtering
    EEG.data = double(EEG.data);

    % Apply zero-phase filtering to each channel
    for i = 1:size(EEG.data, 1)
        EEG.data(i, :) = filtfilt(b, a, EEG.data(i, :));
    end

end
