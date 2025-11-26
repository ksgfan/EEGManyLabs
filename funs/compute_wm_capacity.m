function [K, D_prime] = compute_wm_capacity(EEG, id, current_lab)

    % K - WM capacity
    % S - size of the array
    % H - observed hit rate
    % F - false alarm rate
    
    setSizes = [2, 4, 6];
    
    for setSize = 1 : 3
        
        i_setSize = (EEG.beh.data.trialSetSize == setSizes(setSize));
        
        if strcmp(current_lab, 'Zürich (Prof. Langer)') | ... 
           strcmp(current_lab, 'NorthDakota') | ...
           strcmp(current_lab, 'Münster') 
            KeyCodeA = 39; % on linux
            KeyCodeL = 47; % on linux
    
        elseif strcmp(current_lab, 'Mainz') | strcmp(current_lab, 'Sheffield') | ...
               strcmp(current_lab, 'Ohio') | strcmp(current_lab, 'Reykjavik')| ...
               strcmp(current_lab, 'Dartmouth') | strcmp(current_lab, 'Zürich (Prof. Sauseng)') | ...
               strcmp(current_lab, 'Florida')
            KeyCodeA = 65;
            KeyCodeL = 76; 
        end
    
        % Assign response keys
        if mod(id,2) == 0   % Use subject ID for assignment to ensure counterbalancing
            changeIsL = true;       % L is change, A is no change.
            hits = find(EEG.beh.data.trialIfChange == 1 & EEG.beh.data.allResponses == KeyCodeL & i_setSize); % trial change, response change
            misses = find(EEG.beh.data.trialIfChange == 1 & EEG.beh.data.allResponses == KeyCodeA & i_setSize); % trial change, response no change
            false_alarm = find(EEG.beh.data.trialIfChange == 0 & EEG.beh.data.allResponses == KeyCodeL & i_setSize); % trial no change, response change
            corr_rejection = find(EEG.beh.data.trialIfChange == 0 & EEG.beh.data.allResponses == KeyCodeA & i_setSize ); % trial no change, response no change
        elseif mod(id,2) == 1
            changeIsL = false;      % L is no change, A is change.
            hits = find(EEG.beh.data.trialIfChange == 1 & EEG.beh.data.allResponses == KeyCodeA & i_setSize); % trial change, response change
            misses = find(EEG.beh.data.trialIfChange == 1 & EEG.beh.data.allResponses == KeyCodeL & i_setSize); % trial change, response no change
            false_alarm = find(EEG.beh.data.trialIfChange == 0 & EEG.beh.data.allResponses == KeyCodeA & i_setSize); % trial no change, response change
            corr_rejection = find(EEG.beh.data.trialIfChange == 0 & EEG.beh.data.allResponses == KeyCodeL & i_setSize ); % trial no change, response no change
        end
        
        % Hit rate
        H = numel(hits) / (numel(hits) + numel(misses));
        
        % False alarm rate
        F = numel(false_alarm) / (numel(false_alarm) + numel(corr_rejection)); 
    
        % WM capacity
        K(setSize) = setSizes(setSize) * (H - F);
        
        % d'
        % first, check if H = 1 or F = 0
        if H == 1
            H = 1-1/(2*numel(EEG.beh.data.allResponses));
        end
        if F == 0
            F = 1/(2*numel(EEG.beh.data.allResponses));
        end
            
        D_prime(setSize) = norminv(H) - norminv(F);
    end
end