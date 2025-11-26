function [EEG, EYE] = rename_channels(EEG, EYE, current_lab)

    % Rename channel labels to ensure consistent naming across labs
    if strcmp(current_lab, 'Dartmouth')
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'Fp1')).labels = 'VEOGU'; % U - Upper
    elseif strcmp(current_lab, 'Florida')
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E57')).labels  = 'M1';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E100')).labels = 'M2';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E25')).labels  = 'VEOGU'; % U - Upper
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E127')).labels = 'VEOGL'; % L - Lower
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E32')).labels  = 'HEOGL';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E125')).labels = 'HEOGR';
    
        EYE.chanlocs(ismember({EYE.chanlocs.labels}, 'E57')).labels  = 'M1';
        EYE.chanlocs(ismember({EYE.chanlocs.labels}, 'E100')).labels = 'M2';
        EYE.chanlocs(ismember({EYE.chanlocs.labels}, 'E32')).labels  = 'HEOGL';
        EYE.chanlocs(ismember({EYE.chanlocs.labels}, 'E125')).labels = 'HEOGR';
    
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E52')).labels  = 'P3';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E58')).labels = 'P7';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E70')).labels  = 'O1';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E92')).labels = 'P4';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E96')).labels  = 'P8';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'E83')).labels = 'O2';
    elseif strcmp(current_lab, 'Mainz')
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'VEOG_O')).labels = 'VEOGU'; % O - oben, U - unten (!!!)
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'VEOG_U')).labels = 'VEOGL'; 
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'HEOG_L')).labels = 'HEOGL';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'HEOG_R')).labels = 'HEOGR';
    
        EYE.chanlocs(ismember({EYE.chanlocs.labels}, 'HEOG_L')).labels = 'HEOGL';
        EYE.chanlocs(ismember({EYE.chanlocs.labels}, 'HEOG_R')).labels = 'HEOGR';
    elseif strcmp(current_lab, 'Münster')
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'Fp1')).labels = 'VEOGU'; 
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'IO1')).labels = 'VEOGL';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'Afp9')).labels = 'HEOGL';
        EEG.chanlocs(ismember({EEG.chanlocs.labels}, 'Afp10')).labels = 'HEOGR';
    
        EYE.chanlocs(ismember({EYE.chanlocs.labels}, 'Afp9')).labels = 'HEOGL';
        EYE.chanlocs(ismember({EYE.chanlocs.labels}, 'Afp10')).labels = 'HEOGR';
    elseif strcmp(current_lab, 'Ohio')
        EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'Fp1'))).labels = 'VEOGU'; 
        EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'FT9'))).labels = 'HEOGL';
        EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'FT10'))).labels = 'HEOGR';
        
        EYE.chanlocs(find(ismember({EYE.chanlocs.labels}, 'FT9'))).labels = 'HEOGL';
        EYE.chanlocs(find(ismember({EYE.chanlocs.labels}, 'FT10'))).labels = 'HEOGR';
    elseif strcmp(current_lab, 'Reykjavik')
        EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'Fp1'))).labels = 'VEOGU';
        EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'FT9'))).labels = 'HEOGL';
        EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'FT10'))).labels = 'HEOGR';
    
        EYE.chanlocs(find(ismember({EYE.chanlocs.labels}, 'FT9'))).labels = 'HEOGL';
        EYE.chanlocs(find(ismember({EYE.chanlocs.labels}, 'FT10'))).labels = 'HEOGR';
    elseif strcmp(current_lab, 'Zürich (Prof. Sauseng)')
        EEG.chanlocs(find(ismember({EEG.chanlocs.labels}, 'Fp1'))).labels = 'VEOGU';
    end

end