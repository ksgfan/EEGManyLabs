

%% description

% BrainVision
% reference: M2
% some subjects have electrode 'Fp1', others 'FP1'
% no eye tracker

% To explain more, M1 is left mastoid and M2 is right mastoid. 
% But when we record data we already subtract M2 activity from channels (including M1). 
% This is default setting from the system so M2 works as implicit reference, 
% and that's why we technically don't see data from M2. 
% Then we "re-reference" to the average of M1 and M2 in our analysis, 
% which mathematically is identical to subtracting half of M1 from the channels. 

% we only get 1 HEOG channel data from the system. I believe it's 
% Left - Right, and when I plot left eye movements are positive and right movements are negative.

% s01_cda_3.eeg (and vhdr, vmrk) was manually renamed to s1_cda_3.eeg

% subject '3' has only 705 EEG trials and 720 responses in beh files => recompute from events
% subject '6' has HEOG reversed (right - left). multiply the entire time series by (-1)
% Subject 12 eye and cda data are merged (they're same copies), but because the event codes are not overlapping it shouldn't cause issues in analysis.
% subject 14 has CDA and resting in the same file
% subject 18 has CDA and resting in the same file
% subject 23: 599 EEG trials, 576 responses  => recompute from events

% subject 11, 18, 21 - very high amplitudes
% subject 24 - P8 is very noisy, but the subject is not marked as 'bad' by exclusion criteria 
% subject 25 - P4 and P7 are very noisy, but the subject is not marked as 'bad' by exclusion criteria

%% merge EEG and ET data, and save as mat file
d = dir(strcat(fullfile(rawdatapath, 'data_dartmouth/')));
d = d(~contains({d.name}, {'.', '..', '.DS_Store', 'txt'}));

% run loop
result_folder = fullfile(formatted_data, 'dartmouth');

for sub = 1 : size(d, 1)
    
    id = d(sub).name;
    disp(id)
    
    d_sub = dir(fullfile(d(sub).folder, d(sub).name, '*eeg'));
    
    % subject 1 has 5 files (for each block). merge them
    if strcmp(id, '1')
        EEGconcat = {};
        d_01_cda = dir(fullfile(d(sub).folder, d(sub).name, '*cda*eeg'));
        for block = 1 : size(d_01_cda, 1)
            EEG = pop_fileio(fullfile(d_01_cda(block).folder, d_01_cda(block).name));
            EEGconcat{block} = EEG;
        end
        
        E = EEGconcat{1};
        for ss = 2 : length(EEGconcat)
            E = pop_mergeset(E, EEGconcat{ss},  0);
        end
        tmpEEG = E;
        
        % remove 4 of 5 blocks from d_sub
        d_sub(1:4, :) = [];
    end
                

    mkdir(fullfile(result_folder, id))
    for f = 1 : size(d_sub, 1)
        
        % special case subject 1
        if strcmp(id, '1') & contains(d_sub(f).name, 'cda')
            EEG = tmpEEG; % dont load but use already merged data from workspace
        else
            % load EEG data
            EEG = pop_fileio(fullfile(d_sub(f).folder, d_sub(f).name));
        end
        
        % remove diode
        EEG = pop_select(EEG, 'nochannel', [33]);

        % load chanloc file and add reference
        EEG.data(33, :) = 0;
        EEG.nbchan = 33;
        EEG.chanlocs(33).labels = 'M2';    
        EEG.chanlocs(33).type = 'eeg';  
        locspath = '32channelsGreenwithDiode.loc';
        EEG = pop_chanedit(EEG, 'lookup', locspath);
        
        % 
        EEG.chanlocs(32).type = 'eog'; 

        % % sanity check
        % figure;
        % topoplot([], EEG.chanlocs, 'electrodes', 'labels')

        % info for BIDS format
        EEG.tInfo.OldID = id;
        EEG.tInfo.EEGReference = 'M2';
        EEG.tInfo.InstitutionName = 'Dartmouth College';
        EEG.tInfo.PowerLineFrequency = 60;
        EEG.tInfo.EEGGround = 'Fpz';
        EEG.tInfo.CapManufacturer = 'BrainVision';
        EEG.tInfo.SoftwareFilters = "n/a";
        
        % rename FP1 to Fp1
        fp1 = find(ismember({EEG.chanlocs.labels}, 'FP1'));
        if not(isempty(fp1))
            EEG.chanlocs(fp1).labels = 'Fp1';
        end

        % Remove 'S ' from EEG.event.type
        A = {EEG.event.type};
        for i = 1:numel(A)
            if startsWith(EEG.event(i).type, 'S')
                EEG.event(i).type = strip(num2str(A{i}(2:end)));
            end
        end

        % check for 'empty' event
        if any(strcmp({EEG.event.type}, 'empty'))
            i_empty = strcmp({EEG.event.type}, 'empty');
            [EEG.event(i_empty).type] = deal('999');
            [EEG.event(i_empty).timestamp] = deal([]);
        end
        if any(strcmp({EEG.urevent.type}, 'empty'))
            i_empty = strcmp({EEG.urevent.type}, 'empty');
            [EEG.urevent(i_empty).type] = deal('999');
            [EEG.urevent(i_empty).timestamp] = deal([]);
        end

        % check if it's resting (triggers 20 and 30 are only in Resting EEG)
        if sum(ismember({EEG.event.type}, '20')) > 1
            task = [id, '_Resting_EEG.mat'];
        elseif sum(ismember({EEG.event.type}, '55')) > 2
            task = [id, '_Eye_EEG.mat'];
        else
            task = [id, '_CDA_EEG.mat'];
        end

        % special case subject 6 - reverse HEOG channel
        if strcmp(id, '6')
            i_heog = find(ismember({EEG.chanlocs.labels}, 'HEOG'));
            EEG.data(i_heog, :) = EEG.data(i_heog, :) * (-1);
        end
        
        % special case subject 12
        if strcmp(id, '12') & contains(d_sub(f).name, 'cda')
            task = [id, '_CDA_EEG.mat'];
        elseif strcmp(id, '12') & contains(d_sub(f).name, 'eye')
            task = [id, '_Eye_EEG.mat'];
        end
        
        % special case subject 14 (cda file has a resting state and cda)
        if (strcmp(id, '14') | strcmp(id, '18')) & contains(d_sub(f).name, 'cda')
            task = [id, '_CDA_EEG.mat'];
            EEGorig = EEG;
            i10 = find(ismember({EEGorig.event.type}, '10'));
            i90 = find(ismember({EEGorig.event.type}, '90'));
            EEG = pop_select(EEGorig, 'point', [EEGorig.event(i10).latency, EEGorig.event(i90).latency]);
            save(fullfile(result_folder, id, strcat(id,'_Resting_EEG.mat')), 'EEG');
            EEG = EEGorig;
        end
        
        % save Eye and Resting to a file or prepare to concat CDA blocks
        if not(contains(task, '_CDA_EEG'))
            save(fullfile(result_folder, id, task), 'EEG', '-v7.3')
        else
            % load and merge CDA behavioral data 
            BEHconcat = {};
            beh = [];
            trigger = [];
            d_beh = dir(fullfile(d_sub(f).folder, '*ManyLabs*block*task.mat'));
            for b = 1 : size(d_beh, 1)
                load(fullfile(d_beh(b).folder, d_beh(b).name)); % behavioral
                BEHconcat{b} = beh.data;  
            end
            BEH = [BEHconcat{:}];
            newBEH = struct;
            field_names = fieldnames(BEH);
            for n = 1 : length(field_names)
                tmp = [];
                field_size = size(BEH(1).(field_names{n}));
                for ss = 1 : length(BEHconcat)
                    if field_size(1) > field_size(2)
                        tmp = [tmp; BEH(ss).(field_names{n})]; % concat verticaly
                    elseif field_size(1) < field_size(2)
                        tmp = [tmp, BEH(ss).(field_names{n})]; % concat horizonztally
                    end
                end
                % overwrite the field with merged data
                newBEH.(field_names{n}) = tmp;
            end
            
            
            % subject '3' has less EEG trials than responses. double check
            if strcmp(id, '3') | strcmp(id, '23')
                i_events = find(ismember({EEG.event.type}, {'21', '41', '61'}));
                num_trials = length(i_events); % 705 instead of 720

                % recompute all events
                % first, we need longer segment than contains cues and responses
                tmpEEG = pop_epoch(EEG, {21, 41, 61}, [-0.7 3]);
                events = {tmpEEG.event.type};
                epochs = {tmpEEG.event.epoch};
    
                % init beh data
                allResponses = nan(tmpEEG.trials, 1); % 76, 65
                tmpResponses = nan(tmpEEG.trials, 1);
                allCorrect = nan(tmpEEG.trials, 1); % 0, 1
                trialSetSize = nan(1, tmpEEG.trials); % 2, 4, 6
                trialIfChange = nan(1, tmpEEG.trials); % 0, 1
                trialCuedSide = nan(1, tmpEEG.trials); % 0, 1
    
                % extract data 
                ev_onset = events(ismember(events, '21') | ismember(events, '41') |ismember(events, '61'));
                trialSetSize(ismember(ev_onset, '21')) = 2;
                trialSetSize(ismember(ev_onset, '41')) = 4;
                trialSetSize(ismember(ev_onset, '61')) = 6;
    
                ev_cue = events(ismember(events, '3') | ismember(events, '7'));
                trialCuedSide(ismember(ev_cue, '3')) = 0; % left
                trialCuedSide(ismember(ev_cue, '7')) = 1; % right
    
                ev_resp = events(ismember(events, '76') | ismember(events, '77') | ismember(events, '78') | ismember(events, '79'));
                epochs_of_resp = cell2mat(epochs(ismember(events, '76') | ismember(events, '77') | ismember(events, '78') | ismember(events, '79')));
                tmpResponses(epochs_of_resp) = str2double(ev_resp);
                allCorrect(tmpResponses == 76 | tmpResponses == 77) = 1;
                allCorrect(tmpResponses == 78 | tmpResponses == 79) = 0;
                
                trialIfChange(tmpResponses == 76 | tmpResponses == 79) = 0; % no change
                trialIfChange(tmpResponses == 77 | tmpResponses == 78) = 1; % change
    
                if mod(str2double(id),2) == 0 
                    % L is change, A is no change.
                    allResponses(tmpResponses == 77 | tmpResponses == 79) = 76; % double checked with existing data
                    allResponses(tmpResponses == 76 | tmpResponses == 78) = 65;
                elseif mod(str2double(id),2) == 1
                    % L is no change, A is change.
                    allResponses(tmpResponses == 77 | tmpResponses == 79) = 65;
                    allResponses(tmpResponses == 76 | tmpResponses == 78) = 76;
                end

                % sanity check: old vs new beh data
                san_check = [newBEH.allResponses(1:10), allResponses(1:10)]
                
                % which block is it? => block 4 has 129 instead of 144
                % trials
                i14 = find(ismember({EEG.event.type}, '14'));
                i94 = find(ismember({EEG.event.type}, '94'));
                i_events = find(ismember({EEG.event(i14:i94).type}, {'21', '41', '61'}));
                num_trials = length(i_events) % 129 instead of 144
                

                % now replace BEH data 
                newBEH.allResponses = allResponses;
                newBEH.allCorrect = allCorrect;
                newBEH.trialSetSize = trialSetSize;
                newBEH.trialIfChange = trialIfChange;
                newBEH.trialCuedSide = trialCuedSide;

                
            end

            % append to merged EEG dataset
            EEG.beh = beh; % beh contains also experiment info
            EEG.beh.data = newBEH; % replace with merged data
            EEG.beh.triggers = trigger;

            % save CDA file
            save(fullfile(result_folder, id, task), 'EEG', '-v7.3')
        end
    end
end





