%% Convert data to BIDS format

p = pwd;
funpath = strsplit(p, filesep);
addpath(fullfile(strjoin(funpath(1:end-1), filesep), 'funs'))
initPaths;

%% overview

% data format: 
% - sub-Vog001
%   - eeg
%     - sub-Vog001_task-CDA_channels.tsv
%     - sub-Vog001_task-CDA_eeg.eeg
%     - sub-Vog001_task-CDA_eeg.json
%     - sub-Vog001_task-CDA_eeg.vhdr
%     - sub-Vog001_task-CDA_eeg.vmrk
%     - sub-Vog001_task-CDA_events.json
%     - sub-Vog001_task-CDA_events.tsv
%     - sub-Vog001_task-Eye_channels.tsv
%     - sub-Vog001_task-Eye_eeg.eeg
%     - sub-Vog001_task-Eye_eeg.json
%     - sub-Vog001_task-Eye_eeg.vhdr
%     - sub-Vog001_task-Eye_eeg.vmrk
%     - sub-Vog001_task-Eye_events.json
%     - sub-Vog001_task-Eye_events.tsv
%   - sub-Vog001_scans.tsv with filename and acq_time
% - dataset_description.json
% - participants.json
% - participants.tsv with id, age, gender etc


% This function provides a comprehensive example of using the bids_export
% function. Note that eventually, you may simply use bids_export({file1.set file2.set}) 
% and that all other parameters are highly recommended but optional.
%
% You need the raw data files to run this script
% These can be downloaded from https://openneuro.org/datasets/ds001787
% Then enter the path in this script to re-export the BDF files
%
% Arnaud Delorme - Jan 2019

% Education
% 1: Kein Abschluss
% 2: Hauptschul-/Volksschulabschluss / obligatorische Schule
% 3: Berufslehre, Berufsschule, Berufsmittelschule
% 4: Abitur / Matura
% 5: Fachhochschulabschluss
% 6: Universität / Hochschulabschluss
% 7: Anderer Abschluss


%% create data structure with all subjects

clear data generalInfo pInfo pInfoDesc eInfoDesc README CHANGES stimuli code tInfo chanlocs;

% list folders
d = dir(formatted_data);
d = d(~contains({d.name}, {'.', '..', '.DS_Store'}));

d_demo = dir(demoDataPath);
d_demo = d_demo(~contains({d_demo.name}, {'.', '..', '.DS_Store'}));

% init struct
data = struct();
sub_count = 1;
headers = {'participant_id', 'species', 'age', 'sex', 'handedness', 'tod', 'education_year', 'education_degree', 'replication', 'lab', 'old_id'};
demo_table = table('Size', [0 11], ...           
                  'VariableNames', headers, ...
                  'VariableTypes', {'cell', 'cell', 'cell', 'cell', 'cell', 'cell' , 'cell', 'cell', 'cell', 'cell', 'cell'});
ids = arrayfun(@(x) sprintf('%03d', x), 1:320, 'UniformOutput', false);
labs = {'Dartmouth', 'Florida', 'Mainz', 'Münster', 'NorthDakota', 'Ohio', 'Reykjavik', 'Sheffield', 'Zürich (Prof. Langer)', 'Zürich (Prof. Sauseng)'};


% loop over all labs
for lab = 1 : size(d, 1)
    
    
    % find and read corresponding demographic data
    clear demo
    i_lab = ismember({d_demo.name},  d(lab).name);
    demo = readtable(fullfile(d_demo(i_lab).folder, d_demo(i_lab).name, strcat('demo_', d_demo(i_lab).name, '.xlsx')));
    if isnumeric(demo.TimeOfRecording)
        demo.TimeOfRecording = string(datetime(demo.TimeOfRecording,'ConvertFrom','datenum','Format','HH:mm'));
        if any(ismissing(demo.TimeOfRecording))
            demo.TimeOfRecording(ismissing(demo.TimeOfRecording)) = "NA";
        end
    end

    % list subjects within a lab
    d_lab = dir(fullfile(d(lab).folder, d(lab).name));
    d_lab = d_lab(~contains({d_lab.name}, {'.', '..', '.DS_Store'}));

    % loop over subjects 
    for sub = 1 : size(d_lab, 1)

        d_sub = dir(fullfile(d_lab(sub).folder, d_lab(sub).name));

        % find CDA and Eye
        i_cda = find(contains({d_sub.name}, 'CDA'));
        i_eye = find(contains({d_sub.name}, 'Eye'));
        
        % check if CDA exist
        if(isempty(i_cda))
            continue
        end

        % check if Eye task is available
        if not(isempty(i_eye))
            % create data structure
            data(sub_count).file = {fullfile(d_sub(i_cda).folder, d_sub(i_cda).name)
                                    fullfile(d_sub(i_eye).folder, d_sub(i_eye).name)};
            data(sub_count).task     = {'CDA', 'Eye'};
        else
            % only CDA available
            data(sub_count).file = {fullfile(d_sub(i_cda).folder, d_sub(i_cda).name)};
            data(sub_count).task     = {'CDA'};
        end
        
        % participant information for participants.tsv file
        oldID = strsplit(d_sub(i_cda).name, '_');
        oldID = oldID{1};
        
        % find corresponding row
        if exist('demo') & size(demo,1) > 2
            i_row = demo.ID == str2num(oldID);
            AGE = demo.Age(i_row);
            SEX = demo.Sex{i_row};
            HAND = demo.Handedness{i_row};
            TOD = demo.TimeOfRecording{i_row};
            EDU_YEAR = demo.EducationYears{i_row};
            EDU_DEG = demo.Education(i_row);
        else
            AGE = {};
            SEX = {};
            HAND = {};
            TOD = {};
            EDU_YEAR = {};
            EDU_DEG = {};
        end
        
        % define patterns
        SEX_lower = lower(SEX); 
        female_patterns = ["f","w","frau","female","weiblich","weibl"];
        male_patterns   = ["m","mann","male","männlich","maennlich"];
        other_patterns  = ["other","they", "o"];
        na_patterns = ["na"];
        if ismember(SEX_lower,female_patterns)
            SEX = 'f';
        elseif ismember(SEX_lower,male_patterns)
            SEX = 'm';
        elseif ismember(SEX_lower,other_patterns)
            SEX = 'o';
        elseif ismember(SEX_lower,na_patterns)
            SEX = {};
        elseif isempty(SEX)
            SEX = {};
        else
            error('check me')
        end
        % convert NA to {}
        if strcmp(HAND, 'NA')
            HAND = {};
        end
        if strcmp(TOD, 'NA')
            TOD = {};
        end
        if strcmp(EDU_YEAR, 'NA')
            EDU_YEAR = {};
        end
        if strcmp(EDU_DEG, 'NA')
            EDU_DEG = {};
        end

        if strcmp(AGE, 'NA')
            AGE = {};
        end
        
        % insert to table
        demo_table.participant_id{sub_count} = strcat('Vog', ids{sub_count});
        demo_table.species{sub_count} = 'homo sapiens';
        demo_table.age{sub_count} = AGE;
        demo_table.sex{sub_count} = SEX;
        demo_table.handedness{sub_count} = HAND;
        demo_table.tod{sub_count} = TOD;
        demo_table.education_year{sub_count} = EDU_YEAR;
        demo_table.education_degree{sub_count} = EDU_DEG;
        demo_table.replication{sub_count} = 'Vogel & Machizawa (2004)';
        demo_table.lab{sub_count} = labs{lab};
        demo_table.old_id{sub_count} = oldID;
        
        % update counter
        sub_count = sub_count + 1;

    end
end


pInfo = table2cell(demo_table);
pInfo = [headers; pInfo];
        
%% Descriptions

% general information for dataset_description.json file
% -----------------------------------------------------
generalInfo.Name = 'EEGManyLabs: Repication of Vogel & Machizawa (2004)';
generalInfo.BIDSVersion = '1.8.0';
generalInfo.DatasetType = 'raw';
generalInfo.License = 'CC BY';
generalInfo.Acknowledgements = 'We thank all contributing labs for collecting these data and gnode for hosting the dataset.';
generalInfo.HowToAcknowledge = 'Please cite the EEGManyLabs whitepaper (https://doi.org/10.1016/j.cortex.2021.03.013) and the replication paper (tbd)';
generalInfo.Funding = {
                        'DFG (PA 4005/1-1) provided to Y.G. Pavlov', ...
                        'BBSRC (BB/X008428/1) provided to F. Mushtaq'
                    };
generalInfo.EthicsApprovals = 'Ethics approvals obtained by the replicating labs; see replication paper for details';
generalInfo.ReferencesAndLinks = {'https://osf.io/pbr8c/overview'};
generalInfo.Authors = {'Dawid Strzelczyk', 'Peter E. Clayson', 'Heida Maria Sigurdardottir', ...
'Hélène Devillez', 'Anton Lukashevich', 'Harold A. Rocha', ...
'Yong Hoon Chung', 'Kevin M. Ortego', 'Viola S. Stoermer', 'José C. García Alanis', ...
'Christoph Löffler', 'Anna-Lena Schubert', 'Anna Lena Biel', 'Ariane Tretow', ...
'Weiyong Xu', 'Jarmo Hämäläinen', 'Zitong Lu', 'Yong Min Choi', 'Eva Lout', ...
'Julie D. Golomb', 'Shuangke Jiang', 'Myles Jones', 'Eda Mizrak', 'Claudia C. von Bastian', ...
'Niko A. Busch', 'Maro G. Machizawa', 'William X. Q. Ngiam', 'Edward K. Vogel', ...
'Samuel A. Birkholz', 'Emily Johnson', 'Jeffrey S. Johnson', 'Charline Peylo', 'Larissa Behnke', ...
'Yannik Hilla', 'Paul Sauseng', 'Faisal Mushtaq', 'Yuri G. Pavlov', 'Nicolas Langer'};


% participant column description for participants.json file - check OSF!!!
% ---------------------------------------------------------
pInfoDesc = struct;
pInfoDesc.participant_id.Description = 'Unique participant identifier';
pInfoDesc.species.Description = 'species of the participants';
pInfoDesc.age.Description = 'age of the participant';
pInfoDesc.age.Units       = 'years';
pInfoDesc.sex.Description = 'sex of the participant as reported by the participant';
pInfoDesc.sex.Levels.m = 'male';
pInfoDesc.sex.Levels.f = 'female';
pInfoDesc.sex.Levels.o = 'other';
pInfoDesc.handedness.Description = 'handedness of the participant as determined by Edinburgh Handedness Inventory (EHI)';
pInfoDesc.handedness.Levels.l = 'left';
pInfoDesc.handedness.Levels.r = 'right';
pInfoDesc.handedness.Levels.a = 'ambidextrous';
pInfoDesc.replication.Description = 'EEGManyLabs replication project';
pInfoDesc.replication.Levels.VogelMachizawa2004 = 'Vogel & Machizawa (2004)';
pInfoDesc.tod.Description = 'Time of day; local time of data collection (HH:MM)';
pInfoDesc.lab.Description = 'EEGManyLabs replication Team';
pInfoDesc.old_id.Description = 'Original Participant ID during EEGManyLabs data collection';
pInfoDesc.education_year.Description = 'Total years of education; self report';
pInfoDesc.education_degree.Description = 'Highest educational degree; self report';
pInfoDesc.education_degree.Levels.l1 = 'No formal degree';
pInfoDesc.education_degree.Levels.l2 = 'Lower secondary school';
pInfoDesc.education_degree.Levels.l3 = 'Vocational training';
pInfoDesc.education_degree.Levels.l4 = 'Matura (upper secondary school qualification)';
pInfoDesc.education_degree.Levels.l5 = 'University of Applied Sciences';
pInfoDesc.education_degree.Levels.l6 = 'University degree (academic higher education)';
pInfoDesc.education_degree.Levels.l7 = 'Other qualification';

% event column description for xxx-events.json file (only one such file)
% ----------------------------------------------------------------------
eInfoDesc = struct(); % check bids_export_dawid
eInfo = {}; % check bids_export_dawid

% Content for README file
% -----------------------
README = sprintf( [ 'EEGManyLabs - Replication Raw Dataset - Vogel & Machizawa (2004)\n\n' ...
    ...
    'Raw data for the #EEGManyLabs replication of Vogel & Machizawa seminal study: ' ...
    ' Neural activity predicts individual differences in visual working memory capacity (2004).\n\n'...
    ...
    'BE AWARE - EEG data from multiple different sites. Recordings have different sampling rates, ' ...
    'channel layouts, recording references, and have not been rereferenced.\n\n'...
    ...
    'See also:\n' ...
    'Replication paper: xxx\n' ...
    'OSF page of the replication: https://osf.io/pbr8c/overview\n' ...
    'Replicated publication: Vogel & Machizawa, 2004; https://www.nature.com/articles/nature02447\n' ...
    'EEGManyLabs project: Pavlov et al., 2021; https://doi.org/gjx8qt\n' ...
    ]);
 

% % List of script to run the experiment
% % ------------------------------------
% code = { fullfile( dataPath, 'code', 'run_mw_experiment6.m') mfilename('fullpath') };

% Task information for xxxx-eeg.json file - details in bids_export_dawid.m
% ---------------------------------------
tInfo.PowerLineFrequency = struct();

     
% call to the export function
% ---------------------------
bids_export_dawid(data, 'targetdir', bidspath, 'taskName', 'EEGManyLabs', ...
            'modality', 'eeg', ...
            'exportformat', 'vhdr', 'elecexport', 'on', 'eventexport', 'on', ...
            'gInfo', generalInfo, ...,
            'pInfo', pInfo, 'pInfoDesc', pInfoDesc, ...
            'eInfoDesc', eInfoDesc, 'eInfo', eInfo, 'individualEventsJson', 'on', ...
            'README', README, ...
            'tInfo', tInfo, ...
            'deleteExportDir', 'off');



%