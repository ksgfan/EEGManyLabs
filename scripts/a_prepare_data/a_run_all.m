%% 

% load and inscpect all of the data
% merge CDA blocks
% concat and inspect behavioral data
% load chanlocs
% check reference
% save Eye, CDA and Resting as indivudual mat files

p = pwd;
funpath = strsplit(p, filesep);
addpath(fullfile(strjoin(funpath(1:end-2), filesep),'funs'))
initPaths;


%% run all 'prepare' scripts 

% prepare_dartmouth % OK, ch OK
% 
% prepare_florida % OK, ch OK
% 
% prepare_mainz % OK, ch OK
% 
% prepare_munster % OK, ch OK
% 
% prepare_north_dakota % OK, ch OK
% 
% prepare_ohio % OK, ch OK
% 
% prepare_reykjavik % OK, ch OK
% 
% prepare_sheffield % OK, ch OK

% prepare_zurich % OK, ch OK

% prepare_zurich_sauseng % OK, ch OK








