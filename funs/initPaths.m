%%
clc
clear
restoredefaultpath


%% set paths

prefix = '/Volumes';

% project folder path
projectpath = fullfile(prefix, '/G_PSYPLAFOR_methlab$/EEGManyLabs/');

% raw data parth (as recieved from other labs)
rawdatapath = fullfile(prefix, '/G_PSYPLAFOR_methlab_data$/EEGManyLabs/raw_data/');

% path to all formatted, preprocessed data
datapath = fullfile(projectpath, '/data/');

% formatted data (merged, added chanlocs, etc)
formatted_data = fullfile(datapath, 'formatted_data/'); % result folder for formatted data (beh files merged, chanlocs loaded, etc.)

% path to data in BIDS format
bidspath = fullfile(datapath, 'dataBIDS/');

% path to demographic data
demoDataPath = fullfile(datapath, 'csv_files/demo_data/');

% paths to preprocessed data
preprocessedDirect = fullfile(datapath, 'prepDirect/');
preprocessedAdvanced = fullfile(datapath, 'prepAdvanced/');
preprocessedICA = fullfile(datapath, 'prepICA/');

% path to folder with functions
funspath = fullfile(projectpath, 'funs');
eeglabpath = fullfile(funspath, 'eeglab2025.0.0');
pathEdf2Asc = fullfile(funspath, 'edf2asc-mac');  
chanlocspath = fullfile(funspath, 'chanloc_files');


%% EEGlab
p = pwd;
cd(fullfile(eeglabpath))
eeglab
close()
cd(p)

% make sure that option_single is 0 (keep data in double) and
% option_computeica is 1 (precompute ICA activations)
evalc("pop_editoptions('option_single', 0, 'option_computeica', 1);");

%% add paths
addpath(genpath(fullfile(eeglabpath, '/plugins/Fileio20210601')))
addpath(funspath)
addpath(chanlocspath)


%