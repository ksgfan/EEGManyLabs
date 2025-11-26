rm(list = ls(all.names = TRUE))
#install.packages("tidyverse")
#install.packages("meta")
#install.packages("R.matlab")
#install.packages("DescTools")
library(tidyverse)
library(meta)
library(DescTools)

# paths
project_path = '/Volumes/G_PSYPLAFOR_methlab$/EEGManyLabs'
figure_path = file.path(project_path, 'figures')
setwd(file.path(project_path, 'scripts'))

# pipelines
pipes_labels = c('direct', 'advanced', 'ica', 'keep_all');

# init empty df
df_results = data.frame()

# loop over pipelines and save results
for (pip in 1 : length(pipes_labels)){
  
  # load data
  dat = read.csv(file.path(project_path, 'data/csv_files', paste0('all_stats_', pipes_labels[pip], '.csv')), header = T)
  labs = dat$Lab