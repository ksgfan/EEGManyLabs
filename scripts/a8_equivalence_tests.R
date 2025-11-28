rm(list = ls(all.names = TRUE))

project_path = '/Volumes/G_PSYPLAFOR_methlab$/EEGManyLabs'
figure_path = file.path(project_path, 'figures')
setwd(file.path(project_path, 'scripts'))

library(tidyverse)
library(meta)
library(DescTools)
library(TOSTER)

# pipelines
pipes_labels = c('direct', 'advanced', 'ica', 'keep_all');

# loop over pipelines and save results
for (pip in 1 : length(pipes_labels)){
  
  # load data
  dat = read.csv(file.path(project_path, 'data/csv_files', paste0('all_stats_', pipes_labels[pip], '.csv')), header = T)
  labs = dat$Lab
  
  # equivalence test results
  load(file.path(project_path, 'data/csv_files/', paste0(pipes_labels[pip], '_for_eq_tests.Rda')))
  
  ##############################################################################
  ##############################################################################
  ##### Equivalence Test for H1.3.
  ##############################################################################
  ##############################################################################
  SESOI_d = 0.36
  
  ##### equivalence test for meta-analysis results
  png(file.path(figure_path, paste0(pipes_labels[pip], "_TOST_H1_1.png")), 
      width = 15, height = 7, units = "cm", res = 300)
  # do the test
  tost_H13 = TOSTmeta(
    ES = meta_4_vs_6$TE.random,
    se = meta_4_vs_6$seTE.random,
    low_eqbound_d  = -SESOI_d,
    high_eqbound_d =  SESOI_d,
    alpha = 0.02   # if you want to keep your p < .02 convention
  )
  dev.off()
  
  ##### Now do only 1 equivalence test for H1.3 that treats all data as coming from 1 lab
  # first, concatenate data
  n = sum(dat$num_sub)
  weights = dat$num_subs
  m1 = weighted.mean(dat$setsize4_m, weights) # mean CDA amplitude for set size 4
  sd1 = weighted.mean(dat$setsize4_sd, weights) # sd CDA amplitude for set size 4
  m2 = weighted.mean(dat$setsize6_m, weights) # mean CDA amplitude for set size 6
  sd2 = weighted.mean(dat$setsize6_sd, weights) # sd CDA amplitude for set size 6
  r12	= weighted.mean(dat$corr_4_6, weights) # observed correlation of dependent variable between CDA amplitude set size 4 and 6
  
  png(file.path(figure_path, paste0(pipes_labels[pip], "_TOST_H1_1_collapsed.png")), 
                width = 15, height = 7, units = "cm", res = 300)
  res = TOSTpaired(n=n,m1=m1,m2=m2,sd1=sd1,sd2=sd2,r12 = r12, low_eqbound_dz=-SESOI_d,high_eqbound_dz=SESOI_d)
  dev.off()
  

  ##############################################################################
  ##############################################################################
  ##### Equivalence Test for H2.2.
  ##############################################################################
  ##############################################################################
  SESOI_r = 0.299
  
  ##### equivalence test for meta-analysis results
  png(file.path(figure_path, paste0(pipes_labels[pip], "_TOST_H2_2.png")), 
      width = 15, height = 7, units = "cm", res = 300)
  # do the test
  tost_H22 = TOSTmeta(
    ES = meta_correlation_4_6$TE.random,
    se = meta_correlation_4_6$seTE.random,
    low_eqbound_d  = -SESOI_r,
    high_eqbound_d =  SESOI_r,
    alpha = 0.02   # if you want to keep your p < .02 convention
  )
  dev.off()
  
  ##### Now do only 1 equivalence test for H2.2 that treats all data as coming from 1 lab
  # reverse the correlation to be positive
  dat$wm_corr_4_6_r = dat$wm_corr_4_6_r*(-1)
  n = sum(dat$num_subs)
  weights = dat$num_subs
  obs_r = weighted.mean(dat$wm_corr_4_6_r, weights)
  
  png(file.path(figure_path, paste0(pipes_labels[pip], "_TOST_H2_2_collapsed.png")), 
      width = 15, height = 7, units = "cm", res = 300)
  res = TOSTr(n = n, r = obs_r, low_eqbound_r=-SESOI_r, high_eqbound_r = SESOI_r, alpha=0.05, plot = TRUE, verbose = TRUE)
  dev.off()
  
} # end loop
