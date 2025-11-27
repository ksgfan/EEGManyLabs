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

  ##### Equivalence Test for H1.3.
  cohens_d_h1.3 = c()
  cohens_d_se_h1.3  = c()
  for (lab in 1 : length(labs)){
    n = dat$num_subs[lab] # final sample
    m1 = dat$setsize4_m[lab] # mean CDA amplitude for set size 4
    sd1 = dat$setsize4_sd[lab] # sd CDA amplitude for set size 4
    m2 = dat$setsize6_m[lab] # mean CDA amplitude for set size 6
    sd2 = dat$setsize6_sd[lab] # sd CDA amplitude for set size 6
    r12	= dat$corr_4_6[lab] # observed correlation of dependent variable between CDA amplitude set size 4 and 6
    
    # do the equivalence test
    res = TOSTpaired(n=n,m1=m1,m2=m2,sd1=sd1,sd2=sd2,r12 = r12, low_eqbound_dz=-0.36,high_eqbound_dz=0.36)
    
    # calculate cohen's d from the t-test
    cohens_d_h1.3[lab] = res$TOST_t1 / sqrt(n) # should be always the same? (i.e., t1 or t2)
    cohens_d_se_h1.3[lab] = sqrt(( (2 * (1 - r12)) / n ) + ( cohens_d_h1.3[lab]^2 / (2 * n) ))
  
  }
  
  ##### Now do only 1 equivalence test for H1.3 that treats all data as coming from 1 lab
  
  # first, concatenate data
  n = sum(dat$num_sub)
  weights = dat$num_subs
  m1 = weighted.mean(dat$setsize4_m, weights) # mean CDA amplitude for set size 4
  sd1 = weighted.mean(dat$setsize4_sd, weights) # sd CDA amplitude for set size 4
  m2 = weighted.mean(dat$setsize6_m, weights) # mean CDA amplitude for set size 6
  sd2 = weighted.mean(dat$setsize6_sd, weights) # sd CDA amplitude for set size 6
  r12	= weighted.mean(dat$corr_4_6, weights) # observed correlation of dependent variable between CDA amplitude set size 4 and 6
  
  png(file.path(figure_path, paste0(pipes_labels[pip], "_TOST_H1_1.png")), 
                width = 15, height = 7, units = "cm", res = 300)
  res = TOSTpaired(n=n,m1=m1,m2=m2,sd1=sd1,sd2=sd2,r12 = r12, low_eqbound_dz=-0.36,high_eqbound_dz=0.36)
  dev.off()
  
  ##### Equivalence Test for H2.2.
  # reverse the correlation to be positive
  dat$wm_corr_4_6_r = dat$wm_corr_4_6_r*(-1)
  
  r_h2.2 = c()
  r_se_h2.2 = c()
  p1 = c()
  p2 = c()
  for (lab in 1 : length(labs)){
    n = dat$num_subs[lab] # final sample
    obs_r = dat$wm_corr_4_6_r[lab] # observed correlation between VWM capacity and CDA amplitude increase from 4 to 6
    res = TOSTr(n = n, r = obs_r, low_eqbound_r=-0.299, high_eqbound_r= 0.299, alpha=0.05, plot = TRUE, verbose = TRUE)
  
    # calculate cohen's d from the t-test
    r_h2.2[lab] = res$r
    r_se_h2.2[lab] = (1 - r_h2.2[lab]^2) / sqrt(n-3)
    p1[lab] = res$TOST_p1
    p2[lab] = res$TOST_p2
    
  }
  
  ##### Now do only 1 equivalence test for H2.2 that treats all data as coming from 1 lab
  n = sum(dat$num_subs)
  weights = dat$num_subs
  obs_r = weighted.mean(dat$wm_corr_4_6_r, weights)
  
  png(file.path(figure_path, paste0(pipes_labels[pip], "_TOST_H2_2.png")), 
      width = 15, height = 7, units = "cm", res = 300)
  res = TOSTr(n = n, r = obs_r, low_eqbound_r=-0.299, high_eqbound_r= 0.299, alpha=0.05, plot = TRUE, verbose = TRUE)
  dev.off()
  
  # save the variables to a file
  save(cohens_d_h1.3, cohens_d_se_h1.3, r_h2.2, r_se_h2.2, 
       file = file.path(project_path, 'data/csv_files/', paste0(pipes_labels[pip], '_eq_test_results.Rda')))


} # end loop
