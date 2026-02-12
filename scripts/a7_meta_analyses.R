
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
pipes_labels = c('direct', 'advanced', 'ica');

# init empty df for summary stats
df_results = data.frame()

ci_level = 0.98

lab_labels = c('DART', 'USF', 'JGU', 'WWU', 'NDSU', 'OSU', 'UI' , 'TUOS', 'UZH-MPR', 'UZH-NCN')


# loop over pipelines and save results
for (pip in 1 : length(pipes_labels)){

  # load data
  dat = read.csv(file.path(project_path, 'data/csv_files', paste0('all_stats_', pipes_labels[pip], '.csv')), header = T)
  # dat = dat[dat$Lab != "Florida",] # check the results, when Florida is excluded
  labs = lab_labels
  
  ############ outcome neutral test: contra vs ipsi
  d = dat$on_gz
  d_se = dat$on_gz_se
  
  meta_outcome_neutral = metagen(TE = d,
                      seTE = d_se,
                      studlab = labs,
                      data = NULL,
                      sm = "SMD",
                      common = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      title = "Contra vs. Ipsi",
                      level.ci = ci_level,
                      level = ci_level,
                      level.ma = ci_level,
                      level.predict = ci_level,
                      level.comb = ci_level)
  summary(meta_outcome_neutral)
  # save results
  df_results <- rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'OutcomeNeutral',
      SMD = meta_outcome_neutral$TE.random,
      CI1 = meta_outcome_neutral$lower.random,
      CI2 = meta_outcome_neutral$upper.random,
      t_stat = meta_outcome_neutral$statistic.random,
      p_val = meta_outcome_neutral$pval.random
    )
  )
  
  # plot results
  pdf(file.path(figure_path, paste0(pipes_labels[pip], "_outcome_neutral.pdf")), width = 8, height = 4)
  forest(meta_outcome_neutral, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       # leftlabs = c("Lab", expression("Cohen's d"), "SE", "98% CI"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  png(file.path(figure_path, paste0(pipes_labels[pip], "_outcome_neutral.png")), width = 8, height = 4, units = "in", res = 300)
  forest(meta_outcome_neutral, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  
  ############ correct vs incorrect trials
  d = dat$corrincorr_gz
  d_se = dat$corrincorr_gz_se
  meta_correct_vs_incorrect = metagen(TE = d,
                               seTE = d_se,
                               studlab = labs,
                               data = NULL,
                               sm = "SMD",
                               common = FALSE,
                               random = TRUE,
                               method.tau = "REML",
                               hakn = TRUE,
                               prediction = TRUE,
                               title = "Correct vs. Incorrect",
                               level.ci = ci_level,
                               level = ci_level,
                               level.ma = ci_level,
                               level.predict = ci_level,
                               level.comb = ci_level)
  summary(meta_correct_vs_incorrect)
  # save results
  df_results <- rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'Corr_vs_Incorr',
      SMD = meta_correct_vs_incorrect$TE.random,
      CI1 = meta_correct_vs_incorrect$lower.random,
      CI2 = meta_correct_vs_incorrect$upper.random,
      t_stat = meta_correct_vs_incorrect$statistic.random,
      p_val = meta_correct_vs_incorrect$pval.random
    )
  )
  
  # plot results
  pdf(file.path(figure_path, paste0(pipes_labels[pip], "_correct_vs_incorrect.pdf")), width = 8, height = 4)
  forest(meta_correct_vs_incorrect, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  png(file.path(figure_path, paste0(pipes_labels[pip], "_correct_vs_incorrectt.png")), width = 8, height = 4, units = "in", res = 300)
  forest(meta_correct_vs_incorrect, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  
  ############ set size 2 vs 4
  d = dat$ph_2_4_gz
  d_se = dat$ph_2_4_gz_se
  meta_2_vs_4 = metagen(TE = d,
                         seTE = d_se,
                         studlab = labs,
                         data = NULL,
                         sm = "SMD",
                         common = FALSE,
                         random = TRUE,
                         method.tau = "REML",
                         hakn = TRUE,
                         prediction = TRUE,
                         title = "Set Size 2 vs Set Size 4",
                        level.ci = ci_level,
                        level = ci_level,
                        level.ma = ci_level,
                        level.predict = ci_level,
                        level.comb = ci_level)
  summary(meta_2_vs_4)
  # save results
  df_results <- rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'S2_vs_S4',
      SMD = meta_2_vs_4$TE.random,
      CI1 = meta_2_vs_4$lower.random,
      CI2 = meta_2_vs_4$upper.random,
      t_stat = meta_2_vs_4$statistic.random,
      p_val = meta_2_vs_4$pval.random
    )
  )
  
  # plot results
  pdf(file.path(figure_path, paste0(pipes_labels[pip], "_t_test_2_vs_4.pdf")), width = 8, height = 4)
  forest(meta_2_vs_4, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  png(file.path(figure_path, paste0(pipes_labels[pip], "_t_test_2_vs_4.png")), width = 8, height = 4, units = "in", res = 300)
  forest(meta_2_vs_4, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  
  ############ set size 2 vs 6
  d = dat$ph_2_6_gz
  d_se = dat$ph_2_6_gz_se
  meta_2_vs_6 = metagen(TE = d,
                     seTE = d_se,
                     studlab = labs,
                     data = NULL,
                     sm = "SMD",
                     common = FALSE,
                     random = TRUE,
                     method.tau = "REML",
                     hakn = TRUE,
                     prediction = TRUE,
                     title = "Set Size 2 vs Set Size 6",
                     level.ci = ci_level,
                     level = ci_level,
                     level.ma = ci_level,
                     level.predict = ci_level,
                     level.comb = ci_level)
  summary(meta_2_vs_6)
  # save results
  df_results <- rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'S2_vs_S6',
      SMD = meta_2_vs_6$TE.random,
      CI1 = meta_2_vs_6$lower.random,
      CI2 = meta_2_vs_6$upper.random,
      t_stat = meta_2_vs_6$statistic.random,
      p_val = meta_2_vs_6$pval.random
    )
  )
  
  # plot results
  pdf(file.path(figure_path, paste0(pipes_labels[pip], "_t_test_2_vs_6.pdf")), width = 8, height = 4)
  forest(meta_2_vs_6, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  png(file.path(figure_path, paste0(pipes_labels[pip], "_t_test_2_vs_6.png")), width = 8, height = 4, units = "in", res = 300)
  forest(meta_2_vs_6, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  
  ############ set size 4 vs 6
  d = dat$ph_4_6_gz
  d_se = dat$ph_4_6_gz_se
  meta_4_vs_6 = metagen(TE = d,
                     seTE = d_se,
                     studlab = labs,
                     data = NULL,
                     sm = "SMD",
                     common = FALSE,
                     random = TRUE,
                     method.tau = "REML",
                     hakn = TRUE,
                     prediction = TRUE,
                     title = "Set Size 4 vs Set Size 6",
                     level.ci = ci_level,
                     level = ci_level,
                     level.ma = ci_level,
                     level.predict = ci_level,
                     level.comb = ci_level)
  summary(meta_4_vs_6)
  # save results
  df_results <- rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'S4_vs_S6',
      SMD = meta_4_vs_6$TE.random,
      CI1 = meta_4_vs_6$lower.random,
      CI2 = meta_4_vs_6$upper.random,
      t_stat = meta_4_vs_6$statistic.random,
      p_val = meta_4_vs_6$pval.random
    )
  )
  
  # plot results
  pdf(file.path(figure_path, paste0(pipes_labels[pip], "_t_test_4_vs_6.pdf")), width = 8, height = 4)
  forest(meta_4_vs_6, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  png(file.path(figure_path, paste0(pipes_labels[pip], "_t_test_4_vs_6.png")), width = 8, height = 4, units = "in", res = 300)
  forest(meta_4_vs_6, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       leftcols = c("studlab", "TE", "seTE", "ci"),
       leftlabs = c("Lab", expression(italic("g")["z"]), "SE", "98% CI"),
       rightcols = c("w.random"))
  dev.off()
  
  ############ correlation of amplitude increase with WM capacity 2 to 4
  r = dat$wm_corr_2_4_r * (-1) ### flip sign
  meta_correlation = metacor(
    cor = r,
    n = dat$num_subs,
    studlab = labs,
    data = NULL,
    sm = "ZCOR", # concert to Fischer's z
    common = FALSE,
    random = TRUE,
    method.tau = "REML",
    hakn = TRUE,
    prediction = T,
    backtransf = TRUE, # convert Fishers'z back to correlation
    title = "Amplitude decrease vs WM capacity",
    level.ci = ci_level,
    level = ci_level,
    level.ma = ci_level,
    level.predict = ci_level,
    level.comb = ci_level
  )
  meta_correlation$Pearsons_r <- r
  summary(meta_correlation)
  forest(meta_correlation, 
         prediction = TRUE, 
         print.tau2 = FALSE,
         #leftcols = c("studlab", "TE", "seTE", "ci"),
         leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
         #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
         leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
         rightcols = c("w.random"),
         xlim = c(-0.7, 0.9),
         at = seq(-0.4, 0.8, by = 0.2))
  # save results
  df_results <- rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'Correlation_2-4_vs_VWM',
      SMD = meta_correlation$TE.random,
      CI1 = meta_correlation$lower.random,
      CI2 = meta_correlation$upper.random,
      t_stat = meta_correlation$statistic.random,
      p_val = meta_correlation$pval.random
    )
  )
  
  # plot results
  pdf(file.path(figure_path, paste0(pipes_labels[pip], "_corr_2_to_4.pdf")), width = 8, height = 4)
  forest(meta_correlation, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       #leftcols = c("studlab", "TE", "seTE", "ci"),
       leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
       #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
       leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
       rightcols = c("w.random"),
       xlim = c(-0.7, 0.9),
       at = seq(-0.4, 0.8, by = 0.2))
  dev.off()
  png(file.path(figure_path, paste0(pipes_labels[pip], "_corr_2_to_4.png")), width = 8, height = 4, units = "in", res = 300)
  forest(meta_correlation, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       #leftcols = c("studlab", "TE", "seTE", "ci"),
       leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
       #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
       leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
       rightcols = c("w.random"),
       xlim = c(-0.7, 0.9),
       at = seq(-0.4, 0.8, by = 0.2))
  dev.off()
  
  ############ correlation of amplitude at set size 6 with WMC (K)
  r = dat$wm_corr_6_r * (-1) ### flip sign
  meta_correlation_6 = metacor(
    cor = r,
    n = dat$num_subs,
    studlab = labs,
    data = NULL,
    sm = "ZCOR", # concert to Fischer's z
    common = FALSE,
    random = TRUE,
    method.tau = "REML",
    hakn = TRUE,
    prediction = T,
    backtransf = TRUE, # convert Fishers'z back to correlation
    title = "Amplitude decrease vs WM capacity",
    level.ci = ci_level,
    level = ci_level,
    level.ma = ci_level,
    level.predict = ci_level,
    level.comb = ci_level
  )
  meta_correlation_6$Pearsons_r <- r
  summary(meta_correlation_6)
  forest(meta_correlation_6, 
         prediction = TRUE, 
         print.tau2 = FALSE,
         #leftcols = c("studlab", "TE", "seTE", "ci"),
         leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
         #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
         leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
         rightcols = c("w.random"),
         xlim = c(-0.7, 0.9),
         at = seq(-0.4, 0.8, by = 0.2))
  # save results
  df_results <- rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'Correlation_6_vs_VWM',
      SMD = meta_correlation_6$TE.random,
      CI1 = meta_correlation_6$lower.random,
      CI2 = meta_correlation_6$upper.random,
      t_stat = meta_correlation_6$statistic.random,
      p_val = meta_correlation_6$pval.random
    )
  )
  
  # plot results
  pdf(file.path(figure_path, paste0(pipes_labels[pip], "_corr_6.pdf")), width = 8, height = 4)
  forest(meta_correlation_6, 
         prediction = TRUE, 
         print.tau2 = FALSE,
         #leftcols = c("studlab", "TE", "seTE", "ci"),
         leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
         #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
         leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
         rightcols = c("w.random"),
         xlim = c(-0.7, 0.9),
         at = seq(-0.4, 0.8, by = 0.2))
  dev.off()
  png(file.path(figure_path, paste0(pipes_labels[pip], "_corr_6.png")), width = 8, height = 4, units = "in", res = 300)
  forest(meta_correlation_6, 
         prediction = TRUE, 
         print.tau2 = FALSE,
         #leftcols = c("studlab", "TE", "seTE", "ci"),
         leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
         #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
         leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
         rightcols = c("w.random"),
         xlim = c(-0.7, 0.9),
         at = seq(-0.4, 0.8, by = 0.2))
  dev.off()
  
  ############ correlation of amplitude increase with WM capacity 4 to 6
  r = dat$wm_corr_4_6_r * (-1) ### flip sign
  meta_correlation_4_6 = metacor(
    cor = r,
    n = dat$num_subs,
    studlab = labs,
    data = NULL,
    sm = "ZCOR",
    common = FALSE,
    random = TRUE,
    method.tau = "REML",
    hakn = TRUE,
    prediction = T,
    backtransf = TRUE, 
    title = "Amplitude decrease vs WM capacity",
    level.ci = ci_level,
    level = ci_level,
    level.ma = ci_level,
    level.predict = ci_level,
    level.comb = ci_level
  )
  meta_correlation_4_6$Pearsons_r <- r
  summary(meta_correlation_4_6)
  # save results
  df_results <- rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'Correlation_4-6_vs_VWM',
      SMD = meta_correlation_4_6$TE.random,
      CI1 = meta_correlation_4_6$lower.random,
      CI2 = meta_correlation_4_6$upper.random,
      t_stat = meta_correlation_4_6$statistic.random,
      p_val = meta_correlation_4_6$pval.random
    )
  )
  
  # plot results
  pdf(file.path(figure_path, paste0(pipes_labels[pip], "_corr_4_to_6.pdf")), width = 8, height = 4)
  forest(meta_correlation_4_6, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       #leftcols = c("studlab", "TE", "seTE", "ci"),
       leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
       #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
       leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
       rightcols = c("w.random"),
       xlim = c(-1, 1),
       at = seq(-1, 1, by = 0.5))
  dev.off()
  png(file.path(figure_path, paste0(pipes_labels[pip], "_corr_4_to_6.png")), width = 8, height = 4, units = "in", res = 300)
  forest(meta_correlation_4_6, 
       prediction = TRUE, 
       print.tau2 = FALSE,
       #leftcols = c("studlab", "TE", "seTE", "ci"),
       leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
       #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
       leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
       rightcols = c("w.random"),
       xlim = c(-1, 1),
       at = seq(-1, 1, by = 0.5))
  dev.off()
  
  
  ############ correlation of amplitude increase with WM capacity (d-prime) 2 to 4
  r = dat$wm_d_corr_2_4_r * (-1) ### flip sign
  meta_correlation_d_prime = metacor(
    cor = r,
    n = dat$num_subs,
    studlab = labs,
    data = NULL,
    sm = "ZCOR", # concert to Fischer's z
    common = FALSE,
    random = TRUE,
    method.tau = "REML",
    hakn = TRUE,
    prediction = T,
    backtransf = TRUE, # convert Fishers'z back to correlation
    title = "Amplitude decrease vs WM capacity",
    level.ci = ci_level,
    level = ci_level,
    level.ma = ci_level,
    level.predict = ci_level,
    level.comb = ci_level
  )
  meta_correlation_d_prime$Pearsons_r <- r
  summary(meta_correlation_d_prime)
  forest(meta_correlation_d_prime, 
         prediction = TRUE, 
         print.tau2 = FALSE,
         #leftcols = c("studlab", "TE", "seTE", "ci"),
         leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
         #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
         leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
         rightcols = c("w.random"),
         xlim = c(-0.5, 0.9),
         at = seq(-0.4, 0.8, by = 0.2))
  # save results
  df_results <- rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'Correlation_2-4_vs_VWM_d_prime',
      SMD = meta_correlation_d_prime$TE.random,
      CI1 = meta_correlation_d_prime$lower.random,
      CI2 = meta_correlation_d_prime$upper.random,
      t_stat = meta_correlation_d_prime$statistic.random,
      p_val = meta_correlation_d_prime$pval.random
    )
  )
  
  # plot results
  pdf(file.path(figure_path, paste0(pipes_labels[pip], "_corr_2_to_4_d_prime.pdf")), width = 8, height = 4)
  forest(meta_correlation_d_prime, 
         prediction = TRUE, 
         print.tau2 = FALSE,
         #leftcols = c("studlab", "TE", "seTE", "ci"),
         leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
         #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
         leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
         rightcols = c("w.random"),
         xlim = c(-0.5, 0.9),
         at = seq(-0.4, 0.8, by = 0.2))
  dev.off()
  png(file.path(figure_path, paste0(pipes_labels[pip], "_corr_2_to_4_d_prime.png")), width = 8, height = 4, units = "in", res = 300)
  forest(meta_correlation_d_prime, 
         prediction = TRUE, 
         print.tau2 = FALSE,
         #leftcols = c("studlab", "TE", "seTE", "ci"),
         leftcols = c("studlab", "Pearsons_r", "seTE", "ci"),
         #leftlabs = c("Lab", expression("Fischer's z"), "SE", "98% CI"),
         leftlabs = c("Lab", expression("Pearson's r"), "SE", "98% CI"),
         rightcols = c("w.random"),
         xlim = c(-0.5, 0.9),
         at = seq(-0.4, 0.8, by = 0.2))
  dev.off()
  
  
  
  save(meta_4_vs_6, meta_correlation, meta_correlation_4_6,
       file = file.path(project_path, 'data/csv_files/', paste0(pipes_labels[pip], '_for_eq_tests.Rda')))
  
}

# save result to a file.
save(df_results,
     file = file.path(project_path, 'data/csv_files/meta_analyses_results.Rda'))

df_results[, 3:6] = round(df_results[, 3:6], 2)
