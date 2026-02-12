################################################################################
### Bayesian stats #############################################################
################################################################################

rm(list = ls(all.names = TRUE))
project_path = '/Volumes/G_PSYPLAFOR_methlab$/EEGManyLabs'
setwd(file.path(project_path, 'scripts'))

# install.packages("bayestestR")
library(bayestestR)
library(brms)
library(tidyverse)
#library(loo)

pipes_labels = c('direct', 'advanced', 'ica');

# init empty df for summary stats
df_results = data.frame()

ci_level = 0.98

# loop over pipelines and save results
for (pip in 1 : length(pipes_labels)){
  
  # load data
  dat = read.csv(file.path(project_path, 'data/csv_files', paste0('subject_averages_', pipes_labels[pip], '.csv')), header = T)
  dat$Hand[dat$Hand == "n/a"] = NA
  dat$Sex[dat$Sex == "n/a"] = NA
  # exclude gender 'o' and hand 'a': not enough data
  dat$Hand[dat$Hand == "a"] = NA 
  dat$Sex[dat$Sex == "o"] = NA

  ##############################################################################
  ### Statistical analysis for Hypothesis Outcome Neutral ######################
  ##############################################################################
  
  ##############################################################################
  ### Linear mixed effect model ################################################
  ##############################################################################
  
  # convert to long format
  df_long = dat %>%
    pivot_longer(
      cols = c('Contra', 'Ipsi'),  
      names_to  = "Cluster",
      values_to = "CDA"
    )
  
  # convert vars to factor, scale CDA
  df = df_long %>%
    mutate(
      Cluster = factor(Cluster, levels = c("Contra", "Ipsi")),
      Sex  = factor(Sex, levels = c("f", "m")),
      Hand     = factor(Hand, levels = c("r", "l")), 
      Lab    = factor(Lab),
      ID = factor(ID),
      
      # Scale CDA (Gelman 2007: SD = 0.5)
      CDA_s = as.numeric(scale(CDA)) / 2
    )
  
  # sanity check
  levels(df$Cluster)
  levels(df$Hand)
  levels(df$Hand)
  contrasts(df$Hand)
  contrasts(df$Sex)
  
  # priors
  priors = c(
    prior(student_t(1, 0, 2.5), class = "b"),          # all fixed effects (Cluster, Sex, Hand, Lab)
    prior(student_t(3, 0, 10), class = "Intercept"),   # weak prior for intercept
    prior(student_t(3, 0, 10), class = "sigma")        # residual SD
  )
  
  # model formula with covariates and lab as random effect
  formula_cda <- bf(
    CDA_s ~ Cluster * Sex * Hand + (1|ID) + (1|Lab),
    family = gaussian()
  )
  
  # fit model
  fit_cda = brm(
    formula = formula_cda,
    data = df,
    prior = priors,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    cores = 4,
    sample_prior = "no",
    control = list(adapt_delta = 0.99)
  )
  summary(fit_cda, prob = ci_level)

  # Extract exact posterior
  post = posterior_summary(
    fit_cda,
    pars = "^b_",
    probs = c(0.01, 0.99)
  )
  post
  
  # Mark effects whose 95% CI excludes 0
  res = as.data.frame(post) %>%
    mutate(
      parameter = rownames(post),
      sig98 = if_else(Q1 > 0 | Q99 < 0, "CI excludes 0*", "CI includes 0")
    )
  res
  
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'ON_LMM',
      Param = res$parameter,
      SMD = res$Estimate,
      CI1 = res$Q1,
      CI2 = res$Q99,
      sig98 = res$sig98
    )
  )
  
  ##############################################################################
  ### Sequential updating for Outcome Neutral ##################################
  ##############################################################################
  
  # convert to long format
  df_long = dat %>%
    pivot_longer(
      cols = c('Contra', 'Ipsi'),  
      names_to  = "Cluster",
      values_to = "CDA"
    )
  labs = unique(df_long$Lab)
  
  # formula
  formula_h1_su = bf(
    CDA_s ~ Cluster + (1|ID), # + Sex + Hand #### contrast need at lweast 2 levels! sometimes a lab has only right-handed subject, or e.g. demographics are missing (n/a)
    family = gaussian()
  )
  
  # priors
  priors_initial = c(
    prior(student_t(1, 0, 2.5), class = "b"),
    prior(student_t(3, 0, 10), class = "Intercept"),
    prior(student_t(3, 0, 10), class = "sigma")
  )
  
  # Sequential Bayesian updating 
  fits = list()
  current_priors = priors_initial
  
  for (i in seq_along(labs)) {
    lab_name = labs[i]
    data_lab = df_long[df_long$Lab == lab_name, ]
    
    df_lab = data_lab %>%
      mutate(
        Cluster = factor(Cluster),
        Sex  = factor(Sex, levels = c("f", "m")),
        Hand     = factor(Hand, levels = c("r", "l")), 
        Lab    = factor(Lab),
        ID = factor(ID),
        
        # Scale CDA (Gelman 2007: SD = 0.5)
        CDA_s = as.numeric(scale(CDA)) / 2
      )
    
    # print progress
    cat("Fitting model for lab:", lab_name, "\n")
    
    fit = brm(
      formula = formula_h1_su,
      data    = df_lab,
      prior   = current_priors,
      chains  = 4,
      iter    = 4000,
      warmup  = 2000,
      cores   = 4,
      sample_prior = "yes", 
      control = list(adapt_delta = 0.99)
    )
    fits[[lab_name]] = fit
    
    # update prioirs
    post_summ = posterior_summary(fit, probs = c(0.01, 0.99))
    post_b = post_summ[grepl("^b_", rownames(post_summ)), ]
    
    current_priors = c(
      set_prior(
        paste0("normal(", post_b["b_Intercept",   "Estimate"], ", ", post_b["b_Intercept",   "Est.Error"], ")"),
        class = "Intercept"
      ),
      set_prior(
        paste0("normal(", post_b["b_ClusterIpsi", "Estimate"], ", ", post_b["b_ClusterIpsi", "Est.Error"], ")"),
        class = "b",
        coef  = "ClusterIpsi"
      ),
      set_prior("student_t(3, 0, 10)", class = "sigma")
    )
  }
  
  # save to a file
  saveRDS(fits, file.path(project_path, 'data/csv_files/', paste0(pipes_labels[pip], 'ON_sequential_results.rds')))
  
  # final posterir
  final_fit = fits[[length(fits)]]
  summary(final_fit)
  
  # Extract fixed-effect summaries (Est, Error, CI)
  post = posterior_summary(
    final_fit,                         
    pars = "^b_",
    probs = c(0.01, 0.99)
  )
  
  # Convert to dataframe
  res = as.data.frame(post)
  
  # Add parameter labels
  res$parameter = rownames(post)
  
  # Flag whether the 95% CI excludes zero
  res$sig98 = ifelse(
    res$Q1 > 0 | res$Q99 < 0,
    "CI excludes 0*",
    "CI includes 0"
  )
  res
  
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'ON_Seq',
      Param = res$parameter,
      SMD = res$Estimate,
      CI1 = res$Q1,
      CI2 = res$Q99,
      sig98 = res$sig98
    )
  )
  
  ################################################################################
  ### Outcome Neutral end ########################################################
  ################################################################################
  
  ################################################################################
  ### Statistical analysis for Hypothesis #1 #####################################
  ################################################################################
  
  ################################################################################
  ### Linear mixed effect model ##################################################
  ################################################################################
  
  # convert to long format
  df_long = dat %>%
    pivot_longer(
      cols = c('CDA2', 'CDA4', 'CDA6'),  
      names_to  = "SetSize",
      values_to = "CDA"
    )
  
  # convert vars to factor, scale CDA
  df = df_long %>%
    mutate(
      SetSize = factor(SetSize),
      Sex  = factor(Sex, levels = c("f", "m")),
      Hand     = factor(Hand, levels = c("r", "l")), 
      Lab    = factor(Lab),
      ID = factor(ID),
      
      # Scale CDA (Gelman 2007: SD = 0.5)
      CDA_s = as.numeric(scale(CDA)) / 2
    )
  
  # priors
  priors = c(
    prior(student_t(1, 0, 2.5), class = "b"),          # all fixed effects (SetSize, Sex, Hand, Lab)
    prior(student_t(3, 0, 10), class = "Intercept"),   # weak prior for intercept
    prior(student_t(3, 0, 10), class = "sigma")        # residual SD
  )
  
  # model formula with covariates and lab as random effect
  formula_cda = bf(
    CDA_s ~ SetSize * Sex * Hand + (1|ID) + (1|Lab),
    family = gaussian()
  )
  
  # fit model
  fit_cda = brm(
    formula = formula_cda,
    data = df,
    prior = priors,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    cores = 4,
    sample_prior = "no",
    control = list(adapt_delta = 0.99)
  )
  summary(fit_cda, prob = ci_level)
  # contrast 4 vs 6: option 1
  h = hypothesis(fit_cda, "SetSizeCDA6 - SetSizeCDA4 = 0")
  # contrast 4 vs 6: option 2
  #df_relevel = df
  #df_relevel$SetSize <- relevel(df_relevel$SetSize, ref = "CDA4")
  #fit_cda2 <- update(fit_cda, newdata = df_relevel)
  #summary(fit_cda, prob = ci_level)2)
  # equivalent to 'hypothesis' command
  
  
  # Extract exact posterior
  post = posterior_summary(
    fit_cda,
    pars = "^b_",
    probs = c(0.01, 0.99)
  )
  post
  
  # Mark effects whose 95% CI excludes 0
  res = as.data.frame(post) %>%
    mutate(
      parameter = rownames(post),
      sig98 = if_else(Q1 > 0 | Q99 < 0, "CI excludes 0*", "CI includes 0")
    )
  res
  # save results
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'H1_LMM',
      Param = res$parameter,
      SMD = res$Estimate,
      CI1 = res$Q1,
      CI2 = res$Q99,
      sig98 = res$sig98
    )
  )
  # save results set size 4 vs 6.
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'H1_LMM',
      Param = "4vs6",
      SMD = h$hypothesis$Estimate,
      CI1 = h$hypothesis$CI.Lower,
      CI2 = h$hypothesis$CI.Upper,
      sig98 = if_else(h$hypothesis$CI.Lower > 0 | h$hypothesis$CI.Upper < 0, "CI excludes 0*", "CI includes 0")
    )
  )
  
  ################################################################################
  ### Sequential updating for #H1 ################################################
  ################################################################################
  
  # convert to long format
  df_long = dat %>%
    pivot_longer(
      cols = c('CDA2', 'CDA4', 'CDA6'),  
      names_to  = "SetSize",
      values_to = "CDA"
    )
  labs = unique(df_long$Lab)
  
  # formula
  formula_h1_su = bf(
    CDA_s ~ SetSize + (1|ID), # + Sex + Hand #### contrast need at lweast 2 levels! sometimes a lab has only right-handed subject, or e.g. demographics are missing (n/a)
    family = gaussian()
  )
  
  # priors
  priors_initial = c(
    prior(student_t(1, 0, 2.5), class = "b"),
    prior(student_t(3, 0, 10), class = "Intercept"),
    prior(student_t(3, 0, 10), class = "sigma")
  )
  
  # Sequential Bayesian updating 
  fits = list()
  current_priors = priors_initial
  
  for (i in seq_along(labs)) {
    lab_name = labs[i]
    data_lab = df_long[df_long$Lab == lab_name, ]
    
    df_lab = data_lab %>%
      mutate(
        SetSize = factor(SetSize),
        Sex  = factor(Sex, levels = c("f", "m")),
        Hand     = factor(Hand, levels = c("r", "l")), 
        Lab    = factor(Lab),
        ID = factor(ID),
        
        # Scale CDA (Gelman 2007: SD = 0.5)
        CDA_s = as.numeric(scale(CDA)) / 2
      )
    
    # print progress
    cat("Fitting model for lab:", lab_name, "\n")
    
    fit = brm(
      formula = formula_h1_su,
      data    = df_lab,
      prior   = current_priors,
      chains  = 4,
      iter    = 4000,
      warmup  = 2000,
      cores   = 4,
      sample_prior = "yes", 
      control = list(adapt_delta = 0.99)
    )
    fits[[lab_name]] = fit
    
    # update prioirs
    post_summ = posterior_summary(fit, probs = c(0.01, 0.99))
    post_b = post_summ[grepl("^b_", rownames(post_summ)), ]
    
    current_priors = c(
      set_prior(
        paste0("normal(", post_b["b_Intercept",   "Estimate"], ", ", post_b["b_Intercept",   "Est.Error"], ")"),
        class = "Intercept"
      ),
      set_prior(
        paste0("normal(", post_b["b_SetSizeCDA4", "Estimate"], ", ", post_b["b_SetSizeCDA4", "Est.Error"], ")"),
        class = "b",
        coef  = "SetSizeCDA4"
      ),
      set_prior(
        paste0("normal(", post_b["b_SetSizeCDA6", "Estimate"], ", ", post_b["b_SetSizeCDA6", "Est.Error"], ")"),
        class = "b",
        coef  = "SetSizeCDA6"
      ),
      set_prior("student_t(3, 0, 10)", class = "sigma")
    )
  }
  
  # save to a file
  saveRDS(fits, file.path(project_path, 'data/csv_files/', paste0(pipes_labels[pip], 'H1_sequential_results.rds')))
  
  # final posterir
  final_fit = fits[[length(fits)]]
  summary(final_fit)
  
  # contrast 4 vs 6: option 1
  h = hypothesis(final_fit, "SetSizeCDA6 - SetSizeCDA4 = 0")
  
  # Extract fixed-effect summaries (Est, Error, CI)
  post = posterior_summary(
    final_fit,                         
    pars = "^b_",
    probs = c(0.01, 0.99)
  )
  
  # Convert to dataframe
  res = as.data.frame(post)
  
  # Add parameter labels
  res$parameter = rownames(post)
  
  # Flag whether the 95% CI excludes zero
  res$sig98 = ifelse(
    res$Q1 > 0 | res$Q99 < 0,
    "CI excludes 0*",
    "CI includes 0"
  )
  res
  # save results
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'H1_Seq',
      Param = res$parameter,
      SMD = res$Estimate,
      CI1 = res$Q1,
      CI2 = res$Q99,
      sig98 = res$sig98
    )
  )
  # save results set size 4 vs 6.
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'H1_Seq',
      Param = "4vs6",
      SMD = h$hypothesis$Estimate,
      CI1 = h$hypothesis$CI.Lower,
      CI2 = h$hypothesis$CI.Upper,
      sig98 = if_else(h$hypothesis$CI.Lower > 0 | h$hypothesis$CI.Upper < 0, "CI excludes 0*", "CI includes 0")
    )
  )
  
  
  ################################################################################
  ### H1 end #####################################################################
  ################################################################################
  
  
  ################################################################################
  ### Statistical analysis for Hypothesis #2.1 ####################################
  ################################################################################
  
  ################################################################################
  ### Linear mixed effect model ##################################################
  ################################################################################
  
  # compute difference and change sign
  dat$CDA_24 = (dat$CDA4 - dat$CDA2)*(-1)
  
  # sanity check
  cor.test(dat$K, dat$CDA_24) 
  
  # scale
  df = dat %>%
    mutate(
      K_s  = as.numeric(scale(K)) / 2,
      CDA_24_s  = as.numeric(scale(CDA_24)) / 2,
      Sex  = factor(Sex, levels = c("f", "m")),
      Hand = factor(Hand, levels = c("r", "l")),
      Lab  = factor(Lab),
      ID = factor(ID),
    )
  
  # priors
  priors_h2 = c(
    # Informative prior for the K effect based on r = .78
    prior(student_t(3, 0.78, 0.25), class = "b", coef = "K_s"),
    
    # Weak priors for other fixed effects
    prior(student_t(1, 0, 2.5), class = "b", coef = "Sexm"),
    prior(student_t(1, 0, 2.5), class = "b", coef = "Handl"),
    
    # Intercept and sigma
    prior(student_t(3, 0, 10), class = "Intercept"),
    prior(student_t(3, 0, 10), class = "sigma")
  )
  
  
  # model with lab as random effect
  formula_h2_mixed = bf(
    CDA_24_s ~ K_s * Sex * Hand + (1 | Lab),
    family = gaussian()
  )
  fit_h2_mixed = brm(
    formula = formula_h2_mixed,
    data    = df,
    prior   = priors_h2,
    chains  = 4,
    iter    = 4000,
    warmup  = 2000,
    cores   = 4,
    control = list(adapt_delta = 0.99),
    sample_prior = "no"
  )
  # show summary
  summary(fit_h2_mixed)
  post_h2_mixed = posterior_summary(
    fit_h2_mixed,
    pars = "^b_",
    probs = c(0.01, 0.99)
  )
  post_h2_mixed
  res_h2_mixed = as.data.frame(post_h2_mixed) %>%
    mutate(
      parameter = rownames(post_h2_mixed),
      sig98 = if_else(Q1 > 0 | Q99 < 0,
                      "CI excludes 0*",
                      "CI includes 0")
    )
  res_h2_mixed
  
  # save results
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'H2.1_LMM',
      Param = res_h2_mixed$parameter,
      SMD = res_h2_mixed$Estimate,
      CI1 = res_h2_mixed$Q1,
      CI2 = res_h2_mixed$Q99,
      sig98 = res_h2_mixed$sig98
    )
  )
  
  ################################################################################
  ### Sequential updating for #H2.1 ##############################################
  ################################################################################
  
  # compute difference and change sign
  dat$CDA_24 = (dat$CDA4 - dat$CDA2)*(-1)
  
  labs = unique(dat$Lab)
  
  # formula
  formula_h2_su = bf(
    CDA_24_s ~ K_s, # + Sex + Hand #### contrast need at lweast 2 levels!
    family = gaussian()
  )
  
  # priors
  priors_initial = c(
    prior(student_t(1, 0, 2.5), class = "b"),
    prior(student_t(3, 0, 10), class = "Intercept"),
    prior(student_t(3, 0, 10), class = "sigma")
  )
  
  # Sequential Bayesian updating 
  fits_h2 = list()
  current_priors = priors_initial
  
  for (i in seq_along(labs)) {
    lab_name = labs[i]
    data_lab = dat[dat$Lab == lab_name, ]
    
    df_lab = data_lab %>%
      mutate(
        Sex  = factor(Sex, levels = c("f", "m")),
        Hand     = factor(Hand, levels = c("r", "l")), 
        Lab    = factor(Lab),
        ID = factor(ID),
        
        # Scale CDA (Gelman 2007: SD = 0.5)
        CDA_24_s = as.numeric(scale(CDA_24)) / 2,
        K_s = as.numeric(scale(K)) / 2
      )
    
    # print progress
    cat("Fitting model for lab:", lab_name, "\n")
    
    fit = brm(
      formula = formula_h2_su,
      data    = df_lab,
      prior   = current_priors,
      chains  = 4,
      iter    = 4000,
      warmup  = 2000,
      cores   = 4,
      sample_prior = "yes", 
      control = list(adapt_delta = 0.99)
    )
    fits_h2[[lab_name]] = fit
    
    # update prioirs
    post_summ = posterior_summary(fit, probs = c(0.01, 0.99))
    post_b = post_summ[grepl("^b_", rownames(post_summ)), ]
    
    current_priors = c(
      set_prior(
        paste0("normal(", post_b["b_Intercept",   "Estimate"], ", ", post_b["b_Intercept",   "Est.Error"], ")"),
        class = "Intercept"
      ),
      set_prior(
        paste0("normal(", post_b["b_K_s", "Estimate"], ", ", post_b["b_K_s", "Est.Error"], ")"),
        class = "b",
        coef  = "K_s"
      ),
      set_prior("student_t(3, 0, 10)", class = "sigma")
    )
  }
  
  # save to a file
  saveRDS(fits_h2, file.path(project_path, 'data/csv_files/', paste0(pipes_labels[pip], 'H2.1_sequential_results.rds')))
  
  # final posterir
  final_fit_h2 = fits_h2[[length(fits_h2)]]
  summary(final_fit_h2)
  
  # Extract fixed-effect summaries (Est, Error, CI)
  post_h2 = posterior_summary(
    final_fit_h2,                         
    pars = "^b_",
    probs = c(0.01, 0.99)
  )
  
  # Convert to dataframe
  res_h2 = as.data.frame(post_h2)
  
  # Add parameter labels
  res_h2$parameter = rownames(post_h2)
  
  # Flag whether the 95% CI excludes zero
  res_h2$sig98 = ifelse(
    res_h2$Q1 > 0 | res_h2$Q99 < 0,
    "CI excludes 0*",
    "CI includes 0"
  )
  res_h2
  
  # save results
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'H2.1_Seq',
      Param = res_h2$parameter,
      SMD = res_h2$Estimate,
      CI1 = res_h2$Q1,
      CI2 = res_h2$Q99,
      sig98 = res_h2$sig98
    )
  )
  
  ################################################################################
  ### H2.1 end ###################################################################
  ################################################################################
  
  ################################################################################
  ### Statistical analysis for Hypothesis #2.2 ####################################
  ################################################################################
  
  ################################################################################
  ### Linear mixed effect model ##################################################
  ################################################################################
  
  # compute difference and change sign
  dat$CDA_46 = (dat$CDA6 - dat$CDA4)*(-1)
  
  # sanity check
  cor.test(dat$K, dat$CDA_46) 
  
  # scale
  df = dat %>%
    mutate(
      K_s  = as.numeric(scale(K)) / 2,
      CDA_46_s  = as.numeric(scale(CDA_46)) / 2,
      Sex  = factor(Sex, levels = c("f", "m")),
      Hand = factor(Hand, levels = c("r", "l")),
      Lab  = factor(Lab),
      ID = factor(ID),
    )
  
  # priors
  priors_h2 = c(
    # Weak prior (exact correlation was not reported)
    prior(student_t(3, 0, 1), class = "b", coef = "K_s"),
    
    # Weak priors for other fixed effects
    prior(student_t(1, 0, 2.5), class = "b", coef = "Sexm"),
    prior(student_t(1, 0, 2.5), class = "b", coef = "Handl"),
    
    # Intercept and sigma
    prior(student_t(3, 0, 10), class = "Intercept"),
    prior(student_t(3, 0, 10), class = "sigma")
  )
  
  
  # model with lab as random effect
  formula_h2_mixed = bf(
    CDA_46_s ~ K_s * Sex * Hand + (1 | Lab),
    family = gaussian()
  )
  fit_h2_mixed = brm(
    formula = formula_h2_mixed,
    data    = df,
    prior   = priors_h2,
    chains  = 4,
    iter    = 4000,
    warmup  = 2000,
    cores   = 4,
    control = list(adapt_delta = 0.99),
    sample_prior = "no"
  )
  # show summary
  summary(fit_h2_mixed)
  post_h2_mixed = posterior_summary(
    fit_h2_mixed,
    pars = "^b_",
    probs = c(0.01, 0.99)
  )
  post_h2_mixed
  res_h2_mixed = as.data.frame(post_h2_mixed) %>%
    mutate(
      parameter = rownames(post_h2_mixed),
      sig98 = if_else(Q1 > 0 | Q99 < 0,
                      "CI excludes 0*",
                      "CI includes 0")
    )
  res_h2_mixed
  
  # save results
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'H2.2_LMM',
      Param = res_h2_mixed$parameter,
      SMD = res_h2_mixed$Estimate,
      CI1 = res_h2_mixed$Q1,
      CI2 = res_h2_mixed$Q99,
      sig98 = res_h2_mixed$sig98
    )
  )
  
  ################################################################################
  ### Sequential updating for #H2.2 ##############################################
  ################################################################################
  
  # compute difference and change sign
  dat$CDA_46 = (dat$CDA6 - dat$CDA4)*(-1)
  
  labs = unique(dat$Lab)
  
  # formula
  formula_h2_su = bf(
    CDA_46_s ~ K_s, # + Sex + Hand #### contrast need at lweast 2 levels!
    family = gaussian()
  )
  
  # priors
  priors_initial = c(
    prior(student_t(1, 0, 2.5), class = "b"),
    prior(student_t(3, 0, 10), class = "Intercept"),
    prior(student_t(3, 0, 10), class = "sigma")
  )
  
  # Sequential Bayesian updating 
  fits_h2 = list()
  current_priors = priors_initial
  
  for (i in seq_along(labs)) {
    lab_name = labs[i]
    data_lab = dat[dat$Lab == lab_name, ]
    
    df_lab = data_lab %>%
      mutate(
        Sex  = factor(Sex, levels = c("f", "m")),
        Hand     = factor(Hand, levels = c("r", "l")), 
        Lab    = factor(Lab),
        ID = factor(ID),
        
        # Scale CDA (Gelman 2007: SD = 0.5)
        CDA_46_s = as.numeric(scale(CDA_46)) / 2,
        K_s = as.numeric(scale(K)) / 2
      )
    
    # print progress
    cat("Fitting model for lab:", lab_name, "\n")
    
    fit = brm(
      formula = formula_h2_su,
      data    = df_lab,
      prior   = current_priors,
      chains  = 4,
      iter    = 4000,
      warmup  = 2000,
      cores   = 4,
      sample_prior = "yes", 
      control = list(adapt_delta = 0.99)
    )
    fits_h2[[lab_name]] = fit
    
    # update prioirs
    post_summ = posterior_summary(fit, probs = c(0.01, 0.99))
    post_b = post_summ[grepl("^b_", rownames(post_summ)), ]
    
    current_priors = c(
      set_prior(
        paste0("normal(", post_b["b_Intercept",   "Estimate"], ", ", post_b["b_Intercept",   "Est.Error"], ")"),
        class = "Intercept"
      ),
      set_prior(
        paste0("normal(", post_b["b_K_s", "Estimate"], ", ", post_b["b_K_s", "Est.Error"], ")"),
        class = "b",
        coef  = "K_s"
      ),
      set_prior("student_t(3, 0, 10)", class = "sigma")
    )
  }
  
  # save to a file
  saveRDS(fits_h2, file.path(project_path, 'data/csv_files/', paste0(pipes_labels[pip], 'H2.2_sequential_results.rds')))
  
  # final posterir
  final_fit_h2 = fits_h2[[length(fits_h2)]]
  summary(final_fit_h2)
  
  # Extract fixed-effect summaries (Est, Error, CI)
  post_h2 = posterior_summary(
    final_fit_h2,                         
    pars = "^b_",
    probs = c(0.01, 0.99)
  )
  
  # Convert to dataframe
  res_h2 = as.data.frame(post_h2)
  
  # Add parameter labels
  res_h2$parameter = rownames(post_h2)
  
  # Flag whether the 95% CI excludes zero
  res_h2$sig98 = ifelse(
    res_h2$Q1 > 0 | res_h2$Q99 < 0,
    "CI excludes 0*",
    "CI includes 0"
  )
  res_h2
  
  # save results
  df_results = rbind(
    df_results,
    data.frame(
      Pipeline = pipes_labels[pip],
      TestType = 'H2.2_Seq',
      Param = res_h2$parameter,
      SMD = res_h2$Estimate,
      CI1 = res_h2$Q1,
      CI2 = res_h2$Q99,
      sig98 = res_h2$sig98
    )
  )
  
  ################################################################################
  ### H2.2 end ###################################################################
  ################################################################################
  
  

} # loop end

# save result to a file.
save(df_results,
     file = file.path(project_path, 'data/csv_files/bayesian_analyses_results.Rda'))

df_results[, 4:6] = round(df_results[, 4:6], 2)

