rm(list = ls(all.names = TRUE))

library(officer)
library(flextable)

# paths
project_path = '/Volumes/G_PSYPLAFOR_methlab$/EEGManyLabs'
wordfiles_path = file.path(project_path, 'data', 'word_files')
setwd(file.path(project_path, 'scripts'))

# pipelines
pipes_labels = c('direct', 'advanced', 'ica', 'keep_all');

# fun for p-vals formatting
format_p <- function(x) {
  out <- ifelse(
    x < 0.001,
    formatC(x, format = "e", digits = 2),   # scientific for very small p
    formatC(x, format = "f", digits = 3)    # fixed, 3 decimals
  )
  out
}


# loop over pipelines and save results
for (pip in 1 : length(pipes_labels)){
  
  # load data
  dat = read.csv(file.path(project_path, 'data/csv_files', paste0('all_stats_', pipes_labels[pip], '.csv')), header = T)
  # names(dat)
  
  ##############################################################################
  # outcome neutral
  ##############################################################################
  tab = dat[ , c("Lab", "on_m_contra", "on_sd_contra", "on_m_ipsi", "on_sd_ipsi", 
                 "on_t_stat", "on_df", "on_p", "on_gz") ]
  # round all values except p-vals 
  is_num <- sapply(tab, is.numeric)
  p_cols <- grepl("_p$", names(tab))
  is_num[p_cols] <- FALSE
  tab[ , is_num] <- round(tab[ , is_num], 2)
  tab[ , p_cols ] <- format_p(tab[ , p_cols ])
  # save to word
  ft <- flextable(tab)
  ft <- set_header_labels(
                      ft,
                        Lab = 'Lab',
                        on_m_contra = "M",
                        on_sd_contra = 'SD',
                        on_m_ipsi = "M",
                        on_sd_ipsi = 'SD',
                        on_t_stat   = "t-value",
                        on_df  = "df",
                        on_p      = "p-value",
                        on_gz = "Hedges' gz"
                        
  )
  ft <- autofit(ft)
  ft <- padding(ft, padding.top = 0, padding.bottom = 0, part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- font(ft, fontname = "Arial", part = "all")
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- height(ft, height = 0.5, part = "body")
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft)
  print(doc, target = file.path(wordfiles_path, paste0(pipes_labels[pip], "_outcome_neutral.docx")))
  
  
  ##############################################################################
  # set size 2 vs 4
  ##############################################################################
  tab = dat[ , c("Lab", "setsize2_m", "setsize2_sd", "setsize4_m", "setsize4_sd", 
                 "ph_2_4_t_stat", "ph_2_4_df", "ph_2_4_p", "ph_2_4_gz") ]
  # round all values except p-vals 
  is_num <- sapply(tab, is.numeric)
  p_cols <- grepl("_p$", names(tab))
  is_num[p_cols] <- FALSE
  tab[ , is_num] <- round(tab[ , is_num], 2)
  tab[ , p_cols ] <- format_p(tab[ , p_cols ])
  # save to word
  ft <- flextable(tab)
  ft <- set_header_labels(
    ft,
    Lab = 'Lab',
    setsize2_m = "M",
    setsize2_sd = 'SD',
    setsize4_m = "M",
    setsize4_sd = 'SD',
    ph_2_4_t_stat   = "t-value",
    ph_2_4_df  = "df",
    ph_2_4_p      = "p-value",
    ph_2_4_gz = "Hedges' gz"
    
  )
  ft <- autofit(ft)
  ft <- padding(ft, padding.top = 0, padding.bottom = 0, part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- font(ft, fontname = "Arial", part = "all")
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- height(ft, height = 0.5, part = "body")
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft)
  print(doc, target = file.path(wordfiles_path, paste0(pipes_labels[pip], "_t_test_2_4.docx")))
  
  
  ##############################################################################
  # set size 2 vs 6
  ##############################################################################
  tab = dat[ , c("Lab", "setsize2_m", "setsize2_sd", "setsize6_m", "setsize6_sd", 
                 "ph_2_6_t_stat", "ph_2_6_df", "ph_2_6_p", "ph_2_6_gz") ]
  # round all values except p-vals 
  is_num <- sapply(tab, is.numeric)
  p_cols <- grepl("_p$", names(tab))
  is_num[p_cols] <- FALSE
  tab[ , is_num] <- round(tab[ , is_num], 2)
  tab[ , p_cols ] <- format_p(tab[ , p_cols ])
  # save to word
  ft <- flextable(tab)
  ft <- set_header_labels(
    ft,
    Lab = 'Lab',
    setsize2_m = "M",
    setsize2_sd = 'SD',
    setsize6_m = "M",
    setsize6_sd = 'SD',
    ph_2_6_t_stat   = "t-value",
    ph_2_6_df  = "df",
    ph_2_6_p      = "p-value",
    ph_2_6_gz = "Hedges' gz"
    
  )
  ft <- autofit(ft)
  ft <- padding(ft, padding.top = 0, padding.bottom = 0, part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- font(ft, fontname = "Arial", part = "all")
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- height(ft, height = 0.5, part = "body")
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft)
  print(doc, target = file.path(wordfiles_path, paste0(pipes_labels[pip], "_t_test_2_6.docx")))
  
  
  ##############################################################################
  # correlation of amplitude decrease (2 to 4) and VWM capacity
  ##############################################################################
  tab = dat[ , c("Lab", "wm_corr_2_4_amp_m", "wm_corr_2_4_amp_sd", "wm_corr_2_4_wm_m", "wm_corr_2_4_wm_sd", 
                 "wm_corr_2_4_r", "wm_corr_2_4_p") ]
  
  # make the correlation positive
  tab$wm_corr_2_4_r = tab$wm_corr_2_4_r*(-1)
  
  # round all values except p-vals 
  is_num <- sapply(tab, is.numeric)
  p_cols <- grepl("_p$", names(tab))
  is_num[p_cols] <- FALSE
  tab[ , is_num] <- round(tab[ , is_num], 2)
  tab[ , p_cols ] <- format_p(tab[ , p_cols ])
  # save to word
  ft <- flextable(tab)
  ft <- set_header_labels(
    ft,
    Lab = 'Lab',
    wm_corr_2_4_amp_m = "M",
    wm_corr_2_4_amp_sd = 'SD',
    wm_corr_2_4_wm_m = "M",
    wm_corr_2_4_wm_sd = 'SD',
    wm_corr_2_4_r   = "Pearson's r",
    wm_corr_2_4_p      = "p-value"

  )
  ft <- autofit(ft)
  ft <- padding(ft, padding.top = 0, padding.bottom = 0, part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- font(ft, fontname = "Arial", part = "all")
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- height(ft, height = 0.5, part = "body")
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft)
  print(doc, target = file.path(wordfiles_path, paste0(pipes_labels[pip], "_correlation_2_4_vwm.docx")))
  
} 
# loop end
