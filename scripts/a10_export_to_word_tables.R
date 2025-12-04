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
format_p = function(x) {
  out = ifelse(
    x < 0.001,
    formatC(x, format = "e", digits = 2),   # scientific for very small p
    formatC(x, format = "f", digits = 3)    # fixed, 3 decimals
  )
  out
}

################################################################################
# results tables
################################################################################

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
  is_num = sapply(tab, is.numeric)
  p_cols = grepl("_p$", names(tab))
  is_num[p_cols] = FALSE
  tab[ , is_num] = round(tab[ , is_num], 2)
  tab[ , p_cols ] = format_p(tab[ , p_cols ])
  # save to word
  ft = flextable(tab)
  ft = set_header_labels(
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
  ft = autofit(ft)
  ft = padding(ft, padding.top = 0, padding.bottom = 0, part = "all")
  ft = valign(ft, valign = "center", part = "all")
  ft = font(ft, fontname = "Arial", part = "all")
  ft = fontsize(ft, size = 8, part = "all")
  ft = height(ft, height = 0.5, part = "body")
  doc = read_docx()
  doc = body_add_flextable(doc, ft)
  print(doc, target = file.path(wordfiles_path, paste0(pipes_labels[pip], "_outcome_neutral.docx")))


  ##############################################################################
  # set size 2 vs 4
  ##############################################################################
  tab = dat[ , c("Lab", "setsize2_m", "setsize2_sd", "setsize4_m", "setsize4_sd",
                 "ph_2_4_t_stat", "ph_2_4_df", "ph_2_4_p", "ph_2_4_gz") ]
  # round all values except p-vals
  is_num = sapply(tab, is.numeric)
  p_cols = grepl("_p$", names(tab))
  is_num[p_cols] = FALSE
  tab[ , is_num] = round(tab[ , is_num], 2)
  tab[ , p_cols ] = format_p(tab[ , p_cols ])
  # save to word
  ft = flextable(tab)
  ft = set_header_labels(
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
  ft = autofit(ft)
  ft = padding(ft, padding.top = 0, padding.bottom = 0, part = "all")
  ft = valign(ft, valign = "center", part = "all")
  ft = font(ft, fontname = "Arial", part = "all")
  ft = fontsize(ft, size = 8, part = "all")
  ft = height(ft, height = 0.5, part = "body")
  doc = read_docx()
  doc = body_add_flextable(doc, ft)
  print(doc, target = file.path(wordfiles_path, paste0(pipes_labels[pip], "_t_test_2_4.docx")))


  ##############################################################################
  # set size 2 vs 6
  ##############################################################################
  tab = dat[ , c("Lab", "setsize2_m", "setsize2_sd", "setsize6_m", "setsize6_sd",
                 "ph_2_6_t_stat", "ph_2_6_df", "ph_2_6_p", "ph_2_6_gz") ]
  # round all values except p-vals
  is_num = sapply(tab, is.numeric)
  p_cols = grepl("_p$", names(tab))
  is_num[p_cols] = FALSE
  tab[ , is_num] = round(tab[ , is_num], 2)
  tab[ , p_cols ] = format_p(tab[ , p_cols ])
  # save to word
  ft = flextable(tab)
  ft = set_header_labels(
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
  ft = autofit(ft)
  ft = padding(ft, padding.top = 0, padding.bottom = 0, part = "all")
  ft = valign(ft, valign = "center", part = "all")
  ft = font(ft, fontname = "Arial", part = "all")
  ft = fontsize(ft, size = 8, part = "all")
  ft = height(ft, height = 0.5, part = "body")
  doc = read_docx()
  doc = body_add_flextable(doc, ft)
  print(doc, target = file.path(wordfiles_path, paste0(pipes_labels[pip], "_t_test_2_6.docx")))


  ##############################################################################
  # correlation of amplitude decrease (2 to 4) and VWM capacity
  ##############################################################################
  tab = dat[ , c("Lab", "wm_corr_2_4_amp_m", "wm_corr_2_4_amp_sd", "wm_corr_2_4_wm_m", "wm_corr_2_4_wm_sd",
                 "wm_corr_2_4_r", "wm_corr_2_4_p") ]

  # make the correlation positive
  tab$wm_corr_2_4_r = tab$wm_corr_2_4_r*(-1)

  # round all values except p-vals
  is_num = sapply(tab, is.numeric)
  p_cols = grepl("_p$", names(tab))
  is_num[p_cols] = FALSE
  tab[ , is_num] = round(tab[ , is_num], 2)
  tab[ , p_cols ] = format_p(tab[ , p_cols ])
  # save to word
  ft = flextable(tab)
  ft = set_header_labels(
    ft,
    Lab = 'Lab',
    wm_corr_2_4_amp_m = "M",
    wm_corr_2_4_amp_sd = 'SD',
    wm_corr_2_4_wm_m = "M",
    wm_corr_2_4_wm_sd = 'SD',
    wm_corr_2_4_r   = "Pearson's r",
    wm_corr_2_4_p      = "p-value"

  )
  ft = autofit(ft)
  ft = padding(ft, padding.top = 0, padding.bottom = 0, part = "all")
  ft = valign(ft, valign = "center", part = "all")
  ft = font(ft, fontname = "Arial", part = "all")
  ft = fontsize(ft, size = 8, part = "all")
  ft = height(ft, height = 0.5, part = "body")
  doc = read_docx()
  doc = body_add_flextable(doc, ft)
  print(doc, target = file.path(wordfiles_path, paste0(pipes_labels[pip], "_correlation_2_4_vwm.docx")))

}
# loop end
################################################################################
# results tables end
################################################################################


################################################################################
# rejected trials
################################################################################

all_dat = list()
pipes_labels = c("Direct", "Advanced", "ICA")

# loop over pipelines and save results
for (pip in 1 : length(pipes_labels)){
  
  # load data
  dat = read.csv2(file.path(project_path, 'data/csv_files', paste0('excl_tab_', pipes_labels[pip], '.csv')), header = T, sep = ",")
  # names(dat)
  
  # concat
  dat$Pipeline = pipes_labels[pip]
  all_dat[[pip]] = dat
  
} # loop end

all_dat_bind = bind_rows(all_dat)
colnames(all_dat_bind)
all_dat_bind$final_N = all_dat_bind$N - all_dat_bind$all_excluded_subs_bad_trials_and_performance

table1 = all_dat_bind %>%
  arrange(Lab, factor(Pipeline, levels = pipes_labels)) %>%
  transmute(
    Lab,
    Pipeline,
    Amplitude_M   = avg_excl_AmpThresh,
    Amplitude_SD  = sd_excl_AmpThresh,
    Blocking_M    = avg_excl_blocking,
    Blocking_SD   = sd_excl_blocking,
    VEOG_M        = avg_excl_VEOG,
    VEOG_SD       = sd_excl_VEOG,
    HEOG_M        = avg_excl_HEOG,
    HEOG_SD       = sd_excl_HEOG,
    ET_Blink_M    = avg_excl_ET_Blink,
    ET_Blink_SD   = sd_excl_ET_Blink,
    ET_Saccade_M  = avg_excl_ET_Sacc,
    ET_Saccade_SD = sd_excl_ET_Sacc,
    Final_N      = final_N
  )

# round all values except p-vals
is_num = sapply(table1, is.numeric)
p_cols = grepl("_p$", names(table1))
is_num[p_cols] = FALSE
table1[ , is_num] = round(table1[ , is_num], 2)
table1[ , p_cols ] = format_p(table1[ , p_cols ])
# save to word
ft = flextable(table1)

ft = autofit(ft)
ft = padding(ft, padding.top = 0, padding.bottom = 0, part = "all")
ft = valign(ft, valign = "center", part = "all")
ft = font(ft, fontname = "Arial", part = "all")
ft = fontsize(ft, size = 8, part = "all")
ft = height(ft, height = 0.5, part = "body")
ft = width(ft, width = 0.5)   
sect_landscape = prop_section(
  page_size = page_size(orient = "landscape")
)
doc = read_docx()
doc = body_add_flextable(doc, ft)
print(doc, target = file.path(wordfiles_path, "rejected_trials.docx"))

################################################################################
# rejected trials end
################################################################################