library(metafor)

################################################################################
################################################################################
# do the funnel-plot diagnostic without EEGManyLabs studies
################################################################################
################################################################################

# correlations and samples
r   = c(0.61, 0.48, 0.72, 0.48, 0.62, 0.65, 0.63, 0.52, 0.59, 0.40, 0.80, 0.78, # Luria, 2017, from Fig2.A. 
        0.33, 0.27, 0.24) # Unsworth 2015, Adam 2018, Tröndle, 2024

n   = c(14, 33, 18, 25, 24, 30, 18, 39, 23, 35, 25, 36,
        170, 72, 55)

group = c(
  rep("Luria (2017)", 12),
  rep("Unsworth (2015)", 1),
  rep("Adam (2018)", 1),
  rep("Tröndle (2024)", 1)
)

cols = c(
  "Luria (2017)" = "blue",
  "Unsworth (2015)"   = "#006400", 
  "Adam (2018)"       = "#228B22", 
  "Tröndle (2024)"    = "#32CD32" 
)

# Convert correlations to Fisher's z
dat = escalc(measure = "ZCOR", ri = r, ni = n)

# Random-effects meta-analysis
res_without = rma(yi, vi, data = dat)

# Funnel plot
cex_vec = 0.05 * n
funnel(res_without, 
       xlab = "Fisher's z",
       main = "Funnel Plot", 
       pch = 19, 
       cex = cex_vec,
       col = cols[group])
legend("topright", legend = names(cols), col = cols, pch = 1, bty = "n")


rtest_without = regtest(res_without, model = "rma")
taf_without = trimfill(res_without)

# funnel with missing studies
funnel(taf_without, 
       xlab = "Fisher's z",
       main = "Funnel Plot", 
       pch = 19, 
       cex = cex_vec,
       col = cols[group])
legend("topright", legend = names(cols), col = cols, pch = 1, bty = "n")

# print
res_without
res_without$pval
rtest_without
taf_without
taf_without$pval

################################################################################
################################################################################
# do the funnel-plot diagnostic WITH EEGManyLabs studies
################################################################################
################################################################################

# correlations and samples
r   = c(0.61, 0.48, 0.72, 0.48, 0.62, 0.65, 0.63, 0.52, 0.59, 0.40, 0.80, 0.78, # Luria, 2017, from Fig2.A. 
         0.38, -0.30, 0.42, 0.13,0.23,0.32,0.26,0.38, -0.05, 0.05, # EEGManyLabs
         0.33, 0.27, 0.24) # Unsworth 2015, Adam 2018, Tröndle, 2024

n   = c(14, 33, 18, 25, 24, 30, 18, 39, 23, 35, 25, 36,
         23, 26,11,26,22,20,24,15,30,20,
         170, 72, 55)

group = c(
  rep("Luria (2017)", 12),
  rep("EEGManyLabs", 10),
  rep("Unsworth (2015)", 1),
  rep("Adam (2018)", 1),
  rep("Tröndle (2024)", 1)
)

cols = c(
  "Luria (2017)" = "blue",
  "EEGManyLabs" = "red",
  "Unsworth (2015)"   = "#006400", 
  "Adam (2018)"       = "#228B22", 
  "Tröndle (2024)"    = "#32CD32" 
)

# Convert correlations to Fisher's z
dat = escalc(measure = "ZCOR", ri = r, ni = n)

# Random-effects meta-analysis
res = rma(yi, vi, data = dat)

# Funnel plot
cex_vec = 0.05 * n
funnel(res, 
       xlab = "Fisher's z",
       main = "Funnel Plot", 
       pch = 19, 
       cex = cex_vec,
       col = cols[group])
legend("topright", legend = names(cols), col = cols, pch = 1, bty = "n")

rtest = regtest(res, model = "rma")
taf = trimfill(res)

# funnel with missing studies
funnel(taf, 
       xlab = "Fisher's z",
       main = "Funnel Plot", 
       pch = 19, 
       cex = cex_vec,
       col = cols[group])
legend("topright", legend = names(cols), col = cols, pch = 1, bty = "n")

# print
res
res$pval
rtest
taf
taf$pval




