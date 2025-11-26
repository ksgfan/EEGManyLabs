library("TOSTER")
library("pwr")

##### Compute SESOI
?pwr.t.test
pwr.t.test(power = 0.9, n=70, sig.level=0.02, type= "paired", alternative="two.sided") # n could be also 25 (for each lab)
# The smallest effect size that we have sufficient power (90%) to detect is d = 0.44


##### Hypothesis H2.2

# Compute SESOI
pwr.r.test(n = 250, sig.level = 0.02, power = 0.9,
           alternative = c("two.sided", "less","greater"))



# Alternatively, we could compute the equivalence bounds using the 'small telescope approach'
# https://lakens.github.io/statistical_inferences/09-equivalencetest.html#specifying-the-sesoi-using-the-small-telescopes-approach
pwr::pwr.t.test(
  n = 30, 
  sig.level = 0.02, 
  power = 0.33333, 
  type = "one.sample",
  alternative = "two.sided"
)
# Given the sample size from the original study n = 30, the effect size that gives us 33% power is d = 0.36 - similar to SESOI

# small telescope approach for correlations
pwr::pwr.r.test(
  n = 30, 
  sig.level = 0.02, 
  power = 0.33333, 
  alternative = "greater"
)




