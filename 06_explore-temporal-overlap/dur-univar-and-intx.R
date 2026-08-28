
################################################################################
# Compare models of shared fert at CHR using different measures of temporal overlap
# Check for interactions between temopral overlap and phenotypes
# Heather Kenny-Duddela
# Aug 4, 2025
################################################################################

# load data
# from script adjust-table.R
dyad.ep <- read.csv("input-files/dyad_ep_complete_cases_CHRonly_with fertile.csv")

# libraries
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(bbmle) # for AICtab
library(jtools) # for plotting model predictions
library(interactions) # for plotting interaction effects

## add variables
# quadratic terms
dyad.ep$total.dur.sqr <- dyad.ep$total.dur^2
dyad.ep$dur.fertile.sqr <- dyad.ep$dur.fertile^2
dyad.ep$ci.diff.sqr <- dyad.ep$ci.diff^2

#-------------------------------------------------------------------------------
# compare models with different predictors
#-------------------------------------------------------------------------------

# Binomial----------------------------------------------------------------------

# active dur
bin.dur.active <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                            total.dur, data=dyad.ep, family=binomial)
summary(bin.dur.active)


bin.dur.active.quad <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                            total.dur + total.dur.sqr,
                            data=dyad.ep, family=binomial)
summary(bin.dur.active.quad)


# fertile dur
bin.dur.fertile <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                             dur.fertile, data=dyad.ep, family=binomial)
summary(bin.dur.fertile)


bin.dur.fertile.quad <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                             dur.fertile + dur.fertile.sqr,
                             data=dyad.ep, family=binomial)
summary(bin.dur.fertile.quad)


# ci.diff
bin.dur.cidiff <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                            ci.diff, data=dyad.ep, family=binomial)
summary(bin.dur.cidiff)


bin.dur.cidiff.quad <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                            ci.diff + ci.diff.sqr,
                            data=dyad.ep, family=binomial)
summary(bin.dur.cidiff.quad)


### compare models
bin.dur.aic <- AICtab(bin.dur.active, bin.dur.active.quad, bin.dur.fertile, 
                      bin.dur.fertile.quad, bin.dur.cidiff, bin.dur.cidiff.quad)
bin.dur.aic
# dur.fertile is the best predictor, but only slightly better than 
# active (dAIC 1.7) and ci.diff (dAIC 1.9). No support for quadratic models

write.csv(bin.dur.aic, "output-files/binary univariate dur aic table.csv")


# Count ------------------------------------------------------------------------

# active dur
count.dur.active <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                              total.dur +
                              offset(log(pat.clutch.size)),
                            data=dyad.ep, family=nbinom2)
summary(count.dur.active)


count.dur.active.quad <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                              total.dur + total.dur.sqr +
                              offset(log(pat.clutch.size)),
                            data=dyad.ep, family=nbinom2)
summary(count.dur.active.quad)


# fertile dur
count.dur.fertile <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                              dur.fertile +
                              offset(log(pat.clutch.size)),
                            data=dyad.ep, family=nbinom2)
summary(count.dur.fertile)


count.dur.fertile.quad <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                                   dur.fertile + dur.fertile.sqr +
                                   offset(log(pat.clutch.size)),
                                 data=dyad.ep, family=nbinom2)
summary(count.dur.fertile.quad)


# ci.diff
count.dur.ci.diff <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                               ci.diff +
                               offset(log(pat.clutch.size)),
                             data=dyad.ep, family=nbinom2)
summary(count.dur.ci.diff)


count.dur.ci.diff.quad <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                                    ci.diff + ci.diff.sqr +
                                    offset(log(pat.clutch.size)),
                                  data=dyad.ep, family=nbinom2)
summary(count.dur.ci.diff.quad)


### compare models

count.dur.aic <- AICtab(count.dur.active, count.dur.active.quad, 
                        count.dur.fertile, 
                      count.dur.fertile.quad, count.dur.ci.diff, 
                      count.dur.ci.diff.quad)
count.dur.aic
# dur.active and dur.fertile are essentially the same. dur.active is the top
# but dAIC with dur.fertile is only 0.6. No support for quadratic models

write.csv(count.dur.aic, "output-files/count univariate dur aic table.csv")


# Proportion -------------------------------------------------------------------

# active dur
prop.dur.active <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                             total.dur,
                           weights=pat.clutch.size, 
                           data=dyad.ep, family=binomial)
summary(prop.dur.active)


prop.dur.active.quad <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                             total.dur + total.dur.sqr,
                           weights=pat.clutch.size, 
                           data=dyad.ep, family=binomial)
summary(prop.dur.active.quad)


# fertile dur
prop.dur.fertile <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                             dur.fertile,
                           weights=pat.clutch.size, 
                           data=dyad.ep, family=binomial)
summary(prop.dur.fertile)


prop.dur.fertile.quad <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                  dur.fertile + dur.fertile.sqr,
                                weights=pat.clutch.size, 
                                data=dyad.ep, family=binomial)
summary(prop.dur.fertile.quad)


# ci.diff
prop.dur.ci.diff <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                              ci.diff,
                            weights=pat.clutch.size, 
                            data=dyad.ep, family=binomial)
summary(prop.dur.ci.diff)


prop.dur.ci.diff.quad <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                   ci.diff + ci.diff.sqr,
                                 weights=pat.clutch.size, 
                                 data=dyad.ep, family=binomial)
summary(prop.dur.ci.diff.quad)


### compare models

prop.dur.aic <- AICtab(prop.dur.active, prop.dur.active.quad, 
                        prop.dur.fertile, 
                        prop.dur.fertile.quad, prop.dur.ci.diff, 
                        prop.dur.ci.diff.quad)
prop.dur.aic
# dur.fertile and dur.fertile.quad are the same and top. Marginally better than
# dur.active (dAIC 1.8)

write.csv(prop.dur.aic, "output-files/prop univariate dur aic table.csv")


#-------------------------------------------------------------------------------
# explore dur.fertile as categorical
#-------------------------------------------------------------------------------

# add bins
range(dyad.ep$dur.fertile) # 0-11
hist(dyad.ep$dur.fertile)

### zero should have its own bin
dyad.ep$dur.fert.cat <- NA
dyad.ep$dur.fert.cat[dyad.ep$dur.fertile==0] <- "dur.fert.0"
dyad.ep$dur.fert.cat[dyad.ep$dur.fertile>0 & 
                       dyad.ep$dur.fertile<4] <- "dur.fert.1to3" # 1-3
dyad.ep$dur.fert.cat[dyad.ep$dur.fertile>3 & 
                       dyad.ep$dur.fertile<7] <- "dur.fert.4to6" # 4-6
dyad.ep$dur.fert.cat[dyad.ep$dur.fertile>6 & 
                       dyad.ep$dur.fertile<10] <- "dur.fert.7to9" # 7-9
dyad.ep$dur.fert.cat[dyad.ep$dur.fertile>9] <- "dur.fert.9up" # 10-11

## plot categories 
# binary
ggplot(dyad.ep, aes(x=dur.fert.cat, fill=bin_fert)) +
  geom_bar(position="fill")

# count
ggplot(dyad.ep, aes(x=dur.fert.cat, y=log(shared_fert))) +
  geom_boxplot()

# proportion
ggplot(dyad.ep, aes(x=dur.fert.cat, y=prop.fert)) +
  geom_boxplot()

ggplot(dyad.ep, aes(x=dur.fert.cat, y=log(prop.fert))) +
  geom_boxplot()


### try two level category of 0-3 and 4+
dyad.ep$dur.fert.cat2 <- NA
dyad.ep$dur.fert.cat2[dyad.ep$dur.fertile<4] <- "dur.fert.0to3" # 0-3
dyad.ep$dur.fert.cat2[dyad.ep$dur.fertile>=4] <- "dur.fert.4up" # 4-11

# plot these categories
# binary
ggplot(dyad.ep, aes(x=dur.fert.cat2, fill=bin_fert)) +
  geom_bar(position="fill")

ggplot(dyad.ep, aes(x=dur.fert.cat2, fill=bin_fert)) +
  geom_bar()

# count
ggplot(dyad.ep, aes(x=dur.fert.cat2, y=log(shared_fert))) +
  geom_boxplot()

# proportion
ggplot(dyad.ep, aes(x=dur.fert.cat2, y=log(prop.fert))) +
  geom_boxplot()


### try 3-level category of 0, 1-8, 9-11
dyad.ep$dur.fert.cat3 <- NA
dyad.ep$dur.fert.cat3[dyad.ep$dur.fertile==0] <- "dur.fert.0"
dyad.ep$dur.fert.cat3[dyad.ep$dur.fertile>0 & 
                       dyad.ep$dur.fertile<9] <- "dur.fert.1to8" # 1-8
dyad.ep$dur.fert.cat3[dyad.ep$dur.fertile>=9] <- "dur.fert.9up" # 9-11

# plot
# binary
ggplot(dyad.ep, aes(x=dur.fert.cat3, fill=bin_fert)) +
  geom_bar(position="fill")

# count
ggplot(dyad.ep, aes(x=dur.fert.cat3, y=log(shared_fert+0.001))) +
  geom_boxplot()

# proportion
ggplot(dyad.ep, aes(x=dur.fert.cat3, y=log(prop.fert+0.001))) +
  geom_boxplot()


### 3-level again except 10-11 instead of 9-11
dyad.ep$dur.fert.cat4 <- NA
dyad.ep$dur.fert.cat4[dyad.ep$dur.fertile==0] <- "dur.fert.0"
dyad.ep$dur.fert.cat4[dyad.ep$dur.fertile>0 & 
                        dyad.ep$dur.fertile<10] <- "dur.fert.1to9" # 1-9
dyad.ep$dur.fert.cat4[dyad.ep$dur.fertile>=10] <- "dur.fert.9up" # 10-11

# plot
# proportion
ggplot(dyad.ep, aes(x=dur.fert.cat4, y=log(prop.fert))) +
  geom_boxplot()


# Binary -----------------------------------------------------------------------

bin.dur.cat1 <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                            dur.fert.cat, 
                        data=dyad.ep, family=binomial)
summary(bin.dur.cat1)
# 4to6 is marginally significant


bin.dur.cat2 <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                          dur.fert.cat2, 
                        data=dyad.ep, family=binomial)
summary(bin.dur.cat2)
# 4up is marginally significant


bin.dur.cat3 <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                          dur.fert.cat3, 
                        data=dyad.ep, family=binomial)
summary(bin.dur.cat3)
# none are significant


bin.dur.cat4 <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                          dur.fert.cat4, 
                        data=dyad.ep, family=binomial)
summary(bin.dur.cat4)


# Count ------------------------------------------------------------------------

count.dur.cat1 <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                              dur.fert.cat +
                              offset(log(pat.clutch.size)),
                            data=dyad.ep, family=nbinom2)
summary(count.dur.cat1)
# none are significant, but 1to3 is negative while others are positive

count.dur.cat2 <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                            dur.fert.cat2 +
                            offset(log(pat.clutch.size)),
                          data=dyad.ep, family=nbinom2)
summary(count.dur.cat2)
# not significant, but estimate is positive


count.dur.cat3 <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                            dur.fert.cat3 +
                            offset(log(pat.clutch.size)),
                          data=dyad.ep, family=nbinom2)
summary(count.dur.cat3)
# none are significant


count.dur.cat4 <- glmmTMB(shared_fert ~ (1|female) + (1|male) + 
                            dur.fert.cat4 +
                            offset(log(pat.clutch.size)),
                          data=dyad.ep, family=nbinom2)
summary(count.dur.cat4)
# none are significant, but estimate for 9up is bigger


# Proportion -------------------------------------------------------------------

prop.dur.cat1 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                             dur.fert.cat,
                           weights=pat.clutch.size, 
                           data=dyad.ep, family=binomial)
summary(prop.dur.cat1)
# 9up is significant! 


prop.dur.cat2 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                           dur.fert.cat2,
                         weights=pat.clutch.size, 
                         data=dyad.ep, family=binomial)
summary(prop.dur.cat2)
# 4up is marginally significant


prop.dur.cat3 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                           dur.fert.cat3,
                         weights=pat.clutch.size, 
                         data=dyad.ep, family=binomial)
summary(prop.dur.cat3)
# none are significant


prop.dur.cat4 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                           dur.fert.cat4,
                         weights=pat.clutch.size, 
                         data=dyad.ep, family=binomial)
summary(prop.dur.cat4)
# 9up is significant (10-11 days)!

### CONCLUSIONS: 
# 2-level category is best to use with binary
# alternative 3-level is best to use with proportion
# none are very good for count, but try alternative 3-level

#-------------------------------------------------------------------------------
# explore interactions with phenotype and categorical dur
#-------------------------------------------------------------------------------

# Binomial ---------------------------------------------------------------------

# female throat - not significant
bin.dur.cat2.fthroat <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                          dur.fert.cat2*fem_throat_avg_bright, 
                        data=dyad.ep, family=binomial)
summary(bin.dur.cat2.fthroat)

# female breast - SIGNIFICANT! 
bin.dur.cat2.fbreast <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat2*fem_breast_avg_bright, 
                                data=dyad.ep, family=binomial)
summary(bin.dur.cat2.fbreast)
# the increase in probability of fertilization as female breast gets brighter
# is stronger for the 0to3 category than the 4up
interact_plot(model=bin.dur.cat2.fbreast, pred=fem_breast_avg_bright,
              modx=dur.fert.cat2,
              main.title="Binary intx model categorical dur") 
ggsave("output-files/binary dur fbreast intx plot.png", w=5, h=3)


# female tail - not significant
bin.dur.cat2.ftail <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat2*fem_TS, 
                                data=dyad.ep, family=binomial)
summary(bin.dur.cat2.ftail)


# male throat - not significant
bin.dur.cat2.mthroat <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat2*male_throat_avg_bright, 
                                data=dyad.ep, family=binomial)
summary(bin.dur.cat2.mthroat)


# male breast - not significant
bin.dur.cat2.mbreast <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat2*male_breast_avg_bright, 
                                data=dyad.ep, family=binomial)
summary(bin.dur.cat2.mbreast)


# male tail - not significant
bin.dur.cat2.mtail <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat2*male_TS, 
                                data=dyad.ep, family=binomial)
summary(bin.dur.cat2.mtail)


# socF throat - not significant
bin.dur.cat2.socFthroat <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat2*soc.fem_throat_avg_bright_scaled, 
                                data=dyad.ep, family=binomial)
summary(bin.dur.cat2.socFthroat)


# socF breast - not significant
bin.dur.cat2.socFbreast <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                     dur.fert.cat2*soc.fem_breast_avg_bright_scaled, 
                                   data=dyad.ep, family=binomial)
summary(bin.dur.cat2.socFbreast)


# socF tail - not significant
bin.dur.cat2.socFtail <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                     dur.fert.cat2*soc.fem_TS_scaled, 
                                   data=dyad.ep, family=binomial)
summary(bin.dur.cat2.socFtail)


# socM throat - not significant
bin.dur.cat2.socMthroat <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                     dur.fert.cat2*soc.male_throat_avg_bright_scaled, 
                                   data=dyad.ep, family=binomial)
summary(bin.dur.cat2.socMthroat)


# socM breast - not significant
bin.dur.cat2.socMbreast <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                     dur.fert.cat2*soc.male_breast_avg_bright_scaled, 
                                   data=dyad.ep, family=binomial)
summary(bin.dur.cat2.socMbreast)


# socM tail - not significant
bin.dur.cat2.socMtail <- glmmTMB(bin_num_fert ~ (1|female) + (1|male) +
                                     dur.fert.cat2*soc.male_TS_scaled, 
                                   data=dyad.ep, family=binomial)
summary(bin.dur.cat2.socMtail)


# Count ------------------------------------------------------------------------


# female throat - 1to9 is marginal
count.dur.cat4.fthroat <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat4*fem_throat_avg_bright + 
                                    offset(log(pat.clutch.size)), 
                                data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.fthroat)

interact_plot(model=count.dur.cat4.fthroat, pred=fem_throat_avg_bright,
              modx=dur.fert.cat4, set.offset=5,
              main.title="Count intx model categorical dur") 

# female breast - 1to9 is marginal 
count.dur.cat4.fbreast <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat4*fem_breast_avg_bright + 
                                    offset(log(pat.clutch.size)), 
                                data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.fbreast)

interact_plot(model=count.dur.cat4.fbreast, pred=fem_breast_avg_bright,
              modx=dur.fert.cat4, set.offset=5,
              main.title="Count intx model categorical dur") 



# female tail - CONVERGENCE ISSUE!
count.dur.cat4.ftail <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                dur.fert.cat4*fem_TS + 
                                  offset(log(pat.clutch.size)), 
                              data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.ftail)

# ALSO CONVERGENCE ISSUES HERE
count.dur.cat2.ftail <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat2*fem_TS + 
                                  offset(log(pat.clutch.size)), 
                                data=dyad.ep, family=nbinom2)
summary(count.dur.cat2.ftail)



# male throat - not significant
count.dur.cat4.mthroat <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat4*male_throat_avg_bright + 
                                    offset(log(pat.clutch.size)), 
                                data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.mthroat)


# male breast - not significant
count.dur.cat4.mbreast <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat4*male_breast_avg_bright + 
                                    offset(log(pat.clutch.size)), 
                                data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.mbreast)


# male tail - CONVERGENCE ISSUES
count.dur.cat4.mtail <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                dur.fert.cat4*male_TS + 
                                  offset(log(pat.clutch.size)), 
                              data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.mtail)

# with 2-level categorical, not significant
count.dur.cat2.mtail <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                  dur.fert.cat2*male_TS + 
                                  offset(log(pat.clutch.size)), 
                                data=dyad.ep, family=nbinom2)
summary(count.dur.cat2.mtail)


# socF throat - not significant
count.dur.cat4.socFthroat <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                     dur.fert.cat4*soc.fem_throat_avg_bright_scaled + 
                                       offset(log(pat.clutch.size)), 
                                   data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.socFthroat)


# socF breast - not significant
count.dur.cat4.socFbreast <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                     dur.fert.cat4*soc.fem_breast_avg_bright_scaled +
                                       offset(log(pat.clutch.size)), 
                                   data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.socFbreast)


# socF tail - not significant
count.dur.cat4.socFtail <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                   dur.fert.cat4*soc.fem_TS_scaled + 
                                     offset(log(pat.clutch.size)), 
                                 data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.socFtail)


# socM throat - not significant
count.dur.cat4.socMthroat <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                     dur.fert.cat4*soc.male_throat_avg_bright_scaled + 
                                       offset(log(pat.clutch.size)), 
                                   data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.socMthroat)


# socM breast - CONVERGENCE ISSUES
count.dur.cat4.socMbreast <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                     dur.fert.cat4*soc.male_breast_avg_bright_scaled + 
                                       offset(log(pat.clutch.size)), 
                                   data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.socMbreast)

# not significant
count.dur.cat2.socMbreast <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                       dur.fert.cat2*soc.male_breast_avg_bright_scaled + 
                                       offset(log(pat.clutch.size)), 
                                     data=dyad.ep, family=nbinom2)
summary(count.dur.cat2.socMbreast)


# socM tail - 1to9 marginal
count.dur.cat4.socMtail <- glmmTMB(shared_fert ~ (1|female) + (1|male) +
                                   dur.fert.cat4*soc.male_TS_scaled + 
                                     offset(log(pat.clutch.size)), 
                                 data=dyad.ep, family=nbinom2)
summary(count.dur.cat4.socMtail)

interact_plot(model=count.dur.cat4.socMtail, pred=soc.male_TS_scaled,
              modx=dur.fert.cat4, set.offset=5,
              main.title="Count intx model categorical dur")

# Proportion--------------------------------------------------------------------
# use dur.fert.cat4


# female throat - not significant
prop.dur.cat4.fthroat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                  dur.fert.cat4*fem_throat_avg_bright, 
                                data=dyad.ep, family=binomial, 
                                weights=pat.clutch.size)
summary(prop.dur.cat4.fthroat)

# female breast - not significant
prop.dur.cat4.fbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                  dur.fert.cat4*fem_breast_avg_bright, 
                                data=dyad.ep, family=binomial, 
                                weights=pat.clutch.size)
summary(prop.dur.cat4.fbreast)


# female tail - not significant
prop.dur.cat4.ftail <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                dur.fert.cat4*fem_TS, 
                              data=dyad.ep, family=binomial, 
                              weights=pat.clutch.size)
summary(prop.dur.cat4.ftail)


# male throat - not significant
prop.dur.cat4.mthroat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                  dur.fert.cat4*male_throat_avg_bright, 
                                data=dyad.ep, family=binomial, 
                                weights=pat.clutch.size)
summary(prop.dur.cat4.mthroat)


# male breast - not significant
prop.dur.cat4.mbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                  dur.fert.cat4*male_breast_avg_bright, 
                                data=dyad.ep, family=binomial, 
                                weights=pat.clutch.size)
summary(prop.dur.cat4.mbreast)


# male tail - not significant
prop.dur.cat4.mtail <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                dur.fert.cat4*male_TS, 
                              data=dyad.ep, family=binomial, 
                              weights=pat.clutch.size)
summary(prop.dur.cat4.mtail)


# socF throat - CONVERGENCE ISSUES
prop.dur.cat4.socFthroat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                     dur.fert.cat4*soc.fem_throat_avg_bright_scaled, 
                                   data=dyad.ep, family=binomial, 
                                   weights=pat.clutch.size)
summary(prop.dur.cat4.socFthroat)


# socF breast - marginally significant for both 1to9 and 9up
# negative effect sizes, so as soc female breast gets brighter, proportion shared
# gets smaller and the effect is stronger for higher levels of overlap
prop.dur.cat4.socFbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                     dur.fert.cat4*soc.fem_breast_avg_bright_scaled, 
                                   data=dyad.ep, family=binomial, 
                                   weights=pat.clutch.size)
summary(prop.dur.cat4.socFbreast)



# socF tail - not significant
prop.dur.cat4.socFtail <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                   dur.fert.cat4*soc.fem_TS_scaled, 
                                 data=dyad.ep, family=binomial, 
                                 weights=pat.clutch.size)
summary(prop.dur.cat4.socFtail)


# socM throat - not significant
prop.dur.cat4.socMthroat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                     dur.fert.cat4*soc.male_throat_avg_bright_scaled, 
                                   data=dyad.ep, family=binomial, 
                                   weights=pat.clutch.size)
summary(prop.dur.cat4.socMthroat)


# socM breast - not significant
prop.dur.cat4.socMbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                     dur.fert.cat4*soc.male_breast_avg_bright_scaled, 
                                   data=dyad.ep, family=binomial, 
                                   weights=pat.clutch.size)
summary(prop.dur.cat4.socMbreast)


# socM tail - marginal for 1to9
# negative main effect, as soc male tail gets longer, proportion shared decreases
# positive interaction at 1to9, so increase in tail length at this overlap is more flat
# negative interaction at 9up, so main negative effect is stronger at this overlap
prop.dur.cat4.socMtail <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                   dur.fert.cat4*soc.male_TS_scaled, 
                                 data=dyad.ep, family=binomial, 
                                 weights=pat.clutch.size)
summary(prop.dur.cat4.socMtail)

#-------------------------------------------------------------------------------

dyad.ep.real <- subset(dyad.ep, dyad.ep$bin_fert=="yes")
median(dyad.ep.real$process.dist.m)

all.site <- read.csv('input-files/dyad_ep_complete_cases_with fertile.csv')
all.site.real <- subset(all.site, all.site$bin_fert=="yes")
all.site.sum <- all.site.real %>%
  group_by(barn) %>%
  summarise(median.dist = median(process.dist.m),
            mean.dist = mean(process.dist.m))
