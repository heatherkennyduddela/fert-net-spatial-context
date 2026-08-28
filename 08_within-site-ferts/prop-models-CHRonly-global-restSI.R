
################################################################################
# Proportion models of CO Fert Net 2022, CHR only
# Heather Kenny-Duddela
# June 3, 2025
# updated May 28, 2026
################################################################################

# load data

# from script adjust-table-CHRonly.R
dyad.ep2 <- read.csv("input-files/dyad_ep_complete_cases_CHRonly_with fertile and rest SI.csv")

# load libraries
library(ggplot2)
library(dplyr)
library(bbmle) # for AICtab
library(glmmTMB) # also for fitting mixed effect glm
library(DHARMa) # model diagnostics for mixed effect models
library(jtools)
library(ggpubr) # for multi-panel plots


#-------------------------------------------------------------------------------
# Fit models for EP fert only
#-------------------------------------------------------------------------------

# null
prop.null.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male),
                       data=dyad.ep2, family=binomial, weights=pat.clutch.size)

# dist
prop.dist.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                          dist.scale,
                       data=dyad.ep2, family=binomial, weights=pat.clutch.size)
summary(prop.dist.ep)


# dur 
prop.dur.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                         dur.scale,
                      data=dyad.ep2, family=binomial, weights=pat.clutch.size)
summary(prop.dur.ep)


# rest.SI
prop.restSI.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                         rest.SI.scale,
                       data=dyad.ep2, family=binomial, weights=pat.clutch.size)
summary(prop.restSI.ep)


# female tail
prop.f.tail.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                            fem_TS_scaled,
                         data=dyad.ep2, family=binomial, weights=pat.clutch.size)
summary(prop.f.tail.ep)

# male tail 
prop.m.tail.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                            male_TS_scaled,
                         data=dyad.ep2, family=binomial, weights=pat.clutch.size)
summary(prop.m.tail.ep)

# female throat 
prop.f.throat.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                             fem_throat_avg_bright_scaled,
                           data=dyad.ep2, family=binomial, weights=pat.clutch.size)
summary(prop.f.throat.ep)

# male throat
prop.m.throat.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                             male_throat_avg_bright_scaled,
                           data=dyad.ep2, family=binomial, 
                           weights=pat.clutch.size)
summary(prop.m.throat.ep)

# female breast
prop.f.breast.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                             fem_breast_avg_bright_scaled,
                           data=dyad.ep2, family=binomial, 
                           weights=pat.clutch.size)
summary(prop.f.breast.ep)

# male breast
prop.m.breast.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                             male_breast_avg_bright_scaled,
                           data=dyad.ep2, family=binomial, 
                           weights=pat.clutch.size)
summary(prop.m.breast.ep)


# attempt number 
prop.attempt.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                            attempt2,
                          data=dyad.ep2, family=binomial, 
                          weights=pat.clutch.size)
summary(prop.attempt.ep)


## Social mate traits-----------------------------------------------------------

# social male tail
prop.socM.tail.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                              soc.male_TS_scaled,
                            data=dyad.ep2, family=binomial, 
                            weights=pat.clutch.size)
summary(prop.socM.tail.ep)

# social female tail
prop.socF.tail.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                              soc.fem_TS_scaled,
                            data=dyad.ep2, family=binomial, 
                            weights=pat.clutch.size)
summary(prop.socF.tail.ep)

# social male throat
prop.socM.throat.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                                soc.male_throat_avg_bright_scaled,
                              data=dyad.ep2, family=binomial, 
                              weights=pat.clutch.size)
summary(prop.socM.throat.ep)

# social female throat
prop.socF.throat.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                                soc.fem_throat_avg_bright_scaled,
                              data=dyad.ep2, family=binomial, 
                              weights=pat.clutch.size)
summary(prop.socF.throat.ep)


# social male breast
prop.socM.breast.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                                soc.male_breast_avg_bright_scaled,
                              data=dyad.ep2, family=binomial, 
                              weights=pat.clutch.size)
summary(prop.socM.breast.ep)

# social female breast 
prop.socF.breast.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                                soc.fem_breast_avg_bright_scaled,
                              data=dyad.ep2, family=binomial, 
                              weights=pat.clutch.size)
summary(prop.socF.breast.ep)

## Model comparison table
prop.univar.EP.compare <- as.data.frame(AICtab(prop.null.ep, prop.dist.ep, 
                                              prop.dur.ep,
                                              prop.restSI.ep,
                                              prop.f.tail.ep, 
                                              prop.m.tail.ep, 
                                              prop.f.throat.ep, 
                                              prop.m.throat.ep, 
                                              prop.f.breast.ep, 
                                              prop.m.breast.ep, 
                                              prop.attempt.ep,
                                              prop.socM.tail.ep, 
                                              prop.socF.tail.ep, 
                                              prop.socM.throat.ep, 
                                              prop.socF.throat.ep,
                                              prop.socM.breast.ep, 
                                              prop.socF.breast.ep))

prop.univar.EP.compare
# better than null is socF.throat, socF.breast, m.breast, f.breast, m.tail
# socM.tail, socM.breast, attempt,dist, dur

write.csv(prop.univar.EP.compare, "output-files/prop univar aic table EPfert.csv", 
          row.names=T)

#-------------------------------------------------------------------------------
# Test for interaction between rest.SI and attempt2 (not significant)
#-------------------------------------------------------------------------------

prop.ep.restSI.intx.attempt <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                         rest.SI.scale * attempt2,
                                       data=dyad.ep2, family=binomial,
                                       weights=pat.clutch.size)
summary(prop.ep.restSI.intx.attempt)
# Family: binomial  ( logit )
# Formula:          prop.fert ~ (1 | female) + (1 | male) + rest.SI.scale * attempt2
# Data: dyad.ep2
# Weights: pat.clutch.size
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 286.9     315.2    -137.5     274.9       813 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# female (Intercept) 0.9139   0.956   
# male   (Intercept) 2.0398   1.428   
# Number of obs: 819, groups:  female, 23; male, 22
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                  -5.49884    0.61309  -8.969   <2e-16 ***
#   rest.SI.scale                -0.06932    0.38016  -0.182    0.855    
# attempt2second               -0.73672    0.51036  -1.444    0.149    
# rest.SI.scale:attempt2second  0.05200    0.61411   0.085    0.933    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#-------------------------------------------------------------------------------
# fit global model for EPfert
#-------------------------------------------------------------------------------

# Method for global model is to include all varibles, except those that are 
# highly correlated (see below). Prioritized including terms that were
# significant in univariate models. Compared AIC of these two and they were 
# not meaningfully different. Proceeded with model that has more balance between
# female and male traits. 

# NOTE: due to correlations, male TS and soc.fem breast should not be included
# in the same model. Same goes for female breast and soc.male TS

# model with EVERYTING for reference
# prop.ep.global <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
#                              dist.scale + 
#                              dur.scale +
#                              rest.SI.scale +
#                              attempt2 +
#                              soc.fem_throat_avg_bright_scaled +
#                              soc.fem_breast_avg_bright_scaled +
#                              soc.fem_TS_scaled +
#                              soc.male_throat_avg_bright_scaled + 
#                              soc.male_breast_avg_bright_scaled + 
#                              soc.male_TS_scaled +
#                              male_throat_avg_bright_scaled +
#                              male_breast_avg_bright_scaled +
#                              male_TS_scaled +
#                              fem_throat_avg_bright_scaled +
#                              fem_breast_avg_bright_scaled +
#                              fem_TS_scaled,
#                            data=dyad.ep2, family=binomial, 
#                            weights=pat.clutch.size)

# first keep male TS, exclude socF_breast
# keep female_breast, exclude socM_TS

# try specifying starting values from the logit model fit
prop.ep.global1 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                             dist.scale +
                             dur.scale +
                             rest.SI.scale +
                             attempt2 +
                             soc.fem_throat_avg_bright_scaled +
                             soc.fem_TS_scaled +
                             soc.male_throat_avg_bright_scaled +
                             soc.male_breast_avg_bright_scaled +
                             male_throat_avg_bright_scaled +
                             male_breast_avg_bright_scaled +
                             male_TS_scaled +
                             fem_throat_avg_bright_scaled +
                             fem_breast_avg_bright_scaled +
                             fem_TS_scaled,
                           data=dyad.ep2, family=binomial,
                           weights=pat.clutch.size)
summary(prop.ep.global1)


# next keep male_TS, exclude socF_breast
# keep socM_TS, exclude female_breast
prop.ep.global2 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                             dist.scale +
                             dur.scale +
                             rest.SI.scale +
                             attempt2 +
                             soc.fem_throat_avg_bright_scaled +
                             soc.fem_TS_scaled +
                             soc.male_throat_avg_bright_scaled +
                             soc.male_breast_avg_bright_scaled +
                             soc.male_TS_scaled +
                             male_throat_avg_bright_scaled +
                             male_breast_avg_bright_scaled +
                             male_TS_scaled +
                             fem_throat_avg_bright_scaled +
                             fem_TS_scaled,
                           data=dyad.ep2, family=binomial,
                           weights=pat.clutch.size)
summary(prop.ep.global2)

# keep socF_breast, exclude male_TS
# keep fem_breast, exclude socM_TS
prop.ep.global3 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                             dist.scale +
                             dur.scale +
                             attempt2 +
                             rest.SI.scale +
                             soc.fem_throat_avg_bright_scaled +
                             soc.fem_breast_avg_bright_scaled +
                             soc.fem_TS_scaled +
                             soc.male_throat_avg_bright_scaled +
                             soc.male_breast_avg_bright_scaled +
                             male_throat_avg_bright_scaled +
                             male_breast_avg_bright_scaled +
                             fem_throat_avg_bright_scaled +
                             fem_breast_avg_bright_scaled +
                             fem_TS_scaled,
                           data=dyad.ep2, family=binomial,
                           weights=pat.clutch.size)
summary(prop.ep.global3)

# keep socF_breast, exclude male_TS
# keep socM_TS, exclude fem_breast
prop.ep.global4 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                             dist.scale +
                             dur.scale +
                             attempt2 +
                             rest.SI.scale +
                             soc.fem_throat_avg_bright_scaled +
                             soc.fem_breast_avg_bright_scaled +
                             soc.fem_TS_scaled +
                             soc.male_throat_avg_bright_scaled +
                             soc.male_breast_avg_bright_scaled +
                             soc.male_TS_scaled +
                             male_throat_avg_bright_scaled +
                             male_breast_avg_bright_scaled +
                             fem_throat_avg_bright_scaled +
                             fem_TS_scaled,
                           data=dyad.ep2, family=binomial,
                           weights=pat.clutch.size)
summary(prop.ep.global4)

## Should use model 1 or 3 because they include fem_breast which is significant
# global 3 is very slightly better, but global1 includes male_TS and seems more balanced
# proceed with global1
AICtab(prop.ep.global1, prop.ep.global3)


## save model results for global1
prop.ep.global1.output <- summary(prop.ep.global1)
prop.ep.global1.coeff <- prop.ep.global1.output$coefficients$cond
prop.ep.global1.confint <- confint(prop.ep.global1)
prop.ep.global1.details <- cbind(
  prop.ep.global1.coeff, 
  prop.ep.global1.confint[1:15, 1:2],
  ngrps = c(prop.ep.global1.output$ngrps$cond, rep(NA, 13)),
  grp = c("female","male", rep(NA,13)),
  N = c(prop.ep.global1.output$nobs, rep(NA, 14)))

prop.ep.global1.details <- as.data.frame(prop.ep.global1.details)

# add back-transformed coefficeints and CIs by exponentiating
prop.ep.global1.details$bt.estimate <- exp(
  as.numeric(prop.ep.global1.details$Estimate))

prop.ep.global1.details$bt.lower.ci <- exp(
  as.numeric(prop.ep.global1.details$`2.5 %`))

prop.ep.global1.details$bt.upper.ci <- exp(
  as.numeric(prop.ep.global1.details$`97.5 %`))

write.csv(prop.ep.global1.details, "output-files/full global prop epFert mod.csv")





## Model diagnostics------------------------------------------------------------

# calculate simulated residuals
prop.ep.global1.simResids <- simulateResiduals(prop.ep.global1)

# no obvious issues
plotQQunif(prop.ep.global1.simResids) 
plotResiduals(prop.ep.global1.simResids)


# residuals against dist
plotResiduals(prop.ep.global1.simResids, form=dyad.ep2$dist.scale)

# residuals against dur
plotResiduals(prop.ep.global1.simResids, form=dyad.ep2$dur.scale)

# residuals against rest.SI
plotResiduals(prop.ep.global1.simResids, form=dyad.ep2$rest.SI.scale)

# against soc.fem_throat
plotResiduals(prop.ep.global1.simResids, 
              form=dyad.ep2$soc.fem_throat_avg_bright_scaled)

# against socM_breast
plotResiduals(prop.ep.global1.simResids, 
              form=dyad.ep2$soc.male_breast_avg_bright_scaled)

## try plotting residuals by grouping variables

# female
plotResiduals(prop.ep.global1.simResids, form=dyad.ep2$female)

# male 
plotResiduals(prop.ep.global1.simResids, form=dyad.ep2$male)

#-------------------------------------------------------------------------------
# Test for interactions between significant pheno, dist, and dur variables
#-------------------------------------------------------------------------------

# Significant terms in additive model are: 
# dist
# soc.fem_throat
# soc.male_breast
# male_breast
# fem_breast

## fit models to test for interactions between phenotype and dist

## soc.fem_throat - significant!------------------------------------------------
prop.intx.dist.socFthroat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                       dist.scale +
                                       dur.scale +
                                       rest.SI.scale +
                                       attempt2 +
                                       soc.fem_throat_avg_bright_scaled +
                                       soc.fem_TS_scaled +
                                       soc.male_throat_avg_bright_scaled +
                                       soc.male_breast_avg_bright_scaled +
                                       male_throat_avg_bright_scaled +
                                       male_breast_avg_bright_scaled +
                                       male_TS_scaled +
                                       fem_throat_avg_bright_scaled +
                                       fem_breast_avg_bright_scaled +
                                       fem_TS_scaled +
                                       dist.scale:soc.fem_throat_avg_bright_scaled,
                                     data=dyad.ep2, family=binomial,
                                     weights=pat.clutch.size)
summary(prop.intx.dist.socFthroat)

# save model results
prop.intx.dist.socFthroat.output <- summary(prop.intx.dist.socFthroat)
prop.intx.dist.socFthroat.coeff <- prop.intx.dist.socFthroat.output$coefficients$cond
prop.intx.dist.socFthroat.confint <- confint(prop.intx.dist.socFthroat)
prop.intx.dist.socFthroat.details <- cbind(
  prop.intx.dist.socFthroat.coeff, 
  prop.intx.dist.socFthroat.confint[1:16, 1:2],
  ngrps = c(prop.intx.dist.socFthroat.output$ngrps$cond, rep(NA, 14)),
  grp = c("female","male", rep(NA,14)),
  N = c(prop.intx.dist.socFthroat.output$nobs, rep(NA, 15)))

prop.intx.dist.socFthroat.details <- as.data.frame(prop.intx.dist.socFthroat.details)

# add back-transformed coefficeints and CIs by exponentiating
prop.intx.dist.socFthroat.details$bt.estimate <- exp(
  as.numeric(prop.intx.dist.socFthroat.details$Estimate))

prop.intx.dist.socFthroat.details$bt.lower.ci <- exp(
  as.numeric(prop.intx.dist.socFthroat.details$`2.5 %`))

prop.intx.dist.socFthroat.details$bt.upper.ci <- exp(
  as.numeric(prop.intx.dist.socFthroat.details$`97.5 %`))

str(prop.intx.dist.socFthroat.details)

write.csv(prop.intx.dist.socFthroat.details, 
          "output-files/prop epFert intx dist socFthroat mod.csv")


## soc.male_breast - not significant--------------------------------------------
prop.intx.dist.socMbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                       dist.scale +
                                       dur.scale +
                                       rest.SI.scale +
                                       attempt2 +
                                       soc.fem_throat_avg_bright_scaled +
                                       soc.fem_TS_scaled +
                                       soc.male_throat_avg_bright_scaled +
                                       soc.male_breast_avg_bright_scaled +
                                       male_throat_avg_bright_scaled +
                                       male_breast_avg_bright_scaled +
                                       male_TS_scaled +
                                       fem_throat_avg_bright_scaled +
                                       fem_breast_avg_bright_scaled +
                                       fem_TS_scaled +
                                       dist.scale:soc.male_breast_avg_bright_scaled,
                                     data=dyad.ep2, family=binomial,
                                     weights=pat.clutch.size)
summary(prop.intx.dist.socMbreast)

# save model results
prop.intx.dist.socMbreast.output <- summary(prop.intx.dist.socMbreast)
prop.intx.dist.socMbreast.coeff <- prop.intx.dist.socMbreast.output$coefficients$cond
prop.intx.dist.socMbreast.confint <- confint(prop.intx.dist.socMbreast)
prop.intx.dist.socMbreast.details <- cbind(
  prop.intx.dist.socMbreast.coeff, 
  prop.intx.dist.socMbreast.confint[1:16, 1:2],
  ngrps = c(prop.intx.dist.socMbreast.output$ngrps$cond, rep(NA, 14)),
  grp = c("female","male", rep(NA,14)),
  N = c(prop.intx.dist.socMbreast.output$nobs, rep(NA, 15)))


prop.intx.dist.socMbreast.details <- as.data.frame(prop.intx.dist.socMbreast.details)

# add back-transformed coefficeints and CIs by exponentiating
prop.intx.dist.socMbreast.details$bt.estimate <- exp(
  as.numeric(prop.intx.dist.socMbreast.details$Estimate))

prop.intx.dist.socMbreast.details$bt.lower.ci <- exp(
  as.numeric(prop.intx.dist.socMbreast.details$`2.5 %`))

prop.intx.dist.socMbreast.details$bt.upper.ci <- exp(
  as.numeric(prop.intx.dist.socMbreast.details$`97.5 %`))

str(prop.intx.dist.socMbreast.details)

write.csv(prop.intx.dist.socMbreast.details, 
          "output-files/prop epFert intx dist socMbreast mod.csv")


## male breast - not significant -----------------------------------------------
prop.intx.dist.Mbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                    dist.scale +
                                    dur.scale +
                                    rest.SI.scale +
                                    attempt2 +
                                    soc.fem_throat_avg_bright_scaled +
                                    soc.fem_TS_scaled +
                                    soc.male_throat_avg_bright_scaled +
                                    soc.male_breast_avg_bright_scaled +
                                    male_throat_avg_bright_scaled +
                                    male_breast_avg_bright_scaled +
                                    male_TS_scaled +
                                    fem_throat_avg_bright_scaled +
                                    fem_breast_avg_bright_scaled +
                                    fem_TS_scaled +
                                    dist.scale:male_breast_avg_bright_scaled,
                                     data=dyad.ep2, family=binomial,
                                     weights=pat.clutch.size)
summary(prop.intx.dist.Mbreast)

# save model results
prop.intx.dist.Mbreast.output <- summary(prop.intx.dist.Mbreast)
prop.intx.dist.Mbreast.coeff <- prop.intx.dist.Mbreast.output$coefficients$cond
prop.intx.dist.Mbreast.confint <- confint(prop.intx.dist.Mbreast)
prop.intx.dist.Mbreast.details <- cbind(
  prop.intx.dist.Mbreast.coeff, 
  prop.intx.dist.Mbreast.confint[1:16, 1:2],
  ngrps = c(prop.intx.dist.Mbreast.output$ngrps$cond, rep(NA, 14)),
  grp = c("female","male", rep(NA,14)),
  N = c(prop.intx.dist.Mbreast.output$nobs, rep(NA, 15)))

prop.intx.dist.Mbreast.details <- as.data.frame(prop.intx.dist.Mbreast.details)

# add back-transformed coefficeints and CIs by exponentiating
prop.intx.dist.Mbreast.details$bt.estimate <- exp(
  as.numeric(prop.intx.dist.Mbreast.details$Estimate))

prop.intx.dist.Mbreast.details$bt.lower.ci <- exp(
  as.numeric(prop.intx.dist.Mbreast.details$`2.5 %`))

prop.intx.dist.Mbreast.details$bt.upper.ci <- exp(
  as.numeric(prop.intx.dist.Mbreast.details$`97.5 %`))


write.csv(prop.intx.dist.Mbreast.details, 
          "output-files/prop epFert intx dist Mbreast mod.csv")


## female breast - not significant -----------------------------------------------
prop.intx.dist.Fbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                    dist.scale +
                                    dur.scale +
                                    rest.SI.scale +
                                    attempt2 +
                                    soc.fem_throat_avg_bright_scaled +
                                    soc.fem_TS_scaled +
                                    soc.male_throat_avg_bright_scaled +
                                    soc.male_breast_avg_bright_scaled +
                                    male_throat_avg_bright_scaled +
                                    male_breast_avg_bright_scaled +
                                    male_TS_scaled +
                                    fem_throat_avg_bright_scaled +
                                    fem_breast_avg_bright_scaled +
                                    fem_TS_scaled +
                                    dist.scale:fem_breast_avg_bright_scaled,
                                  data=dyad.ep2, family=binomial,
                                  weights=pat.clutch.size)
summary(prop.intx.dist.Fbreast)

# save model results
prop.intx.dist.Fbreast.output <- summary(prop.intx.dist.Fbreast)
prop.intx.dist.Fbreast.coeff <- prop.intx.dist.Fbreast.output$coefficients$cond
prop.intx.dist.Fbreast.confint <- confint(prop.intx.dist.Fbreast)
prop.intx.dist.Fbreast.details <- cbind(
  prop.intx.dist.Fbreast.coeff, 
  prop.intx.dist.Fbreast.confint[1:16, 1:2],
  ngrps = c(prop.intx.dist.Fbreast.output$ngrps$cond, rep(NA, 14)),
  grp = c("female","male", rep(NA,14)),
  N = c(prop.intx.dist.Fbreast.output$nobs, rep(NA, 15)))

prop.intx.dist.Fbreast.details <- as.data.frame(prop.intx.dist.Fbreast.details)

# add back-transformed coefficeints and CIs by exponentiating
prop.intx.dist.Fbreast.details$bt.estimate <- exp(
  as.numeric(prop.intx.dist.Fbreast.details$Estimate))

prop.intx.dist.Fbreast.details$bt.lower.ci <- exp(
  as.numeric(prop.intx.dist.Fbreast.details$`2.5 %`))

prop.intx.dist.Fbreast.details$bt.upper.ci <- exp(
  as.numeric(prop.intx.dist.Fbreast.details$`97.5 %`))

write.csv(prop.intx.dist.Fbreast.details, 
          "output-files/prop epFert intx dist Fbreast mod.csv")

#-------------------------------------------------------------------------------
# visualize results
#-------------------------------------------------------------------------------

# Significant terms in additive model are: 
# dist
# soc.fem_throat
# soc.male_breast
# male_breast
# fem_breast

# Significant intx terms are:
# dist*soc.fem_throat

# distance---------------------------------------------------------------------- 
# Manually generate newdata for plotting
newdata.prop.epFert.dist <- custom.newdata(dyad.ep2, "epFert", "dist.scale",
                                                   "dist", soc.pair=NA)
newdata.prop.epFert.dist$pat.clutch.size <- 4

# predictions and CIs
fit.prop.epFert.dist <- calculate.glmm.ci(prop.ep.global1, 
                                                  newdata.prop.epFert.dist)

ggplot(fit.prop.epFert.dist, aes(x=dist, y=fit)) +
  geom_point(data=dyad.ep2, 
             aes(x=process.dist.m, y=prop.fert,size=pat.clutch.size), 
             alpha=0.5) +
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="red", alpha=0.3) +
  geom_line(color="red", linewidth=1) +
  ylab("Proportion of shared fert in clutch") +
  xlab("Distance between nests in meters")

# without raw data
dist.pred <- ggplot(fit.prop.epFert.dist, aes(x=dist, y=fit)) +
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="black", alpha=0.3) +
  geom_line(color="black", linewidth=1) +
  ylab("Proportion of clutch with\nshared fertilizations") +
  xlab("Distance between nests in meters") +
  theme_light()

dist.pred

# ggsave("output-files/prop mod dist plot.png", h=3, w=5)


# soc Fem throat intx distance -------------------------------------------------

# new data for SDplus
newdata.dist.SDplus <- custom.newdata(dyad.ep2, "epFert", 
                                      "soc.fem_throat_avg_bright_scaled",
                                      "fem.throat", soc.pair=NA)
newdata.dist.SDplus$pat.clutch.size <- 4
newdata.dist.SDplus$dist.scale <- 1

# new data for SDmean
newdata.dist.mean <- newdata.dist.SDplus
newdata.dist.mean$dist.scale <- 0

# new data for SDminus
newdata.dist.SDminus <- newdata.dist.SDplus
newdata.dist.SDminus$dist.scale <- -1

# calculate model predictions
# FIRST RUN FUNCTION TO CALCULATE GLMM CIS

fit.prop.SDplus <- calculate.glmm.ci(prop.intx.dist.socFthroat,
                                     newdata.dist.SDplus)

fit.prop.mean <- calculate.glmm.ci(prop.intx.dist.socFthroat,
                                   newdata.dist.mean)

fit.prop.SDminus <- calculate.glmm.ci(prop.intx.dist.socFthroat,
                                      newdata.dist.SDminus)

# plot with CI and reversed brightness
ggplot(fit.prop.SDplus, aes(x=fem.throat*-1, 
                            y=fit)) +
  geom_line(color="#FFD393", linewidth=1) +
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="#FFD393", alpha=0.3) + 
  geom_line(data=fit.prop.mean, 
            aes(x=fem.throat*-1, y=fit), 
            color="#FFB64B",linewidth=1) +
  geom_ribbon(data=fit.prop.mean,
              aes(ymax=upper, ymin=lower), fill="#FFB64B", alpha=0.3) +
  geom_line(data=fit.prop.SDminus, 
            aes(x=fem.throat*-1, y=fit), 
            color="#FF9900",linewidth=1) +
  geom_ribbon(data=fit.prop.SDminus,
              aes(ymax=upper, ymin=lower), fill="#FF9900", alpha=0.3) +
  theme_light() +
  xlab("Social female throat brightness") +
  ylab("Proportion of the clutch with shared fertilizations")

ggsave("output-files/prop mod dist socFthroat intx plot CIs reverse.png", h=4, w=3.5)


# plot with CI
fthroat.intx <- ggplot(fit.prop.SDplus, aes(x=fem.throat, 
                            y=fit)) +
  geom_line(color="darkblue", linewidth=1) +
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="darkblue", alpha=0.3) + 
  geom_line(data=fit.prop.mean, 
            aes(x=fem.throat, y=fit), 
            color="blue",linewidth=1) +
  geom_ribbon(data=fit.prop.mean,
              aes(ymax=upper, ymin=lower), fill="blue", alpha=0.3) +
  geom_line(data=fit.prop.SDminus, 
            aes(x=fem.throat, y=fit), 
            color="lightblue",linewidth=1) +
  geom_ribbon(data=fit.prop.SDminus,
              aes(ymax=upper, ymin=lower), fill="lightblue", alpha=0.3) +
  theme_light() +
  xlab("Social female throat brightness") +
  ylab("Proportion of the clutch with shared fertilizations")

fthroat.intx

# ggsave("output-files/prop mod dist socFthroat intx plot CIs reverse.png", h=4, w=3.5)


# plot without CI bands
ggplot(fit.prop.SDplus, aes(x=fem.throat, 
                            y=fit)) +
  geom_line(color="darkblue", linewidth=2) +
  geom_line(data=fit.prop.mean, 
            aes(x=fem.throat, y=fit), 
            color="blue",linewidth=2) +
  geom_line(data=fit.prop.SDminus, 
            aes(x=fem.throat, y=fit), 
            color="lightblue",linewidth=2) +
  theme_light() +
  xlab("Social female throat brightness") +
  ylab("Proportion of the clutch with\nshared fertilizations")

ggsave("output-files/prop mod dist socFthroat intx plot.png", h=3, w=5)

# Without CI and zoomed in
fthroat.intx.noci <- ggplot(fit.prop.SDplus, aes(x=fem.throat, 
                            y=fit)) +
  geom_line(color="darkblue", linewidth=2) +
  geom_line(data=fit.prop.mean, 
            aes(x=fem.throat, y=fit), 
            color="blue",linewidth=2) +
  geom_line(data=fit.prop.SDminus, 
            aes(x=fem.throat, y=fit), 
            color="lightblue",linewidth=2) +
  theme_light() +
  xlab("Social female throat brightness") +
  ylab("Proportion of the clutch with\nshared fertilizations") +
  ylim(0,0.02)

fthroat.intx.noci

# soc Fem throat intx distance on X-axis ---------------------------------------

# new data for SDplus
newdata.SFthroat.SDplus <- custom.newdata(dyad.ep2, "epFert", 
                                      "dist.scale",
                                      "dist", soc.pair=NA)




newdata.SFthroat.SDplus$pat.clutch.size <- 4
newdata.SFthroat.SDplus$soc.fem_throat_avg_bright_scaled <- 1

# new data for SDmean
newdata.SFthroat.mean <- newdata.SFthroat.SDplus
newdata.SFthroat.mean$soc.fem_throat_avg_bright_scaled <- 0

# new data for SDminus
newdata.SFthroat.SDminus <- newdata.SFthroat.SDplus
newdata.SFthroat.SDminus$soc.fem_throat_avg_bright_scaled <- -1

# calculate model predictions
# FIRST RUN FUNCTION TO CALCULATE GLMM CIS

fit.prop2.SDplus <- calculate.glmm.ci(prop.intx.dist.socFthroat,
                                     newdata.SFthroat.SDplus)

fit.prop2.mean <- calculate.glmm.ci(prop.intx.dist.socFthroat,
                                   newdata.SFthroat.mean)

fit.prop2.SDminus <- calculate.glmm.ci(prop.intx.dist.socFthroat,
                                      newdata.SFthroat.SDminus)

# plot with CI
ggplot(fit.prop2.SDplus, aes(x=dist,  y=fit)) +
  geom_line(color="darkblue", linewidth=1) +
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="darkblue", alpha=0.3) + 
  geom_line(data=fit.prop2.mean, 
            aes(x=dist, y=fit), 
            color="blue",linewidth=1) +
  geom_ribbon(data=fit.prop2.mean,
              aes(ymax=upper, ymin=lower), fill="blue", alpha=0.3) +
  geom_line(data=fit.prop2.SDminus, 
            aes(x=dist, y=fit), 
            color="lightblue",linewidth=1) +
  geom_ribbon(data=fit.prop2.SDminus,
              aes(ymax=upper, ymin=lower), fill="lightblue", alpha=0.3) +
  theme_light() +
  xlab("Nest distance (m)") +
  ylab("Proportion of the clutch with shared fertilizations")


# without CI
fdist.intx.noci <- ggplot(fit.prop2.SDplus, aes(x=dist,  y=fit)) +
  geom_line(color="burlywood2", linewidth=2) +
  geom_line(data=fit.prop2.mean, 
            aes(x=dist, y=fit), 
            color="chocolate",linewidth=2) +
  geom_line(data=fit.prop2.SDminus, 
            aes(x=dist, y=fit), 
            color="coral4",linewidth=2) +
  theme_light() +
  xlab("Nest distance (m)") +
  ylab("Proportion of the clutch\nwith shared fertilizations")

fdist.intx.noci

mean(dyad.ep2$fem_throat_avg_bright) # 15.3
sd(dyad.ep2$fem_throat_avg_bright) # 2.57



# soc male breast---------------------------------------------------------------

newdata.prop.epFert.socMbreast <- custom.newdata(dyad.ep2, "epFert", 
                                                 "soc.male_breast_avg_bright_scaled",
                                           "male.breast", soc.pair=NA)
newdata.prop.epFert.socMbreast$pat.clutch.size <- 4

# predictions and CIs
fit.prop.epFert.socMbreast <- calculate.glmm.ci(prop.ep.global1, 
                                          newdata.prop.epFert.socMbreast)

# without raw data
socMbreast <- ggplot(fit.prop.epFert.socMbreast, aes(x=male.breast, y=fit)) +
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="black", alpha=0.3) +
  geom_line(color="black", linewidth=1) +
  ylab("Proportion of clutch with\nshared fertilizations") +
  xlab("Social male breast brightness") +
  theme_light()

socMbreast

# ggsave("output-files/prop mod socMbreast plot.png", h=3, w=5)




# male breast-------------------------------------------------------------------

newdata.prop.epFert.Mbreast <- custom.newdata(dyad.ep2, "epFert", 
                                                 "male_breast_avg_bright_scaled",
                                                 "male.breast", soc.pair=NA)
newdata.prop.epFert.Mbreast$pat.clutch.size <- 4

# predictions and CIs
fit.prop.epFert.Mbreast <- calculate.glmm.ci(prop.ep.global1, 
                                                newdata.prop.epFert.Mbreast)

# without raw data
Mbreast <- ggplot(fit.prop.epFert.Mbreast, aes(x=male.breast, y=fit)) +
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="black", alpha=0.3) +
  geom_line(color="black", linewidth=1) +
  ylab("Proportion of clutch with\nshared fertilizations") +
  xlab("Focal male breast brightness") +
  theme_light()

Mbreast

# ggsave("output-files/prop mod male breast plot.png", h=3, w=5)


# female breast-----------------------------------------------------------------

newdata.prop.epFert.Fbreast <- custom.newdata(dyad.ep2, "epFert", 
                                              "fem_breast_avg_bright_scaled",
                                              "fem.breast", soc.pair=NA)
newdata.prop.epFert.Fbreast$pat.clutch.size <- 4

# predictions and CIs
fit.prop.epFert.Fbreast <- calculate.glmm.ci(prop.ep.global1, 
                                             newdata.prop.epFert.Fbreast)

# without raw data
Fbreast <- ggplot(fit.prop.epFert.Fbreast, aes(x=fem.breast, y=fit)) +
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="black", alpha=0.3) +
  geom_line(color="black", linewidth=1) +
  ylab("Proportion of clutch with\nshared fertilizations") +
  xlab("Focal female breast brightness") +
  theme_light()

Fbreast

# ggsave("output-files/prop mod female breast plot.png", h=3, w=5)


# Combine plots ----------------------------------------------------------------

# first combine phenotype plots
# pheno.plot <- ggarrange((Fbreast + ylab(NULL) + ggtitle(" ")),  
#                         Mbreast + ylab(NULL)+ ggtitle(NULL), 
#                         fthroat.intx + ylab(NULL)+ ggtitle(NULL), 
#                         socMbreast + ylab(NULL)+ ggtitle(NULL), 
#           nrow=2, ncol=2, align = "hv",
#           labels = c("B","C","D","E"), hjust = -1)
# 
# pheno.plot

add.plot <- ggarrange(dist.pred + ylab(NULL), 
                       Fbreast + ylab(NULL) + ggtitle(" "),  
                      Mbreast + ylab(NULL)+ ggtitle(NULL), 
                      socMbreast + ylab(NULL)+ ggtitle(NULL), 
                      nrow=2, ncol=2, align = "hv",
                      labels = c("AUTO"), hjust = -1)
add.plot

# add intx
dist.pheno.plot <- ggarrange(add.plot,
                             fthroat.intx + ylab(NULL) + ggtitle(" "),
                             ncol=1, nrow=2,
                             labels=c(" ", "E"), hjust = -1,
                             heights=c(3,2))
dist.pheno.plot

full <- annotate_figure(
  dist.pheno.plot,
  left = text_grob("Proportion of clutch with shared fertilizations", 
                   rot = 90))

full

ggsave("output-files/main and intx effects plots.png", h=8, w=6)


# combined intx plots for sup mat ----------------------------------------------

sup.intx <- ggarrange(fthroat.intx.noci + ggtitle(" "), 
                      fdist.intx.noci + ylab(NULL) + ggtitle(" "),
                      nrow=1, ncol=2, labels = c("AUTO"))
sup.intx

ggsave("output-files/intx for sup mat.png", h=3, w=7)
