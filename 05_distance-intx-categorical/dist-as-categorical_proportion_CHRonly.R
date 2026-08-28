
################################################################################
# Explore distance as a categorical variable, CO Fert Net 2022, proportion outcome
# Heather Kenny-Duddela
# May 20, 2025
################################################################################

## load data

# load data
# from script adjust-table.R
dyad <- read.csv("output-files/dyad-complete-cases-CHRonly.csv")
# from script adjust-table.R
dyad.ep2 <- read.csv("output-files/dyad_ep_complete_cases_CHRonly.csv")


## libraries
library(ggplot2)
library(dplyr)
library(bbmle) # for AICtab
library(glmmTMB) # also for fitting mixed effect glm
library(DHARMa) # model diagnostics for mixed effect models]

# format variables
dyad$soc.pair <- factor(dyad$soc.pair, levels=c(1,0))

dyad$dist.cat <- factor(dyad$dist.cat, levels=c("dist.nest", "dist.01",
                                                "dist.05","dist.10","dist.20",
                                                "dist.20+"))

dyad$dist.cat2 <- factor(dyad$dist.cat2, levels=c("dist.nest",
                                                  "dist.05","dist.10","dist.20",
                                                  "dist.20+"))

dyad$dist.cat3 <- factor(dyad$dist.cat3, levels=c("dist.nest",
                                                  "dist.05","dist.10","dist.20",
                                                  "dist.20+"))

dyad.ep2$dist.cat <- factor(dyad.ep2$dist.cat, levels=c("dist.nest", "dist.01",
                                                "dist.05","dist.10","dist.20",
                                                "dist.20+"))

dyad.ep2$dist.cat2 <- factor(dyad.ep2$dist.cat2, levels=c("dist.nest",
                                                  "dist.05","dist.10","dist.20",
                                                  "dist.20+"))

dyad.ep2$dist.cat3 <- factor(dyad.ep2$dist.cat3, levels=c("dist.nest",
                                                  "dist.05","dist.10","dist.20",
                                                  "dist.20+"))


#-------------------------------------------------------------------------------
# proportion allFert, explore distance categories
#-------------------------------------------------------------------------------

# global model with continuous distance - dist is significant
prop.allFert.global <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                dist.scale + dur.scale + soc.pair +
                                fem_throat_avg_bright_scaled + 
                                fem_breast_avg_bright_scaled +
                                fem_TS_scaled +
                                male_throat_avg_bright_scaled +
                                male_breast_avg_bright_scaled + 
                                male_TS_scaled,
                               weights=pat.clutch.size, 
                              data=dyad, family=binomial)

summary(prop.allFert.global)



# don't include soc_pair as variable here because it is completely correlated
# with the dist.nest category

# fist dist category - DOES NOT CONVERGE
prop.allFert.global.dcat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                dist.cat + dur.scale +
                                fem_throat_avg_bright_scaled + 
                                fem_breast_avg_bright_scaled +
                                fem_TS_scaled +
                                male_throat_avg_bright_scaled +
                                male_breast_avg_bright_scaled + 
                                male_TS_scaled,
                              weights=pat.clutch.size, data=dyad, 
                              family=binomial)

summary(prop.allFert.global.dcat)



# fit model again with cat2 - all categories significant
prop.allFert.global.dcat2 <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                                     dist.cat2 + dur.scale +
                                     fem_throat_avg_bright_scaled + 
                                     fem_breast_avg_bright_scaled +
                                     fem_TS_scaled +
                                     male_throat_avg_bright_scaled +
                                     male_breast_avg_bright_scaled + 
                                     male_TS_scaled,
                                   weights=pat.clutch.size, data=dyad, 
                                   family=binomial)

summary(prop.allFert.global.dcat2)
# all categories sig and similar (-5.1 TO -6.5)


# fit model again with cat3 - similar results
prop.allFert.global.dcat3 <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                                      dist.cat3 + dur.scale +
                                      fem_throat_avg_bright_scaled + 
                                      fem_breast_avg_bright_scaled +
                                      fem_TS_scaled +
                                      male_throat_avg_bright_scaled +
                                      male_breast_avg_bright_scaled + 
                                      male_TS_scaled,
                                    weights=pat.clutch.size, data=dyad, 
                                    family=binomial)

summary(prop.allFert.global.dcat3)
# similar to dcat2


# reduced model with just significant trait variables
prop.allFert.maleTS.dcat3 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                      dist.cat3 + dur.scale +
                                       fem_throat_avg_bright_scaled +
                                      male_TS_scaled,
                                    weights=pat.clutch.size, data=dyad, 
                                    family=binomial)

summary(prop.allFert.maleTS.dcat3)

# reduced model with continuous distance - dist is significant
prop.allFert.maleTS.dist <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                      soc.pair +
                                      dist.scale + dur.scale +
                                      fem_throat_avg_bright_scaled +
                                      male_TS_scaled,
                                    weights=pat.clutch.size, data=dyad, 
                                    family=binomial)

summary(prop.allFert.maleTS.dist)


#-------------------------------------------------------------------------------
# proportion allFert, interaction between dist and trait variables and dur
#-------------------------------------------------------------------------------

# fem_throat - marginal interaction with dist20 but not the other categories
prop.allFert.dcat3.fthroat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                        dur.scale +
                                      dist.cat3*fem_throat_avg_bright_scaled,
                                    weights=pat.clutch.size, data=dyad, 
                                    family=binomial)

summary(prop.allFert.dcat3.fthroat)

# fem_breast - significant for dist10 and dist20+, but not dist20...
prop.allFert.dcat3.fbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                        dur.scale +
                                       dist.cat3*fem_breast_avg_bright_scaled,
                                     weights=pat.clutch.size, data=dyad, 
                                     family=binomial)

summary(prop.allFert.dcat3.fbreast)

# fem_TS - marginal interaction with dist10 but not others
prop.allFert.dcat3.fTS <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                    dur.scale +
                                       dist.cat3*fem_TS_scaled,
                                     weights=pat.clutch.size, data=dyad, 
                                  family=binomial)

summary(prop.allFert.dcat3.fTS)

# male_throat - no significant interactions
prop.allFert.dcat3.mthroat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                        dur.scale +
                                   dist.cat3*male_throat_avg_bright_scaled,
                                 weights=pat.clutch.size, data=dyad, 
                                 family=binomial)
summary(prop.allFert.dcat3.mthroat)

# male_breast - significant with dist10, marginal with dist20
prop.allFert.dcat3.mbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                        dur.scale +
                                       dist.cat3*male_breast_avg_bright_scaled,
                                     weights=pat.clutch.size, data=dyad, 
                                     family=binomial)
summary(prop.allFert.dcat3.mbreast)


# dur - not significant
prop.allFert.dcat3.dur <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                    dur.scale*dist.cat3,
                                     weights=pat.clutch.size, data=dyad, 
                                  family=binomial)
summary(prop.allFert.dcat3.dur)

#-------------------------------------------------------------------------------
# proportion EP, explore distance categories
#-------------------------------------------------------------------------------


# global model with continuous distance - dist IS significant
prop.global.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                           dist.scale + dur.scale +
                           male_breast_avg_bright_scaled +
                           soc.fem_throat_avg_bright_scaled +
                           soc.male_TS_scaled,
                         weights=pat.clutch.size, data=dyad.ep2, 
                         family=binomial)

summary(prop.global.ep)

# with dist.cat - dist is NOT significant
prop.global.ep.dcat <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                           dist.cat + dur.scale +
                           male_breast_avg_bright_scaled +
                           soc.fem_throat_avg_bright_scaled +
                           soc.male_TS_scaled,
                         weights=pat.clutch.size, data=dyad.ep2, 
                         family=binomial)

summary(prop.global.ep.dcat)
# none of the categories are significant, and the effect sizes are all very 
# similar (14 to 15)

# with dist.cat3 - dist20 sig but not dist20+
prop.global.ep.dcat3 <- glmmTMB(prop.fert ~ (1|female) + (1|male) +  
                                dist.cat3 + dur.scale +
                                male_breast_avg_bright_scaled +
                                soc.fem_throat_avg_bright_scaled +
                                soc.male_TS_scaled,
                              weights=pat.clutch.size, data=dyad.ep2, 
                              family=binomial)

summary(prop.global.ep.dcat3)


## Larger global model - dist is significant
prop.global.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                           dist.scale + dur.scale +
                           soc.fem_throat_avg_bright_scaled +
                           soc.male_throat_avg_bright_scaled +
                           soc.fem_breast_avg_bright_scaled +
                           soc.male_TS_scaled +
                           soc.fem_TS_scaled +
                           male_TS_scaled +
                           male_breast_avg_bright_scaled +
                           male_throat_avg_bright_scaled +
                           fem_breast_avg_bright_scaled,
                         weights=pat.clutch.size, data=dyad.ep2, 
                         family=binomial)
summary(prop.global.ep)

# with cat3 - dist20 is significant
prop.global.ep.dcat3 <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                           dist.cat3 + dur.scale +
                           soc.fem_throat_avg_bright_scaled +
                           soc.male_throat_avg_bright_scaled +
                           soc.fem_breast_avg_bright_scaled +
                           soc.male_TS_scaled +
                           soc.fem_TS_scaled +
                           male_TS_scaled +
                           male_breast_avg_bright_scaled +
                           male_throat_avg_bright_scaled +
                           fem_breast_avg_bright_scaled,
                         weights=pat.clutch.size, data=dyad.ep2, 
                         family=binomial)
summary(prop.global.ep.dcat3)

## Univariate models

# continuous distance - significant
prop.uni.ep <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                            dist.scale,
                          weights=pat.clutch.size, data=dyad.ep2, 
                          family=binomial)
summary(prop.uni.ep)

# dist cat - none significant
prop.uni.ep.dcat <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                         dist.cat,
                       weights=pat.clutch.size, data=dyad.ep2, 
                       family=binomial)
summary(prop.uni.ep.dcat)

# dist.cat2 - dist20 significant, all are negative but strongest at 20
prop.uni.ep.dcat2 <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                              dist.cat2,
                            weights=pat.clutch.size, data=dyad.ep2, 
                            family=binomial)
summary(prop.uni.ep.dcat2)

# dist.cat3 - dist20 is significant
prop.uni.ep.dcat3 <- glmmTMB(prop.fert ~ (1|female) + (1|male) + 
                               dist.cat3,
                             weights=pat.clutch.size, data=dyad.ep2, 
                             family=binomial)
summary(prop.uni.ep.dcat3)

## Compare models - continuous distance is top, dcat3 is next best
AICtab(prop.uni.ep, prop.uni.ep.dcat, prop.uni.ep.dcat2, prop.uni.ep.dcat3)


#-------------------------------------------------------------------------------
# proportion EP, interaction between dist and trait variables and dur
#-------------------------------------------------------------------------------

# fem_throat - none significant
prop.EPfert.dcat3.fthroat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                       dist.cat3*fem_throat_avg_bright_scaled,
                                     weights=pat.clutch.size, data=dyad.ep2, 
                                     family=binomial)

summary(prop.EPfert.dcat3.fthroat)

# fem_breast - no significant interactions
prop.EPfert.dcat3.fbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                       dist.cat3*fem_breast_avg_bright_scaled,
                                     weights=pat.clutch.size, data=dyad.ep2, 
                                     family=binomial)

summary(prop.EPfert.dcat3.fbreast)

# fem_TS - no significant interactions
prop.EPfert.dcat3.fTS <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                   dist.cat3*fem_TS_scaled,
                                 weights=pat.clutch.size, data=dyad.ep2, 
                                 family=binomial)

summary(prop.EPfert.dcat3.fTS)

# male_throat - no significant interactions
prop.EPfert.dcat3.mthroat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                       dist.cat3*male_throat_avg_bright_scaled,
                                     weights=pat.clutch.size, data=dyad.ep2, 
                                     family=binomial)
summary(prop.EPfert.dcat3.mthroat)

# male_breast - no significant interactions
prop.EPfert.dcat3.mbreast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                       dist.cat3*male_breast_avg_bright_scaled,
                                     weights=pat.clutch.size, data=dyad.ep2, 
                                     family=binomial)
summary(prop.EPfert.dcat3.mbreast)


# dur - not significant
prop.EPfert.dcat3.dur <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                   dur.scale*dist.cat3,
                                 weights=pat.clutch.size, data=dyad.ep2, 
                                 family=binomial)
summary(prop.EPfert.dcat3.dur)

# soc.fem_throat - intx with dist20 is significant!
prop.EPfert.dcat3.socF_throat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                    dist.cat3*soc.fem_throat_avg_bright_scaled,
                                weights=pat.clutch.size, data=dyad.ep2, 
                                family=binomial)
summary(prop.EPfert.dcat3.socF_throat)

# soc.fem_breast - intx with 20+ is marginal
prop.EPfert.dcat3.socF_breast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                          dist.cat3*soc.fem_breast_avg_bright_scaled,
                                        weights=pat.clutch.size, data=dyad.ep2, 
                                        family=binomial)
summary(prop.EPfert.dcat3.socF_breast)

# soc.fem_TS - no significant intx
prop.EPfert.dcat3.socF_TS <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                          dist.cat3*soc.fem_TS_scaled,
                                        weights=pat.clutch.size, data=dyad.ep2, 
                                     family=binomial)
summary(prop.EPfert.dcat3.socF_TS)

# soc.male_throat - marginal with dist20
prop.EPfert.dcat3.socM_throat <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                      dist.cat3*soc.male_throat_avg_bright_scaled,
                                    weights=pat.clutch.size, data=dyad.ep2, 
                                    family=binomial)
summary(prop.EPfert.dcat3.socM_throat)

# soc.male_breast - no significant interactions
prop.EPfert.dcat3.socM_breast <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                          dist.cat3*soc.male_breast_avg_bright_scaled,
                                        weights=pat.clutch.size, data=dyad.ep2, 
                                        family=binomial)
summary(prop.EPfert.dcat3.socM_breast)

# soc.male_TS - significant with dist20
prop.EPfert.dcat3.socM_TS <- glmmTMB(prop.fert ~ (1|female) + (1|male) +
                                          dist.cat3*soc.male_TS_scaled,
                                        weights=pat.clutch.size, data=dyad.ep2, 
                                     family=binomial)
summary(prop.EPfert.dcat3.socM_TS)





