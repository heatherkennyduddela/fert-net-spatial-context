
################################################################################
# Additional adjustments to tables for final models
# Heather Kenny-Duddela
# May 21, 2025
# updated May 28, 2026
################################################################################

# libraries
library(tidyverse)

# load data
# from script add-phenotype-CHRonly.R
dyad <- read.csv("input-files/fert dist plumage-CHRonly_with fertile and rest SI.csv")
# from script summary-stats-and-corr.R
# pair <- read.csv("input-files/pair pheno.csv")

# from script summary-stats-and-corr-CHRonly.R
pair.chr <- read.csv("input-files/pair pheno-CHRonly.csv")

#-------------------------------------------------------------------------------
# Check difference between pair and pair.chr
# 
# check.pair <- left_join(select(pair.chr,female, male,
#                                fem_TS_scaled, 
#                                fem_throat_avg_bright_scaled,
#                                fem_breast_avg_bright_scaled, 
#                                male_TS_scaled,
#                                male_throat_avg_bright_scaled,
#                                male_breast_avg_bright_scaled),
#                         select(pair, female, male,
#                                          fem_TS_scaled, 
#                                          fem_throat_avg_bright_scaled,
#                                          fem_breast_avg_bright_scaled, 
#                                          male_TS_scaled,
#                                          male_throat_avg_bright_scaled,
#                                          male_breast_avg_bright_scaled),
#                         by=c("female","male")) 
# 
# check.pair.fem <- select(check.pair, fem_TS_scaled.x, fem_TS_scaled.y,
#                          fem_throat_avg_bright_scaled.x, fem_throat_avg_bright_scaled.y,
#                          fem_breast_avg_bright_scaled.x, fem_breast_avg_bright_scaled.y)

#-------------------------------------------------------------------------------

## Add scaled distance and duration
dyad$dist.scale <- scale(dyad$process.dist.m)
dyad$dur.scale <- scale(dyad$dur.fertile)
colnames(dyad)[72:73]

# add bin_num_fert to have it in zeros and ones
dyad$bin_num_fert <- ifelse(dyad$bin_fert=="yes", 1, 0)

# add site size category, only 2 because not enough data to properly estimate
# effects for small or medium sites
dyad$size.cat <- "other"
dyad$size.cat[which(dyad$barn=="CHR")] <- "large"


# add 2-level attempt column instead of 3 level
dyad$attempt2 <- ifelse(dyad$attempt1==1, "first", "second")

# make soc.pair==1 the intercept condition
dyad$soc.pair <- factor(dyad$soc.pair, levels=c(1,0))


## Add distance as categorical

## create distance bins, 0.08, 0.08-1, 1-5, 5-10, 10-20 (CHR and MB),
# 20+ (CHR only)

ggplot(dyad, aes(x=process.dist.m)) + geom_histogram()

# cat1
dyad$dist.cat <- NA
dyad$dist.cat[which(dyad$process.dist.m==0.08)] <- "dist.nest"
dyad$dist.cat[which(dyad$process.dist.m>0.08 & 
                      dyad$process.dist.m<=1)] <- "dist.01"
dyad$dist.cat[which(dyad$process.dist.m>1 & 
                      dyad$process.dist.m<=5)] <- "dist.05"
dyad$dist.cat[which(dyad$process.dist.m>5 & 
                      dyad$process.dist.m<=10)] <- "dist.10"
dyad$dist.cat[which(dyad$process.dist.m>10 & 
                      dyad$process.dist.m<=20)] <- "dist.20"
dyad$dist.cat[which(dyad$process.dist.m>20)] <- "dist.20+"

dyad$dist.cat <- factor(dyad$dist.cat, levels=c("dist.nest", "dist.01",
                                                "dist.05","dist.10","dist.20",
                                                "dist.20+"))

# cat 2
## Try second dist.cat, 0.08, 0.08-5, 5-10, 10-20, 20+

dyad$dist.cat2 <- NA
dyad$dist.cat2[which(dyad$process.dist.m==0.08)] <- "dist.nest"
dyad$dist.cat2[which(dyad$process.dist.m>0.08 & 
                       dyad$process.dist.m<=5)] <- "dist.05"
dyad$dist.cat2[which(dyad$process.dist.m>5 & 
                       dyad$process.dist.m<=10)] <- "dist.10"
dyad$dist.cat2[which(dyad$process.dist.m>10 & 
                       dyad$process.dist.m<=20)] <- "dist.20"
dyad$dist.cat2[which(dyad$process.dist.m>20)] <- "dist.20+"

dyad$dist.cat2 <- factor(dyad$dist.cat2, levels=c("dist.nest",
                                                  "dist.05","dist.10","dist.20",
                                                  "dist.20+"))

# cat3
## Try third dist cat with nest, 0.08-10, 10-20, 20+

dyad$dist.cat3 <- NA
dyad$dist.cat3[which(dyad$process.dist.m==0.08)] <- "dist.nest"
dyad$dist.cat3[which(dyad$process.dist.m>0.08 & 
                       dyad$process.dist.m<=10)] <- "dist.10"
dyad$dist.cat3[which(dyad$process.dist.m>10 & 
                       dyad$process.dist.m<=20)] <- "dist.20"
dyad$dist.cat3[which(dyad$process.dist.m>20)] <- "dist.20+"

dyad$dist.cat3 <- factor(dyad$dist.cat3, levels=c("dist.nest",
                                                  "dist.05","dist.10","dist.20",
                                                  "dist.20+"))



# add column for proportion shared fert
dyad$prop.fert <- dyad$shared_fert/dyad$pat.clutch.size



# keep only relevant trait variables (tail, throat bright, breast bright)
dyad2 <- dyad %>% select(!contains("belly", ignore.case=T)) %>%
  select(!contains("vent")) %>%
  select(!contains("chroma")) %>%
  select(!contains("hue")) %>%
  select(!contains("total_")) %>%
  select(!contains("mean"))


# save full dyad with CHR only
write.csv(dyad2, "output-files/dyad-table-CHRonly_with fertile and rest SI.csv", 
          row.names=F)


# complete cases for CHRonly
dyad3 <- dyad2[complete.cases(dyad2[,20:33]), ]
#re-calculate scaled dist and dur
dyad3$dist.scale <- scale(dyad3$process.dist.m)
dyad3$dur.scale <- scale(dyad3$dur.fertile)
dyad3$rest.SI.scale <- scale(dyad3$rest.SI.chr)
write.csv(dyad3, 
          "output-files/dyad-complete-cases-CHRonly_with fertile and rest SI.csv", 
          row.names=F)

dyad4 <- dyad3



#-------------------------------------------------------------------------------
# Adjustments for ep only table
#-------------------------------------------------------------------------------

# make ep only table
dyad.ep <- subset(dyad4, dyad4$soc.pair==0)

# add social pair phenotype info
fem.mate <- select(pair.chr, female, male_TS_scaled:male_vent_avg_bright_scaled)
colnames(fem.mate)[2:6] <- paste("soc.", colnames(fem.mate)[2:6], sep="")

dyad.ep.fmate <- left_join(dyad.ep, fem.mate, by="female")

# Male 2850-57582 was simultaneously socially paired with 2 females, 2850-57746 
# at CHR-6, and 2870-64909 at CHR-148. Nest 6 became active about 1 week before
# nest 148, so keep only female 57746 as the social female
male.mate <- select(pair.chr, male, female, fem_TS_scaled:fem_vent_avg_bright_scaled)
male.mate2 <- male.mate[-which(male.mate$male=="2850-57582" & 
                                 male.mate$female=="2870-64909"), ]
# Same issue for male 2850-57806, socially paired with 2850-57568 and 2850-57735.
# Zero WP shared with 57568, and 2 WP shared for each attempt with 57735. 
# Keep 57735 as social female
male.mate3 <- male.mate2[-which(male.mate$male=="2850-57806" & 
                                  male.mate$female=="2850-57568"), ]

colnames(male.mate3)[3:7] <- paste("soc.", colnames(male.mate3)[3:7], sep="")

dyad.ep2 <- left_join(dyad.ep.fmate, male.mate3[ ,-2], by="male")


# keep only relevant trait variables (tail, throat bright, breast bright)
dyad.ep3 <- dyad.ep2 %>% select(!contains("belly", ignore.case=T)) %>%
  select(!contains("vent"))



## ep table for CHR only
write.csv(dyad.ep3, "output-files/dyad_ep_CHRonly_with fertile and rest SI.csv", 
          row.names=F)


# complete cases for CHRonly
dyad.ep4 <- dyad.ep3[complete.cases(dyad.ep3[,c(20:33,44:49)]), ]
# re-calculate scaled dist and dur and rest.SI
dyad.ep4$dist.scale <- scale(dyad.ep4$process.dist.m)
dyad.ep4$dur.scale <- scale(dyad.ep4$dur.fertile)
dyad.ep4$rest.SI.scale <- scale(dyad.ep4$rest.SI.chr) # correct set for scaling, previous variable included the wrong set of observations
write.csv(dyad.ep4, 
          "output-files/dyad_ep_complete_cases_CHRonly_with fertile and rest SI.csv",
          row.names=F)
