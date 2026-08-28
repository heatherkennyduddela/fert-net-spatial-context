#-------------------------------------------------------------------------------
# Calculate number of total unknown sires for 2021  
# Heather Kenny-Duddela
# May 23, 2026
#-------------------------------------------------------------------------------


# libraries
library(tidyverse)
library(ggplot2)

### Set working directory
setwd("~/CU Boulder/BARS fieldwork/2022 Field and Lab/Ch 3 analyses/02_paternity-assignment/2021")


### load data tables for 2021

# full kinship table with metadata
kin2.21 <- read.csv("01_output-files/kin_2021_all_with_labels.csv")

# table of offspring with assigned parents
assigned <- read.csv("02_output-files/kin_2021_parent_offspring_assigned_wide.csv")


#-------------------------------------------------------------------------------

################################################################################
#-------------------------------------------------------------------------------
# Pseudo code:
# 1) Get list of true half-sibs based on parentage assignment
#   a) loop through table of offspring with assigned parents and identify
#   kids with same mom but different dad, or same dad but different mom
#   b) For each set of kids with one parent in common, save their IDs and then
#   extract kinship relationships for that set of kids
#-------------------------------------------------------------------------------

# remove kids with unknown mom or dad
assigned.noNA <- subset(assigned, !is.na(assigned$Band_dad) &
                          !is.na(assigned$Band_mom))

storage.hs2 <- as.data.frame(matrix(nrow=1,ncol=6,NA))
colnames(storage.hs2) <- c("Ind1","Ind2","mom_same","dad_same", "mom_band", "dad_band")

# loop through kids
for (i in 1:length(assigned.noNA$Ind2)) {
  # loop through moms and dads
  for (j in 1:length(assigned.noNA$Band_mom)) {
    # case of same dad, different mom
    if(assigned.noNA$Band_dad[i] == assigned.noNA$Band_dad[j] &
       assigned.noNA$Band_mom[i] != assigned.noNA$Band_mom[j]) {
      storage <- c(assigned.noNA$Ind2[i], 
                   assigned.noNA$Ind2[j],
                   F,T, NA, assigned.noNA$Band_dad[i])
      storage.hs2 <- rbind(storage.hs2, storage)
    }
    # case of same mom, different dad
    if(assigned.noNA$Band_dad[i] != assigned.noNA$Band_dad[j] &
       assigned.noNA$Band_mom[i] == assigned.noNA$Band_mom[j]) {
      storage <- c(assigned.noNA$Ind2[i], 
                   assigned.noNA$Ind2[j],
                   T,F, assigned.noNA$Band_mom[i], NA)
      storage.hs2 <- rbind(storage.hs2, storage)
    }
    # case of same mom, same dad (full sibs)
    if(assigned.noNA$Band_dad[i] == assigned.noNA$Band_dad[j] &
       assigned.noNA$Band_mom[i] == assigned.noNA$Band_mom[j]) {
      storage <- c(assigned.noNA$Ind2[i], 
                   assigned.noNA$Ind2[j],
                   T,T, assigned.noNA$Band_mom[i], assigned.noNA$Band_dad[i])
      storage.hs2 <- rbind(storage.hs2, storage)
    }
    # case of diff mom, diff dad (unrelated)
    if(assigned.noNA$Band_dad[i] != assigned.noNA$Band_dad[j] &
       assigned.noNA$Band_mom[i] != assigned.noNA$Band_mom[j]) {
      storage <- c(assigned.noNA$Ind2[i], 
                   assigned.noNA$Ind2[j],
                   F,F, NA,NA)
      storage.hs2 <- rbind(storage.hs2, storage)
    }
  }
}

storage.hs2 <- storage.hs2[-1,]

# list of half sibs includes duplicate pairs (ex: AB and BA)
# do a full join to pull out half sib pairs from the kinship list

kin.hs2 <- inner_join(storage.hs2, kin2.21, by=c("Ind1","Ind2"))

# add column for sib relationship
kin.hs2$sib <- "unrelated"
kin.hs2$sib[which(kin.hs2$mom_same==T & kin.hs2$dad_same==T)] <- "full"
kin.hs2$sib[which(kin.hs2$mom_same==T & kin.hs2$dad_same==F |
                    kin.hs2$mom_same==F & kin.hs2$dad_same==T )] <- "half"

#-------------------------------------------------------------------------------
# Get kinship values for kids with unknown dads
#-------------------------------------------------------------------------------

# get list of kids with unknown dads
dad.unk <- subset(assigned, is.na(assigned$pi_HAT_dad))

# pull out kinship values for unrelated to unrelated
unk2unk <- kin2.21[which(kin2.21$Ind1 %in% dad.unk$Ind2 &
                           kin2.21$Ind2 %in% dad.unk$Ind2), ]

# plot unknown sibs over known relationships
ggplot(kin.hs2, aes(x=k0_hat, y=pi_HAT, color=sib)) +
  geom_point(shape=2) +
  geom_point(data=unk2unk, aes(x=k0_hat, y=pi_HAT), color="black")


# use pi_hat values for known half-sibs and unrelated to classify unknowns

# wide cutoff uses max and min of known half-sibs
max(subset(kin.hs2$pi_HAT, kin.hs2$sib=="half")) # 0.343
min(subset(kin.hs2$pi_HAT, kin.hs2$sib=="half")) # 0.088

# narrow cutoff uses min of full and max of unrelated
min(subset(kin.hs2$pi_HAT, kin.hs2$sib=="full")) # 0.227
max(subset(kin.hs2$pi_HAT, kin.hs2$sib=="unrelated")) # 0.143

# make overlay plot again with cutoff lines
sib_categories_2021 <- ggplot(kin.hs2, aes(x=k0_hat, y=pi_HAT, color=sib)) +
  geom_point(shape=2) +
  geom_point(data=unk2unk, aes(x=k0_hat, y=pi_HAT), color="black")+
  geom_hline(yintercept=0.343, linetype="dashed") +
  geom_hline(yintercept=0.088, linetype="dashed") +
  geom_hline(yintercept=0.227, linetype="solid") +
  geom_hline(yintercept=0.143, linetype="solid")

sib_categories_2021


# save plot
ggsave("03_output-files/unknown kinship overlayed with known sibs cutoff lines 2021.png", h=5, w=7)

save(sib_categories_2021, file="03_output-files/sib_categories_2021_plot.Rdata")


# add classifications to unknowns for wide and narrow cutoffs

# wide
unk2unk$sib_wide <- NA
unk2unk$sib_wide[unk2unk$pi_HAT>0.343] <- "full"
unk2unk$sib_wide[unk2unk$pi_HAT<=0.343 &
                   unk2unk$pi_HAT>=0.088] <- "half"
unk2unk$sib_wide[unk2unk$pi_HAT<=0.088] <- "unrelated"

# narrow
unk2unk$sib_narrow <- NA
unk2unk$sib_narrow[unk2unk$pi_HAT>0.227] <- "full"
unk2unk$sib_narrow[unk2unk$pi_HAT<=0.227 &
                     unk2unk$pi_HAT>=0.143] <- "half"
unk2unk$sib_narrow[unk2unk$pi_HAT<=0.143] <- "unrelated"


#-------------------------------------------------------------------------------
# Calculate total number of unknown sires
#-------------------------------------------------------------------------------

# first categorize matching or not matching familyIDs
unk2unk$famID_check <- ifelse(unk2unk$FamilyID_ind1 == unk2unk$FamilyID_ind2, T, F)

# check that all full sibs have same familyID
sum(unk2unk$sib_wide=="full" & unk2unk$famID_check==T) # 47
sum(unk2unk$sib_wide=="full") # 47

# check that unrelated have different familyIDs
sum(unk2unk$sib_wide=="unrelated" & unk2unk$famID_check==F) # 588
sum(unk2unk$sib_wide=="unrelated") # 588

# subset only full and half sibs
unk_fh <- subset(unk2unk, unk2unk$sib_wide != "unrelated")

# Same dad, different moms -----------------------------------------------------
# pull out cases of half sibs where family ID is different
# this means the kids have the same dad but different moms
unk_h_famID_diff <- subset(unk_fh, unk_fh$famID_check==F & unk_fh$sib_wide=="half")

# One kid from Cooks_22_1 (family Cooks-23) with one kid from Cooks_31_collected
# (family Cooks-31)
Cooks23_Cooks31_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Cooks-23" & 
                              unk2unk$FamilyID_ind2=="Cooks-31" |
                              unk2unk$FamilyID_ind1=="Cooks-31" &
                              unk2unk$FamilyID_ind2=="Cooks-23" |
                              unk2unk$FamilyID_ind1=="Cooks-23" &
                              unk2unk$FamilyID_ind2=="Cooks-23" |
                              unk2unk$FamilyID_ind1=="Cooks-31" &
                              unk2unk$FamilyID_ind2=="Cooks-31")
# Three of four kids from 22_1 (Cooks-23) are full sibs with each other and
# unrelated to 31_collected or Hepp. One kid from 22_1 is half sibs with the
# others in this clutch, and also half/full sibs with the kid from 31_collected. 
# This kid from 22_1 (82971) must be half sibs and not full with 31_collected
# because they have different moms. 2640-97370 is mom to 22_1 and 2850-57558
# is mom to 31_collected and Hepp. Male uns-1 sired 3/4 from 22_1. Male uns-2
# sired 1/4 from 22_1 AND the one from 31_collected. Kid from Hepp is half
# sibs with kid from 31_collected, and they have the same mom. This means male
# uns-3 sired the kid from Hepp. 


# Schaaps_97_2 (family Schaaps-108) with Schaaps-54_collected (family Schaaps-54)
Schaaps108_Schaaps54_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Schaaps-108" & 
                                  unk2unk$FamilyID_ind2=="Schaaps-54" |
                                  unk2unk$FamilyID_ind1=="Schaaps-54" &
                                  unk2unk$FamilyID_ind2=="Schaaps-108" |
                                  unk2unk$FamilyID_ind1=="Schaaps-108" &
                                  unk2unk$FamilyID_ind2=="Schaaps-108" |
                                  unk2unk$FamilyID_ind1=="Schaaps-54" &
                                  unk2unk$FamilyID_ind2=="Schaaps-54")
# All eggs from 54_collected (family Schaaps-54) are full sibs with mom 2850-57603 
# and have dad male uns-1. All are half sibs with the kid from 33_2 
# (family Schaaps-54) who has the same mom (57603) 
# but a different dad (male uns-2). E59 from 54_collected is half/full from the 
# wide/narrow categories, but they must be half with 33_2 because they are 
# full sibs with the rest of the eggs from 54_collected and 54_collected and
# 33_2 have the same dam. The kid from 97_2 (family Schaaps-108) is half sibs 
# with all the eggs from 54_collected and has a different mom (2850-57553)
# so the same dad (male uns-1). 


# same mom, different dads -----------------------------------------------------
# pull out cases of half sibs where family ID is same
# this means the kids have the same mom but different dads
unk_h_famID_same <- subset(unk_fh, unk_fh$famID_check==T & unk_fh$sib_wide=="half")

# BlueCloud-09
BlueCloud9_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="BlueCloud-09" &
                         unk2unk$FamilyID_ind2=="BlueCloud-09")
# Some uncertainty with E29 from 9_collected being called half (wide) with 83405
# from 21_1 but not the rest of that clutch, and E29 being called half with E30 
# but not the rest of that clutch. Logically it makes the most sense for there
# to be only one sire here. 

# Cooks-27
Cooks27_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Cooks-27" &
                             unk2unk$FamilyID_ind2=="Cooks-27")
# Kids within 27_collected are all full sibs, and kids within 29_2 are all
# full sibs. Uncertainty about relationships between clutches. Wide category
# is consistent saying all between-clutch relationships are half sib. Narrow
# category is inconsistent saying some between-clutch relationships are full
# and some are half. Wide makes more logical sense, that each clutch has a different
# unsampled sire. Males: uns-1 and uns-2

# Cooks-31, already checked above. Two different dads for kids from 31_collected
# and Hepp. 

# Cooks-23 also already checked above. One kid from 22_1 (82971) has a different
# dad than the other three in that clutch. Male uns-2 sired one kid in 22_1 and
# one kid in 31_collected. 

# Schaaps-54 also already checked above. Kids from 54_collected and 33_2 have
# the same mom but different dads. 

# Schaaps-71
Schaaps71_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Schaaps-71" &
                          unk2unk$FamilyID_ind2=="Schaaps-71")
# All kids are from 134_2. Most supported scenario is 83001, 83002, and 83005
# are full sibs, and they are half sibs to 83003. This implies two unsampled sires.
# Uncertainty about relationship between 001 and 002 (half or full), and between
# 003 and 005 (half or full). However, 005 is definitely full with 001 and 002, 
# and 001 and 002 are definitely half with 003. 001 and 002 can't both be full with 
# 005 but be half with each other. And 003 can't be full with 005.

# Schaaps-80
Schaaps80_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Schaaps-80" &
                            unk2unk$FamilyID_ind2=="Schaaps-80")
# The two kids from 80_2 are half sibs, so there were two unsampled sires for 
# this clutch. 

# Check UrbanFarm for good measure, should be all full sibs
UrbanFarm2_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="UrbanFarm-02" &
                             unk2unk$FamilyID_ind2=="UrbanFarm-02")
# All full sibs within and between clutches = just one sire

# determine number of unique family IDs-----------------------------------------
unique_fam <- unique(c(unk2unk$FamilyID_ind1, unk2unk$FamilyID_ind2)) # 10
length(unique_fam)
sort(unique_fam)

# Assume one uns sire for each dam, except cases discussed above
# Cooks-23 had two sires in a single clutch (+1)
# One of the sires from Cooks-23 was shared with Cooks-31 (-1)
# Cooks-31 had two different sires for her two clutches (+1)
# Schaaps-54 had two sires for her two clutches (+1)
# One of the sires from Schaaps-54 was shared with Schaaps-108 (-1)
# Blue Cloud-9 had only one sire (+0)
# Cooks-27 had two sires in two separate clutches (+1)
# Schaaps-71 had two sires in a single clutch (+1)
# Schaaps-80 had two sires in a single clutch (+1)

# Total: 10 +1 -1 +1 +1 -1 +0 +1 +1 +1 = 14


#-------------------------------------------------------------------------------
# Calculate total RS for each unsampled sire
#-------------------------------------------------------------------------------

# summarize unknown kids by clutch and famID
unk.sum <- dad.unk %>% group_by(clutch_id_ind2, FamilyID_mom) %>%
  summarise(n = n())

# make column for unsampled male ID
unk.sum$uns_male <- paste(unk.sum$FamilyID_mom,"uns1", sep="_")

# add extra rows for Cooks-23, Schaaps-71, and Schaaps-80 which each had
# multiple sires witin a single clutch
unk.sum2 <- rbind(unk.sum, unk.sum[c(5,9,12), ])

# update Schaaps-71 uns1 has 3 chicks, uns2 has 1
unk.sum2$n[9] <- 3
unk.sum2$n[17] <- 1
unk.sum2$uns_male[17] <- "Schaaps-71_uns2"

# update Schaaps-80 uns1 with 1 chick, and uns2 with 1 chick
unk.sum2$n[12] <- 1
unk.sum2$n[18] <- 1
unk.sum2$uns_male[18] <- "Schaaps-80_uns2"

# Sire from Cooks_31_collected is shared with 22_1 (one kid in 22_1)
# Update 22_1 uns-1 with 3 chicks, and update 22_1 with 1 chick and
# sire Cooks-31_uns1. Update Hepp with uns-2 
unk.sum2$n[5] <- 3
unk.sum2$n[16] <- 1
unk.sum2$uns_male[16] <- "Cooks-31_uns1"
unk.sum2$uns_male[7] <- "Cooks-31_uns2"

# update Schaaps-54 Schaaps_33_2 to uns-2
unk.sum2$uns_male[which(unk.sum2$clutch_id_ind2=="Schaaps_33_2")] <- "Schaaps-54_uns2"

# update Schaaps-108 (clutch 97_2) to share Schaaps-54_uns1
unk.sum2$uns_male[which(unk.sum2$clutch_id_ind2=="Schaaps_97_2")] <- "Schaaps-54_uns1"

# update Cooks-27 second clutch (Cooks_29_2) with uns-2
unk.sum2$uns_male[which(unk.sum2$clutch_id_ind2=="Cooks_29_2")] <- "Cooks-27_uns2"


# sum offspring by unsampled male ID
uns_male_RS <- unk.sum2 %>% group_by(uns_male) %>%
  summarise(n_chick = sum(n))

# save table
write.csv(uns_male_RS, "03_output-files/unsampled male RS 2021.csv", row.names=F)


#-------------------------------------------------------------------------------
# update assigned table with uns_male IDs
#-------------------------------------------------------------------------------

# add uns_male IDs to table where each row is a kid
# need to manually add the cases where there were two sires within a single clutch
# Cooks-23, Schaaps-71, and Schaaps-80
# Add everyone except these uns2 sires
dad.unk2 <- left_join(dad.unk, unk.sum2[1:15,c(1,4)], by=c("clutch_id_ind2"))

# Cooks-23 kid 82971 is half sib to others in this clutch. They actually get
# the sire from Cooks_31_collected
dad.unk2$uns_male[which(dad.unk2$Ind2=="CO_82971")] <- "Cooks-31_uns1"

# For Schaaps-71,  83003 is the one with a different dad, uns2
dad.unk2$uns_male[which(dad.unk2$Ind2=="CO_83003")] <- "Schaaps-71_uns2"

# For Schaaps-80, arbitrarily give CO_83046 uns2
dad.unk2$uns_male[which(dad.unk2$Ind2=="CO_83046")] <- "Schaaps-80_uns2"



# update Band_dad 
dad.unk2$Band_dad <- dad.unk2$uns_male

# update genetic family
dad.unk2$genetic_fam <- paste(dad.unk2$Band_dad, dad.unk2$Band_mom, sep="_")

## bind updated unk with rest of assigned
assigned.notUNS <- subset(assigned, !is.na(assigned$Band_dad))

assigned.with.uns <- rbind(assigned.notUNS, dad.unk2[,-29])

# save table
write.csv(assigned.with.uns, 
          "03_output-files/kin_2021_parent_offspring_assigned_withUNS.csv",
          row.names=F)





