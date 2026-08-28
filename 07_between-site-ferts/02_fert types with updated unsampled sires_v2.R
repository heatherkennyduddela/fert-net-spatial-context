
################################################################################
# Update unsampled sires for all years to include known between-site ferts
# at sites where all nesting males were banded. Where unbanded males existed, 
# assign these to be either the same as the unsampled males or different males,
# to calculate the maximum and minimum possible number of different-structure
# fertilizations. 

# Heather Kenny-Duddela
# May 25, 2026
################################################################################

# Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr) # multipanel plots

# load files

# updated parentage assignment with labeled unsampled males for each year
# From the scripts calculating the number of unsampled sires for each year
po21 <- read.csv("input-files/kin_2021_parent_offspring_assigned_withUNS.csv")
po22 <- read.csv("input-files/kin_2022_parent_offspring_assigned_withUNS.csv")
po23 <- read.csv('input-files/kin_2023_parent_offspring_assigned_withUNS.csv')

# all adults analyzed for paternity, from paternity assignment with lcMLkin scripts 
# for each year
kin21 <- read.csv("input-files/kin_2021_adults_all.csv")
kin22 <- read.csv("input-files/kin_2022_adults_all.csv")
kin23 <- read.csv("input-files/kin_2023_adults_all.csv")

#-------------------------------------------------------------------------------
# sum total number of females and males sampled for paternity
#-------------------------------------------------------------------------------

colnames(kin23)[11] <- "Type_ind1"
colnames(kin23)[18] <- "Type_ind2"


ind1.all <- rbind(select(kin21, Ind1, Type_ind1),
                  select(kin22, Ind1, Type_ind1),
                  select(kin23, Ind1, Type_ind1))

ind2.all <- rbind(select(kin21, Ind2, Type_ind2),
                  select(kin22, Ind2, Type_ind2),
                  select(kin23, Ind2, Type_ind2))

colnames(ind1.all) <- c("Ind","Type")
colnames(ind2.all) <- c("Ind","Type")

all.adult <- rbind(ind1.all, ind2.all)

adult.tally <- all.adult %>%
  group_by(Type) %>%
  summarise(n.adult = length(unique(Ind)))



#-------------------------------------------------------------------------------
# Add fert types, including ep_unsamp_diff where possible
#-------------------------------------------------------------------------------

# 2021 -------------------------------------------------------------------------

# summarise at clutch level
shared.fert.clutch21 <- po21 %>% 
  group_by(clutch_id_ind2, Band_dad, Band_mom, 
           Site_dad, Site_mom, FamilyID_dad, FamilyID_mom, Brood_ind2, 
           site_cat) %>%
  summarise(fert = n())

# number of clutches (56)
length(unique(shared.fert.clutch21$clutch_id_ind2))

shared.fert.clutch21$fert_type <- NA
shared.fert.clutch21$fert_type[which(shared.fert.clutch21$FamilyID_dad == 
                                     shared.fert.clutch21$FamilyID_mom)] <- "wp"
shared.fert.clutch21$fert_type[which(shared.fert.clutch21$FamilyID_dad != 
                                     shared.fert.clutch21$FamilyID_mom &
                                     shared.fert.clutch21$Site_dad == 
                                     shared.fert.clutch21$Site_mom)] <- "ep_same"
shared.fert.clutch21$fert_type[which(shared.fert.clutch21$FamilyID_dad != 
                                     shared.fert.clutch21$FamilyID_mom &
                                     shared.fert.clutch21$Site_dad != 
                                     shared.fert.clutch21$Site_mom)] <- "ep_diff"
shared.fert.clutch21$fert_type[which(is.na(shared.fert.clutch21$FamilyID_dad))] <- "ep_unsampled"



# update shared.fert.clutch to re-classify ferts at the same site but 
# between different structures
# Blue Cloud-19 to Blue Cloud-14
# Cooks-23 to Hepp

shared.fert.clutch21.update <- shared.fert.clutch21

shared.fert.clutch21.update$fert_type[which(
  shared.fert.clutch21.update$clutch_id_ind2=="BlueCloud-14_collected" &
    shared.fert.clutch21.update$FamilyID_dad=="BlueCloud-19")] <- "ep_diff"

shared.fert.clutch21.update$fert_type[which(
  shared.fert.clutch21.update$clutch_id_ind2=="Hepp (near Cook)_1_1" &
    shared.fert.clutch21.update$FamilyID_dad=="Cooks-23")] <- "ep_diff"

# add new category for sites where we banded all known breeding males, and
# unsampled sires are therefore from outside the site: 
# Cathys/Mary Anns, Cooks, Hepp, Make Believe, Struthers, Reinarz 

# Only Cooks has fert originally classified as ep_unsampled
shared.fert.clutch21.update$fert_type[which(
  shared.fert.clutch21.update$Site_mom=="Cooks" &
    shared.fert.clutch21.update$fert_type=="ep_unsampled")] <- "ep_unsamp_diff"



# 2022 -------------------------------------------------------------------------

# summarise at clutch level
shared.fert.clutch22 <- po22 %>% 
  group_by(clutch_id_ind2, Band_dad, Band_mom, 
           Site_dad, Site_mom, FamilyID_dad, FamilyID_mom, brood, 
           site_type) %>%
  summarise(fert = n())

# number of clutches (92)
length(unique(shared.fert.clutch22$clutch_id_ind2))

shared.fert.clutch22$fert_type <- NA
shared.fert.clutch22$fert_type[which(shared.fert.clutch22$FamilyID_dad == 
                                       shared.fert.clutch22$FamilyID_mom)] <- "wp"
shared.fert.clutch22$fert_type[which(shared.fert.clutch22$FamilyID_dad != 
                                       shared.fert.clutch22$FamilyID_mom &
                                       shared.fert.clutch22$Site_dad == 
                                       shared.fert.clutch22$Site_mom)] <- "ep_same"
shared.fert.clutch22$fert_type[which(shared.fert.clutch22$FamilyID_dad != 
                                       shared.fert.clutch22$FamilyID_mom &
                                       shared.fert.clutch22$Site_dad != 
                                       shared.fert.clutch22$Site_mom)] <- "ep_diff"
shared.fert.clutch22$fert_type[which(is.na(shared.fert.clutch22$FamilyID_dad))] <- "ep_unsampled"

# update shared.fert.clutch to re-classify ferts at the same site but 
# between different structures
# For 2022, there were active nests in the following outbuildings
# Blue cloud middle shed (19), wash room (24), gazebo (26), and trailer (27)
# Make Believe shed (62)
# Urban Farm sheep shed (2,5), Urban Farm Red Barn (4, 7)

## Check if any of these birds were involved in EP ferts 
# visual inspection of ep.all.22, Blue Cloud is the only one

# BlueCloud-19 dad to BlueCloud-18 mom, clutch 21_2
# Pair from 19 moved to wash room nest (24) after the shed nest failed. They 
# had already initiated a nest at 24 by the time 21 started laying. Distance 
# for this EP should be wash room to main  barn

shared.fert.clutch22.update <- shared.fert.clutch22
shared.fert.clutch22.update$fert_type[which(shared.fert.clutch22.update$FamilyID_dad=="BlueCloud-19" &
                              shared.fert.clutch22.update$FamilyID_mom=="BlueCloud-18")] <- "ep_diff"

# Determine which sites had all males banded, and therefore unk ferts should be 
# changed to ep_unsamp_diff

# Sites with unk_ep were, with 100% banding listed after dash: 
# Cooks - yes (one UNB female from paternity, all males banded)
# CHR - no (4 UNB males)
# Dome House - yes (solitary)
# Make Believe - yes
# Speedwell - yes (solitary)
# Struthers - yes
# Urban Farm Girlz - yes

# Change dad_unk at these sites to "ep_unsamp_diff"
shared.fert.clutch22.update$fert_type[which(is.na(shared.fert.clutch22.update$Site_dad) &
                              shared.fert.clutch22.update$Site_mom == "Cooks")] <- "ep_unsamp_diff"
shared.fert.clutch22.update$fert_type[which(is.na(shared.fert.clutch22.update$Site_dad) &
                              shared.fert.clutch22.update$Site_mom == "Dome House")] <- "ep_unsamp_diff"
shared.fert.clutch22.update$fert_type[which(is.na(shared.fert.clutch22.update$Site_dad) &
                              shared.fert.clutch22.update$Site_mom == "Make Believe")] <- "ep_unsamp_diff"
shared.fert.clutch22.update$fert_type[which(is.na(shared.fert.clutch22.update$Site_dad) &
                              shared.fert.clutch22.update$Site_mom == "Speedwell")] <- "ep_unsamp_diff"
shared.fert.clutch22.update$fert_type[which(is.na(shared.fert.clutch22.update$Site_dad) &
                              shared.fert.clutch22.update$Site_mom == "Struthers")] <- "ep_unsamp_diff"
shared.fert.clutch22.update$fert_type[which(is.na(shared.fert.clutch22.update$Site_dad) &
                              shared.fert.clutch22.update$Site_mom == "Urban Farm Girlz")] <- "ep_unsamp_diff"
shared.fert.clutch22.update$fert_type[which(is.na(shared.fert.clutch22.update$Site_dad) &
                                              shared.fert.clutch22.update$Site_mom == "Cooks")] <- "ep_unsamp_diff"


# 2023 -------------------------------------------------------------------------

# summarise at clutch level
shared.fert.clutch23 <- po23 %>% 
  group_by(clutchID_ind2, band_dad, band_mom, 
           site_dad, site_mom, FamilyID_dad, FamilyID_mom, brood_ind2, 
           site_type) %>%
  summarise(fert = n())

# number of clutches (60)
length(unique(shared.fert.clutch23$clutchID_ind2))

shared.fert.clutch23$fert_type <- NA
shared.fert.clutch23$fert_type[which(shared.fert.clutch23$FamilyID_dad == 
                                       shared.fert.clutch23$FamilyID_mom)] <- "wp"
shared.fert.clutch23$fert_type[which(shared.fert.clutch23$FamilyID_dad != 
                                       shared.fert.clutch23$FamilyID_mom &
                                       shared.fert.clutch23$site_dad == 
                                       shared.fert.clutch23$site_mom)] <- "ep_same"
shared.fert.clutch23$fert_type[which(shared.fert.clutch23$FamilyID_dad != 
                                       shared.fert.clutch23$FamilyID_mom &
                                       shared.fert.clutch23$site_dad != 
                                       shared.fert.clutch23$site_mom)] <- "ep_diff"
shared.fert.clutch23$fert_type[which(is.na(shared.fert.clutch23$FamilyID_dad))] <- "ep_unsampled"


# update shared.fert.clutch to re-classify ferts at the same site but 
# between different structures
# For 2023, the following outbuilding nests were active and sampled for paternity: 
# Blue Cloud: 28 (north shed) and 19 (middle shed) were active but not sampled
# Cooks: 15 (house) was active but not sampled
# Make Believe: no shed nests were active
#--------------------------
# Struthers: 23 (round shed)
# Urban Farm Girlz: 8, 5, 2 (sheep shed); 4, 7 (Red barn)

# Change Urban Farm 7 to Urban Farm 8 to ep_diff
shared.fert.clutch23.update <- shared.fert.clutch23
shared.fert.clutch23.update$fert_type[which(shared.fert.clutch23.update$FamilyID_dad=="UrbanFarm-07" &
                              shared.fert.clutch23.update$FamilyID_mom=="UrbanFarm-08")] <- "ep_unsamp_diff"

# Identify unsampled sires at sites with 100% of males banded-------------------

# For 2023, sites with 100% of adult males banded: 
# Blue Cloud - no
# *Boyer - yes
# *Cathys - yes (solitary)
# Cooks - no
# *Jays - yes (solitary)
# Make Beleive - no
# *Nixon - yes (solitary)
# *Plumbumpy - yes (solitary)
# *Red Wagon - yes (solitary)
# *Speedwell - yes (solitary)
# *Struthers - yes
# Urban Farm Girlz - no (?)

# change Boyer, Jays, and Speedwell unsampled 
shared.fert.clutch23.update$fert_type[which(shared.fert.clutch23.update$site_mom=="Boyer" &
                              is.na(shared.fert.clutch23.update$site_dad)|
                              shared.fert.clutch23.update$site_mom=="Jays" &
                              is.na(shared.fert.clutch23.update$site_dad)|
                              shared.fert.clutch23.update$site_mom=="Speedwell"&
                              is.na(shared.fert.clutch23.update$site_dad))] <- "ep_unsamp_diff"

# combine years ----------------------------------------------------------------

# harmonize names
colnames(shared.fert.clutch21.update)
colnames(shared.fert.clutch22.update)
colnames(shared.fert.clutch23.update)

colnames(shared.fert.clutch22.update) <- colnames(shared.fert.clutch21.update)
colnames(shared.fert.clutch23.update) <- colnames(shared.fert.clutch21.update)

# add year identifier
shared.fert.clutch21.update$year <- 2021
shared.fert.clutch22.update$year <- 2022
shared.fert.clutch23.update$year <- 2023

# harmonize object class for brood
shared.fert.clutch21.update$Brood_ind2 <- as.character(shared.fert.clutch21.update$Brood_ind2)
shared.fert.clutch22.update$Brood_ind2 <- as.character(shared.fert.clutch22.update$Brood_ind2)
shared.fert.clutch23.update$Brood_ind2 <- as.character(shared.fert.clutch23.update$Brood_ind2)

# # combine tables
# all.years.clutch <- rbind(shared.fert.clutch21.update, 
#                           shared.fert.clutch22.update,
#                           shared.fert.clutch23.update)
# 
# # update fert type for CHR-101 and Cooks-13 which had unsampled moms
# # For CHR_101_2 the dad familyID matches with the mom -> wp
# # For Cooks_13_1 the dad familyID matches with the mom -> wp
# all.years.clutch$fert_type[which(is.na(all.years.clutch$Band_mom))] <- "wp"
# 
# # save table
# write.csv(all.years.clutch, 
#           file="02_output-files/fert type by clutch unsamp_diff update.csv",
#           row.names = F)



#-------------------------------------------------------------------------------
# Identify males that did not have their social nest kids included in the
# paternity analysis for each year
#-------------------------------------------------------------------------------


### Identify males from each year with zero fert

## 2021

male.21 <- subset(kin21, kin21$Type_ind1=="dad" & kin21$Type_ind2=="dad")
male.21.unique <- unique(male.21[,8:11]) # 37
# figure out number that had social nests monitored
sum(male.21.unique$FamilyID_ind1 != "none") # 29 with nests monitored

# List of dads that did not have social nest kids included in paternity analysis
male.21.nokid <- male.21.unique[
  which(male.21.unique$FamilyID_ind1 %in% po21$FamilyID_ind2 == F), ]

# List of dads with zero fertilizations
nofert.21 <- male.21.unique[
  which(male.21.unique$Band_ind1 %in% po21$Band_dad == F), ]

# for these males with familyIDs, check if any offspring were included in
# the paternity analysis. Only candidate is MakeBelieve-61
# Searched in po21 for MakeBelieve-61 and there were no kids. 

# In 2021, none of the monitored males with social offspring included in the 
# paternity analysis (27 of the 29 monitored) had zero success


## 2022

male.22 <- subset(kin22, kin22$Type_ind1=="dad" & kin22$Type_ind2=="dad")
male.22.unique <- unique(male.22[,8:11]) # 104
# figure out number that had social nests monitored
sum(male.22.unique$FamilyID_ind1 != "none") # 73 with nests monitored

# List of dads that did not have social nest kids included in paternity analysis
male.22.nokid <- male.22.unique[
  which(male.22.unique$FamilyID_ind1 %in% po22$FamilyID_ind2 == F), ]

# List of dads with zero ferts
nofert.22 <- male.22.unique[
  which(male.22.unique$Band_ind1 %in% po22$Band_dad == F), ]

# for these males with familyIDs, check if any offspring from his social nest
# were included in the paternity analysis.

# BlueCloud-26, no offspring analyzed
# Boyer-11, no offspring analyzed
# CHR-058, 4 analyzed
# CHR-082, no offspring analyzed
# CHR-085, no offspring analyzed
# CHR-149, no offspring analyzed
# Cooks-06, 5 analyzed
# Jays-02, no offspring analyzed
# Ericas-01, no offspring analyzed
# MakeBelieve-62, no offspring analyzed


## 2023
male.23 <- subset(kin23, kin23$type_ind1=="dad" & kin23$type_ind2=="dad")
male.23.unique <- unique(male.23[,8:12]) # 73
# figure out number that had social nests monitored
sum(male.23.unique$FamilyID_ind1 != "none") # 37 with nests monitored

# List of dads that did not have social nest kids included in paternity analysis
male.23.nokid <- male.23.unique[
  which(male.23.unique$FamilyID_ind1 %in% po23$FamilyID_ind2 == F), ]
# Remove Boyer-15 dad because he did have  social nests included. He wasn't 
# confirmed at nest 13 originally but he was later banded and the paternity
# makes sense for him to be the social dad for family Boyer-13. 
male.23.nokid2 <- subset(male.23.nokid, male.23.nokid$FamilyID_ind1 != "Boyer-15")

# List of dads with zero ferts
nofert.23 <- male.23.unique[
  which(male.23.unique$band_ind1 %in% po23$band_dad == F), ]

# for these males with familyIDs, check if any offspring from his social nest
# were included in the paternity analysis.

# Cooks-18, IDobs at this nest unconfirmed. A male was caught off the nest 
# at night but he is not the genetic dad of any of the collected eggs
# Cooks-4.2, 3 chicks analyzed
# Make Believe-17, 2870-65018 was seen near 17 after eggs collected and ID not 
# confirmed. He is not the dad of any of the eggs
# Make Believe-34, 2870-65017 caught with 2850-57938 who is the mom of nest 34, 
# but IDobs were not confirmed. He is not the dad of any collected eggs
# UrbanFarm-07, 2850-57976 was caught at nest 7 but ID never confirmed. 4 chicks
# analyzed from nest 7 and dad of all was actually male thought to be at
# nest 5 (2850-57976) and he also sired offspring in nest 8.2...

# WP and EP same-site fert assignment are too uncertain at Urban Farm Girlz 2023 
# so remove these from the analysis of fert types! Can still keep in the between-site
# distances. 

# Will also remove Cooks-18, MakeBelieve-17, and MakeBelieve-34 from the analysis
# of fert types because again social dad IDs are too uncertain. Will consider
# social males at these nests to be unsampled (UNB) for the purposes of determining
# unsampled sire fert types (same site or diff site)

#-------------------------------------------------------------------------------
# Save list of dads from all years with no social kids analyzed
#-------------------------------------------------------------------------------

# add year label
male.21.nokid$year <- 2021
male.22.nokid$year <- 2022
male.23.nokid$nest_ind1 <- 2023
colnames(male.23.nokid) <- c("Band_ind1", "Site_ind1", "year", "Type_ind1",
                             "FamilyID_ind1")

# combine tables
male.all.nokid <- rbind(male.21.nokid, male.22.nokid, male.23.nokid)

#-------------------------------#
# Table of dads with uncertain IDobs that should NOT be analyzed for fert types
male.uncertain <- subset(male.23.unique, 
                         male.23.unique$FamilyID_ind1 == "Cooks-18" |
                           male.23.unique$FamilyID_ind1 == "MakeBelieve-17" |
                           male.23.unique$FamilyID_ind1 == "MakeBelieve-34" |
                           male.23.unique$site_ind1 == "Urban Farm Girlz")
# add year label
male.uncertain$year <- 2023

# Table of moms with no kids analyzed for paternity ----------------------------

# 2021

fem.21 <- subset(kin21, kin21$Type_ind1 =="mom" & kin21$Type_ind2=="mom")
fem.21.unique <- unique(fem.21[,8:12]) 

# List of moms that did not have social nest kids included in paternity analysis
fem.21.nokid <- fem.21.unique[
  which(fem.21.unique$FamilyID_ind1 %in% po21$FamilyID_ind2 == F), ]
# none for 2021


# 2022

fem.22 <- subset(kin22, kin22$Type_ind1 =="mom" & kin22$Type_ind2=="mom")
fem.22.unique <- unique(fem.22[,8:12]) 

# List of moms that did not have social nest kids included in paternity analysis
fem.22.nokid <- fem.22.unique[
  which(fem.22.unique$FamilyID_ind1 %in% po22$FamilyID_ind2 == F), ]


# 2023

fem.23 <- subset(kin23, kin23$type_ind1 =="mom" & kin23$type_ind2=="mom")
fem.23.unique <- unique(fem.23[,8:12]) 

# List of moms that did not have social nest kids included in paternity analysis
fem.23.nokid <- fem.23.unique[
  which(fem.23.unique$FamilyID_ind1 %in% po23$FamilyID_ind2 == F), ]
# exclude female from Boyer-15 who was not confirmed for nest 13 but is the mom
# just don't add this table to the others


# Combine tables, but 2022 is the only year with legit females to check
# just add year label
fem.22.nokid$Nest.._ind1 <- 2022
colnames(fem.22.nokid)[5] <- "year"


# add females to the male table
all.nokid <- rbind(male.all.nokid, fem.22.nokid)

# also add females to the uncertain table for Urban Farm, Cooks-18, MakeBelieve-17, and MakeBelieve-34
fem.uncertain <- subset(fem.23.unique, 
                        fem.23.unique$FamilyID_ind1 == "Cooks-18" |
                          fem.23.unique$FamilyID_ind1 == "MakeBelieve-17" |
                          fem.23.unique$FamilyID_ind1 == "MakeBelieve-34" |
                          fem.23.unique$site_ind1 == "Urban Farm Girlz")
# add year label
fem.uncertain$year <- 2023

# combine
all.uncertain <- rbind(male.uncertain, fem.uncertain)

# save both tables
write.csv(all.nokid, "02_output-files/birds with no kids analyzed all years.csv",
          row.names=F)

write.csv(all.uncertain, "02_output-files/birds with uncertain IDobs.csv",
          row.names=F)

#-------------------------------------------------------------------------------
# Known UNB males at each site
#-------------------------------------------------------------------------------

# Blue Cloud 2021 - 1 UNB male, nest 17 none analyzed
# Blue Cloud 2022 - 1 UNB male, nest 27 none analyzed
# Blue Cloud 2023 - 2 UNB male (nest 22 and 19), 3 eggs analyzed from nest 22

# CHR 2022 - 4 UNB male
# - CHR 7 and 9, female 2850-57584. None analyzed
# - CHR 77, female 2850-57696. None analyzed
# - CHR 76 and 12, female 2850-57737. None analyzed
# - CHR-145, female 2870-65004. None analyzed

# Cooks 2023 - potentially 3 UNB males: nest 18, and male ID at 2.2 was unconfirmed,
# and another UNB seen but not tied to a nest. 
# Nest 18 had collected clutch analyzed, 2.2 and other UNB did not have offspring analyzed

# Make Believe 2023 - up to 2 UNB males (nest 34 and 17), both had clutches analyzed
# One egg analyzed for 17, 3 egg analyzed for 34

# Schaaps 2021 - 2 UNB male (nest 131 and 75), none analyzed

# Urban Farm Girlz 2021 - 1 UNB male, none analyzed

##################### Other notes ############################
# Urban Farm Girlz 2023 - sparse IDbos, estimate 2 UNB males (nests 8 and 4), 
# one with his social nest analyzed


#-------------------------------------------------------------------------------
# Number of unsample sires compared to the number of known UNB males by site
#-------------------------------------------------------------------------------

# 2021

# create table with number of unsampled sires for each site and year
clutch21.unsamp  <- subset(shared.fert.clutch21, 
                           is.na(shared.fert.clutch21$Site_dad))

# add site label
clutch21.unsamp$Site_unsamp <- sub("-.*", "", clutch21.unsamp$Band_dad)


unsamp21 <- clutch21.unsamp %>% group_by(Site_unsamp) %>%
  summarise(n.clutch = n(),
            n.sire = length(unique(Band_dad)),
            n.fert = sum(fert))

# Manually add number of unbanded
unsamp21$n.unb <- c(1, 0, 2, 1)

# add year label
unsamp21$year <- 2021


# 2022

# create table with number of unsampled sires for each site and year
clutch22.unsamp  <- subset(shared.fert.clutch22, 
                           is.na(shared.fert.clutch22$Site_dad))

# add site label
clutch22.unsamp$Site_unsamp <- sub("-.*", "", clutch22.unsamp$Band_dad)


unsamp22 <- clutch22.unsamp %>% group_by(Site_unsamp) %>%
  summarise(n.clutch = n(),
            n.sire = length(unique(Band_dad)),
            n.fert = sum(fert))

# Manually add number of unbanded
unsamp22$n.unb <- c(4, 0, 0, 0, 0, 0, 0)

# add year label
unsamp22$year <- 2022


# 2023

# create table with number of unsampled sires for each site and year
clutch23.unsamp  <- subset(shared.fert.clutch23, 
                           is.na(shared.fert.clutch23$site_dad))

# add site label
clutch23.unsamp$Site_unsamp <- sub("-.*", "", clutch23.unsamp$band_dad)


unsamp23 <- clutch23.unsamp %>% group_by(Site_unsamp) %>%
  summarise(n.clutch = n(),
            n.sire = length(unique(band_dad)),
            n.fert = sum(fert))

# Manually add number of unbanded
unsamp23$n.unb <- c(2, 0, 3, 2, 0, 2)

# add year label
unsamp23$year <- 2023


# combine tables
unsamp.all <- rbind(unsamp21, unsamp22, unsamp23)

# subset only those with known UNB males
unsamp.unb <- subset(unsamp.all, unsamp.all$n.unb > 0)


# save both tables
write.csv(unsamp.all, 
          "02_output-files/unsampled and unbanded sires all years.csv",
          row.names = F)

write.csv(unsamp.unb, 
          "02_output-files/unsampled and unbanded sires with UNB.csv",
          row.names=F)

#-------------------------------------------------------------------------------
# Use number of UNB and number of unsampled sires to calculate the max and min
# number of diff-site ferts
#-------------------------------------------------------------------------------

# For minimizing different-site ferts, the male(s) with the most fertilizations
# should be UNB males from that site. 

# summarise fertilizations by male and site

clutch21.unsamp$year <- 2021
clutch22.unsamp$year <- 2022
clutch23.unsamp$year <- 2023
colnames(clutch23.unsamp)[2] <- "Band_dad"

clutch.unsamp.all <- rbind(clutch21.unsamp, clutch22.unsamp, clutch23.unsamp)

fert.by.unsamp <- clutch.unsamp.all %>% group_by(year, Site_unsamp, Band_dad) %>%
  summarise(n.fert = sum(fert))

# order by number of fert in descending order
fert.by.unsamp.sort <- sort_by(fert.by.unsamp, 
                                 list(fert.by.unsamp$year,
                                      fert.by.unsamp$Site_unsamp,
                                      -fert.by.unsamp$n.fert))

# keep only rows where site had UNB

# add site by year label
fert.by.unsamp.sort$site.year <- paste(fert.by.unsamp.sort$Site_unsamp, 
                                       fert.by.unsamp.sort$year, 
                                       sep="_")
unsamp.unb$site.year <- paste(unsamp.unb$Site_unsamp,
                              unsamp.unb$year,
                              sep="_")

# summarize at the year level
uns.unb.year <- unsamp.unb %>% group_by(year) %>%
  summarise(tot.sire = sum(n.sire),
            tot.unb = sum(n.unb))

# subset
fert.by.unsamp.sort.unb <- subset(fert.by.unsamp.sort, 
                                  fert.by.unsamp.sort$site.year %in%
                                    unsamp.unb$site.year)

# Min diff ---------------------------------------------------------------------

clutch21.min.diff <- shared.fert.clutch21.update
clutch22.min.diff <- shared.fert.clutch22.update
clutch23.min.diff <- shared.fert.clutch23.update

# Loop through each site and assign UNB to sires -------------------------------

for (i in 1:length(unique(unsamp.unb$site.year))) {
  # use different data for each year
  
  # 2021 ---------------------------#
  if (unsamp.unb$year[i] == 2021) {
    
    # get number of unbanded for that site in that year
    n.unb <- unsamp.unb$n.unb[i]
    
    # subset unsampled sires for that site and year
    sires <- subset(fert.by.unsamp.sort.unb, 
                    fert.by.unsamp.sort.unb$site.year==
                      unsamp.unb$site.year[i])
    
    # take the first n.unb number of rows from the sires table
    sires.unb <- sires[1:n.unb, ]
    
    # loop through each of these sires
    for (j in 1:length(sires.unb$Band_dad)) {
      
      # change fert type label for these sires to "ep_same"
      clutch21.min.diff$fert_type[which(
        clutch21.min.diff$Band_dad==sires.unb$Band_dad[j])] <- "ep_same"
    }
    
    # take rest of unsampled sires (if any)
    sires.unsamp <- sires[-(1:n.unb), ]
    
    # loop through each of the remaining unsampled sires and 
    # change fert type label to "ep_diff"
    for (k in 1:length(sires.unsamp$Band_dad)) {
      
      clutch21.min.diff$fert_type[which(
        clutch21.min.diff$Band_dad==sires.unsamp$Band_dad[k])] <- "ep_diff"
    }
    
    # change any remaining "ep_unsamp_diff" to "ep_diff"
    clutch21.min.diff$fert_type[which(
      clutch21.min.diff$fert_type=="ep_unsamp_diff")] <- "ep_diff"
    
  }
  
  
  # 2022 ---------------------------#
  if (unsamp.unb$year[i] == 2022) {
    
    # get number of unbanded for that site in that year
    n.unb <- unsamp.unb$n.unb[i]
    
    # subset unsampled sires for that site and year
    sires <- subset(fert.by.unsamp.sort.unb, 
                    fert.by.unsamp.sort.unb$site.year==
                      unsamp.unb$site.year[i])
    
    # take the first n.unb number of rows from the sires table
    sires.unb <- sires[1:n.unb, ]
    
    # loop through each of these sires
    for (j in 1:length(sires.unb$Band_dad)) {
      
      # change fert type label for these sires to "ep_same"
      clutch22.min.diff$fert_type[which(
        clutch22.min.diff$Band_dad==sires.unb$Band_dad[j])] <- "ep_same"
    }
    
    # take rest of unsampled sires (if any)
    sires.unsamp <- sires[-(1:n.unb), ]
    
    # loop through each of the remaining unsampled sires and 
    # change fert type label to "ep_diff"
    for (k in 1:length(sires.unsamp$Band_dad)) {
      
      clutch22.min.diff$fert_type[which(
        clutch22.min.diff$Band_dad==sires.unsamp$Band_dad[k])] <- "ep_diff"
    }
    
    # change any remaining "ep_unsamp_diff" to "ep_diff"
    clutch22.min.diff$fert_type[which(
      clutch22.min.diff$fert_type=="ep_unsamp_diff")] <- "ep_diff"
    
  }
  
  
  # 2023 ---------------------------#
  if (unsamp.unb$year[i] == 2023) {
    
    # get number of unbanded for that site in that year
    n.unb <- unsamp.unb$n.unb[i]
    
    # subset unsampled sires for that site and year
    sires <- subset(fert.by.unsamp.sort.unb, 
                    fert.by.unsamp.sort.unb$site.year==
                      unsamp.unb$site.year[i])
    
    # take the first n.unb number of rows from the sires table
    sires.unb <- sires[1:n.unb, ]
    
    # loop through each of these sires
    for (j in 1:length(sires.unb$Band_dad)) {
      
      # change fert type label for these sires to "ep_same"
      clutch23.min.diff$fert_type[which(
        clutch23.min.diff$Band_dad==sires.unb$Band_dad[j])] <- "ep_same"
    }
    
    # take rest of unsampled sires (if any)
    sires.unsamp <- sires[-(1:n.unb), ]
    
    # loop through each of the remaining unsampled sires and 
    # change fert type label to "ep_diff"
    for (k in 1:length(sires.unsamp$Band_dad)) {
      
      clutch23.min.diff$fert_type[which(
        clutch23.min.diff$Band_dad==sires.unsamp$Band_dad[k])] <- "ep_diff"
    }
    
    # change any remaining "ep_unsamp_diff" to "ep_diff"
    clutch23.min.diff$fert_type[which(
      clutch23.min.diff$fert_type=="ep_unsamp_diff")] <- "ep_diff"
    
  }
}

# For the unsampled male that had ferts at both Boyer and UrbanFarmGirlz,
# change the fert type to ep_diff for the ferts at Boyer
clutch23.min.diff$fert_type[which(
  clutch23.min.diff$clutch_id_ind2=="Boyer_3.2_1" &
    is.na(clutch23.min.diff$Site_dad))] <- "ep_diff"


# Max diff ---------------------------------------------------------------------

# For maximizing different-site ferts, all unsampled sires should be 
# different site

clutch21.max.diff <- shared.fert.clutch21.update
clutch22.max.diff <- shared.fert.clutch22.update
clutch23.max.diff <- shared.fert.clutch23.update

# change fert labels
clutch21.max.diff$fert_type[which(clutch21.max.diff$Band_dad %in%
                                    fert.by.unsamp$Band_dad)] <- "ep_diff"

clutch22.max.diff$fert_type[which(clutch22.max.diff$Band_dad %in%
                                    fert.by.unsamp$Band_dad)] <- "ep_diff"

clutch23.max.diff$fert_type[which(clutch23.max.diff$Band_dad %in%
                                    fert.by.unsamp$Band_dad)] <- "ep_diff"

# make bar plots for min.diff and max.diff -------------------------------------

# combine years
clutch.all.min.diff <- rbind(clutch21.min.diff, clutch22.min.diff, 
                             clutch23.min.diff)

clutch.all.max.diff <- rbind(clutch21.max.diff, clutch22.max.diff, 
                             clutch23.max.diff)

# correct NA labels for nests with unsampled moms
# CHR 101 social dad is 2850-57733 => WP
# Cooks 13 social dad is 2850-57675 => wp
clutch.all.min.diff$fert_type[which(
  clutch.all.min.diff$clutch_id_ind2=="CHR_101_2" |
    clutch.all.min.diff$clutch_id_ind2=="Cook's_13_1" )] <- "wp"

clutch.all.max.diff$fert_type[which(
  clutch.all.max.diff$clutch_id_ind2=="CHR_101_2" |
    clutch.all.max.diff$clutch_id_ind2=="Cook's_13_1" )] <- "wp"

# Remove birds identified above that should not be included because there
# was uncertainty for the IDobs and identity of social dads

# WP and EP same-site fert assignment are too uncertain at Urban Farm Girlz 2023 
# so remove these from the analysis of fert types! Can still keep in the between-site
# distances. 
# Will also remove Cooks-18, MakeBelieve-17, and MakeBelieve-34 from the analysis
# of fert types because again social dad IDs are too uncertain in 2023. 

clutch.all.min.diff2 <- subset(
  clutch.all.min.diff, 
    clutch.all.min.diff$clutch_id_ind2 != "Make Believe_17_collected" &
    clutch.all.min.diff$clutch_id_ind2 != "Make Believe_34_collected" &
    clutch.all.min.diff$clutch_id_ind2 != "Cooks_18_collected" &
    clutch.all.min.diff$clutch_id_ind2 != "Urban Farm Girlz_5_1" &
    clutch.all.min.diff$clutch_id_ind2 != "Urban Farm Girlz_7_2" &
    clutch.all.min.diff$clutch_id_ind2 != "Urban Farm Girlz_8.2_2" &
    clutch.all.min.diff$clutch_id_ind2 != "Urban Farm Girlz_8_1" )

clutch.all.max.diff2 <- subset(
  clutch.all.max.diff, 
  clutch.all.max.diff$clutch_id_ind2 != "Make Believe_17_collected" &
    clutch.all.max.diff$clutch_id_ind2 != "Make Believe_34_collected" &
    clutch.all.max.diff$clutch_id_ind2 != "Cooks_18_collected" &
    clutch.all.max.diff$clutch_id_ind2 != "Urban Farm Girlz_5_1" &
    clutch.all.max.diff$clutch_id_ind2 != "Urban Farm Girlz_7_2" &
    clutch.all.max.diff$clutch_id_ind2 != "Urban Farm Girlz_8.2_2" &
    clutch.all.max.diff$clutch_id_ind2 != "Urban Farm Girlz_8_1" )

# also make plots for "raw" fert types with unsampled as a fert category
clutch.all.update <- rbind(shared.fert.clutch21.update,
                           shared.fert.clutch22.update,
                           shared.fert.clutch23.update)

clutch.all.update2 <- subset(
  clutch.all.update, 
  clutch.all.update$clutch_id_ind2 != "Make Believe_17_collected" &
    clutch.all.update$clutch_id_ind2 != "Make Believe_34_collected" &
    clutch.all.update$clutch_id_ind2 != "Cooks_18_collected" &
    clutch.all.update$clutch_id_ind2 != "Urban Farm Girlz_5_1" &
    clutch.all.update$clutch_id_ind2 != "Urban Farm Girlz_7_2" &
    clutch.all.update$clutch_id_ind2 != "Urban Farm Girlz_8.2_2" &
    clutch.all.update$clutch_id_ind2 != "Urban Farm Girlz_8_1" )

clutch.all.update2$fert_type[which(
  clutch.all.update2$clutch_id_ind2=="CHR_101_2" |
    clutch.all.update2$clutch_id_ind2=="Cook's_13_1" )] <- "wp"

# change "ep_unsamp_diff" to "ep_diff"
clutch.all.update2$fert_type[which(
  clutch.all.update2$fert_type=="ep_unsamp_diff")] <- "ep_diff"



# save tables
write.csv(clutch.all.min.diff2,
          "02_output-files/fert type min diff site.csv", row.names = F)

write.csv(clutch.all.max.diff2,
          "02_output-files/fert type max diff site.csv", row.names = F)

write.csv(clutch.all.update2,
          "02_output-files/fert type by clutch unsamp_diff update.csv", 
          row.names = F)


## plot for min diff----------#
# 
# # adjust factor order
# clutch.all.min.diff2$fert_type <- factor(clutch.all.min.diff2$fert_type, 
#                                     levels=c("ep_diff",
#                                              "ep_same","wp"))
# 
# min.diff.bar <- ggplot(clutch.all.min.diff2, aes(x=year, y=fert, fill=fert_type)) + 
#   geom_col() +
#   scale_fill_manual(values=c("#FDE725FF","#31688EFF","#440154FF"))+
#   theme_light() +
#   ylab("Number of fertilizations") + xlab(NULL) +
#   ggtitle("With minimum possible 'different structure'")
# 
# min.diff.bar
# 
# # ggsave("02_output-files/fert type by sire bar plot all years min diff site.png")
# 
# 
# ## plot for max diff----------#
# 
# # adjust factor order
# clutch.all.max.diff2$fert_type <- factor(clutch.all.max.diff2$fert_type, 
#                                         levels=c("ep_diff",
#                                                  "ep_same","wp"))
# 
# max.diff.bar <- ggplot(clutch.all.max.diff2, aes(x=year, y=fert, fill=fert_type)) + 
#   geom_col() +
#   scale_fill_manual(values=c("#FDE725FF","#31688EFF","#440154FF"))+
#   theme_light() +
#   ylab("Number of fertilizations") + xlab("Year") +
#   ggtitle("With maximum possible 'different structure'")
# 
# max.diff.bar
# 
# # ggsave("02_output-files/fert type by sire bar plot all years max diff site.png")
# 
# 
## plot for ferts with unsampled----------#



# adjust factor order
clutch.all.update2$fert_type <- factor(clutch.all.update2$fert_type,
                                    levels=c("ep_unsampled","ep_diff",
                                             "ep_same","wp"))

raw.bar <- ggplot(clutch.all.update2, aes(x=year, y=fert, fill=fert_type)) +
  geom_col() +
  scale_fill_manual(values=c("#35B779FF","#FDE725FF","#31688EFF","#440154FF"))+
  theme_light() +
  ylab("Number of fertilizations") + xlab("Year")

raw.bar

ggsave("02_output-files/fert types with unclassified.png", w=4, h=4)


# ## combine plots-----------------------------#
# ggarrange(plotlist = list(raw.bar, min.diff.bar, max.diff.bar),
#           ncol=1, nrow=3,
#           common.legend = T,legend = "bottom",
#           labels=c("A","B","C"))
# 
# ggsave("02_output-files/fert type bar plots multi.png", h=6, w=3, scale=1.5)


# Subset each scenario and year combination for plotting
# keep only EP ferts


# update factor orders
clutch.all.update2$fert_type <- factor(clutch.all.update2$fert_type,
                                       levels=c("ep_diff","ep_unsampled",
                                                "ep_same","wp"))

clutch.all.min.diff2$fert_type <- factor(clutch.all.min.diff2$fert_type,
                                    levels=c("ep_diff",
                                             "ep_same","wp"))

clutch.all.max.diff2$fert_type <- factor(clutch.all.max.diff2$fert_type,
                                        levels=c("ep_diff",
                                                 "ep_same","wp"))

# add scenario label
clutch.all.update2$scenario <- "Raw"
clutch.all.min.diff2$scenario <- "Min"
clutch.all.max.diff2$scenario <- "Max"

# add genetic fam for troubleshooting
clutch.all.update2$genetic_fam <- paste(clutch.all.update2$Band_dad, 
                                        clutch.all.update2$Band_mom, sep="_")

clutch.all.min.diff2$genetic_fam <- paste(clutch.all.min.diff2$Band_dad, 
                                        clutch.all.min.diff2$Band_mom, sep="_")

clutch.all.max.diff2$genetic_fam <- paste(clutch.all.max.diff2$Band_dad, 
                                          clutch.all.max.diff2$Band_mom, 
                                          sep="_")


# min
min21 <- subset(clutch.all.min.diff2, clutch.all.min.diff2$year == 2021 &
                  clutch.all.min.diff2$fert_type != "wp")

min22 <- subset(clutch.all.min.diff2, clutch.all.min.diff2$year == 2022 &
                  clutch.all.min.diff2$fert_type != "wp")

min23 <- subset(clutch.all.min.diff2, clutch.all.min.diff2$year == 2023 &
                  clutch.all.min.diff2$fert_type != "wp")

# max
max21 <- subset(clutch.all.max.diff2, clutch.all.max.diff2$year == 2021 &
                  clutch.all.max.diff2$fert_type != "wp")

max22 <- subset(clutch.all.max.diff2, clutch.all.max.diff2$year == 2022 &
                  clutch.all.max.diff2$fert_type != "wp")

max23 <- subset(clutch.all.max.diff2, clutch.all.max.diff2$year == 2023 &
                  clutch.all.max.diff2$fert_type != "wp")

# raw
raw21 <- subset(clutch.all.update2, clutch.all.update2$year == 2021 &
                  clutch.all.update2$fert_type != "wp")

raw22 <- subset(clutch.all.update2, clutch.all.update2$year == 2022 &
                  clutch.all.update2$fert_type != "wp")

raw23 <- subset(clutch.all.update2, clutch.all.update2$year == 2023 &
                  clutch.all.update2$fert_type != "wp")



# Plot for 2021 --------------------------------# 

plot21.df <- rbind(raw21, min21, max21)
plot21.df$scenario <- factor(plot21.df$scenario,
                                         levels=c("Raw", "Min","Max"))

fert21 <- ggplot(plot21.df, aes(x=scenario, y=fert, fill=fert_type)) + 
    geom_col() +
    scale_fill_manual(values=c("#FDE725FF","#35B779FF","#31688EFF"))+
    theme_light() +
    ylab("Number of fertilizations") +
    ggtitle("2021") + xlab(NULL)


# Plot for 2022 --------------------------------# 

plot22.df <- rbind(raw22, min22, max22)
plot22.df$scenario <- factor(plot22.df$scenario,
                             levels=c("Raw", "Min","Max"))

fert22 <- ggplot(plot22.df, aes(x=scenario, y=fert, fill=fert_type)) + 
  geom_col() +
  scale_fill_manual(values=c("#FDE725FF","#35B779FF","#31688EFF"))+
  theme_light() +
  ylab("Number of fertilizations") +
  ggtitle("2022") + xlab(NULL)


# Plot for 2023 --------------------------------# 

plot23.df <- rbind(raw23, min23, max23)
plot23.df$scenario <- factor(plot23.df$scenario,
                             levels=c("Raw", "Min","Max"))

fert23 <- ggplot(plot23.df, aes(x=scenario, y=fert, fill=fert_type)) + 
  geom_col() +
  scale_fill_manual(values=c("#FDE725FF","#35B779FF","#31688EFF"))+
  theme_light() +
  ylab("Number of fertilizations") +
  ggtitle("2023") + xlab(NULL)


# plot all years together------------------------# 

plot.years.df <- rbind(plot21.df, plot22.df, plot23.df)

fert.all.scenario <- ggplot(plot.years.df, aes(x=scenario, y=fert, fill=fert_type)) + 
  geom_col() +
  scale_fill_manual(values=c("#FDE725FF","#35B779FF","#31688EFF"))+
  theme_light() + ylim(0,120) +
  ylab("Number of fertilizations") +
  ggtitle("Reclassified fertilizations for each year") + facet_grid(.~year) +
  xlab("Scenario")
  
ggsave("02_output-files/fert type bar plots EP reclassify.png", h=3, w=6)

#-------------------------------------------------------------------------------

## Number of fert types in each year -------------------#
fert.type.table <- clutch.all.update2 %>% group_by(fert_type, year) %>%
  summarise(num.fert = sum(fert)) %>%
  group_by(year) %>%
  mutate(tot.fert = sum(num.fert),
         prop.fert = (num.fert/tot.fert)*100) 




#-------------------------------------------------------------------------------
# Chi-squared tests for ep fert type
#-------------------------------------------------------------------------------

## for "raw" fert types -------------------------------------------------------#
clutch.ep.update2 <- subset(clutch.all.update2, 
                            clutch.all.update2$fert_type != "wp")

# all years
ep.sum.year <- clutch.ep.update2 %>% group_by(fert_type) %>%
  summarise(n.fert = sum(fert))
# fert_type    n.fert
# <fct>         <int>
# 1 ep_unsampled     51
# 2 ep_diff          71
# 3 ep_same         101



clutch.ep.sum <- clutch.ep.update2 %>% group_by(year, fert_type) %>%
  summarise(n.fert = sum(fert))

clutch.ep.sum
# year fert_type    n.fert
# <dbl> <fct>         <int>
# 1  2021 ep_unsampled     23
# 2  2021 ep_diff          21
# 3  2021 ep_same          24
# 4  2022 ep_unsampled     12
# 5  2022 ep_diff          35
# 6  2022 ep_same          59
# 7  2023 ep_unsampled     16
# 8  2023 ep_diff          15
# 9  2023 ep_same          18

# 2021
ep.obs.21 <- c(rep("ep_same", 24), rep("ep_diff", 21))
ep.obs.21.table <- table(ep.obs.21)
chisq.test(ep.obs.21.table)
# Chi-squared test for given probabilities
# 
# data:  ep.obs.21.table
# X-squared = 0.2, df = 1, p-value = 0.6547


# 2022
ep.obs.22 <- c(rep("ep_same", 59), rep("ep_diff", 35))
ep.obs.22.table <- table(ep.obs.22)
chisq.test(ep.obs.22.table)
# Chi-squared test for given probabilities
# 
# data:  ep.obs.22.table
# X-squared = 6.1277, df = 1, p-value = 0.01331


# 2023
ep.obs.23 <- c(rep("ep_same", 18), rep("ep_diff", 15))
ep.obs.23.table <- table(ep.obs.23)
chisq.test(ep.obs.23.table)
# Chi-squared test for given probabilities
# 
# data:  ep.obs.23.table
# X-squared = 0.27273, df = 1, p-value = 0.6015


## for min diff site ----------------------------------------------------------#

min.diff.ep <- subset(clutch.all.min.diff2, 
                               clutch.all.min.diff2$fert_type != "wp")

min.ep.sum <- min.diff.ep %>% group_by(year, fert_type) %>%
  summarise(n.fert = sum(fert))

min.ep.sum
# year fert_type n.fert
# <dbl> <fct>      <int>
# 1  2021 ep_diff       26
# 2  2021 ep_same       42
# 3  2022 ep_diff       38
# 4  2022 ep_same       68
# 5  2023 ep_diff       19
# 6  2023 ep_same       31

# 2021
ep.min.21 <- c(rep("ep_same", 42), rep("ep_diff", 26))
ep.min.21.table <- table(ep.min.21)
chisq.test(ep.min.21.table)
# Chi-squared test for given probabilities
# 
# data:  ep.min.21.table
# X-squared = 3.7647, df = 1, p-value = 0.05235


# 2022
ep.min.22 <- c(rep("ep_same", 68), rep("ep_diff", 38))
ep.min.22.table <- table(ep.min.22)
chisq.test(ep.min.22.table)
# Chi-squared test for given probabilities
# 
# data:  ep.min.22.table
# X-squared = 8.4906, df = 1, p-value = 0.00357


# 2023
ep.min.23 <- c(rep("ep_same", 31), rep("ep_diff", 19))
ep.min.23.table <- table(ep.min.23)
chisq.test(ep.min.23.table)
# Chi-squared test for given probabilities
# 
# data:  ep.min.23.table
# X-squared = 2.88, df = 1, p-value = 0.08969


## for max diff site ----------------------------------------------------------#

max.diff.ep <- subset(clutch.all.max.diff2, 
                      clutch.all.max.diff2$fert_type != "wp")

max.ep.sum <- max.diff.ep %>% group_by(year, fert_type) %>%
  summarise(n.fert = sum(fert))

max.ep.sum
# year fert_type n.fert
# <dbl> <fct>      <int>
# 1  2021 ep_diff       44
# 2  2021 ep_same       24
# 3  2022 ep_diff       47
# 4  2022 ep_same       59
# 5  2023 ep_diff       31
# 6  2023 ep_same       18

# 2021
ep.max.21 <- c(rep("ep_same", 24), rep("ep_diff", 44))
ep.max.21.table <- table(ep.max.21)
chisq.test(ep.max.21.table)
# Chi-squared test for given probabilities
# 
# data:  ep.max.21.table
# X-squared = 5.8824, df = 1, p-value = 0.01529


# 2022
ep.max.22 <- c(rep("ep_same", 59), rep("ep_diff", 47))
ep.max.22.table <- table(ep.max.22)
chisq.test(ep.max.22.table)
# Chi-squared test for given probabilities
# 
# data:  ep.max.22.table
# X-squared = 1.3585, df = 1, p-value = 0.2438


# 2023
ep.max.23 <- c(rep("ep_same", 18), rep("ep_diff", 31))
ep.max.23.table <- table(ep.max.23)
chisq.test(ep.max.23.table)
# Chi-squared test for given probabilities
# 
# data:  ep.max.23.table
# X-squared = 3.449, df = 1, p-value = 0.06329



