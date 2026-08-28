
#-------------------------------------------------------------------------------
# Calculate total number of unsampled sires for 2023
# Heather Kenny-Duddela
# May 21, 2026
#-------------------------------------------------------------------------------

# libraries
library(tidyverse)
library(ggplot2)

# set working directory
setwd("~/CU Boulder/BARS fieldwork/2022 Field and Lab/Ch 3 analyses/02_paternity-assignment/2023")


### load data tables for 2023

# full kinship table with metadata
load("01_output-files/kin2_23.Rdata")

# table of offspring with assigned parents
assigned <- read.csv("02_output-files/kin_2023_po_wide.csv")



#-------------------------------------------------------------------------------
# data table modifications
#-------------------------------------------------------------------------------

# Noticed an error in labeling of familyID for kids from Urban Farm Girlz_8_1
# This was a tricky nest because the female (2850-57629) nested here for her
# first brood, but then moved to nest 5 for her second brood. A different female
# took over nest 8.2 later in the season (2850-57764). Update familyID for the
# kids from clutch 8_1 from UrbanFarm-08 to UrbanFarm-05. 

# This has already been updated in the assigned table. 

kin2.23$FamilyID_ind2[which(kin2.23$clutchID_ind2=="Urban Farm Girlz_8_1")] <-
  "UrbanFarm-05"

kin2.23$FamilyID_ind1[which(kin2.23$clutchID_ind1=="Urban Farm Girlz_8_1")] <-
  "UrbanFarm-05"


# One of the Boyer females also had an incorrect family ID which was already 
# corrected in the assigned table but not in the kin2.23 table. Female 2850-57688 
# had offspring in nests 13 (eggs) and 15 (chicks), but not in nest 5 as was
# assumed from IDobs.Change her family ID from Boyer-05 to Boyer-13 

kin2.23$FamilyID_ind2[which(kin2.23$band_ind2=="2850-57688")] <-
  "Boyer-13"

kin2.23$FamilyID_ind1[which(kin2.23$band_ind2=="2850-57688")] <-
  "Boyer-13"

# Another Boyer female 2850-57685 had offspring in nests 5 (eggs), 13, and 16
# Kids from clutch Boyer_16_2 and the female need to have their family ID changed from 
# Boyer-13 to Boyer-05 in kin2.23. Already corrected in assigned table. 

kin2.23$FamilyID_ind2[which(kin2.23$band_ind2=="2850-57685")] <-
  "Boyer-05"

kin2.23$FamilyID_ind1[which(kin2.23$band_ind2=="2850-57685")] <-
  "Boyer-05"

kin2.23$FamilyID_ind2[which(kin2.23$clutchID_ind2=="Boyer_16_2")] <-
  "Boyer-05"

kin2.23$FamilyID_ind1[which(kin2.23$clutchID_ind1=="Boyer_16_2")] <-
  "Boyer-05"

#-------------------------------------------------------------------------------

# remove kids with unknown mom or dad
assigned.noNA <- subset(assigned, !is.na(assigned$band_dad) &
                          !is.na(assigned$band_mom))

storage.hs2 <- as.data.frame(matrix(nrow=1,ncol=6,NA))
colnames(storage.hs2) <- c("Ind1","Ind2","mom_same","dad_same", "mom_band", "dad_band")

# loop through kids
for (i in 1:length(assigned.noNA$Ind2)) {
  # loop through moms and dads
  for (j in 1:length(assigned.noNA$band_mom)) {
    # case of same dad, different mom
    if(assigned.noNA$band_dad[i] == assigned.noNA$band_dad[j] &
       assigned.noNA$band_mom[i] != assigned.noNA$band_mom[j]) {
      storage <- c(assigned.noNA$Ind2[i], 
                   assigned.noNA$Ind2[j],
                   F,T, NA, assigned.noNA$band_dad[i])
      storage.hs2 <- rbind(storage.hs2, storage)
    }
    # case of same mom, different dad
    if(assigned.noNA$band_dad[i] != assigned.noNA$band_dad[j] &
       assigned.noNA$band_mom[i] == assigned.noNA$band_mom[j]) {
      storage <- c(assigned.noNA$Ind2[i], 
                   assigned.noNA$Ind2[j],
                   T,F, assigned.noNA$band_mom[i], NA)
      storage.hs2 <- rbind(storage.hs2, storage)
    }
    # case of same mom, same dad (full sibs)
    if(assigned.noNA$band_dad[i] == assigned.noNA$band_dad[j] &
       assigned.noNA$band_mom[i] == assigned.noNA$band_mom[j]) {
      storage <- c(assigned.noNA$Ind2[i], 
                   assigned.noNA$Ind2[j],
                   T,T, assigned.noNA$band_mom[i], assigned.noNA$band_dad[i])
      storage.hs2 <- rbind(storage.hs2, storage)
    }
    # case of diff mom, diff dad (unrelated)
    if(assigned.noNA$band_dad[i] != assigned.noNA$band_dad[j] &
       assigned.noNA$band_mom[i] != assigned.noNA$band_mom[j]) {
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

kin.hs2 <- inner_join(storage.hs2, kin2.23, by=c("Ind1","Ind2"))

# add column for sib relationship
kin.hs2$sib <- "unrelated"
kin.hs2$sib[which(kin.hs2$mom_same==T & kin.hs2$dad_same==T)] <- "full"
kin.hs2$sib[which(kin.hs2$mom_same==T & kin.hs2$dad_same==F |
                    kin.hs2$mom_same==F & kin.hs2$dad_same==T )] <- "half"

# scatterplot with colored points
ggplot(kin.hs2, aes(x=k0_hat, y=pi_HAT, color=sib)) +
  geom_point(shape=2)


#-------------------------------------------------------------------------------
# Get kinship values for kids with unknown dads
#-------------------------------------------------------------------------------

# get list of kids with unknown dads
dad.unk <- subset(assigned, is.na(assigned$pi_HAT_dad))

# pull out kinship values for unrelated to unrelated
unk2unk <- kin2.23[which(kin2.23$Ind1 %in% dad.unk$Ind2 &
                           kin2.23$Ind2 %in% dad.unk$Ind2), ]

# plot unknown sibs over known relationships
ggplot(kin.hs2, aes(x=k0_hat, y=pi_HAT, color=sib)) +
  geom_point(shape=2) +
  geom_point(data=unk2unk, aes(x=k0_hat, y=pi_HAT), color="black")


# use pi_hat values for known half-sibs and unrelated to classify unknowns

# wide cutoff uses max and min of known half-sibs
max(subset(kin.hs2$pi_HAT, kin.hs2$sib=="half")) # 0.321
min(subset(kin.hs2$pi_HAT, kin.hs2$sib=="half")) # 0.044

# narrow cutoff uses min of full and max of unrelated
min(subset(kin.hs2$pi_HAT, kin.hs2$sib=="full")) # 0.253
max(subset(kin.hs2$pi_HAT, kin.hs2$sib=="unrelated")) # 0.104

# make overlay plot again with cutoff lines
sib_categories_2023 <- ggplot(kin.hs2, aes(x=k0_hat, y=pi_HAT, color=sib)) +
  geom_point(shape=2) +
  geom_point(data=unk2unk, aes(x=k0_hat, y=pi_HAT), color="black")+
  geom_hline(yintercept=0.321, linetype="dashed") +
  geom_hline(yintercept=0.044, linetype="dashed") +
  geom_hline(yintercept=0.253, linetype="solid") +
  geom_hline(yintercept=0.104, linetype="solid")

# save plot
ggsave("04_output-files/unknown kinship overlayed with known sibs cutoff lines 2023.png", h=5, w=7)

save(sib_categories_2023, file="04_output-files/sib_categories_2023_plot.Rdata")


# add classifications to unknowns for wide and narrow cutoffs

# wide
unk2unk$sib_wide <- NA
unk2unk$sib_wide[unk2unk$pi_HAT>0.321] <- "full"
unk2unk$sib_wide[unk2unk$pi_HAT<=0.321 &
                   unk2unk$pi_HAT>=0.044] <- "half"
unk2unk$sib_wide[unk2unk$pi_HAT<=0.044] <- "unrelated"

# narrow
unk2unk$sib_narrow <- NA
unk2unk$sib_narrow[unk2unk$pi_HAT>0.253] <- "full"
unk2unk$sib_narrow[unk2unk$pi_HAT<=0.253 &
                     unk2unk$pi_HAT>=0.104] <- "half"
unk2unk$sib_narrow[unk2unk$pi_HAT<=0.104] <- "unrelated"

#-------------------------------------------------------------------------------
# Calculate total number of unknown sires
#-------------------------------------------------------------------------------

# first categorize matching or not matching familyIDs
unk2unk$famID_check <- ifelse(unk2unk$FamilyID_ind1 == unk2unk$FamilyID_ind2, T, F)

# check that all full sibs have same familyID
sum(unk2unk$sib_wide=="full" & unk2unk$famID_check==T) # 25
sum(unk2unk$sib_wide=="full") # 25

# check that unrelated have different familyIDs
sum(unk2unk$sib_wide=="unrelated" & unk2unk$famID_check==F) # 518
sum(unk2unk$sib_wide=="unrelated") # 518

# subset only full and half sibs
unk_fh <- subset(unk2unk, unk2unk$sib_wide != "unrelated")

# Same dad but different moms---------------------------------------------------

# pull out cases of half sibs where family ID is different
unk_h_famID_diff <- subset(unk_fh, unk_fh$famID_check==F & unk_fh$sib_wide=="half")

## check these cases of same dad and different mom
# Boyer 2.3_2 and Boyer 16_2
# Cooks 17.2 and Cooks 29
# Boyer 3.2 and UrbanFarm 8.2
# Cooks 2 and Cooks 18

Boyer2.3_Boyer16_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Boyer-02" & 
                              unk2unk$FamilyID_ind2=="Boyer-05" |
                              unk2unk$FamilyID_ind1=="Boyer-05" &
                              unk2unk$FamilyID_ind2=="Boyer-02" |
                              unk2unk$FamilyID_ind1=="Boyer-02" &
                              unk2unk$FamilyID_ind2=="Boyer-02" |
                              unk2unk$FamilyID_ind1=="Boyer-05" &
                              unk2unk$FamilyID_ind2=="Boyer-05")
# kids in clutch 2.3 (family 2) are all full sibs, and they are all half sibs with the one
# kid from clutch 16 (family 05). The same unsampled male sired chicks in both nests. 

Cooks17.2_Cooks29_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Cooks-15" & 
                                   unk2unk$FamilyID_ind2=="Cooks-29" |
                                   unk2unk$FamilyID_ind1=="Cooks-29" &
                                   unk2unk$FamilyID_ind2=="Cooks-15" |
                                   unk2unk$FamilyID_ind1=="Cooks-15" &
                                   unk2unk$FamilyID_ind2=="Cooks-15" |
                                   unk2unk$FamilyID_ind1=="Cooks-29" &
                                   unk2unk$FamilyID_ind2=="Cooks-29")
# kids from 17.2 are all full sibs. They are all half sibs with one kid from
# 17.1 because they have the same mom (2850-57948) but different dads, UNS-1 for
# 17.1 and UNS-2 for 17.2. Kids from clutch 17.2 are also half sibs with one
# egg from 29. This egg has mom 2850-57853 and must share the dad UNS-2 from 17.2. 

Boyer3.2_UrbanFarm8.2_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Boyer-03" & 
                                    unk2unk$FamilyID_ind2=="UrbanFarm-08" |
                                    unk2unk$FamilyID_ind1=="UrbanFarm-08" &
                                    unk2unk$FamilyID_ind2=="Boyer-03" |
                                    unk2unk$FamilyID_ind1=="Boyer-03" &
                                    unk2unk$FamilyID_ind2=="Boyer-03" |
                                    unk2unk$FamilyID_ind1=="UrbanFarm-08" &
                                    unk2unk$FamilyID_ind2=="UrbanFarm-08")
# the two kids from Boyer 3.2 are full sibs. They are both half sibs with the 
# one kid from Urban Farm 8.2. The narrow sibling categories say that one of the
# Boyer kids is a full sib with the UF kid, while the other Boyer kid is a half
# sib with the UF kid. This doesn't make sense because both wide and narrow say
# that the Boyer kids are full siblings. Go with the inference from the wide. 
# A single male sired kids in both nests!

Cooks2_Cooks18_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Cooks-02" & 
                                    unk2unk$FamilyID_ind2=="Cooks-18" |
                                    unk2unk$FamilyID_ind1=="Cooks-18" &
                                    unk2unk$FamilyID_ind2=="Cooks-02" |
                                    unk2unk$FamilyID_ind1=="Cooks-02" &
                                    unk2unk$FamilyID_ind2=="Cooks-02" |
                                    unk2unk$FamilyID_ind1=="Cooks-18" &
                                    unk2unk$FamilyID_ind2=="Cooks-18")
# the four Cooks 18 kids are all full sibs, and they are each half sibs with 
# the one kid from Cooks 2. The same male sired everyone. 


# Same mom but different dads --------------------------------------------------

# pull out cases of half sibs where family ID is same
unk_h_famID_same <- subset(unk_fh, unk_fh$famID_check==T & unk_fh$sib_wide=="half")

# Three kids from Cooks 17.2 with one kid from 17.1
# Three kids from Blue Cloud 18.2 with one other kid from Blue Cloud 18.2

Cooks15_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="Cooks-15" & 
                                 unk2unk$FamilyID_ind2=="Cooks-15" )
# One male for each clutch, so two total

BlueCloud18_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="BlueCloud-18" & 
                              unk2unk$FamilyID_ind2=="BlueCloud-18" )
# three of the kids are all full sibs, and are half sibs with the fourth kid. 
# Two different unsampled males for this clutch


# Determine number of family IDs and total unique unsampled males---------------

unique_fam <- unique(c(unk2unk$FamilyID_ind1, unk2unk$FamilyID_ind2))
length(unique_fam) # 16

# There is one unk sire per family, except for the cases discussed above
# Boyer-2 and Boyer-05 had the same sire (-1)
# Cooks-15 had two sires (+1), but one sire was shared with Cooks-29 (-1)
# Boyer-3 and UrbanFarm-8 shared the same sire (-1)
# Cooks-18 and  Cooks-2 shared the same sire (-1)
# BlueCloud-18 had two sires (+1)

# 16 -1 +1 -1 -1 -1 +1 = 14 total sires

#-------------------------------------------------------------------------------
# Calculate total RS for each unsampled sire
#-------------------------------------------------------------------------------

# summarize unknown kids by clutch and famID
unk.sum <- dad.unk %>% group_by(clutchID_ind2, FamilyID_mom, FamilyID_ind2) %>%
  summarise(n = n())

# make column for unsampled male ID
unk.sum$uns_male <- paste(unk.sum$FamilyID_mom,"uns1", sep="_")

# change male ID for Boyer-05 which shared with Boyer-2
unk.sum$uns_male[which(unk.sum$clutchID_ind2=="Boyer_16_2")] <- "Boyer-02_uns1"

# change male ID for Cooks 17.2 to be shared with Cooks-29
unk.sum$uns_male[which(unk.sum$clutchID_ind2=="Cooks_17.2_2")] <- "Cooks-29_uns1"

# change male ID for Boyer 3.2 to share with UrbanFarm-08
unk.sum$uns_male[which(unk.sum$clutchID_ind2=="Boyer_3.2_1")] <- "UrbanFarm-08_uns1"

# change male ID for Cooks-18 to share with Cooks-02
unk.sum$uns_male[which(unk.sum$clutchID_ind2=="Cooks_18_collected")] <- "Cooks-02_uns1"

# Add additional row for BlueCloud-18 because she had two sires
# first sire had 3 offspring, second had 1
unk.sum$n[which(unk.sum$clutchID_ind2=="Blue Cloud_18.2_1")] <- 3
unk.sum2 <- rbind(unk.sum, unk.sum[1, ])
unk.sum2[18,4] <- 1
unk.sum2[18,5] <- "BlueCloud-18_uns2"


# sum offspring by unsampled male ID
uns_male_RS <- unk.sum2 %>% group_by(uns_male) %>%
  summarise(n_chick = sum(n))

# save table
write.csv(uns_male_RS, "04_output-files/unsampled male RS 2023.csv", row.names=F)

#-------------------------------------------------------------------------------
# update assigned table with uns_male IDs
#-------------------------------------------------------------------------------

# add uns_male IDs to table where each row is a kid
# need to manually add the second Blue Cloud sire, otherwise everyone from
# that clutch will get two sires
dad.unk2 <- left_join(dad.unk, unk.sum2[1:17,c(1,5)], by=c("clutchID_ind2"))
# 33919 is the half sibling to the other 3 in this clutch, so give them the
# second sire as the dad
dad.unk2$uns_male[which(dad.unk2$band_ind2=="2881-33919")] <- "BlueCloud-18_uns2"

# update Band_dad 
dad.unk2$band_dad <- dad.unk2$uns_male

# update genetic family
dad.unk2$genetic_fam <- paste(dad.unk2$band_dad, dad.unk2$band_mom, sep="_")

## bind updated unk with rest of assigned
assigned.notUNS <- subset(assigned, !is.na(assigned$band_dad))

assigned.with.uns <- rbind(assigned.notUNS, dad.unk2[,-33])

# save table
write.csv(assigned.with.uns, 
          "04_output-files/kin_2023_parent_offspring_assigned_withUNS.csv",
          row.names=F)




