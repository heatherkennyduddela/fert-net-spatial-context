
#-------------------------------------------------------------------------------
# Calculate number of total unknown sires for 2022
# Heather Kenny-Duddela
# May 23, 2026
#-------------------------------------------------------------------------------


# libraries
library(tidyverse)
library(ggplot2)


### load data tables for 2022

# full kinship table with metadata
# from 01_paternity 2022 assignment with lcMLkin output
load("input-files/kin2_22.Rdata")

# table of offspring with assigned parents
# from 02_summarise total RS and fertilization types
assigned <- read.csv("input-files/kin_2022_parent_offspring_assigned_wide.csv")


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

kin.hs2 <- inner_join(storage.hs2, kin2.22, by=c("Ind1","Ind2"))

# add column for sib relationship
kin.hs2$sib <- "unrelated"
kin.hs2$sib[which(kin.hs2$mom_same==T & kin.hs2$dad_same==T)] <- "full"
kin.hs2$sib[which(kin.hs2$mom_same==T & kin.hs2$dad_same==F |
                    kin.hs2$mom_same==F & kin.hs2$dad_same==T )] <- "half"


# plot sib and half sib relatedness --------------------------------------------
plot(kin.hs2$k0_hat, kin.hs2$pi_HAT)

# histograms of pi_hat for full and half
ggplot(kin.hs2, aes(x=pi_HAT, fill=sib)) + geom_histogram() +
  facet_grid(sib~.)

# scatterplot with colored points
ggplot(kin.hs2, aes(x=k0_hat, y=pi_HAT, color=sib)) +
  geom_point(shape=2)

# zoom in on overlapping area for full and half sibs
ggplot(kin.hs2, aes(x=k0_hat, y=pi_HAT, color=sib)) +
  geom_point(shape=2) + ylim(0.25, 0.35) +
  xlim(0.4,0.7)

# identify area of relatedness overlap -----------------------------------------

# based on visual assessment of plot
kin.hs.close <- subset(kin.hs2, kin.hs2$pi_HAT<0.325 & kin.hs2$pi_HAT>0.25 &
                          kin.hs2$k0_hat>0.5 & kin.hs2$k0_hat<0.65)


#-------------------------------------------------------------------------------
# Get kinship values for kids with unknown dads
#-------------------------------------------------------------------------------

# get list of kids with unknown dads
dad.unk <- subset(assigned, is.na(assigned$pi_HAT_dad))

# pull out kinship values for unrelated to unrelated
unk2unk <- kin2.22[which(kin2.22$Ind1 %in% dad.unk$Ind2 &
                           kin2.22$Ind2 %in% dad.unk$Ind2), ]

# plot unknown sibs over known relationships
ggplot(kin.hs2, aes(x=k0_hat, y=pi_HAT, color=sib)) +
  geom_point(shape=2) +
  geom_point(data=unk2unk, aes(x=k0_hat, y=pi_HAT), color="black")


# use pi_hat values for known half-sibs and unrelated to classify unknowns

# wide cutoff uses max and min of known half-sibs
max(subset(kin.hs2$pi_HAT, kin.hs2$sib=="half")) # 0.311
min(subset(kin.hs2$pi_HAT, kin.hs2$sib=="half")) # 0.077

# narrow cutoff uses min of full and max of unrelated
min(subset(kin.hs2$pi_HAT, kin.hs2$sib=="full")) # 0.244
max(subset(kin.hs2$pi_HAT, kin.hs2$sib=="unrelated")) # 0.156

# make overlay plot again with cutoff lines
ggplot(kin.hs2, aes(x=k0_hat, y=pi_HAT, color=sib)) +
  geom_point(shape=2) +
  geom_point(data=unk2unk, aes(x=k0_hat, y=pi_HAT), color="black")+
  geom_hline(yintercept=0.311, linetype="dashed") +
  geom_hline(yintercept=0.077, linetype="dashed") +
  geom_hline(yintercept=0.244, linetype="solid") +
  geom_hline(yintercept=0.156, linetype="solid")
  

# save plot
ggsave("generated-files/unknown kinship overlayed with known sibs cutoff lines 2022.png", h=5, w=7)


# add classifications to unknowns for wide and narrow cutoffs

# wide
unk2unk$sib_wide <- NA
unk2unk$sib_wide[unk2unk$pi_HAT>0.311] <- "full"
unk2unk$sib_wide[unk2unk$pi_HAT<=0.311 &
                   unk2unk$pi_HAT>=0.077] <- "half"
unk2unk$sib_wide[unk2unk$pi_HAT<=0.077] <- "unrelated"

# narrow
unk2unk$sib_narrow <- NA
unk2unk$sib_narrow[unk2unk$pi_HAT>0.244] <- "full"
unk2unk$sib_narrow[unk2unk$pi_HAT<=0.244 &
                   unk2unk$pi_HAT>=0.156] <- "half"
unk2unk$sib_narrow[unk2unk$pi_HAT<=0.156] <- "unrelated"


#-------------------------------------------------------------------------------
# Calculate total number of unknown sires
#-------------------------------------------------------------------------------

# first categorize matching or not matching familyIDs
unk2unk$famID_check <- ifelse(unk2unk$FamilyID_ind1 == unk2unk$FamilyID_ind2, T, F)

# check that all full sibs have same familyID
sum(unk2unk$sib_wide=="full" & unk2unk$famID_check==T) # 31
sum(unk2unk$sib_wide=="full") # 31

# check that unrelated have different familyIDs
sum(unk2unk$sib_wide=="unrelated" & unk2unk$famID_check==F) # 663
sum(unk2unk$sib_wide=="unrelated") # 663

# subset only full and half sibs
unk_fh <- subset(unk2unk, unk2unk$sib_wide != "unrelated")

# pull out cases of half sibs where family ID is different
# this means the kids have the same dad but different moms
unk_h_famID_diff <- subset(unk_fh, unk_fh$famID_check==F & unk_fh$sib_wide=="half")

# A single unknown dad sired offspring from familyID CHR-74 (clutch CHR-6_2) and
# familyID CHR-93 (clutch CHR-93_2)

# check these two families compared to each other
CHR74_CHR93_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="CHR-074" & 
                              unk2unk$FamilyID_ind2=="CHR-093" |
                              unk2unk$FamilyID_ind1=="CHR-093" &
                              unk2unk$FamilyID_ind2=="CHR-074" |
                              unk2unk$FamilyID_ind1=="CHR-074" &
                              unk2unk$FamilyID_ind2=="CHR-074" |
                              unk2unk$FamilyID_ind1=="CHR-093" &
                              unk2unk$FamilyID_ind2=="CHR-093")
# This shows that the female from family CHR-074 mated with one unk dad for her
# first clutch (74_1) and a different unk dad for her second clutch (6_2).
# The same unk male from clutch 6_2 also sired one offspring with the female 
# from family CHR-093 in clutch 93_2. 


# pull out cases of half sibs where family ID is same
# this means the kids have the same mom but different dads
unk_h_famID_same <- subset(unk_fh, unk_fh$famID_check==T & unk_fh$sib_wide=="half")

## From wide cutoffs

# Two different unk dads sired offspring in family CHR-74, one in clutch 
# CHR-74_1 and one in clutch CHR_6_2

# Two different unk dads sired offspring in family CHR-112, one in clutch 
# CHR_112_1 and one in clutch CHR_112.2_2. This makes sense for the wide 
# cutoff sib relationships, but not with the narrow cutoff relationships.
# The narrow and wide cutoffs agree about full sib relationships for individuals
# from the same clutch. Across clutches with the same mom, the wide cutoffs
# show only half-sib relationships which makes sense. Across clutches for the
# narrow cutoff some individuals are being called unrelated, or full sibs but
# these relationships are not consistent between full siblings. For example, 
# 83300 is from clutch 1 and is full sibs with 83413 who is also from clutch 1.
# 83300 is called unrelated to 83443 from clutch 2, but 83413 is called a half
# sib with 83443. 

# Check family CHR_112
CHR112_check <- subset(unk2unk, unk2unk$FamilyID_ind1=="CHR-112" &
                         unk2unk$FamilyID_ind2=="CHR-112")


# determine number of unique family IDs
unique_fam <- unique(c(unk2unk$FamilyID_ind1, unk2unk$FamilyID_ind2)) # 16
unique_fam

# There is one unk sire per family, except for the cases discussed above
# CHR-74 had two sires, but one of those was also the sire for CHR-93
# CHR-112 also had two sires

# 16 + 1 -1 +1 = 17

#-------------------------------------------------------------------------------
# Calculate total RS for each unsampled sire
#-------------------------------------------------------------------------------

# summarize unknown kids by clutch and famID
unk.sum <- dad.unk %>% group_by(clutch_id_ind2, FamilyID_mom) %>%
  summarise(n = n())

# make column for unsampled male ID
unk.sum$uns_male <- paste(unk.sum$FamilyID_mom,"uns1", sep="_")

# change male ID for CHR-112 which had two differen males, one for each clutch
unk.sum$uns_male[which(unk.sum$clutch_id_ind2=="CHR_112.2_2")] <- "CHR-112_uns2"

# change male ID for CHR_6_2 to have an additional sire
unk.sum$uns_male[which(unk.sum$clutch_id_ind2=="CHR_6_2")] <- "CHR-074_uns2"

# change male ID for CHR-93 which had same sire as CHR_6_2
unk.sum$uns_male[which(unk.sum$clutch_id_ind2=="CHR_93_2")] <- "CHR-074_uns2"

# sum offspring by unsampled male ID
uns_male_RS <- unk.sum %>% group_by(uns_male) %>%
  summarise(n_chick = sum(n))

# save table
write.csv(uns_male_RS, "generated-files/unsampled male RS 2022.csv", row.names=F)

#-------------------------------------------------------------------------------
# update assigned table with uns_male IDs
#-------------------------------------------------------------------------------

# add uns_male IDs to table where each row is a kid
dad.unk2 <- left_join(dad.unk, unk.sum[,c(1,4)], by=c("clutch_id_ind2"))

# update Band_dad 
dad.unk2$Band_dad <- dad.unk2$uns_male

# update genetic family
dad.unk2$genetic_fam <- paste(dad.unk2$Band_dad, dad.unk2$Band_mom, sep="_")

## bind updated unk with rest of assigned
assigned.notUNS <- subset(assigned, !is.na(assigned$Band_dad))

assigned.with.uns <- rbind(assigned.notUNS, dad.unk2[,-30])

# save table
write.csv(assigned.with.uns, 
          "generated-files/kin_2022_parent_offspring_assigned_withUNS.csv",
          row.names=F)



