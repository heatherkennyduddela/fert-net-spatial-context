
################################################################################
# Genetic parentage assignemnt with lcMLkin output for CO 2023
# Heather Kenny-Duddela
# May 14, 2025
################################################################################

# For 2023, 92,430 SNPs analysed for kinship analysis in lcMLkin
# Only autosomes, sex-linked scaffolds removed from dataset
# BARS genome is 1.2GB total size

## Load data

# kinship results from lcMLkin
kin23 <- read.table("input-files/CO-2023batch2.thin10k.GL-recalc.relate", header=T)

# family table from 2023
fam.23 <- read.csv("input-files/Family IDs for paternity 2023.csv")

# load libraries
library(tidyverse) # for dplyr, ggplot2

#-------------------------------------------------------------------------------
# Check example individuals to confirm data structure
#-------------------------------------------------------------------------------

# k0_hat = prob no alleles shared that are IBD (unrelated)
# k1_hat = prob 1 allele shared that is IBD
# k2_hat = prob 2 alleles shared that are IBD
# pi_hat = coefficient of relatedness (r)

# coefficient of relatedness (r) = 2*phi = (k1)/2 + k2

## Plot k0 vs. pi_hat for all birds
plot(kin23$k0_hat, kin23$pi_HAT)


# plot for test individual, female from Blue Cloud
kin23.57928 <- subset(kin23, kin23$Ind1=="CO_57928" |
                        kin23$Ind2=="CO_57928")

plot(kin23.57928$k0_hat, kin23.57928$pi_HAT)

# try calculating r for the test individual
# see that the calculation for r matches pi_hat
kin23.57928$r <- kin23.57928$k1_hat/2 + kin23.57928$k2_hat
cor.test(kin23.57928$r, kin23.57928$pi_HAT)
# Pearson's product-moment correlation
# 
# data:  kin23.57928$r and kin23.57928$pi_HAT
# t = 2076.9, df = 320, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9999538 0.9999702
# sample estimates:
#       cor 
# 0.9999629


## social male for test female
kin23.57659 <- subset(kin23, kin23$Ind1=="CO_57659" |
                        kin23$Ind2=="CO_57659")

plot(kin23.57659$k0_hat, kin23.57659$pi_HAT)

## kid of test female
kin23.33940 <- subset(kin23, kin23$Ind1=="CO_33940" |
                        kin23$Ind2=="CO_33940")

plot(kin23.33940$k0_hat, kin23.33940$pi_HAT)


#-------------------------------------------------------------------------------
# Combine family table info with kinship
#-------------------------------------------------------------------------------

# check for duplicate combos of band and familyID
# no duplicates found
fam.duplicate <- anyDuplicated(fam.23[, c(1,5)])

# keep only columns that are relevant to adults and kids
seq.cols <- select(fam.23, band, site, nest, type, FamilyID, brood, seqID)

# add clutch ID
seq.cols$clutchID <- paste(seq.cols$site, seq.cols$nest, seq.cols$brood,
                           sep="_")

# add suffix to specify ind1 metadata
colnames(seq.cols) <- paste(colnames(seq.cols),"ind1",sep="_")
# change colname for seqID to Ind1 to match with kin23
colnames(seq.cols)[7] <- "Ind1"

# combine kinship data with metadata for ind1
kin1.23 <- left_join(kin23, seq.cols, by="Ind1")
kin23$Ind1[46486]

# check that merge worked properly 
test <- subset(kin1.23, kin1.23$Ind1=="CO_57659" | kin1.23$Ind2=="CO_57659")


## Repeat for adding metadata for ind2
# keep only columns that are relevant to adults and kids
seq.cols2 <- select(fam.23, band, site, nest, type, FamilyID, brood, seqID)

# add clutch ID
seq.cols2$clutchID <- paste(seq.cols2$site, seq.cols2$nest, seq.cols2$brood,
                           sep="_")

# add suffix to specify ind1 metadata
colnames(seq.cols2) <- paste(colnames(seq.cols2),"ind2",sep="_")
# change colname for seqID to Ind1 to match with kin23
colnames(seq.cols2)[7] <- "Ind2"

# combine kinship data with metadata for ind2
kin2.23 <- left_join(kin1.23, seq.cols2, by="Ind2")

# check for NAs - fixed, seqID and band had extra space in fam table
# kin2.ind1na <- filter(kin2.23, is.na(kin2.23$band_ind1))
# unique(kin2.ind1na$Ind1) # "CO_57948"
# kin2.ind2na <- filter(kin2.23, is.na(kin2.23$band_ind2))
# unique(kin2.ind2na$Ind2) # "CO_57948"

# save kinship table with metadata
# this has all pairwise comparisons, including parent-kid, parent-parent, kid-kid
write.csv(kin2.23, "01_output-files/kin_2023_all_with_labels.csv", row.names=F)
save(kin2.23, file="01_output-files/kin2_23.Rdata")

#-------------------------------------------------------------------------------
# Make table of all possible kid-parent combinations
#-------------------------------------------------------------------------------

kin23.po <- subset(kin2.23, 
                     kin2.23$type_ind1=="kid" & kin2.23$type_ind2=="mom" |
                     kin2.23$type_ind1=="kid" & kin2.23$type_ind2=="dad" |
                     kin2.23$type_ind1=="egg" & kin2.23$type_ind2=="mom" |
                     kin2.23$type_ind1=="egg" & kin2.23$type_ind2=="dad" |
                     kin2.23$type_ind2=="kid" & kin2.23$type_ind1=="mom" |
                     kin2.23$type_ind2=="kid" & kin2.23$type_ind1=="dad" |
                     kin2.23$type_ind2=="egg" & kin2.23$type_ind1=="mom" |
                     kin2.23$type_ind2=="egg" & kin2.23$type_ind1=="dad")

## restructure so that ind1 are all parents and ind2 are all kids

## cases where ind1 is kid or egg
kin23.po.ind1kid <- subset(kin23.po, kin23.po$type_ind1=="kid" |
                             kin23.po$type_ind1=="egg")
# keep only ind1 columns
kin23.po.ind1kid.before <- kin23.po.ind1kid[,1:14]
# change ind1 to ind2
colnames(kin23.po.ind1kid.before) <- gsub("ind1", "ind2", 
                                          colnames(kin23.po.ind1kid.before))
colnames(kin23.po.ind1kid.before)[1:2] <- c("Ind2","Ind1")
# reorder cols
kin23.po.ind1kid.after <- relocate(kin23.po.ind1kid.before, Ind2, .after=Ind1)

# keep only ind2 columns
kin23.po.ind2par.before <- kin23.po.ind1kid[,c(1:7, 15:21)]
# change ind2 to ind1
colnames(kin23.po.ind2par.before) <- gsub("ind2", "ind1", 
                                          colnames(kin23.po.ind2par.before))
colnames(kin23.po.ind2par.before)[1:2] <- c("Ind2","Ind1")
# reorder cols
kin23.po.ind2par.after <- relocate(kin23.po.ind2par.before, Ind2, .after=Ind1)

# combine ind1 and ind2
kin23.po.ind1kid.swap <- left_join(kin23.po.ind2par.after, 
                                   kin23.po.ind1kid.after, 
                                   by=c("Ind1", "Ind2", "k0_hat", "k1_hat",
                                        "k2_hat", "pi_HAT", "nbSNP"))

# Add the swapped table to the "good" cases from the original table

# cases where ind1 is parents
kin23.po.ind1par <- subset(kin23.po, kin23.po$type_ind1=="mom" |
                             kin23.po$type_ind1=="dad")
# rowbind the two tables
kin23.po.org <- rbind(kin23.po.ind1par, kin23.po.ind1kid.swap)

# check that organization worked
unique(kin23.po.org$type_ind1)
unique(kin23.po.org$type_ind2)

## check number of kids, moms, and dads
# 216 kids and eggs
length(unique(kin23.po.org$Ind2)) 
# 72 dads
length(unique(subset(kin23.po.org$Ind1, kin23.po.org$type_ind1=="dad")))
# 35 moms
length(unique(subset(kin23.po.org$Ind1, kin23.po.org$type_ind1=="mom")))


# save table of all parent-offspring kinship
write.csv(kin23.po.org, file="01_output-files/kin_2023_parent_offspring_all.csv", 
          row.names=F)

#-------------------------------------------------------------------------------
# Make table of adult only relationships
#-------------------------------------------------------------------------------

kin.adults <- subset(
  kin2.23,
    kin2.23$type_ind1=="mom" & kin2.23$type_ind2=="dad" |
    kin2.23$type_ind1=="dad" & kin2.23$type_ind2=="mom" |
    kin2.23$type_ind1=="mom" & kin2.23$type_ind2=="mom" |
    kin2.23$type_ind1=="dad" & kin2.23$type_ind2=="dad")

# 20 cases of close relatives between adults (pi_hat > 0.1)

# save file
write.csv(kin.adults, "01_output-files/kin_2023_adults_all.csv", row.names=F)

#-------------------------------------------------------------------------------
# Get relatedness of matched mom-kid pairs, use for assigning genetic sires
#-------------------------------------------------------------------------------

kin23.match <- subset(kin23.po.org, 
                      kin23.po.org$FamilyID_ind1==kin23.po.org$FamilyID_ind2)

# filter just moms
kin23.match.mom <- subset(kin23.match, kin23.match$type_ind1=="mom")

# two cases of social mom not being the genetic mom: 
# Urban Farm Girlz 8, brood 1
# Boyer-13 collected eggs, Boyer-05 collected eggs


# keep only high kinship moms
# Now the relatedness range is 0.389-0.473
kin23.match.mom2 <- subset(kin23.match.mom, kin23.match.mom$pi_HAT > 0.1)
range(kin23.match.mom2$pi_HAT)

# save minimum relatedness value for matched moms
mom.kid.min <- min(kin23.match.mom2$pi_HAT)


## use min relatedness cutoff to look at only assigned parents
kin23.assign <- subset(kin23.po.org, kin23.po.org$pi_HAT >= mom.kid.min)

# check how many kids have both, one, and no parents assigned
assign.par <- kin23.assign %>% group_by(Ind2)%>%
  summarise(mom.assign = length(which(type_ind1=="mom")),
            dad.assign = length(which(type_ind1=="dad")))

# 14 kids have 2 moms (!), the rest have just 1
# 37 kids have no dad assigned, 22 have 4, 2 have 2, the rest have just 1

# Check if there are any kids from kin23.po that are not in kin23.assign
# There should be 216 kids total, and 216 are in assign.par
kid216.index <- kin23.po.org$Ind2 %in% kin23.assign$Ind2
# check where index==F, empty dataframe
kid216.check <- kin23.po.org[which(kid216.index==F), ]


## Check kids who have 2 mom assigned

# add column for mom or dad assigned
kin23.assign.count <- kin23.assign %>% group_by(Ind2) %>%
  mutate(mom.assign = length(which(type_ind1=="mom")),
         dad.assign = length(which(type_ind1=="dad")))

double.mom <- subset(kin23.assign.count, kin23.assign.count$mom.assign==2)
length(unique(double.mom$band_ind2)) # 14
unique(double.mom$FamilyID_ind2) 
# "Boyer-15" Boyer-13" "MakeBelieve-70.3" "MakeBelieve-31"

# these are duplicates, with no true multiple maternal assignments
# The same female was at 

## check confusing Boyer mom assignments
# genetic nests were 13 (eggs) and then 15 (chicks)
# Should remove the nest 5 ID from this female, which will reduce duplicate rows
Boyer.mom.57688 <- subset(kin23.assign.count, kin23.assign.count$Ind1=="CO_57688")

# genetic nests were 5 (eggs) and then 13 (1 kid) and 16 (4 kids)
Boyer.mom.57685 <- subset(kin23.assign.count, kin23.assign.count$Ind1=="CO_57685")


## Check confusing Make Believe genetics
MB.mom <- subset(kin23.assign.count, kin23.assign.count$site_ind2=="Make Believe" &
               kin23.assign.count$type_ind1=="mom")
MB.dad <- subset(kin23.assign.count, kin23.assign.count$site_ind2=="Make Believe" &
                   kin23.assign.count$type_ind1=="dad")
# Will need to manually remove duplicates for nest 31 and 70.3 because genetic mom
# is the same at these two, but social and genetic dads are different


## check cases of multiple dads
# Boyer and Make Believe dads are duplicates
# One case of two different males assigned = in nest Struthers-24, male from
# Struthers-24 assigned (57552) AND male from CHR assigned (34050)!
double.dad <- subset(kin23.assign.count, kin23.assign.count$dad.assign>1 &
                     kin23.assign.count$type_ind1=="dad")

# Boyer
Boyer.dad.57686 <- subset(kin23.assign.count, kin23.assign.count$Ind1=="CO_57686")

# Check relatedness between the Struthers and CHR males
double.dad.check <- subset(kin.adults, kin.adults$Ind1=="CO_57552" & 
                             kin.adults$Ind2=="CO_34050" |
                             kin.adults$Ind1=="CO_34050" &
                             kin.adults$Ind2=="CO_57552")
# Their relatedness is high (0.399) so they are either full sibs or father/son
# 57552 in adult data starting 2021, 34050 first there starting 2023
# 34050 was hatched at Struthers in 2022, and his dad is 57552. This makes him 
# full-sibs with the chicks from 2023! 

# check assignment for full Struthers-24 family
# The Struthers-24 pair had multiple chicks together, and not all of those have
# the CHR male assigned as a possible sire...
struthers24 <- subset(kin23.assign.count, 
                      kin23.assign.count$FamilyID_ind2=="Struthers-24")

# check CHR male's relatedness to all Struthers-24 kids
# range of close relation is 0.324-0.487, with kinship to parents 0.391 and 0.395
# This suggest he may be the true sire for 83151 (0.439) and 83149 (0.487),
# although overall range of mom-kid relatedness values is 0.395 - 0.475. 
CHRsire.struthers24 <- subset(kin2.23, kin2.23$Ind1=="CO_34050" &
                                kin2.23$FamilyID_ind2=="Struthers-24" |
                                kin2.23$FamilyID_ind1=="Struthers-24" &
                                kin2.23$Ind2=="CO_34050")

#-------------------------------------------------------------------------------
# Remove duplicated parental assignments
#-------------------------------------------------------------------------------

## Change Boyer family IDs to be more accurate

# familyID should be 13 for clutches 13_collected, and 15_1. 
# social female 57688, social male 57686
kin23.assign2 <- kin23.assign.count

kin23.assign2$FamilyID_ind2[which(
  kin23.assign.count$clutchID_ind2=="Boyer_13_collected" |
  kin23.assign.count$clutchID_ind2=="Boyer_15_1")] <- "Boyer-13"

kin23.assign2$FamilyID_ind1[which(
  kin23.assign.count$Ind1=="CO_57688")] <- "Boyer-13"

# familyID should be 5 for clutches 5_collected, 13.2_1, and 16_2. 
# social female 57685, social male ALSO 57686...
kin23.assign2$FamilyID_ind2[which(
  kin23.assign.count$clutchID_ind2=="Boyer_5_collected" |
    kin23.assign.count$clutchID_ind2=="Boyer_13.2_1" |
    kin23.assign.count$clutchID_ind2=="Boyer_16_2")] <- "Boyer-05"

kin23.assign2$FamilyID_ind1[which(
  kin23.assign.count$Ind1=="CO_57685")] <- "Boyer-05"

boyer13 <- subset(kin23.assign2, kin23.assign2$FamilyID_ind2=="Boyer-13")
boyer05 <- subset(kin23.assign2, kin23.assign2$FamilyID_ind2=="Boyer-05")

# within Boyer-05, change familyID of male 57686 to Boyer-05
kin23.assign2$FamilyID_ind1[
  which(kin23.assign2$Ind1=="CO_57686" & 
          kin23.assign2$FamilyID_ind2=="Boyer-05")] <- "Boyer-05"

boyer05 <- subset(kin23.assign2, kin23.assign2$FamilyID_ind2=="Boyer-05")

# within Boyer-13, change familyID of male 57686 to Boyer-13
kin23.assign2$FamilyID_ind1[
  which(kin23.assign2$Ind1=="CO_57686" & 
          kin23.assign2$FamilyID_ind2=="Boyer-13")] <- "Boyer-13"

boyer13 <- subset(kin23.assign2, kin23.assign2$FamilyID_ind2=="Boyer-13")

# male 57686 also sired chicks in nest 2. Harmonize name for those chicks
kin23.assign2$FamilyID_ind1[
  which(kin23.assign2$Ind1=="CO_57686" & 
          kin23.assign2$FamilyID_ind2=="Boyer-02")] <- "Boyer-13"


## identify and remove replicated duplicates
kin23.assign2$parent.dup <- duplicated(kin23.assign2 %>% 
                                        select(Ind1, Ind2, FamilyID_ind1,
                                               FamilyID_ind2))

kin23.assign3 <- subset(kin23.assign2, kin23.assign2$parent.dup==F)

# count new number of assigned parents
kin23.assign3 <- kin23.assign3 %>% group_by(Ind2) %>%
  mutate(mom.assign2 = length(which(type_ind1=="mom")),
         dad.assign2 = length(which(type_ind1=="dad")))

## Harmonize problematic Make Believe familyIDs
# Male 21085 has family ID both 70 and 70.3 
# Female 57878 has family ID both 31 and 70.3
# female 57596 has family ID only 70
# male 57877 has family ID only 31

MB.doubles <- subset(kin23.assign3, kin23.assign3$mom.assign2>1 &
                       kin23.assign3$site_ind2=="Make Believe" |
                       kin23.assign3$dad.assign2>1 &
                       kin23.assign3$site_ind2=="Make Believe")

# should identify cases where parent band is duplicated, and keep row where 
# familyID of parent matches familyID of kid
dup.index <- kin23.assign3[which(duplicated(kin23.assign3 %>% select(
  Ind1, Ind2))), ]

kin23.assign3$MB.dup <- ifelse(kin23.assign3$Ind1 %in% dup.index$Ind1 &
                                 kin23.assign3$Ind2 %in% dup.index$Ind2, 
                               T, F)
# remove non-matched rows
kin23.assign4 <- kin23.assign3[-which(kin23.assign3$MB.dup==T &
                                        kin23.assign3$FamilyID_ind1 !=
                                        kin23.assign3$FamilyID_ind2), ]

# Manually checked removed rows vs. kept
# dup.check.removed <- kin23.assign3[which(kin23.assign3$MB.dup==T &
#                                             kin23.assign3$FamilyID_ind1 !=
#                                             kin23.assign3$FamilyID_ind2), ]
# 
# dup.check.keep <- subset(kin23.assign4, kin23.assign4$MB.dup==T)
# 
# write.csv(dup.check.removed, "01_output-files/dup.check.removed.csv")
# write.csv(dup.check.keep, "01_output-files/dup.check.keep.csv")

# calculate new number of assigned parents
# Now it's just the three kids from Struthers who have multiple dads assigned
# no kids have multiple moms assigned
kin23.assign4 <- kin23.assign4 %>% group_by(Ind2) %>%
  mutate(mom.assign3 = length(which(type_ind1=="mom")),
         dad.assign3 = length(which(type_ind1=="dad")))

# Check for mis-matched moms
mom.mismatch <- subset(kin23.assign4, kin23.assign4$FamilyID_ind1 !=
                         kin23.assign4$FamilyID_ind2 &
                         kin23.assign4$type_ind1=="mom")
# one case at Urban Farm Girlz
ufg.check <- subset(kin23.assign4, kin23.assign4$site_ind2=="Urban Farm Girlz")

# Incorrect ID obs for female at clutch 8_1
# Correct familyIDs should be: 
# UrbanFarm-07: clutch 7_2, female 57608, male 57643 (change for male)
# UrbanFarm-05: clutches 8_1 and 5_1, female 57629. Male UNK for clutch 8_1,
# male 57765 for clutch 5_1. Clutch 5_1 should actually be named 5_2 (change)
# UrbanFarm-08: clutch 8.2_2, female 57764, male 57976. (change)
# Male 57976 had DNA extracted but is missing from sequence data!

kin23.assign5 <- kin23.assign4

kin23.assign5$FamilyID_ind1[which(
  kin23.assign4$Ind1=="CO_57643")] <- "UrbanFarm-07"

kin23.assign5$FamilyID_ind1[which(
  kin23.assign4$Ind1=="CO_57765")] <- "UrbanFarm-05"

kin23.assign5$FamilyID_ind2[which(
  kin23.assign4$clutchID_ind2=="Urban Farm Girlz_8_1")] <- "UrbanFarm-05"

kin23.assign5$brood_ind2[which(
  kin23.assign4$clutchID_ind2=="Urban Farm Girlz_5_1")] <- "2"


# remove filtering columns that are no longer needed
kin23.assign6 <- kin23.assign5 %>% select(!c(mom.assign, dad.assign, parent.dup,
                                             mom.assign2, dad.assign2, MB.dup))


## save duplicate-free parent-offspring assignments
# save table of assigned parents
write.csv(kin23.assign6, file="01_output-files/kin_2023_assigned_parents.csv", 
          row.names=F)


#-------------------------------------------------------------------------------
# Make plot of mom-offspring vs. dad-offspring relatedness
#-------------------------------------------------------------------------------

## make plot to compare distribution of assigned mom-offspring relatedness to
# distribution of assigned dad-offspring relatedness

# make table for plotting
parent.offspring.plot <- select(kin23.assign6, pi_HAT, type_ind1)
parent.offspring.plot$parent <- ifelse(parent.offspring.plot$type_ind1=="mom",
                                       "Mothers", "Sires")

# set colors to create a legend
cols <- c("Mothers"="maroon", "Sires"="orange")

ggplot(parent.offspring.plot, aes(x=pi_HAT, fill=parent)) +
  geom_histogram(color="black") +
  facet_grid(parent~.) +
  ggtitle("Relatedness values of parents with offspring") +
  xlab("Estimated kinship coefficients between parents and offspring") +
  scale_fill_manual(name="Legend",values=cols)

ggsave("01_output-files/kinship values parents offspring.png", h=3, w=5)

t.test(pi_HAT ~ parent, parent.offspring.plot)
# Welch Two Sample t-test
# 
# data:  pi_HAT by parent
# t = 1.1561, df = 362.99, p-value = 0.2484
# alternative hypothesis: true difference in means between group Mothers and group Sires is not equal to 0
# 95 percent confidence interval:
#   -0.001269540  0.004891512
# sample estimates:
#   mean in group Mothers   mean in group Sires 
# 0.4367731             0.4349622


#-------------------------------------------------------------------------------
# Save table of kids with unknown sires
#-------------------------------------------------------------------------------

dad.unk23 <- kin23.assign7 %>% filter(dad.assign3==0)

# 34 unassigned offspring are from 16 families
length(unique(dad.unk23$FamilyID_ind2))

# save table of kids with unassigned dads
save(dad.unk23, file="01_output-files/kids_unk_dads_2023.Rdata")


#-------------------------------------------------------------------------------
# Calculate summary stats by clutch 1 and 2
#-------------------------------------------------------------------------------

# total offspring - 216
length(unique(kin23.assign6$Ind2))
# total clutches - 60
length(unique(kin23.assign6$clutchID_ind2))

unique(kin23.assign6$brood_ind2)
collected <- kin23.assign6 %>% filter(brood_ind2=="collected")
brood1 <- kin23.assign6 %>% filter(brood_ind2=="1")
brood2 <- kin23.assign6 %>% filter(brood_ind2=="2")

# number of offspring per brood
length(unique(collected$Ind2)) # 77
length(unique(brood1$Ind2)) # 93
length(unique(brood2$Ind2)) # 46

# number of clutches per brood
length(unique(collected$clutchID_ind2)) # 19
length(unique(brood1$clutchID_ind2)) # 26
length(unique(brood2$clutchID_ind2)) # 15

# assigned sires
length(unique(subset(collected$Ind1, collected$type_ind1=="dad"))) # 14
length(unique(subset(brood1$Ind1, brood1$type_ind1=="dad"))) # 26
length(unique(subset(brood2$Ind1, brood2$type_ind1=="dad"))) # 12

# assigned dams
length(unique(subset(collected$Ind1, collected$type_ind1=="mom"))) # 19
length(unique(subset(brood1$Ind1, brood1$type_ind1=="mom"))) # 26
length(unique(subset(brood2$Ind1, brood2$type_ind1=="mom"))) # 15

# number of kids with unassigned sires
dad.unk23.collected <- subset(dad.unk23, dad.unk23$brood_ind2=="collected")
length(unique(dad.unk23.collected$Ind2)) # 7

dad.unk23.brood1 <- subset(dad.unk23, dad.unk23$brood_ind2=="1")
length(unique(dad.unk23.brood1$Ind2)) # 17

dad.unk23.brood2 <- subset(dad.unk23, dad.unk23$brood_ind2=="2")
length(unique(dad.unk23.brood2$Ind2)) # 10
