

# Heather Kenny-Duddela
# Feb 16, 2023
# Updated Aug 4, 2025



# read in lcMLkin output file for updated 2021 data
kin <- read.table("input-files/CO-2021.thin10k.GL-recalc_update2023.relate",
                  header=T)

# read in families table
fam.21 <- read.csv("input-files/fam_clutch21_updated broods 2023 with ids.csv")


# load libraries
library(tidyverse)
library(ggplot2)


#-------------------------------------------------------------------------------

# # add seqID
# fam.21$Seq.ID[fam.21$New2023=="Yes" & 
#                 fam.21$Type!="egg"] <- paste("CO-",substr(
#                   subset(fam.21$Band, fam.21$New2023=="Yes" &
#                            fam.21$Type!="egg"),6,11), sep="")
# 
# # add clutch ID
# fam.21$clutch_id[fam.21$New2023=="Yes" & fam.21$Type=="kid" |
#                    fam.21$New2023=="Yes" & 
#                    fam.21$Type=="egg"] <- paste(
#                      fam.21$Site[fam.21$New2023=="Yes" & fam.21$Type=="kid" |
#                                    fam.21$New2023=="Yes" & fam.21$Type=="egg"],
#                      fam.21$Nest[fam.21$New2023=="Yes" & fam.21$Type=="kid" |
#                                    fam.21$New2023=="Yes" & fam.21$Type=="egg"],
#                      fam.21$Brood[fam.21$New2023=="Yes" & fam.21$Type=="kid" |
#                                     fam.21$New2023=="Yes" & fam.21$Type=="egg"],
#                      sep="_")
# 
# # replace "egg" with "collected" for old clutch.IDs
# fam.21$clutch_id <- gsub("egg","collected", fam.21$clutch_id)
# 
# write.csv(fam.21, "input-files/fam_clutch21_updated broods 2023 with ids.csv",
#           row.names=F)

#-------------------------------------------------------------------------------



# k0_hat = prob no alleles shared that are IBD (unrelated)
# k1_hat = prob 1 allele shared that is IBD
# k2_hat = prob 2 alleles shared that are IBD
# pi_hat = coefficient of relatedness (r)

# coefficient of relatedness (r) = 2*phi = (k1)/2 + k2

# pull out one focal individual
# Blue Cloud 09
kin.83404 <- subset(kin, kin$Ind1=="CO_83404")

# try calculating r for the test individual
# see that the calculation for r matches pi_hat
kin.83404$r <- kin.83404$k1_hat/2 + kin.83404$k2_hat

# plot k0 vs. pi_hat
plot(kin.83404$k0_hat, kin.83404$pi_HAT)


# make plot for all individuals
plot(kin$k0_hat, kin$pi_HAT)


# test birds from 2023 update
# Schaaps-110

# social mom is genetic mom, and social dad is genetic dad!
kin.83218 <- subset(kin, kin$Ind1 == "CO_83218" |
                      kin$Ind2 == "CO_83218")
plot(kin.83218$k0_hat, kin.83218$pi_HAT)

#-------------------------------------------------------------------------------
# # dad 57701 is not in updated kin, but was in original kin21. 
# # Figure out which other samples need to be included
# 
# # load original kin21
# original <- read.table("input-files/troubleshooting/CO-2021.thin10k.GL-recalc.relate",
#                        header=T)
# # load 2021old sample list used for updated lcMLkin analysis
# old <- read.csv("input-files/troubleshooting/sample_list_2021_old.csv")
# 
# 
# ## check files from original that are not in old
# 
# # first make list of unique seqIDs
# original.unique <- unique(c(original$Ind1, original$Ind2))
# 
# missing.old <- as.data.frame(original.unique[which(original.unique %in% old$sample == F)])
# colnames(missing.old) <- "missing.samples"
# write.csv(missing.old, "input-files/troubleshooting/missing from 2021 old.csv")
#-------------------------------------------------------------------------------


# seems like good relatedness cut-off values are 
# <0.9 for k0 and 
# >0.1 for pi_hat


# Add band, site, family, sex, type to kinship table for ind1 and ind2
# this way we can filter to only parent-offspring relationships
# and can identify social vs. extra-pair dams and sires

# replace "-" in seqID with "_"
fam.21$Seq.ID.test <- gsub("-","_",fam.21$Seq.ID)


# relevant info from fam.21
# band, site, nest, type, familyID, brood, clutch_id, seq.ID.test
seq.cols <- fam.21[,c(2,3,5,6,9,10,11,14)]

colnames(seq.cols) <- paste(colnames(seq.cols),"ind1",sep="_")
colnames(seq.cols)[8] <- "Ind1"

kin1 <- left_join(kin, seq.cols, by="Ind1")

seq.cols2 <- fam.21[,c(2,3,5,6,9,10,11,14)]

colnames(seq.cols2) <- paste(colnames(seq.cols2),"ind2",sep="_")
colnames(seq.cols2)[8] <- "Ind2"

kin2 <- left_join(kin1, seq.cols2, by="Ind2")

# save kinship table with all pairwise comparisons and labels
write.csv(kin2, "01_output-files/kin_2021_all_with_labels.csv", row.names=F)

#-------------------------------------------------------------------------------
# make table of all possible kid-parent combinations
#-------------------------------------------------------------------------------

# keep all possible kid-parent combinations from the table
# this object is the same size as kin.kids.par
kin.po <- subset(kin2, kin2$Type_ind1=="kid" & kin2$Type_ind2=="mom" |
                     kin2$Type_ind1=="kid" & kin2$Type_ind2=="dad" |
                     kin2$Type_ind1=="egg" & kin2$Type_ind2=="mom" |
                     kin2$Type_ind1=="egg" & kin2$Type_ind2=="dad" |
                     kin2$Type_ind1=="mom" & kin2$Type_ind2=="egg" |
                     kin2$Type_ind1=="mom" & kin2$Type_ind2=="kid" |
                     kin2$Type_ind1=="dad" & kin2$Type_ind2=="egg" |
                     kin2$Type_ind1=="dad" & kin2$Type_ind2=="kid")

## restructure so that ind1 are all parents and ind2 are all kids

## cases where ind1 is kid or egg
kin.po.ind1kid <- subset(kin.po, kin.po$Type_ind1=="kid" |
                             kin.po$Type_ind1=="egg")
# keep only ind1 columns
kin.po.ind1kid.before <- kin.po.ind1kid[,1:14]
# change ind1 to ind2
colnames(kin.po.ind1kid.before) <- gsub("ind1", "ind2", 
                                          colnames(kin.po.ind1kid.before))
colnames(kin.po.ind1kid.before)[1:2] <- c("Ind2","Ind1")
# reorder cols
kin.po.ind1kid.after <- relocate(kin.po.ind1kid.before, Ind2, .after=Ind1)

# keep only ind2 columns
kin.po.ind2par.before <- kin.po.ind1kid[,c(1:7, 15:21)]
# change ind2 to ind1
colnames(kin.po.ind2par.before) <- gsub("ind2", "ind1", 
                                          colnames(kin.po.ind2par.before))
colnames(kin.po.ind2par.before)[1:2] <- c("Ind2","Ind1")
# reorder cols
kin.po.ind2par.after <- relocate(kin.po.ind2par.before, Ind2, .after=Ind1)

# combine ind1 and ind2
kin.po.ind1kid.swap <- left_join(kin.po.ind2par.after, 
                                   kin.po.ind1kid.after, 
                                   by=c("Ind1", "Ind2", "k0_hat", "k1_hat",
                                        "k2_hat", "pi_HAT", "nbSNP"))

# Add the swapped table to the "good" cases from the original table

# cases where ind1 is parents
kin.po.ind1par <- subset(kin.po, kin.po$Type_ind1=="mom" |
                             kin.po$Type_ind1=="dad")
# rowbind the two tables
kin.po.org <- rbind(kin.po.ind1par, kin.po.ind1kid.swap)

# check that organization worked
unique(kin.po.org$Type_ind1)
unique(kin.po.org$Type_ind2)

## check number of kids, moms, and dads
# 224 kids and eggs
length(unique(kin.po.org$Ind2)) 
# 38 dads
length(unique(subset(kin.po.org$Ind1, kin.po.org$Type_ind1=="dad")))
# 27 moms
length(unique(subset(kin.po.org$Ind1, kin.po.org$Type_ind1=="mom")))


# save table of all parent-offspring kinship
write.csv(kin.po.org, file="01_output-files/kin_2021_parent_offspring_all.csv", 
          row.names=F)

#-------------------------------------------------------------------------------
# make table of adult only relationships
#-------------------------------------------------------------------------------

kin.adults <- subset(
  kin2,
  kin2$Type_ind1=="mom" & kin2$Type_ind2=="dad" |
    kin2$Type_ind1=="dad" & kin2$Type_ind2=="mom" |
    kin2$Type_ind1=="mom" & kin2$Type_ind2=="mom" |
    kin2$Type_ind1=="dad" & kin2$Type_ind2=="dad")

# 3 cases of close relatives between adults (pi_hat > 0.1)

# save file
write.csv(kin.adults, "01_output-files/kin_2021_adults_all.csv", row.names=F)

#-------------------------------------------------------------------------------
# Get relatedness of matched mom-kid pairs, use for assigning genetic sires
#-------------------------------------------------------------------------------

kin.match <- subset(kin.po.org, 
                      kin.po.org$FamilyID_ind1==kin.po.org$FamilyID_ind2)

# filter just moms (pi_hat range 0.315-0.571)
kin.match.mom <- subset(kin.match, kin.match$Type_ind1=="mom")

# save minimum relatedness value for matched moms
mom.kid.min <- min(kin.match.mom$pi_HAT)

## use min relatedness cutoff to look at only assigned parents
kin.assign <- subset(kin.po.org, kin.po.org$pi_HAT >= mom.kid.min)

# check how many kids have both, one, and no parents assigned
assign.par <- kin.assign %>% group_by(Ind2)%>%
  summarise(mom.assign = length(which(Type_ind1=="mom")),
            dad.assign = length(which(Type_ind1=="dad")))
# all kids have exactly 1 mom assigned. 37 kids have no dad, the rest have
# exactly 1 dad

# 222 kids in assign.par, and should be 224 total
# check for kids from kin.po.org that are not in assign.par
unique(kin.po.org$Ind2[which(kin.po.org$Ind2 %in% assign.par$Ind2 == F)])
# two missing kids: "CO_39180" "CO_39179" from Blue Cloud-17
# no blood was collected from the mom, and the dad was UNB

## save parent-offspring assignments
write.csv(kin.assign, "01_output-files/kin_2021_assigned_parents.csv", 
          row.names = F)

#-------------------------------------------------------------------------------
# Make plot of mom-offspring vs. dad-offspring relatedness
#-------------------------------------------------------------------------------

## make plot to compare distribution of assigned mom-offspring relatedness to
# distribution of assigned dad-offspring relatedness

# make table for plotting
parent.offspring.plot <- select(kin.assign, pi_HAT, Type_ind1)
parent.offspring.plot$parent <- ifelse(parent.offspring.plot$Type_ind1=="mom",
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
# t = 0.94114, df = 399.85, p-value = 0.3472
# alternative hypothesis: true difference in means between group Mothers and group Sires is not equal to 0
# 95 percent confidence interval:
#   -0.002602467  0.007382647
# sample estimates:
#   mean in group Mothers   mean in group Sires 
# 0.4159414             0.4135514

#-------------------------------------------------------------------------------
# Save table of kids with unknown sires
#-------------------------------------------------------------------------------

dad0 <- subset(assign.par, assign.par$dad.assign==0)

dad.unk <- kin.assign %>% filter(Ind2 %in% dad0$Ind2)

# 37 unassigned offspring are from 10 families (excluding 2 from Blue Cloud 17)
length(unique(dad.unk$FamilyID_ind2))

# save table of kids with unassigned dads
save(dad.unk, file="01_output-files/kids_unk_dads_2021.Rdata")


# Check extra-pair kids 
kid.ep <- subset(kin.assign, kin.assign$FamilyID_ind2 != 
                   kin.assign$FamilyID_ind1 &
                   kin.assign$Type_ind1=="dad")




