
# Heather Kenny-Duddela
# Feb 24, 2023

# Paternity assignment using lcMLkin output
# For 2022 data

# read in lcMLkin output file
# approx 90,000 SNP used, filtered for LD at 10kb windows
# Only autosomes, sex-linked scaffolds removed from dataset
# BARS genome is 1.2GB total size

kin22 <- read.table("input-files/CO-2022.thin10K.GL-recalc.relate", header=T)

# read in families table with 1st and 2nd brood info

fam.22 <- read.csv("input-files/fam_clutch_2022_updated broods.csv")



# load libraries
library(tidyverse)
library(ggplot2)



# k0_hat = prob no alleles shared that are IBD (unrelated)
# k1_hat = prob 1 allele shared that is IBD
# k2_hat = prob 2 alleles shared that are IBD
# pi_hat = coefficient of relatedness (r)

# coefficient of relatedness (r) = 2*phi = (k1)/2 + k2

# pull out one focal individual
# CHR 004
kin22.83277 <- subset(kin22, kin22$Ind1=="CO_83277" |
                        kin22$Ind2=="CO_83277")

# try calculating r for the test individual
# see that the calculation for r matches pi_hat
kin22.83277$r <- kin22.83277$k1_hat/2 + kin22.83277$k2_hat

# plot k0 vs. pi_hat
plot(kin22.83277$k0_hat, kin22.83277$pi_HAT)


# make plot for all individuals
plot(kin22$k0_hat, kin22$pi_HAT)

# more test birds
# blue Cloud 22
kin22.83079 <- subset(kin22, kin22$Ind1=="CO_83079" |
                        kin22$Ind2=="CO_83079")
plot(kin22.83079$k0_hat, kin22.83079$pi_HAT)

# cooks 29 
kin22.83069 <- subset(kin22, kin22$Ind1=="CO_83069" |
                        kin22$Ind2=="CO_83069")
plot(kin22.83069$k0_hat, kin22.83069$pi_HAT )


# seems like good relatedness cut-off values are 
# <0.9 for k0 and 
# >0.1 for pi_hat


# Add band, site, family, sex, type, clutch to kinship table for ind1 and ind2



# replace "-" in seqID with "_"
fam.22$Seq.ID.test <- gsub("-","_",fam.22$Seq.ID)


# remove duplicates from fam22, or where no famID is assigned
fam.duplicate <- anyDuplicated(fam.22[,c(3,7)])

# pull out relevant metadata to add to kinship data
seq.cols <- fam.22[,c(3,4,6,7,19:22)]

seq.cols[duplicated(seq.cols$Ind1),]

# two males are duplicated because they were attending two nests during
# first and/or second broods. I think this accounts for the difference
# in number of rows between kin22 and kin22.1
#seq.dup <- subset(seq.cols, seq.cols$Ind1=="CO_57806" |
#                    seq.cols$Ind1=="CO_57582")

# add suffix to specify ind1 metadata
colnames(seq.cols) <- paste(colnames(seq.cols),"ind1",sep="_")
colnames(seq.cols)[8] <- "Ind1"

# combine kinship data with metadata for ind1
kin1.22 <- left_join(kin22, seq.cols, by="Ind1")

# check whether merge worked properly
# note that only kids have nest, clutch, and brood metadata. Adults don't have this
# because site and family ID are sufficient to specify adults

test <- subset(kin1.22, kin1.22$Ind1=="CO_57806" | kin1.22$Ind2=="CO_57806")

test2 <- subset(kin1.22, kin1.22$Ind1=="CO_57649" | kin1.22$Ind2=="CO_57649")
test3 <- subset(kin1.22, kin1.22$Ind1=="CO_83428" | kin1.22$Ind2=="CO_83428")

test4 <- subset(kin1.22, kin1.22$Ind1=="CO_57582" | kin1.22$Ind2=="CO_57582")


# repeat process of adding metadata for ind2

seq.cols2 <- fam.22[,c(3,4,6,7,19:22)]

colnames(seq.cols2) <- paste(colnames(seq.cols2),"ind2",sep="_")
colnames(seq.cols2)[8] <- "Ind2"

kin2.22 <- left_join(kin1.22, seq.cols2, by="Ind2")

test5 <- subset(kin2.22, kin2.22$Ind1=="CO_57806" | kin2.22$Ind2=="CO_57806")

# save kinship table with metadata
# this has all pairwise comparisons, including parent-kid, parent-parent, kid-kid
write.csv(kin2.22, "generated-files/kin_2022_all_with_labels.csv", row.names=F)
save(kin2.22, file="generated-files/kin2_22.Rdata")

# make list of kids that were sequenced from 2022
kids.sequenced <- subset(kin2.22, kin2.22$Type_ind1=="kid" &
                           kin2.22$Type_ind2=="kid")
kids.sequenced2 <- kids.sequenced[, c(1,8:14)]
kids.sequenced3 <- kids.sequenced[, c(2,15:21)]
colnames(kids.sequenced2) <- c("SeqID", "Band","Site", "Type", "FamilyID",
                               "Nest","Brood","clutchID")
colnames(kids.sequenced3) <- c("SeqID", "Band","Site", "Type", "FamilyID",
                               "Nest","Brood","clutchID")
kids.sequenced4 <- rbind(kids.sequenced2, kids.sequenced3)
kids.sequenced5 <- kids.sequenced4[-(which(duplicated(kids.sequenced4$Band))), ]
# save
write.csv(kids.sequenced5, "generated-files/list of sequenced nestlings 2022.csv", row.names = F)



# Pull out all possible kid-parent combinations from the table
# due to structure of table, ind1 is all mom/dad and ind2 is all kid
kin22.po <- subset(kin2.22, kin2.22$Type_ind1=="kid" & kin2.22$Type_ind2=="mom" |
                   kin2.22$Type_ind1=="kid" & kin2.22$Type_ind2=="dad" |
                   kin2.22$Type_ind1=="egg" & kin2.22$Type_ind2=="mom" |
                   kin2.22$Type_ind1=="egg" & kin2.22$Type_ind2=="dad" |
                   kin2.22$Type_ind1=="mom" & kin2.22$Type_ind2=="egg" |
                   kin2.22$Type_ind1=="mom" & kin2.22$Type_ind2=="kid" |
                   kin2.22$Type_ind1=="dad" & kin2.22$Type_ind2=="egg" |
                   kin2.22$Type_ind1=="dad" & kin2.22$Type_ind2=="kid")

# save table of all parent-offspring kinship
write.csv(kin22.po, file="generated-files/kin_2022_parent_offspring_all.csv", row.names=F)

# pull out test dads
# Looks like 325 kids total
length(unique(kin22.po$Ind2))

# number of unique males (103)
length(unique(subset(kin22.po$Ind1, kin22.po$Type_ind1=="dad")))

# number of unique females (60)
length(unique(subset(kin22.po$Ind1, kin22.po$Type_ind1=="mom")))

# CHR-43, all offspring are within-pair, he is dad to all kids in the family
dad.57567 <- subset(kin22.po, kin22.po$Ind1=="CO_57567")
plot(dad.57567$k0_hat, dad.57567$pi_HAT)

dad.97548 <- subset(kin22.po, kin22.po$Ind1=="CO_97548")
plot(dad.97548$k0_hat, dad.97548$pi_HAT)


# filter so k0<0.9 and pi_hat > 0.1
kin22.close <- subset(kin22.po, kin22.po$k0_hat<0.9 & kin22.po$pi_HAT>0.1)

# save file
write.csv(kin22.close, file="generated-files/kin_2022_parent_offspring_close_relatives.csv")

# check that test dads retain the expected offspring
dad.assign.57567 <- subset(kin22.close, kin22.close$Ind1=="CO_57567")

dad.assign.97548 <- subset(kin22.close, kin22.close$Ind1=="CO_97548")


# Identify cases where parent familyID does not match kid familyID (extra-pair sire)
index22 <- which(kin22.close$FamilyID_ind1 != kin22.close$FamilyID_ind2)

# look at mismatched rows for EP offspring
# 197 pairs total
mismatch22 <- kin22.close[index22,]

# some mismatched parents have pi_hat values lower than 0.3
# I think this is because they are grandparents or aunts/uncles instead of parents
plot(mismatch22$k0_hat, mismatch22$pi_HAT)


# try subsetting to test specific sib/PO relationships between different adults
# first get table of just moms and dads
adults22 <- subset(kin2.22, kin2.22$Type_ind1=="mom" & kin2.22$Type_ind2=="dad" |
                    kin2.22$Type_ind1=="mom" & kin2.22$Type_ind2=="mom" |
                    kin2.22$Type_ind1=="dad" & kin2.22$Type_ind2=="mom" |
                    kin2.22$Type_ind1=="dad" & kin2.22$Type_ind2=="dad")

# save file
write.csv(adults22, file="generated-files/kin_2022_adults_all.csv", row.names=F)

# now filter for close relatives
adults22.close <- subset(adults22, adults22$k0_hat<0.9 & adults22$pi_HAT>0.1)

# check instances of mismatched sites for po relations (different site EP)
# 91 kids total!
mismatch.site22 <- kin22.close[which(kin22.close$Site_ind1 != kin22.close$Site_ind2), ]

# look at matched rows (within-pair sires)
match22 <- kin22.close[which(kin22.close$FamilyID_ind1 == kin22.close$FamilyID_ind2), ]

# all matched parents have pi_hat values greater than 0.3
plot(match22$k0_hat, match22$pi_HAT)

# check range of matched moms only
range(subset(match22$pi_HAT, match22$Type_ind1=="mom"))# 0.384 to 0.442

# save minimum relatedness value for matched moms
mom.kid.min <- min(subset(match22$pi_HAT, match22$Type_ind1=="mom"))

# check how many kids have both, one, and no parents assigned
assign.par22 <- kin22.close %>% group_by(Ind2) %>%
  summarise(mom.assign = length(which(Type_ind1=="mom")),
            dad.assign = length(which(Type_ind1 == "dad")))

# check which kid from kin22.po is not in kin22.close
kid325.index <- kin22.po$Ind2 %in% kin22.close$Ind2
# check where index==F
kin325.check <- kin22.po[which(kid325.index==FALSE),]
# only one kid is missing, and this was a fledgling who was caught at CHR
# but was never assigned to a nest. True kid sample size is 324


# check cases where >1 dads are assigned to the same kid
dad.2.index22 <- which(kin22.close$Ind2 %in% assign.par22$Ind2[which(assign.par22$dad.assign>1)])

# some "dads" must be other close relatives
dad.2.22 <- kin22.close[dad.2.index22, ]

hist(dad.2.22$pi_HAT)



# filter close relatives to keep only highest of each mom and dad
# group  by each kid and parent type (mom and dad)
kin.highest.par22 <- kin22.close %>%
  group_by(Ind2, Type_ind1) %>%
  summarise(par.hi = Ind1[which(pi_HAT == max(pi_HAT))],
            pi_HAT_max = pi_HAT[which(pi_HAT == max(pi_HAT))])

kin.hi <- kin22.close %>%
  group_by(Ind2, Type_ind1) %>%
  summarise(par.hi = Ind1[which(pi_HAT == max(pi_HAT))])

colnames(kin.hi)[3] <- "Ind1"

# add metadata back to list of highest value moms and dads
kin.hi.info <- left_join(kin.hi[,c(1,3)], kin22.close, by=c("Ind1","Ind2"))

# make sure that using max dad is not cutting out close calls that we should pay attention to
# kin.cutoff.dad has fewer rows than kin.hi.info
kin.cutoff.dad <- kin22.close %>%
  group_by(Ind2, Type_ind1, Band_ind2) %>%
  summarise(Ind1 = Ind1[which(pi_HAT >= mom.kid.min)],
            pi_HAT = pi_HAT[which(pi_HAT >= mom.kid.min)])

# check which rows are missing from kin.cutoff.dad
kin.missing.kid <- kin.hi.info$Ind2 %in% kin.cutoff.dad$Ind2
index.missing.kid <- which(kin.missing.kid != TRUE)
# index missing kid is empty so all kids are accounted for

# check for missing parents
kin.missing.par <- kin.hi.info$Ind1 %in% kin.cutoff.dad$Ind1
index.missing.par <- which(kin.missing.par != TRUE)
# row 288 of kin.hi.info, this is dad CO_57665 from Cooks, and his pi_hat is 0.249
kin.hi.info[288,]

## Use kin.cutoff.dad because it is more conservative

# check number of moms and dads kept for each kid
# there are still some kids with 2 dads assigned
kin.cutoff.summary <- kin.cutoff.dad %>%
  group_by(Ind2, Type_ind1) %>%
  summarise(num.par = n() )

kin.cutoff.summary2 <- kin.cutoff.dad %>%
  group_by(Ind2, Band_ind2) %>%
  summarise(mom.assign = length(which(Type_ind1=="mom")),
            dad.assign = length(which(Type_ind1 == "dad")))

# check number of kids with unassigned dads = 38 kids
length(which(kin.cutoff.summary2$dad.assign==0))

# number of kids with unassigned moms = 8
# one clutch is Cooks-13 (kids 34060 to 34063) where the mom ID was suspected but not confirmed
# other clutch is CHR-011 (83142 to 83145) where the mom was UNB and never caught
length(which(kin.cutoff.summary2$mom.assign==0))

# pull out kids with 2 dads after new filtering
# these are only cases where 2 dads are kept but they are the same individual
# These are the males at CHR that were attending two nests, so they are in the data set twice 
# to keep both sets of family IDs
kin.cutoff.dad2 <- kin.cutoff.dad[kin.cutoff.dad$Ind2 %in% subset(kin.cutoff.summary$Ind2, kin.cutoff.summary$num.par>1) ,]
kin.cutoff.dad2.only <- subset(kin.cutoff.dad2, kin.cutoff.dad2$Type_ind1=="dad")

# add metadata to table of assigned parents using the cutoff value
kin.cutoff.info <- left_join(kin.cutoff.dad[,c(1,4)], kin22.close, by=c("Ind1","Ind2"))


# check for mismatched moms
mom.mismatch <- subset(kin.cutoff.info, kin.cutoff.info$FamilyID_ind1 != kin.cutoff.info$FamilyID_ind2 &
                         kin.cutoff.info$Type_ind1 == "mom")

# check total number of moms
mom.tot <- subset(kin.cutoff.info, kin.cutoff.info$Type_ind1=="mom")

# the only "mismatched" moms are at Boyer where the mom family ID includes both 5 and 7, 
# but kid family ID just includes 5
mom.match <- subset(kin.cutoff.info, kin.cutoff.info$FamilyID_ind1 == kin.cutoff.info$FamilyID_ind2 &
                      kin.cutoff.info$Type_ind1=="mom")




# save table of assigned parents
write.csv(kin.cutoff.info, file="generated-files/kin_2022_assigned_parents.csv", row.names=F)
# save as Rdata
save(kin.cutoff.info, file="generated-files/kin_cutoff_info.Rdata")

## make plot to compare distribution of assigned mom-offspring relatedness to
# distribution of assigned dad-offspring relatedness

# make table for plotting
parent.offspring.plot <- select(kin.cutoff.info, pi_HAT, Type_ind1)
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

ggsave("generated-files/kinship values parents offspring.png", h=3, w=5)

t.test(pi_HAT ~ parent, parent.offspring.plot)
# Welch Two Sample t-test
# 
# data:  pi_HAT by parent
# t = 2.6228, df = 632.99, p-value = 0.008932
# alternative hypothesis: true difference in means between group Mothers and group Sires is not equal to 0
# 95 percent confidence interval:
#   0.0005203446 0.0036211577
# sample estimates:
#   mean in group Mothers   mean in group Sires 
# 0.4171899             0.4151191 




# check again for extra-pair sires after filtering extra dads
kin.ep <- kin.cutoff.info[which(kin.cutoff.info$FamilyID_ind1 != kin.cutoff.info$FamilyID_ind2), ]

# save list of ep dads
write.csv(kin.ep, file="generated-files/kin_2022_ep_dads.csv", row.names=F)

# check between site ep fertilizations
kin.ep.site <- kin.cutoff.info[which(kin.cutoff.info$FamilyID_ind1 != kin.cutoff.info$FamilyID_ind2 &
                                   kin.cutoff.info$Site_ind1 != kin.cutoff.info$Site_ind2), ]

# save between site ep dads
write.csv(kin.ep.site, file="generated-files/kin_2022_ep_dads_diff_site.csv", row.names=F)

# filter to only keep matches >0.3
# this is no different from the current kin.ep.site object
kin.ep.site.strict <- subset(kin.ep.site, kin.ep.site$pi_HAT>0.3)


# check number of known assigned dads and moms
# assigned dads = 62
length(unique(subset(kin.cutoff.info$Band_ind1, kin.cutoff.info$Type_ind1=="dad")))
# assigned moms = 57
length(unique(subset(kin.cutoff.info$Band_ind1, kin.cutoff.info$Type_ind1=="mom")))


### check number of unknown dads

# pull out all kids from families table
kid.bands <- subset(fam.22, fam.22$Type=="kid")
length(unique(kid.bands$Band)) # 326 kids!

# only 324 kids in the kin.cutoff.summary2 - figure out who is missing
missing.kids.index <- which(kid.bands$Band %in% kin.cutoff.summary2$Band_ind2)
missing.kids <- kid.bands[-(missing.kids.index),]

# check for missing kids in the kinship table
# highest pi_hat is 0.092 with mom from Speedwell - no relatives

# okay because this was a juvenile caught at CHR
kin22.64911 <- subset(kin2.22, kin2.22$Ind1=="CO_64911" | kin2.22$Ind2=="CO_64911")

# other missing sample was too degraded so not included in kinship analysis
kin22.T017 <- subset(kin2.22, kin2.22$Ind1=="CO_T017" | kin2.22$Ind2=="CO_T017")

# all other kids are represented in the kin.cutoff.summary2 table, and have at least mom or dad assigned


length(which(kin.cutoff.summary2$dad.assign==0))

# pull out offspring with unassigned dads
dad.unk22.index <- which(kin22.close$Ind2 %in% 
                           subset(kin.cutoff.summary2$Ind2, kin.cutoff.summary2$dad.assign<1))

dad.unk22 <- kin22.close[dad.unk22.index, ]

# the 38 unassigned offspring are from 16 different families
length(unique(dad.unk22$Ind2))
length(unique(dad.unk22$FamilyID_ind2))

# save table of kids with unassigned, unknown dads
save(dad.unk22, file="generated-files/kids_unk_dads_2022.Rdata")



################################################################################
# calculate summary stats by clutch 1 and 2

# total offspring, assigned sires, unassigned sires, unique clutches
# assigned dams, unassigned dams

length(unique(kin.cutoff.info$Ind2))# 324
length(unique(kin.cutoff.info$clutch_id_ind2)) # 92

unique(kin.cutoff.info$Brood.._ind2)

brood1 <- subset(kin.cutoff.info, kin.cutoff.info$Brood.._ind2=="1")
brood2 <- subset(kin.cutoff.info, kin.cutoff.info$Brood.._ind2=="2")
brood3 <- subset(kin.cutoff.info, kin.cutoff.info$Brood.._ind2=="3")


# number offspring
length(unique(brood1$Ind2)) # 176
length(unique(brood2$Ind2)) # 146
length(unique(brood3$Ind2)) # 2

# number assigned sires
length(unique(subset(brood1$Ind1, brood1$Type_ind1=="dad"))) # 48
length(unique(subset(brood2$Ind1, brood2$Type_ind1=="dad"))) # 42
length(unique(subset(brood3$Ind1, brood3$Type_ind1=="dad"))) # 1

# number assigned dams
length(unique(subset(brood1$Ind1, brood1$Type_ind1=="mom"))) # 49
length(unique(subset(brood2$Ind1, brood2$Type_ind1=="mom"))) # 40
length(unique(subset(brood3$Ind1, brood3$Type_ind1=="mom"))) # 1

# number of clutches
length(unique(brood1$clutch_id_ind2)) # 50
length(unique(brood2$clutch_id_ind2)) # 41
length(unique(brood3$clutch_id_ind2)) # 1

# number unassigned sires
dad.unk22.brood1 <- subset(dad.unk22, dad.unk22$Brood.._ind2=="1")
length(unique(dad.unk22.brood1$Ind2)) # 16

dad.unk22.brood2 <- subset(dad.unk22, dad.unk22$Brood.._ind2=="2")
length(unique(dad.unk22.brood2$Ind2)) # 22






