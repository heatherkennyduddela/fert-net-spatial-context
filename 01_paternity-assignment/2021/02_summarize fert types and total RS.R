
################################################################################
# Summarise fertiliztion types and total RS for females and males in 2021
# with updated samples from 2023 sequencing effort (added Schaaps)
# Heather Kenny-Duddela
# Aug 21, 2025
################################################################################

# load data
po.assign <- read.csv("01_output-files/kin_2021_assigned_parents.csv")

# load libraries
library(tidyverse)

#-------------------------------------------------------------------------------
# switch table to wide format where each row is a kid and there are columns
# for mom and dad
#-------------------------------------------------------------------------------

# switch to wide format where each row is a kid, and there are cols for dad and mom
# first just do dads
po.dad.wide <- subset(po.assign, po.assign$Type_ind1=="dad")
# pull out only cols that need to get "dad" labels added
po.dad.wide2 <- po.dad.wide[,c(2,3:11)]
colnames(po.dad.wide2)[2:10] <- c("k0_hat_dad", "k1_hat_dad", "k2_hat_dad", 
                                  "pi_HAT_dad",
                                  "nbSNP_dad", "Band_dad","Site_dad", 
                                  "Type_dad" , "FamilyID_dad")

# do same for moms
po.mom.wide <- subset(po.assign, po.assign$Type_ind1=="mom")
po.mom.wide2 <- po.mom.wide[,c(2,3:11)]
colnames(po.mom.wide2)[2:10] <- c("k0_hat_mom", "k1_hat_mom", "k2_hat_mom", 
                                  "pi_HAT_mom",
                                  "nbSNP_mom", "Band_mom","Site_mom", 
                                  "Type_mom" , "FamilyID_mom")

# pull out columns with relevant kid info, and no duplicated kids
po.kids <- po.assign[-which(duplicated(po.assign$Ind2)),c(2,15:21)]

# add dad columns to kid table
po.wide.d <- left_join(po.kids, po.dad.wide2, by="Ind2")

# add mom columns to kid table
po.wide <- left_join(po.wide.d, po.mom.wide2, by="Ind2")

# create column to specify genetic family, kids with same genetic mom and dad
po.wide$genetic_fam <- paste(po.wide$Band_dad, po.wide$Band_mom, sep="_")

# move genetic family label near the front of the table
po.wide <- po.wide[,c(1, 27, 2:26)]

# check for duplicated rows due to some males having multiple family IDs
# identify duplicated rows across kid ID, kid family, and dad band
po.wide.dup <- po.wide[duplicated(po.wide[,c(1,6,15)]),]

# add column for categorical site size
po.wide$site_cat <- NA
po.wide$site_cat[which(po.wide$Site_ind2=="Cathys" |
                         po.wide$Site_ind2=="Reinarz" |
                         po.wide$Site_ind2=="Struthers" |
                         po.wide$Site_ind2=="Urban Farm Girlz")] <- "small"
po.wide$site_cat[which(po.wide$Site_ind2=="Schaaps" |
                         po.wide$Site_ind2=="Blue Cloud" |
                         po.wide$Site_ind2=="Make Believe" |
                         po.wide$Site_ind2=="Cooks")] <- "medium"
po.wide$site_cat[which(po.wide$Site_ind2=="Hepp (near Cook)")] <- "solitary"



## save wide format table
write.csv(po.wide, "02_output-files/kin_2021_parent_offspring_assigned_wide.csv",
          row.names = F)

#-------------------------------------------------------------------------------
# Count number of offspring produced by each genetic family
#-------------------------------------------------------------------------------

# count up number of offspring produced by each genetic family
shared.fert.list <- po.wide %>% 
  group_by(genetic_fam, Band_dad, Band_mom, 
           Site_dad, Site_mom, FamilyID_dad, FamilyID_mom, site_cat) %>%
  summarise(fert = n())

hist(shared.fert.list$fert)

# add labels to specify fertilization type: within-pair (wp), extra-pair same site (ep_same), 
# extra-pair different site (ep_diff), unknown dad (dad_unk), unknown mom (mom_unk)
# based on matching mom and dad FamilyIDs and Site 

shared.fert.list$fert_type <- NA
shared.fert.list$fert_type[which(shared.fert.list$FamilyID_dad == 
                                   shared.fert.list$FamilyID_mom)] <- "wp"
shared.fert.list$fert_type[which(shared.fert.list$FamilyID_dad != 
                                   shared.fert.list$FamilyID_mom &
                                   shared.fert.list$Site_dad == 
                                   shared.fert.list$Site_mom)] <- "ep_same"
shared.fert.list$fert_type[which(shared.fert.list$FamilyID_dad != 
                                   shared.fert.list$FamilyID_mom &
                                   shared.fert.list$Site_dad != 
                                   shared.fert.list$Site_mom)] <- "ep_diff"
shared.fert.list$fert_type[which(is.na(shared.fert.list$FamilyID_dad))] <- "ep_unsampled"

# order factors
shared.fert.list$fert_type <- factor(shared.fert.list$fert_type, 
                                  levels=c("wp","ep_same","ep_unsampled",
                                           "ep_diff"))

# add column for year
shared.fert.list$year <- 2021


# save shared fertilizations list
write.csv(shared.fert.list, 
          file="02_output-files/fert_2021_by_genetic_family.csv", 
          row.names=F)

# summarise number and proportion of each offspring type
total.fert <- sum(shared.fert.list$fert)

fert.type.summary <- shared.fert.list %>%
  group_by(fert_type) %>%
  summarise(num.fert = sum(fert),
            prop.fert = num.fert/total.fert)


# table of just between-site ferts (Blue Cloud 18 male to Make  Believe 
# 69 and Struthers 20)
ep.diff <- subset(shared.fert.list, shared.fert.list$fert_type=="ep_diff")

# check for ferts between out buildings within the same site

# Blue Cloud (shed/washroom male to main barn 14)
ep.BlueCloud <- subset(shared.fert.list, shared.fert.list$fert_type=="ep_same" &
                         shared.fert.list$Site_mom=="Blue Cloud")

# Struthers (none)
ep.Struthers <- subset(shared.fert.list, shared.fert.list$fert_type=="ep_same" &
                         shared.fert.list$Site_mom=="Struthers")

# Make Believe (no ferts from shed male to main barn)
ep.MakeBelieve <- subset(shared.fert.list, shared.fert.list$fert_type=="ep_same" &
                           shared.fert.list$Site_mom=="Make Believe")

# Urban Farm (none)
ep.UFG <- subset(shared.fert.list, shared.fert.list$fert_type=="ep_same" &
                   shared.fert.list$Site_mom=="Urban Farm Girlz")

#-------------------------------------------------------------------------------
# Visualize fertilization types
#-------------------------------------------------------------------------------

## By female ID
ggplot(shared.fert.list, aes(x=reorder(Band_mom, fert),
                             y=fert, fill=fert_type)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() + ggtitle("Number of fertilization types by female ID") +
  xlab("Nest site") + ylab("Number of\n fertilizations")

# by site size
ggplot(shared.fert.list, aes(x=site_cat,
                             y=fert, fill=fert_type)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() + ggtitle("Number of fertilization types by site size") +
  xlab("Site size category") + ylab("Proprtion of\n fertilizations")

# by site
ggplot(shared.fert.list, aes(x=Site_mom,
                             y=fert, fill=fert_type)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() + ggtitle("Number of fertilization types by site") +
  xlab("Nest site") + ylab("Number of\n fertilizations")

# by year
ggplot(shared.fert.list, aes(x=as.factor(year),
                             y=fert, fill=fert_type)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d() + ggtitle("Number of fertilization types") +
  xlab("Year") + ylab("Number of\n fertilizations")


#-------------------------------------------------------------------------------
# Count number of offspring for each genetic family by clutch
#-------------------------------------------------------------------------------


# summarise at clutch level
shared.fert.clutch <- po.wide %>% 
  group_by(clutch_id_ind2, Band_dad, Band_mom, 
           Site_dad, Site_mom, FamilyID_dad, FamilyID_mom, Brood_ind2, 
           site_cat) %>%
  summarise(fert = n())

# number of clutches (56)
length(unique(shared.fert.clutch$clutch_id_ind2))

shared.fert.clutch$fert_type <- NA
shared.fert.clutch$fert_type[which(shared.fert.clutch$FamilyID_dad == 
                                   shared.fert.clutch$FamilyID_mom)] <- "wp"
shared.fert.clutch$fert_type[which(shared.fert.clutch$FamilyID_dad != 
                                   shared.fert.clutch$FamilyID_mom &
                                   shared.fert.clutch$Site_dad == 
                                   shared.fert.clutch$Site_mom)] <- "ep_same"
shared.fert.clutch$fert_type[which(shared.fert.clutch$FamilyID_dad != 
                                   shared.fert.clutch$FamilyID_mom &
                                   shared.fert.clutch$Site_dad != 
                                   shared.fert.clutch$Site_mom)] <- "ep_diff"
shared.fert.clutch$fert_type[which(is.na(shared.fert.clutch$FamilyID_dad))] <- "ep_unsampled"

# order factors
shared.fert.clutch$fert_type <- factor(shared.fert.clutch$fert_type, 
                                     levels=c("wp","ep_same","ep_unsampled",
                                              "ep_diff"))

# update shared.fert.clutc to re-classify ferts at the same site but 
# between different out buildings
# Blue Cloud-19 to Blue Cloud-14
# Cooks-23 to Hepp

shared.fert.clutch.update <- shared.fert.clutch

shared.fert.clutch.update$fert_type[which(
  shared.fert.clutch.update$clutch_id_ind2=="BlueCloud-14_collected" &
    shared.fert.clutch.update$FamilyID_dad=="BlueCloud-19")] <- "ep_diff"

shared.fert.clutch.update$fert_type[which(
  shared.fert.clutch.update$clutch_id_ind2=="Hepp (near Cook)_1_1" &
    shared.fert.clutch.update$FamilyID_dad=="Cooks-23")] <- "ep_diff"

# add new category for sites where we banded all known breeding males, and
# unsampled sires are therefore from outside the site: 
# Cooks, Make Believe, Struthers
# Only Cooks has fert originally classified as ep_unsampled
shared.fert.clutch.update$fert_type <- as.character(
  shared.fert.clutch.update$fert_type)

shared.fert.clutch.update$fert_type[which(
  shared.fert.clutch.update$Site_mom=="Cooks" &
    shared.fert.clutch.update$fert_type=="ep_unsampled")] <- "ep_unsamp_diff"


## save ferts by clutch
write.csv(shared.fert.clutch.update, "02_output-files/fert_2021_by_clutchID.csv",
          row.names = F)

# summarise number and proportion of each offspring type
total.fert.clutch.update <- sum(shared.fert.clutch.update$fert)

# summary of fert types classified at the clutch level
fert.type.clutch <- shared.fert.clutch.update %>%
  group_by(fert_type) %>%
  summarise(num.fert = sum(fert),
            prop.fert = num.fert/total.fert)

# by year
shared.fert.clutch.update$year <- 2021

ggplot(shared.fert.clutch.update, aes(x=as.factor(year),
                             y=fert, fill=fert_type)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d() + ggtitle("Number of fertilization types") +
  xlab("Year") + ylab("Number of\n fertilizations")



