
################################################################################
# script for summarizing total RS and multiple mating by males and females
# Identifies genetic families and full sib relationships
# Makes plots of fertilization types based on sites size, brood, etc
# 
# Heather Kenny-Duddela
# July 14, 2023

# libraries
library(tidyverse)
library(ggplot2)

# load data - table of assigned parents
# from object in paternity 2022 assignment with lcMLkin called kin.cutoff.info
load("input-files/kin_cutoff_info.Rdata")




#-------------------------------------------------------------------------------
# Assess variation in relatedness value for full sibs
# Identify full sibs by having same mom and dad assigned
#-------------------------------------------------------------------------------


# switch to wide format where each row is a kid, and there are cols for dad and mom
# first just do dads
po.dad.wide <- subset(kin.cutoff.info, kin.cutoff.info$Type_ind1=="dad")
# pull out only cols that need to get "dad" labels added
po.dad.wide2 <- po.dad.wide[,c(1,3:11)]
colnames(po.dad.wide2)[2:10] <- c("k0_hat_dad", "k1_hat_dad", "k2_hat_dad", "pi_HAT_dad",
                                  "nbSNP_dad", "Band_dad","Site_dad", "Type_dad" , "FamilyID_dad")

# do same for moms
po.mom.wide <- subset(kin.cutoff.info, kin.cutoff.info$Type_ind1=="mom")
po.mom.wide2 <- po.mom.wide[,c(1,3:11)]
colnames(po.mom.wide2)[2:10] <- c("k0_hat_mom", "k1_hat_mom", "k2_hat_mom", "pi_HAT_mom",
                                  "nbSNP_mom", "Band_mom","Site_mom", "Type_mom" , "FamilyID_mom")

# pull out columns with relevant kid info, and no duplicated kids
po.kids <- kin.cutoff.info[-which(duplicated(kin.cutoff.info$Ind2)),c(1,15:21)]

# add dad columns to kid table
po.wide.d <- left_join(po.kids, po.dad.wide2, by="Ind2")

# add mom columns to kid table
po.wide <- left_join(po.wide.d, po.mom.wide2, by="Ind2")

# create column to specify genetic family, kids with same genetic mom and dad
po.wide$genetic_fam <- paste(po.wide$Band_dad, po.wide$Band_mom, sep="_")

# move genetic family label near the front of the table
po.wide <- po.wide[,c(1, 27, 2:26)]

# remove duplicate rows due to some males having multiple family IDs
# identify duplicated rows across kid ID, kid family, and dad band
po.wide.dup <- po.wide[duplicated(po.wide[,c(1,6,15)]),]
# get all rows with duplicated kids, not just ones identified by duplicated
po.wide.dup2 <- po.wide[(po.wide$Ind2 %in% po.wide.dup$Ind2), ]
# make subset of kid table with duplicated kids and dads removed
po.wide.sub <- po.wide[-(which(po.wide$Ind2 %in% po.wide.dup$Ind2)), ]
# Reduce the duplicated kids table so that each kid-dad pair has only a single row
# For most cases there is a FamilyID_dad that matches the kid FamilyID, so preferentially take
# these matching FamilyIDs instead of non-matching ones. These are true within-pair offspring. 
# For CO_83279 there is not a matching FamilyID_dad, so this is an EP kid from a dad that was
# attending two other nests. Here just take the first occurrence, which will not be marked as 
# T from duplicated()
po.wide.dup3 <- subset(po.wide.dup2, po.wide.dup2$Ind2=="CO_83279" &
                         !duplicated(po.wide.dup2[,c(1,6,15)]) # keep first row of CO_83279
                       |
                         po.wide.dup2$Ind2!="CO_83279") # Keep all other kids
# Keep other kid rows where kid FamilyID matches dad FamilyID, and also keep the one CO_83279 row
po.wide.dup4 <- po.wide.dup3[which(po.wide.dup3$FamilyID_ind2 == po.wide.dup3$FamilyID_dad |
                                     po.wide.dup3$Ind2=="CO_83279"), ]
# Remove duplicate rows where kid FamilyID matches dad FamilyID (because these rows are true duplicates)
po.wide.dup5 <- po.wide.dup4[-which(duplicated(po.wide.dup4[,c(1,6,15)])),]

# add keeper rows back to the main dataframe
po.wide.c <- rbind(po.wide.sub, po.wide.dup5)


################################################################################

# add column for categorical site size
po.wide.c$site_type <- NA

po.wide.c$site_type[which(po.wide.c$Site_ind2 == "Cathy's" |
                            po.wide.c$Site_ind2 == "Karen's" |
                            po.wide.c$Site_ind2 == "Dome House" |
                            po.wide.c$Site_ind2 == "Marte's" |
                            po.wide.c$Site_ind2 == "Mary Ann's" |
                            po.wide.c$Site_ind2 == "Speedwell")] <- "solitary"

po.wide.c$site_type[which(po.wide.c$Site_ind2 == "Urban Farm Girlz" |
                            po.wide.c$Site_ind2 == "McCauley" |
                            po.wide.c$Site_ind2 == "Struthers")] <- "small"

po.wide.c$site_type[which(po.wide.c$Site_ind2 == "Boyer" |
                            po.wide.c$Site_ind2 == "Blue Cloud" |
                            po.wide.c$Site_ind2 == "Make Believe" |
                            po.wide.c$Site_ind2 == "Cooks")] <- "medium"          

po.wide.c$site_type[which(po.wide.c$Site_ind2 == "CHR")] <- "large"

# reorder factors
po.wide.c$site_type <- factor(po.wide.c$site_type,
                              levels=c("solitary", "small","medium","large"))

# rename brood column
colnames(po.wide.c)[8] <- "brood"

# save wide file
write.csv(po.wide.c, file="generated-files/kin_2022_parent_offspring_assigned_wide.csv")


################################################################################


# count up number of offspring produced by each genetic family
shared.fert.list <- po.wide.c %>% 
  group_by(genetic_fam, Band_dad, Band_mom, 
           Site_dad, Site_mom, FamilyID_dad, FamilyID_mom, site_type) %>%
  summarise(fert = n())

hist(shared.fert.list$fert)

# add labels to specify fertilization type: within-pair (wp), extra-pair same site (ep_same), 
# extra-pair different site (ep_diff), unknown dad (dad_unk), unknown mom (mom_unk)
# based on matching mom and dad FamilyIDs and Site 

shared.fert.list$fert_type <- NA
shared.fert.list$fert_type[which(shared.fert.list$FamilyID_dad == 
                                   shared.fert.list$FamilyID_mom)] <- "wp"
shared.fert.list$fert_type[which(shared.fert.list$FamilyID_dad != shared.fert.list$FamilyID_mom &
                                   shared.fert.list$Site_dad == shared.fert.list$Site_mom)] <- "ep_same"
shared.fert.list$fert_type[which(shared.fert.list$FamilyID_dad != shared.fert.list$FamilyID_mom &
                                   shared.fert.list$Site_dad != shared.fert.list$Site_mom)] <- "ep_diff"
shared.fert.list$fert_type[which(is.na(shared.fert.list$FamilyID_mom))] <- "mom_unk"
shared.fert.list$fert_type[which(is.na(shared.fert.list$FamilyID_dad))] <- "dad_unk"


# save shared fertilizations list
write.csv(shared.fert.list, file="generated-files/fert_2022_by_genetic_family.csv", row.names=F)



# summarise number and proportion of each offspring type
total.fert <- sum(shared.fert.list$fert)

fert.type.summary <- shared.fert.list %>%
  group_by(fert_type) %>%
  summarise(num.fert = sum(fert),
            prop.fert = num.fert/total.fert)

# summarize fertilization types within mom sites
fert.site.sum <- shared.fert.list %>%
  group_by(Site_mom, fert_type, site_type) %>%
  summarise(num.fert = sum(fert))

# reorder site factors to be in increasing site size
fert.site.sum$Site_mom <- factor(fert.site.sum$Site_mom,
                                 levels=c("Cathy's", "Karen's", "Dome House", "Marte's",
                                          "Mary Ann's", "Speedwell", "Urban Farm Girlz",
                                          "McCauley","Struthers","Boyer","Blue Cloud",
                                          "Make Believe", "Cooks", "CHR"))

fert.site.sum$fert_type <- factor(fert.site.sum$fert_type, levels=c("wp","ep_same","mom_unk",
                                                                    "dad_unk","ep_diff"))

ggplot(subset(fert.site.sum, !is.na(fert.site.sum$Site_mom)), 
       aes(x=Site_mom, y=num.fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() + ggtitle("Number of fertilization types by site") +
  xlab("Nest site (mom, arranged by size)") + ylab("Number of\n fertilizations")

ggsave("generated-files/fert type number by site 2022.png", width=5, height=3)


ggplot(subset(fert.site.sum, !is.na(fert.site.sum$Site_mom)), 
       aes(x=Site_mom, y=num.fert, fill=fert_type)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d()+ ggtitle("Proportion of fertilization types by site") +
  xlab("Nest site (mom, arranged by size)") + ylab("Proportion of\n fertilizations")

# ggsave("fert type proportion by site 2022.png", width=5, height=3)


#-------------------------------------------------------------------------------
# summarise by site size
#-------------------------------------------------------------------------------

ggplot(subset(fert.site.sum, !is.na(fert.site.sum$Site_mom)), 
       aes(x=site_type, y=num.fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() + ggtitle("Number of fertilization types by site size") +
  xlab("Site size type") + ylab("Number of fertilizations")

ggsave("generated-files/fert number by site size.png", height=3, width=5)

ggplot(subset(fert.site.sum, !is.na(fert.site.sum$Site_mom)), 
       aes(x=site_type, y=num.fert, fill=fert_type)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() + ggtitle("Number of fertilization types by site size") +
  xlab("Site size type") + ylab("Proportion of fertilizations")

# ggsave("fert proportion by site size.png", height=3, width=5)




#-------------------------------------------------------------------------------
# at the nest/clutch level, look at proportion of outside site dads
# get a sense for whether patterns are usually lots of EP at few nests, or 
# some EP at most nests
#-------------------------------------------------------------------------------


# summarise at clutch level
shared.fert.clutch <- po.wide.c %>% 
  group_by(clutch_id_ind2, Band_dad, Band_mom, 
           Site_dad, Site_mom, FamilyID_dad, FamilyID_mom, brood, site_type) %>%
  summarise(fert = n())

# number of clutches
length(unique(shared.fert.clutch$clutch_id_ind2))

# assign fert categories
shared.fert.clutch$fert_type <- NA
shared.fert.clutch$fert_type[which(shared.fert.clutch$FamilyID_dad == 
                                     shared.fert.clutch$FamilyID_mom)] <- "wp"
shared.fert.clutch$fert_type[which(shared.fert.clutch$FamilyID_dad != shared.fert.clutch$FamilyID_mom &
                                     shared.fert.clutch$Site_dad == shared.fert.clutch$Site_mom)] <- "ep_same"
shared.fert.clutch$fert_type[which(shared.fert.clutch$FamilyID_dad != shared.fert.clutch$FamilyID_mom &
                                     shared.fert.clutch$Site_dad != shared.fert.clutch$Site_mom)] <- "ep_diff"
shared.fert.clutch$fert_type[which(is.na(shared.fert.clutch$FamilyID_mom))] <- "mom_unk"
shared.fert.clutch$fert_type[which(is.na(shared.fert.clutch$FamilyID_dad))] <- "dad_unk"


# summarise number and proportion of each offspring type
total.fert.clutch <- sum(shared.fert.clutch$fert)

# number of females where we have both 1st and 2nd broods
fem.both.broods <- shared.fert.clutch %>%
  group_by(FamilyID_mom, brood) %>%
  summarise(chicks = n())

fem.both.brood2 <- subset(fem.both.broods, fem.both.broods$brood == "1" |
                            fem.both.broods$brood == "2")

fem.both.brood3 <- fem.both.brood2 %>%
  group_by(FamilyID_mom) %>%
  summarise(brood_count = n())

# 33 females had multiple broods
length(subset(fem.both.brood3$FamilyID_mom, fem.both.brood3$brood_count==2))  

# proportion of within-pair offspring per clutch
shared.fert.clutch2 <- shared.fert.clutch %>% 
  group_by(clutch_id_ind2) %>%
  mutate(clutch_size = sum(fert))

shared.fert.wp <- subset(shared.fert.clutch2, shared.fert.clutch2$fert_type=="wp")
shared.fert.wp <- shared.fert.wp %>%
  mutate(prop.wp = fert/clutch_size)

shared.fert.wp2 <- left_join(shared.fert.clutch2, shared.fert.wp[,c(1,13)], by="clutch_id_ind2")
shared.fert.wp2$prop.wp[which(is.na(shared.fert.wp2$prop.wp))] <- 0

ggplot(shared.fert.wp2, aes(x=prop.wp)) + 
  geom_histogram(fill="#440154FF", color="black", bins = 16) +
  xlab("Proportion of within-pair offspring") + ylab("Count") + 
  ggtitle("Distirbution of within-pair proportions per clutch")

ggsave("generated-files/prop within pair by clutch.png", height=4, width=6)


# plot across all nests

# change factor order to match previous plots
shared.fert.clutch$fert_type <- factor(shared.fert.clutch$fert_type, levels=c("wp","ep_same","mom_unk",
                                                                              "dad_unk","ep_diff"))

# reorder site factors to be in increasing site size
shared.fert.clutch$Site_mom <- factor(shared.fert.clutch$Site_mom,
                                      levels=c("Cathy's", "Karen's", "Dome House", "Marte's",
                                               "Mary Ann's", "Speedwell", "Urban Farm Girlz",
                                               "McCauley","Struthers","Boyer","Blue Cloud",
                                               "Make Believe", "Cooks", "CHR"))

ggplot(shared.fert.clutch, aes(x=clutch_id_ind2, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Clutch ID") + ylab("Number of\n fertilizations")

ggplot(shared.fert.clutch, aes(x=clutch_id_ind2, y=fert, fill=fert_type)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d()+ ggtitle("Proportion of fertilization types by clutch") +
  xlab("Clutch ID") + ylab("Proportion of\n fertilizations")

# facet by brood 1 and 2
ggplot(shared.fert.clutch, aes(x=clutch_id_ind2, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Clutch ID") + ylab("Number of\n fertilizations") + facet_wrap("brood")

# brood 1 and 2 by site
ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1" &
                !is.na(shared.fert.clutch$Site_mom)|
                shared.fert.clutch$brood=="2" & 
                !is.na(shared.fert.clutch$Site_mom)), aes(x=Site_mom, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Site") + ylab("Number of\n fertilizations") + facet_grid(brood~.) +
  ggtitle("Fertilizations types by site for brood 1 and 2")

ggsave("generated-files/fert by site for broods 1 and 2.png", width=5, height=4)

ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1"|
                shared.fert.clutch$brood=="2"), aes(x=Site_mom, y=fert, fill=fert_type)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Site") + ylab("Number of\n fertilizations") + facet_grid(brood~.) +
  ggtitle("Fertilizations types by site for brood 1 and 2")

# brood 1 and 2 by mom ID
ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1" &
                !is.na(shared.fert.clutch$Site_mom)|
                shared.fert.clutch$brood=="2" & 
                !is.na(shared.fert.clutch$Site_mom)), 
       aes(x=Band_mom, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Mom Band") + ylab("Number of\n fertilizations") + facet_grid(brood~.) +
  ggtitle("Fertilizations types by female for brood 1 and 2")

ggsave("generated-files/fert by female for brood 1 and 2.png", width=7, height=4)

ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1"|
                shared.fert.clutch$brood=="2"), aes(x=Band_mom, y=fert, fill=fert_type)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Mom Band") + ylab("Number of\n fertilizations") + facet_grid(brood~.) +
  ggtitle("Fertilizations types by female for brood 1 and 2")

# brood 1 and 2 by dad ID
ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1" &
                !is.na(shared.fert.clutch$Site_mom)|
                shared.fert.clutch$brood=="2" & 
                !is.na(shared.fert.clutch$Site_mom)), 
       aes(x=Band_dad, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Dad Band") + ylab("Number of\n fertilizations") + facet_grid(brood~.) +
  ggtitle("Fertilizations types by male for brood 1 and 2")

# ggsave("fert by male for brood 1 and 2.png", width=7, height=4)

# remove unknown dads
ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1" &
                !is.na(shared.fert.clutch$Site_mom) &
                !is.na(shared.fert.clutch$Band_dad)
              |
                shared.fert.clutch$brood=="2" & 
                !is.na(shared.fert.clutch$Site_mom) &
                !is.na(shared.fert.clutch$Band_dad)), 
       aes(x=Band_dad, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Dad Band") + ylab("Number of\n fertilizations") + facet_grid(brood~.) +
  ggtitle("Fertilizations types by male for brood 1 and 2")

ggsave("generated-files/fert by male for brood 1 and 2 no NA.png", width=7, height=4)


# summarize based on brood 1 and 2
ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1" &
                !is.na(shared.fert.clutch$Site_mom)|
                shared.fert.clutch$brood=="2" & 
                !is.na(shared.fert.clutch$Site_mom)), 
       aes(x=brood, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("brood") + ylab("Number of fertilizations") +
  ggtitle("Fertilizations types for brood 1 and 2")

# ggsave("fert types by brood.png", height=5, width=5)

ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1" &
                !is.na(shared.fert.clutch$Site_mom)|
                shared.fert.clutch$brood=="2" & 
                !is.na(shared.fert.clutch$Site_mom)), 
       aes(x=brood, y=fert, fill=fert_type)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("brood") + ylab("Proportion of fertilizations") +
  ggtitle("Fertilizations types for brood 1 and 2")

ggsave("generated-files/fert proportion by brood.png", height=5, width=5)


# by brood 1 and 2 across site sizes
ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1" &
                !is.na(shared.fert.clutch$Site_mom)|
                shared.fert.clutch$brood=="2" & 
                !is.na(shared.fert.clutch$Site_mom)), 
       aes(x=site_type, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() + facet_grid(brood~.) +
  xlab("Site size type") + ylab("Number of fertilizations") +
  ggtitle("Fertilizations types for brood 1 and 2 by site type")

ggsave("generated-files/fert number by brood and site type.png", width=6, height=5)

ggplot(subset(shared.fert.clutch, shared.fert.clutch$brood=="1" &
                !is.na(shared.fert.clutch$Site_mom)|
                shared.fert.clutch$brood=="2" & 
                !is.na(shared.fert.clutch$Site_mom)), 
       aes(x=site_type, y=fert, fill=fert_type)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() + facet_grid(brood~.) +
  xlab("Site size type") + ylab("Proportion of fertilizations") +
  ggtitle("Fertilizations types for brood 1 and 2 by site type")

# ggsave("fert prop by brood and site type.png", width=6, height=5)

#-------------------------------------------------------------------------------
# pull out full-sib relatedness values based on genetic fams

# get list of unique genetic families with mom and dad IDs known

known.gen.fam <- subset(po.wide, !is.na(po.wide$FamilyID_dad) & !is.na(po.wide$FamilyID_mom))
num.gen.fam <- unique(known.gen.fam$genetic_fam)

# load full kinship table (not just parent-offspring relationships)
load("input-files/kin2_22.Rdata")

# loop through the different families

# storage for full genetic sibs
sib.storage <- kin2.22[1,]
sib.storage[,] <- NA
sib.storage$genetic_fam <- NA

for (i in 1:length(num.gen.fam)) {
  fam <- subset(known.gen.fam, known.gen.fam$genetic_fam == num.gen.fam[i])
  if (length(unique(fam$Ind2)) > 1) {
    gen.sib <- subset(kin2.22, kin2.22$Ind1 %in% fam$Ind2 &
                        kin2.22$Ind2 %in% fam$Ind2)
    gen.sib$genetic_fam <- num.gen.fam[i]
    sib.storage <- rbind(sib.storage, gen.sib)}
}

# remove row of NAs
sib.storage <- subset(sib.storage, !is.na(sib.storage$Ind1))

# plot relatedness values for full sibs
plot(sib.storage$k0_hat, sib.storage$pi_HAT)

ggplot(sib.storage, aes(k0_hat, pi_HAT)) + geom_point(shape=1) +
  ylab("kinship value (min=0.244, max=0.572)") +
  xlab("Est probability of unrelated (min=0.206, max=0.637)") +
  ggtitle("Range of kinship and k0-hat values for full sibs")

ggsave("generated-files/full sib kinship values.png")

# max and min cutoff values
range(sib.storage$k0_hat) # k0 from 0.206 to 0.637
range(sib.storage$pi_HAT) # pi_hat from 0.244 to 0.572

full.sib.pi.cutoff <- min(sib.storage$pi_HAT)

#-------------------------------------------------------------------------------
# add site sizes to shared.fert.clutch2
#-------------------------------------------------------------------------------

# Blue Cloud - 7
# Cooks - 9
# CHR - 32
# Make Believe - 7
# Struthers - 3
# McCauley - 2
# Urban Farm Girlz - 2
# Dome House - 1
# Marte's - 1
# Boyer - 5
# Cathy's - 1
# Mary Ann's - 1

# for dads
shared.fert.clutch2$site_size_dad <- NA
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Blue Cloud")] <- 7
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Starlight")] <- 3
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Jay's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Boyer")] <- 5
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Cathy's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="CHR")] <- 32
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Mayas")] <- 9
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Cooks")] <- 9
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Reinarz")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Dome House")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Karen's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Make Believe")] <- 7
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Marte's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Mary Ann's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="McCauley")] <- 2
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Speedwell")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Struthers")] <- 3
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Urban Farm Girlz")] <- 2

# for moms
shared.fert.clutch2$site_size_mom <- NA
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Blue Cloud")] <- 7
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Starlight")] <- 3
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Jay's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Boyer")] <- 5
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Cathy's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="CHR")] <- 32
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Mayas")] <- 9
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Cooks")] <- 9
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Reinarz")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Dome House")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Karen's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Make Believe")] <- 7
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Marte's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Mary Ann's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="McCauley")] <- 2
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Speedwell")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Struthers")] <- 3
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Urban Farm Girlz")] <- 2

#-------------------------------------------------------------------------------
# Calculate total reproductive success for each male
# first subset males where we monitored their social nest
#-------------------------------------------------------------------------------

# get males where social nest was monitored
# this should be cases where familyID is not "none"

load("input-files/fam_clutch.RData") # generated by script "add clutch id to family info"
fam.22 <- fam.clutch

males.sampled <- subset(fam.22, fam.22$Type=="dad")

# where male familyID is represented in kid FamilyIDs
males.monitored <- subset(males.sampled, males.sampled$FamilyID %in%
                            kin.cutoff.info$FamilyID_ind2)



# pull out monitored males from shared fert list (59 male nests monitored)
shared.fert.monitored <- subset(shared.fert.clutch2, shared.fert.clutch2$Band_dad %in%
                                  males.monitored$Band)

# summarise across monitored males
fert.male.monitored <- shared.fert.monitored %>%
  group_by(Band_dad, Site_dad, site_size_dad) %>%
  summarise(num.chicks = sum(fert))

ggplot(fert.male.monitored, aes(x=reorder(Band_dad, num.chicks), y=num.chicks)) + 
  geom_bar(fill="#440154FF", color="black", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("Total reproductive success for each male") +
  xlab("Male band number") + ylab("Total offspring")

ggsave("generated-files/Male total RS 2022.png", h=4, w=6)

# color by site
ggplot(fert.male.monitored, aes(x=reorder(Band_dad, num.chicks), y=num.chicks)) + 
  geom_bar(aes(fill=Site_dad), color="black", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("Total reproductive success for each male") +
  xlab("Male band number")

# facet by site
ggplot(fert.male.monitored, aes(x=reorder(Band_dad, num.chicks), y=num.chicks)) + 
  geom_bar(fill="#440154FF", color="black", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("Total reproductive success for each male") +
  xlab("Male band number") + facet_wrap("Site_dad")

# check correlation with site size
ggplot(fert.male.monitored, aes(x=site_size_dad, y=num.chicks)) + geom_point(alpha=0.5) + 
  geom_smooth(method=lm, se=F) + xlab("Number of breeding pairs at the male's nesting site") +
  ylab("Total offspring") + 
  ggtitle("Relationship between site size and male total RS\nSpearman Rho=0.24, p=0.079")

ggsave("generated-files/site size and male total RS.png")

cor.test(fert.male.monitored$site_size_dad, fert.male.monitored$num.chicks, method="spearman")


## summarise males by fertilization type
fert.male.monitored.type <- shared.fert.monitored %>%
  group_by(Band_dad, FamilyID_dad, Site_dad, fert_type, site_size_dad) %>%
  summarise(num.chicks = sum(fert))

fert.male.monitored.type <- fert.male.monitored.type %>%
  group_by(Band_dad) %>%
  mutate(tot.chicks = sum(num.chicks))

# add column for just ep and wp
fert.male.monitored.type$wp_ep <- NA
fert.male.monitored.type$wp_ep[which(fert.male.monitored.type$fert_type=="wp")] <- "wp" 
fert.male.monitored.type$wp_ep[which(fert.male.monitored.type$fert_type!="wp")] <- "ep" 

# change factor order to match previous plots
fert.male.monitored.type$fert_type <- factor(fert.male.monitored.type$fert_type, 
                                             levels=c("dad_unk", "mom_unk","ep_diff","ep_same", "wp"))

# reorder site factors to be in increasing site size
fert.male.monitored.type$Site_dad <- factor(fert.male.monitored.type$Site_dad,
                                            levels=c("Cathy's", "Karen's", "Dome House", "Marte's",
                                                     "Mary Ann's", "Speedwell", "Urban Farm Girlz",
                                                     "McCauley","Struthers","Boyer","Blue Cloud",
                                                     "Make Believe", "Cooks", "CHR"))

# reorder site factors to be in increasing site size
fert.male.monitored$Site_dad <- factor(fert.male.monitored$Site_dad,
                                       levels=c("Cathy's", "Karen's", "Dome House", "Marte's",
                                                "Mary Ann's", "Speedwell", "Urban Farm Girlz",
                                                "McCauley","Struthers","Boyer","Blue Cloud",
                                                "Make Believe", "Cooks", "CHR"))


# color by fert type
ggplot(fert.male.monitored.type, aes(x=reorder(Band_dad, tot.chicks), 
                                     y=num.chicks, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values = c("#39568CFF","#FDE725FF", "#55C667FF", "#440154FF"))+
  ggtitle("Total reproductive success for each male by type") +
  xlab("Male band number") + ylab("Total offspring")

ggsave("generated-files/Male total RS by fertilization type 2022.png", h=4, w=7)

# color by ep and wp
ggplot(fert.male.monitored.type, aes(x=reorder(Band_dad, tot.chicks), 
                                     y=num.chicks, fill=wp_ep)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values = c("#55C667FF", "#440154FF"))+
  ggtitle("Total reproductive success for each male by ep and wp") +
  xlab("Male band number") + ylab("Total offspring")

ggsave("generated-files/Male total RS by wp and ep 2022.png", h=4, w=7)

# arrange by site
fert.male.monitored %>% arrange(across(.cols=c("Site_dad","num.chicks"))) %>%
  rowid_to_column %>%
  ggplot() +
  geom_bar(aes(x=reorder(Band_dad, rowid), y=num.chicks, fill=Site_dad), stat="identity") +
  theme(axis.text.x = element_text(angle=90))

# arrange by site and color by fer type
fert.male.monitored.type %>% arrange(across(.cols=c("Site_dad","tot.chicks"))) %>%
  rowid_to_column %>%
  ggplot() +
  geom_bar(aes(x=reorder(Band_dad, rowid), y=num.chicks, fill=fert_type), 
           position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values = c("#39568CFF","#FDE725FF", "#55C667FF", "#440154FF"))

ggplot(fert.male.monitored.type, aes(x=reorder(Band_dad, Site_dad), 
                                     y=num.chicks, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values = c("#39568CFF","#FDE725FF", "#55C667FF", "#440154FF"))+
  ggtitle("Total reproductive success for each male by type") +
  xlab("Male band number") + ylab("Total offspring")


# check relationship between site size and ep wp for males

# ep only
ggplot(fert.male.monitored, aes(x=site_size_dad, y=num.chicks)) + geom_point(alpha=0.5) +
  geom_smooth(data=subset(fert.male.monitored.type, fert.male.monitored.type$fert_type=="wp"), 
              aes(x=site_size_dad, y=num.chicks), method=lm, color="blue", se=F) +
  geom_smooth(data=subset(fert.male.monitored.type, fert.male.monitored.type$fert_type=="ep_diff" |
                            fert.male.monitored.type$fert_type=="ep_same"), 
              aes(x=site_size_dad, y=num.chicks), method=lm, color="red", se=F) +
  ggtitle("Relationship between site size and male RS\n (EP=red, WP=blue)") +
  ylab("Total offspring") + xlab("Number of breeding pairs at male's nesting site")

ggsave("generated-files/site size and male RS by EP and WP.png")



##################################################################################

### summarise by female

fert.female <- shared.fert.clutch2 %>%
  group_by(Band_mom, FamilyID_mom, Site_mom, site_size_mom) %>%
  summarise(num.chicks = sum(fert))


ggplot(subset(fert.female, !is.na(fert.female$Band_mom)), 
       aes(x=reorder(Band_mom, num.chicks), y=num.chicks)) + 
  geom_bar(fill="#440154FF", color="black", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("Total reproductive success for each female") +
  xlab("Female band number") + ylab("Total offspring") + ylim(0,17)

ggsave("generated-files/female total RS 2022.png", h=4, w=6)

#  include fertilization type
fert.female.type <- shared.fert.clutch2 %>%
  group_by(Band_mom, FamilyID_mom, Site_mom, fert_type, site_size_mom) %>%
  summarise(num.chicks = sum(fert))

fert.female.type <- fert.female.type %>%
  group_by(Band_mom) %>%
  mutate(tot.chicks = sum(num.chicks))

# add column for just wp or ep
fert.female.type$wp_ep <- NA
fert.female.type$wp_ep[which(fert.female.type$fert_type=="wp")] <- "wp"
fert.female.type$wp_ep[which(fert.female.type$fert_type!="wp")] <- "ep"

# change factor order to match previous plots
fert.female.type$fert_type <- factor(fert.female.type$fert_type, 
                                     levels=c("dad_unk", "mom_unk","ep_diff","ep_same", "wp"))

ggplot(subset(fert.female.type, !is.na(fert.female.type$Band_mom)), 
       aes(x=reorder(Band_mom, tot.chicks), y=num.chicks, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values = c("#39568CFF","#FDE725FF","#55C667FF","#440154FF"))+
  ggtitle("Total reproductive success for each female by type") +
  ylab("Total offspring") + xlab("Female band number") + ylim(0,17)

ggsave("generated-files/Female total RS by fertilization type 2022.png", h=4, w=7)

# plot with just two colors for ep and wp
ggplot(subset(fert.female.type, !is.na(fert.female.type$Band_mom)), 
       aes(x=reorder(Band_mom, tot.chicks), y=num.chicks, fill=wp_ep)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values = c("#55C667FF","#440154FF"))+
  ggtitle("Total reproductive success for each female by wp and ep") +
  ylab("Total offspring") + xlab("Female band number") + ylim(0,17)

ggsave("generated-files/Female total RS by wp and ep 2022.png", h=4, w=7)

## check correlation between female RS and site size ##

ggplot(fert.female, aes(x=site_size_mom, y=num.chicks)) + geom_point(alpha=0.5) + 
  geom_smooth(method=lm, se=F) + xlab("Number of breeding pairs at the female's nesting site") +
  ylab("Total offspring") + 
  ggtitle("Relationship between site size and female total RS\nSpearman Rho=0.288, p=0.029")

ggsave("generated-files/site size and female total RS.png")

cor.test(fert.female$site_size_mom, fert.female$num.chicks, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fert.female$site_size_mom and fert.female$num.chicks
# S = 21944, p-value = 0.02934
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2888229 

# trends for within pair and extra-pair

# both ep and wp for females
ggplot(fert.female, aes(x=site_size_mom, y=num.chicks)) + geom_point(alpha=0.5) +
  geom_smooth(data=subset(fert.female.type, fert.female.type$fert_type=="wp"), 
              aes(x=site_size_mom, y=num.chicks), method=lm, color="blue", se=F) +
  geom_smooth(data=subset(fert.female.type, fert.female.type$fert_type!="wp"), 
              aes(x=site_size_mom, y=num.chicks), method=lm, color="red", se=F) +
  ggtitle("Relationship between site size and female RS\n (EP=red, WP=blue)") +
  ylab("Total offspring") + xlab("Number of breeding pairs at female's nesting site")

ggsave("generated-files/site size and female RS by EP and WP.png")

# just wp
ggplot(subset(fert.female.type, fert.female.type$fert_type=="wp"), 
       aes(x=site_size_mom, y=num.chicks)) +
  geom_point(alpha=0.5) + geom_smooth(method=lm) +
  ggtitle("Relationship between site size and female WP offspring") +
  ylab("Number of within-pair offspring") +
  xlab("Number of breeding pairs at the female's nesting site")

# just non-wp
ggplot(subset(fert.female.type, fert.female.type$fert_type!="wp"), 
       aes(x=site_size_mom, y=num.chicks)) +
  geom_point(alpha=0.5) + geom_smooth(method=lm) +
  ggtitle("Relationship between site size and female EP offspring") +
  ylab("Number of extra-pair offspring") +
  xlab("Number of breeding pairs at the female's nesting site")


############# calculate proportion of males and females that mate multiply ############

# proportion of males that mate multiply 
fert.male.monitored.type$prop.wp <- fert.male.monitored.type$num.chicks/fert.male.monitored.type$tot.chicks

fert.male.wp <- subset(fert.male.monitored.type, fert.male.monitored.type$wp_ep=="wp")

# 33 out of 54 males with only wp offspring
fert.male.wp.only <- subset(fert.male.wp, fert.male.wp$prop.wp==1)

# 0.611 proportion of males with only WP.
prop.wp.male <- length(unique(fert.male.wp.only$Band_dad))/length(unique(fert.male.monitored.type$Band_dad))

# check whether any males only had EP and no WP
# 21 males had ep offspring
fert.male.ep <- subset(fert.male.monitored.type, fert.male.monitored.type$wp_ep=="ep")

# only 17 of these had both, so missing 4 males
both <- sum(fert.male.ep.only$Band_dad %in% fert.male.wp$Band_dad)

# try merging EP and WP
fert.male.both <- left_join(fert.male.ep, fert.male.wp, by="Band_dad")
# CHR-11 male (2850-57733) should have 50% WP but mom was not sampled for DNA
# Cooks-13 male (2850-57675) should have only WP but mom was not sampled for DNA
# CHR-045 male (2850-57653) had only one EP kid and zero WP

# add one to number of males with only WP = 33 +1 = 34
# proportion males that did not engage in EP is 34/54 = 0.629
# proportion of males that DID engage in EP is 20/54 = 0.37



### females

fert.female.type$prop.wp <- fert.female.type$num.chicks/fert.female.type$tot.chicks

fert.female.wp <- subset(fert.female.type, fert.female.type$fert_type=="wp")

fert.female.wp.only <- subset(fert.female.type, fert.female.type$prop.wp==1)

# proportion of females that had only WP = 0.413 (24/58)
# proportion of females that did engage in EP is 0.586 (34/58)
prop.wp.female <- length(unique(fert.female.wp.only$Band_mom))/length(unique(fert.female.type$Band_mom))

####################################################################################
# Save key data tables

# save fertilization types by clutch ID
write.csv(shared.fert.clutch2, file="generated-files/fert_2022_by_clutchID.csv", row.names=F)








