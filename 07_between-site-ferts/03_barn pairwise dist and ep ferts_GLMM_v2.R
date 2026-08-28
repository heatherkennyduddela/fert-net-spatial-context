
################################################################################
# Script to add between-site distance for EP ferts from all 3 years
# Heather Kenny-Duddela
# Aug 27, 2025
################################################################################

# load data

# generated in directory: 
# CU Boulder/BARS fieldwork/2022 Field and Lab/Nest maps and Sketch Up/pairwise dist from ArcGIS
dist <- read.csv("input-files/distances between Barns with ids.csv")

## 2021
# adults analyzed for paternity 2021 - from 2021 paternity with 2023 update directory
adult.kin <- read.csv("input-files/kin_2021_adults_all.csv")
# fertilizations types for 2021 by clutch
kin.assign <- read.csv("input-files/fert_2021_by_clutchID.csv")
# parent-offspring assignments 2021
po21 <- read.csv("input-files/kin_2021_parent_offspring_assigned_withUNS.csv")


## 2022
# fertilizations types for 2022 by clutch
assign22 <- read.csv("input-files/fert_2022_by_clutchID.csv")
# adults analyzed for paternity 2022 - from directory:
# 2022 paternity analysis v2/01_paternity 2022 assignment with lcMLkin output
adult22 <- read.csv("input-files/kin_2022_adults_all.csv")
# parent-offspring assignments 2022
po22 <- read.csv("input-files/kin_2022_parent_offspring_assigned_withUNS.csv")


## 2023
# fertilization types for 2023 by clutch - from directory: 
# 2023 Field and lab\Paternity analysis 2023\paternity 2023 with lcmlkin output\02_output-files
assign23 <- read.csv("input-files/fert_2023_by_clutchID.csv")
# adults analyzed for paternity 2023 - from directory:
# Paternity analysis 2023\paternity 2023 with lcmlkin output\01_output-files
adult23 <- read.csv("input-files/kin_2023_adults_all.csv")
# parent-offspring assignments 2023
po23 <- read.csv("input-files/kin_2023_parent_offspring_assigned_withUNS.csv")


# libraries
library(tidyverse)
library(ggplot2)
library(glmmTMB) # also for fitting mixed effect glm
library(DHARMa) # model diagnostics for mixed effect models
library(car)

#-------------------------------------------------------------------------------
# Add distances to between-site EP ferts 2021
#-------------------------------------------------------------------------------

# all ep
ep.all <- subset(kin.assign, kin.assign$fert_type=="ep_diff" |
                   kin.assign$fert_type=="ep_unsamp_diff" |
                   kin.assign$fert_type=="ep_same" |
                   kin.assign$fert_type=="ep_unsampled")

# change ep_unsamp_diff to ep_unsampled
ep.all.update <- ep.all
ep.all.update$fert_type[which(ep.all.update$fert_type=="ep_unsamp_diff")] <- "ep_unsampled"

# sum(ep.all.update$fert) # 70
# 
# table(ep.all.update$fert_type)
# 
# sum(subset(ep.all.update$fert, ep.all.update$fert_type=="ep_same" )) # 24
# sum(subset(ep.all.update$fert, ep.all.update$fert_type=="ep_unsampled" )) # 37
# sum(subset(ep.all.update$fert, ep.all.update$fert_type=="ep_diff" )) # 9
# 
# ## Chi-squared test
# 
# ep.all.obs <- c(rep("ep_same",24), rep("ep_unsampled", 37), rep("ep_diff", 9))
# ep.all.table <- table(ep.all.obs)
# 
# chisq.test(ep.all.table)
# Chi-squared test for given probabilities
# 
# data:  ep.all.table
# X-squared = 16.829, df = 2, p-value = 0.0002217

## add distances for known outside site sires-----------------------------------

# only known-distance sires
ep.diff <- subset(kin.assign, kin.assign$fert_type=="ep_diff")

# add building labels to match dist
ep.diff$building_dad <- NA
ep.diff$building_dad[which(ep.diff$Site_dad=="Cooks")] <- "Cooks"
ep.diff$building_dad[which(ep.diff$Site_dad=="Blue Cloud")] <- "Blue Cloud main"
ep.diff$building_dad[which(
  ep.diff$FamilyID_dad=="BlueCloud-19")] <- "Blue Cloud Middle Shed"

ep.diff$building_mom <- NA
ep.diff$building_mom[which(ep.diff$Site_mom=="Blue Cloud")] <- "Blue Cloud main"
ep.diff$building_mom[which(ep.diff$Site_mom=="Struthers")] <- "Struthers main"
ep.diff$building_mom[which(ep.diff$Site_mom=="Cooks")] <- "Hepp awning"
ep.diff$building_mom[which(ep.diff$Site_mom=="Make Believe")] <- "Make Believe"

# change colnames in dist
colnames(dist)[3:4] <- c("building_dad","building_mom")

# add distances
ep.diff.dist <- left_join(ep.diff, dist[,c(1,3,4)], 
                          by=c("building_dad","building_mom"))

#-------------------------------------------------------------------------------
# Determine number of available distances based on adults sampled 2021
#-------------------------------------------------------------------------------

# only keep female-male rows
adult.kin2 <- subset(adult.kin, adult.kin$Type_ind1=="mom" &
                       adult.kin$Type_ind2=="dad" |
                       adult.kin$Type_ind1=="dad" &
                       adult.kin$Type_ind2=="mom")

## make all dads ind1 and all moms ind2

adult.dad.mom <- subset(adult.kin, adult.kin$Type_ind1=="dad" &
                          adult.kin$Type_ind2=="mom")

adult.mom.dad <- subset(adult.kin, adult.kin$Type_ind1=="mom" &
                          adult.kin$Type_ind2=="dad")

# pull out ind1 info columns and change to ind2
adult.mom.ind1 <- adult.mom.dad[,1:11]
colnames(adult.mom.ind1)[8:11] <- gsub("ind1","ind2", colnames(adult.mom.ind1)[8:11])
colnames(adult.mom.ind1)[1:2] <- c("Ind2","Ind1")

# pull out ind2 cols and change to ind1
adult.dad.ind2 <- adult.mom.dad[,c(1:7,15:18)]
colnames(adult.dad.ind2)[8:11] <- gsub("ind2","ind1", colnames(adult.dad.ind2)[8:11])
colnames(adult.dad.ind2)[1:2] <- c("Ind2","Ind1")

# combine back together
adult.fixed <- full_join(adult.dad.ind2, adult.mom.ind1, 
                         by=c("Ind1","Ind2","k0_hat","k1_hat","k2_hat",
                              "pi_HAT","nbSNP"))

# add in rows that were correct to start with
adult.fixed2 <- rbind(adult.fixed, adult.dad.mom[,-c(12:14,19:21)])

## Match buildings with sites for males
# "Alishas" = "Alishas est"
# "Blue Cloud" = "Blue Cloud main" , "Blue Cloud wash room" , "Blue Cloud Middle Shed"
# "Cathys" = "Cathys est"
# "Cooks" = "Cooks" , "Hepp awning"
# "Jims" = ?? (missing)
# "Make Believe" = "Make Believe" , "Make Believe Shed"
# "Reinarz" = "Reinarz shed est"
# "Schaaps" = "Schaaps"
# "Struthers" = "Struthers main"
# "Urban Farm Girlz" = "Urban Farm Girlz Sheep Shed" , "Urban Farm Girlz Red Barn"

# Note that some males had more than one location through the season: 
# BlueCloud-19 male started in "Blue Cloud Middle Shed" and then moved to 
# "Blue Cloud wash room"
# Cooks-31 male was in "Cooks" and then moved to "Hepp awning"
# Cathys-3 male was in the trailer and then moved to Mary Anns

## First add simple building matches for dads-----------------------------------
adult.fixed2$building_dad <- NA
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Alishas")] <- "Alishas est"
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Blue Cloud")] <- "Blue Cloud main"
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Cathys")] <- "Cathys est"
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Cooks")] <- "Cooks"
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Make Believe")] <- "Make Believe"
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Reinarz")] <- "Reinarz shed est"
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Schaaps")] <- "Schaaps"
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Struthers")] <- "Struthers main"
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Urban Farm Girlz")] <- "Urban Farm Girlz Sheep Shed"
adult.fixed2$building_dad[which(
  adult.fixed2$Site_ind1=="Jims")] <- "Jims Green Barn est"

## add corrections for out buildings
# Make Believe 61 should be in shed, not main barn
# Urban Farm 2 and 5 in sheep shed, 4 and 7 in red barn
# Cathys-3 was in trailer and then at Mary Ann's
adult.fixed2$building_dad[which(
  adult.fixed2$FamilyID_ind1=="MakeBelieve-61")] <- "Make Believe Shed"
adult.fixed2$building_dad[which(
  adult.fixed2$FamilyID_ind1=="UrbanFarm-04")] <- "Urban Farm Girlz Red Barn"
adult.fixed2$building_dad[which(
  adult.fixed2$FamilyID_ind1=="Cathys-03")] <- "Mary Anns est"

# Blue Cloud 19 birds should occur twice, once in middle shed and once in wash room
# First change current set of dads to middle shed
adult.fixed2$building_dad[which(
  adult.fixed2$FamilyID_ind1=="BlueCloud-19")] <- "Blue Cloud Middle Shed"
# next duplicate these rows and change building to washroom
bc19.dad <- subset(adult.fixed2, adult.fixed2$FamilyID_ind1=="BlueCloud-19")
bc19.dad$building_dad <- "Blue Cloud wash room"

# Also need to add duplicate location for Cooks-31 male
co31.dad <- subset(adult.fixed2, adult.fixed2$FamilyID_ind1=="Cooks-31")
co31.dad$building_dad <- "Hepp awning"

# add duplicate location for Cathys-3 in the trailer est
ca03.dad <- subset(adult.fixed2, adult.fixed2$FamilyID_ind1=="Cathys-03")
ca03.dad$building_dad <- "Cathys trailer est"


# combine the duplicate locations with the main table
adult.fixed3 <- rbind(adult.fixed2, bc19.dad, co31.dad, ca03.dad)

## Same process for adding mom buildings----------------------------------------
adult.fixed3$building_mom <- NA
adult.fixed3$building_mom[which(
  adult.fixed3$Site_ind2=="Alishas")] <- "Alishas est"
adult.fixed3$building_mom[which(
  adult.fixed3$Site_ind2=="Blue Cloud")] <- "Blue Cloud main"
adult.fixed3$building_mom[which(
  adult.fixed3$Site_ind2=="Cathys")] <- "Cathys est"
adult.fixed3$building_mom[which(
  adult.fixed3$Site_ind2=="Cooks")] <- "Cooks"
adult.fixed3$building_mom[which(
  adult.fixed3$Site_ind2=="Make Believe")] <- "Make Believe"
adult.fixed3$building_mom[which(
  adult.fixed3$Site_ind2=="Reinarz")] <- "Reinarz shed est"
adult.fixed3$building_mom[which(
  adult.fixed3$Site_ind2=="Schaaps")] <- "Schaaps"
adult.fixed3$building_mom[which(
  adult.fixed3$Site_ind2=="Struthers")] <- "Struthers main"
adult.fixed3$building_mom[which(
  adult.fixed3$Site_ind2=="Urban Farm Girlz")] <- "Urban Farm Girlz Sheep Shed"


# corrections for out buildings
# Note: no clutches sampled from MB-61 so no female there
# Urban Farm 2 and 5 in sheep shed, 4 and 7 in red barn
# Cathys-3 was in trailer and then at Mary Ann's
adult.fixed3$building_mom[which(
  adult.fixed3$FamilyID_ind2=="UrbanFarm-04")] <- "Urban Farm Girlz Red Barn"
adult.fixed3$building_mom[which(
  adult.fixed3$FamilyID_ind2=="Cathys-03")] <- "Mary Anns est"

# Blue Cloud 19 birds should occur twice, once in middle shed and once in wash room
# First change current set of moms to middle shed
adult.fixed3$building_mom[which(
  adult.fixed3$FamilyID_ind2=="BlueCloud-19")] <- "Blue Cloud Middle Shed"
# next duplicate these rows and change building to washroom
bc19.mom <- subset(adult.fixed3, adult.fixed3$FamilyID_ind2=="BlueCloud-19")
bc19.mom$building_mom <- "Blue Cloud wash room"

# Also need to add duplicate location for Cooks-31 male
co31.mom <- subset(adult.fixed3, adult.fixed3$FamilyID_ind2=="Cooks-31")
co31.mom$building_mom <- "Hepp awning"

# add duplicate location for Cathys-3 in the trailer est
ca03.mom <- subset(adult.fixed3, adult.fixed3$FamilyID_ind2=="Cathys-03")
ca03.mom$building_mom <- "Cathys trailer est"


# combine the duplicate locations with the main table
adult.fixed4 <- rbind(adult.fixed3, bc19.mom, co31.mom, ca03.mom)

## Remove birds that are in the same building (social pairs and same-site EP)
adult.fixed5 <- subset(adult.fixed4, adult.fixed4$building_dad !=
                         adult.fixed4$building_mom)

## add distances
adult.fixed.dist <- left_join(adult.fixed5, dist[,c(1,3,4)], 
                          by=c("building_dad","building_mom"))

## histogram of available and realized between-site fert distances 2021
ggplot(adult.fixed.dist, aes(x=dist)) + 
  geom_histogram(binwidth=500, fill="lightblue", color="black") +
  geom_histogram(data=ep.diff.dist, aes(x=dist),
                 binwidth=500, fill="darkblue", color="black") +
  xlab("Distance in meters") +
  ggtitle("Available and realized between-site fertilizations 2021")

ggsave("output-files/available and realized fert dist 2021 histogram.png",
       h=5, w=5)

## density plots
ggplot(adult.fixed.dist, aes(x=dist)) + 
  geom_density(fill="lightblue", color="black") +
  geom_density(data=ep.diff.dist, aes(x=dist),
                 fill="darkblue", color="black", alpha=0.4) +
  xlab("Distance in meters") +
  geom_rug(data=ep.diff.dist, aes(x=dist))

ggsave("output-files/available and realized fert dist 2021 density.png",
       h=2, w=3)


# Add total analyzed for paternity 2021 ----------------------------------------

# calculate total offspring analyzed for paternity for each clutch
clutch.size <- po21 %>% group_by(Band_mom) %>%
  summarise(pat.clutch.size = n())

length(unique(clutch.size$Band_mom))

# add pat.clutch.size to ep table
ep.diff.dist2 <- left_join(ep.diff.dist, clutch.size, by="Band_mom")

# update clutch size for Hepp clutch to 2
ep.diff.dist2$pat.clutch.size[which(
  ep.diff.dist2$clutch_id_ind2=="Hepp (near Cook)_1_1")] <- 2


## add pat.clutch.size to available distances table

# change column name
colnames(adult.fixed.dist)[12] <- "Band_mom"

adult.fixed.dist2 <- left_join(adult.fixed.dist, clutch.size, by="Band_mom")

# update clutch sizes for females that had multiple locations
# BlueCloud-19 only clutch in washroom was analyzed. Remove rows where she 
# is in the middle shed (2850-57613)
adult.fixed.dist3 <- filter(adult.fixed.dist2, 
                            !(Band_mom == "2850-57613" & building_mom == 
                                "Blue Cloud Middle Shed"))

# Cooks-31 clutch was 4 in original nest, then 2 at Hepp
adult.fixed.dist3$pat.clutch.size[which(
  adult.fixed.dist3$building_mom=="Cooks" &
    adult.fixed.dist3$FamilyID_ind2=="Cooks-31")] <- 4

adult.fixed.dist3$pat.clutch.size[which(
  adult.fixed.dist3$building_mom=="Hepp awning" &
    adult.fixed.dist3$FamilyID_ind2=="Cooks-31")] <- 2

# Cathys-03 clutch was 4 for collected and 4 at Mary Anns
adult.fixed.dist3$pat.clutch.size[which(
  adult.fixed.dist3$Band_mom=="2850-57550")] <- 4


# Add proportion shared fert to all distances table ----------------------------

colnames(adult.fixed.dist3)[8] <- "Band_dad"

# first add number shared fert
adult.dist.fert <- left_join(adult.fixed.dist3, 
                             select(ep.diff.dist2, Band_dad, Band_mom, 
                                    building_dad, building_mom, fert),
                             by=c("Band_mom","Band_dad","building_mom",
                                  "building_dad"))

adult.dist.fert$fert[which(is.na(adult.dist.fert$fert))] <- 0

adult.dist.fert$prop.fert <- adult.dist.fert$fert / adult.dist.fert$pat.clutch.size

# add scaled distance
adult.dist.fert$dist.scale <- scale(adult.dist.fert$dist)


################################################################################
############################### 2022 ###########################################
################################################################################

#-------------------------------------------------------------------------------
# Add distances to between-site EP ferts for 2022
#-------------------------------------------------------------------------------

ep.all.22 <- subset(assign22, assign22$fert_type != "wp" &
                      assign22$fert_type != "mom_unk")

sum(ep.all.22$fert) # 106

table(ep.all.22$fert_type)
sum(subset(ep.all.22$fert, ep.all.22$fert_type=="ep_same")) # 60
sum(subset(ep.all.22$fert, ep.all.22$fert_type=="ep_diff")) # 8
sum(subset(ep.all.22$fert, ep.all.22$fert_type=="dad_unk")) # 38


# Identify same property, between structure EP ferts----------------------------

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

ep22.update <- ep.all.22
ep22.update$fert_type[which(ep22.update$FamilyID_dad=="BlueCloud-19" &
                              ep22.update$FamilyID_mom=="BlueCloud-18")] <- "ep_diff"



# Add distances for known sire identity between-site ---------------------------

ep22.diff <- subset(ep22.update, ep22.update$fert_type=="ep_diff" &
                      !is.na(ep22.update$Site_dad))

# add dad_buildings
ep22.diff$building_dad <- NA
ep22.diff$building_dad[which(ep22.diff$Site_dad=="Starlight")] <- "Starlight est"
ep22.diff$building_dad[which(ep22.diff$Site_dad=="Blue Cloud")] <- "Blue Cloud wash room"
ep22.diff$building_dad[which(ep22.diff$Site_dad=="Jay's")] <- "Jays est"
ep22.diff$building_dad[which(ep22.diff$Site_dad=="Mayas")] <- "Mayas est"
ep22.diff$building_dad[which(ep22.diff$Site_dad=="Reinarz")] <- "Reinarz Garage est"
ep22.diff$building_dad[which(ep22.diff$Site_dad=="Dome House")] <- "Dome House est"

# add mom buildings
ep22.diff$building_mom <- NA
ep22.diff$building_mom[which(ep22.diff$Site_mom=="Blue Cloud")] <- "Blue Cloud main"
ep22.diff$building_mom[which(ep22.diff$Site_mom=="CHR")] <- "CHR"
ep22.diff$building_mom[which(ep22.diff$Site_mom=="Cooks")] <- "Cooks"
ep22.diff$building_mom[which(
  ep22.diff$Site_mom=="Urban Farm Girlz")] <- "Urban Farm Girlz Sheep Shed"

# add distances
ep22.diff.dist <- left_join(ep22.diff, dist[,c(1,3,4)], 
                          by=c("building_dad","building_mom"))

#-------------------------------------------------------------------------------
# Determine number of available distances based on adults sampled 2022
#-------------------------------------------------------------------------------

# only keep dad-mom and mom-dad rows
adult22.2 <- subset(adult22, adult22$Type_ind1=="dad" & adult22$Type_ind2=="mom" |
                      adult22$Type_ind1=="mom" & adult22$Type_ind2=="dad")

# in this table, all dads are already Ind1, and all moms are Ind2
unique(adult22.2$Type_ind1)
unique(adult22.2$Type_ind2)


# view list of sites from 2022
sort(unique(c(adult22.2$Site_ind1, adult22.2$Site_ind2)))
sort(unique(c(adult22.2$Site_ind1)))

# view list of building
sort(unique(c(dist$building_dad, dist$building_mom)))

## Match sites from 2022 with buildings
# "Blue Cloud" = "Blue Cloud main",  "Blue Cloud Gazebo",  "Blue Cloud wash room", "Blue Cloud Middle Shed", "Blue Cloud trailer est"
# "Boyer" = "Boyer"            
# "Cargill"  = "Cargill est"        
# "Cathy's" = "Cathys est"        
# "CHR" = "CHR"             
# "Cooks"  = "Cooks"          
# "Dome House"  = "Dome House est"     
# "Erica's" = "Ericas"        
# "Grizz" = "Grizz est"
# "Jame's Ditch"  = "James Ditch est"   
# "Jay's"  = "Jays est"          
# "Karen's"  = "Karens est"       
# "Make Believe" = "Make Believe", "Make Believe shed"
# "Marte's" = "Martes est"         
# "Mary Ann's" = "Mary Anns est"      
# "Mayas" = "Mayas est"          
# "McCauley" = "McCauley est"
# "Plumbumpy" = "Plumbumpy est"
# "Quick"  = "Quick est"          
# "Reinarz"  = "Reinarz Garage est"       
# "Schaaps" = "Schaaps"         
# "Sol y Sombra"   = "Sol y Sombra est"  
# "Speedwell"  = "Speedwell est"      
# "Starlight"  = "Starlight est"     
# "Struthers" = "Struthers main", "Struthers round", "Struthers back shed"       
# "Urban Farm Girlz" = "Urban Farm Girlz Red Barn", "Urban Farm Girlz Sheep Shed"


## First add simple building matches for dads-----------------------------------

adult22.2$building_dad <- NA
adult22.2$building_dad[which(adult22.2$Site_ind1=="Blue Cloud")] <- "Blue Cloud main"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Boyer")] <- "Boyer"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Cargill")] <- "Cargill est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Cathy's")] <- "Cathys est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="CHR")] <- "CHR"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Cooks")] <- "Cooks"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Dome House")] <- "Dome House est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Erica's")] <- "Ericas"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Grizz")] <- "Grizz est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Jame's Ditch")] <- "James Ditch est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Jay's")] <- "Jays est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Karen's")] <- "Karens est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Make Believe")] <- "Make Believe"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Marte's")] <- "Martes est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Mary Ann's")] <- "Mary Anns est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Mayas")] <- "Mayas est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="McCauley")] <- "McCauley est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Plumbumpy")] <- "Plumbumpy est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Quick")] <- "Quick est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Reinarz")] <- "Reinarz Garage est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Schaaps")] <- "Schaaps"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Sol y Sombra")] <- "Sol y Sombra est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Speedwell")] <- "Speedwell est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Starlight")] <- "Starlight est"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Struthers")] <- "Struthers main"
adult22.2$building_dad[which(adult22.2$Site_ind1=="Urban Farm Girlz")] <- "Urban Farm Girlz Sheep Shed"

# Blue cloud middle shed (19), wash room (24), gazebo (26), and trailer (27)
# Make Believe shed (62)
# Urban Farm sheep shed (2,5), Urban Farm Red Barn (4, 7)
# Struthers 22 (back shed) and 23 (round shed)

# Blue Cloud 19 birds should occur twice, once in middle shed and once in wash room
# First change current set of dads to middle shed
adult22.2$building_dad[which(
  adult22.2$FamilyID_ind1=="BlueCloud-19")] <- "Blue Cloud Middle Shed"
# next duplicate these rows and change building to washroom
bc22.19.dad <- subset(adult22.2, adult22.2$FamilyID_ind1=="BlueCloud-19")
bc22.19.dad$building_dad <- "Blue Cloud wash room"

# Blue Cloud gazebo only needs to be listed once
adult22.2$building_dad[which(
  adult22.2$FamilyID_ind1=="BlueCloud-26")] <- "Blue Cloud Gazebo"

# Make Believe 62 and 62.2 are in shed
adult22.2$building_dad[which(
  adult22.2$FamilyID_ind1=="MakeBelieve-62" |
    adult22.2$FamilyID_ind1=="MakeBelieve-62.2")] <- "Make Believe Shed"

# Struthers 23 is goat shed, and Struthers 3 should be duplicated at 22 in back shed
adult22.2$building_dad[which(
  adult22.2$FamilyID_ind1=="Struthers-23")] <- "Struthers round"

st22.03.dad <- subset(adult22.2, adult22.2$FamilyID_ind1=="Struthers-03")
st22.03.dad$building_dad <- "Struthers back shed"

# combine the duplicate locations with the main table
adult22.fixed <- rbind(adult22.2, bc22.19.dad, st22.03.dad)


## same process for females ----------------------------------------------------

adult22.fixed$building_mom <- NA
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Blue Cloud")] <- "Blue Cloud main"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Boyer")] <- "Boyer"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Cargill")] <- "Cargill est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Cathy's")] <- "Cathys est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="CHR")] <- "CHR"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Cooks")] <- "Cooks"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Dome House")] <- "Dome House est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Erica's")] <- "Ericas"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Grizz")] <- "Grizz est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Jame's Ditch")] <- "James Ditch est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Jay's")] <- "Jays est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Karen's")] <- "Karens est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Make Believe")] <- "Make Believe"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Marte's")] <- "Martes est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Mary Ann's")] <- "Mary Anns est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Mayas")] <- "Mayas est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="McCauley")] <- "McCauley est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Plumbumpy")] <- "Plumbumpy est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Quick")] <- "Quick est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Reinarz")] <- "Reinarz Garage est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Schaaps")] <- "Schaaps"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Sol y Sombra")] <- "Sol y Sombra est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Speedwell")] <- "Speedwell est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Starlight")] <- "Starlight est"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Struthers")] <- "Struthers main"
adult22.fixed$building_mom[which(adult22.fixed$Site_ind2=="Urban Farm Girlz")] <- "Urban Farm Girlz Sheep Shed"

# Blue Cloud 19 birds should occur twice, once in middle shed and once in wash room
# First change current set of moms to middle shed
adult22.fixed$building_mom[which(
  adult22.fixed$FamilyID_ind2=="BlueCloud-19")] <- "Blue Cloud Middle Shed"
# next duplicate these rows and change building to washroom
bc22.19.mom <- subset(adult22.fixed, adult22.fixed$FamilyID_ind2=="BlueCloud-19")
bc22.19.mom$building_mom <- "Blue Cloud wash room"


# combine the duplicate locations with the main table
adult22.fixed2 <- rbind(adult22.fixed, bc22.19.mom)


## Remove birds that are in the same building (social pairs and same-site EP)
adult22.fixed3 <- subset(adult22.fixed2, adult22.fixed2$building_dad !=
                         adult22.fixed2$building_mom)

# remove moms with no clutches analyzed for paternity

# Blue Cloud gazebo didn't have any offspring included in paternity, so 
# remove this female
# Make Believe 62 and 62.2 didn't have any offspring analyzed so remove female
# Struthers 23 is goat shed, but no offspring analyzed for this female, remove
# Struthers 3 only had offspring analyzed from 
# the main barn nest, so don't include 22 for this female
# Cargill female remove
# CHR-149 female remove

adult22.fixed3 <- subset(adult22.fixed3, 
                         adult22.fixed3$FamilyID_ind2 != "BlueCloud-26" &
                           adult22.fixed3$FamilyID_ind2 != "MakeBelieve-62" &
                           adult22.fixed3$FamilyID_ind2 != "MakeBelieve-62.2" &
                           adult22.fixed3$FamilyID_ind2 != "Struthers-23" &
                           adult22.fixed3$FamilyID_ind2 != "CHR-149" &
                           adult22.fixed3$FamilyID_ind2 != "none")


## add distances
adult22.fixed.dist <- left_join(adult22.fixed3, dist[,c(1,3,4)], 
                              by=c("building_dad","building_mom"))

## histogram of available and realized between-site fert distances 2021
ggplot(adult22.fixed.dist, aes(x=dist)) + 
  geom_histogram(binwidth=500, fill="lightblue", color="black") +
  geom_histogram(data=ep22.diff.dist, aes(x=dist),
                 binwidth=500, fill="darkblue", color="black") +
  xlab("Distance in meters") +
  ggtitle("Available and realized between-site fertilizations 2022")

ggsave("output-files/available and realized fert dist 2022 histogram.png",
       h=5, w=5)

## density plots
ggplot(adult22.fixed.dist, aes(x=dist)) + 
  geom_density(fill="lightblue", color="black") +
  geom_density(data=ep22.diff.dist, aes(x=dist),
               fill="darkblue", color="black", alpha=0.4) +
  xlab("Distance in meters") +
  geom_rug(data=ep22.diff.dist, aes(x=dist))

ggsave("output-files/available and realized fert dist 2022 density.png",
       h=2, w=3)

#-------------------------------------------------------------------------------
# Add total analyzed for paternity 2022 
#-------------------------------------------------------------------------------

# first add  band labels for the two unsampled females, Cooks-13 and CHR-011
po22$Band_mom[which(po22$FamilyID_ind2=="Cooks-13")] <- "Cooks-13_UNS"
po22$Band_mom[which(po22$FamilyID_ind2=="CHR-011")] <- "CHR-011_UNS"

# calculate total offspring analyzed for paternity for each clutch
clutch.size22 <- po22 %>% group_by(Band_mom) %>%
  summarise(pat.clutch.size = n())

length(unique(clutch.size22$Band_mom))

# update band label for Cooks-13 female because the female who was
# assigned to this nest wasn't actually the genetic mom
adult22.fixed.dist$Band_ind2[which(
  adult22.fixed.dist$FamilyID_ind2=="Cooks-13")] <- "Cooks-13_UNS"


## add pat.clutch.size to available distances table

# change column name
colnames(adult22.fixed.dist)[15] <- "Band_mom"

adult22.fixed.dist2 <- left_join(adult22.fixed.dist, clutch.size22, by="Band_mom")


# update clutch sizes for females that had multiple locations
# Blue Cloud-19 had 1 chick from shed and 5 from washroom

adult22.fixed.dist2$pat.clutch.size[which(
  adult22.fixed.dist2$FamilyID_ind2=="BlueCloud-19" & 
    adult22.fixed.dist2$building_mom=="Blue Cloud Middle Shed")] <- 1

adult22.fixed.dist2$pat.clutch.size[which(
  adult22.fixed.dist2$FamilyID_ind2=="BlueCloud-19" & 
    adult22.fixed.dist2$building_mom=="Blue Cloud wash room")] <- 5


# Add proportion shared fert to all distances table ----------------------------

colnames(adult22.fixed.dist2)[8] <- "Band_dad"

# first add number shared fert
adult22.dist.fert <- left_join(adult22.fixed.dist2, 
                             select(ep22.diff.dist, Band_dad, Band_mom, 
                                    building_dad, building_mom, fert),
                             by=c("Band_mom","Band_dad","building_mom",
                                  "building_dad"))

adult22.dist.fert$fert[which(is.na(adult22.dist.fert$fert))] <- 0

adult22.dist.fert$prop.fert <- adult22.dist.fert$fert / adult22.dist.fert$pat.clutch.size

# add scaled distance
adult22.dist.fert$dist.scale <- scale(adult22.dist.fert$dist)




################################################################################
############################# 2023 #############################################
################################################################################


#-------------------------------------------------------------------------------
# Add distances to between-site EP ferts 2023
#-------------------------------------------------------------------------------

ep.all.23 <- subset(assign23, assign23$fert_type != "wp")

sum(ep.all.23$fert) # 64

table(ep.all.23$fert_type)
sum(subset(ep.all.23$fert, ep.all.23$fert_type=="ep_same")) # 23
sum(subset(ep.all.23$fert, ep.all.23$fert_type=="ep_diff")) # 7
sum(subset(ep.all.23$fert, ep.all.23$fert_type=="ep_unsampled")) # 34

# Identify same property, between structure EP ferts----------------------------

# For 2023, the following outbuilding nests were active and sampled for paternity: 
# Blue Cloud: 28 (north shed) and 19 (middle shed) were active but not sampled
# Cooks: 15 (house) was active but not sampled
# Make Believe: no shed nests were active
#--------------------------
# Struthers: 23 (round shed)
# Urban Farm Girlz: 8, 5, 2 (sheep shed); 4, 7 (Red barn)

# Change Urban Farm 7 to Urban Farm 8 to ep_diff
ep23.update <- ep.all.23
ep23.update$fert_type[which(ep23.update$FamilyID_dad=="UrbanFarm-07" &
                              ep23.update$FamilyID_mom=="UrbanFarm-08")] <- "ep_diff"



## Add distances to known-sire between-site ferts-------------------------------

ep23.diff <- subset(ep23.update, ep23.update$fert_type=="ep_diff" &
                      !is.na(ep23.update$site_dad))

# add dad_buildings
ep23.diff$building_dad <- NA
ep23.diff$building_dad[which(
  ep23.diff$site_dad=="Blue Cloud")] <- "Blue Cloud main"
ep23.diff$building_dad[which(
  ep23.diff$FamilyID_dad=="UrbanFarm-07")] <- "Urban Farm Girlz Red Barn"
ep23.diff$building_dad[which(
  ep23.diff$site_dad=="Cathys")] <- "Cathys est"
ep23.diff$building_dad[which(
  ep23.diff$site_dad=="Green Mill")] <- "Green Mill est"

# add mom buildings
ep23.diff$building_mom <- NA
ep23.diff$building_mom[which(
  ep23.diff$site_ind2=="Struthers")] <- "Struthers main"
ep23.diff$building_mom[which(
  ep23.diff$FamilyID_mom=="UrbanFarm-08")] <- "Urban Farm Girlz Sheep Shed"
ep23.diff$building_mom[which(
  ep23.diff$site_ind2=="Boyer")] <- "Boyer"
ep23.diff$building_mom[which(
  ep23.diff$site_ind2=="Jays")] <- "Jays est"

# add distances
ep23.diff.dist <- left_join(ep23.diff, dist[,c(1,3,4)], 
                            by=c("building_dad","building_mom"))

# Remove Urban Farm Girlz dads because IDobs were too uncertain even to assign
# dads to one structure or the other

ep23.diff.dist2 <- subset(ep23.diff.dist, 
                          ep23.diff.dist$site_dad != "Urban Farm Girlz")



#-------------------------------------------------------------------------------
# Determine number of available distances based on adults sampled 2023
#-------------------------------------------------------------------------------

# only keep dad-mom and mom-dad rows
adult23.2 <- subset(adult23, adult23$type_ind1=="dad" & adult23$type_ind2=="mom" |
                      adult23$type_ind1=="mom" & adult23$type_ind2=="dad")

## make all dads ind1 and all moms ind2

adult23.dad.mom <- subset(adult23.2, adult23.2$type_ind1=="dad" &
                          adult23.2$type_ind2=="mom")

adult23.mom.dad <- subset(adult23.2, adult23.2$type_ind1=="mom" &
                          adult23.2$type_ind2=="dad")

# pull out ind1 info columns and change to ind2
adult23.mom.ind1 <- adult23.mom.dad[,c(1:9, 11:12)]
colnames(adult23.mom.ind1)[8:11] <- gsub("ind1","ind2", colnames(adult23.mom.ind1)[8:11])
colnames(adult23.mom.ind1)[1:2] <- c("Ind2","Ind1")

# pull out ind2 cols and change to ind1
adult23.dad.ind2 <- adult23.mom.dad[,c(1:7,15:16, 18:19)]
colnames(adult23.dad.ind2)[8:11] <- gsub("ind2","ind1", colnames(adult23.dad.ind2)[8:11])
colnames(adult23.dad.ind2)[1:2] <- c("Ind2","Ind1")

# combine back together
adult23.fixed <- full_join(adult23.dad.ind2, adult23.mom.ind1, 
                         by=c("Ind1","Ind2","k0_hat","k1_hat","k2_hat",
                              "pi_HAT","nbSNP"))

# add in rows that were correct to start with
adult23.fixed2 <- rbind(adult23.fixed, adult23.dad.mom[,-c(10,13,14,17,20,21)])

# remove duplicate rows caused by males who attended more than one social nest
adult23.fixed3 <- adult23.fixed2[-which(duplicated(adult23.fixed2[,1:2])), ]

## First add simple building matches for dads-----------------------------------

sort(unique(adult23$site_ind1))

adult23.fixed3$building_dad <- NA
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Blue Cloud")] <- "Blue Cloud main"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Boyer")] <- "Boyer"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Cathys")] <- "Cathys est"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="CHR")] <- "CHR"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Cooks")] <- "Cooks"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Green Mill")] <- "Green Mill est"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Jays")] <- "Jays est"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Karens")] <- "Karens est"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Make Believe")] <- "Make Believe"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Plumbumpy")] <- "Plumbumpy est"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Quick")] <- "Quick est"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Red Wagon")] <- "Red Wagon porch est"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Speedwell")] <- "Speedwell est"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Struthers")] <- "Struthers main"
adult23.fixed3$building_dad[which(adult23.fixed3$site_ind1=="Urban Farm Girlz")] <- "Urban Farm Girlz Sheep Shed"

# update outbuildings for active nests
# Blue Cloud-28 in north shed, pair then moved to Nixon later so will need duplicate rows
# Blue Cloud males 57926, 57927, and 57930 never assigned to nests, but caught in main barn
# Struthers 23 in round shed
# UrbanFarm-07 in Red Barn
adult23.fixed3$building_dad[which(adult23.fixed3$FamilyID_ind1=="Nixon-01")] <- "Blue Cloud pony shed"
adult23.fixed3$building_dad[which(adult23.fixed3$FamilyID_ind1=="Struthers-23")] <- "Struthers round"
adult23.fixed3$building_dad[which(adult23.fixed3$FamilyID_ind1=="UrbanFarm-07")] <- "Urban Farm Girlz Red Barn"

# make duplicate rows for Nixon-01
nixon.dad <- subset(adult23.fixed3, adult23.fixed3$FamilyID_ind1=="Nixon-01")
nixon.dad$building_dad <- "Nixon est"

# add back to dads
adult23.fixed4 <- rbind(adult23.fixed3, nixon.dad)

## Same process for building_mom------------------------------------------------

sort(unique(adult23.fixed4$site_ind2))

adult23.fixed4$building_mom <- NA
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Blue Cloud")] <- "Blue Cloud main"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Boyer")] <- "Boyer"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Cathys")] <- "Cathys est"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="CHR")] <- "CHR"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Cooks")] <- "Cooks"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Jays")] <- "Jays est"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Karens")] <- "Karens est"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Make Believe")] <- "Make Believe"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Plumbumpy")] <- "Plumbumpy est"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Quick")] <- "Quick est"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Red Wagon")] <- "Red Wagon porch est"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Speedwell")] <- "Speedwell est"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Struthers")] <- "Struthers main"
adult23.fixed4$building_mom[which(adult23.fixed4$site_ind2=="Urban Farm Girlz")] <- "Urban Farm Girlz Sheep Shed"

# update outbuildings for active nests
# Blue Cloud-28 in north shed, pair then moved to Nixon later. Only Nixon nest
# analyzed so only include this location for the female 2850-57946

# Struthers 23 in round shed
# Urban Farm 7 in red barn
adult23.fixed4$building_mom[which(adult23.fixed4$FamilyID_ind2=="Nixon-01")] <- "Nixon est"
adult23.fixed4$building_mom[which(adult23.fixed4$FamilyID_ind2=="Struthers-23")] <- "Struthers round"
adult23.fixed4$building_mom[which(adult23.fixed4$FamilyID_ind2=="UrbanFarm-07")] <- "Urban Farm Girlz Red Barn"


# Remove Urban Farm Girlz dads because IDobs were too uncertain even to assign
# dads to one structure or the other
adult23.fixed5 <- subset(adult23.fixed4, 
                         adult23.fixed4$site_ind1 != "Urban Farm Girlz")


## Remove birds that are in the same building (social pairs and same-site EP)
adult23.fixed6 <- subset(adult23.fixed5, adult23.fixed5$building_dad !=
                           adult23.fixed5$building_mom)



## add distances----------------------------------------------------------------
adult23.fixed.dist <- left_join(adult23.fixed6, dist[,c(1,3,4)], 
                                by=c("building_dad","building_mom"))

## histogram of available and realized between-site fert distances 2023
ggplot(adult23.fixed.dist, aes(x=dist)) + 
  geom_histogram(binwidth=500, fill="lightblue", color="black") +
  geom_histogram(data=ep23.diff.dist, aes(x=dist),
                 binwidth=500, fill="darkblue", color="black") +
  xlab("Distance in meters") +
  ggtitle("Available and realized between-site fertilizations 2023")

ggsave("output-files/available and realized fert dist 2023 histogram.png",
       h=5, w=5)

## density plots
ggplot(adult23.fixed.dist, aes(x=dist)) + 
  geom_density(fill="lightblue", color="black") +
  geom_density(data=ep23.diff.dist, aes(x=dist),
               fill="darkblue", color="black", alpha=0.4) +
  xlab("Distance in meters") +
  geom_rug(data=ep23.diff.dist, aes(x=dist))

ggsave("output-files/available and realized fert dist 2023 density.png",
       h=2, w=3)

#-------------------------------------------------------------------------------
# Add total analyzed for paternity 2023 
#-------------------------------------------------------------------------------

# calculate total offspring analyzed for paternity for each clutch
clutch.size23 <- po23 %>% group_by(band_mom) %>%
  summarise(pat.clutch.size = n())

length(unique(clutch.size23$band_mom))


## add pat.clutch.size to available distances table

# change column name
colnames(adult23.fixed.dist)[12] <- "band_mom"

adult23.fixed.dist2 <- left_join(adult23.fixed.dist, clutch.size23, by="band_mom")


# Add proportion shared fert to all distances table ----------------------------

colnames(adult23.fixed.dist2)[8] <- "band_dad"

# first add number shared fert
adult23.dist.fert <- left_join(adult23.fixed.dist2, 
                               select(ep23.diff.dist2, band_dad, band_mom, 
                                      building_dad, building_mom, fert),
                               by=c("band_mom","band_dad","building_mom",
                                    "building_dad"))


adult23.dist.fert$fert[which(is.na(adult23.dist.fert$fert))] <- 0

adult23.dist.fert$prop.fert <- adult23.dist.fert$fert / adult23.dist.fert$pat.clutch.size

# add scaled distance
adult23.dist.fert$dist.scale <- scale(adult23.dist.fert$dist)

################################################################################
###################### Combine years ###########################################
################################################################################

# combine years-----------------------------------------------------------------

adult.dist.fert$year <- 2021
adult22.dist.fert$year <- 2022
adult23.dist.fert$year <- 2023

colnames(adult23.dist.fert)[c(8,12)] <- c("Band_dad","Band_mom")

all.dist.fert <- rbind(select(adult.dist.fert, Band_dad, Band_mom, 
                              building_dad, building_mom, dist, pat.clutch.size,
                              fert, prop.fert, year), 
                       select(adult22.dist.fert, Band_dad, Band_mom, 
                              building_dad, building_mom, dist, pat.clutch.size,
                              fert, prop.fert, year),
                       select(adult23.dist.fert, Band_dad, Band_mom, 
                              building_dad, building_mom, dist, pat.clutch.size,
                              fert, prop.fert, year))

# save final table used for model
write.csv(all.dist.fert, "03_output-files/all_dist_fert_between_site_final.csv", row.names = F)


# Calculate bin cutoff----------------------------------------------------------

# Check average min and max nest distances for available
all.summary <- all.dist.fert %>% group_by(Band_mom) %>%
  summarise(min.neighbor = min(dist),
            max.neighbor = max(dist),
            mean.neighbor = mean(dist),
            sd.neighbor = sd(dist))

# look at averages for the summary stats
all.summary2 <- data.frame(min.neighbor.mean = mean(all.summary$min.neighbor),
                              max.neighbor.mean = mean(all.summary$max.neighbor),
                              mean.neighbor.mean = mean(all.summary$mean.neighbor),
                              sd.neighbor.mean = mean(all.summary$sd.neighbor))

# breaks at 
# 0-half of min.neighbor.mean: 0-546
# next break at min.neighbor.mean: 546-1092
# next break at (min.neighbor.mean + sd.neighbor.mean): 1092-4052
# next break at mean.neighbor.mean: 4052-6305
# next break at max.neighbor.mean: 6305-11323
# last break larger than max: 11323+

# add breaks to dataframe
all.dist.fert$dist.bin6 <- NA
all.dist.fert$dist.bin6[which(all.dist.fert$dist <= 546)] <- "0-546"
all.dist.fert$dist.bin6[which(
  all.dist.fert$dist > 546 & all.dist.fert$dist <= 1092)] <- "546-1092"
all.dist.fert$dist.bin6[which(
  all.dist.fert$dist > 1092 & all.dist.fert$dist <= 4052)] <- "1092-4052"
all.dist.fert$dist.bin6[which(
  all.dist.fert$dist > 4052 & all.dist.fert$dist <= 6305)] <- "4052-6305"
all.dist.fert$dist.bin6[which(
  all.dist.fert$dist > 6305 & all.dist.fert$dist <= 11323)] <- "6305-11323"
all.dist.fert$dist.bin6[which(
  all.dist.fert$dist > 11323)] <- "11323+"

all.dist.fert$dist.bin6 <- factor(all.dist.fert$dist.bin6, 
                                      levels=c("0-546", "546-1092",
                                               "1092-4052", "4052-6305",
                                               "6305-11323", "11323+"))

# first three bins have ferts, collapse others
all.dist.fert$dist.bin3 <- NA
all.dist.fert$dist.bin3[which(all.dist.fert$dist <= 546)] <- "0-546"
all.dist.fert$dist.bin3[which(
  all.dist.fert$dist > 546 & all.dist.fert$dist <= 1092)] <- "546-1092"
all.dist.fert$dist.bin3[which(
  all.dist.fert$dist > 1092)] <- "1092+"

all.dist.fert$dist.bin3 <- factor(all.dist.fert$dist.bin3, 
                                  levels=c("0-546", "546-1092",
                                           "1092+"))

# histogram with bins and realized ferts ---------------------------------------
# column for only realized fert dist
all.dist.fert$realized.dist <- ifelse(all.dist.fert$prop.fert > 0, 
                                          all.dist.fert$dist, NA)

# histogram with dot plot
ggplot(all.dist.fert, aes(x=dist)) + 
  geom_histogram(fill="lightblue", color="black", binwidth = 546,
                 origin=0) +
  geom_dotplot(aes(x=realized.dist), 
               method="histodot", binwidth = 546, origin=0,
               alpha=0.6) +
  geom_vline(xintercept=546, color="red") +
  geom_vline(xintercept = 1092, color="red") +
  geom_vline(xintercept = 4052, linetype="dashed", color="red") +
  geom_vline(xintercept = 6305, linetype="dashed", color="red") +
  geom_vline(xintercept = 11323, linetype="dashed", color="red") +
  xlab("Distance between structures in meters") +
  theme_light()

ggsave("03_output-files/hist-dot-all.png", h=4, w=7)


#-------------------------------------------------------------------------------
# Fit GLMM of proportion shared ferts between sites all years
#-------------------------------------------------------------------------------


# 3 bins -----------------------------------------------------------------------
prop.ep.all.bin3 <- glmmTMB(prop.fert ~ (1|Band_mom) + (1|Band_dad) +
                              dist.bin3 + as.factor(year),
                            data=all.dist.fert, family=binomial,
                            weights=pat.clutch.size)

summary(prop.ep.all.bin3)
# Family: binomial  ( logit )
# Formula:          
#   prop.fert ~ (1 | Band_mom) + (1 | Band_dad) + dist.bin3 + as.factor(year)
# Data: all.dist.fert
# Weights: pat.clutch.size
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 223.1     272.5    -104.5     209.1      8670 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance Std.Dev.
# Band_mom (Intercept) 48.71    6.979   
# Band_dad (Intercept) 47.47    6.890   
# Number of obs: 8677, groups:  Band_mom, 100; Band_dad, 169
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -14.3337     2.9259  -4.899 9.63e-07 ***
#   dist.bin3546-1092    -0.3903     1.1824  -0.330   0.7413    
# dist.bin31092+       -5.3177     1.1886  -4.474 7.68e-06 ***
#   as.factor(year)2022  -0.7315     1.0471  -0.699   0.4848    
# as.factor(year)2023  -4.8616     2.1436  -2.268   0.0233 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

modsum.confint <- as.data.frame(confint(prop.ep.all.bin3))
modsum.bt <- exp(modsum.confint)
colnames(modsum.bt) <- c("lower.bt","upper.bt","estimate.bt")

modsum.save <- cbind(modsum.confint, modsum.bt)

write.csv(modsum.save, "03_output-files/binned dist mod sum all.csv")


car::Anova(prop.ep.all.bin3)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: prop.fert
# Chisq Df Pr(>Chisq)    
# dist.bin3       22.5549  2  1.265e-05 ***
#   as.factor(year)  5.1621  2     0.0757 .  


# model diagnostics-------------------------------------------------------------

# calculate simulated residuals
prop.all.simResids <- simulateResiduals(prop.ep.all.bin3)

# no obvious issues
plotQQunif(prop.all.simResids) 
plotResiduals(prop.all.simResids)

# residuals against dist
plotResiduals(prop.all.simResids, form=all.dist.fert$dist.bin3)

## try plotting residuals by grouping variables

# female
plotResiduals(prop.all.simResids, form=all.dist.fert$Band_mom)

# male 
plotResiduals(prop.all.simResids, form=all.dist.fert$Band_dad)

#-------------------------------------------------------------------------------
# Use log link instead of logit
#-------------------------------------------------------------------------------

# Model would not run with default starting values. 
# Specify starting values from the logit model coefficients

log_starts <- list(
  theta = c(sqrt(48), sqrt(47)),
  beta = c(-14, -0.4, -5.3, -0.7, -4.9)
)

prop.ep.all.log <- glmmTMB(prop.fert ~ (1|Band_mom) + (1|Band_dad) +
                             dist.bin3 + as.factor(year),
                           data=all.dist.fert, family=binomial(link="log"),
                           weights=pat.clutch.size, 
                           start = log_starts)

summary(prop.ep.all.log)


#-------------------------------------------------------------------------------
# Calculate mean available and realized distances for all years
#-------------------------------------------------------------------------------

mean(all.dist.fert$dist) # 6015
mean(all.dist.fert$realized.dist, na.rm=T) # 1011


# For each year separately, for shared fert = Y and shared fert = N ------------

# 2021 
mean(subset(adult.dist.fert$dist, adult.dist.fert$fert > 0)) # yes = 1456.003
mean(subset(adult.dist.fert$dist, adult.dist.fert$fert == 0)) # no = 8122.588

# 2022
mean(subset(adult22.dist.fert$dist, adult22.dist.fert$fert > 0)) # yes = 724.478
mean(subset(adult22.dist.fert$dist, adult22.dist.fert$fert == 0)) # no = 5811.282

# 2023
mean(subset(adult23.dist.fert$dist, adult23.dist.fert$fert > 0)) # yes = 993.5637
mean(subset(adult23.dist.fert$dist, adult23.dist.fert$fert == 0)) # no = 5555.975



#-------------------------------------------------------------------------------
# Calculate number of clutches, moms, and sires sampled at each site and year
#-------------------------------------------------------------------------------

# 2021 -------------------------------------------------------------------------

## Number of females and males sampled
adults21.sampled.male <- adult.fixed.dist3 %>% 
  group_by(building_dad) %>%
  summarise(n.males = length(unique(Band_dad)))
colnames(adults21.sampled.male)[1] <- "site"

adults21.sampled.fem <- adult.fixed.dist3 %>% 
  group_by(building_mom) %>%
  summarise(n.females = length(unique(Band_mom)))
colnames(adults21.sampled.fem)[1] <- "site"

adults21.sampled <- left_join(adults21.sampled.male, adults21.sampled.fem,
                              by="site")
# replace NAs with zeros
adults21.sampled$n.females[which(is.na(adults21.sampled$n.females))] <- 0


## number of clutches sampled
clutch21.sampled <- kin.assign %>%
  group_by(Site_mom) %>%
  summarise(clutch.per.site = length(unique(clutch_id_ind2)))

# make corrections for Mary Anns, Hepp, and Urban Farm, Blue Cloud
clutch21.sampled2 <- rbind(clutch21.sampled, 
                           c("Mary Anns est",1),
                           c("Hepp awning",1),
                           c("Urban Farm Girlz Sheep Shed",3),
                           c("Urban Farm Girlz Red Barn", 1),
                           c("Cathys trailer est",1),
                           c("Cathys est",1),
                           c("Blue Cloud main",4),
                           c("Blue Cloud wash room", 1),
                           c("Reinarz shed est", 1),
                           c("Struthers main", 3))

clutch21.sampled2$clutch.per.site[which(
  clutch21.sampled2$Site_mom=="Cooks")] <- 11


## Combine clutches with adults
colnames(clutch21.sampled2)[1] <- "site"
sampled21 <- left_join(adults21.sampled, clutch21.sampled2, 
                       by="site")
sampled21$clutch.per.site[which(is.na(sampled21$clutch.per.site))] <- 0

# add sample type column
sampled21$sample.type <- ifelse(sampled21$clutch.per.site>0,
                                "focal nests","sires only")

# save file
write.csv(sampled21, 
          "03_output-files/clutches females males sampled in 2021 by site.csv",
          row.names = F)

# 2022 -------------------------------------------------------------------------

# number of females and males sampled by structure
adult22.sampled.male <- adult22.fixed.dist2 %>%
  group_by(building_dad) %>%
  summarise(n.male = length(unique(Band_dad)))
colnames(adult22.sampled.male)[1] <- "building"

adult22.sampled.fem <- adult22.fixed.dist2 %>%
  group_by(building_mom) %>%
  summarise(n.female = length(unique(Band_mom)))
colnames(adult22.sampled.fem)[1] <- "building"

adult22.sampled <- full_join(adult22.sampled.male, adult22.sampled.fem,
                             by="building")

# replace NAs with zeros
adult22.sampled$n.female[which(is.na(adult22.sampled$n.female))] <- 0

## number of clutches sampled

# fix NAs in assign22
assign22$Site_mom[which(assign22$clutch_id_ind2=="CHR_101_2")] <- "CHR"
assign22$Site_mom[which(assign22$clutch_id_ind2=="Cook's_13_1")] <- "Cooks"

clutch22.sampled <- assign22 %>%
  group_by(Site_mom) %>%
  summarise(clutch.per.site = length(unique(clutch_id_ind2)))

# add correct building names and outbuildings
clutch22.sampled2 <- rbind(clutch22.sampled, 
                           c("Blue Cloud Middle Shed", 1),
                           c("Blue Cloud wash room", 1),
                           c("Blue Cloud main", 8),
                           c("Cathys est", 1),
                           c("Dome House est", 2),
                           c("Karens est", 1),
                           c("Martes est", 2),
                           c("Mary Anns est", 1),
                           c("McCauley est", 3),
                           c("Speedwell est", 1),
                           c("Struthers main", 3),
                           c("Urban Farm Girlz Sheep Shed", 3))

## Combine clutches with adults
colnames(clutch22.sampled2)[1] <- "building"
sampled22 <- left_join(adult22.sampled, clutch22.sampled2, 
                       by="building")
sampled22$clutch.per.site[which(is.na(sampled22$clutch.per.site))] <- 0

# add sample type column
sampled22$sample.type <- ifelse(sampled22$clutch.per.site>0,
                                "focal nests","sires only")

# save file
write.csv(sampled22, 
          "03_output-files/clutches females males sampled in 2022 by site.csv",
          row.names = F)


# 2023 -------------------------------------------------------------------------

# number of females and males sampled by structure
adult23.sampled.male <- adult23.fixed.dist2 %>%
  group_by(building_dad) %>%
  summarise(n.male = length(unique(band_dad)))
colnames(adult23.sampled.male)[1] <- "building"

adult23.sampled.fem <- adult23.fixed.dist %>%
  group_by(building_mom) %>%
  summarise(n.female = length(unique(band_mom)))
colnames(adult23.sampled.fem)[1] <- "building"

adult23.sampled <- full_join(adult23.sampled.male, adult23.sampled.fem,
                             by="building")

# replace NAs with zeros
adult23.sampled$n.female[which(is.na(adult23.sampled$n.female))] <- 0
adult23.sampled$n.male[which(is.na(adult23.sampled$n.male))] <- 0
 
## number of clutches sampled
clutch23.sampled <- assign23 %>%
  group_by(site_ind2) %>%
  summarise(clutch.per.site = length(unique(clutchID_ind2)))

# update with proper building names and outbuilding clutches
clutch23.sampled2 <- rbind(clutch23.sampled, 
                           c("Blue Cloud main", 7),
                           c("Nixon est", 1),
                           c("Cathys est", 2),
                           c("Jays est", 2),
                           c("Karens est", 1),
                           c("Plumbumpy est", 1),
                           c("Red Wagon porch est", 1),
                           c("Speedwell est", 2),
                           c("Struthers main", 4),
                           c("Struthers round", 1),
                           c("Urban Farm Girlz Red Barn", 1),
                           c("Urban Farm Girlz Sheep Shed", 3))


## Combine clutches with adult
colnames(clutch23.sampled2)[1] <- "building"
sampled23 <- left_join(adult23.sampled, clutch23.sampled2,
                       by="building")
sampled23$clutch.per.site[which(is.na(sampled23$clutch.per.site))] <- 0

# add sample type column
sampled23$sample.type <- ifelse(sampled23$clutch.per.site>0,
                                "focal nests","sires only")

# save file
write.csv(sampled23, 
          "03_output-files/clutches females males sampled in 2023 by site.csv",
          row.names = F)


