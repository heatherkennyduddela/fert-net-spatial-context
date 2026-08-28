
################################################################################
# Script to determine temporal overlap of nest activity 2022 CO
# Heather Kenny-Duddela
# Feb 21, 2025
################################################################################

# libraries
library(tidyverse) # included ggplot2 and lubridate for formatting dates
library(reshape2) # for switching between matrix and edge list

## load data
fert.female3 <- read.csv("input-files/fert start end each female each nest.csv")
kin.ci.nest <- read.csv("input-files/paternity with ci dates and diff 2022.csv")

# matrix of fertile days from 03_timing and synchrony > 02_calculate breeding
# synchrony index.R
mat <- read.csv("input-files/season matrix multi sites.csv")

# distance file, generated in directory:
# CU Boulder/BARS fieldwork/2022 Field and Lab/Nest maps and Sketch Up/pairwise dist from ArcGIS
dist <- read.csv("input-files/distances with BarnNest ids.csv") 

# family IDs table, generated manually
dads <- read.csv("input-files/fam_clutch_2022_updated broods.csv")

# date format
fert.female3$fert.start1 <- ymd(fert.female3$fert.start1)
fert.female3$fert.end1 <- ymd(fert.female3$fert.end1)
fert.female3$fert.start2 <- ymd(fert.female3$fert.start2)
fert.female3$fert.end2 <- ymd(fert.female3$fert.end2)
fert.female3$fert.start3 <- ymd(fert.female3$fert.start3)
fert.female3$fert.end3 <- ymd(fert.female3$fert.end3)

#-------------------------------------------------------------------------------
# We want to create a spatial relationship matrix where each row is a clutch
# (rather than each row being a female-male pair). Basically, each female-male
# pair will have two rows, one for each clutch, assuming both were present during
# 1st and 2nd clutches. Some birds were only present during 2nd clutches so will 
# not be included in 1st clutch pairs. Females that only laid a 1st clutch will not
# be included in the 2nd clutch pairs. 

# Because birds can switch nests between attempts, having each row be a clutch
# allows for different spatial relationships during 1st and 2nd broods.
#-------------------------------------------------------------------------------

## add dates of nest activity for classifying spatial relationships
# Fertile period was defined as 7 days before the CI until the day the 
# penultimate egg was laid. 
# Nest is active from start of fertile period until one day before the start
# of the fertile period for the female's next brood. For 3rd broods, active
# period ends on estimated nestling day 18, fert end + 32 days

fert.female3$active.start1 <- fert.female3$fert.start1
fert.female3$active.end1 <- fert.female3$fert.start2 - 1

fert.female3$active.start2 <- fert.female3$fert.start2
fert.female3$active.end2 <- fert.female3$fert.start3 - 1

fert.female3$active.start3 <- fert.female3$fert.start3
fert.female3$active.end3 <- fert.female3$fert.end3 + 32

# For birds without 2nd or 3rd broods, active end on nestling day 18 (estimated) 
# Estimate as fert end + 1 + 13 days + 18 days
fert.female3$active.end1[which(is.na(fert.female3$active.end1))] <- 
  fert.female3$fert.end1[which(is.na(fert.female3$active.end1))] + 32

fert.female3$active.end2[which(is.na(fert.female3$active.end2) &
                                 !is.na(fert.female3$active.start2))] <- 
  fert.female3$fert.end2[which(is.na(fert.female3$active.end2) &
                                 !is.na(fert.female3$active.start2))] + 32


# Reshape so there is one row per nest attempt ---------------------------------
#-------------------------------------------------------------------------------

nest1 <- select(fert.female3, band, site, nest1, fert.start1, fert.end1,
                active.start1, active.end1)

nest2 <- select(fert.female3, band, site, nest2, fert.start2, fert.end2,
                active.start2, active.end2)
# remove NA rows
nest2 <- filter(nest2, !is.na(fert.start2))
# make nests integer
nest2$nest2 <- as.integer(nest2$nest2)

nest3 <- select(fert.female3, band, site, nest3, fert.start3, fert.end3,
                active.start3, active.end3)
# remove NA rows
nest3 <- filter(nest3, !is.na(fert.start3))
# make nests integer
nest3$nest3 <- as.integer(nest3$nest3)


## harmonize column names
colnames(nest1)[3:7] <- c("nest","fert.start","fert.end",
                          "active.start", "active.end")
colnames(nest2)[3:7] <- c("nest","fert.start","fert.end",
                          "active.start", "active.end")
colnames(nest3)[3:7] <- c("nest","fert.start","fert.end",
                          "active.start", "active.end")
# add column for attempt
nest1$attempt <- 1
nest2$attempt <- 2
nest3$attempt <- 3

## combine again
nest.all <- rbind(nest1, nest2, nest3)

# CI is fert start + 7 days
nest.all$ci <- nest.all$fert.start + 7

# add BarnNest ID
nest.all$BarnNest <- paste(nest.all$site, nest.all$nest, sep="-")

#-------------------------------------------------------------------------------
# Add social males to each nest attempt
#-------------------------------------------------------------------------------

## for this we will use the family table

# first keep only relevant columns
dads2 <- select(dads, Band, Site, Type, FamilyID)

# filter only moms and dads, and check for duplicated bands
dads3 <- filter(dads2, Type=="mom" | Type=="dad")
dads3$dup <- duplicated(dads3$Band)

## Deal with dads that were attending two social nests
# 2850-57806 was at both CHR 96 and CHR 37
# 2850-57582 was at both CHR 148 and CHR 74
# For both cases, check which nest initiated earlier, and call this the male's 
# location in the barn. 
# Mom for CHR 148 is 2870-64909, mom for CHR 74 is 2850-57746
# Mom for CHR 96 is 2850-57568, mom for CHR 37 is 2850-57735

# CI for CHR 148 is 2022-07-30, CI for CHR 74 is 2022-05-30 and for 2nd brood
# of that female (CHR 6) is 2022-07-21. Nest 6 should be the 2nd brood location
# for male 2850-57582. 

# CI for CHR 96 is 2022-05-18 and 2022-07-10. CI for nest 37 is 2022-05-18
# and 2022-07-09. Will set nest 96 as location for brood1 and 37 as location
# for brood2 for male 2850-57806


### for non-duplicated dads, we can link them using the mom IDs

# change table structure to have one column per family
# this is messy because some dads don't have familyIDs, and all get lumped
# together under "none". Try using this and see if it works
dads4 <- pivot_wider(dads3, names_from = Type, values_from = Band)
# convert any list rows to strings
dads5 <- dads4 %>% mutate(dad = sapply(dad, toString),
                          mom = sapply(mom, toString))

dads6 <- dads5 %>% filter(nchar(dad) <= 10 & nchar(dad) > 0)

## link dad bands to mom bands in nest.all

# change nest.all colnames
colnames(nest.all)[1] <- "mom"

# now add dads
nest.all2 <- left_join(nest.all, dads6[,4:5], by="mom")

# Some dads are missing because they were unbanded:
# - CHR 7 and 9, female 2850-57584
# - CHR 77, female 2850-57696
# - CHR 76 and 12, female 2850-57737
# - Blue Cloud-27, female 2850-57886
# - CHR-145, female 2870-65004
# Add in unique labels for these UNB males
nest.all2$dad[which(nest.all2$mom=="2850-57584")] <- "unb-CHR-7and9"
nest.all2$dad[which(nest.all2$mom=="2850-57696")] <- "unb-CHR-77"
nest.all2$dad[which(nest.all2$mom=="2850-57737")] <- "unb-CHR-76and12"
nest.all2$dad[which(nest.all2$mom=="2850-57886")] <- "unb-BC-27"
nest.all2$dad[which(nest.all2$mom=="2870-65004")] <- "unb-CHR-145"


# Some dads are missing but are banded and can be added: 
# - CHR 96, 2850-57806 (use nest 96 as his location for brood1)
# - CHR 82 and 115, female 57655, male 2850-57745
# - Sol y Sombra, female 2850-57641, male 2850-57642
# - Struthers-23, female 2850-57683, male 2850-57684
# - Erica's, female 2850-57741, male 2811-00632
# - Blue Cloud-26, female 2850-57885, male 2870-64905
# - CHR-11 and CHR-101, female 2870-64908, male 2850-57733
# - CHR-148, female 2870-64909, male 2850-57582 who also attended nest 6

# One nest had two different dads for brood1 and brood2
# - Make Believe 62, female 57876, brood1 male 2850-57680, brood2 male 2850-57750

## add in missing dads, and make sure double nesters are correct: 

# add male from CHR 96, female there is 2850-57568
nest.all2$dad[which(nest.all2$mom=="2850-57568")] <- "2850-57806"

# add male from 82 and 115 (82 is location for brood1)
nest.all2$dad[which(nest.all2$mom=="2850-57655")] <- "2850-57745"

# Sol y Sombra male
nest.all2$dad[which(nest.all2$mom=="2850-57641")] <- "2850-57642"

# Struthers 23 male
nest.all2$dad[which(nest.all2$mom=="2850-57683")] <- "2850-57684"

# Erica's male
nest.all2$dad[which(nest.all2$mom=="2850-57741")] <- "2811-00632"

# Blue Cloud-26 male
nest.all2$dad[which(nest.all2$mom=="2850-57885")] <- "2870-64905"

# CHR 11 and 101
nest.all2$dad[which(nest.all2$mom=="2870-64908")] <- "2850-57733"

# CHR-148 male
nest.all2$dad[which(nest.all2$mom=="2870-64909")] <- "2850-57582"

# Make Believe-62 males
nest.all2$dad[which(nest.all2$mom=="2850-57876" & 
                      nest.all2$attempt==1)] <- "2850-57680"
nest.all2$dad[which(nest.all2$mom=="2850-57876" & 
                      nest.all2$attempt==2)] <- "2850-57750"



### Identify nests with overlap of fert and active------------------------------
#-------------------------------------------------------------------------------

## Temporal overlap should occur if active end is not > fert start OR
# if active start is not > fert end.
# In other words, active start OR active end needs to be between fert start and end
# OR active start < fert start & active end > fert end

## Make a matrix for all the nests

# first make unique ID for each nest with BarnNest and band
nest.all2$matID <- paste(nest.all2$BarnNest, nest.all2$mom, 
                        nest.all2$dad, nest.all2$attempt, sep="_")

# make matrix for duration of overlap
mat.temporal.dur <- matrix(NA, nrow=120, ncol=120)

# save this table
write.csv(nest.all2, "output-files/nest all table.csv", row.names = F)


## loop to identify which nests had temporal overlap----------------------------


# for each nest in nest.all2 (i is focal female ID)
for (i in 1:120) {
  # loop through the other nests (n is male ID) and check whether there 
  # is temporal overlap
  for (n in 1:120) {
    
    # check if active start is in fertile period
    start.in <- nest.all2$active.start[n] >= nest.all2$fert.start[i] & 
      nest.all2$active.start[n] <= nest.all2$fert.end[i]
    # check if active end is in fertile period
    end.in <- nest.all2$active.end[n] >= nest.all2$fert.start[i] & 
      nest.all2$active.end[n] <= nest.all2$fert.end[i]
    # check if fertile period is inside active period
    fert.in <- nest.all2$active.start[n] <= nest.all2$fert.start[i] &
      nest.all2$active.end[n] >= nest.all2$fert.end[i]
    
    ## calculate duration based on type of overlap
    # add 1 to account for same day overlap 
    # (ex: 2025-06-12 minus 2025-06-12 = 0, but should be 1 day overlap)
    
    # active start in fertile period
    if (start.in==T) {
      dur <- nest.all2$fert.end[i] - nest.all2$active.start[n] +1
    } else {
      # active end in fertile period
      if (end.in == T) {
        dur <- nest.all2$active.end[n] - nest.all2$fert.start[i] + 1
      } else {
        # fertile period inside active period
        if (fert.in==T) {
          dur <- nest.all2$fert.end[i] - nest.all2$fert.start[i] + 1
        } else {
          dur <- 0
        }
        
      }
    }
    
    # fill dur matrix with value
    mat.temporal.dur[i, n] <- dur
    
  }
}

# try converting to dataframe first
mat.temp.df <- as.data.frame(mat.temporal.dur)
colnames(mat.temp.df) <- nest.all2$matID
rownames(mat.temp.df) <- nest.all2$matID
mat.temp.df$nest.mom <- row.names(mat.temp.df)

# try using pivot_longer
temporal.el <- pivot_longer(mat.temp.df, !nest.mom , names_to = c("nest.dad"), 
                          values_to = "dur")
# add column for binary overlap
temporal.el$overlap <- ifelse(temporal.el$dur > 0, 1, 0)


### Another loop to calculate overlap in fertile periods (rather than active period)

# make matrix for duration of overlap
mat.fertile.dur <- matrix(NA, nrow=120, ncol=120)

# for each nest in nest.all2 (i is focal female ID)
for (i in 1:120) {
  # loop through the other nests (n is male ID) and check whether there 
  # is temporal overlap
  for (n in 1:120) {
    
    # check if fert1 start is in fert2 period
    start.in <- nest.all2$fert.start[n] >= nest.all2$fert.start[i] & 
      nest.all2$fert.start[n] <= nest.all2$fert.end[i]
    # check if fert1 end is in fert2 period
    end.in <- nest.all2$fert.end[n] >= nest.all2$fert.start[i] & 
      nest.all2$fert.end[n] <= nest.all2$fert.end[i]
    # check if fert1 period is inside fert2 period
    fert.in <- nest.all2$fert.start[n] <= nest.all2$fert.start[i] &
      nest.all2$fert.end[n] >= nest.all2$fert.end[i]
    
    ## calculate duration based on type of overlap
    # add 1 to account for same day overlap 
    # (ex: 2025-06-12 minus 2025-06-12 = 0, but should be 1 day overlap)
    
    # fert start in fert period
    if (start.in==T) {
      dur <- nest.all2$fert.end[i] - nest.all2$fert.start[n] +1
    } else {
      # fert end in fert period
      if (end.in == T) {
        dur <- nest.all2$fert.end[n] - nest.all2$fert.start[i] + 1
      } else {
        # fert period inside fert period
        if (fert.in==T) {
          dur <- nest.all2$fert.end[i] - nest.all2$fert.start[i] + 1
        } else {
          dur <- 0
        }
        
      }
    }
    
    # fill dur matrix with value
    mat.fertile.dur[i, n] <- dur
    
  }
}

# try converting to dataframe first
mat.fert.df <- as.data.frame(mat.fertile.dur)
colnames(mat.fert.df) <- nest.all2$matID
rownames(mat.fert.df) <- nest.all2$matID
mat.fert.df$nest.mom <- row.names(mat.fert.df)

# try using pivot_longer
fertile.el <- pivot_longer(mat.fert.df, !nest.mom , names_to = c("nest.dad"), 
                            values_to = "dur.fertile")
# add column for binary overlap
fertile.el$overlap.fertile <- ifelse(fertile.el$dur.fertile > 0, 1, 0)



## Second loop to calculate differences in CI dates-----------------------------

# make matrix for duration of overlap
mat.ci.diff <- matrix(NA, nrow=120, ncol=120)

# for each nest in nest.all2 (i is focal female ID)
for (i in 1:120) {
  # loop through the other nests (n is male ID) and check whether there 
  # is temporal overlap
  for (n in 1:120) {
    # value should be negative if male CI was sooner than female CI
    # subtract female CI from male CI
    ci.diff <- nest.all2$ci[n] - nest.all2$ci[i] 
    
    # fill diff matrix with value
    mat.ci.diff[i, n] <- ci.diff
  }
}

# try converting to dataframe first
mat.ci.df <- as.data.frame(mat.ci.diff)
colnames(mat.ci.df) <- nest.all2$matID
rownames(mat.ci.df) <- nest.all2$matID
mat.ci.df$nest.mom <- row.names(mat.ci.df)

# try using pivot_longer
ci.el <- pivot_longer(mat.ci.df, !nest.mom , names_to = c("nest.dad"), 
                            values_to = "ci.diff")


## Loop to calculate rest of barn synchrony index for CHR only -----------------

# subset only CHR nests from dyad table
nest.chr <- subset(nest.all2, nest.all2$site=="CHR")

# subset only CHR birds from matrix of fertile days
mat.chr <- subset(mat, mat$site=="CHR" | is.na(mat$site))

# convert julian days to yyyy-mm-dd in matrix of fertile days
date.format <- parse_date_time(x = paste(2022, as.character(mat.chr[1, 4:104])),
                               orders = "yj")
# add to mat.chr
mat.chr2 <- as.data.frame(mat.chr[-1, ])
colnames(mat.chr2)[4:104] <- as.character(date.format)
# convert zeros and ones to numeric
mat.num <- as.matrix(mat.chr2[,4:104])
# add labels back
mat.num2 <- cbind(mat.chr2[,1:3], mat.num)

# check that colSums works
colSums(mat.num2[,4:104])


# add number of fertile days to dyad table
# Not using fertile days from mat because that includes multiple nesting 
# attempts for each female. 
nest.chr$fert.start <- ymd(nest.chr$fert.start)
nest.chr$fert.end <- ymd(nest.chr$fert.end)
nest.chr$num.fert.days <- nest.chr$fert.end - nest.chr$fert.start


## Loop for calculating rest.SI for each dyad

# This is based off of Kempanears 1993 synchrony index, and uses the inner
# term of that metric. 0% means completely asynchronous, and 100% means 
# completely synchronous
# For each dyad, of the days that the focal female is fertile, count up the number of
# other females that were also fertile that day excluding the other female in
# that dyad. Sum over all of the fertile days for focal female p
# divide by: the number of fertile days for female p times the total number of
# other females excluding both females in the dyad. 

# Within each female, the rest.SI value will only vary by dyad if the dyad partner
# is fertile at the same time. Even then, it should't change the rest.SI 
# value very much

# make matrix for SI
mat.rest.SI <- matrix(NA, nrow=59, ncol=59)

# for each nest in nest.all2 (i is focal female ID)
for (i in 1:59) {
  # loop through the other nests (n is dyad partner) and calculate overlap SI 
  # for focal female excluding dyad partner
  for (n in 1:59) {
    
    # fertile dates for focal female
    dates <- mat.num[ ,which(colnames(mat.num)==nest.chr$fert.start[i]) : 
                        which(colnames(mat.num)==nest.chr$fert.end[i])]
    # row of focal female in mat.num
    focal.row <- which(mat.num2$band == nest.chr$mom[i])
    # row of dyad partner in mat.num
    dyad.row <- which(mat.num2$band == nest.chr$mom[n])
    # number of other females fertile on those days
    # exclude focal female and dyad partner rows
    other.fem.by.day <- colSums(dates[-c(focal.row, dyad.row),])
    # sum across all fertile days
    other.fem.tot <- sum(other.fem.by.day)
    # divide by max possible synchrony (all other females fertile on all days)
    rest.SI <- other.fem.tot/(length(mat.num2$band)*as.numeric(nest.chr$num.fert.days[i]))
    # save result
    mat.rest.SI[i,n] <- rest.SI
  }
}

# convert to dataframe first
rest.SI.df <- as.data.frame(mat.rest.SI)
colnames(rest.SI.df) <- nest.chr$matID
rownames(rest.SI.df) <- nest.chr$matID
rest.SI.df$nest.mom <- row.names(rest.SI.df)

# pivot_longer to convert to an edge list
rest.SI.el <- pivot_longer(rest.SI.df, !nest.mom , names_to = c("nest.dad"), 
                           values_to = "rest.SI")



# add columns for nest and band-------------------------------------------------

# combine tables
temporal.el2 <- left_join(temporal.el, fertile.el, 
                          by=c("nest.mom","nest.dad")) %>%
  left_join(ci.el, by=c("nest.mom","nest.dad")) %>%
  left_join(rest.SI.el, by=c("nest.mom","nest.dad"))

# split up labels columns
temporal.el.ids <- temporal.el2 %>%
  separate(nest.mom, c("location1", "mom1", "dad1", "attempt1"), "_",
           remove = FALSE) %>%
  separate(nest.dad, c("location2", "mom2", "dad2", "attempt2"), "_",
           remove = FALSE) %>%
  select(location1, location2, mom1, dad2, attempt1, attempt2,
         dur, overlap, dur.fertile, overlap.fertile, ci.diff, rest.SI)

temporal.el.ids$pairwise <- paste(temporal.el.ids$mom1, temporal.el.ids$dad2,
                                  sep="_")


## check pairwise and repeated nest locations by pulling out example
example <- filter(
  temporal.el.ids,
  location1=="CHR-27" & location2=="Make Believe-73" |
    location1=="CHR-29" & location2=="Make Believe-73" |
    location1=="Make Believe-73" & location2=="CHR-27" |
    location1=="Make Believe-73" & location2=="CHR-29")


## Merge with nest distance table-----------------------------------------------

# harmonize names
colnames(dist)[6:7] <- c("location1","location2")

# make join columns
dist$locations <- paste(dist$location1, dist$location2, sep="_")
temporal.el.ids$locations <- paste(temporal.el.ids$location1, 
                                   temporal.el.ids$location2, sep="_")
# join tables
# note that dist only includes distances for nests within the same barn
el.time.dist <- left_join(temporal.el.ids, dist[,c(1,4,5,10)],by="locations")

# for social pairs, add distance of 0.08 (8cm)
el.time.dist$distance.m[which(el.time.dist$location1==el.time.dist$location2)] <- 0.08

# remove NA rows for distance
el.time.dist.sameBarn <- subset(el.time.dist, !is.na(el.time.dist$distance.m))

# including zeros for dur means we have social pairs being compared to their own 
# different nesting attempts. 
# example: 2850-57586_2850-57669. Fix this by adding a column indicating social pairs
# at the pairwise level

# add column to indicate social pairs based on pairwise label
# pull out set of pairwise where location1=location2
soc.pair <- subset(el.time.dist.sameBarn, el.time.dist.sameBarn$location1==
                     el.time.dist.sameBarn$location2)
# add column
el.time.dist.sameBarn$soc.pair <- ifelse(el.time.dist.sameBarn$pairwise %in% 
                                           soc.pair$pairwise, 1, 0)

# set all social pair distances to 0.08, to address cases where pairs are being
# compared to their own previous or later nesting attempts in different nest 
# locations
el.dur <- el.time.dist.sameBarn
el.dur$distance.m[which(el.dur$soc.pair==1)] <- 0.08

# summarize by pairwise and attempt1 (from mom nest) to check number of rows
el.attempt.summary <- el.dur %>%
  group_by(pairwise, attempt1) %>%
  summarise(replicates=n(),
            max.dist=max(distance.m),
            min.dist=min(distance.m),
            diff=max.dist - min.dist,
            ci.diff.min=min(ci.diff))

# explore cases where replicates > 2. Main issue is NA dads and 2850-57806 male
# that was attending 2 nests simultaneously. 

# CI for CHR 96 is 2022-05-18 and 2022-07-10. CI for nest 37 is 2022-05-18
# and 2022-07-09. Will set nest 96 as location for brood1 and 37 as location
# for brood2 for male 2850-57806

# Remove rows for 2850-57806 with location 37 and attempt 1, and 
# with location 96 and attempt 2

# Male 2850-57582 also attended nest 6 and nest 148 simultaneously, but nest 6
# initiated first, so make that his location and remove 148 from location 2

# identify rows to remove
el.dur$remove <- ifelse(
  el.dur$location2=="CHR-37" & el.dur$attempt2==1 |
    el.dur$location2=="CHR-96" & el.dur$attempt2==2 |
    el.dur$location2=="CHR-148", 1, 0)

# remove the rows
el.dur1.2 <- subset(el.dur, el.dur$remove==0)

# make summary again
el.attempt.summary2 <- el.dur1.2 %>%
  group_by(pairwise, attempt1) %>%
  summarise(replicates=n(),
            max.dist=max(distance.m),
            min.dist=min(distance.m),
            diff=max.dist - min.dist, 
            ci.diff.min=min(ci.diff))

# Check rows for replicates=3
replicates3 <- filter(el.attempt.summary2, replicates==3)
# split out mom and dad
replicates3.id <- replicates3 %>%
  separate(pairwise, c("mom.band","dad.band"), sep="_")
# Three problematic dads to investigate: 
# 2850-57633, Cooks male with 3 nesting attempts
# 2850-57669, CHR male with 3 attempts
# 2850-57731, another CHR male with 3 attempts
# These all seem real so can stay

#-------------------------------------------------------------------------------
# Remove "fringe" nests in different structures
#-------------------------------------------------------------------------------

# Some sites had additional out buildings where birds nested, but which were 
# physically separate from the main nest colony building. We will exclude these
# for the analysis of fine-scale spatial effects.

# Most inclusive would be to only exclude BlueCloud-24, which is >100m from other nests
# Logic is that satellite birds were mostly observed on wires with birds from main 
# structure, except in the case of BC-24. Max dist is then 71m

# Most restrictive would be to only include nests within the same continuous structure,
# so exclude BlueCloud-24, 19, 26 and 27, plus Struthers-22 and 23, and MakeBelieve-62
# This reduces the range to max of ~31m. Histogram supports this strategy.

el.dur2 <- filter(el.dur1.2, 
                  location1 != "Blue Cloud-24" & location2 != "Blue Cloud-24" &
                    location1 != "Blue Cloud-19" & location2 != "Blue Cloud-19" & 
                    location1 != "Blue Cloud-26" & location2 != "Blue Cloud-26" &
                    location1 != "Blue Cloud-27" & location2 != "Blue Cloud-27" &
                    location1 != "Struthers-22" & location2 != "Struthers-22" &
                    location1 != "Struthers-23" & location2 != "Struthers-23" &
                    location1 != "Make Believe-62" & location2 != "Make Believe-62")

#-------------------------------------------------------------------------------
# Summarize so one distance value for each pairwise - for dur.active first
#-------------------------------------------------------------------------------

## For positive dur, can use time weighted distances----------------------------

el.pos.dur2 <- subset(el.dur2, el.dur2$dur>0)

# add column for distance times dur
el.pos.dur2$distXdur <- el.pos.dur2$distance.m *
  el.pos.dur2$dur

## summarize and add mean daily distance, weighted by duration of overlap

# for each pairwise, summarize number of replicates and max and min distance

el.attempt.summary3 <- el.pos.dur2 %>%
  group_by(pairwise, attempt1, location1, soc.pair) %>%
  summarise(replicates=n(),
            max.dist=max(distance.m),
            min.dist=min(distance.m),
            diff=max.dist - min.dist,
            max.dur=max(dur),
            min.dur=min(dur),
            dur.diff=max.dur-min.dur,
            total.dur=sum(dur),
            total.weight.dist=sum(distXdur),
            mean.daily.dist=total.weight.dist/total.dur,
            diff.min.weight=round( (mean.daily.dist - min.dist), digits=2) )


#-------------------------------------------------------------------------------
## For zero dur, need to choose active nest that is closest in time to focal

el.zero.dur <- subset(el.dur2, el.dur2$dur==0)

# check range for min ci diff
# smallest absolute value is positive 9

el.zero.summary <- el.zero.dur %>%
  group_by(pairwise, attempt1, location1, soc.pair) %>%
  summarise(replicates=n(),
            max.dist=max(distance.m),
            min.dist=min(distance.m),
            diff=max.dist - min.dist,
            min.ci.diff=min(ci.diff),
            max.ci.diff=max(ci.diff),
            total.dur=sum(dur))

# check cases where min is negative and max is positive
zero.min.max <- subset(el.zero.summary, el.zero.summary$min.ci.diff<0 &
                         el.zero.summary$max.ci.diff>0)
# for the cases where there is a difference in spatial distance, the positive CI
# difference has a much smaller absolute value, at least 2-3 times smaller than
# the negative CI difference. Good to go with the distance that corresponds to
# the smallest absolute value difference in CI.

# add column for ci.diff absolute value
el.zero.dur$ci.diff.abs <- abs(el.zero.dur$ci.diff)

# loop to identify smallest ci.diff for each attempt and take that distance
el.zero.summary$distance.m <- NA

for (i in 1:length(el.zero.summary$pairwise)) {
  
  pair <- el.zero.summary$pairwise[i]
  attempt <- el.zero.summary$attempt1[i]
  
  table <- subset(el.zero.dur, el.zero.dur$pairwise==pair &
                    el.zero.dur$attempt1==attempt)
  
  min.ci.diff <- table[which(table$ci.diff.abs==min(table$ci.diff.abs)), ]
  
  distance <- min.ci.diff$distance.m
  
  el.zero.summary$distance.m[i] <- distance
}

# Combine pos and zero dur tables back together---------------------------------

# rename distance columns
colnames(el.zero.summary)[12] <- "process.dist.m"
colnames(el.attempt.summary3)[14] <- "process.dist.m"

# merge, only keeping relevant columns
el.attempt.dist1 <- rbind(select(el.attempt.summary3, pairwise, attempt1,location1,
                                soc.pair, total.dur, process.dist.m),
                         select(el.zero.summary, pairwise, attempt1, location1, 
                                soc.pair,total.dur, process.dist.m))

ggplot(el.attempt.dist1, aes(x=process.dist.m)) + geom_histogram(binwidth = 5)

# check that there are not duplicates within pairwise and attempt1
dup.check <- el.attempt.dist1 %>% group_by(pairwise, attempt1) %>%
  summarise(n=n(),
            max.dur=max(total.dur))

# there are no cases where rep=2 but max dur is zero
which(dup.check$n==2 & dup.check$max.dur==0)

# most of the pairwise-attempt have 2 replicates instead of 1
# When this happens, keep the row with the non-zero dur

list.pairwise <- unique(el.attempt.dist1$pairwise)

storage <- as.data.frame(matrix(nrow=1, ncol=6, NA))
colnames(storage) <- colnames(el.attempt.dist1)

for (i in 1:length(list.pairwise)) {
  pair <- subset(el.attempt.dist1, el.attempt.dist1$pairwise==list.pairwise[i])
  attempts <- unique(pair$attempt1)
  for (j in 1:length(attempts)) {
    table <- subset(pair, pair$attempt1==attempts[j])
    keep <- table[which(table$total.dur==max(table$total.dur)), ]
    storage <- rbind(storage, keep)
  }
}

# remove NA row
el.attempt.dist2 <- storage[-1, ]

#-------------------------------------------------------------------------------
# Summarize so one distance value for each pairwise - for dur.fertile
#-------------------------------------------------------------------------------

# check replicates and range of fert.dur and rest.SI
el.fert.summary <- el.dur2 %>%
  group_by(pairwise, attempt1, location1, soc.pair) %>%
  summarise(replicates=n(),
            max.dist=max(distance.m),
            min.dist=min(distance.m),
            max.dur.fert=max(dur.fertile),
            min.dur.fert=min(dur.fertile),
            max.rest.SI = max(rest.SI, na.rm=T),
            min.rest.SI = min(rest.SI, na.rm=T))

# any non-zero min.dur.fert cases have only one replicate and min=max
# for each pairwise-attempt keep row with non-zero dur
# if multiple rows with zero dur, keep one with closest CI diff

list.pairwise.fertile <- unique(el.fert.summary$pairwise)

storage <- as.data.frame(matrix(nrow=1, ncol=19, NA))
colnames(storage) <- colnames(el.dur2)

for (i in 1:length(list.pairwise)) {
  pair <- subset(el.dur2, el.dur2$pairwise==list.pairwise[i])
  attempts <- unique(pair$attempt1)
  for (j in 1:length(attempts)) {
    table <- subset(pair, pair$attempt1==attempts[j])
    keep <- table[which(table$dur.fertile==max(table$dur.fertile)), ]
    keep2 <- keep[which(keep$ci.diff==min(abs(keep$ci.diff)) |
                    keep$ci.diff==-min(abs(keep$ci.diff))), ]
    storage <- rbind(storage, keep2)
  }
}

# remove NA row
el.fert.dist <- storage[-1, ]

el.fert.summary2 <- el.fert.dist %>%
  group_by(pairwise, attempt1, location1, soc.pair) %>%
  summarise(replicates=n(),
            max.dist=max(distance.m),
            min.dist=min(distance.m),
            max.dur.fert=max(dur.fertile),
            min.dur.fert=min(dur.fertile),
            early.ci.diff=min(ci.diff),
            close.ci.diff=min(abs(ci.diff)),
            diff.ci.diff=abs(early.ci.diff) - abs(close.ci.diff),
            max.rest.SI = max(rest.SI),
            min.rest.SI = min(rest.SI))

# check for cases where max and min are different for dur.fert
# all are the same, so just keep max
el.fert.summary2$diff <- el.fert.summary2$max.dur.fert==el.fert.summary2$min.dur.fert

# check cases where max and min are different for rest.SI
# all are the same
el.fert.summary2$diff.rest.SI <- el.fert.summary2$max.rest.SI==el.fert.summary2$min.rest.SI

# combine el.attempt.dist with dur.fert
# close ci and early ci are the same, so just keep early
el.attempt.dist <- left_join(el.attempt.dist2, el.fert.summary2[, c(1:4,8,10, 13)], 
                             by=c("pairwise","attempt1","location1","soc.pair"))
colnames(el.attempt.dist)[7:9] <- c("dur.fertile","ci.diff","rest.SI.chr")


#-------------------------------------------------------------------------------
# Add paternity data to filter female clutches that were not analyzed for paternity
#-------------------------------------------------------------------------------

# Can't just use genetic_fam colunm in kin.ci.nest to merge, because this would 
# lose dyads that are valid but didn't share any offspring. Also need to make
# sure that the clutch attempt labels merge correctly. 

# add BarnNest column to kin.ci.nest
kin.ci.nest$BarnNest <- paste(kin.ci.nest$Site_ind2, round(kin.ci.nest$nest, digits=0), sep="-")

# identify nests that are in kin.ci.nest
nest.kin.index <- el.attempt.dist$location1 %in% kin.ci.nest$BarnNest

# subset only indexed nests
# Note, can't have location2 because sometimes for a single female clutch, the
# focal male was in more than one location. Hence the mean daily dist calculation.
el.attempt.dist2 <- el.attempt.dist[which(nest.kin.index==T), ]

ggplot(el.attempt.dist2, aes(x=process.dist.m)) + geom_histogram(binwidth = 2)


## Check that females within the same site have the same number of male comparisons
el.attempt.dist3 <- el.attempt.dist2 %>%
  separate_wider_delim(pairwise, names=c("female","male"), delim="_", cols_remove=F) %>%
  separate_wider_delim(location1, names=c("barn","nest"), delim="-", cols_remove = F)

fem.check <- el.attempt.dist3 %>%
  group_by(barn, female, attempt1) %>%
  summarise(n.male = length(unique(male)))

male.check <- el.attempt.dist3 %>%
  group_by(barn, male) %>%
  summarise(n.female = length(unique(female)))

#-------------------------------------------------------------------------------
# Add shared fertilizations at the clutch level
#-------------------------------------------------------------------------------

# add clutch ID to el.attempt for merging with fert
el.attempt.dist3$clutch_id_merge <- paste(el.attempt.dist3$barn, 
                                         el.attempt.dist3$nest,
                                         el.attempt.dist3$attempt1, sep="_")
# correct Mary Ann's_1_1 should be changed to Mary Ann's_1_2
el.attempt.dist3$clutch_id_merge[which(
  el.attempt.dist3$clutch_id_merge=="Mary Ann's_1_1")] <- "Mary Ann's_1_2"
# Cathy's_2_3 should be changed to 2_2 in attempt.dist
el.attempt.dist3$clutch_id_merge[which(
  el.attempt.dist3$clutch_id_merge=="Cathy's_2_3")] <- "Cathy's_2_2"
# CHR_148_1 should be changed to CHR_148_2 in attempt.dist because of timing
el.attempt.dist3$clutch_id_merge[which(
  el.attempt.dist3$clutch_id_merge=="CHR_148_1")] <- "CHR_148_2"

# add rounded nest column to fert
kin.ci.nest$nest_round <- round(kin.ci.nest$nest)

# add new clutch ID to fert
kin.ci.nest$clutch_id_merge <- paste(kin.ci.nest$Site_ind2,
                                     kin.ci.nest$nest_round,
                                     kin.ci.nest$brood, sep="_")

# calculate total offspring analyzed for paternity for each clutch
clutch.size <- kin.ci.nest %>% group_by(clutch_id_merge) %>%
  summarise(pat.clutch.size = n())
# Speedwell has strange duplicated columns, and I'm not sure why. So does
# CHR_101_2, and Marte's 1 and Marte's 2

# remove duplicated rows
kin.ci.nest2 <- kin.ci.nest[-which(kin.ci.nest$notes=="likely 1st clutch in same nest as 2nd brood but not monitored until later in the season. Saw 4 fledglings flying around during adult captures"), ]

kin.ci.nest3 <- kin.ci.nest2[-which(kin.ci.nest2$FamilyID_ind2=="Speedwell-01" &
                                     kin.ci.nest2$ci_dad=="unknown"), ]

kin.ci.nest4 <- kin.ci.nest3[-which(kin.ci.nest3$clutch_id_ind2=="CHR_101_2" &
                                      kin.ci.nest3$ci_dad=="2022-08-10"), ]

kin.ci.nest5 <- kin.ci.nest4[-c(59,61,63,280,282,330), ]

clutch.size2 <- kin.ci.nest5 %>% group_by(clutch_id_merge) %>%
  summarise(pat.clutch.size = n())
# Marte's 1 is still wrong and I don't know why the specified rows weren't removed!
# Carry on because Marte's will be taken out anyways

# add clutch size to kin.ci.nest5
kin.ci.nest6 <- left_join(kin.ci.nest5, clutch.size2, by="clutch_id_merge")


# summarise fert by clutch ID
kin.clutch <- kin.ci.nest6 %>%
  group_by(clutch_id_merge, genetic_fam, fert_type, pat.clutch.size) %>%
  summarise(shared_fert = n())

# re-order genetic fam to be mom_dad instead of dad_mom
kin.clutch2 <- kin.clutch %>%
  separate_wider_delim(genetic_fam, names=c("dad","mom"), delim="_", cols_remove = F)

# Correct CHR_101_NA, with mom band still left as NA, should be 2870-64908
kin.clutch2$mom[which(kin.clutch2$mom=="NA" & 
                        kin.clutch2$clutch_id_merge=="CHR_101_NA")] <- "2870-64908"
kin.clutch2$clutch_id_merge[which(kin.clutch2$clutch_id_merge=="CHR_101_NA")] <- "CHR_101_2"

# make reordered pairwise column
kin.clutch2$pairwise <- paste(kin.clutch2$mom, kin.clutch2$dad, sep="_")

# first check overlap by clutch ID
el.attempt.dist3$clutch_match <- ifelse(el.attempt.dist3$clutch_id_merge %in% 
                                          kin.clutch2$clutch_id_merge, 1, 0)

# summarise clutch IDs that are not in fert
clutch.no.fert <- el.attempt.dist3 %>%
  filter(clutch_match==0)
clutch.no.fert2 <- sort(unique(clutch.no.fert$clutch_id_merge))
# view list of 13 clutches
clutch.no.fert2
## manually check for name typos
# CHR_78_2 is not in the fert, unsure why because nest fledged and nestlings were banded
# Karens_1_1 nest failed before nestlings reached D12
# Make Beleive_8_2 never hatched

# keep only clutches that are matched in fert
el.attempt.dist4 <- filter(el.attempt.dist3, clutch_match==1)

# add shared ferts by merging tables
el.attempt.dist5 <- left_join(el.attempt.dist4, select(kin.clutch2[,-4], pairwise, 
                                                       clutch_id_merge, fert_type,
                                                       shared_fert),
                              by=c("clutch_id_merge", "pairwise"))

# add clutch size by clutch_id but not pairwise. All observations should have
# clutch size, not just those that share ferts
el.attempt.dist5.2 <- left_join(el.attempt.dist5, clutch.size2, 
                                by="clutch_id_merge")

# fill in NAs with zeros
el.attempt.dist5.2$shared_fert[which(is.na(el.attempt.dist5.2$shared_fert))] <- 0

# add column for binary shared fert (0=no, 1=yes)
el.attempt.dist5.2$bin_fert <- ifelse(el.attempt.dist5.2$shared_fert>0, "yes", "no")

# remove cases of social pairs where dur=0, because these are comparing a pair's 
# current attempt to later or earlier attempts which doesn't make sense
el.attempt.dist5.2$remove <- ifelse(el.attempt.dist5.2$soc.pair==1 & 
                                    el.attempt.dist5.2$total.dur==0, 1, 0)

el.attempt.dist6 <- filter(el.attempt.dist5.2, remove==0)

## remove solitary sites: Dome House, Mary Ann's, Karen's, Cathy's
el.attempt.dist7 <- filter(el.attempt.dist6, barn!="Dome House" &
                              barn!="Mary Ann's" & barn!="Karen's" &
                              barn!="Cathy's")

# remove remove column
el.attempt.dist8 <- select(el.attempt.dist7, -remove)



# save final table
write.csv(el.attempt.dist8, 
          "output-files/edge list distance and fert_with fertile and rest SI.csv", 
          row.names = F)




