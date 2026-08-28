
################################################################################
# Summary stats for CO Fert Net 2022, CHR only
# Heather Kenny-Duddela
# June 3, 2025
# updated Aug 18, 2026
################################################################################

# libraries
library(dplyr)



# load data

# from script adjust-table-CHRonly.R
dat <- read.csv("input-files/dyad_ep_complete_cases_CHRonly_with fertile and rest SI.csv")

pheno <- read.csv("input-files/pheno-CHR-for-unscaling.csv")


#-------------------------------------------------------------------------------
# summary stat for all CHR birds (used for scaling)
#-------------------------------------------------------------------------------

# remove "F?"
pheno$sex[which(pheno$sex == "F?")] <- "F"

pheno.num <- select(pheno, band, sex, TS, TS_scaled,
                    throat_avg_bright, throat_avg_bright_scaled,
                    breast_avg_bright, breast_avg_bright_scaled)

## Females

pheno.fem <- filter(pheno.num, sex=="F") %>%
  select(-band, -sex)

# make summary stat table (empty)
summary.stats.fem.all <- as.data.frame(matrix(NA, nrow=5, ncol=6))
colnames(summary.stats.fem.all) <- colnames(pheno.fem)
# sample size
summary.stats.fem.all[1,] <- colSums(!is.na(pheno.fem))
# mean
summary.stats.fem.all[2,] <- apply(pheno.fem, 2 ,mean,na.rm=T)
# sd
summary.stats.fem.all[3,] <- apply(pheno.fem, 2, sd, na.rm=T)
# min
summary.stats.fem.all[4,] <- apply(pheno.fem, 2, min, na.rm=T)
# max
summary.stats.fem.all[5,] <- apply(pheno.fem, 2, max, na.rm=T)
# transpose
t.summary.stats.fem.all <- as.data.frame(t(as.matrix(summary.stats.fem.all)))
# add column labels
colnames(t.summary.stats.fem.all) <- c("N","mean","sd","min","max")



## males

pheno.male <- filter(pheno.num, sex=="M") %>%
  select(-band, -sex)

# make summary stat table (empty)
summary.stats.male.all <- as.data.frame(matrix(NA, nrow=5, ncol=6))
colnames(summary.stats.male.all) <- colnames(pheno.male)
# sample size
summary.stats.male.all[1,] <- colSums(!is.na(pheno.male))
# mean
summary.stats.male.all[2,] <- apply(pheno.male, 2 ,mean,na.rm=T)
# sd
summary.stats.male.all[3,] <- apply(pheno.male, 2, sd, na.rm=T)
# min
summary.stats.male.all[4,] <- apply(pheno.male, 2, min, na.rm=T)
# max
summary.stats.male.all[5,] <- apply(pheno.male, 2, max, na.rm=T)
# transpose
t.summary.stats.male.all <- as.data.frame(t(as.matrix(summary.stats.male.all)))
# add column labels
colnames(t.summary.stats.male.all) <- c("N","mean","sd","min","max")


#-------------------------------------------------------------------------------
# summary stats at dyad level
#-------------------------------------------------------------------------------


# filter for numeric columns
numeric <- dat[ ,c(4, 10, 11, 13, 17,18, 42)]
# binary attempt
numeric$brood2 <- ifelse(numeric$attempt1 == 2, 1, 0)

# make summary stat table (empty)
summary.stats.dyad <- as.data.frame(matrix(NA, nrow=3, ncol=8))
colnames(summary.stats.dyad) <- colnames(numeric)
# sample size
summary.stats.dyad[1,] <- colSums(!is.na(numeric))
# mean
summary.stats.dyad[2,] <- apply(numeric, 2 ,mean,na.rm=T)
# sd
summary.stats.dyad[3,] <- apply(numeric, 2, sd, na.rm=T)
# transpose
t.summary.stats.dyad <- as.data.frame(t(as.matrix(summary.stats.dyad)))
# add column labels
colnames(t.summary.stats.dyad) <- c("N","mean","sd")
# remove attempt1
summary.stats.dyad.final <- t.summary.stats.dyad[-1, ]
# export table
write.csv(summary.stats.dyad.final, "output-files/summary stats local dyad.csv")


#-------------------------------------------------------------------------------
# Make table of complete cases used for modeling
#-------------------------------------------------------------------------------

# unscale social partner phenotypes
# females
dat$soc.fem_breast_avg_bright_unscale <- dat$soc.fem_breast_avg_bright_scaled*t.summary.stats.fem.all$sd[5] + t.summary.stats.fem.all$mean[5]

range(dat$soc.fem_breast_avg_bright_unscale)

dat$soc.fem_throat_avg_bright_unscale <- dat$soc.fem_throat_avg_bright_scaled*t.summary.stats.fem.all$sd[3] + t.summary.stats.fem.all$mean[3]

range(dat$soc.fem_throat_avg_bright_unscale)

dat$soc.fem_TS_unscale <- dat$soc.fem_TS_scaled*t.summary.stats.fem.all$sd[1] + t.summary.stats.fem.all$mean[1]

range(dat$soc.fem_TS_unscale)

# males
dat$soc.male_breast_avg_bright_unscale <- dat$soc.male_breast_avg_bright_scaled*t.summary.stats.male.all$sd[5] + t.summary.stats.male.all$mean[5]

range(dat$soc.male_breast_avg_bright_unscale)

dat$soc.male_throat_avg_bright_unscale <- dat$soc.male_throat_avg_bright_scaled*t.summary.stats.male.all$sd[3] + t.summary.stats.male.all$mean[3]

range(dat$soc.male_throat_avg_bright_unscale)

dat$soc.male_TS_unscale <- dat$soc.male_TS_scaled*t.summary.stats.male.all$sd[1] + t.summary.stats.male.all$mean[1]

range(dat$soc.male_TS_unscale)



# split females and males
# females
fem <- select(dat, female, fem_throat_avg_bright:fem_TS,  
              soc.male_TS_unscale:soc.male_breast_avg_bright_unscale)

fem2 <- fem[-which(duplicated(fem$female)), ]


# males
male <- select(dat, male, male_throat_avg_bright:male_TS, 
               soc.fem_TS_unscale:soc.fem_breast_avg_bright_unscale)

male2 <- male[-which(duplicated(male$male)), ]


#-------------------------------------------------------------------------------
# summary stats for each sex
#-------------------------------------------------------------------------------

## females

# filter for numeric columns
fem.numeric <- fem2[,-1 ]

# make summary stat table (empty)
summary.stats.fem <- as.data.frame(matrix(NA, nrow=3, ncol=6))
colnames(summary.stats.fem) <- colnames(fem.numeric)
# sample size
summary.stats.fem[1,] <- colSums(!is.na(fem.numeric))
# mean
summary.stats.fem[2,] <- apply(fem.numeric, 2 ,mean,na.rm=T)
# sd
summary.stats.fem[3,] <- apply(fem.numeric, 2, sd, na.rm=T)
# transpose
t.summary.stats.fem <- as.data.frame(t(as.matrix(summary.stats.fem)))
# add column labels
colnames(t.summary.stats.fem) <- c("N","mean","sd")


## males

# filter for numeric columns
male.numeric <- male2[,-1 ]

# make summary stat table (empty)
summary.stats.male <- as.data.frame(matrix(NA, nrow=3, ncol=6))
colnames(summary.stats.male) <- colnames(male.numeric)
# sample size
summary.stats.male[1,] <- colSums(!is.na(male.numeric))
# mean
summary.stats.male[2,] <- apply(male.numeric, 2 ,mean,na.rm=T)
# sd
summary.stats.male[3,] <- apply(male.numeric, 2, sd, na.rm=T)
# transpose
t.summary.stats.male <- as.data.frame(t(as.matrix(summary.stats.male)))
# add column labels
colnames(t.summary.stats.male) <- c("N","mean","sd")


#-------------------------------------------------------------------------------
# combine tables
#-------------------------------------------------------------------------------

# designate types
t.summary.stats.male.all$type <- "all males"
t.summary.stats.fem.all$type <- "all females"
t.summary.stats.fem$type <- "modeled females"
t.summary.stats.male$type <- "modeled males"

all <- rbind(t.summary.stats.male.all[,-c(4,5)], t.summary.stats.fem.all[,-c(4,5)],
             t.summary.stats.male, t.summary.stats.fem)

# save table
write.csv(all, "output-files/summary stats local node.csv")


