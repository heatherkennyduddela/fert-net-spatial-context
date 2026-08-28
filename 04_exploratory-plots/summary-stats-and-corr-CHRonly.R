
################################################################################
# Summary stats and correlations for ch 3 analysis
# Heather Kenny-Duddela
# May 1, 2025
################################################################################

# load data
dat <- read.csv("output-files/fert dist plumage-CHRonly_with fertile and rest SI.csv")

# load libraries
library(tidyverse)
library(Hmisc)

#-------------------------------------------------------------------------------
# Summary stats at dyad level
#-------------------------------------------------------------------------------

# filter for numeric columns
numeric <- dat[ ,c(10, 11, 13, 17,18)]
# add proportion fert
numeric$prop.fert <- numeric$shared_fert/numeric$pat.clutch.size

# make summary stat table (empty)
summary.stats.dyad <- as.data.frame(matrix(NA, nrow=3, ncol=6))
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
# export table
write.csv(t.summary.stats.dyad, "output-files/summary stats dyad.csv")


#-------------------------------------------------------------------------------
# Correlation between dist and dur and rest.SI
#-------------------------------------------------------------------------------

# dist and dur
ggplot(dat, aes(x=process.dist.m, y=total.dur)) + geom_point()

cor.test(dat$process.dist.m, dat$total.dur, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat$process.dist.m and dat$total.dur
# S = 197975589, p-value = 0.3979
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.02611347 


# dur and rest.SI
ggplot(dat, aes(x=rest.SI.chr)) + geom_histogram()

ggplot(dat, aes(x=rest.SI.chr, y=total.dur)) + geom_point()

cor.test(dat$rest.SI.chr, dat$total.dur, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat$rest.SI.chr and dat$total.dur
# S = 180555611, p-value = 0.0376
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.0641748


# rest.SI and dist
ggplot(dat, aes(x=rest.SI.chr, y=process.dist.m)) + geom_point()

cor.test(dat$rest.SI.chr, dat$process.dist.m, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat$rest.SI.chr and dat$process.dist.m
# S = 188978975, p-value = 0.5066
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.02051625



#-------------------------------------------------------------------------------
# Make social node table
#-------------------------------------------------------------------------------

# keep only social pairs
soc <- subset(dat, dat$soc.pair==1)

# keep only phenotype columns
soc.pheno <- select(soc, female, male, fem_mass_1:male_vent_avg_bright_scaled,
                    pat.clutch.size)

# keep only unique rows
pair.list <- which(duplicated(soc.pheno[,1:2]))
soc.node <- soc.pheno[-pair.list, ]

# save pair table
write.csv(soc.node, "output-files/pair pheno-CHRonly.csv", row.names=F)


#-------------------------------------------------------------------------------
# summary stats for soc.node
#-------------------------------------------------------------------------------


# filter for numeric columns
node.numeric <- soc.node[ ,c(4,8,)]
# add proportion fert
numeric$prop.fert <- numeric$shared_fert/numeric$pat.clutch.size

# make summary stat table (empty)
summary.stats.dyad <- as.data.frame(matrix(NA, nrow=3, ncol=6))
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
# export table
write.csv(t.summary.stats.dyad, "output-files/summary stats dyad.csv")


#-------------------------------------------------------------------------------
# Make correlation table for soc.node
#-------------------------------------------------------------------------------



### Females---------------------------------------------------------------------
fem <- select(soc.node, fem_mass_1, fem_throat_avg_bright, fem_breast_avg_bright,
              fem_belly_avg_bright, fem_vent_avg_bright, fem_TS, fem_mean_rwl,
              pat.clutch.size)

# calculate correlations
fem.table <- rcorr(as.matrix(fem), type="spearman")
# pull out correlation values
fem.R <- fem.table$r
# pull out p-vals
fem.p <- fem.table$P

## Define notions for significance levels; spacing is important.
mystars <- ifelse(fem.p < .0001, "****", 
                  ifelse(fem.p < .001, "*** ",
                         ifelse(fem.p < .01, "**  ",
                                ifelse(fem.p < .05, "*   ", "    "))))

## truncate the correlation matrix to two decimal
fem.R2 <- format(round(cbind(rep(-1.11, ncol(fem)), fem.R), 2))[,-1]

## build a new matrix that includes the correlations with their appropriate stars
Rnew <- matrix(paste(fem.R2, mystars, sep=""), ncol=ncol(fem))
diag(Rnew) <- paste(diag(fem.R), " ", sep="")
rownames(Rnew) <- colnames(fem)
colnames(Rnew) <- paste(colnames(fem), "", sep="")

# Hide upper triangle
fem.upper<-Rnew
fem.upper[upper.tri(Rnew)]<-""
fem.upper<-as.data.frame(fem.upper)

write.csv(fem.upper, "output-files/female pheno corr table-CHRonly.csv")


### Males---------------------------------------------------------------------
male <- select(soc.node, male_mass_1, male_throat_avg_bright, male_breast_avg_bright,
              male_belly_avg_bright, male_vent_avg_bright, male_TS, male_mean_rwl)

# calculate correlations
male.table <- rcorr(as.matrix(male), type="spearman")
# pull out correlation values
male.R <- male.table$r
# pull out p-vals
male.p <- male.table$P

## Define notions for significance levels; spacing is important.
mystars <- ifelse(male.p < .0001, "****", 
                  ifelse(male.p < .001, "*** ",
                         ifelse(male.p < .01, "**  ",
                                ifelse(male.p < .05, "*   ", "    "))))

## truncate the correlation matrix to two decimal
male.R2 <- format(round(cbind(rep(-1.11, ncol(male)), male.R), 2))[,-1]

## build a new matrix that includes the correlations with their appropriate stars
male.Rnew <- matrix(paste(male.R2, mystars, sep=""), ncol=ncol(male))
diag(male.Rnew) <- paste(diag(male.R), " ", sep="")
rownames(male.Rnew) <- colnames(male)
colnames(male.Rnew) <- paste(colnames(male), "", sep="")

# Hide upper triangle
male.upper<-male.Rnew
male.upper[upper.tri(male.Rnew)]<-""
male.upper<-as.data.frame(male.upper)

write.csv(male.upper, "output-files/male pheno corr table-CHRonly.csv")


### Between social pairs--------------------------------------------------------
pair <- select(soc.node, male_mass_1, male_throat_avg_bright, male_breast_avg_bright,
               male_belly_avg_bright, male_vent_avg_bright, male_TS, 
               male_mean_rwl,
               fem_mass_1, fem_throat_avg_bright, fem_breast_avg_bright,
               fem_belly_avg_bright, fem_vent_avg_bright, fem_TS, fem_mean_rwl,
               pat.clutch.size)

# calculate correlations
pair.table <- rcorr(as.matrix(pair), type="spearman")
# pull out correlation values
pair.R <- pair.table$r
# pull out p-vals
pair.p <- pair.table$P

## Define notions for significance levels; spacing is important.
mystars <- ifelse(pair.p < .0001, "****", 
                  ifelse(pair.p < .001, "*** ",
                         ifelse(pair.p < .01, "**  ",
                                ifelse(pair.p < .05, "*   ", "    "))))

## truncate the correlation matrix to two decimal
pair.R2 <- format(round(cbind(rep(-1.11, ncol(pair)), pair.R), 2))[,-1]

## build a new matrix that includes the correlations with their appropriate stars
pair.Rnew <- matrix(paste(pair.R2, mystars, sep=""), ncol=ncol(pair))
diag(pair.Rnew) <- paste(diag(pair.R), " ", sep="")
rownames(pair.Rnew) <- colnames(pair)
colnames(pair.Rnew) <- paste(colnames(pair), "", sep="")

# Hide upper triangle
pair.upper<-pair.Rnew
pair.upper[upper.tri(pair.Rnew)]<-""
pair.upper<-as.data.frame(pair.upper)

write.csv(pair.upper, "output-files/pair pheno corr table-CHRonly.csv")

#-------------------------------------------------------------------------------
# Check association of phenotypes with dist and dur
#-------------------------------------------------------------------------------

# visualize dist and clutch size
ggplot(dat, aes(x=as.factor(pat.clutch.size), y=process.dist.m)) +
  geom_boxplot() + geom_point(position=position_jitter(h=0, w=0.1),
                              shape=1) +
  xlab("No. offspring analyzed for paternity") +
  ylab("Distance between nests (m)")
ggsave("output-files/boxplot dist vs clutch size.png")

# visualize dur and clutch size
# these are correlated because fertile period depends on clutch size
ggplot(dat, aes(x=as.factor(pat.clutch.size), y=total.dur)) +
  geom_boxplot() + geom_point(position=position_jitter(h=0, w=0.1),
                              shape=1)+
  xlab("No. offspring analyzed for paternity") +
  ylab("No. days of overlap")
ggsave("output-files/boxplot dur vs clutch size.png")







