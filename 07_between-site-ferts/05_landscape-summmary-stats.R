################################################################################
# summary stats for landscape scale variables
# Heather Kenny-Duddela
# Aug 20, 2026
################################################################################

# load libraries
library(dplyr)

# load data
dat <- read.csv("03_output-files/all_dist_fert_between_site_final.csv")

#-------------------------------------------------------------------------------

# keep numeric columns only
dat2 <- dat[,5:8]

# make summary stat table (empty)
summary.stats.btwn <- as.data.frame(matrix(NA, nrow=3, ncol=4))
colnames(summary.stats.btwn) <- colnames(dat2)
# sample size
summary.stats.btwn[1,] <- colSums(!is.na(dat2))
# mean
summary.stats.btwn[2,] <- apply(dat2, 2 ,mean,na.rm=T)
# sd
summary.stats.btwn[3,] <- apply(dat2, 2, sd, na.rm=T)
# transpose
t.summary.stats.btwn <- as.data.frame(t(as.matrix(summary.stats.btwn)))
# add column labels
colnames(t.summary.stats.btwn) <- c("N","mean","sd")


# save
write.csv(t.summary.stats.btwn, "05_output-files/summary_stat_landscape.csv")



