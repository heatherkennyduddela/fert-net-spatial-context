
################################################################################
# For CHR in 2022, calculate "rest of barn" synchrony index for each dyad
# Heather Kenny-Duddela
# May 25, 2026
################################################################################

# Use Kempenaers 1993 breeding synchrony index calculation, inner term only
# average percent of females that are [simultaneously] fertile per day during 
# the breeding season. 0% means completely asynchronous, and 100% means 
# completely synchronous


# For each dyad, of the days that the focal female is fertile, count up the number of
# other females that were also fertile that day excluding the other female in
# that dyad. Sum over all of the fertile days for focal female p
# divide by: the number of fertile days for female p times (total num other females -1)

# NOTE: this script has been incorporated into the make-full-data-table script
# and does not need to be run separately

# load files and libraries------------------------------------------------------

## libraries
library(tidyverse)
library(lubridate)
library(ggplot2)

## files

# Matrix of fertile days for females at non-solitary sites in 2022
mat <- read.csv("02_output-files/season matrix multi sites.csv")

# table of pairwise clutches from make-full-data-table
nest.all2 <- read.csv("input-files/nest all table.csv")

#-------------------------------------------------------------------------------
# calculate "rest of barn SI" for each focal nesting attemtp at CHR
#-------------------------------------------------------------------------------

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



### Loop to calculate rest of barn SI for each dyad (by nesting attempt)

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

# split up labels columns
rest.SI.el.ids <- rest.SI.el %>%
  separate(nest.mom, c("location1", "mom1", "dad1", "attempt1"), "_",
           remove = FALSE) %>%
  separate(nest.dad, c("location2", "mom2", "dad2", "attempt2"), "_",
           remove = FALSE) %>%
  select(location1, location2, mom1, dad2, attempt1, attempt2,
         rest.SI)

rest.SI.el.ids$pairwise <- paste(rest.SI.el.ids$mom1, rest.SI.el.ids$dad2,
                                  sep="_")

# save table
write.csv(rest.SI.el.ids, file="03_output-files/rest of barn synchrony CHR.csv",
          row.names=F)


