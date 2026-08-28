
################################################################################
# Function for calculating custom new data for raw or scaled predictors
# Heather Kenny-Duddela
# March 16, 2025
################################################################################

# libraries
library(tidyselect)

# old.data [original dataframe]
# data.type ["allFert", "epFert"]
# xvar [name of the x-variable to be used for plotting]
# xunscale [indicator for unscaled variable: "fem.TS", "fem.throat","fem.breast", 
#           "male.TS", "male.throat", "male.breast", "dist", "dur", "rest.SI.chr"]
# soc.pair [0, 1] for allFert only, 1 means predict for social pairs, 0 means non-social pairs

# Note that phenotypes were scaled using all of the birds present at CHR, 37 females
# and 32 males. This is the same set of birds included in the pheno object

# dist and dur and rest.SI were scaled using only the complete cases included in the file
# "dyad_ep_CHRonly_with fertile.csv" which is the data file used for fitting
# the final models. 



# For troubleshooting

# old.data <- dyad.ep2
# data.type <- "epFert"
# xvar <- "dist.scale"
# xunscale <- "dist"
# soc.pair <- NA


custom.newdata <- function(old.data, data.type, xvar, xunscale, soc.pair) {
  
  # load phenotype data, not in dyad format
  pheno <- read.csv("input-files/pheno-CHR-for-unscaling.csv")
  
  pheno.fem <- pheno %>% filter(sex=="F"|sex=="F?") %>%
    select(TS, throat_avg_bright, breast_avg_bright)
  
  stat.fem <- data.frame(stat = c("mean","sd"),
                         fem.TS = c(mean(pheno.fem$TS), sd(pheno.fem$TS)),
                         fem.throat = c(mean(pheno.fem$throat_avg_bright),
                                        sd(pheno.fem$throat_avg_bright)),
                         fem.breast = c(mean(pheno.fem$breast_avg_bright),
                                        sd(pheno.fem$breast_avg_bright)))
  
  pheno.male <- pheno %>% filter(sex=="M"|sex=="M?") %>%
    select(TS, throat_avg_bright, breast_avg_bright) 
  
  stat.male <- data.frame(stat = c("mean","sd"),
                         male.TS = c(mean(pheno.male$TS), sd(pheno.male$TS)),
                         male.throat = c(mean(pheno.male$throat_avg_bright, na.rm=T),
                                        sd(pheno.male$throat_avg_bright, na.rm=T)),
                         male.breast = c(mean(pheno.male$breast_avg_bright),
                                        sd(pheno.male$breast_avg_bright)))
  
  stat.dd <- data.frame(stat=c("mean","sd"),
                        dist = c(mean(old.data$process.dist.m),
                                 sd(old.data$process.dist.m)),
                        dur = c(mean(old.data$dur.fertile),
                                sd(old.data$dur.fertile)),
                        rest.SI = c(mean(old.data$rest.SI.chr),
                                sd(old.data$rest.SI.chr)))
  
  pheno.both <- cbind(stat.fem, stat.male[,2:4], stat.dd[,2:3])
  
  # Pull columns from allFert or epFert table
  if (data.type=="allFert") {
    list <- colnames(select(old.data, attempt2, soc.pair, dur.fertile,
                            dur.scale, process.dist.m, dist.scale, 
                            rest.SI.scale, rest.SI.chr,
                            shared_fert, pat.clutch.size,
                            bin_num_fert, prop.fert,
                            fem_throat_avg_bright:fem_breast_avg_bright_scaled, 
                            male_throat_avg_bright:male_breast_avg_bright_scaled))
  } else {
    list <- colnames(select(old.data, attempt2, soc.pair, dur.fertile, 
                            dur.scale, process.dist.m, dist.scale,
                            rest.SI.scale, rest.SI.chr,
                            shared_fert, bin_num_fert, prop.fert,
                            pat.clutch.size, 
                            fem_throat_avg_bright:fem_breast_avg_bright_scaled,
                            male_throat_avg_bright:male_breast_avg_bright_scaled,
                            soc.male_TS_scaled:soc.fem_breast_avg_bright_scaled))
  }
  
  # vector for range of xvar
  seq.xvar <- seq(range(old.data[ ,which(colnames(old.data)==xvar)], 
                        na.rm=T)[1], 
                      range(old.data[ ,which(colnames(old.data)==xvar)], 
                            na.rm=T)[2], 0.001)
  xvar.n <- length(seq.xvar)  
  
  ## Different variable names for allFert and epFert
  if (data.type=="allFert") {
    
    
    ## For scaled variables
      # generate new data
      newdata <- 
        data.frame(
          attempt2 = rep("first", xvar.n),
          soc.pair = rep(soc.pair, xvar.n),
          dur.scale = rep(0, xvar.n),
          dist.scale = rep(0, xvar.n),
          rest.SI.scale = rep(0,xvar.n),
          fem_TS_scaled = rep(mean(old.data$fem_TS_scaled), xvar.n),
          fem_throat_avg_bright_scaled = rep(mean(old.data$fem_throat_avg_bright_scaled), xvar.n),
          fem_breast_avg_bright_scaled = rep(mean(old.data$fem_breast_avg_bright_scaled), xvar.n),
          male_TS_scaled = rep(mean(old.data$male_TS_scaled), xvar.n),
          male_throat_avg_bright_scaled = rep(mean(old.data$male_throat_avg_bright_scaled), xvar.n),
          male_breast_avg_bright_scaled = rep(mean(old.data$male_breast_avg_bright_scaled), xvar.n))
      
      # add sequence for xvar
      newdata[, which(colnames(newdata)==xvar)] <- seq.xvar
    
      # un-scale the predictor
      newdata$unscale <- (newdata[, which(colnames(newdata)==xvar)] * 
                            pheno.both[2, 
                                       which(colnames(pheno.both)==xunscale)] +
                            pheno.both[1, 
                                       which(colnames(pheno.both)==xunscale)])
      # informative colname
      colnames(newdata)[12] <- xunscale
      
   
  
  # For epFert data table    
  } else {
    
      ## For scaled variables
      # generate new data
      newdata <- 
        data.frame(
          attempt2 = rep("first", xvar.n),
          dur.scale = rep(0, xvar.n),
          dist.scale = rep(0, xvar.n),
          rest.SI.scale = rep(0, xvar.n),
          fem_TS_scaled = rep(mean(old.data$fem_TS_scaled), xvar.n),
          fem_throat_avg_bright_scaled = rep(mean(old.data$fem_throat_avg_bright_scaled), xvar.n),
          fem_breast_avg_bright_scaled = rep(mean(old.data$fem_breast_avg_bright_scaled), xvar.n),
          male_TS_scaled = rep(mean(old.data$male_TS_scaled), xvar.n),
          male_throat_avg_bright_scaled = rep(mean(old.data$male_throat_avg_bright_scaled), xvar.n),
          male_breast_avg_bright_scaled = rep(mean(old.data$male_breast_avg_bright_scaled), xvar.n),
          soc.male_TS_scaled = rep(mean(old.data$soc.male_TS_scaled), xvar.n),
          soc.male_throat_avg_bright_scaled = rep(mean(old.data$soc.male_throat_avg_bright_scaled), xvar.n),
          soc.male_breast_avg_bright_scaled = rep(mean(old.data$soc.male_breast_avg_bright_scaled), xvar.n),
          soc.fem_TS_scaled = rep(mean(old.data$soc.fem_TS_scaled), xvar.n),
          soc.fem_throat_avg_bright_scaled = rep(mean(old.data$soc.fem_throat_avg_bright_scaled), xvar.n),
          soc.fem_breast_avg_bright_scaled = rep(mean(old.data$soc.fem_breast_avg_bright_scaled), xvar.n))
      
      # add sequence for xvar
      newdata[, which(colnames(newdata)==xvar)] <- seq.xvar
      
      # un-scale the predictor
      newdata$unscale <- (newdata[, which(colnames(newdata)==xvar)] * 
                            pheno.both[2, 
                                       which(colnames(pheno.both)==xunscale)] +
                            pheno.both[1, 
                                       which(colnames(pheno.both)==xunscale)])
      # informative colname
      colnames(newdata)[17] <- xunscale
  } 
  return(newdata)
} 
