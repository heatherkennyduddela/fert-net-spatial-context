
################################################################################
# script for combining full sib/half sib plots from each year
# Heather Kenny-Duddela
# June 12, 2026
################################################################################

# load data
load("input-files/sib_categories_2021_plot.Rdata")
load("input-files/sib_categories_2022_plot.Rdata")
load("input-files/sib_categories_2023_plot.Rdata")

# libraries
library(ggplot2)
library(ggpubr)

#-------------------------------------------------------------------------------
# make combined plot
#-------------------------------------------------------------------------------

new21 <- sib_categories_2021 +
  scale_color_manual(values=c("darkred","orange", "#666666")) +
  theme_light() +
  ggtitle("2021") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Kinship coefficient (pi_HAT)") + xlab(NULL)

new22 <- sib_categories_2022 +
  scale_color_manual(values=c("darkred","orange", "#666666")) +
  theme_light() +
  ggtitle("2022") +
  ylab(NULL) + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Probablity of unrelated (k0_hat)")

new23 <- sib_categories_2023 +
  scale_color_manual(values=c("darkred","orange", "#666666")) +
  theme_light() +
  ggtitle("2023") +
  ylab(NULL) + theme(plot.title = element_text(hjust = 0.5)) +
  xlab(NULL)


ggarrange(new21, new22, new23, 
          common.legend = T, labels=c("A","B","C"),
          ncol=3, nrow=1, align = "hv") 

ggsave("04_output-files/combined sib categories plots.png", h=3.5, w=9)
