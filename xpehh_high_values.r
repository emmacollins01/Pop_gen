library(dplyr)

setwd(('/mnt/storage12/emma/selection/'))
all <- read.csv("all_xpehh_significant.csv", sep = ",", header = TRUE)
head(all)
all$XPEHH <- as.numeric(all$XPEHH)

high <- all %>% filter(XPEHH > 5 | XPEHH < -5)
#high <- all[all$XPEHH >5 | all$XPEHH < -5,]
head(high)

high_pr <- high[high$Country1 == "Puerto Rico" | high$Country2 == "Puerto Rico",]
head(high_pr)
high_pr <- high_pr[!(duplicated(paste(high_pr$Chromosome,high_pr$Position))),]

top <- high_pr %>% top_n(XPEHH, n = 20)
