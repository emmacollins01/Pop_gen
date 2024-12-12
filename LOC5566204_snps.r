
## Investigate GABA alpha gene LOC556204

##### Script to provide counts of number of SNPs in Illumina amplicons Puerto Rico Aedes aegypti data #####
rm(list = ls())
## load packages
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)

########################################################################
############################# ADD DATA #############################
#######################################################################

setwd("/mnt/storage12/emma/PR_snps/")

# all missense snps
subset <- read.delim("LOC5566204.txt", sep = "\t", header = FALSE)
nrow(subset)

## add column headers
colnames(subset) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELE", "DP", "AD", "ANN")
head(subset)

## remove all positions without a call
subset <- subset(subset, GT != "./.")
nrow(subset)

###################################################
################ CHECK HOW MANY ###################

length(unique(subset$POS))
# 64,579 removing no calls ./.
length(unique(subset$POS[subset$GT != "0/0"]))
# same if you remove reference calls

###################################################
################ Add metadata #####################
metadata <- read.csv("/mnt/storage12/emma/PR_combine/all_sample_country_metadata.csv")
subset$Country <- metadata$Country[match(subset$SAMPLE, metadata$ID)]
head(subset)

###################################################
################ Seperate gene, csq etc #####################


for (i in 1:nrow(subset)){
    #missense <- grep("missense_variant", unlist(str_split(subset$ANN[i], ",")))
    #missense_index <- missense[i]
    #to_extract <- unlist(str_split(subset$ANN[i], ","))[missense_index]
    subset$impact[i] <- unlist(str_split(subset$ANN[i], "\\|"))[2]
    subset$effect[i] <- unlist(str_split(subset$ANN[i], "\\|"))[3]
    subset$gene[i] <- unlist(str_split(subset$ANN[i], "\\|"))[5]
    subset$csq[i] <- unlist(str_split(subset$ANN[i], "\\|"))[11]
}

unique(subset$impact)
unique(subset$effect)
unique(subset$gene)
unique(subset$csq)

##
subset$DP <- as.numeric(subset$DP)
missense_df <- subset(subset, impact == "missense_variant")
missens_pos <- subset(missense_df, GT != "0/0")
length(unique(missens_pos$POS))
unique(missens_pos$POS)


list_pos <- c("28116470", "28116530", "28116545", "28116641", "28116656", "28116683", "28116926", "28116932", "28138163", "281158804")


sig_pos <- subset(subset, POS %in% list_pos)
sig_pos_dp <- subset(sig_pos, DP >= 10)
length(unique(sig_pos_dp$POS))
table(sig_pos_dp$POS, sig_pos_dp$Country)
unique(sig_pos_dp$impact)
