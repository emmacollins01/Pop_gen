
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
subset <- read.delim("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.missense_only.txt")
subset <- read.delim("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest.txt")
## add column headers
colnames(subset) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELE", "DP", "AD", "ANN")
## remove all positions without a call
subset <- subset(subset, GT != "./.")

## Add metadata
metadata <- read.csv("/mnt/storage12/emma/PR_combine/all_sample_country_metadata.csv")
subset$Country <- metadata$Country[match(subset$SAMPLE, metadata$ID)]
head(subset)

# snps per country
country_snp <- subset %>% group_by(Country) %>% summarise(n = n_distinct(CHROM,POS))

# snps per country
df_mut_sum <- subset %>% group_by(Country) %>% summarise(count = n_distinct(CHROM, POS))
df_mut_sum <- subset %>% group_by(Country) %>% summarise(count = n())
df_mut_sum2 <- df_mut_sum %>% group_by(Country) %>% summarise(n = n())
head(df_mut_sum)

subset$ANN[1]
str_split(subset$ANN[1], "\\|")
unlist(str_split(subset$ANN[1], "\\|"))[10]

str_split(subset$ANN[1], ",")

high <- subset[grep("HIGH", subset$ANN),]


subset$transcript_no <- NA

for (i in 1:nrow(subset)){

    subset$transcript_no[i] <- length(str_split(subset$ANN[i], ",")) 


}

## MISSENSE
missense <- df[grep("MODIFIER", df$ANN),]

write.csv(missense, "all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.missense.only.txt")