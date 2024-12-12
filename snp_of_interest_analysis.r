
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
subset <- read.delim("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.missense_only.txt")
nrow(subset)
# subset to genes of interest
subset <- read.delim("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest.txt")
# subset to expanded list of genes of interest
subset <- read.delim("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann_of_interest_expanded.txt")

## add column headers
colnames(subset) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELE", "DP", "AD", "ANN")
## remove all positions without a call
subset <- subset(subset, GT != "./.")

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

# snps per country
country_snp <- subset %>% group_by(Country) %>% summarise(n = n_distinct(CHROM,POS))

test <- subset %>% filter(GT != "0/0") %>% group_by(Country) %>% summarise(n = n_distinct(CHROM,POS))
ggplot() +
geom_bar(data = test, aes(x = Country, y = n), stat = "identity")

# snps per country
df_mut_sum <- subset %>% filter(GT != "0/0") %>% group_by(Country) %>% summarise(count = n_distinct(CHROM, POS))
df_mut_sum2 <- subset %>% filter(GT != "0/0") %>% group_by(Country) %>% summarise(count = n())
head(df_mut_sum)

subset$ANN[1]
str_split(subset$ANN[1], "\\|")
unlist(str_split(subset$ANN[1], "\\|"))[10]

str_split(subset$ANN[1], ",")

# dont do this because it brings out any although some have high impact some have moderate
#high <- subset[grep("HIGH", subset$ANN),]

## EXTRACT MISSENSE ONLY ###

missense <- subset[grep("missense_variant", subset$ANN),]
missense_noref <- missense %>% filter(GT != "0/0")
length(unique(missense$POS))
length(unique(missense_noref$POS))

mis_high <- missense[grep("HIGH", missense$ANN),]
length(unique(mis_high$POS[mis_high$GT != "0/0"]))
length(unique(mis_high$Country[mis_high$GT != "0/0"]))

unique(mis_high$Country)

missense_country <- missense %>% filter(GT != "0/0") %>% group_by(Country) %>% summarise(n = n_distinct(CHROM, POS))

missense_country$CHROM[missense_country$CHROM == "NC_035107.1"] <- "1"
missense_country$CHROM[missense_country$CHROM == "NC_035108.1"] <- "2"
missense_country$CHROM[missense_country$CHROM == "NC_035109.1"] <- "3"

ggplot() + 
geom_bar(data = missense_country, aes(x = Country, y = n, fill = CHROM), stat = "identity") +
#facet_grid(~CHROM) +
scale_fill_manual(values = c("#40476D", "#749cb4",  "#A6C48A", "#51A3A3","#EE964B"), name = "Chromosome") +
ylab("Number of non-synonymous SNPs") +
xlab("") +
theme_minimal() +
theme(axis.text.x = element_text(size = 28, angle = 65, hjust =0.7),
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 30),
        legend.position = "bottom",
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28))

####### extract from annotation ####### 
high <- missense

for (i in 1:nrow(high)){
    missense <- grep("missense_variant", unlist(str_split(high$ANN[i], ",")))
    missense_index <- missense[1]
    to_extract <- unlist(str_split(high$ANN[i], ","))[missense_index]
    high$impact[i] <- unlist(str_split(to_extract, "\\|"))[2]
    high$effect[i] <- unlist(str_split(to_extract, "\\|"))[3]
    high$gene[i] <- unlist(str_split(to_extract, "\\|"))[5]
    high$csq[i] <- unlist(str_split(to_extract, "\\|"))[11]
}

write.csv(high, "/mnt/storage12/emma/PR_snps/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest_gene_column.csv")

for (i in 1:nrow(high)){
    high$gene_name <- str_extract(high$ANN[i], "LOC[0-9]+")
}

#write.csv(high, "/mnt/storage12/emma/PR_snps/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest_gene_column_LOC.csv")

high <- read.csv("/mnt/storage12/emma/PR_snps/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest_gene_column_LOC.csv")

length(unique(high$gene_name))
unique(high$gene_name)
unique(high$impact)


## add query to find how many samples have each gene
high_pos <- high %>% group_by(CHROM, POS, gene_name) %>% summarise(n =n())
high_gene <- high_pos %>% group_by(gene_name) %>% summarise(n = n())

query <- read.csv("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest.txt", sep = "\t", header = FALSE)
head(query)
colnames(query) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELES", "DP", "AD", "ANN")
query <- query[query$POS %in%high$POS,]
query <- query[query$GT != "./.",]
unique(query$GT)
query$GT[query$GT == "0|1"] <- "0/1"
query$GT[query$GT == "0|0"] <- "0/0" 
query$GT[query$GT == "0|2"] <- "0/2"
query$GT[query$GT == "2|2"] <- "2/2"
query$GT[query$GT == "1|1"] <- "1/1"

query_group_sum <- query%>% group_by(CHROM, POS, GT, REF, ALT) %>% summarise(n = n())


library(dplyr)

# Step 1: Summarise to get count `n` for each unique combination
query_group_sum <- query %>%
  group_by(CHROM, POS, GT, REF, ALT, gene_name) %>%
  summarise(n = n(), .groups = 'drop')

# Step 2: Calculate the total count per position
query_group_sum <- query_group_sum %>%
  group_by(CHROM, POS) %>%
  mutate(total_per_position = sum(n)) %>%
  ungroup()

# Step 3: Calculate the total count per genotype
query_group_sum <- query_group_sum %>%
  group_by(GT) %>%
  mutate(total_per_GT = sum(n)) %>%
  ungroup()

# Print the result
print(query_group_sum)


# Merge the data frames based on the POS column
query_group <- merge(query_group_sum, high[, c("POS", "gene_name")], by = "POS", all.x = TRUE)

# number of genes
length(unique(query_group$gene_name[query_group$GT != "0/0"]))
head(query_group)
#qur <- query_group %>% group_by(CHROM, POS, GT, REF, ALT, gene_name) %>% summarise(n = n()) %>% mutate(total = sum(n))
qur <- query_group_sum %>% group_by(CHROM, POS, GT, REF, ALT, n) %>% mutate(total = sum(n), .groups ='drop')
# Sample pipeline to find total by GT and total per position
qur <- query_group %>%
  group_by(CHROM, POS, GT, REF, ALT, gene_name) %>%
  summarise(n = n()) %>% # Get the count for each combination
  group_by(CHROM, POS, REF, ALT) %>%
  mutate(total_per_position = sum(n)) 
  
  %>% # Calculate the total count per position
  ungroup() %>%
  group_by(GT) %>%
  mutate(total_per_GT = sum(n)) %>% # Calculate the total count per GT
  ungroup()

# Print the result
print(qur)


# Print the updated query_group
print(query_group)


unique(query_group$gene_name)

gene <- as.data.frame(high$gene)
gene[1]
grep("GENE_exon-XM_021839123.1-1", gff$V9)
length(unique(gene$'high$gene'))
colnames(gene) <- c("gene_name")
gene_sum <- gene %>% group_by(gene_name) %>% summarise(n = n())


####### match exon to gene ####### 
gff <- read.csv("/mnt/storage12/emma/Reference_files/genes.gff", header = FALSE, comment.char = '#', sep = "\t")
head(gff)

for (i in 1:nrow(gene)){


}

for (i in 1:nrow(high)){
    exon_name <- str_sub(high$gene[i], 11, nchar(high$gene[i]))
    line <- gff$V9[grep(exon_name, gff$V9)][1]
    gene_id <- unlist(str_split(line, ";"))[5]
    high$gene_id[i] <- str_sub(gene_id, 6, nchar(gene_id))
}







## match to fst values

fst <- read.csv("/mnt/storage12/emma/PR_fst/fst_output/all_countries_fst_2024-03-19.csv")
head(fst)

fst_high <- subset(fst, WEIR_AND_COCKERHAM_FST > 0.9)

for (i in 1:nrow(fst_high)){
    fst_high$ANN[i] <- high$ANN[match(paste0(fst_high$CHROM, fst_high$POS), paste0(high$CHROM, high$POS))]
}



subset$transcript_no <- NA

for (i in 1:nrow(subset)){

    subset$transcript_no[i] <- length(str_split(subset$ANN[i], ",")) 


}

## MISSENSE
missense <- df[grep("MODIFIER", df$ANN),]

write.csv(missense, "all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.missense.only.txt")



## check GSTE2 Leu111 mutation

test <- subset(subset, POS %in% c(351634048, 351634049))


testt <- subset(subset, POS %in% c(351634752, 351634753))



