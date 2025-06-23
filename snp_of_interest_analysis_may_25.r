
##### Script to compare SNPs in Puerto Rico Aedes aegypti  #####
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


# subset to genes of interest
subset <- read.delim("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest.txt")
nrow(subset)
# subset to expanded list of genes of interest
#subset <- read.delim("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann_of_interest_expanded.txt")

## add column headers
colnames(subset) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELE", "DP", "AD", "ANN")
## remove all positions without a call
subset <- subset(subset, GT != "./.")

###################################################
################ CHECK HOW MANY ###################

## only missense SNPs
high <- read.csv("/mnt/storage12/emma/PR_snps/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest_gene_column_LOC.csv")


length(unique(high$gene_name))
unique(high$impact)

high_qual <- subset(high, DP > 10)

## add query to find how many samples have each gene
high_pos <- high %>% group_by(CHROM, POS, gene_name) %>% summarise(n =n())
high_pos <- high_qual %>% filter(GT != "0/0") %>% group_by(CHROM, POS, gene_name) %>% summarise(n_per_snp =n())
## filter so only common snps included
high_pos <- high_pos %>% filter(n_per_snp > 5)
## number of snps per gene - note does not reflect frequency
high_pos <- high_pos %>% group_by(CHROM, gene_name) %>% summarise(n_per_gene = n())
high_gene <- high_pos2 %>% group_by(gene_name) %>% summarise(n = n())

## Genes of interest regions
gff_genes <- read.delim("/mnt/storage12/emma/Reference_files/of_interest_no_contig.bed", header = FALSE)
interest <- read.delim("/mnt/storage12/emma/Reference_files/genes_of_interest.txt", header = FALSE)

gff_genes$gene_name <- str_extract(gff_genes$V4, "Name=[^;]+")
gff_genes$gene_name <- str_replace(gff_genes$gene_name, "Name=", "") 
length(unique(gff_genes$gene_name))

gff <- subset(gff_genes, gene_name %in% interest$V1)
length(gff$gene_name)

gff$length <- gff$V3 - gff$V2


### see how many snps per gene

high_pos$range <- gff$length[match(high_pos$gene_name, gff$gene_name)]
table(high_pos$range)

high_pos_rm <- high_pos[!is.na(high_pos$range),]

high_pos_rm$snp_ratio <- (high_pos_rm$n_per_gene/high_pos_rm$range)*100

### check how many genes of each type there are eg/ are there lots of cuticle proteins

descrip <- read.delim("/mnt/storage12/emma/Reference_files/of_interest_no_contig_with_descriptions.bed", header = FALSE)
descrip$gene_name <- str_extract(descrip$V4, "Name=[^;]+")
descrip$gene_name <- str_replace(descrip$gene_name, "Name=", "")

descrip$V5[descrip$V5 == "Description not found"] <- "Voltage-gated sodium channel"
high_pos_rm$descript <- descrip$V5[match(high_pos_rm$gene_name, descrip$gene_name)]

table(high_pos_rm$descript)


high_pos_rm$descript_simple <- NA
high_pos_rm$descript_simple[grep("cytochrome", high_pos_rm$descript)] <- "P450"
high_pos_rm$descript_simple[grep("cuticle", high_pos_rm$descript)] <- "Cuticle"
high_pos_rm$descript_simple[grep("G-", high_pos_rm$descript)] <- "G- protien coupled receptor"
high_pos_rm$descript_simple[grep("glutathione", high_pos_rm$descript)] <- "GST"
high_pos_rm$descript_simple[grep("esterase", high_pos_rm$descript)] <- "EST"
high_pos_rm$descript_simple[grep("acetylcholinesterase", high_pos_rm$descript)] <- "ACE"
high_pos_rm$descript_simple[grep("UDP", high_pos_rm$descript)] <- "UGT"
high_pos_rm$descript_simple[grep("glucuronosyltransferase", high_pos_rm$descript)] <- "UGT"
high_pos_rm$descript_simple[grep("GABA", high_pos_rm$descript)] <- "GABA"
high_pos_rm$descript_simple[grep("gamma-aminobutyric", high_pos_rm$descript)] <- "GABA"
high_pos_rm$descript_simple[grep("salivary", high_pos_rm$descript)] <- "Salivary"
high_pos_rm$descript_simple[grep("voltage-gated", high_pos_rm$descript)] <- "VGSC"

table(high_pos_rm$descript)
table(high_pos_rm$descript_simple)
## Majority are P450s, then EST, then cuticle

write.csv(high_pos_rm, "/mnt/storage12/emma/PR_snps/snps_matched_descript_of_interest.csv")



#####################################################


top_snp_genes <- high_pos_rm[high_pos_rm$snp_ratio > 1,]

ggplot(top_snp_genes, aes(x = gene_name, y = order(snp_ratio), fill = descript_simple)) +
  geom_boxplot()

library(forcats)

ggplot(top_snp_genes, aes(x = fct_reorder(gene_name, snp_ratio, .fun = median), 
                          y = snp_ratio, 
                          fill = descript_simple)) +
  geom_boxplot() +
  xlab("Gene name (ordered by SNP ratio)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(top_snp_genes, aes(x = fct_reorder(gene_name, snp_ratio, .desc = TRUE), 
                          y = snp_ratio, 
                          fill = descript_simple)) +
  geom_bar(stat = "identity") +
  xlab("Gene name (ordered by SNP ratio)") +
  ylab("SNP Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###################################################
################ Add metadata #####################
metadata <- read.csv("/mnt/storage12/emma/PR_combine/metadata_country_year.csv")
high$Country <- metadata$Country[match(high$SAMPLE, metadata$ID)]
high$year <- metadata$year[match(high$SAMPLE, metadata$ID)]
head(high)

################ By Country #####################
high_meta_qual <- subset(high, DP > 10)

## add query to find how many samples have each gene
high_pos <- high %>% group_by(CHROM, POS, gene_name) %>% summarise(n =n())
high_meta_pos <- high_meta_qual %>% filter(GT != "0/0") %>% group_by(CHROM, POS, Country, gene_name) %>% summarise(n_per_snp =n())
## filter so only common snps included
high_meta_pos <- high_meta_pos %>% filter(n_per_snp > 5)
## number of snps per gene - note does not reflect frequency
high_meta_pos <- high_meta_pos %>% group_by(CHROM, Country, gene_name) %>% summarise(n_per_gene = n())
high_gene <- high_pos2 %>% group_by(gene_name) %>% summarise(n = n())


### see how many snps per country and gene

high_meta_pos$range <- gff$length[match(high_meta_pos$gene_name, gff$gene_name)]
table(high_meta_pos$range)

high_meta_pos_rm <- high_meta_pos[!is.na(high_meta_pos$range),]

high_meta_pos_rm$snp_ratio <- (high_meta_pos_rm$n_per_gene/high_meta_pos_rm$range)*100

high_meta_pos_rm$descript <- high_pos_rm$descript[match(high_meta_pos_rm$gene_name, high_pos_rm$gene_name)]
high_meta_pos_rm$descript_simple <- high_pos_rm$descript_simple[match(high_meta_pos_rm$gene_name, high_pos_rm$gene_name)]


top_snp_meta_genes <- high_meta_pos_rm[high_meta_pos_rm$snp_ratio > 1,]

ggplot(top_snp_genes, aes(x = gene_name, y = order(snp_ratio), fill = descript_simple)) +
  geom_boxplot()

library(forcats)

grouped <- high_meta_pos_rm %>% group_by(Country, descript_simple) %>% summarise(mean = mean(snp_ratio))
ggplot(grouped, aes(x = fct_reorder(Country, mean), 
                          y = mean, 
                          fill = descript_simple)) +
                          geom_bar(stat = "identity") +
                          xlab("Gene name (ordered by SNP ratio)") +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(top_snp_genes, aes(x = fct_reorder(gene_name, snp_ratio, .fun = median), 
                          y = snp_ratio, 
                          fill = descript_simple)) +
                          geom_boxplot() +
                          xlab("Gene name (ordered by SNP ratio)") +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(top_snp_genes, aes(x = fct_reorder(gene_name, snp_ratio, .desc = TRUE), 
                          y = snp_ratio, 
                          fill = descript_simple)) +
                          geom_bar(stat = "identity") +
                          xlab("Gene name (ordered by SNP ratio)") +
                          ylab("SNP Ratio") +
                          theme_minimal() +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1))


























# snps per country
country_snp <- subset %>% group_by(Country) %>% summarise(n = n_distinct(CHROM,POS))
year_snp <- subset %>% group_by(year) %>% summarise(n = n_distinct(CHROM,POS))

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

unique(mis_high$Country[mis_high$GT != "0/0"])

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

high_qual <- subset(high DP > 10)


## add query to find how many samples have each gene
high_pos <- high %>% group_by(CHROM, POS, gene_name) %>% summarise(n =n())
high_pos <- high_qual %>% filter(GT != "0/0") %>% group_by(CHROM, POS, gene_name) %>% summarise(n =n())
high_gene <- high_pos %>% group_by(gene_name) %>% summarise(n = n())

query <- read.csv("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest.txt", sep = "\t", header = FALSE)
head(query)

colnames(query) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELES", "DP", "AD", "ANN")
query <- high
query <- query[query$POS %in% high_pos$POS,]
query <- query[query$GT != "./.",]
unique(query$GT)
query$GT[query$GT == "0|1"] <- "0/1"
query$GT[query$GT == "0|0"] <- "0/0" 
query$GT[query$GT == "0|2"] <- "0/2"
query$GT[query$GT == "2|2"] <- "2/2"
query$GT[query$GT == "1|1"] <- "1/1"

query_group_sum <- query%>% group_by(CHROM, POS, GT, REF, ALT) %>% summarise(n = n())

percent_genotype <- query_group_sum %>%
  group_by(CHROM, POS) %>%
  mutate(TOTAL = sum(n)) %>%
  group_by(CHROM, POS, GT) %>%
  summarise(GT_COUNT = sum(n), TOTAL = first(TOTAL), .groups = "drop") %>%
  mutate(PERCENT = round((GT_COUNT / TOTAL) * 100, 2)) %>%
  arrange(CHROM, POS, GT)


top_muts <- percent_genotype %>% filter(GT != "0/0") %>% top_n(10, PERCENT)

query_group <- merge(top_muts, high[, c("POS", "gene_name")], by = "POS", all.x = TRUE)
query_group <- query_group[duplicated(query_group$gene_name),]
library(dplyr)

# Step 1: Summarise to get count `n` for each unique combination
query_group_sum <- high %>%
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



