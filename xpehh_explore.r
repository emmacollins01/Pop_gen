rm(list = ls())
library(dplyr)

xp <- read.csv("/mnt/storage12/emma/selection/all_xpehh_significant.csv")
#xp <- read.csv("/mnt/storage12/emma/selection/new_GT_array/all_xpehh_significant.csv")
#xp <- read.csv("/mnt/storage12/emma/selection/MAF_test/all_xpehh_significant_MAF.csv")


## number of rows = 31969750

length(unique(xp$Position))
head(xp)

summary(xp$XPEHH)
# range -11.2 to 11.2

unique(xp$Chromosome)
# chrom names = 35107 etc

### check in IR genes

vgsc <- xp[(xp$Position > 315926360 & xp$Position < 316405639)  & xp$Chromosome == "35109",]
unique(vgsc$Chromosome)
rdl <- xp[(xp$Position > 41628484 & xp$Position < 41861946) & xp$Chromosome == "35108",]
ace <- xp[xp$Position > 161486025 & xp$Position < 161871375 & xp$Chromosome == "35109",]
gste <- xp[xp$Position > 351633367 & xp$Position < 351634957 & xp$Chromosome == "35109",]


pr <- subset(xp, Country1 == "Puerto Rico")
length(unique(pr$Position))
vgsc2 <- pr[(pr$Position > 315926360 & pr$Position < 316405639)  & pr$Chromosome == "35109",]
rdl2 <- pr[(pr$Position > 41628484 & pr$Position < 41861946) & pr$Chromosome == "35108",]
ace2 <- pr[pr$Position > 161486025 & pr$Position < 161871375 & pr$Chromosome == "35109",]
gste2 <- pr[pr$Position > 351633367 & pr$Position < 351634957 & pr$Chromosome == "35109",]
length(unique(paste(pr$Chromosome, pr$Position)))


xp_high <- xp %>% filter(XPEHH > 5 | XPEHH < -5)

summary(xp_high$XPEHH)
length(unique(xp_high$Position))

vgsc3 <- xp_high[(xp_high$Position > 315926360 & xp_high$Position < 316405639)  & xp_high$Chromosome == "35109",]
rdl3 <- xp_high[(xp_high$Position > 41628484 & xp_high$Position < 41861946) & xp_high$Chromosome == "35108",]
ace3 <- xp_high[xp_high$Position > 161486025 & xp_high$Position < 161871375 & xp_high$Chromosome == "35109",]
gste3 <- xp_high[xp_high$Position > 351633367 & xp_high$Position < 351634957 & xp_high$Chromosome == "35109",]


## Modal position- across group comparisons
xp_group <- xp %>% filter(XPEHH > 5 | XPEHH < -5) %>% 
                group_by(Chromosome, Position) %>% 
                summarise(n = n())


# look at highest value
test <- xp[xp$Position == "313184797",]
test$XPEHH <- abs(test$XPEHH)
class(test$XPEHH)
test2 <- test[!duplicated(round(test$XPEHH, 7)),]
xp$uni <- paste(xp$Country1, xp$Country2)

## check sign direction if needed
uganda_buk <- read.csv("/mnt/storage12/emma/selection/df_significant_xpehh_threshold_Uganda_vs_Burkina Faso.csv")
head(uganda_buk)
uganda_buk[uganda_buk$Position == "313184797",]
buk_uganda <- read.csv("/mnt/storage12/emma/selection/df_significant_xpehh_threshold_Burkina Faso_vs_Uganda.csv")
buk_uganda[buk_uganda$Position == "313184797",]

head(xp_group)
table(xp_group$n)


xp$Position <- as.numeric(xp$Position)

xp_high <- xp %>% filter(XPEHH > 5 | XPEHH < -5) 
xp_high_top <- xp_high %>% top_n(XPEHH, n = 20)
xp_high_bottom <- xp_high %>% top_n(XPEHH, n = -20)

xp$combo <- paste(xp$Country1, "_", xp$Country2)
library(ggplot2)
ggplot() +
geom_boxplot(data = xp, aes(x = combo,  y = XPEHH)) +
theme(axis.text = element_text(angle = 90))


PR <- xp_high %>% filter(xp_high$Country1 == "Puerto Rico" | xp_high$Country2 == "Puerto Rico")
nrow(PR)
PR <- PR[!duplicated(round(abs(PR$XPEHH), 10)),]
summary(PR$XPEHH)
table(PR$Country1)
table(PR$Country2)
table(PR$Country1, PR$Country2)
### Add gene to location

gff <- read.csv("/mnt/storage12/emma/Reference_files/genes_only.gff", sep = "\t", header = FALSE)
#gff <- read.csv("/mnt/storage12/emma/Reference_files/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3", sep = "\t", header = FALSE, skip = 2312)
#gff <- read.csv("/mnt/storage12/emma/Reference_files/Aedes_aegypti_lvpagwg.AaegL5.47.gff3", sep = "\t", header = FALSE, skip = 2314)

head(gff)

colnames(gff) <- c("Chromosome", "Ref", "Region", "Start", "End", "x1", "strand", "v2", "ANN")
gff_gene <- gff[gff$Region == "gene",]
gff_gene <- gff
gff_gene <- gff_gene[gff_gene$Chromosome %in% c("035107", "035108", "035109", "035159"),]
xp$Chromosome[xp$Chromosome == 35107] <- 035107
xp$Chromosome[xp$Chromosome == 35108] <- 035108
xp$Chromosome[xp$Chromosome == 35109] <- 035109
xp$Chromosome[xp$Chromosome == 35159] <- 035159

class(xp$Chromosome)
class(gff_gene$Chromosome)
gff_gene$Chromosome <-as.numeric(gff_gene$Chromosome)

unique(gff_gene$Chromosome)
chrom1 <- gff_gene[gff_gene$Chromosome == "35107",]
chrom2 <- gff_gene[gff_gene$Chromosome == "35108",]

xp_pos <- xp %>% select(Chromosome,Position)
xp_pos$uni <- paste0(xp_pos$Chromosome, xp_pos$Position)
length(unique(xp_pos$uni))
xp_pos2 <- xp_pos[!duplicated(xp_pos$uni),]


##########################################
# Merge the data frames based on chromosome and position range
gff_gene <- gff_gene[gff_gene$Region != "region",]

merged_df <- xp_pos2 %>%
  left_join(gff_gene, by = "Chromosome", relationship = "many-to-many") %>%
  filter(Position >= Start & Position <= End) %>%
  select(Chromosome, Position, ANN)

gene_na <- merged_df[is.na(merged_df$gene_name),]

# Create a function to extract and concatenate gene names
extract_genes <- function(ann) {
  genes <- str_extract_all(ann, "LOC[0-9]+")[[1]]  # Extract all matches
  if (length(genes) == 0) {
    return(NA)  # Return NA if no matches found
  } else {
    return(paste(unique(genes), collapse = "; "))  # Concatenate unique matches
  }
}

extract_product <- function(ann) {
  product <- str_extract_all(ann, "product=[^;]*;")[[1]]  # Extract all matches
  if (length(product) == 0) {
    return(NA)  # Return NA if no matches found
  } else {
    return(paste(unique(product), collapse = "; "))  # Concatenate unique matches
  }
}

# Apply the function to each row in the ANN column and store in gene_name column
merged_df <- merged_df %>%
  mutate(gene_name = sapply(ANN, extract_genes))

merged_df <- merged_df %>%
  mutate(product_name = sapply(ANN, extract_product))

head(merged_df)
length(unique(merged_df$gene_name))
# 1498
test_na <- merged_df[is.na(merged_df$gene_name),]

head(merged_df)
head(xp)
merged_df$uni_pos <- paste(merged_df$Chromosome, merged_df$Position)
xp$uni_pos <- paste(xp$Chromosome, xp$Position)

# gene name
xp$gene_name <- merged_df$gene_name[match(paste(xp$Chromosome, xp$Position), paste(merged_df$Chromosome, merged_df$Position))]
length(unique(xp$gene_name))
xp_na <- xp[is.na(xp$gene_name),]
head(xp_na)
# product name
merged_df_sub <- merged_df[complete.cases(merged_df$product_name),]
xp$product_name <- merged_df_sub$product_name[match(paste(xp$Chromosome, xp$Position), paste(merged_df_sub$Chromosome, merged_df_sub$Position))]
xp_na <- xp[is.na(xp$product_name),]
head(xp_na)

##########################################

xp_sig <- xp %>% filter(XPEHH > 5 | XPEHH < -5) 
xp_sig <- xp_sig[!duplicated(round(xp_sig$XPEHH), 7),]

length(unique(xp_sig$Position))

grep("LOC110678408", gff_gene$ANN)
grep("LOC5579970", gff_gene$ANN)

descrip <- read.csv("/mnt/storage12/emma/Reference_files/gene_id_product_uniq.txt", sep = "\t", header = FALSE)

xp_sig$descr <- descrip$V2[match(xp_sig$gene_name, descrip$V1)]
xp_sig_grp <- xp_sig %>% group_by(gene_name) %>% summarise(n = n())
xp_sig_grp_country <- xp_sig %>% group_by(Chromosome, gene_name, Country1, Country2, descr) %>% summarise(n = n())
xp_sig_grp_country2 <- xp_sig_grp_country %>% group_by(Chromosome, gene_name, descr) %>% summarise(no_comparisons = n())


xp_top <- xp_sig %>% top_n(abs(XPEHH), n = 20)
xp_bottom <- xp_sig %>% top_n(XPEHH, n = -10)

##########################################

### Look at genes of interest

interest <- read.table("/mnt/storage12/emma/Reference_files/genes_of_interest.txt")
interest <- read.table("/mnt/storage12/emma/Reference_files/enrichment/genes_name_only_expanded.txt")


xp_interest <- xp_sig[xp_sig$gene_name %in% interest$V1,]
length(unique(xp_interest$gene_name))
unique(xp_interest$gene_name)
xp_interest_na <- xp_interest[is.na(xp_interest$gene_name),]

summary(xp_interest$XPEHH)
xp_interest_top <- xp_interest %>% top_n(XPEHH, n = 10)

all_gene_xp <- unique(xp_sig[9])
write.csv(all_gene_xp, "selection/xpehh_significant_genes.csv")

xp_interest_grp <- xp_interest %>% group_by(Chromosome, Country1, Country2, gene_name, descr) %>% summarise(n =n())
length(unique(xp_interest_grp$gene_name))

xp_interest_grp2 <- xp_interest_grp %>% group_by(Chromosome, gene_name, descr) %>% summarise(n =n())

##########################################
### Other file with more genes

regions_gff <- read.table("/mnt/storage12/emma/Reference_files/region_of_interest.gff", sep = "\t", header = FALSE)

head(regions_gff)
regions_gff <- regions_gff[regions_gff$V1 %in% c("035107", "035108", "035109", "013159"),]
for ( i in 1:nrow(regions_gff)){
    regions_gff$gene[i] <- unique(unlist(str_extract_all(regions_gff$V9[i], "LOC[0-9]+")))
}

xp_regions <- xp_sig[xp_sig$gene_name %in% regions_gff$gene,]
xp_regions_na <- xp_regions[is.na(xp_regions$gene_name),]
length(unique(xp_regions$gene_name))

xp_regions_top <- xp_regions %>% top_n(XPEHH, n = 20)
xp_regions_top <- xp_regions %>% top_n(abs(XPEHH), n = 20)

table(xp_regions_top$Country1, xp_regions_top$Country2)

unique(xp_regions_top$descr)

xp_regions_grp_pos <- xp_regions %>% group_by(Chromosome, Position, gene_name, descr) %>% summarise(n= n())
xp_regions_grp <- xp_regions_grp_pos %>% group_by(gene_name, descr) %>% summarise(total=n())

##########################################
# Puerto Rico

PR_all <- xp %>% filter(Country1 == "Puerto Rico" | Country2 == "Puerto Rico")

PR <- xp_sig %>% filter(Country1 == "Puerto Rico" | Country2 == "Puerto Rico")

PR_regions <- xp_regions %>% filter(Country1 == "Puerto Rico" | Country2 == "Puerto Rico")

## look for genes found in ihs

# check genes 
grep("LOC5564970", xp_sig$gene_name)

ihs_comp <- xp_sig[xp_sig$gene_name == "LOC5566204",]

ihs_comp2 <- xp_sig[xp_sig$gene_name == "LOC5576390",]

# high affinity cAMP-specific and IBMX-insensitive 3',5'-cyclic phosphodiesterase 8 
sv_comp <- xp_sig[xp_sig$gene_name == "LOC5577718",]
# high affinity cGMP-specific 3',5'-cyclic phosphodiesterase 9A 
sv_comp2 <- xp_sig[xp_sig$gene_name == "LOC5572215",]
# no receptor potential A 
sv_comp3 <- xp_sig[xp_sig$gene_name == "LOC5571895",]
