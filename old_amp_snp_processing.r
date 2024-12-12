
##### Script to provide counts of number of SNPs in Illumina amplicons Puerto Rico Aedes aegypti data #####
rm(list = ls())
## load packages
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(tidyr)

########################################################################
############################# ADD DATA #############################
#######################################################################

### Puerto Rico data 
df <- read.delim("/mnt/storage12/emma/PR_snps/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.IR_genes.txt", sep = "\t", header = FALSE)

## add column headers
colnames(df) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELE", "DP", "AD", "ANN")
## remove all positions without a call
df <- subset(df, GT != "./.")

df$GT[df$GT == "0|1"] <- "0/1"
df$GT[df$GT == "1|1"] <- "1/1"

## Add metadata
metadata <- read.csv("/mnt/storage12/emma/PR_combine/all_sample_country_metadata.csv")
df$Country <- metadata$Country[match(df$SAMPLE, metadata$ID)]
df$Region <- metadata$Region[match(df$SAMPLE, metadata$ID)]

head(df)
## Add snp locations
SNP_locs <- read.delim("/mnt/storage12/emma/PR_snps/IR_mut_pos.txt", header = TRUE)

# snps per country
country_snp <- df %>% group_by(Country) %>% summarise(n = n_distinct(CHROM,POS))

unique(df$POS)

## match to gene
df$gene <- NA
for (i in 1:nrow(df)){

  if(grepl("LOC5567355", df$ANN[i])){
    df$gene[i] <- "VGSC"
  }
  else if(grepl("LOC5570466", df$ANN[i])){
    df$gene[i] <- "RDL"
  }
  else if(grepl("LOC5578456", df$ANN[i])){
    df$gene[i] <- "ACE"
  }
  else if(grepl("LOC110676855", df$ANN[i])){
    df$gene[i] <- "GSTE"
  }
  else{
    df$gene[i] <- "other"
  }
}

unique(df$gene)
table(df$gene, df$Country)

missense <- df[grep("missense_variant", df$ANN),]
missense_no_ref <- missense[missense$GT != "0/0",]
length(unique(missense$POS))
length(unique(missense_no_ref$POS))
gene_mis <- missense %>% group_by(gene) %>% summarise(n = n_distinct(POS))
gene_mis
gene_mis <- missense %>% group_by(gene, Country) %>% summarise(n = n_distinct(POS))
gene_mis
unique(missense$POS)

ggplot(data = gene_mis, aes(x = gene, y =n,  fill = Country))+
  geom_bar(stat = "identity", position = position_dodge(width =1))


## allele frequency per country and region
missense_region <- missense %>% group_by(CHROM, POS, GT, Region) %>% summarise(count = n())
missense_region_sum <- missense %>% group_by(CHROM, POS, Region) %>% summarise(total = n())
missense_region$total <- missense_region_sum$total[match(paste0(missense_region$POS, missense_region$Region), paste0(missense_region_sum$POS, missense_region_sum$Region))]

missense_region_wide <- spread(missense_region, GT, count)
missense_region_wide$'0/0'[is.na(missense_region_wide$'0/0')] <- 0
missense_region_wide$'0/1'[is.na(missense_region_wide$'0/1')] <- 0
missense_region_wide$'1/1'[is.na(missense_region_wide$'1/1')] <- 0

missense_region_wide$af <- (missense_region_wide$'0/1'+(missense_region_wide$'1/1'*2)) / (missense_region_wide$total*2)

missense_wider <- missense_region_wide %>% select(-total, -'0/0', -'0/1', -'1/1')
missense_wider <- spread(missense_wider, Region, af)

write.csv(missense_wider, "../af_per_region_table.csv")


## allele frequency for Puerto Rico
missense_country <- missense %>% group_by(CHROM, POS, GT, Country) %>% summarise(count = n())
missense_country_sum <- missense %>% group_by(CHROM, POS, Country) %>% summarise(total = n())
missense_country$total <- missense_country_sum$total[match(paste0(missense_country$POS, missense_country$Country), paste0(missense_country_sum$POS, missense_country_sum$Country))]

missense_country_wide <- spread(missense_country, GT, count)
missense_country_wide$'0/0'[is.na(missense_country_wide$'0/0')] <- 0
missense_country_wide$'0/1'[is.na(missense_country_wide$'0/1')] <- 0
missense_country_wide$'1/1'[is.na(missense_country_wide$'1/1')] <- 0

missense_country_wide$af <- (missense_country_wide$'0/1'+(missense_country_wide$'1/1'*2)) / (missense_country_wide$total*2)

missense_wider <- missense_country_wide %>% select(-total, -'0/0', -'0/1', -'1/1')
missense_wider <- spread(missense_wider, Country, af)

write.csv(missense_wider, "../af_per_country_table.csv")

############ synonymous #############

synon <- df[grep("synonymous", df$ANN),]
length(unique(synon$POS))

table(synon$gene)
synon_gene <- synon %>% group_by(gene) %>% summarise(n = n_distinct(POS))
synon_gene <- synon %>% group_by(gene, Country) %>% summarise(n = n_distinct(POS))



###################################################
################### ONLY IR SNPS #################
###################################################
# IR snps per country
df_mut <- df[df$POS %in% SNP_locs$POS,]
df_mut <- df
for (i in 1:nrow(df_mut)){
  df_mut$effect[i] <- unlist(strsplit(df_mut$ANN[i], "\\|"))[[2]]
  df_mut$csq[i] <- df_mut
}
df_mut <- df_mut[df_mut$effect != "synonymous_variant",]
unique(df_mut$POS)
df_mut_sum <- df_mut %>% group_by(SAMPLE, Country) %>% summarise(count = n())
df_mut_sum <- df_mut_sum %>% group_by(Country) %>% summarise(n = n())
head(df_mut_sum)

unique(df_mut$effect)
table(df_mut$gene[df_mut$effect == "missense_variant" | df_mut$effect == "missense_variant&splice_region_variant"])

gene_mis <- df_mut[df_mut$effect == "missense_variant"| df_mut$effect == "missense_variant&splice_region_variant",]
gene_mis2 <- gene_mis %>% group_by(as.factor(gene)) %>% summarise(n = n_distinct(POS))

sum <- gene_mis %>% group_by(CHROM, POS, REF, ALT, GT, ALLELE, )

# number of IR mutations per sample 0-4
df_mut_sample <- df_mut %>% filter(GT != "0/0") %>% dplyr::group_by(SAMPLE, Country) %>% summarise(count = n_distinct(POS))
# number of IR mutations per country 0-4
df_mut_country <- df_mut %>% filter(GT != "0/0") %>% group_by(Country) %>% summarise(count = n_distinct(POS))

## MISSENSE
missense <- df[grep("MODIFIER", df$ANN),]

########################################################################
############################# PLOT IR SNPS  #############################
#######################################################################

### plots for Country and insecticide summary
ggplot(data = df_mut_sum, aes(x = reorder(Country, n), y = n)) +
  geom_bar(fill = "#465c7a", stat = "identity") +
  xlab("Country") +
  ylab("Frequency") +
  geom_text(aes(x = reorder(Country, n), y = n, label = n), size = 12, vjust = -0.2) +
  theme_classic() +
  theme(axis.text = element_text(size = 30, angle = 90),
        axis.title = element_text(size = 30))

df_mut_sum2 <- df_mut %>% 
                group_by(Country, POS, GT) %>% 
                summarise(count = n(), .groups = "drop")

df_mut_sum2 <- df_mut_sum2[df_mut_sum2$GT !=  "0/0" & df_mut_sum2$GT != "./.",]

df_mut_sum3 <- df_mut_sum2 %>% group_by(Country, POS) %>% summarise(n = sum(count))
df_mut_sum4 <- df_mut_sum2 %>% group_by(Country, POS) %>% summarise(prop = sum(proportion))
df_total <- df_mut_sum3 %>% group_by(Country) %>% summarise(total = sum(n))
df_mut_sum2$total <- df_total$total[match(df_mut_sum2$Country, df_total$Country)]
df_mut_sum2$proportion <- signif((df_mut_sum2$count/df_mut_sum2$total)*100, 3)

# test they are all 100%
test <- df_mut_sum2 %>% group_by(Country) %>% summarise(sum = sum(proportion))

mycols <- c("#ffaa5c", "#db7376", "#52a884", "#bc80bd", "#465c7a")
mycols <- c("#40476D", "#EE964B", "#749cb4",  "#A6C48A",  "#51A3A3")
mycols <- c("#40476D", "#749cb4",  "#A6C48A", "#bc80bd")
#mycols <- unikn::usecol(c("#CDD3D5", "#75B8C8", "#197278", "#03045E", "#922D50"), n = 5)
#mycols <- c("#40476D", "#678D58", "#A6C48A","#52a884","#bc80bd","#51A3A3", "#7D387D", "#bc80bd")

total_only <- df_mut_sum2 %>% group_by(Country) %>% summarise(total = sum(count))


plot1 <- ggplot(data = df_mut_sum2, aes(x = reorder(Country, count), y = count)) +
  geom_bar(stat = "identity", aes(fill = as.factor(POS))) +
  xlab("") +
  ylab("Frequency of insecticide resistance mutations") +
    scale_fill_manual(values = mycols) +
  #facet_wrap(~POS, ncol = 5) +
  #geom_text(aes(x = Country, label = total), size = 8, position = position_stack(height = 1)) +
  #geom_text(aes(label = after_stat(y), group = POS), stat = 'summary', fun = sum, vjust = -1, size = 14) +
  stat_summary(fun = sum, geom = "text", aes(label = after_stat(y)), vjust = -1, size = 8) +
  theme_classic() +
  labs(fill = "Mutation Position") +
  theme(axis.text = element_text(size = 24, angle = 90),
        axis.title = element_text(size = 30),
        strip.text = element_text(size = 24),
        #legend.position = "none",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

plot1

plot2 <- ggplot(data = df_mut_sum2, aes(x = reorder(Country, count), y = proportion)) +
  geom_bar(stat = "identity", aes(fill = as.factor(POS))) +
  xlab("") +
  ylab("Frequency") +
  labs(fill = "Mutation Position") +
  scale_y_continuous(limits = c(0,100.1)) +
  scale_fill_manual(values = mycols) +
  #geom_text(aes(x = reorder(Country, n), y = signif(proportion,3), label = proportion), 
   #         size = 8, position = position_stack(0.4)) +
 #geom_text(data = total_only, aes(x = Country, label = total), 
  #          size = 8, position = position_stack(vjust = 1.1)) +
  theme_classic() +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text = element_text(size = 24, angle = 90),
        axis.title = element_text(size = 30),
        strip.text = element_text(size = 24),
        legend.position = 'right',
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

plot2

library(cowplot)
library(patchwork)
plot_grid(plot1, plot2)
combined_plot <- plot1 + plot2 + plot_layout(ncol= 2, nrow = 1, widths = c(10.5,11))
combined_plot

dim(df)
length(unique(df$POS))

df_mut_prop  <- df %>% 
                filter(POS %in% c(315939224, 315983763, 316080722,  41847790)) %>%
                group_by(Country, POS, GT) %>% 
                mutate(total = n()) %>% 
                summarise(count = n(), .groups = "drop")


df_total <- df_mut_prop %>% group_by(Country) %>% summarise(total = sum(count))
df_mut_prop$total <- df_total$total[match(df_mut_prop$Country, df_total$Country)]
df_mut_prop$proportion <- signif((df_mut_prop$count/df_mut_prop$total)*100, 3)


## make pie chart to show proportions
mycols2 <- c("#40476D", "#678D58","#bc80bd","#51A3A3", "#7D387D", "#bc80bd")
mycols2 <- mycols
plot3 <- ggplot(df_mut_prop, aes(x="", y=proportion, fill = as.factor(POS))) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = mycols) +
  coord_polar("y", start=0) +
  theme_void() +
  facet_wrap(~Country, nrow = 4) +
  #ggtitle("Proportion of each insecticide resistance mutation per country") +
  guides(fill = guide_legend(title = "No. of insecticide \nresistance mutations"), size = 18, ncol = 1) +
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 18)
  ) #,panel.background = element_rect(fill = 'grey55')) +
  #geom_text(aes(x = 1.6,y = sum, label =sum), position = position_stack(vjust = .5), color = "black", size = 4)


combined_plot <- plot1 + plot3 + plot_layout(ncol= 2, nrow = 1, widths = c(10.5,11), heights = c(6,15))
combined_plot



### calculate allele frequency 
library(tidyr)
df_mut$GT[df_mut$GT == "0|1"] <- "0/1"
df_mut$GT[df_mut$GT == "1|1"] <- "1/1"
af <- df_mut %>% 
                group_by(Country, POS, GT) %>% 
                summarise(count = n(), .groups = "drop")

af <- af[af$GT != "./.",]
af2 <- spread(af, GT, count)

af2$'0/0'[is.na(af2$'0/0')] <- 0
af2$'0/1'[is.na(af2$'0/1')] <- 0
af2$'1/1'[is.na(af2$'1/1')] <- 0


af2$total <- (af2$'0/0' + af2$'0/1' + af2$'1/1')
af2$af <- ((af2$'0/1' + (af2$'1/1'*2)) / (af2$total*2))*100
mycols <- c("#40476D", "#678D58", "#A6C48A","#52a884","#bc80bd","#51A3A3", "#7D387D", "#c7522a")

ggplot(data = af2, aes(x = as.factor(Country), y = af)) +
  geom_bar(stat = "identity", aes(fill = as.factor(POS)), position = "dodge") +
  xlab("") +
  ylab("Allele Frequency") +
    scale_fill_manual(values = mycols) +
  facet_wrap(~Country, scales = "free_x", ncol = 8) +
  #geom_text(aes(x = Country, label = total), size = 8, position = position_stack(height = 1)) +
  #geom_text(aes(label = after_stat(y), group = POS), stat = 'summary', fun = sum, vjust = -1, size = 14) +
  #stat_summary(fun = sum, geom = "text", aes(label = after_stat(y)), vjust = -1, size = 8) +
  theme_bw() +
  labs(fill = "Mutation Position") +
  theme(axis.text = element_text(size = 24, angle = 90),
        axis.title = element_text(size = 30),
        strip.text = element_text(size = 24),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))


ggplot(af2, aes(x="", y=af, fill = as.factor(POS))) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = mycols) +
  coord_polar("y", start=0) +
  theme_void() +
  facet_wrap(~Country, nrow = 4) +
  #ggtitle("Proportion of each insecticide resistance mutation per country") +
  guides(fill = guide_legend(title = "No. of insecticide \nresistance mutations"), size = 18, ncol = 1) +
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 18)
  ) #,panel.background = element_rect(fill = 'grey55')) +
  #geom_text(aes(x = 1.6,y = sum, label =sum), position = position_stack(vjust = .5), color = "black", size = 4)


## group by number of mutations per sample

summary <- df_mut %>% 
            filter(GT != "0/0") %>%
            group_by(SAMPLE, Country) %>% 
            summarise(count_mut = length(unique(POS)), .groups = "drop")


summary2 <- summary %>% group_by(Country, count_mut) %>% 
        summarise(sum_mut = n()) %>% 
        mutate(total = sum(sum_mut), prop = sum_mut/total)


ggplot(summary2, aes(x="", y=prop, fill = as.factor(count_mut))) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = mycols) +
  coord_polar("y", start=0) +
  theme_void() +
  facet_wrap(~Country, nrow = 4) +
  ggtitle("Proportion of samples with each number \nof insecticide resistance mutations") +
  guides(fill = guide_legend(title = "No. of insecticide \nresistance mutations"), ncol = 1) +
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"),
    title = element_text(size =22),
    #legend.title = element_text(size = 22),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 18))






####################################################################################################







# make depth numeric and only include > 10
df$DP <- as.numeric(df$DP)
df <- subset(df, DP >10)

## check how many samples
length(unique(df$SAMPLE))
length(unique(df$Country))
length(unique(df$POS))

## only keep snps where AR is at least 0.05 or AD > 2
df_AD <- df[grep("AD", colnames(df))]
df_AD_split <- (str_split_fixed(as.character(df_AD$AD), ",", max(lengths(gregexpr("([^,]+)", df_AD$AD)))))
df_AD_split <- as.data.frame(df_AD_split)
no_ADs <- seq(1, max(lengths(gregexpr("([^,]+)", df_AD$AD))))
colnames(df_AD_split) <- rep(paste0("AD", no_ADs))
df_bind <- cbind(df, df_AD_split)
length(unique(paste(df$POS, df$Sample_file)))

rm( df_AD, df_AD_split, df)

df_bind$DP <- as.numeric(as.character(df_bind$DP))
df_bind$AD1 <- as.numeric(as.character(df_bind$AD1))
df_bind$AD2 <- as.numeric(as.character(df_bind$AD2))

df_bind$AD1[is.na(df_bind$AD1)] <- 0
df_bind$AD2[is.na(df_bind$AD2)] <- 0


for (i in 1:nrow(df_bind)){
  df_bind$DP_sum[i] <- sum(df_bind[i,grep("AD[0-9]", colnames(df_bind))])
}

df_bind$AR_1 <- df_bind$AD1/df_bind$DP_sum
df_bind$AR_2 <- df_bind$AD2/df_bind$DP_sum

df_bind$GT <- gsub("\\|", "/", df_bind$GT)

list_rows <- list()
for (i in 1:nrow(df_bind)){
  
  if(unlist(str_split(df_bind$GT[i], "/"))[1] != unlist(str_split(df_bind$GT[i], "/"))[2]){
    
    minimum_ad <- min(ads[1:2])
    
    ars <- df_bind[i,grep("AR_[0-9]", colnames(df_bind))][order(-unlist(df_bind[i,grep("AR_[0-9]", colnames(df_bind))]))]
    minimum_ar <- min(ars[1:2])
    
    
    if(minimum_ad <2 | minimum_ar < 0.05 | df_bind$DP_sum[i] < 20){
      print(df_bind$POS[i])
      print(df_bind$SAMPLE[i])
      list_rows <- c(list_rows, i) 
    }
  }
}

length(list_rows)

df_bind_filt <- df_bind[-unlist(list_rows),]
length(unique(df_bind$POS))

df <- df_bind_filt



# check how many have calls at each position
df_grouped <- df %>% group_by(CHROM, POS) %>% summarise(n = n())
# remove where only have < 50 samples at each position
df_grouped_high <- subset(df_grouped, n > 50)
df <- subset(df, POS %in% df_grouped_high$POS)

## remove SNPS where only homozygous ref called
pos_list <- list()

for (i in unique(df$POS)){
  
  df_i <- subset(df, POS == i)
  
  if(length(unique(df_i$GT)) == 1 & unique(df_i$GT)[1] == "0/0"){
    
    pos_list <- c(pos_list, i)
    
  }
}

length(pos_list)

df <- df[!(df$POS %in% pos_list),]

length(unique(df$POS))
# 49 SNPs left

## Add the metadata to the SNPs dataframe
df$Country <- metadata$Country[match(df$SAMPLE, metadata$ID)]
df$Region <- metadata$Region[match(df$SAMPLE, metadata$ID)]
df$Subregion <- metadata$Subregion[match(df$SAMPLE, metadata$ID)]

# check which Countrys are included
unique(df$Country)



samples_kept <- as.data.frame(unique(df$SAMPLE))
length(unique(df$SAMPLE))
length(unique(df$POS))
write.csv(samples_kept, paste0("/mnt/storage12/emma/PR_snps/samples_kept", Sys.Date(), ".csv"))

snp_pos <- as.data.frame(unique(df$POS))
write.csv(snp_pos, paste0("/mnt/storage12/emma/PR_snps/snp_pos_kept", Sys.Date(), ".csv"))

write.csv(df, paste0("/mnt/storage12/emma/PR_snps/filtered_snps", Sys.Date(), ".csv"))


###################################################################################################
###################################### stopping point 1 #########################################
###################################################################################################

## Summary of number of SNPS
length(unique(df$SAMPLE))
# now 212 samples
length(unique(df$POS))
# 49

df_sub <- df
df_sub <- read.csv("/mnt/storage12/emma/PR_snps/filtered_snps2024-03-01.csv")
## Summary by Country of SNPS
snps <- df_sub %>% group_by(SAMPLE, Country) %>% summarise(count = n())
snp_loc <- snps %>% group_by(Country) %>% summarise(n = n())

SNP_locs <- read.delim("/mnt/storage12/emma/PR_snps/IR_mut_pos.txt", header = TRUE)
mut_pos <- SNP_locs$POS
snps_mut <- df_sub[df_sub$POS %in% mut_pos,]
snps <- df_sub %>% group_by(SAMPLE, Country) %>% summarise(count = n())
snp_loc <- snps %>% group_by(Country) %>% summarise(n = n())

### plots for Country and insecticide summary
mycols <-  unikn::usecol(c("#ffaa5c", "#db7376", "#52a884", "#bc80bd", "#465c7a"), n = 8)
ggplot(data = snp_loc, aes(x = reorder(Country, n), y = n)) +
  geom_bar(fill = mycols, stat = "identity") +
  xlab("Country") +
  ylab("Frequency") +
  geom_text(aes(x = reorder(Country, n), y = n, label = n), size = 12, vjust = -0.2) +
  theme_classic() +
  theme(axis.text = element_text(size = 30, angle = 90),
        axis.title = element_text(size = 30))

# check adds up
sum(snp_loc$n)

## Summary of exposure to insecticides
insect <- df %>% group_by(SAMPLE, insecticide) %>% summarise(count = n())
insect_sum <- insect %>% group_by(insecticide) %>% summarise(n = n())
### plots for Country and insecticide summary
cols3 <- unikn::usecol(c("#CDD3D5", "#75B8C8", "#197278", "#03045E", "#922D50"), n = 5)


ggplot(data = insect_sum, aes(x = "", y = n, fill = insecticide)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols3) +
  theme_void() +
  theme(legend.title = element_text(size = 30), 
        legend.text = element_text(size = 28),
        legend.position = "right") +
  labs(fill = "Insecticide")+
  geom_text(aes(x = 1.6,y = n, label = n), position = position_stack(vjust = .5), color = "black", size = 16)




## translation
df2 <- read.delim("/Users/lsh1513859/OneDrive - London School of Hygiene and Tropical Medicine/vector_genomics2/PuertoRico/Output_10_01_24_amplicon_combine_no_samclip2/combined_genotyped_filtered_formatted.snps.trans_norm.txt", header = FALSE)
colnames(df2) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELE", "DP", "AD", "CSQ1")
length(unique(df2$SAMPLE))
length(unique(df2$POS))
# will match to other data frame so doesnt matter that there are extra samples in this one
# will only keep samples in first df

df2_sub <- df2
df2_sub <- df2_sub[paste0(df2_sub$SAMPLE, df2_sub$POS) %in% paste0(df_sub$SAMPLE, df_sub$POS),]

## check does this match number of samples and SNPs after filtering
# should be 178 and 57
length(unique(df2_sub$SAMPLE))
# 178
length(unique(df2_sub$POS))
# 57

## read in normalised non-translated file 
# subset based on filtering above
df3 <- read.delim("/Users/lsh1513859/OneDrive - London School of Hygiene and Tropical Medicine/vector_genomics2/PuertoRico/Output_10_01_24_amplicon_combine_no_samclip2/combined_genotyped_filtered_formatted.snps.norm.txt", header = FALSE)
colnames(df3) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELE", "DP", "AD")

df3 <- df3[paste0(df3$SAMPLE, df3$POS) %in% paste0(df_sub$SAMPLE, df_sub$POS),]




## add consequence to the SNP file so you get reference calls and consequence in one file
df_grouped <- df3 %>% dplyr::group_by(CHROM, POS, GT, ALT) %>% dplyr::summarise(n = n())
df_grouped$CSQ1 <- df2_sub$CSQ1[match(paste(df_grouped$POS, df_grouped$GT, df_grouped$ALT), paste(df2_sub$POS, df2_sub$GT, df2_sub$ALT))]
#df_grouped$CSQ2 <- df2_sub$CSQ2[match(paste(df_grouped$POS, df_grouped$GT), paste(df2_sub$POS, df2_sub$GT))]



# synonymous positions
synon_pos <- df_grouped$POS[grep("synonymous", df_grouped$CSQ2)]
synon_pos2 <- df_grouped$POS[grep("synonymous", df_grouped$CSQ1)]
unique(c(synon_pos, synon_pos2))
synon <- df_grouped[df_grouped$POS %in% synon_pos2,]

synon_long <- synon %>% select(CHROM, POS, GT, ALT, n) %>% tidyr::spread(GT, n)
# make table nice for paper
synon_long[is.na(synon_long)] <- 0
synon_long[,4:6] <- apply(synon_long[,4:6],2, as.numeric)
synon_long$total <- rowSums(synon_long[,4:6])

synon_long <- subset(synon_long, POS %in% c)

# splice positions
splice_pos <- df_grouped$POS[grep("splice", df_grouped$CSQ2)]
splice_pos2 <- df_grouped$POS[grep("splice", df_grouped$CSQ1)]
unique(c(splice_pos, splice_pos2))

# select only missense mutations
missense_pos <- df_grouped$POS[grep("missense", df_grouped$CSQ2)]
missense_pos2 <- df_grouped$POS[grep("missense", df_grouped$CSQ1)]
unique(c(missense_pos, missense_pos2))
missense <- df_grouped[df_grouped$POS %in% missense_pos2,]

# only include SNPs called in more than 50 samples
missense_good <- missense %>% group_by(POS) %>% dplyr::summarise(sum = sum(n))
missense_good <- subset(missense_good, sum > 50)

# count number of each without reference included
csq_only <- subset(df2_sub, GT != "0/0")
a <- unique(csq_only$POS[grep("missense", csq_only$CSQ1)])
b <- unique(csq_only$POS[grep("splice", csq_only$CSQ1)])
c <- unique(csq_only$POS[grep("synonymous", csq_only$CSQ1)])

positions <- c(a,b, c)
unique(csq_only$POS[!(csq_only$POS %in% positions)])

test2 <- csq_only %>% group_by(POS, CSQ1) %>% summarise(n =n())

CSQ <- strsplit(missense$CSQ1,split=',', fixed=TRUE)

missense_long <- missense %>% select(CHROM, POS, GT, ALT, n) %>% tidyr::spread(GT, n)
# make table nice for paper
missense_long[is.na(missense_long)] <- 0
missense_long[,4:6] <- apply(missense_long[,4:6],2, as.numeric)
missense_long$total <- rowSums(missense_long[,4:6])

str_split(missense$CSQ1, "\\|")[[6]]

df_missense <- missense %>%
  mutate(allele_info = strsplit(missense$CSQ1, ",")) %>%
  rowwise() %>%
  mutate(info_6 = sapply(allele_info, function(x) unlist(strsplit(x[6], "\\|"))[2]))

## check which mutation are in which Countrys
df_missense <- subset(df, POS %in% missense_good$POS)
df_missense <- df_missense %>% 
  group_by(POS, GT, ALLELE, Country) %>% 
  summarise(n = n())


df_missense$IR <- NA
IR_pos <- c("41847790", "315939224", "31598763", "315998453", "316080722")
df_missense$IR[df_missense$POS %in% IR_pos] <- "Yes"
df_missense$IR[!df_missense$POS %in% IR_pos] <- "No"

##################### stopping point before plots 

library(viridis)
library(unikn)

stop2<- df_missense

df_missense$POS <- as.factor(df_missense$POS)

## bar chart of proportion of each SNP in each Country
ggplot() + 
  geom_bar(data = df_missense, aes(x = reorder(Country, IR), y = n, fill = GT), 
           position = "fill", stat = "identity") +
  xlab("Position") +
  ylab("Percentage") +
  geom_text(aes(label=df$CSQ1), vjust=0) +
  scale_fill_viridis(discrete = TRUE) +
  theme_classic() +
  facet_wrap(~POS) +
  theme(axis.text.x = element_text(angle = 90))

## proportion of each IR SNP in each Country

IR_missense <- df_missense[df_missense$IR == "Yes",]
IR_missense_total <- IR_missense %>% group_by(POS, Country) %>% summarise(total = sum(n))
IR_missense$total <- IR_missense_total$total[match(paste(IR_missense$POS, IR_missense$Country), paste(IR_missense_total$POS, IR_missense_total$Country))]
IR_missense$proportion <- (IR_missense$n/IR_missense$total)*100

source("/Users/lsh1513859/Documents/git/Vector_Genomics/Puerto_Rico/PR_colour_palette.R")
cols <- colours(number = 4)

ggplot(IR_missense, aes(x="", y=proportion, fill = GT)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = cols) +
  coord_polar("y", start=0) +
  facet_grid(cols = vars(POS), rows = vars(Country)) +
  theme_void() +
  theme(strip.text = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) +
  #ggtitle("Global proportion of samples with multiple insecticide resistance mutations") +
  guides(fill = guide_legend(title = "Genotype", ncol = 1)) #+
#theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))  #,panel.background = element_rect(fill = 'grey55')) +
#geom_text(aes(x = 1.6,y = sum, label =sum), position = position_stack(vjust = .5), color = "black", size = 4)


# Allele pie charts for 4 IR SNPs
#cols <- unikn::usecol(c("goldenrod2", "coral2", "darkmagenta", "dodgerblue"), n = 9)
cols <- colours(number = 9)
ggplot(IR_missense, aes(x="", y=proportion, fill = ALLELE)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = cols) +
  coord_polar("y", start=0) +
  facet_grid(cols = vars(POS), rows = vars(Country)) +
  theme_void() +
  theme(strip.text = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) +
  #ggtitle("Global proportion of samples with multiple insecticide resistance mutations") +
  guides(fill = guide_legend(title = "Alleles", ncol = 1)) #+
#theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))  #,panel.background = element_rect(fill = 'grey55')) +
#geom_text(aes(x = 1.6,y = sum, label =sum), position = position_stack(vjust = .5), color = "black", size = 4)

IR_missense2 <- IR_missense %>% dplyr::group_by(POS, ALLELE) %>% dplyr::summarise(n = sum(n)) %>% mutate(proportion = n/sum(n))
cols <- colours(number = 9)
ggplot(IR_missense2, aes(x="", y=proportion, fill = ALLELE)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = cols) +
  coord_polar("y", start=0) +
  facet_grid(cols = vars(POS)) +
  theme_void() +
  theme(strip.text = element_text(size = 16), legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), legend.position = "bottom", legend.key.size = unit(1, "cm")) +
  guides(fill = guide_legend(title = "Alleles", ncol = 9))

IR_missense3 <- IR_missense %>% dplyr::group_by(POS, GT) %>% dplyr::summarise(n = sum(n)) %>% mutate(proportion = n/sum(n))

ggplot() + 
  geom_bar(data = IR_missense3, aes(x = POS, y = proportion, fill = GT), stat = "identity") +
  scale_fill_manual(values = cols) +
  xlab("Position") +
  ylab("Proportion") +
  theme_classic() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Alleles", ncol = 9))


############## stopping point

stop3 <- IR_missense

### calculate allele frequency

IR_missense$ALT <- df$ALT[match(paste(IR_missense$POS, IR_missense$GT), paste(df$POS, df$GT))]

df_IR <- subset(df, POS %in% IR_missense$POS)

#df_AF <- df_IR %>% select(POS, REF, ALT, GT, ALLELE) %>% group_by(POS, REF, ALT, GT, ALLELE) %>% summarise(n = n())
#df_AF$allele1 <- NA
#df_AF$allele2 <- NA
#for (i in 1:nrow(df_AF)){
# df_AF$allele1[i] <- unlist(str_split(df_AF$ALLELE[i], "/"))[[1]]
# df_AF$allele2[i] <- unlist(str_split(df_AF$ALLELE[i], "/"))[[2]]
#}

# Extract allele1 and allele2 from the ALLELE column
df_AF <- df_IR %>%
  select(POS, REF, ALT, GT, ALLELE) %>%
  group_by(POS, REF, ALT, ALLELE) %>%
  summarise(n = n()) %>%
  mutate(allele1 = str_split(ALLELE, "/") %>% sapply(function(x) x[1]),
         allele2 = str_split(ALLELE, "/") %>% sapply(function(x) x[2]))
df_AF$AF <- NA
df_AF$total <- NA


## Start loop to calculate allele frequency for each position
for ( i in unique(df_IR$POS)){
  
  df_IR_sub <- df_IR[df_IR$POS == i,]
  no_alleles <- length(unlist(str_split(df_IR_sub$ALT[1], ",")))
  alleles <- unlist(str_split(df_IR_sub$ALT[1], ","))
  df_IR_sub$allele1 <- NA
  df_IR_sub$allele2 <- NA
  
  for (j in 1:nrow(df_IR_sub)){
    df_IR_sub$allele1[j] <- unlist(str_split(df_IR_sub$ALLELE[j], "/"))[[1]]
    df_IR_sub$allele2[j] <- unlist(str_split(df_IR_sub$ALLELE[j], "/"))[[2]]
    
    alleles_to_match <- df_IR_sub$ALLELE[j]
    df_IR_sub$AF <- NA
    
    ## Need to calculate allele frequency by adding number of rows with each allele
    
  }
  for (k in unique(c(df_IR_sub$allele1, df_IR_sub$allele2))){
    
    AF <- (nrow(df_IR_sub[df_IR_sub$allele1 == k,])) + (nrow(df_IR_sub[df_IR_sub$allele2 == k,])) 
    total <- nrow(df_IR_sub)*2
    
    
    df_AF$total[(df_AF$POS == i & 
                   (df_AF$allele1 == k | df_AF$allele2 == k) & 
                   (df_AF$allele1 != df_AF$REF[1] | df_AF$allele2 != df_AF$REF[1]))] <- total
    
    df_AF$AF[(df_AF$POS == i & 
                (df_AF$allele1 == k | df_AF$allele2 == k) & 
                (df_AF$allele1 != df_AF$REF[1] | df_AF$allele2 != df_AF$REF[1]))] <- AF
    
  }
}



df_AF2 <- df_AF %>% dplyr::group_by(POS, REF, ALLELE, allele2) %>% reframe(AF_prop = AF/total)
df_AF2 <- df_AF2[!duplicated(df_AF2$AF_prop),]
test <- df_AF2 %>% group_by(POS) %>% summarise(sum = sum(AF_prop))

df_AF3 <- subset(df_AF2, allele2 != REF)
### do chi squared between groups

rdl <- subset(df_AF2, POS == "41847790")

table(rdl$Country, rdl$MAF)
chisq.test(rdl$Country, rdl$)


### first attempt AF
IR_missense$MAF <- NA
IR_missense$MAF_prop <- NA

for (i in unique(IR_missense$POS)){
  for (j in unique(IR_missense$Country)){
    
    IR_sub <- IR_missense[IR_missense$POS == i & IR_missense$Country == j,]
    
    ref <- 0
    ref_freq <- 0
    homo_alt <- 0
    homo_alt_freq <- 0
    hetero <- 0
    hetero_freq <- 0
    
    length(unique(str_split_i(unique(IR_sub$GT), "/", 1)))
    if (length(unq))  
      
      
      for (n in 1:nrow(IR_sub)){
        
        if(substring(IR_sub$GT[n], 1,1) == 0 & substring(IR_sub$GT[n], 3,3) == 0){
          ref <- IR_sub$n[n]
          ref_freq <- ref*2
        }
        if(substring(IR_sub$GT[n], 1,1) != 0 & substring(IR_sub$GT[n], 3,3) != 0 &
           substring(IR_sub$GT[n], 1,1) == substring(IR_sub$GT[n], 3,3)){
          homo_alt <- IR_sub$n[n]
          homo_alt_freq <- homo_alt*2
        }
        if(substring(IR_sub$GT[n], 1,1) != substring(IR_sub$GT[n], 3,3)){
          hetero <- IR_sub$n[n]
          hetero_freq <- hetero*2
        } 
      }
    if (exists("hetero_freq") & exists("homo_alt_freq") & exists("ref_freq")) {
      MAF <- (hetero + homo_alt_freq)
      MAF_prop <- ((hetero + homo_alt_freq) / (hetero_freq + homo_alt_freq + ref_freq))  
      MAF2 <- ((hetero + homo_alt_freq) / (IR_sub$total[n]*2))
    }
    # if (exists("hetero_freq") & exists("homo_alt_freq")) {
    # MAF <- ((hetero_freq + homo_alt_freq) / (hetero_freq + homo_alt_freq))             
    # }
    # if (exists("hetero_freq") & exists("ref_freq")) {
    #  MAF <- ((hetero_freq ) / (hetero_freq + ref_freq))           
    # }
    #if (exists("homo_alt_freq") & exists("ref_freq")) {
    # MAF <- ((homo_alt_freq) / (homo_alt_freq + ref_freq))         
    # }
    print(i)
    print(j)
    print(MAF)
    print(MAF2)
    
    IR_missense$MAF[IR_missense$POS == i & IR_missense$Country == j] <- MAF
    IR_missense$MAF_prop[IR_missense$POS == i & IR_missense$Country == j] <- MAF_prop
  }
  
  
}



### do chi squared between groups

rdl <- subset(IR_missense, POS == "41847790")

table(rdl$Country, rdl$MAF)
chisq.test(rdl$Country, rdl$)

# compare Country and mortality % with both insecticides

bio <- read.csv("../../../OneDrive - London School of Hygiene and Tropical Medicine/Puerto_Rico_field_work/Results/IR_results_r_2.csv")


data_frame <- read.csv("https://goo.gl/j6lRXD")  #Reading CSV
table(data_frame$treatment, data_frame$improvement)

delta_1 <- subset(bio, insecticide == "Deltamethrin")
table(delta_1$Country, delta_1$concentration)
chisq.test(delta_1$Country, delta_1$mortality)

mal_1 <- subset(bio, concentration == 1 & insecticide == "Malathion")
chisq.test(mal_1$Country, mal_1$mortality)

# compare minor allele frequency and Country


########################## Synonymous SNPs #############################################

synon_pos <- df_grouped$POS[grep("synonymous", df_grouped$CSQ2)]
synon_pos <- c(synon_pos, df_grouped$POS[grep("synonymous", df_grouped$CSQ1)])
synon <- df_grouped[df_grouped$POS %in% synon_pos,]

# only include SNPs called in more than 50 samples
synon_good <- synon %>% group_by(POS) %>% dplyr::summarise(sum = sum(n))
synon_good <- subset(synon_good, sum > 50)

CSQ <- strsplit(synon$CSQ1,split=',', fixed=TRUE)
CSQ <- strsplit(synon$CSQ2,split=',', fixed=TRUE)

synon_wide <- synon[,1:4]
synon_wide <- spread(synon_wide, GT, n)



## meteadata - counts for each category
meta <- read.csv("../../../OneDrive - London School of Hygiene and Tropical Medicine/vector_genomics2/Reference_files/sample_metadata.csv")

meta_sub <- subset(meta, sample %in% good_list2$SAMPLE)
meta_sub <- meta

meta_sub_group <- meta_sub %>% group_by(insecticide) %>% summarise(count = n())

ggplot() +
  geom_bar(data = meta_sub_group, aes(x = insecticide, y = count), stat = "identity")

ggplot(data = meta_sub_group, aes(y = count, x = insecticide, fill = insecticide)) +
  geom_bar(stat = "identity", position = "stack")

# pie chart of insecticide exposure
cols <- unikn::usecol(c("goldenrod2", "coral2", "darkmagenta", "darkgreen", "darkorchid3", "dodgerblue"), n = 8)

ggplot(meta_sub_group, aes(x="", y=count, fill = insecticide)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = cols) +
  coord_polar("y", start=0) +
  theme_void() +
  guides(fill=guide_legend(title="Insecticide")) +
  theme(legend.position = "right", legend.key.size = unit(1, "cm"), legend.text = element_text(size = 32), legend.title = element_text(size = 32)) + #,panel.background = element_rect(fill = 'grey55')) +
  geom_text(aes(x = 1.6,y = count, label = count), position = position_stack(vjust = .5), color = "black", size = 16)


meta_sub_group2 <- meta_sub %>% group_by(Country) %>% summarise(count = n())

plot <- ggplot() +
  geom_bar(data = meta_sub_group2, aes(x = reorder(Country, count), y = count), fill = "grey65", stat = "identity") +
  ylab("Frequency") +
  xlab("Country") +
  theme_minimal()


plot + theme(axis.text = element_text(size = 24), 
             axis.title = element_text(size = 24))

df2_sub$Country <- meta$Country[match(df2_sub$SAMPLE, meta$sample)]


## SNPs by Country

missense_loc <- subset(df, POS %in% missense_pos | POS %in% synon_pos)
loc_snp <- missense_loc %>% group_by(POS, GT, Country) %>% dplyr::summarise(n= n())
loc_snp_wide <- spread(loc_snp, GT, n)
loc_snp_wide$`0/0`[is.na(loc_snp_wide$`0/0`)] <- 0
loc_snp_wide$`0/1`[is.na(loc_snp_wide$`0/1`)] <- 0
loc_snp_wide$`1/1`[is.na(loc_snp_wide$`1/1`)] <- 0
loc_snp_wide$`1/2`[is.na(loc_snp_wide$`1/2`)] <- 0
#loc_snp_wide$total <- loc_snp_wide$`0/0` + loc_snp_wide$`0/1` +loc_snp_wide$`1/1` + loc_snp_wide$`1/2`

loc_snp_wide$POS <- as.factor(loc_snp_wide$POS)

## Number of each GT for each SNP
snp_counts <- missense_loc %>% group_by(POS, GT) %>% dplyr::summarise(n=n())
snp_counts <- spread(snp_counts, GT, n)

## proportion of GT for each Country
loc_snp_wide2 <- loc_snp_wide %>% 
  group_by(POS, Country) %>% 
  mutate(sum = as.character(rowSums(select(cur_data(), is.numeric)))) %>%
  summarise_if(is.numeric, ~ . / as.numeric(sum))

# which missense SNPs found in all Countrys
loc_snp2 <- missense_loc %>% group_by(POS, Country) %>% dplyr::summarise(n= n())
loc_snp3 <- loc_snp2 %>% group_by(POS) %>% dplyr::summarise(n= n())



