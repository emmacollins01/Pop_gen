
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

## Add metadata
metadata <- read.csv("/mnt/storage12/emma/PR_combine/all_sample_country_metadata.csv")
df$Country <- metadata$Country[match(df$SAMPLE, metadata$ID)]
head(df)
## Add snp locations
SNP_locs <- read.delim("/mnt/storage12/emma/PR_snps/IR_mut_pos.txt", header = TRUE)

# snps per country
country_snp <- df %>% group_by(Country) %>% summarise(n = n_distinct(CHROM,POS))

unique(df$POS)

# IR snps per country
df_mut <- df[df$POS %in% SNP_locs$POS,]
for (i in 1:nrow(df_mut)){
  df_mut$effect[i] <- unlist(strsplit(df_mut$ANN[i], "\\|"))[[2]]
}
df_mut <- df_mut[df_mut$effect != "synonymous_variant",]
unique(df_mut$POS)

# number of IR mutations per sample 0-4
df_mut_sample <- df_mut %>% filter(GT != "0/0") %>% dplyr::group_by(SAMPLE, Country) %>% summarise(count = n_distinct(POS))
# number of IR mutations per country 0-4
df_mut_country <- df_mut %>% filter(GT != "0/0") %>% group_by(Country) %>% summarise(count = n_distinct(POS))

df_mut$csq <- NA
df_mut$csq <- SNP_locs$Ref_Mutation_Pos[match(df_mut$POS, SNP_locs$POS)]

df_mut$csq_bin <- NA


for (i in 1:nrow(df_mut)){
  if(df_mut$GT[i] == "0/0"){
    df_mut$csq_bin[i] <- 0
} 
  if(df_mut$GT[i] != "0/0") {
    df_mut$csq_bin[i] <- 1
  }
}


unique(df_mut$csq_bin)

pie <- df_mut %>% select(SAMPLE, Country, csq, csq_bin)
wide <- spread(pie, csq, csq_bin)
wide$A302S[is.na(wide$A302S)] <- 0
wide$F1534C[is.na(wide$F1534C)] <- 0
wide$V1016I[is.na(wide$V1016I)] <- 0
wide$V410L[is.na(wide$V410L)] <- 0

long <- gather(wide, snp, presence, A302S:V410L)

long_sum <- long %>% group_by(SAMPLE,Country) %>% reframe(sum = sum(presence))
long_sum2 <- long_sum %>% group_by(Country, sum) %>% reframe(n=n())
long_sum2$sum <- as.factor(long_sum2$sum)

seq <- c(0,1,2,3,4)
new_df <- data.frame("Country" = rep(unique(long_sum2$Country), 5), "sum" = rep(seq, 8))
new_df$n <- long_sum2$n[match(paste0(new_df$Country, new_df$sum), paste0(long_sum2$Country, long_sum2$sum))]
new_df$n[is.na(new_df$n)] <- 0

long_sum2 <- new_df

ftable(crossTable)  # print crosstabs
summary(crossTable) # chi-square tests)

#cols <- c("#40476D", "#EE964B", "#749cb4",  "#A6C48A",  "#51A3A3")
#cols <- c("#40476D", "#ffaa5c","#bc80bd","#51A3A3", "#7D387D", "#bc80bd")
cols <- c("#40476D", "#749cb4",  "#A6C48A",  "#51A3A3","#EE964B" )
#cols <- c("#ffaa5c", "#db7376", "#52a884", "#bc80bd", "#465c7a")

plot_list <- list()
for (i in unique(long_sum2$Country)){
  
  plot <- ggplot(long_sum2[long_sum2$Country == i,], aes(x="", y=n, fill = factor(sum))) +
    geom_bar(stat="identity", width=1) +
    #scale_fill_manual(values = cols) +
    scale_fill_manual(values = c("#40476D", "#749cb4",  "#A6C48A", "#51A3A3","#EE964B" ), breaks = c(0,1,2,3,4)) +
    coord_polar("y", start=0) +
    facet_wrap(~Country) +
    theme_void() +
    guides(fill = guide_legend(title = "No. of insecticide \nresistance mutations")) +
    theme(legend.position = "none", strip.text.x = element_text(size = 30)) +
    geom_text(aes(x = 1.8,y = n, label =n), size = 10, position = position_stack(vjust = .5), color = "black")
  
  print(plot)
  plot_list[[i]] <- plot
}

plot_list[[1]]


library(patchwork)

combined_plot <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] +
    plot_list[[5]] + plot_list[[6]] + plot_list[[7]] + plot_list[[8]] + plot_layout(ncol = 2)
combined_plot


## Just for legend
ggplot(long_sum2[long_sum2$Country == "Burkina Faso",], aes(x="", y=n, fill = factor(sum))) +
    geom_bar(stat="identity", width=1) +
    scale_fill_manual(values = c("#40476D", "#749cb4",  "#A6C48A", "#51A3A3","#EE964B" ), breaks = c(0,1,2,3,4)) +
    coord_polar("y", start=0) +
    #facet_wrap(~Country) +
    #theme_void() +
    guides(fill = guide_legend(title = "No. of \ninsecticide \nresistance \nmutations")) +
    theme(legend.position = "bottom", 
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.key.size = unit(3, "lines"), 
    strip.text.x = element_text(size = 30)) +
    geom_text(aes(x = 1.8,y = n, label =n), size = 10, position = position_stack(vjust = .5), color = "black")
  



