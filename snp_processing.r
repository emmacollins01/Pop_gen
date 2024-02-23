
setwd("PR_snps")

library(vcfR)
library(tidyr)
library(dplyr)

vcf <- read.vcfR("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode_IR_genes_ann.vcf.gz")

metadata <- read.csv("../PR_combine/all_sample_country_metadata.csv")
metadata <- as.data.frame(metadata)

gt <- as.data.frame(extract.gt(vcf, element = "GT"))
gt_snp <- gt[row.names(gt) == "NC_035108.1_41847790",]

gt_snp_long <- pivot_longer(gt_snp, cols = BAY_F1:SRR11196652)

gt_snp_long <- as.data.frame(gt_snp_long)

gt_snp_long$country <- metadata$Country[match(gt_snp_long$name, metadata$ID)] 


rdl_sum <- gt_snp_long %>% group_by(value, country) %>% dplyr::summarise(n = n())

rdl_sum <- rdl_sum[order(rdl_sum$country),]
