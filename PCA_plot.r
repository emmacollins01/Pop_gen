
##############################################################
#####################  Read in Packages  #####################
##############################################################

options(scipen=999) #this ensures analysis is in scientific notation
#install.packages("data.table")
require(data.table) #allows you to read in dataframes
require(ape) #popgen package
require(qqman) #popgen package
library(ggplot2) #for plotting
require("colorRamps") #for colouring
library("RColorBrewer") #colouring
#library(tidyverse)

#setwd("~/Desktop/pop_gen/")
## run on server!!
## cmd + option + enter

##############################################################
#######################  Read in Data  #######################
##############################################################

#mat_bin <- read.delim("emmadata.2024_02_02.genotyped.norm.mat.bin")
#mat_bin <- read.delim("/mnt/storage12/emma/PR_combine/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.mat.bin")
#mat_bin <- read.delim("/mnt/storage12/emma/PR_combine/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.indv_miss0.5.recode.mat.bin")
mat_bin <- read.delim("/mnt/storage12/emma/PR_combine/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.mat.bin")
mat_bin <- read.delim("/mnt/storage12/emma/consensus/PR_only_WGS.mat.bin")
mat_bin <- as.data.frame(mat_bin)

# Read in metadata
metadata<-read.delim("/mnt/storage12/emma/PR_combine/all_sample_country_metadata.csv", sep=",", header=T)
metadata<-as.data.frame(metadata)

### You want to check that both the SNP matrix and the metadata have the same samples ###
metadata <- metadata[(metadata$ID) %in% colnames(mat_bin),]
mat_bin_samples <- mat_bin[,colnames(mat_bin )%in%metadata$ID]



##############################################################
########################  Check NAs  #########################
##############################################################

na_count <-sapply(snp_na, function(y) sum(length(which(is.na(y)))))
na_count<-as.data.frame(na_count)
na_count<- na_count[4:nrow(na_count),]
na_count<-cbind(metadata$ID, metadata$Country, na_count)
na_count<-as.data.frame(na_count)

library(ggplot2)
plot <- ggplot() +
  geom_boxplot(data = na_count, aes(x = as.factor(V2), y = as.numeric(na_count), fill = V2)) +
  xlab("")+
  ylab("NA count") +
  ggtitle("all_aedes_norm_filt_miss0.5_mac3_minQ30") +
  theme_minimal()
plot
ggsave("/mnt/storage12/emma/PR_combine/all_aedes_norm_filt_miss0.5_mac3_minQ30_NA_check.png", plot)
dev.off()

rm(na_count)

##############################################################
##################  Imputing NAs option  #####################
##############################################################
snp_imputed_freq<-t(snp_na)

for(i in 1:ncol(snp_imputed_freq)){
  snp_imputed_freq[is.na(snp_imputed_freq[,i]), i] <- as.numeric(names(which.max(table(snp_imputed_freq[,i]))))
}

snp_imputed_freq_t<-t(snp_imputed_freq)
dim(snp_imputed_freq_t)
write.csv(snp_imputed_freq_t, "all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.mat.bin.imput", row.names = FALSE)

snp_imputed_freq_t <- read.delim("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.mat.bin.imput", sep = ",", header = TRUE)
dim(snp_imputed_freq_t)
## if has indexed row names need to remove
snp_imputed_freq_t <- snp_imputed_freq_t[,2:ncol(snp_imputed_freq_t)]

na_count_imp <-sapply(snp_imputed_freq_t, function(y) sum(length(which(is.na(y)))))
na_count_imp<-as.data.frame(na_count_imp)
write.csv(na_count_imp, "na_count_imp_df")

#na_count_imp<- na_count_imp[4:nrow(na_count_imp),]
na_count_imp<-cbind(metadata$ID, metadata$Country, na_count_imp)
na_count_imp<-as.data.frame(na_count_imp)
colnames(na_count_imp) <- c("ID", "Country", "count")
rownames(na_count_imp) <- NULL

library(ggplot2)
plot <- ggplot() +
  geom_point(data = na_count_imp, aes(x = Country, y = count))
ggsave("NA_check_imp.png", plot)
dev.off()

##############################################################
###############  Creating the distance matrix  ###############
##############################################################

PR_meta <- read.csv("/mnt/storage12/emma/PR_combine/all_samples_PR_resistance.csv")
subset_SRR <- unique(PR_meta$ID) 

subset <- mat_bin_samples[,colnames(mat_bin_samples) %in% subset_SRR]

dist_dat<-dist(t(mat_bin_samples), method = "manhattan")
dist_dat<-dist(t(mat_bin), method = "manhattan")
dist_dat_PR<-dist(t(subset), method = "manhattan")
#dist_dat<-dist(t(snp_imputed_freq_t), method = "manhattan")
#write.csv(dist_dat, "all_aedes_norm_ filt_miss0.5_mac3_minQ30.vcf.gz.recode.dist.mat")

#Formatting for a PCA
cmd_dat<-cmdscale(dist_dat, k = 10, eig = TRUE, x.ret = TRUE)
cmd_dat_PR<-cmdscale(dist_dat_PR, k = 10, eig = TRUE, x.ret = TRUE)
#write.csv(cmd_dat, "all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.dist.mat.pca")


# work out eigenvalue for axis label
calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

vars <- calc_variance_explained(cmd_dat)
vars_PR <- calc_variance_explained(cmd_dat_PR)
PC1_lab <- vars["PC1"]
PC2_lab <- vars["PC2"]
PC3_lab <- vars["PC3"]

PC1_lab <- vars_PR["PC1"]
PC2_lab <- vars_PR["PC2"]
PC3_lab <- vars_PR["PC3"]

##############################################################
##########################  PCA plot  ########################
##############################################################


df <- as.data.frame(cmd_dat$points, stringsAsFactors = F)
colnames(df) <- gsub("V", "PC", colnames(df))
df$Country <- metadata$Country
df$Region <- metadata$Region
df$Subregion <- metadata$Subregion

# filtered version
library(stringr)
df_PR <- as.data.frame(cmd_dat_PR$points, stringsAsFactors = F)
colnames(df_PR) <- gsub("V", "PC", colnames(df_PR))
df_PR$Country <- PR_meta$Country
df_PR$Region <- PR_meta$Region
df_PR$Subregion <- PR_meta$Subregion
for (i in 1:nrow(df_PR)){
    df_PR$Subregion_simple[i] <- unlist(strsplit(df_PR$Subregion[i], ":"))[[2]]
}



#DEFINE COLOUR TO USE
#col_ramp <- brewer.pal(8, "Pastel1")
col_ramp <- c("#40476D", "#678D58","#EE964B", "#F95738",  "#A6C48A",  "#51A3A3", "#7D387D", "#bc80bd")
#col_ramp <- c("#bc80bd", "#EE964B", "#F95738", "#7D387D", "#40476D", "#678D58", "#A6C48A", "#51A3A3")
#col_ramp <- c("#db7376", "#52a884", "#bc80bd", "#465c7a", "#ffaa5c",  "#A6C48A",  "#51A3A3", "#7D387D")
colour_by <- "Country"
colour_by <- "Subregion"
colour_by <- "Subregion_simple"


## change data input if using filtered version
####### PC1 vs PC2 #######
plot <- ggplot() + 
    geom_point(data = df, aes(x = PC1, y = PC2, colour = !!sym(colour_by)), size = 5) +
    xlab(paste0("PC1 ", PC1_lab, "%")) +
    ylab(paste0("PC2 ", PC2_lab, "%")) +
    scale_colour_manual(values = col_ramp, name = "Subregion") +
    theme_classic() +
    theme(legend.position = "bottom",
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 24),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24))

plot
tiff("/mnt/storage12/emma/PR_combine/PCA_plots/PCA_all_samples_PC1_PC2_lmiss.tiff",height=1200,width=1200)
dev.off()
tiff("/mnt/storage12/emma/PR_combine/PCA_plots/PCA_PR_samples_PC1_PC2_lmiss.tiff",height=1200,width=1200)


####### PC1 vs PC3 #######
plot <- ggplot() + 
    geom_point(data = df_PR, aes(x = PC1, y = PC3, colour = !!sym(colour_by)), size = 4) +
    xlab(paste0("PC1 ", PC1_lab, "%")) +
    ylab(paste0("PC3 ", PC3_lab, "%")) +
    scale_colour_manual(values = col_ramp) +
    theme_classic() +
    theme(legend.position = "bottom",
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 24),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24))

plot
tiff("/mnt/storage12/emma/PR_combine/PCA_plots/PCA_all_samples_PC1_PC3_lmiss.tiff",height=1200,width=1200)
dev.off()


####### PC2 vs PC3 #######
plot <- ggplot() + 
    geom_point(data = df, aes(x = PC2, y = PC3, colour = !!sym(colour_by)), size = 4) +
    xlab(paste0("PC2 ", PC2_lab, "%")) +
    ylab(paste0("PC3 ", PC3_lab, "%")) +
    scale_colour_manual(values = col_ramp) +
    theme_classic() +
    theme(legend.position = "bottom",
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 24),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24))

plot
tiff("/mnt/storage12/emma/PR_combine/PCA_plots/PCA_all_samples_PC2_PC3_lmiss.tiff",height=1200,width=1200)
dev.off()



###### test remove sample ########
which.max(df$PC2)
dim(df)
df_rm <- df[c(1:183, 185:248),]
plot <- ggplot() + 
    geom_point(data = df_rm, aes(x = PC1, y = PC2, colour = !!sym(colour_by)), size = 4) +
    xlab(paste0("PC1 ", PC1_lab, "%")) +
    ylab(paste0("PC2 ", PC2_lab, "%")) +
    scale_colour_manual(values = col_ramp) +
    theme_classic() +
    theme(legend.position = "bottom",
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 24),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24))

plot

which.max(df_rm$PC2)
dim(df_rm)
df_rm <- df_rm[c(1:111,113:245),]

####### PC1 vs PC2 #######
plot <- ggplot() + 
    geom_point(data = df_rm, aes(x = PC1, y = PC2, colour = !!sym(colour_by)), size = 4) +
    xlab(paste0("PC1 ", PC1_lab, "%")) +
    ylab(paste0("PC2 ", PC2_lab, "%")) +
    scale_colour_manual(values = col_ramp) +
    theme_classic() +
    theme(legend.position = "bottom",
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 24),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24))

plot
tiff("/mnt/storage12/emma/PR_combine/PCA_plots/PCA_all_samples_PC1_PC2_lmiss_no_outlier.tiff",height=1200,width=1200)
dev.off()

####### PC1 vs PC3 #######
plot <- ggplot() + 
    geom_point(data = df_rm, aes(x = PC1, y = PC3, colour = !!sym(colour_by)), size = 4) +
    xlab(paste0("PC1 ", PC1_lab, "%")) +
    ylab(paste0("PC3 ", PC3_lab, "%")) +
    scale_colour_manual(values = col_ramp) +
    theme_classic() +
    theme(legend.position = "bottom",
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 24),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24))

plot
tiff("/mnt/storage12/emma/PR_combine/PCA_plots/PCA_all_samples_PC1_PC3_lmiss_no_outlier.tiff",height=1200,width=1200)
dev.off()

####### PC2 vs PC3 #######
plot <- ggplot() + 
    geom_point(data = df_rm, aes(x = PC2, y = PC3, colour = !!sym(colour_by)), size = 4) +
    xlab(paste0("PC2 ", PC2_lab, "%")) +
    ylab(paste0("PC3 ", PC3_lab, "%")) +
    scale_colour_manual(values = col_ramp) +
    theme_classic() +
    theme(legend.position = "bottom",
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 24),
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24))

plot
tiff("/mnt/storage12/emma/PR_combine/PCA_plots/PCA_all_samples_PC2_PC3_lmiss_no_outlier.tiff",height=1200,width=1200)
dev.off()


#### NJ TREE ####
require("ape")

tree_dat<-nj(dist_dat)
write.tree(tree_dat, file="all_samples_lmiss.tree")


tree_dat<-nj(dist_dat_PR)
write.tree(tree_dat, file="PR_lmiss.tree")

## open in iTOL








##old plot script

plot(cmd_dat$points[,1], cmd_dat$points[,2],pch=16,xlab="PC1",ylab="PC2",cex=2.5,cex.axis=2)
plot(cmd_dat$points[,1], 
        cmd_dat$points[,2],
        pch=16, 
        col=col_ramp[as.factor(metadata$Country)],
        xlab=paste0("PC1 ", PC1_lab, "%"),
        ylab="PC2",cex=2.5,cex.axis=2)

legend("bottomleft",title="Legend",legend=unique(metadata_rem$Country),cex=1,col=col_ramp[unique(as.factor(metadata_rem$Country))],pch=16)

#TO SAVE PC1 vs PC2
tiff("./PCA_all_samples_PC1_PC2_imput.tiff",height=1200,width=1200)
plot(cmd_dat$points[,1], cmd_dat$points[,2], col = col_ramp[as.factor(metadata$Country)],pch=16,xlab="PC1",ylab="PC2",cex=2.5,cex.axis=2)
legend("bottomleft",title="Legend",legend=unique(metadata$Country),cex=1,col=col_ramp[unique(as.factor(metadata$Country))],pch=16)
dev.off()

#TO SAVE PC1 vs PC2
tiff("./PCA_all_samples_PC1_PC3_imput2.tiff",height=1200,width=1200)
plot(cmd_dat$points[,1], cmd_dat$points[,3], col = as.factor(metadata$Country),pch=16,xlab="PC1",ylab="PC3",cex=2.5,cex.axis=2)
legend("bottomright",title="Legend",col = unique(as.factor(metadata$Country)), legend=unique(metadata$Country),cex=1.5,pch=16)
dev.off()

#TO SAVE PC2 vs PC3
tiff("./PCA_all_samples_PC2_PC3_imput2.tiff",height=1200,width=1200)
plot(cmd_dat$points[,2], cmd_dat$points[,3], col = as.factor(metadata$Country),pch=16,xlab="PC2",ylab="PC3",cex=2.5,cex.axis=2)
legend("bottomright",title="Legend",col = unique(as.factor(metadata$Country)), legend=unique(metadata$Country),cex=1.5,pch=16)
dev.off()

## test remove sample
cmd_dat_df <- as.data.frame(cmd_dat$points)

which.max(cmd_dat_df[,1])
dim(cmd_dat_df)
cmd_dat_df <- cmd_dat_df[c(1:130, 132:248),]
which.max(cmd_dat_df[,1])
cmd_dat_df <- cmd_dat_df[c(1:137, 139:247),]

calc_variance_explained <- function(pc_points) {
  vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
  names(vars) <- paste0("PC", seq_len(length(vars)))
  vars
}

vars <- calc_variance_explained(cmd_dat)
library(ggplot2)
tiff("./PCA_all_samples_PC1_PC2_test2.tiff",height=1200,width=1200)
plot(cmd_dat_df$V1, cmd_dat_df$V2, col = col_ramp[as.factor(metadata$Country)],pch=16,xlab="PC1",ylab="PC3",cex=2.5,cex.axis=2)
dev.off()

ggplot(data = cmd_dat_df, aes(x = V1, y = V2)) +
  geom_point(aes(size = 10, colour = as.factor(metadata$Country))) +
  labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
       y = paste0("PC2", " (", vars["PC2"], "%)")) +
  theme_classic() +
  theme(legend.position = "bottom")

dev.off()


#### NJ TREE ####
workdir <- "/mnt/storage12/emma/NJ_tree_files" # Working directory with plink files
prefix <- "PR_WGS" # Prefix for plink files
metadata <- "/mnt/storage12/emma/PR_combine/all_sample_country_metadata.csv" # File path to metadata

#### DIST#
dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))
met <- read.table(metadata, sep = ",", stringsAsFactors = FALSE, header = TRUE)

desc <- id %>% left_join(met, by = c("V1" = "ID"))

dist_m <- as.matrix(dist)

# Export dist_m to .newick to make neighbour joining tree
tree <- nj(dist_m)
write.tree(phy = tree, file = file.path(workdir, paste0(prefix, ".newick.tree")))

require("ape")

tree_dat<-nj(dist_dat)
write.tree(tree_dat, file="/mnt/storage12/emma/NJ_tree/PR_WGS.tree")

## open in iTOL


