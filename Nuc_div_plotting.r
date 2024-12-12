
## Start by making file of samples names per population
# then file with all of these files eg/ ls SRR.txt > countries.txt
# then you can run separate_vcfs.sh to split multisample vcf into separate populations
# OR 
#for i in $(ls *SRR.txt); do bcftools view -S $i all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.vcf.gz -Oz -o ${i}.vcf.gz; done

## to make nucleotide diversity metrics
## Run vcftools - run on 10kb windows
## These files are in input here
#for i in $(ls *SRR.vcf.gz);do vcftools --gzvcf $i --window-pi 10000 --out ${i}_nuc_div_window_10kb;done
library(dplyr)
setwd("/mnt/storage12/emma/nuc_div/")
rm(list = ls())
nuc <- read.csv("all_nuc_div_window_10kb.windowed.pi", sep = "\t", header = TRUE)

x <- nuc %>% slice_max(PI, n = 20)
rm(nuc, x)


files <- list.files(pattern = "*z_nuc_div_window_10kb.windowed.pi")
files <- list.files(pattern = "*z_nuc_div_window_20kb.windowed.pi")
files <- list.files(pattern = "*z_nuc_div_window_100kb.windowed.pi")
files <- list.files(pattern = "*z_nuc_div_window_300kb.windowed.pi")
files <- list.files(pattern = "*z_nuc_div_window_500kb.windowed.pi")
files <- list.files(pattern = "*z_nuc_div_window_1000kb.windowed.pi")
files <- list.files(pattern = "*z_nuc_div_window_5000kb.windowed.pi")

countries <- c()
dataframes <- c()

# loop through the files
for (i in 1:length(files)){
  # save the clean filename as a char var because it will be called again
  #fnClean = substr(files[i], start = 1, stop = nchar(files[i])-10)
  ## CHANGE THIS LINE DEPENDING ON HOW FILES ARE NAMED
  fnClean = unlist(strsplit(files[i], "_subset_SRR\\."))[1]
  dataframes <- c(dataframes, fnClean)
  # create a cmd as a string to be parsed and evaluated on-the-fly
  # the point here is that you can use the 'fnClean' var in the string
  # without knowing what it is - assign is expecting an inline string
  # ...not a string saved as a var, so it can't be used in this case
  loadFileCMD = paste0(fnClean,' = read.delim(files[i], stringsAsFactors = FALSE)')
  print(loadFileCMD) # check the cmd
  eval(parse(text=loadFileCMD))
  
  # create another string command to be evaluated to insert the file name
  # to the 'FileFrom' field
  countries <- c(countries, fnClean)
  if(nrow(eval(parse(text=loadFileCMD)) >0)){
   ## CHANGE THIS LINE DEPENDING ON HOW FILES ARE NAMED - WLL GO ON COUNTRY_NAME COLUMN
    country_name <- unlist(strsplit(files[i], "_subset_SRR\\."))[1]
    addFnCMD = paste0(fnClean,'$FileFrom = "', country_name, '"')
    print(addFnCMD) # check the cmd
    eval(parse(text=addFnCMD)) }
}

# check added column for country
head(Kenya)

table(all$country_name, all$PI)

## combine all dataframes together - now they have country column
dataframes <- Filter(function(x) is(x, "data.frame"), mget(ls()))
all <- do.call(rbind, dataframes)
head(all)

summary(all$PI)

library(ggplot2)

## Subset and change names of chromosomes if needed!
nuc <- subset(all, CHROM %in% c(35107, 35108, 35109, 35159))
nuc$CHROM[nuc$CHROM == 35107] <- "1"
nuc$CHROM[nuc$CHROM == 35108] <- "2"
nuc$CHROM[nuc$CHROM == 35109] <- "3"
nuc$CHROM[nuc$CHROM == 35159] <- "MT"

nuc$FileFrom[nuc$FileFrom == "Burkina_Faso"] <- "Burkina Faso"
nuc$FileFrom[nuc$FileFrom == "Puerto_Rico"] <- "Puerto Rico"
colnames(nuc)
options(scipen = 0)
options(scipen = 999)
library(unikn)
mycols <-  unikn::usecol(c("#ffaa5c", "#db7376", "#52a884", "#bc80bd", "#465c7a"), n = 8)

summary <- nuc %>% group_by(FileFrom, CHROM) %>% summarise(mean = mean(PI))
summary
summary2 <- nuc %>% group_by(FileFrom) %>% summarise(mean = mean(PI))
summary3 <- nuc %>% group_by(CHROM) %>% summarise(mean = mean(PI))

summary(nuc$PI)

pr_plot <- ggplot() + 
geom_line(data = nuc[nuc$FileFrom == "Puerto Rico" & nuc$CHROM != "MT",], aes(x = BIN_START, y = PI), alpha = 0.8, colour = "#709C94", size = 1.2) +
facet_wrap(~CHROM, scale = "free_x", nrow = 3) +
xlab("Position") +
ylab("Nucleotide Diversity (π)") +
theme_classic() +
theme(axis.text = element_text(size = 22),
axis.title = element_text(size = 22),
strip.text = element_text(size =22))

x <- nuc %>% filter(FileFrom == "Puerto Rico") %>% slice_max(PI, n = 20)
x
options(scipen)
all_plot <- ggplot() + 
geom_line(data = nuc[nuc$CHROM != "MT",], aes(x = BIN_START, y = PI, colour = FileFrom), alpha = 0.8, size = 1.2) +
facet_wrap(~CHROM, scale = "free_x", nrow = 3) +
scale_colour_manual(values = mycols) +
xlab("Position") +
ylab("Nucleotide Diversity (π)") +
theme_classic() +
guides(color = guide_legend(override.aes = list(size=5))) +
theme(
legend.title = element_blank(),
axis.title = element_text(size = 22),
strip.background = element_rect(fill = "#a8bac3"),
legend.position = "bottom",
axis.text = element_text(size = 18, angle = 45, hjust =1, vjust = -1),
legend.text = element_text(size = 20),
strip.text = element_text(size =22),
legend.key.width = unit(5,"line"))


ggplot() + 
geom_line(data = nuc[nuc$CHROM != "MT",], aes(x = BIN_START, y = PI, colour = FileFrom), alpha = 0.8) +
facet_wrap(~CHROM, scale = "free_x", nrow = 3) +
scale_fill_manual(values = mycols) +
xlab("Position") +
ylab("Nucleotide Diversity (π)") +
theme_classic() +
guides(color = guide_legend(override.aes = list(size=5))) +
theme(axis.text = element_text(size = 22),
strip.background = element_rect(fill = "#a8bac3"),
axis.title = element_text(size = 22),
legend.position = "bottom",
legend.text = element_text(size = 22),
strip.text = element_text(size =22),
legend.key.width = unit(5,"line"))

all_plot
library(patchwork)
both <- pr_plot + all_plot
both
### BOXPLOTS ###
mycols <-  unikn::usecol(c("#ffaa5c", "#db7376", "#52a884", "#bc80bd", "#465c7a"), n = 8)
#mycols <- unikn::usecol(c("#CDD3D5", "#75B8C8", "#197278", "#03045E", "#922D50"), n = 8)
#mycols <- c("#c7522a", "#e5c185", "#f0daa5", "#fbf2c4", "#b8cdab", "#74a892", "#008585", "#004343")
    plot1 <- ggplot() +
        geom_boxplot(data = nuc[nuc$CHROM != "MT",], aes(x = FileFrom, y = PI, fill = FileFrom), outlier.shape = NA) +
        scale_fill_manual(values = mycols) +
        scale_y_continuous(limits = c(0, 0.0015)) +
        xlab("") +
        ylab("Nucleotide DIversity (π)") +
        ylim(c(0,0.00065))+
        ggtitle("Nucleotide diversity across each chromsome \n (window size 100kbp)") +
        guides(fill = guide_legend(title = "Country")) +
        theme_classic() +
        theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = "None",
        axis.text = element_text(size = 18, angle = 45, hjust =1, vjust = -1),
        axis.title = element_text(size = 18),
        strip.background = element_rect(fill = "#a8bac3"),
        title = element_text(size = 20)) +
        facet_wrap(~CHROM, ncol = 4)

   plot2 <-  ggplot() +
        geom_boxplot(data = nuc[nuc$CHROM != "MT",], aes(x = FileFrom, y = PI, fill = FileFrom), outlier.shape = NA) +
        scale_fill_manual(values = mycols) +
        scale_y_continuous(limits = c(0, 0.0015)) +
        xlab("") +
        #ylab("Nucleotide DIversity (π)") +
        ylim(c(0,0.00065))+
        ggtitle("Nucleotide diversity for each population \n(window size 100kb)") +
        guides(fill = guide_legend(title = "Country")) +
        theme_classic() +
        theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = "None",
        axis.text = element_text(size = 18, angle = 45, hjust =1),
        axis.title = element_text(size = 18),
        title = element_text(size = 20)) 
        #facet_wrap(~CHROM, ncol = 4)

library(cowplot)
library(patchwork)
nuc_div_plot <- plot2 + plot1 + plot_layout(nrow = 2)
nuc_div_plot

## scatter grid per chromosome and country
ggplot() + 
        geom_line(data = nuc[nuc$CHROM != "MT",], aes(x = BIN_START, y = PI, alpha = 0.4)) +
        facet_grid(vars(FileFrom), vars(CHROM), scales = "free_x") +
        #geom_hline(yintercept = mean(sub$PI), colour = "#a8bac3", size = 1.5) +
        xlab("") +
        ylab("Nucleotide Diversity (π)") +
        #ggtitle(paste0("Nuclotide diversity across each chromsomes for ", sub$FileFrom[1])) +
        #ggtitle(paste0(nuc$FileFrom[1])) +
        theme_classic() +
        theme(strip.text = element_text(size = 16),
        axis.text = element_text(size = 14, angle = 45, hjust =1),
        axis.title = element_text(size = 18),
        strip.background = element_rect(fill = "#a8bac3"),
        legend.position = "none",
        title = element_text(size = 18)) 
    

## plotting by country

plot_list <- list()
for (i in unique(nuc$FileFrom)){

    sub = nuc[nuc$FileFrom == i,]

    plot <- ggplot() + 
        geom_point(data = sub[sub$CHROM != "MT",], aes(x = BIN_START, y = PI, alpha = 0.4)) +
        facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
        geom_hline(yintercept = mean(sub$PI), colour = "#a8bac3", size = 1.5) +
        xlab("") +
        ylab("Pi Value") +
        #ggtitle(paste0("Nuclotide diversity across each chromsomes for ", sub$FileFrom[1])) +
        ggtitle(paste0(sub$FileFrom[1])) +
        theme_minimal() +
        theme(strip.text = element_text(size = 16),
        axis.text = element_text(size = 16, angle = 45, hjust =1),
        axis.title = element_text(size = 18),
        legend.position = "none",
        title = element_text(size = 18)) 
    
    ggplot() + 
        geom_line(data = sub[sub$CHROM != "MT",], aes(x = BIN_START, y = PI)) +
        facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
        geom_hline(yintercept = mean(sub$PI), colour = "#a8bac3", size =1.5) +
        xlab("Position") +
        ylab("Pi Value") +
        ggtitle(paste0("Nuclotide diversity across each chromsomes for ", sub$FileFrom[1])) +
        theme_minimal() +
        theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 16, angle = 45, hjust =1),
        axis.title = element_text(size = 16),
        title = element_text(size = 18)) 

    print(i)
    print(plot)
    plot_list[[i]] <- plot
    #ggsave(paste0(i, "_nuc_div_plot.jpeg"), plot)
    }

all_plots <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]] + plot_list[[6]] + plot_list[[7]] + plot_list[[8]]
all_plots


## to plot individually
ggplot() + 
    geom_point(data = nuc, aes(x = BIN_START, y = PI)) +
    facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
    geom_hline(yintercept = 0.0006, colour = "firebrick") +
    xlab("Position") +
    ylab("Pi Value") +
    ggtitle("Nuclotide diversity across each chromsomes") +
    theme_classic() +
    theme(strip.text = element_text(size = 20),
    axis.text = element_text(size = 16, angle = 45, hjust =1),
    axis.title = element_text(size = 18),
    title = element_text(size = 20))

mycols <- c("#c7522a", "#e5c185", "#f0daa5", "#fbf2c4", "#b8cdab", "#74a892", "#008585", "#004343")
ggplot() +
    geom_boxplot(data = nuc, aes(x = FileFrom, y = log(PI), fill = FileFrom)) +
    scale_fill_manual(values = mycols) +
    xlab("Country") +
    ylab("Nucleotide Diversity (π)") +
    ggtitle("π across each chromsomes") +
    guides(fill = guide_legend(title = "Country")) +
    theme_classic() +
    theme(strip.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.position = "None",
    axis.text = element_text(size = 18, angle = 45, hjust =1),
    axis.title = element_text(size = 18),
    title = element_text(size = 20))



## make into a table

table = data.frame(country = unique(nuc$FileFrom), mean = NA, row.names = NULL, check.rows = FALSE, check.names = TRUE)

for (i in unique(nuc$FileFrom)){

    sub = nuc[nuc$FileFrom == i,]
    sub_mean <- mean(sub$PI)

    table$mean[table$country == i] <-sub_mean
}

head(table)
write.csv(table, "Nucleotide_diversity_summary.csv")