
## Start by making file of samples names per population
# then file with all of these files eg/ ls SRR.txt > countries.txt
# then you can run separate_vcfs.sh to split multisample vcf into separate populations
# OR 
#for i in $(ls *SRR.txt); do bcftools view -S $i all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.vcf.gz -Oz -o ${i}.vcf.gz; done

## to make nucleotide diversity metrics
## Run vcftools - run on 10kb windows
## These files are in input here
#for i in $(ls *SRR.vcf.gz);do vcftools --gzvcf $i --window-pi 10000 --out ${i}_nuc_div_window_10kb;done

setwd("/mnt/storage12/emma/nuc_div/")
rm(list = ls())
#nuc <- read.csv("all_nuc_div_window_10kb.windowed.pi", sep = "\t", header = TRUE)

files <- list.files(pattern = "*tajima_window_10kb.Tajima.D")
files <- list.files(pattern = "*tajima_window_20kb.Tajima.D")
files <- list.files(pattern = "*tajima_window_100kb.Tajima.D")
files <- list.files(pattern = "*tajima_window_300kb.Tajima.D")
files <- list.files(pattern = "*tajima_window_500kb.Tajima.D")
files <- list.files(pattern = "*tajima_window_5000kb.Tajima.D")

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

## combine all dataframes together - now they have country column
dataframes <- Filter(function(x) is(x, "data.frame"), mget(ls()))
all <- do.call(rbind, dataframes)
head(all)

## data set to write to csv
summary(all$TajimaD)
sig <- all[abs(all$TajimaD) >= 2,]
sig <- subset(all, TajimaD > 2)
sig2 <- subset(all, TajimaD < -2)
sig <- rbind(sig, sig2)
sig <- sig[complete.cases(sig),]
nrow(sig)
head(sig)
summary(sig$TajimaD)
summary(sig$N_SNPS)
write.csv(sig, "all_sig_2_tajimasD_100kb.csv", row.names = FALSE)
sig <- read.csv("/mnt/storage12/emma/nuc_div/all_sig_2_tajimasD_100kb.csv")
head(sig)
summary(sig$TajimaD)
summary(all$TajimaD)
table(sig$FileFrom)
length(unique(paste(sig$CHROM, sig$BIN_START)))
library(ggplot2)

test <- sig %>% group_by(CHROM, BIN_START) %>% summarise(n= n())
test2 <- sig %>% group_by(FileFrom) %>% summarise(count = n_distinct(unique(paste(CHROM, BIN_START))))

## Subset and change names of chromosomes if needed!
nuc <- subset(all, CHROM %in% c(35107, 35108, 35109, 35159))
nuc$CHROM[nuc$CHROM == 35107] <- "1"
nuc$CHROM[nuc$CHROM == 35108] <- "2"
nuc$CHROM[nuc$CHROM == 35109] <- "3"
nuc$CHROM[nuc$CHROM == 35159] <- "MT"

# overall means
summary(nuc$TajimaD)
library(dplyr)
nuc_median <- nuc[!is.na(nuc$TajimaD),] %>% group_by(FileFrom) %>% summarise(median = median(TajimaD))
nuc_median <- nuc[!is.na(nuc$TajimaD),] %>% group_by(FileFrom, CHROM) %>% summarise(median = median(TajimaD))


summary(nuc$TajimaD[nuc$FileFrom == "Puerto_Rico"])
PR_taj <- nuc[nuc$FileFrom == "Puerto_Rico" & !is.na(nuc$TajimaD),]
PR_taj_median <- PR_taj %>% group_by(CHROM) %>% summarise(median = median(TajimaD))

colnames(nuc)
options(scipen = 0)
mycols <-  unikn::usecol(c("#ffaa5c", "#db7376", "#52a884", "#bc80bd", "#465c7a"), n = 8)

library(ggplot2)
all_plot <- ggplot() + 
geom_line(data = nuc[nuc$CHROM != "MT",], aes(x = BIN_START, y = TajimaD, colour = FileFrom), alpha = 1, size = 1.22) +
facet_wrap(~CHROM, scale = "free_x", nrow = 3) +
scale_colour_manual(values = mycols) +
xlab("Position") +
ylab("Tajima's D") +
theme_classic() +
guides(color = guide_legend(override.aes = list(size=5))) +
theme(axis.text = element_text(size = 22),
legend.title = element_blank(),
axis.title = element_text(size = 22),
legend.position = "bottom",
legend.text = element_text(size = 22),
strip.text = element_text(size =22),
legend.key.width = unit(5,"line"))

all_plot
ggsave("Tajima_5MB_all.png", all_plot)



###################################
for (i in unique(nuc$FileFrom)){

    sub = nuc[nuc$FileFrom == i,]

    plot <- ggplot() + 
        geom_point(data = sub, aes(x = BIN_START, y = TajimaD), alpha = 0.5) +
        facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
        #geom_hline(yintercept = mean(sub$PI), colour = "firebrick") +
        xlab("Position") +
        ylab("Tajimas D") +
        ggtitle(paste0("Tajima's D across each chromsomes for ", sub$FileFrom[1])) +
        theme_minimal() +
        theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 16, angle = 45, hjust =1),
        axis.title = element_text(size = 18),
        title = element_text(size = 20)) 
    print(i)
    print(plot)
    ggsave(paste0(i, "_tajima_plot.jpeg"), plot)

    }

## scatter grid per chromosome and country
ggplot() + 
        geom_line(data = nuc[nuc$CHROM != "MT",], aes(x = BIN_START, y = TajimaD, alpha = 0.4)) +
        facet_grid(vars(FileFrom), vars(CHROM), scales = "free_x") +
        #geom_hline(yintercept = 2.5, colour = "#a8bac3", size = 1.5, alpha = 0.5) +
        #geom_hline(yintercept = -2.5, colour = "#a8bac3", size = 1.5, alpha = 0.5) +
        xlab("") +
        ylab("Tajima's D Value") +
        #ggtitle(paste0("Nuclotide diversity across each chromsomes for ", sub$FileFrom[1])) +
        #ggtitle(paste0(nuc$FileFrom[1])) +
        theme_classic() +
        theme(strip.text = element_text(size = 16),
        axis.text = element_text(size = 14, angle = 45, hjust =1),
        axis.title = element_text(size = 18),
        strip.background = element_rect(fill = "#a8bac3"),
        legend.position = "none",
        title = element_text(size = 18))


## boxplot to compare all countries
#mycols <- c("#db7376", "#52a884", "#bc80bd", "#465c7a", "#ffaa5c", "#480355", "#97B2D8", "#D1CFE2")
#mycols <- c("#c7522a", "#e5c185", "#f0daa5", "#fbf2c4", "#b8cdab", "#74a892", "#008585", "#004343")
mycols <-  unikn::usecol(c("#ffaa5c", "#db7376", "#52a884", "#bc80bd", "#465c7a"), n = 8)

ggplot() +
    geom_boxplot(data = nuc, aes(x = FileFrom, y = TajimaD, fill = FileFrom)) +
    scale_fill_manual(values = mycols) +
    xlab("") +
    ylab("Tajimas D") +
    geom_hline(yintercept = 2, linetype = "dashed") +
    geom_hline(yintercept = -2, linetype = "dashed") +
    ggtitle("Tajima's D across each chromsomes") +
    guides(fill = guide_legend(title = "Country")) +
    theme_classic() +
    theme(strip.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.position = "None",
    axis.text = element_text(size = 18, angle = 45, hjust =1),
    axis.title = element_text(size = 24),
    title = element_text(size = 20))

## to plot individually
ggplot() + 
    geom_point(data = nuc[nuc$FileFrom == "Puerto_Rico",], aes(x = BIN_START, y = TajimaD)) +
    facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
    #geom_hline(yintercept = 0.0006, colour = "firebrick") +
    xlab("Position") +
    ylab("Tajimas D") +
    ggtitle("Tajima's D across each chromsomes") +
    theme_classic() +
    theme(strip.text = element_text(size = 20),
    axis.text = element_text(size = 16, angle = 45, hjust =1),
    axis.title = element_text(size = 18),
    title = element_text(size = 20))

ggplot() +
    geom_boxplot(data = nuc[nuc$FileFrom == "Puerto_Rico" &nuc$CHROM != "MT",], aes(x = FileFrom, y = TajimaD, fill = FileFrom)) +
    scale_fill_manual(values = "#709C94") +
    facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
    ylab("Tajimas D") +
    xlab("") +
    geom_hline(yintercept = 2, linetype = "dashed") +
    geom_hline(yintercept = -2, linetype = "dashed") +
    ggtitle("Tajima's D per chromsomes for Puerto Rico") +
    guides(fill = guide_legend(title = "Country")) +
    theme_classic() +
    theme(strip.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.position = "None",
    axis.text.x=element_blank(),
    axis.text = element_text(size = 18, angle = 45, hjust =1),
    axis.title.y = element_text(size = 24),
    title = element_text(size = 20))

## make into a table

table = data.frame(country = unique(nuc$FileFrom), mean = NA, row.names = NULL, check.rows = FALSE, check.names = TRUE)

for (i in unique(nuc$FileFrom)){

    sub = nuc[nuc$FileFrom == i,]
    sub <- sub[sub$N_SNPS != 0,]
    sub_mean <- mean(sub$TajimaD, rm.na = TRUE)

    table$mean[table$country == i] <-sub_mean
}

head(table)
#write.csv(table, "Tajima_summary_20kbp.csv")

nuc2 <- nuc[!(is.nan(nuc$TajimaD)),]
high_taj <- nuc2 %>% top_n(TajimaD, n = 50)
low_taj <- nuc2 %>% slice_min(TajimaD, n = 50)

top_pos <- nuc2%>% filter(TajimaD > 2.5) %>% group_by(BIN_START) %>% summarise(n = n())
low_pos <- nuc2%>% filter(TajimaD < -2) %>% group_by(BIN_START) %>% summarise(n = n())

high_taj_2.5 <- nuc %>% filter(TajimaD >2.5 | TajimaD < -2.5)

length(unique(top_pos$BIN_START))
length(unique(low_pos$BIN_START))

table(high_taj_2.5$FileFrom)
table(high_taj_2.5$CHROM)
table(high_taj_2.5$CHROM, high_taj_2.5$FileFrom)

## most common
nuc_group <- nuc %>% filter(TajimaD > 2.5) %>% group_by(CHROM, BIN_START) %>% summarise(n= n()) %>% top_n(n, n = 20)
nuc_group2 <- nuc %>% filter(TajimaD < -2.5) %>% group_by(CHROM, BIN_START) %>% summarise(n= n()) %>% top_n(n, n = 20)

### PR
chroms <- c(35107, 35018, 35109, 35159)
PR_sub <- subset(Puerto_Rico, CHROM %in% chroms)
PR_high <- subset(PR_sub, abs(TajimaD) > 2.5)
