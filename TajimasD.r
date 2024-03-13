
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

library(ggplot2)

## Subset and change names of chromosomes if needed!
nuc <- subset(all, CHROM %in% c(35107, 35108, 35109, 35159))
nuc$CHROM[nuc$CHROM == 35107] <- "1"
nuc$CHROM[nuc$CHROM == 35108] <- "2"
nuc$CHROM[nuc$CHROM == 35109] <- "3"
nuc$CHROM[nuc$CHROM == 35159] <- "MT"

colnames(nuc)
options(scipen = 999)

for (i in unique(nuc$FileFrom)){

    sub = nuc[nuc$FileFrom == i,]

    plot <- ggplot() + 
        geom_point(data = sub, aes(x = BIN_START, y = TajimaD)) +
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

## to plot individually
ggplot() + 
    geom_point(data = nuc, aes(x = BIN_START, y = TajimaD)) +
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

## make into a table

table = data.frame(country = unique(nuc$FileFrom), mean = NA, row.names = NULL, check.rows = FALSE, check.names = TRUE)

for (i in unique(nuc$FileFrom)){

    sub = nuc[nuc$FileFrom == i,]
    sub <- sub[sub$N_SNPS != 0,]
    sub_mean <- mean(sub$TajimaD, rm.na = TRUE)

    table$mean[table$country == i] <-sub_mean
}

head(table)
write.csv(table, "Tajima_summary.csv")
