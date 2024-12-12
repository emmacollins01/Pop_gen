## LD ##
library(tidyr)
library(gaston)
library(tidyverse)
library(dplyr)

setwd("/mnt/storage12/emma/PR_LD/IR_gene/")
setwd("/mnt/storage12/emma/PR_LD/of_interest/")

#a<- read.delim("American_Samoa_LD.hap.ld", stringsAsFactors = FALSE)
#b<- read.delim("Brazil_LD.hap.ld", stringsAsFactors = FALSE)

files <- list.files(pattern = "*hap.ld")
countries <- c()
# loop through the files
for (i in 1:length(files)){
  # save the clean filename as a char var because it will be called again
  fnClean = substr(files[i], start = 1, stop = nchar(files[i])-10)
  
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
    addFnCMD = paste0(fnClean,'$FileFrom = files[i]')
    print(addFnCMD) # check the cmd
    eval(parse(text=addFnCMD)) }
}

countries

# no record in Australia
all <- rbind(Brazil, Burkina_Faso, Kenya,
            Mexico, Puerto_Rico,
            Thailand, Uganda, USA)

rm(American_Samoa, Brazil, Burkina_Faso, Costa_Rica, East_Kenya, East_Kenya_Outlier, 
   French_Polynesia, Gabon, Ghana, LAB, Madagascar, Mauritius, Mexico, Nigeria,
   S._Africa, Senegal, Thailand, Uganda, USA, West_Kenya)

ld <- all


ld <- subset(ld, R.2 != "NaN")
ld <- ld[!is.na(ld$R.2),]
#decide if needs to be subset
#summary(ld$R.2)
#ld <- subset(ld, R.2 > 0.0001)

chrom2 <- subset(ld, CHR == "NC_035108.1")
chrom3 <- subset(ld, CHR == "NC_035109.1")
ld <- chrom2
ld <- chrom3


ld_og <- ld

for (i in 1:length(unique(ld_og$FileFrom))){
  
  country_i = unique(ld_og$FileFrom)[i]
  ld_i <- subset(ld_og, FileFrom == country_i)
  ld_country <- ld_i
  print(country_i)
  
  a <- c(unique(ld_i$POS1), (unique(ld_i$POS2)))
  no_snps <- length(unique(a))
  
  #for plot name
  file = unique(substr(ld_i$FileFrom, start = 1, stop = nchar(ld_i$FileFrom)-10))
  region = paste0(unique(ld_i$CHR)[1], "_", unique(ld_i$CHR)[2])
  
  #subset to only important columns
  ld_i <- ld_i[,c(2,3,5)]
  colnames(ld_i) <- c("POS1", "POS2", "R.2")
  ld_i <- ld_i[!duplicated(ld_i),]
  ld_i <- spread(ld_i, POS2, R.2)
  
  ld_i_snps <- unique(ld_i$POS1)
  ld_i_snps2 <- unique(colnames(ld_i))
  ld_i_snps2 <- ld_i_snps2[2:length(ld_i_snps2)]
  all_snps <- unique(c(ld_i_snps, ld_i_snps2))
  
  df <- data.frame(matrix(NA, nrow = no_snps, ncol = (no_snps+1)))
  df$X1 <- all_snps
  all_snps2 <- paste0("X", all_snps)
  all_snps3 <- c("POS", all_snps2)
  colnames(df) <- all_snps3
  
  # Need to set this so its per country
  ld <- ld_country
  
  for (i in unique(all_snps)){
    for (j in unique(all_snps)){
      
      i
      i_x <- paste0("X", i)
      j 
      if(length(ld$CHR[ld$POS1 == i & ld$POS2 == j])>0){
        df[df$POS == j, i_x] <- ld$R.2[ld$POS1 == i & ld$POS2 == j]
      }
      else if(length(ld$CHR[ld$POS1 == j & ld$POS2 == i])>0){
        df[df$POS == j, i_x] <- ld$R.2[ld$POS1 == j & ld$POS2 == i]
      }
      
    }
    
  }
  df[is.na(df)] <- 0
  ld_i <- df
  ld_i <- as.matrix(ld_i)
  snps <- ld_i[,-1]
  #snps <- print(snps, quote=FALSE, na.print="NA")
  rownames(snps) <- ld_i[,1]
  snp_pos <- as.numeric(unique(ld_i[,1]))
  snps <- matrix(as.numeric(snps), ncol = ncol(snps))
  colnames(snps) <- snp_pos
  rownames(snps) <- snp_pos
  LD.plot(snps, snp_pos, max.dist = Inf, depth = nrow(snps), pdf.file = paste0(file, "_", region, '.pdf'))
  
}


####### Table of highest values
summary(all$R.2)
# mean 0.03
summary(all$POS1)

regions <- read.csv("../../../OneDrive - London School of Hygiene and Tropical Medicine/NCBI_Aedes_analysis/Region_metadata.csv")

all$FileFrom <- substr(all$FileFrom, start = 1, stop = nchar(all$FileFrom)-10)

all$region1 <- regions$Region1[match(all$FileFrom, regions$Country)]
all$region2 <- regions$Region2[match(all$FileFrom, regions$Country)]
all$Continent <- regions$Region3[match(all$FileFrom, regions$Country)]

library(ggpubr)
plot_data <- all
plot_data <- subset(plot_data, FileFrom != "LAB")
plot_data <- plot_data[!is.na(plot_data$Continent),]
#RDL
rdl <- ggplot() + 
  geom_point(data = plot_data[plot_data$POS1 >=41630421 & plot_data$POS1 <= 41859336 & plot_data$Continent != "NA",], aes(x = POS1, y = R.2, colour = Continent))+
  theme_bw() +
  xlab("Position") +
  ggtitle("GABA")
# VGSC
vgsc <- ggplot() + 
  geom_point(data = plot_data[plot_data$POS1 >=315931877 & plot_data$POS1 <= 316404055,], aes(x = POS1, y = R.2, colour = Continent))+
  theme_bw() + 
  xlab("Position") +
  ggtitle("VGSC")
#ACE
ace <- ggplot() + 
  geom_point(data = plot_data[plot_data$POS1 >=161486025 & plot_data$POS1 <= 161487543,], aes(x = POS1, y = R.2, colour = Continent))+
  theme_bw() +
  xlab("Position") +
  ggtitle("ACE")
# GSTE2
gste2 <- ggplot() + 
  geom_point(data = plot_data[plot_data$POS1 >=351633367 & plot_data$POS1 <= 351633799,], aes(x = POS1, y = R.2, colour = Continent))+
  theme_bw() +  
  xlab("Position") +
  ggtitle("GSTE2") 


ggarrange(rdl, vgsc, ace, gste2)

all_high <- subset(all, R.2 > 0.5)

summary(all_high$R2)
table(all_high$POS1)
table(all_high$POS2)

missense <- read.csv("../../../OneDrive - London School of Hygiene and Tropical Medicine/NCBI_Aedes_analysis/missense_pos_v2.csv")

all_high$csq1 <- missense$Reference.organism.position[match(all_high$POS1, missense$POS)]
all_high$csq2 <- missense$Reference.organism.position[match(all_high$POS2, missense$POS)]

all_high$FileFrom <- substr(all_high$FileFrom, start = 1, stop = nchar(all_high$FileFrom)-10)

all_high$Gene1 <- missense$GENE[match(all_high$POS1, missense$POS)]
all_high$Gene2 <- missense$GENE[match(all_high$POS2, missense$POS)]

all_high$CHR <- missense$CHROM[match(all_high$POS1, missense$POS)]
all_high$R.2 <- round(all_high$R.2, digits = 2)
all_high$D<- round(all_high$D, digits = 2)
all_high$Dprime <- round(all_high$Dprime, digits = 2)

colnames(all_high) <- c("Chr", "Pos1", "Pos2", "n_chr", "R2", "D", "Dprime", "Country", "Csq1", "Csq2", "Gene1", "Gene2")

write.csv(all_high, "Ld_over_0.5_table.csv", index = FALSE)

length(unique(paste(all_high$Pos1, all_high$Pos2)))
all_high$uni <- paste(all_high$Pos1, all_high$Pos2)
test <- all_high %>% dplyr::group_by(uni) %>% dplyr::summarise(n= n())


table(all_high$Gene1)
table(all_high$Gene2)
table(all_high$Country)


