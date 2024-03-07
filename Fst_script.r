
library(data.table)
library(qqman)

## https://lshtm-1.gitbook.io/sequence-analysis/r-scripts/fst

#nonbinary SNP matrix
snp<-read.delim("/mnt/storage12/emma/PR_combine/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.noniupac.mat",sep="\t",header=T)
snp<-read.delim("/mnt/storage12/emma/PR_combine/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.noniupac.mat",sep="\t",header=T)
snp<-as.data.frame(snp)

#Binary SNP matrix
snp_01<-read.delim("/mnt/storage12/emma/PR_combine/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.mat.bin",sep="\t",header=T)
snp_01<-read.delim("/mnt/storage12/emma/PR_combine/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.mat.bin",sep="\t",header=T)
snp_01<-as.data.frame(snp_01)

#Metadata
metadata<-fread("/mnt/storage12/emma/PR_combine/all_sample_country_metadata.csv",sep=",",header=T)
metadata2<-fread("/mnt/storage12/emma/PR_combine/all_samples_PR_resistance.csv",sep=",",header=T)
metadata<-as.data.frame(metadata)
metadata2<-as.data.frame(metadata2)
head(metadata)

snp_calls<-snp[,-c(1:3)]
snp_calls_01<-snp_01[,-c(1:3)]

legend<-snp[,c(1:3)]

fsts <- function(haps,lege,pop1,pop2){
  
  ## remove monomorphics ##
  haps1 <- haps[,c(pop1,pop2)]
  lege1 <- data.frame(pos=lege)

  ## two populations  
  p1 <- length(pop1)
  p2 <- length(pop2)

  haps1 <- data.frame(lapply(haps1[,1:ncol(haps1)], as.numeric))
  #haps1 <- data.frame(lapply(haps1[,(p1+1):(p2+p1)], as.numeric))

  m1 <- apply(haps1[,1:p1],1,mean,na.rm=T)
  m2 <- apply(haps1[,(p1+1):(p1+p2)],1,mean,na.rm=T)
  pbar <- (p1*m1 + p2*m2) / (p1 + p2)
  nbar <- (p1 + p2) / 2
  
  fst <-  ( p1*(m1-pbar)^2 + p2*(m2-pbar)^2 )/(p1+p2)/pbar/(1-pbar) ## should be n-1
  fst[is.na(fst)] <- 0
  pp1 <- signif(m1,3)
  pp2 <- signif(m2,3)
  tt <- cbind(lege1,pp1,pp2,fst)
  colnames(tt) <-c(colnames(lege1),"Ppop1","Ppop2","fst")
  tt
}


# Create a list for unique chromosome IDs
list_chr<-unique(snp[,1])

#A population could be a group of countries
pop_1<-unique(sort(c(which(metadata_rem$Country=="Cameroon"), which(metadata_rem$Country=="Gabon"), which(metadata_rem$Country=="Kenya"), which(metadata_rem$Country=="Nigeria"),
                     which(metadata_rem$Country=="South_Sudan"),which(metadata_rem$Country=="Uganda"), which(metadata_rem$Country=="Liberia"), which(metadata_rem$Country=="Malawi"))))

#Or just one country
pop_1<-which(metadata$Country=="Uganda") ## Change this to Pop2 of your populations
pop_2<-which(metadata$Country=="Thailand") ## Change this to Pop2 of your populations

pop_1<-which(metadata2$Resistance=="Most_Resistant") ## Change this to Pop2 of your populations
pop_2<-which(metadata2$Resistance=="Least_Resistant") ## Change this to Pop2 of your populations


## use function to create table of SNPs
fst_mut_pop <- fsts(snp_calls_01,legend,pop_1,pop_2)

## get tt by running sections within function - may work now making value numeric
#fst_mut_pop <- tt

Fst_values<-fst_mut_pop$fst
snp_chr_factor<-legend[,1]

for(i in 1:length(list_chr)){
  snp_ite<-list_chr[i]
  snp_chr_factor[snp_chr_factor==list_chr[i]]<-i
}


# Removing letters and symbols so can make numeric
fst_mut_pop$pos.chr <- gsub("^[A-Za-z]+_", "", fst_mut_pop$pos.chr)
unique(fst_mut_pop$pos.chr)

## remove contig to make more simple
fst_mut_pop1<- fst_mut_pop[(fst_mut_pop$pos.chr == "035107.1" |
  fst_mut_pop$pos.chr == "035108.1" |
  fst_mut_pop$pos.chr == "035109.1" |
  fst_mut_pop$pos.chr == "035159.1"),]

## remove contig to make more simple
fst_mut_pop1<- fst_mut_pop[(fst_mut_pop$pos.chr == "35107" |
  fst_mut_pop$pos.chr == "35108" |
  fst_mut_pop$pos.chr == "35109" |
  fst_mut_pop$pos.chr == "35159"),]

fst_mut_pop1$pos.chr <- as.numeric(fst_mut_pop1$pos.chr)

as.data.frame(table(fst_mut_pop1$pos.chr))
fst_mut_pop1$snp <- paste0("snp_", 1:nrow(fst_mut_pop1))

fst_mut_pop1 <- na.omit(fst_mut_pop1)

snpofinterest <- str(paste0("snp_", 150000:151000))
#manhattan(Pop_qqman_fst,logp=FALSE,main="Fst Pop1 vs. Pop2",ylab="Fst",las=1,cex.axis=1.4,cex.lab=2,suggestiveline=0.5,genomewideline=FALSE)
manhattan(fst_mut_pop1[1:100000,], highlight = snpofinterest, 
  chr = "pos.chr", bp = "pos.pos", p = "fst", snp = "snp", 
  logp=FALSE, main="Fst Pop1 vs. Pop2", ylab="Fst",
  las=1, cex.axis=1.4, cex.lab=2, suggestiveline=0.8, genomewideline=FALSE)

manhattan(fst_mut_pop1,
  chr = "pos.chr", bp = "pos.pos", p = "fst", snp = "snp",
  logp=FALSE, main="Fst Most Resistant vs. Least Resistant", ylab="Fst", ylim = c(0,1.1),
  las=1, cex.axis=1.4, cex.lab=2, suggestiveline=0.8, genomewideline=FALSE)

unique(fst_mut_pop$fst)
colnames(fst_mut_pop)

library(ggplot2)
options(scipen =999)
plot <- ggplot()+
geom_point(data = fst_mut_pop1, aes(x = pos.pos, y = fst, alpha = 0.5, fill = fst_mut_pop1$pos.chr), shape = 2)+
facet_wrap(~pos.chr, ncol = 4) +
theme_classic() 
plot + geom_hline(yintercept = 0.8, linetype = "dashed", colour = "firebrick")

ggplot() +
geom_bar(data = fst_mut_pop1, aes(x = pos.pos,))

# if you want to save
tiff("/mnt/storage12/emma/PR_fst/PCA_most_vs_least_resistant.tiff",height=1200,width=1800)
manhattan(fst_mut_pop1, chr = "pos.chr",bp = "pos.pos", 
p = "fst", snp = "pos.ref", logp=FALSE,
main="Fst Most Resistant vs Least Resistant",
ylab="Fst",las=1,cex.axis=1.4,cex.lab=2,
ylim = c(0,1.1),
suggestiveline=0.8,genomewideline=FALSE)
#manhattan(Pop_qqman_fst,logp=FALSE,main="Fst Pop1 vs. Pop2",ylab="Fst",las=1,cex.axis=1.4,cex.lab=2,suggestiveline=0.5,genomewideline=FALSE)
dev.off()

