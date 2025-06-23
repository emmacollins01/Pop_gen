

# Variants in genes of interest that PASS only 

all <- read.csv("/mnt/storage12/emma/delly/June24/merged.sample_filt.BNDfiltered.pass.txt", header = FALSE, sep = "\t")

head(all)
colnames(all) <- c("SAMPLE", "CHROM","POS", "REF", "ALT", "?", "GT", "CHANGE", "ANN")

all$uni <- paste(all$CHROM, all$POS)

length(unique(all$uni))
length(unique(all$POS))

for (i in 1:nrow(all)){
    print(i)
    all$gene[i] <- paste(unique(unlist(str_extract_all(all$ANN[i], "LOC\\d+"))), collapse = ", ")


}

nrow(all)
length(unique(all$gene))

library(stringr)
library(dplyr)
del <- read.csv("/mnt/storage12/emma/delly/June24/merged.sample_filt.BNDfiltered.pass.of_interest_expanded.txt", header = FALSE, sep = "\t")

head(del)
colnames(del) <- c("SAMPLE", "CHROM","POS", "REF", "ALT", "?", "GT", "CHANGE", "ANN")

del$uni <- paste(del$CHROM, del$POS)

length(unique(del$uni))
length(unique(del$POS))

for (i in 1:nrow(del)){
    del$gene[i] <- unique(unlist(str_extract_all(del$ANN[i], "LOC\\d+")))


}

nrow(del)
length(unique(del$gene))


high <- del[grep("HIGH", del$ANN),]
mis <- del[grep("missense", del$ANN),]

length(unique(high$POS))
summary(as.factor(high$CHROM))

high2 <- del %>% group_by(CHROM) %>% summarise(n= n_distinct(POS)) 

high_sum <- del %>% filter(GT != "0/0") %>% group_by(CHROM, POS, REF, ALT) %>% summarise(n = n_distinct(SAMPLE))
top <- high_sum %>% top_n(n, n = 10)
high_sum <- high %>% filter(GT != "0/0") %>% group_by(CHROM, POS, REF, ALT) %>% summarise(n = n_distinct(SAMPLE))

mode <- high_sum[high_sum$n == 33,]

test <- subset(del, POS %in% mode$POS)

unique(high_sum$POS[high_sum$n == 33])
high_sum2 <- high_sum %>% filter(n ==33) %>% group_by(CHROM) %>% summarise(n= n_distinct(POS)) 

for (i in 1:nrow(high)){
    high$type[i] <- unlist(strsplit(high$ANN[i], "\\|"))[[2]]
}
unique(high$type)
unique(high$ALT)

common <- high_sum$POS[high_sum$n == 33]
mode <- high[high$POS %in% common,]
for (i in 1:nrow(mode)){
    mode$impact[i] <- unlist(strsplit(mode$ANN[i], "\\|"))[[1]]
    mode$effect[i] <- unlist(strsplit(mode$ANN[i], "\\|"))[[2]]
    mode$gene[i] <- unlist(strsplit(mode$ANN[i], "\\|"))[[5]]

}
unique(mode$impact)
unique(mode$effect)
unique(mode$gene)

mode2 <- mode %>% select(CHROM, POS, ALT, REF, impact, effect, gene) %>% group_by(CHROM, POS, ALT, REF, impact, effect, gene) %>% summarise(n =n())

mode_high <- mode[grep("HIGH", mode$ANN),]
unique(mode_high$POS)


mode_dup <- mode[grep("DUP", mode$ANN),]
unique(mode_dup$POS)


### positions in most common genes
test <- del %>% group_by(POS, gene) %>% summarise(no_pos = n())
test2 <- test %>% group_by(gene) %>% summarise(no_svs =n())
genes <- c("LOC5572215", "LOC23687794", "LOC5565389", "LOC5577718", "LOC23687733", "LOC5571510", "LOC23687755", "LOC5573499", "LOC5571895", "LOC5572902")
subset <- del[del$gene %in% genes,]

write.csv(subset, "/mnt/storage12/emma/delly/June24/top_10_genes.csv", row.names = FALSE)

### INV ###

inv <- vcf[vcf$type == "SVTYPE=INV",]
unique(inv$V2[inv$V1 == "NC_035107.1"])


### VGSC

vgsc <- high[grep("LOC5567355", high$ANN),]
vgsc_sum <- vgsc %>% filter(GT != "0/0")



dup <- del[grep("DUP", del$ANN),]
dup_high <- dup[grep("HIGH", dup$ANN),]
dup_high_filt <- dup_high %>% filter(GT != "0/0" & GT != "./.")
unique(dup_high_filt$POS)

for (i in 1:nrow(dup_high_filt)){
    dup_high_filt$type[i] <- unlist(strsplit(dup_high_filt$ANN[i], "\\|"))[[2]]
}

unique(dup_high_filt$type)


vcf <- read.csv("/mnt/storage12/emma/delly/merged.sample_filt_of_interest_ann_pass_noheader.vcf.gz",header = FALSE, sep = "\t")

for (i in 1:nrow(vcf)){
    vcf$type[i] <- unlist(strsplit(vcf$V8[i], "\\;"))[[2]]
}

unique(vcf$type)
summary(as.factor(vcf$type))

#### VGSC ####

vgsc <- vcf[grep("LOC5567355", vcf$V8),]

for (i in 1:nrow(vgsc)){
    vgsc$type[i] <- unlist(strsplit(vgsc$V8[i], "\\|"))[[2]]
    vgsc$impact[i] <- unlist(strsplit(vgsc$V8[i], "\\|"))[[3]]
}

vgsc_high <- vgsc[vgsc$impact== "HIGH",]
unique(vgsc_high$type)
unique(vgsc_high$V2)


### RDL ####

rdl <- vcf[grep("LOC5570466", vcf$V8),]

for (i in 1:nrow(rdl)){
    rdl$type[i] <- unlist(strsplit(rdl$V8[i], "\\|"))[[2]]
    rdl$impact[i] <- unlist(strsplit(rdl$V8[i], "\\|"))[[3]]
}

rdl_high <- rdl[rdl$impact== "HIGH",]
unique(rdl_high$type)
unique(rdl_high$V2)



### ACE - LOC5578456
ace <- vcf[grep("LOC5578456", vcf$V8),]

for (i in 1:nrow(ace)){
    ace$type[i] <- unlist(strsplit(ace$V8[i], "\\|"))[[2]]
    ace$impact[i] <- unlist(strsplit(ace$V8[i], "\\|"))[[3]]
}

ace_high <- ace[ace$impact== "HIGH",]
unique(ace_high$type)
unique(ace_high$V2)

### GSTE - LOC110676855
gste <- vcf[grep("LOC110676855", vcf$V8),]

for (i in 1:nrow(gste)){
    gste$type[i] <- unlist(strsplit(gste$V8[i], "\\|"))[[2]]
    gste$impact[i] <- unlist(strsplit(gste$V8[i], "\\|"))[[3]]
}

gste_high <- ace[gste$impact== "HIGH",]
unique(ace_high$type)
unique(ace_high$V2)


### P450 DUPS?

cyp <- vcf[grep("LOC5565578", vcf$ANN),]

