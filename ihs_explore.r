
library(dplyr)

ihs <- read.table("/mnt/storage12/emma/ihs/all_ihs_edit.tsv", header = TRUE, sep = "\t", row.names = NULL)
head(ihs)
summary(ihs$iHS.Value)
ihs_group <- ihs %>% group_by(Position,Chromosome) %>% summarise(n = n())
summary(ihs_group$n)

ir_genes <- c("LOC5567355", "LOC5570466", "LOC110676855", "LOC5578456")
ihs_sub <- ihs[grep(ir_genes, ihs$Gff_Annotation),]

ihs_sub <- ihs[grep("LOC5567355", ihs$Gff_Annotation),]
ihs_sub <- ihs[grep("LOC5570466", ihs$Gff_Annotation),]
ihs_sub <- ihs[grep("LOC110676855", ihs$Gff_Annotation),]
ihs_sub <- ihs[grep("LOC5578456", ihs$Gff_Annotation),]


## IR gene
vgsc <- ihs[ihs$Position > 315926360 & ihs$Position < 316405639,]
rdl <- ihs[ihs$Position > 41628484 & ihs$Position < 41861946,]
ace <- ihs[ihs$Position > 161486025 & ihs$Position < 161871375,]
gste <- ihs[ihs$Position > 351633367 & ihs$Position < 351634957,]


## Genes of interest
interest <- read.table("/mnt/storage12/emma/Reference_files/genes_of_interest.txt")

library(stringr)
list <- c()
position <- c()
for (j in 1:nrow(ihs)){
    gene <- unique(unlist(str_extract_all(ihs$Gff_Annotation[j], "LOC[0-9]+")))
    #interest_gene <- str_extract_all(gene, "LOC[0-9]+") 
    for (i in unique(gene)){
        if(i %in% interest$V1){
        list <- c(list, i)
        position <- c(position, ihs$Position[j])
        }
    }
}

top_ihs <- ihs %>% top_n(iHS.Value, n = 20)

top_ihs_interest <- ihs %>% filter(!is.na(gene)) %>% top_n(iHS.Value, n = 20)
## why is LOC5566204 not in interest list

list
#grep("LOC5578894", interest$V1)
#ihs$Position[grep("LOC5566204", ihs$Gff_Annotation)]
#ihs[ihs$Position == "232928340",]

ihs_ir_genes <- ihs[ihs$Position %in% position,]

ihs$gene <- NA
ihs$product <- NA

for (j in 1:nrow(ihs)){
    gene <- unique(unlist(str_extract_all(ihs$Gff_Annotation[j], "LOC[0-9]+")))
    product <- unique(unlist(str_extract_all(ihs$Gff_Annotation[j], "product=[^;]+")))
    ihs$gene[j] <- gene
    ihs$product[j] <- product
}


interest_ihs <- ihs[ihs$gene %in% interest$V1,]

pr <- ihs[ihs$Country == "Puerto_Rico",]
pr <- pr[,c(1:4,6:7)]




##############################################

library(dplyr)

ihs <- read.table("/mnt/storage12/emma/ihs/all_ihs_edit.tsv", header = TRUE, sep = "\t", row.names = NULL)
head(ihs)

length(unique(ihs$Position))
table(ihs$Chromosome)
summary(ihs$iHS.Value)

ihs$gene <- NA
ihs$product <- NA

# Loop over rows to extract gene and product information
for (j in 1:nrow(ihs)) {
  gene_matches <- unique(unlist(str_extract_all(ihs$Gff_Annotation[j], "LOC[0-9]+")))
  product_matches <- unique(unlist(str_extract_all(ihs$Gff_Annotation[j], "product=[^;]+")))
  
  # Assign non-empty matches to the data frame
  if (length(gene_matches) > 0) {
    ihs$gene[j] <- paste(gene_matches, collapse = "; ")
  } else {
    ihs$gene[j] <- NA
  }
  
  if (length(product_matches) > 0) {
    ihs$product[j] <- paste(product_matches, collapse = "; ")
  } else {
    ihs$product[j] <- NA
  }
}


length(unique(ihs$gene))
# 183 genes

top_ihs <- ihs %>% top_n(iHS.Value, n = 10)
top_ihs2 <- top_ihs[,c(1:4,6:7)]

interest_bed <- read.table("/mnt/storage12/emma/Reference_files/of_interest_no_contig_with_descriptions.bed", sep = "\t", header = FALSE)
head(interest_bed)

## add gene name column to interest file
interest_bed$gene <- NA
for ( j in 1:nrow(interest_bed)){
    interest_bed$gene[j] <- unique(unlist(str_extract_all(interest_bed$V4[j], "LOC[0-9]+")))
}
head(interest_bed)


ihs_interest2<- ihs[ihs$gene %in% interest_bed$gene,]
#ihs_interest2<- subset(ihs, gene %in% interest_bed$gene)
unique(interest_bed$gene)

## LOC5566204 gene missing from interest_bed - because in gff it doesnt have a gene line and I subset by genes

grep("LOC5566204", interest_bed$gene)
ihs[grep("LOC5566204", ihs$gene),]


# Function to check if any gene in ihs$gene is in interest_bed$gene
is_gene_of_interest <- function(gene_string, interest_genes) {
  genes <- unlist(str_split(gene_string, ";\\s*"))
  any(genes %in% interest_genes)
}

# Apply the function to subset the data frame
ihs_interest2 <- ihs %>% filter(sapply(gene, is_gene_of_interest, interest_genes = interest_bed$gene))
length(unique(ihs_interest2$gene))


### Looks at Puerto Rico
pr <- ihs[ihs$Country == "Puerto_Rico",]
pr <- pr[,c(1:4,6:7)]
length(unique(pr$gene))

## LOC5566204 gene missing from interest_bed - because in gff it doesnt have a gene line and I subset by genes
## Fixed by using regions of interest found from original grep rather than subset version
regions_gff <- read.table("/mnt/storage12/emma/Reference_files/region_of_interest.gff", sep = "\t", header = FALSE)

head(regions_gff)
regions_gff <- regions_gff[regions_gff$V1 %in% c("035107", "035108", "035109", "013159"),]
for ( i in 1:nrow(regions_gff)){
    regions_gff$gene[i] <- unique(unlist(str_extract_all(regions_gff$V9[i], "LOC[0-9]+")))
}

unique(regions_gff$gene)
length(unique(regions_gff$gene))
#643
# 616 no contigs


# Apply the function to subset the data frame
ihs_interest3 <- ihs %>% filter(sapply(gene, is_gene_of_interest, interest_genes = regions_gff$gene))
length(unique(ihs_interest2$gene))




###############################################################
# JUNE 2024
###############################################################

library(dplyr)
## pavlue
#ihs_sig_p <- read.table("/mnt/storage12/emma/ihs/all_ihs_significant_4.csv", header = TRUE, sep = ",", row.names = NULL)
#head(ihs_sig_p)

#ihs_all_p <- read.table("/mnt/storage12/emma/ihs/all_ihs_all_values.csv", header = TRUE, sep = ",", row.names = NULL)
#head(ihs_all_p)

ihs_all <- read.table("/mnt/storage12/emma/ihs/all_ihs_value_and_pvalue.csv", header = TRUE, sep = ",", row.names = NULL)
unique(ihs_all$country)
head(ihs_all)
ihs_group <- ihs_all %>% dplyr::group_by(chrom,pos) %>% summarise(n= n())
summary(ihs_group$n)
ihs_group_top <- ihs_group[ihs_group$n >7,]
length(unique(ihs_group_top$pos))

ihs_sig2 <- ihs_all[ihs_all$pval >4,]
length(unique(ihs_sig2$pos))
table(ihs_sig2$chrom)

ihs_group_sig <- ihs_all %>% filter(pval >4) %>% dplyr::group_by(chrom,pos) %>% summarise(n= n())

ihs_group_top <- ihs_group_sig[ihs_group_sig$n >1,]
nrow(ihs_group_top)



