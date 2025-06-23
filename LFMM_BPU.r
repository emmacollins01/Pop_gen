
########################################
######## R CODE: latent factor mixed models
########################################
library("LEA")
setwd("/mnt/storage12/emma/LFMM/test_run/")


# Aedes aegypti
# F1534C -  an example
## first, run sNMF to find appropriate K
#vcf2geno(input.file = "/mnt/storage12/emma/LFMM/PR_merge_chromonly_phase.vcf", output.file = "admix.geno", force = TRUE)
#vcf2geno(input.file = "/mnt/storage12/emma/LFMM/PR_merged_chromonly_phase.vcf", output.file = "admix.geno", force = TRUE)
vcf2geno(input.file = "/mnt/storage12/emma/LFMM/test_run/BPU_phase.vcf", output.file = "admix.geno", force = TRUE)
#vcf2geno(input.file = "/mnt/storage12/emma/LFMM/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.snps_phase.vcf", output.file = "admix.geno", force = TRUE)

## seems unecessary 
#admix_pop=read.delim("vssc395ids.txt",row.names = 1,header = TRUE)

## Takes awhile to run - RUN IN SCREEN - screen -r 890059.pts-17.plum-s12
project = NULL
project = snmf(input.file = "admix.geno", K = 1:25, entropy = TRUE,
               repetitions = 10, CPU = 16, project = "new")


project = load.snmfProject("admix.snmfProject")
#plot(project, col = "blue", pch = 19, cex = 1.2)
png("admixture_plot.png", width = 800, height = 600)
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()
## need to adjust K??
best = which.min(cross.entropy(project, K = 4))

#######################################
############# BASH CODE ###############
#######################################

# creating the .env file for a specific position
## Going to use mutation F1534C - position 315939224
# the env file is a vector file which contains information on the number of copies of a resistance allele in each sample, 0, 1 or 2.
# first, this command extracts the genotypes at position from the VCF and saves them into a file called genotypes.txt.
cp BPU_phase.vcf BPU_phase_copy.vcf
bgzip BPU_phase_copy.vcf
tabix BPU_phase_copy.vcf.gz
bcftools query -r 035109:315939224-315939224 -f '[%GT\t]\n' BPU_phase_copy.vcf.gz > F1534C_BPU_genotypes.txt
bcftools query -r 035109:316080722-316080722 -f '[%GT\t]\n' BPU_phase_copy.vcf.gz > V410L_BPU_genotypes.txt

# then convert the genotype codes into numeric allele counts
sed -i 's/0|0/0/g; s/0|1/1/g; s/1|0/1/g; s/1|1/2/g' F1534C_BPU_genotypes.txt
# Replace tabs (or spaces) with newlines to ensure each value is on its own line
tr '\t' '\n' < F1534C_BPU_genotypes.txt > F1534C_BPU_new6.env
# Verify the resulting .env file
head F1534C_BPU.env
# this file  is then a vector where each row corresponds to a sample, 
#and the values represent the number of copies of the resistance allele for that SNP (position 2422652)


########################################
######## R CODE: latent factor mixed models
########################################

## second, run lfmm using this K. the *.env file 
# is a vector of values 0-2 indicating number of copies of the resistance allele
library("LEA")
F1534C_BPU_new5.env <- read.table("/mnt/storage12/emma/LFMM/test_run/F1534C_BPU_new6.env", quote="\"", comment.char="")
V410L_BPU.env <- read.table("/mnt/storage12/emma/LFMM/test_run_V410L/V410L_BPU_4.env", quote="\"", comment.char="")

vssc.lfmm = vcf2lfmm("/mnt/storage12/emma/LFMM/test_run/BPU_phase.vcf")
vssc_v2.lfmm = vcf2lfmm("/mnt/storage12/emma/LFMM/test_run_V410L/BPU_phase.vcf")

## RUN IN SCREEN - 1442179.pts-21.plum-s12	(17/12/24 16:35:01)	(Detached)
F1534C.lfmm = lfmm("/mnt/storage12/emma/LFMM/test_run/BPU_phase.lfmm", 
                    "/mnt/storage12/emma/LFMM/test_run/F1534C_BPU_new6.env", 
                    K = 3, CPU = 16, rep = 10, project="force", 
                    iterations = 1000, burnin = 50, all = FALSE, missing.data = TRUE)

V410L.lfmm = lfmm("/mnt/storage12/emma/LFMM/test_run_V410L/BPU_phase.lfmm", 
                    "/mnt/storage12/emma/LFMM/test_run_V410L/V410L_BPU_4.env", 
                    K = 3, CPU = 16, rep = 10, project="force", 
                    iterations = 1000, burnin = 50, all = FALSE, missing.data = TRUE)

#project = load.lfmmProject("PR_merged_chromonly_phase_F1534C.lfmmProject")
project = load.lfmmProject("/mnt/storage12/emma/LFMM/test_run/BPU_phase_F1534C_BPU_new6.lfmmProject")
project = load.lfmmProject("/mnt/storage12/emma/LFMM/test_run_V410L/BPU_phase_V410L_BPU_4.lfmmProject")

#Record z-scores from the 10 runs in the zs matrix
zs_GTC = z.scores(F1534C.lfmm, K = 3, d = 1)
zs_GTC = z.scores(V410L.lfmm, K = 3, d = 1)

#Combine z-scores using the median
zs.median = apply(zs_GTC, MARGIN = 1, median)

#Compute the GIF
lambda = median(zs.median^2)/qchisq(0.5, df = 1)
lambda

# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values, col = "red")

## FDR control: Benjamini-Hochberg at level q
## L = number of loci
L = 641034
#fdr level q
q = 0.01
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]
candidates.bh

## save outputs
write.csv(adj.p.values, "adj_k3_V410L_BPU")
saveRDS(candidates.bh, file = "candidates_bh_V410L.rds")

pval <- read.csv("adj_k3_V410L_BPU")

# load candidates.bh with:
#big_candidates.bh <- readRDS("big_candidates_bh.rds")
## run in bash 
## bcftools query -f '%CHROM\t%POS\n' BPU_phase.vcf > BPU_chrom_pos.tsv
#vcf_pos <- read.delim("BPU_chrom_pos.tsv", sep = "\t", header = FALSE)
vcf_pos <- read.delim("BPU_phase.vcfsnp", sep = " ", header = FALSE)
head(vcf_pos)
colnames(vcf_pos) <- c("CHROM", "POS", "Q", "REF", "ALT", "?", "QUAL", "GT")
sig_pos <- vcf_pos[candidates.bh,]

write.csv(sig_pos, "BPU_V410L_sig_chrom_pos.csv")

## all points 
## check if nrow in vcf_pos = adj pvals
nrow(vcf_pos) == nrow(pval)
vcf_pos$p_value <- pval$pval

# Add p-values to these significant positions
adj_sig <- pval[candidates.bh,]
colnames(adj_sig) <- c("index", "pval")
# Add them as a new column
sig_pos$p_value <- adj_sig$pval

write.csv(sig_pos, "BPU_V410L_pos_pval.csv")


library(ggplot2)
library(dplyr)

ggplot()+
geom_point(data = sig_pos, aes(x = POS, y = p_value), colour = "red") +
#geom_point(data = vcf_pos, aes(x = POS, y = p_value), colour = "black", alpha = 0.5) +
facet_wrap(~CHROM, nrow = 3) +
theme_minimal()


library(dplyr)

sig_pos <- read.csv("/mnt/storage12/emma/LFMM/test_run/BPU_pos_pval.csv")
sig_pos <- read.csv("/mnt/storage12/emma/LFMM/test_run_V410L/BPU_V410L_pos_pval.csv")
sig_pos <- read.csv("/mnt/storage12/emma/LFMM/test_run/BPU_sig_chrom_pos.csv")

## gaba_alpha NC_035107.1 (28030989..28161637)
sig_pos$POS <- as.numeric(sig_pos$POS)
subset <- sig_pos[sig_pos$CHROM == "35107" & sig_pos$POS > "28030989",]
subset2 <- subset[subset$POS < "28161637",]

test <- sig_pos[sig_pos$CHROM == "35107",]
summary(test$POS)
test2 <- test[test$POS > "28030989" & test$POS < "28161637",]

# vgsc 3:315926360-316405639

vgsc <- sig_pos[sig_pos$CHROM == 35109 & sig_pos$POS > 315926360 & sig_pos$POS < 316405639,]


top_10 <- sig_pos %>% top_n(10, p_value)

## Add gene to sig pos 

# Required libraries
library(readr)
library(dplyr)

# Input and output file paths
input_csv <- "/mnt/storage12/emma/LFMM/test_run_V410L/BPU_V410L_pos_pval.csv"
output_file <- "/mnt/storage12/emma/LFMM/test_run_V410L/V410L_sig_pos_output_genes.txt"

data <- read_csv(input_csv)


# Initialize an empty data frame to store results
result_df <- data.frame(
  Pos = numeric(),
  Region = character(),
  Gene = character(),
  Product = character(),
  Pvalue = numeric(),
  stringsAsFactors = FALSE
)

# Iterate through each row
for (i in seq_len(nrow(data))) {
  row <- data[i, ]
  chrom <- row$CHROM
  pos <- row$POS
  bin_start <- as.integer(row$POS)
  bin_end <- bin_start
  pavlue <- row$p_value
  
  # Format CHROM to match the GTF file format
  formatted_chrom <- sprintf("NC_0%s.1", chrom)
  
  # Construct the tabix command
  region <- sprintf("%s:%d-%d", formatted_chrom, bin_start, bin_end)
  command <- sprintf("tabix /mnt/storage12/emma/Reference_files/genes_sorted.gff.gz %s", region)
  
  # Run the tabix command and capture the output
  output <- tryCatch({
    system(command, intern = TRUE)  # Run command and capture output
  }, error = function(e) {
    message <- sprintf("Error processing region %s: %s", region, e$message)
    warning(message)
    return(NULL)
  })
  
  if (!is.null(output) && length(output) > 0) {
    for (line in output) {
      fields <- strsplit(line, "\t")[[1]]
      if (length(fields) > 8) {
        attributes <- fields[9]
        match <- regmatches(attributes, regexpr("gene=([^;]+)", attributes))
        gene_name <- ifelse(length(match) > 0, sub("gene=", "", match), "N/A")
        match <- regmatches(attributes, regexpr("product=([^;]+)", attributes))
        product <- ifelse(length(match) > 0, sub("product=", "", match), "N/A")

        # Add the processed line to the data frame
        result_df <- bind_rows(
          result_df,
          data.frame(
            Region = region,
            Pos = pos,
            Gene = gene_name,
            Product = product,
            Pvalue = pavlue,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}


test <- result_df[!duplicated(result_df),]
test2 <- test[test$Gene != "N/A",]

unique(test2$Gene)


gff <- read.csv("/mnt/storage12/emma/Reference_files/genes_only.gff", sep = "\t", header = FALSE)
head(gff)
library(stringr)
# function to extract product from gff
extract_product <- function(ann) {
  product <- str_extract_all(ann, "product=[^;]*;")[[1]]  # Extract all matches
  if (length(product) == 0) {
    return(NA)  # Return NA if no matches found
  } else {
    return(paste(unique(product), collapse = "; "))  # Concatenate unique matches
  }
}

# loop through and extract product to match to gene name

for (j in test2$Gene){
    for (i in 1:length(gff$V9[grep(test2$Gene[j], gff$V9)])){
        ann <- gff$V9[grep(test2$Gene[j], gff$V9)][i]
        ann2 <- str_extract_all(ann, "product=[^;]*;")
        test2$Product[test2$Gene == j] <- ann2[1]
    }
}


ann <- gff$V9[grep(test2$Gene[1], gff$V9)]
ann2 <- list()
for (i in unlist(ann)){
  print(i)
  ann2 <- str_extract_all(ann, "product=[^;]*;")
}
ann$product <- ann2$






extract_product(ann)


