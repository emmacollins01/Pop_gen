
library(rehh)
library(ggplot2)

all <- data2haplohh(hap_file = "/mnt/storage12/emma/selection/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.main_chrom_only.phased.vcf.gz",
                   polarize_vcf = FALSE,
                   chr.name = "35107")


meta <-read.csv("/mnt/storage12/emma/selection/all_sample_country_metadata.csv")

test <- subset(all, min_af = 0.05)

# read in data for each species


# house
chroms <- c("35017", "35108", "35109")
all_PR <- data.frame()
for (i in chroms){
    chr_name <-i
    PR_hh <- data2haplohh(hap_file = "/mnt/storage12/emma/selection/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.main_chrom_only.phased_PuertoRico.vcf.gz",
                   polarize_vcf = FALSE,
                   chr.name = chr_name)
    all_PR <- rbind(PR_hh, all_PR)
}

for(i in chroms) {
  # haplotype file name for each chromosome
  hap_file = paste("hap_chr_", i, ".cgu", sep = "")
  # create internal representation
  hh <- data2haplohh(hap_file = hap_file,
                     chr.name = i,
                     allele_coding = "map")

PR_hh <- data2haplohh(hap_file = "/mnt/storage12/emma/selection/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.main_chrom_only.phased_PuertoRico.vcf.gz",
                   polarize_vcf = FALSE,
                   chr.name = "35108")
# bactrianus
Mexico_hh <- data2haplohh(hap_file = "/mnt/storage12/emma/selection/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.main_chrom_only.phased_Mexico.vcf.gz",
                         polarize_vcf = FALSE,
                         chr.name = "35109")

Kenya_hh <- data2haplohh(hap_file = "/mnt/storage12/emma/selection/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.main_chrom_only.phased_Kenya.vcf.gz",
                         polarize_vcf = FALSE,
                         chr.name = "35109")

PR_hh_f <- subset(PR_hh, min_maf = 0.05)
Mexico_hh_f <- subset(Mexico_hh, min_maf = 0.05)
Kenya_hh_f <- subset(Kenya_hh, min_maf = 0.05)

## iHS

# perform scans
PR_scan <- scan_hh(PR_hh_f, polarized = FALSE)
Mexico_scan <- scan_hh(Mexico_hh_f, polarized = FALSE)
Kenya_scan <- scan_hh(Kenya_hh_f, polarized = FALSE)

# perform iHS on house
PR_ihs <- ihh2ihs(PR_scan, freqbin = 1)
head(PR_ihs$ihs)
summary(PR_ihs$ihs)

head(PR_ihs$frequency.class)

Mexico_ihs <- ihh2ihs(Mexico_scan, freqbin = 1)
head(Mexico_ihs$frequency.class)

Kenya_ihs <- ihh2ihs(Kenya_scan, freqbin = 1)
head(Kenya_ihs$frequency.class)

ggplot(PR_ihs$ihs, aes(POSITION, IHS)) + geom_point()

hit <- house_bac %>% arrange(desc(LOGPVALUE)) %>% top_n(1)
top <- PR_ihs$ihs %>% arrange(desc(LOGPVAULE)) %>% top_n(10)

top <- PR_ihs$ihs[PR_ihs$ihs$LOGPVALUE >6,]

ggplot() +
geom_point(data =PR_ihs$ihs, aes(POSITION, IHS), colour = "red") +
geom_point(data =Mexico_ihs$ihs, aes(POSITION, IHS), colour = "blue") +
geom_point(data =Kenya_ihs$ihs, aes(POSITION, IHS), colour = "darkgreen") 

ggplot() +
geom_point(data =PR_ihs$ihs, aes(POSITION, LOGPVALUE), colour = "red") +
geom_point(data =Mexico_ihs$ihs, aes(POSITION, LOGPVALUE), colour = "blue") +
geom_point(data =Kenya_ihs$ihs, aes(POSITION, LOGPVALUE), colour = "darkgreen") +
ylim(c(6,10))


# plot
ggplot(PR_ihs$ihs, aes(POSITION, LOGPVALUE)) + 
geom_point() +
ylim(c(6,10))

ggplot() +
geom_point(data =PR_ihs$ihs, aes(POSITION, LOGPVALUE), colour = "red") +
geom_point(data =Mexico_ihs$ihs, aes(POSITION, LOGPVALUE), colour = "blue") +
geom_point(data =Kenya_ihs$ihs, aes(POSITION, LOGPVALUE), colour = "darkgreen")

# perform xp-ehh
PR_mex <- ies2xpehh(PR_scan, Mexico_scan,
                       popname1 = "PuertoRico", popname2 = "Mexico",
                       include_freq = T)
## Scan of pop1 contains 189364 markers.
## Scan of pop2 contains 180301 markers.
## Merged data contains 162632 markers.


PR_ken <- ies2xpehh(PR_scan, Kenya_scan,
                       popname1 = "PuertoRico", popname2 = "Kenya",
                       include_freq = T)
## Scan of pop1 contains 189364 markers.
## Scan of pop2 contains 145753 markers.
## Merged data contains 130590 markers.

Mex_ken <- ies2xpehh(Mexico_scan, Kenya_scan,
                       popname1 = "PuertoRico", popname2 = "Kenya",
                       include_freq = T)
#Scan of pop1 contains 180301 markers.
#Scan of pop2 contains 145753 markers.
#Merged data contains 122708 markers.

ggplot(PR_mex, aes(POSITION, XPEHH_PuertoRico_Mexico)) + geom_point() +
    theme(axis.text = element_text(size = 20)) +
    ylim(c(-6, 6))

ggplot(PR_ken, aes(POSITION, XPEHH_PuertoRico_Kenya)) + geom_point() +
    theme(axis.text = element_text(size = 20)) +
    ylim(c(-6,6))

ggplot(PR_mex, aes(POSITION, LOGPVALUE)) + geom_point() +
    theme(axis.text = element_text(size = 20))

ggplot(PR_ken, aes(POSITION, LOGPVALUE)) + geom_point() +
    theme(axis.text = element_text(size = 20))


ggplot() +
geom_point(data = PR_mex, aes(POSITION, LOGPVALUE), colour = "blue") +
geom_point(data = PR_ken, aes(POSITION, LOGPVALUE), colour = "red")


# find the highest hit
hit <- PR_mex %>% arrange(desc(LOGPVALUE)) %>% top_n(1)
# get SNP position
x <- hit$POSITION
x

marker_id_PR <- which(PR_hh_f@position == x)
marker_id_Mex <- which(Mexico_hh_f@POSITION == x)
