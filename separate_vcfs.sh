#!/bin/bash
 
# Define input VCF filename
input_vcf="all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.vcf.gz"
 
# Loop through each population in fst_populations.txt
cat countries.txt | xargs -I {} -P 10 sh -c "
    population='{}'
    output_vcf=\"\${population%.*}.vcf.gz\"
    # Extract samples for the population
    bcftools view -S \"\${population}\" \"${input_vcf}\" -o \"\${output_vcf}\"
    echo \"Population \${population} extracted to \${output_vcf}\""
