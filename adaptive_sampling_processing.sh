
## convert fast5 to pod5
ls fast5/*.fast5 | parallel -j 10 pod5 convert fast5 --output pod5/ --one-to-one ./fast5/ {}

## basecall - pod5 to fastq
conda activate dorado
dorado basecaller hac pod5/ > base_perc.bam


## mapping 
conda activate minion
minimap2 -ax map-ont ../Reference_files/GCF_002204515.2_AaegL5.0_genomic.2023.fa base_perc.fastq.gz | samtools view -S -b - | samtools sort - -o base_perc.bam | samtools index base_perc.bam

