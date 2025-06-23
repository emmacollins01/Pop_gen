

## index individual vcfs
tabix Uganda_subset_SRR.vcf.gz

## Combine all vcfs

bcftools merge Puerto_Rico_subset_SRR.vcf.gz Uganda_subset_SRR.vcf.gz Brazil_subset_SRR.vcf.gz -Oz -o BPU.vcf.gz

# index merged vcf
tabix BPU.vcf.gz

## phase vcf
java -Xmx50G -jar beagle.01Mar24.d36.jar gt=BPU.vcf.gz out=BPU_phase

