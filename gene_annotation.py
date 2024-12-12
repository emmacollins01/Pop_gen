#%%
import pandas as pd
import re

#%%
subset = pd.read_csv("/mnt/storage12/emma/PR_snps/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest.txt", sep = "\t", header = None)

#%%
# all - be cautious!!
subset = pd.read_csv("/mnt/storage12/emma/PR_snps/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.missense_only.txt", sep = "\t", header = None)

#%%
subset.columns = ["SAMPLE", "CHROM", "POS", "REF", "ALT", "QUAL", "GT", "ALLELE", "DP", "AD", "ANN"]

# Remove positions without a call
subset = subset[subset['GT'] != './.']
# %%

#%%
# Assuming `subset` is a pandas DataFrame with an 'ANN' column containing annotation strings

# Filter rows containing 'missense_variant' in the 'ANN' column
missense = subset[subset['ANN'].str.contains("missense_variant")]

# Get unique positions from 'missense'
unique_pos = missense['POS'].nunique()

# Extract information from annotation
high = missense.copy()

#%% 
high_test = missense[missense['ANN'].str.contains("HIGH")]

#%%
for i, row in high.iterrows():
    annotations = row['ANN'].split(',')
    missense_index = next((i for i, ann in enumerate(annotations) if "missense_variant" in ann), None)
    
    if missense_index is not None:
        to_extract = annotations[missense_index]
        parts = to_extract.split('|')
        high.at[i, 'impact'] = parts[1] if len(parts) > 1 else ''
        high.at[i, 'effect'] = parts[2] if len(parts) > 2 else ''
        high.at[i, 'gene'] = parts[4] if len(parts) > 4 else ''
        high.at[i, 'csq'] = parts[10] if len(parts) > 10 else ''
        if 'LOC' in parts:
            gene_name = parts.split('|')[-4]
            if gene_name.startswith('LOC'):
                gene_name.append(gene_name)
                high.at[i, 'gene'] = gene_name
# Now `high` DataFrame contains the extracted information

#%%
high.to_csv("/mnt/storage12/emma/PR_snps/of_interest_snps_only_gene_annotated.csv")


#%%
# Read GFF file
gff = pd.read_csv("/mnt/storage12/emma/Reference_files/genes.gff", header=None, comment='#', sep="\t")

# Assuming `high` is a pandas DataFrame containing the 'gene' column
#%%
for i, row in high.iterrows():
    exon_name = row['gene'][10:]  # Extract exon name
    line = gff[gff[8].str.contains(exon_name)][8].iloc[0]  # Find the line containing exon_name
    gene_id = line.split(';')[4]  # Extract gene ID from the line
    high.at[i, 'gene_id'] = gene_id[5:]  # Extract gene ID value

# `high` now contains the extracted gene IDs
#%%
high.to_csv("/mnt/storage12/emma/PR_snps/of_interest_snps_only_gene_annotated_LOC.csv")
# %%

