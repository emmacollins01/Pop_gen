import pandas as pd
import re

#%%
subset = pd.read_csv("all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.lmiss.recode.ann.of.interest.txt")


#%%
# Assuming `subset` is a pandas DataFrame with an 'ANN' column containing annotation strings

# Filter rows containing 'missense_variant' in the 'ANN' column
missense = subset[subset['ANN'].str.contains("missense_variant")]

# Get unique positions from 'missense'
unique_pos = missense['POS'].nunique()

# Extract information from annotation
high = missense.copy()

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

# Now `high` DataFrame contains the extracted information