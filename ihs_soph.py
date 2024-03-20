######################## SELECTION STATISTICS #########################

## SELECT SAMPLE POPULATION TO WORK WITH

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.
# Note that the zarr file should have been made from a phased vcf
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# %%
import os
import numpy as np
import allel
import zarr
import pandas as pd
import gffutils
import tqdm

# %% set wd
os.chdir('/mnt/storage12/emma/ihs/')
os.getcwd()

# %%
# convert phased, filtered, VCF file to zarr file
# already converted to zarr
#allel.vcf_to_zarr('../selection/all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.main_chrom_only.phased.vcf.gz', 'all_aedes_lmiss_main_chrom_only_phased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('all_aedes_lmiss_main_chrom_only_phased.zarr', mode='r')
#
# callset.tree(expand=True)

# %%
## convert zarr file to genotype array
gt = allel.GenotypeDaskArray(callset['calldata/GT'])
print(gt.shape)

# %%
## import metadata
df_samples=pd.read_csv('all_sample_country_metadata.csv',sep=',',usecols=['ID','Country','Region','Subregion'])
df_samples.head()
df_samples.groupby(by=['Country']).count

# %%
## working with Guinea-Bissau samples
sample_ids = callset['samples'][:]
# Get sample identifiers for Cameroon samples from df_samples
gb_sample_ids = df_samples[df_samples['Country'] == 'Puerto Rico']['ID'].values
# Find indices of these samples in the genotype array
gb_indices = np.array([np.where(sample_ids == id)[0][0] for id in gb_sample_ids if id in sample_ids])
# Verify the indices are within the correct range
print("Max index:", gb_indices.max(), "Sample array size:", len(sample_ids))
# Select genotypes for Cameroon samples using the indices
gt_gb_samples = gt.take(gb_indices, axis=1)
gt_gb_samples

# %%
## select variants that are segregating within gb_samples as only these will be informative
## also some selection tests don't support multiallelic variants, so just keep biallelics
## for this pipeline the VCF is already filtered so should be no biallelic SNPs anyway

ac_gb = gt_gb_samples.count_alleles(max_allele=8).compute()
gb_seg_variants = ac_gb.is_segregating() & ac_gb.is_biallelic_01()
ac_gb_seg = ac_gb.compress(gb_seg_variants, axis=0)
gt_gb_seg = gt_gb_samples.compress(gb_seg_variants, axis = 0)
gt_gb_seg

# %%
## this is from a phased VCF so we can convert this genotype array to haplotype array

h_gb_seg = gt_gb_seg.to_haplotypes().compute()
h_gb_seg

# %%
# we need variant positions
pos = callset['variants/POS'][:]
pos_gb_seg = pos.compress(gb_seg_variants, axis=0)
pos_gb_seg

# %%
# some variants in 1000 genomes project have multiple variants at the same genomic position, 
# which causes problems for some selection tests in scikit-allel. 
# Let's check if there any of these.
count_multiple_variants = np.count_nonzero(np.diff(pos_gb_seg == 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    #sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %%
# compute raw iHS

ihs_gb_raw = allel.ihs(h_gb_seg, pos_gb_seg, use_threads=True, include_edges=True)
ihs_gb_raw
print("Raw iHS computed")

# %%

#%matplotlib inline
import matplotlib.pyplot as plt
from datetime import datetime

# %% view raw iHS values as a histogram
# ~np.isnan(ihs_gb_std[0]) is used to filter out NaN values
fig, ax = plt.subplots()
ax.hist(ihs_gb_raw[~np.isnan(ihs_gb_raw)], bins=20)
ax.set_xlabel('Raw IHS')
ax.set_ylabel('Frequency (no. variants)');

# %% Standardize iHS

ihs_gb_std = allel.standardize_by_allele_count(ihs_gb_raw, ac_gb_seg[:, 1])
ihs_gb_std
print("Standardized iHS computed")

# %% 

# Here we deviate from the Jupyter notebook and use ihs_res_std[0]
# ~np.isnan(ihs_gb_std[0]) is used to filter out NaN values
fig, ax = plt.subplots()
ax.hist(ihs_gb_std[0][~np.isnan(ihs_gb_std[0])], bins=20)
ax.set_xlabel('Standardised IHS')
ax.set_ylabel('Frequency (no. variants)');

# Save the figure as a file (e.g., PNG) in the current working directory
filename = f'standardised_ihs_histogram.png'
plt.savefig(filename)

# show the plot (optional, could # this out)
plt.show()

# %% note that iHS has been calculated with unpolarized data using ~np.isnan, so only the magnitude of iHS
# is informative, not the sign.

# plot over the genome
# np.abs is converting all iHS vales to their absoltue values before plotting. This means that if ihs_gb_std[0]
# contains any negative valyes, those values will be made positive. It is plotting the magnitude of iHS without considering the
# direction of selection, which the sign of iHS could indicate

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos_gb_seg, np.abs(ihs_gb_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
ax.set_xlabel('Genomic position (bp) chromosome agnostic')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(-2, 9);

# Save the figure as a file (e.g., PNG) in the current working directory
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
filename = f'ihs_manhattan_{timestamp}.png'
plt.savefig(filename)

print("iHS plotted")

# %% find the index of the variant with the highest iHS value
idx_hit_max = np.nanargmax(ihs_gb_std[0])

# %% genomic position of top hit
pos_gb_seg[idx_hit_max]
print(f'Genomic position with highest iHS value (chr agnostic):', pos_gb_seg[idx_hit_max])

# %% Visualise EHH decay around top hit with a Voight plot
# pull out haplotypes for region around top hit
flank_size = 2000
h_hit = h_gb_seg[idx_hit_max - flank_size:idx_hit_max + flank_size + 1]
h_hit

# %%
fig = allel.fig_voight_painting(h_hit[:, h_gb_seg[idx_hit_max] == 0], index=flank_size, height_factor=0.02)
fig.suptitle('Reference allele', y=1);

# %% 
fig = allel.fig_voight_painting(h_hit[:, h_gb_seg[idx_hit_max] == 1], index=flank_size, height_factor=0.02)
fig.suptitle('Alternate allele', y=1);

print("EHH decay computed")

# %% Plot iHS manhattan plot with all chromosomes

# length of gb_seg_variants = 6785669, it is a numpy.ndarray
chromosomes = callset['variants/CHROM'][:]
chrom_gb_seg = chromosomes.compress(gb_seg_variants, axis = 0)
chrom_gb_seg
# length of chrom_gb_seg = 5673213, it is a numpy.ndarray
pos = callset['variants/POS'][:]
pos_gb_seg = pos.compress(gb_seg_variants, axis=0)
pos_gb_seg
# length of pos_gb_seg = 5673213, it is a numpy.ndarray

# define chromosome lengths and colours 
chromosome_lengths = {
    '035159': 16790,
    '035107': 310827022,
    '035108': 474425716,
    '035109': 409777670
}

# %% Calculate cumulative offsets for each chromosome
cumulative_lengths = {}
cumulative_length = 0
for chrom, length in chromosome_lengths.items():
    cumulative_lengths[chrom] = cumulative_length
    cumulative_length += length

# %% Plot iHS
# `pos`, `chromosomes`, and `ihs_gb_std[0]` arrays are already defined and aligned
# Define the threshold
import matplotlib.patches as mpatches

threshold = 5

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Define colors for each chromosome (for illustration)
chromosome_colours = {
    '035107': '#3d348b', '035108': '#f18701', '035109': '#f7b801', '035159': '#f35b04'
}

# Create a list to hold the legend patches
legend_patches=[]

# Filter out 'Y_unplaced' or any other chromosomes not in your chromosome_lengths dictionary
filtered_chroms = [chrom for chrom in sorted(set(chrom_gb_seg)) if chrom in chromosome_lengths]

# Iterate through each chromosome to plot its variants
for chrom in filtered_chroms:
    # Create mask for the current chromosome
    mask = chrom_gb_seg == chrom
    
    # Apply the chromosome mask and filter out NaN values simultaneously
    chrom_positions = pos_gb_seg[mask]
    chrom_ihs_values = ihs_gb_std[0][mask]
    non_nan_mask = ~np.isnan(chrom_ihs_values)
    
    # Make sure to apply the non-NaN mask to both the positions and iHS values
    chrom_positions_no_nan = chrom_positions[non_nan_mask]
    chrom_ihs_values_no_nan = chrom_ihs_values[non_nan_mask]
    
    # Adjust positions for visualization if needed
    adjusted_positions = chrom_positions_no_nan + cumulative_lengths[chrom]

    # Now create threshold masks based on the non-NaN iHS values
    below_threshold_mask = chrom_ihs_values_no_nan < threshold
    above_threshold_mask = chrom_ihs_values_no_nan >= threshold
    
    # Plot points below and above the threshold
    ax.scatter(adjusted_positions[below_threshold_mask], 
               chrom_ihs_values_no_nan[below_threshold_mask], 
               color=chromosome_colours[chrom], alpha=0.1, s=10)
    ax.scatter(adjusted_positions[above_threshold_mask], 
               chrom_ihs_values_no_nan[above_threshold_mask], 
               color=chromosome_colours[chrom], alpha=1.0, s=10)
    patch = mpatches.Patch(color=chromosome_colours[chrom], label=chrom)
    legend_patches.append(patch)

legend_patches.append(mpatches.Patch(color='black', label='Significance Threshold'))
ax.legend(handles=legend_patches, title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')

ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('$|iHS|$')
ax.axhline(y=threshold, color='black', linestyle='--')
plt.tight_layout()
plt.show()

print("iHS plotted with all chromosomes")


###########################################

# %% Remove NaN iHS values, and use the same mask to filter pos and chrom.
ihs_vals = ihs_gb_std[0] # the standardised ihs values are saved in ihs_gb_std[0]
mask_non_nan = ~np.isnan(ihs_vals) #remove the ihs_gb_std nan values
ihs_vals_withoutnan = ihs_vals[mask_non_nan] # save these values as vals_withoutnan
pos_withoutnan = pos_gb_seg[mask_non_nan] # use the same mask to get the corresponding positions of vals_withoutnan
chrom_withoutnan = chrom_gb_seg[mask_non_nan] # use the same mask to get the corresponding chromosomes of vals_withoutnan
print("Filtered out NaN values")

# %% Sort these by putting them in ascending order, ready for the calculate_empirical_p_value function
sorted_indices = np.argsort(ihs_vals_withoutnan) # np.argsort returns an array of teh corresponding indices
sorted_ihs = np.sort(ihs_vals_withoutnan) 
print("Put values into ascending order")

# %% Compute empirical p-values, then can log transform these and plot.

def calculate_empirical_p_value(val, sorted_vals,l):
    """
    Calculate the empirical p-value for an observed value in a sorted list of values.

    Parameters:
    - sorted_values: A list of values, sorted in ascending order.
    - observed_value: The observed value for which to calculate the p-value.
    - l: The length of the list of values

    Returns:
    - The empirical p-value.
    """
    return (l-np.where(sorted_vals>=val)[0][0])/l

# %% Calculate p-values for non-NaN iHS scores

print("Starting to calculate p-values of iHS")

# Multithread the p-value calculation because otherwise it is mega slow
import joblib
import tqdm

len_ihs = len(sorted_ihs)
pvals = []
from joblib import Parallel, delayed

parallel = Parallel(n_jobs=15, return_as='generator')
pvals = [r for r in tqdm.tqdm(parallel(delayed(calculate_empirical_p_value)(val,sorted_ihs,len_ihs) for val in sorted_ihs),total=len(sorted_ihs))]

# %% P-values are currently in the order of sorted iHS, need to reorder them back 
# based on the position of the original ihs value in ihs_vals_withoutnan array

reordered_pvals = np.empty(len(pvals))
reordered_pvals[sorted_indices] = pvals # reordered_pvals will now have the p-values in the original order of the ihs values

# %% Take log of p-values
neg_log_pvals = (-1*(np.log10(reordered_pvals)))
print("Computed log10 of p-values")

# %% Plotting adjusted p-values
print("Plotting neg log10 of p-values")

# Define colors for each chromosome (for illustration)
chromosome_colours = {
    '035107': '#3d348b', '035108': '#f18701', '035109': '#f7b801', '035159': '#7678ed'
}

# Plotting

threshold = 4

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))


# Create a list to hold the legend patches
legend_patches=[]

# Filter out 'Y_unplaced' or any other chromosomes not in your chromosome_lengths dictionary
filtered_chroms = ['035107', '035108', '035109', '035159']

# Iterate through each chromosome to plot its variants
for chrom in filtered_chroms:
    # Create mask for the current chromosome
    mask = chrom_withoutnan == chrom
    
    # Apply the chromosome mask
    chrom_positions = pos_withoutnan[mask]
    chrom_p_values = neg_log_pvals[mask]
    
    # Adjust positions for visualization if needed
    adjusted_positions = chrom_positions + cumulative_lengths[chrom]

    # Now create threshold masks based on the non-NaN iHS values
    below_threshold_mask = chrom_p_values < threshold
    above_threshold_mask = chrom_p_values >= threshold
    
    # Plot points below and above the threshold
    ax.scatter(adjusted_positions[below_threshold_mask], 
               chrom_p_values[below_threshold_mask], 
               color=chromosome_colours[chrom], alpha=0.1, s=10)
    ax.scatter(adjusted_positions[above_threshold_mask], 
               chrom_p_values[above_threshold_mask], 
               color=chromosome_colours[chrom], alpha=1.0, s=10)
    patch = mpatches.Patch(color=chromosome_colours[chrom], label=chrom)
    legend_patches.append(patch)

legend_patches.append(mpatches.Patch(color='black', label='Significance Threshold'))
ax.legend(handles=legend_patches, title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')

ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('$- log10 pvalue |iHS|$')
ax.axhline(y=threshold, color='black', linestyle='--')
plt.tight_layout()
plt.show()

print("iHS p-values plotted with all chromosomes")

## %% list all positions with iHS value over certain threshold

mask_pvalues_above_threshold = neg_log_pvals >= threshold
significant_ihs_positions = pos_withoutnan[mask_pvalues_above_threshold]
significant_ihs_chromosomes = chrom_withoutnan[mask_pvalues_above_threshold]
significant_ihs_values = neg_log_pvals[mask_pvalues_above_threshold]

print("iHS p-values above threshold identified")

## Save positions and corresponding iHS values above the threshold to a text file

with open(f"PuertoRico_significant_iHS_threshold_{threshold}.txt", "w") as file:
    for chrom, position, ihs_value in zip(significant_ihs_chromosomes, significant_ihs_positions, significant_ihs_values):
        file.write(f"{chrom}\t{position}\t{ihs_value}\n")

# %% bring in the gff file to understand where each of these variants is

print("Using GFF file to bring in annotations for these positions")

# Parameters
input_file_name = f"PuertoRico_significant_iHS_threshold_{threshold}.txt"
output_file_name = f"PuertoRico_significant_iHS_threshold_{threshold}_GFF_annotated.txt"
gff_file = '/mnt/storage12/emma/Reference_files/genes_only.gff'

# Function to find and format the GFF line(s) that overlap a given position
def find_overlapping_gff_lines(chromosome, position, gff_file):
    overlapping_lines = []
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#') or line.strip() == "":
                continue  # Skip header and empty lines
            parts = line.split('\t')
            if parts[0] == chromosome and int(parts[3]) <= position <= int(parts[4]):
                formatted_line = ",".join(parts).strip()
                overlapping_lines.append(formatted_line)
    return overlapping_lines

# Open the output file to write the annotated positions
with open(output_file_name, "w") as outfile:
    # Write the header line
    outfile.write("Chromosome\tPosition\tiHS Value\tGff_Annotation\n")

    # Open the file containing significant iHS positions to read
    with open(input_file_name, "r") as infile:
        for line in infile:
            parts = line.strip().split("\t")
            chromosome, position, ihs_value = parts[0], int(parts[1]), parts[2]
            
            # Find overlapping GFF lines for the position
            overlapping_gff_lines = find_overlapping_gff_lines(chromosome, position, gff_file)
            
            # Join all overlapping GFF lines into a single string
            gff_annotation = "; ".join(overlapping_gff_lines)
            
            # Write to the output file
            outfile.write(f"{chromosome}\t{position}\t{ihs_value}\t{gff_annotation}\n")

print(f"iHS significant values identified and GFF annotations written here: {output_file_name}")

#%%






























#################################### OLD ###########################################




# %% Compute empirical p-values, then can log transform these and plot.
# the standardised ihs values are saved in ihs_gb_std[0]

# write function
def calculate_empirical_p_value(sorted_values, observed_value):
    """
    Calculate the empirical p-value for an observed value in a sorted list of values.

    Parameters:
    - sorted_values: A list of values, sorted in ascending order.
    - observed_value: The observed value for which to calculate the p-value.

    Returns:
    - The empirical p-value.
    """
    # Count how many values are as extreme or more extreme than the observed value
    # For a one-tailed test, this would be how many values are greater than or equal to the observed value
    extreme_values_count = sum(value >= observed_value for value in sorted_values)

    # Calculate the empirical p-value
    p_value = extreme_values_count / len(sorted_values)

    return p_value

# %% Filter out NaN values from the original iHS scores and keep their original positions
vals = ihs_gb_std[0]
mask_non_nan = ~np.isnan(vals) #remove the ihs_gb_std nan values
vals_withoutnan = vals[mask_non_nan] # save these values as vals_withoutnan
pos_withoutnan = pos_gb_seg[mask_non_nan] # use the same mask to get the corresponding positions of vals_withoutnan
chrom_withoutnan = chrom_gb_seg[mask_non_nan] # use the same mask to get the corresponding chromosomes of vals_withoutnan
print("Filtered out NaN values")

# %% Calculate p-values for non-NaN iHS scores
sorted_ihs = np.sort(vals_withoutnan) #sort these by putting them in ascending order, ready for the calculate_empirical_p_value function
print("Put values into ascending order")

print("Starting to calculate p-values of iHS")

# %% Multithread the p-value calculation because otherwise it is mega slow
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def calculate_pval_concurrently(vals_withoutnan, sorted_ihs):
    with ThreadPoolExecutor(max_workers=25) as executor:
        # Submit all tasks and keep track of futures
        futures = {executor.submit(calculate_empirical_p_value, sorted_ihs, ihs): ihs for ihs in vals_withoutnan}
        # Prepare to collect p-values
        pvals = []
        # Iterate over futures as they complete
        for future in tqdm(as_completed(futures), total=len(futures), desc='Computing p-values'):
            # Collect result from each future
            pvals.append(future.result())
    return pvals

# Use the function to calculate p-values
pvals = calculate_pval_concurrently(vals_withoutnan, sorted_ihs)


print("P-values computed")

# %% save pvals for next time
# output p-values retain their original order corresponding to vals_withoutnan
pvals_df = pd.DataFrame(pvals, columns=['p_value'])
pvals_df.to_csv('pvals_paralleled_screen.csv', index=False)

print("Saved p-values to dataframe")

# %% load pvals back in
#pvals_loaded_df = pd.read_csv('pvals_paralleled_screen.csv')
#pvals_loaded = pvals_loaded_df['p_value'].values

# %% save pvals_loaded as pvals so don't have to change the rest of the script
#pvals = pvals_loaded

# %% Adjust p-values to avoid log10(0)
adjusted_pvals = np.maximum(pvals, 1e-10)

print("Adjusted p-values to avoid log10(0)")

# %% Taking the log10
log_pvals = np.log10(adjusted_pvals)

print("Computed log10 of p-values")

print("Plotting log10 of p-values")

# %% Define colors for each chromosome (for illustration)
chromosome_colours = {
    '2L': '#3d348b', '2R': '#f18701', '3L': '#f7b801', '3R': '#7678ed', 'anop_mito': '#f35b04', 'anop_X': '#119DA4'
}

# %% Plotting
plt.figure(figsize=(10, 6))

# Iterate through each chromosome to plot its variants with colors
for chrom in sorted(set(chrom_withoutnan)):
    mask = chrom_withoutnan == chrom
    plt.scatter(pos_withoutnan[mask], log_pvals[mask], color=chromosome_colours[chrom], alpha=0.6, label=chrom)

plt.xlabel('Genomic Position (bp)')
plt.ylabel('-log10(P-value)')
plt.title('Manhattan Plot of iHS Scores Colored by Chromosome')
plt.grid(True)

# Optional: Create a legend if you want to identify chromosomes by color
plt.legend(title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.show()





















# Hashing out below to run script to here in screen
#
#
## %% list all positions with iHS value over certain threshold
#threshold = 5
#mask_above_threshold = ihs_gb_std[0] >= threshold
#positions_above_threshold = pos_gb_seg[mask_above_threshold]
#ihs_values_above_threshold = ihs_gb_std[0][mask_above_threshold]
#
#chromosomes = callset['variants/CHROM'][:]
#chromosomes_above_threshold = chromosomes[gb_seg_variants][mask_above_threshold]
#
## Save positions and corresponding iHS values above the threshold to a text file
## Adjusted file writing to include chromosome, position, and iHS value
#with open(f"gb_ihs_positions_above_threshold_{threshold}.txt", "w") as file:
#    for chrom, position, ihs_value in zip(chromosomes_above_threshold, positions_above_threshold, ihs_values_above_threshold):
#        file.write(f"{chrom}\t{position}\t{ihs_value}\n")
#
## %% bring in the gff file to understand where each of these variants is
#
#gff_file = '/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3'
#db_file = 'annotations.db'
## create database from gff file (already created annotations.db now)
## db = gffutils.create_db(gff_file, dbfn=db_file, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
#
## Querythe database to identify genes associated with the significant positions
## Initialize an empty list to store the positions
#positions_above_threshold = []
#
## Open and read the file
#with open("gb_ihs_positions_above_threshold_5.txt", "r") as file:
#    for line in file:
#        parts = line.strip().split("\t")
#        # Assuming the format is: chromosome, position, iHS value
#        chromosome, position = parts[0], int(parts[1])
#        positions_above_threshold.append((chromosome, position))
#
## Connect to the database
#db = gffutils.FeatureDB(db_file)
#
## Open a new file to write the annotated positions
#with open("gb_ihs_positions_above_threshold_5_gff_annotated.txt", "w") as outfile:
#    # Write the header line
#    outfile.write("Chromosome\tPosition\tiHS Value\tOverlapping Genes\n")
#    
#    # Read the original file again to keep the iHS values this time
#    with open("gb_ihs_positions_above_threshold_5.txt", "r") as infile:
#        for line in infile:
#            parts = line.strip().split("\t")
#            chromosome, position, ihs_value = parts[0], parts[1], parts[2]
#            overlapping_genes = db.region(seqid=chromosome, start=int(position), end=int(position))
#            
#            # Collect all overlapping gene IDs
#            gene_ids = [gene.id for gene in overlapping_genes]
#            
#            # Write to outfile
#            outfile.write(f"{chromosome}\t{position}\t{ihs_value}\t{','.join(gene_ids)}\n")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
## %% ########### Cross-population extended haplotype homozygosity (XPEHH) ###########
#
### Compute the unstandardized cross-population extended haplotype homozygosity score (XPEHH) for each variant.
### allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
## create h1 and h2, selecting all variants instead of segregating variants only, which is what we did in iHS
#
#
## %%
### VCF is phased so we can convert genotype arrays made earlier to haplotype array
### Create arrays needed for Cameroon samples
#sample_ids = callset['samples'][:]
## Get sample identifiers for Cameroon samples from df_samples
#cam_sample_ids = df_samples[df_samples['country'] == 'Cameroon']['sample'].values
## Find indices of these samples in the genotype array
#cam_indices = np.array([np.where(sample_ids == id)[0][0] for id in cam_sample_ids if id in sample_ids])
## Verify the indices are within the correct range
#print("Max index:", cam_indices.max(), "Sample array size:", len(sample_ids))
## Select genotypes for Cameroon samples using the indices
#gt_cam_samples = gt.take(cam_indices, axis=1)
#
## %% select variants that are segregating within cam_samples as only these will be informative
### also some selection tests don't support multiallelic variants, so just keep biallelics
#
#ac_cam = gt_cam_samples.count_alleles(max_allele=8).compute()
#cam_seg_variants = ac_cam.is_segregating() & ac_cam.is_biallelic_01()
#ac_cam_seg = ac_cam.compress(cam_seg_variants, axis=0)
#gt_cam_seg = gt_cam_samples.compress(cam_seg_variants, axis = 0)
#gt_cam_seg
#
## %% this is from a phased VCF so we can convert this genotype array to haplotype array
#
#h_cam_seg = gt_cam_seg.to_haplotypes().compute()
#h_cam_seg
#
## %% we need variant positions
#pos = callset['variants/POS'][:]
#pos_cam_seg = pos.compress(cam_seg_variants, axis=0)
#pos_cam_seg
#
## Let's check if there any of genomic positions with multiple variants.
#count_multiple_variants = np.count_nonzero(np.diff(pos_cam_seg == 0))
#
#if count_multiple_variants == 0:
#    print("No cases where there are multiple variants at the same genomic position, script will continue")
#else:
#    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
#    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line
#
## %% Continue with xp-ehh 
#
#
#
#
#h_sus = gt_sus_samples.to_haplotypes().compute()
#h_sus
#
#h_res = gt_res_samples.to_haplotypes().compute()
#h_res
#
#ac_gt = gt.count_alleles(max_allele=8).compute()
#
## %%
## get variant positions
#
#pos = callset['variants/POS'][:]
#pos
#
## %% look at shapes of arrays
#print("h_sus shape:", h_sus.shape)
#print("h_res shape:", h_res.shape)
#print("pos shape", pos.shape)
#
## %% compute xpehh
## xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=20000, is_accessible=None, use_threads=True)
#
## xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, include_edges=True, use_threads=True)
#xpehh_raw = allel.xpehh(h_sus, h_res, pos, use_threads=True)
#xpehh_raw
#
## %% look for where the biggest signal is
#xpehh_hit_max = np.nanargmax(xpehh_raw)
#xpehh_hit_max
#
## %% genomic position of top hit
#pos[xpehh_hit_max]
#
## %%
#%matplotlib inline
#import matplotlib.pyplot as plt
#
## %%
#
#fig, ax = plt.subplots()
#ax.hist(xpehh_raw[~np.isnan(xpehh_raw)], bins=20)
#ax.set_xlabel('Raw XP-EHH')
#ax.set_ylabel('Frequency (no. variants)');
#
## %% Standardize XP-EHH - do not think that we really need to do this
#
## xpehh_std = allel.standardize_by_allele_count(xpehh_raw, ac_gt[:, 1])
#
## %% 
##fig, ax = plt.subplots()
##ax.hist(xpehh_std[0][~np.isnan(xpehh_std[0])], bins=20)
##ax.set_xlabel('Standardized XP-EHH')
##ax.set_ylabel('Frequency (no. variants)');
#
## %% note that iHS has been calculated with unpolarized data, so only the magnitude of iHS
## is informative, not the sign.
## xpehh_std
#
## %% look at shapes
#
#print("pos shape:", pos.shape)
#print("xpehh_raw shape:", xpehh_raw.shape)
#
#min_pos = pos.min()
#max_pos = pos.max()
#
#print("Minimum Genomic Position:", min_pos)
#print("Maximum Genomic Position:", max_pos)
#
#
## %% plot on manhattan plot
#
#fig, ax = plt.subplots(figsize=(10, 3))
#ax.plot(pos, np.abs(xpehh_raw), linestyle=' ', marker='o', mfc='none', mew=3, mec='k', label='$|XP-EHH|$')
#ax.axhline(y=4, color='red', linestyle='--')
#ax.set_xlabel('Genomic position (bp)')
#ax.set_ylabel('$|XP-EHH|$')
#ax.set_ylim(0, 8)
#ax.set_xlim(0,61542681)
#ax.legend()
#
## %% list all positions with iHS value over certain threshold
#
#xpehh_positions_above_threshold_4 = pos[xpehh_raw >= 4]
#
## Save positions_above_threshold to a text file
#with open("xpehh_positions_above_threshold_4", "w") as file:
#    for position in xpehh_positions_above_threshold_4:
#        file.write(str(position) + "\n")
#
#
#
## %% ################################ TAJIMA'S D #####################################
#
## allel.moving_delta_tajima_d(ac1, ac2, size, start=0, stop=None, step=None
## compute the difference in Tajima's D between two populations in moving windows
#
## create allele counts arrays ac1 and ac2
## genotype arrays were made earlier in this script: gt_sus_samples and gt_res_samples
#
#ac_sus_samples = gt_sus_samples.count_alleles()
#ac_sus_samples
#
#ac_res_samples = gt_res_samples.count_alleles()
#ac_res_samples
#