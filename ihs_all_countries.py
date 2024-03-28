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
df_samples=pd.read_csv('all_sample_country_metadata_underscore.csv',sep=',',usecols=['ID','Country','Region','Subregion'])
df_samples.head()
df_samples.groupby(by=['Country']).count

# %%
## Looping through all populations to calculate ihs for each
sample_ids = callset['samples'][:]


sample_indices = {}  # Dictionary to store sample indices

for country in np.unique(df_samples.Country):
    country_samples = df_samples[df_samples['Country'] == country]['ID'].values
    sample_indices[f"{country}_samples"] = country_samples
    globals()[f"{country}_indicies"] = np.array([np.where(sample_ids == id)[0][0] for id in country_samples if id in sample_ids]) 
    print("Max index:", globals()[f"{country}_indicies"].max(), "Sample array size:", len(sample_ids))
    gt = allel.GenotypeDaskArray(callset['calldata/GT'])
    globals()[f"gt_{country}_samples"] = gt.take(globals()[f"{country}_indicies"], axis=1)

    ## select variants that are segregating within gb_samples as only these will be informative
    ## also some selection tests don't support multiallelic variants, so just keep biallelics
    ## for this pipeline the VCF is already filtered so should be no biallelic SNPs anyway
    globals()[f"ac_{country}"] = globals()[f"gt_{country}_samples"].count_alleles(max_allele=8).compute()
    globals()[f"{country}_seg_variants"] = globals()[f"ac_{country}"].is_segregating() & globals()[f"ac_{country}"].is_biallelic_01()
    globals()[f"ac_{country}_seg"] = globals()[f"ac_{country}"].compress(globals()[f"{country}_seg_variants"], axis=0)
    globals()[f"gt_{country}_seg"] = globals()[f"gt_{country}_samples"].compress(globals()[f"{country}_seg_variants"], axis = 0)
    
    ## this is from a phased VCF so we can convert this genotype array to haplotype array
    globals()[f"h_{country}_seg"] = globals()[f"gt_{country}_seg"].to_haplotypes().compute()


    # we need variant positions
    pos = callset['variants/POS'][:]
    globals()[f"pos_{country}_seg"] = pos.compress(globals()[f"{country}_seg_variants"], axis=0)
    

    # some variants in 1000 genomes project have multiple variants at the same genomic position, 
    # which causes problems for some selection tests in scikit-allel. 
    # Let's check if there any of these.
    count_multiple_variants = np.count_nonzero(np.diff(globals()[f"pos_{country}_seg"] == 0))

    if count_multiple_variants == 0:
        print("No cases where there are multiple variants at the same genomic position, script will continue")
    else:
        print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    #sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

    # compute raw iHS
    globals()[f"ihs_{country}_raw"] = allel.ihs(globals()[f"h_{country}_seg"], globals()[f"pos_{country}_seg"], use_threads=True, include_edges=True)
    print("Raw iHS computed")



#%%
# matplotlib inline
import matplotlib.pyplot as plt
from datetime import datetime

# %% view raw iHS values as a histogram
# ~np.isnan(ihs_gb_std[0]) is used to filter out NaN values
for country in np.unique(df_samples.Country):
    fig, ax = plt.subplots()
    ax.hist(globals()[f"ihs_{country}_raw"][~np.isnan(globals()[f"ihs_{country}_raw"])], bins=20)
    ax.set_xlabel('Raw IHS')
    ax.set_ylabel('Frequency (no. variants)');

# %% Standardize iHS
for country in np.unique(df_samples.Country):
    globals()[f"ihs_{country}_std"] = allel.standardize_by_allele_count(globals()[f"ihs_{country}_raw"], globals()[f"ac_{country}_seg"][:, 1])
    print("Standardized iHS computed")

# %% 

# Here we deviate from the Jupyter notebook and use ihs_res_std[0]
# ~np.isnan(ihs_gb_std[0]) is used to filter out NaN values
for country in np.unique(df_samples.Country):
    fig, ax = plt.subplots()
    ax.hist(globals()[f"ihs_{country}_std"][0][~np.isnan(globals()[f"ihs_{country}_std"][0])], bins=20)
    ax.set_xlabel('Standardised IHS')
    ax.set_ylabel('Frequency (no. variants)');

    # Save the figure as a file (e.g., PNG) in the current working directory
    filename = f'standardised_ihs_{country}_histogram.png'
    #plt.savefig(filename)

    # show the plot (optional, could # this out)
    plt.show()

# %% note that iHS has been calculated with unpolarized data using ~np.isnan, so only the magnitude of iHS
# is informative, not the sign.

# plot over the genome
# np.abs is converting all iHS vales to their absoltue values before plotting. This means that if ihs_gb_std[0]
# contains any negative valyes, those values will be made positive. It is plotting the magnitude of iHS without considering the
# direction of selection, which the sign of iHS could indicate
for country in np.unique(df_samples.Country):
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(globals()[f"pos_{country}_seg"], np.abs(globals()[f"ihs_{country}_std"][0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
    ax.set_xlabel('Genomic position (bp) chromosome agnostic')
    ax.set_ylabel('$|IHS|$')
    ax.set_ylim(-2, 9);

    # Save the figure as a file (e.g., PNG) in the current working directory
    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    filename = f'ihs_manhattan_{timestamp}_{country}.png'
    #plt.savefig(filename)

    print("iHS plotted")

# %% find the index of the variant with the highest iHS value
for country in np.unique(df_samples.Country):
    globals()[f"idx_hit_max_{country}"] = np.nanargmax(globals()[f"ihs_{country}_std"][0])
    pos_max = globals()[f"pos_{country}_seg"][globals()[f"idx_hit_max_{country}"]]
    print(globals()[f"idx_hit_max_{country}"])
    print(f"highest value is at position {pos_max}")


# %% Visualise EHH decay around top hit with a Voight plot
# pull out haplotypes for region around top hit
for country in np.unique(df_samples.Country):
    flank_size = 2000
    globals()[f"h_hit_{country}"] = globals()[f"h_{country}_seg"][globals()[f"idx_hit_max_{country}"] - flank_size:globals()[f"idx_hit_max_{country}"] + flank_size + 1]
    print(globals()[f"h_hit_{country}"])

# %%

for country in np.unique(df_samples.Country):
    fig = allel.fig_voight_painting(globals()[f"h_hit_{country}"][:, globals()[f"h_{country}_seg"][globals()[f"idx_hit_max_{country}"]] == 0], index=flank_size, height_factor=0.02)
    fig.suptitle('Reference allele', y=1)
    print(country);

#%%
for country in np.unique(df_samples.Country):
    fig2 = allel.fig_voight_painting(globals()[f"h_hit_{country}"][:, globals()[f"h_{country}_seg"][globals()[f"idx_hit_max_{country}"]] == 1], index=flank_size, height_factor=0.02)
    fig2.suptitle('Alternate allele', y=1);


# %% Plot iHS
# `pos`, `chromosomes`, and `ihs_gb_std[0]` arrays are already defined and aligned
# Define the threshold
import matplotlib.patches as mpatches

# %% Plot iHS manhattan plot with all chromosomes
for country in np.unique(df_samples.Country):
    # length of gb_seg_variants = 6785669, it is a numpy.ndarray
    chromosomes = callset['variants/CHROM'][:]
    globals()[f"chrom_{country}_seg"] = chromosomes.compress(globals()[f"{country}_seg_variants"], axis = 0)
    
    pos = callset['variants/POS'][:]
    globals()[f"pos_{country}_seg"] = pos.compress(globals()[f"{country}_seg_variants"], axis=0)

    # define chromosome lengths and colours 
    chromosome_lengths = {
        'MT': 16790,
        '1': 310827022,
        '2': 474425716,
        '3': 409777670
    }

    cumulative_lengths = {}
    cumulative_length = 0
    for chrom, length in chromosome_lengths.items():
        cumulative_lengths[chrom] = cumulative_length
        cumulative_length += length

    threshold = 5

    # Set up the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define colors for each chromosome (for illustration)
    chromosome_colours = {
        '1': '#3d348b', '2': '#f18701', '3': '#f7b801', 'MT': '#f35b04'
    }

    # Create a list to hold the legend patches
    legend_patches=[]

    # Filter out 'Y_unplaced' or any other chromosomes not in your chromosome_lengths dictionary
    filtered_chroms = [chrom for chrom in sorted(set(globals()[f"chrom_{country}_seg"])) if chrom in chromosome_lengths]

    # Iterate through each chromosome to plot its variants
    for chrom in filtered_chroms:
        # Create mask for the current chromosome
        mask = globals()[f"chrom_{country}_seg"] == chrom
    
        # Apply the chromosome mask and filter out NaN values simultaneously
        chrom_positions = globals()[f"pos_{country}_seg"][mask]
        chrom_ihs_values = globals()[f"ihs_{country}_std"][0][mask]
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
    plt.title(f'iHS {country}', size =12)
    ax.axhline(y=threshold, color='black', linestyle='--')
    ax.axhline(y=-5, color = 'black', linestyle = '--')
    plt.ylim(-6.5, 6.5)
    filename = f'iHS_value_200_iterations_{country}.png'
    #plt.savefig(filename)
    plt.tight_layout()
    plt.show()

    print("iHS plotted with all chromosomes")


########################################### REMOVE NAS ################################

# %% Remove NaN iHS values, and use the same mask to filter pos and chrom.

for country in np.unique(df_samples.Country):
    globals()[f"ihs_{country}_vals"] = globals()[f"ihs_{country}_std"][0] # the standardised ihs values are saved in ihs_gb_std[0]
    globals()[f"mask_non_nan_{country}"] = ~np.isnan(globals()[f"ihs_{country}_vals"]) #remove the ihs_gb_std nan values
    globals()[f"ihs_vals_withoutnan_{country}"] = globals()[f"ihs_{country}_vals"][globals()[f"mask_non_nan_{country}"]] # save these values as vals_withoutnan
    globals()[f"pos_withoutnan_{country}"] = globals()[f"pos_{country}_seg"][globals()[f"mask_non_nan_{country}"]] # use the same mask to get the corresponding positions of vals_withoutnan
    globals()[f"chrom_withoutnan_{country}"] = globals()[f"chrom_{country}_seg"][globals()[f"mask_non_nan_{country}"]] # use the same mask to get the corresponding chromosomes of vals_withoutnan
    print("Filtered out NaN values")

    #Sort these by putting them in ascending order, ready for the calculate_empirical_p_value function
    globals()[f"sorted_indices_{country}"] = np.argsort(globals()[f"ihs_vals_withoutnan_{country}"]) # np.argsort returns an array of teh corresponding indices
    globals()[f"sorted_ihs_{country}"] = np.sort(globals()[f"ihs_vals_withoutnan_{country}"]) 
    print("Put values into ascending order")

# %% Compute empirical p-values, then can log transform these and plot.
################################# CALCULATE P VALS ################################

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
for country in np.unique(df_samples.Country):
    print("Starting to calculate p-values of iHS")

    # Multithread the p-value calculation because otherwise it is mega slow
    import joblib
    import tqdm

    len_ihs = len(globals()[f"sorted_ihs_{country}"])
    pvals = []
    from joblib import Parallel, delayed

    parallel = Parallel(n_jobs=15, return_as='generator')
    globals()[f"pvals_{country}"] = [r for r in tqdm.tqdm(parallel(delayed(calculate_empirical_p_value)(val,globals()[f"sorted_ihs_{country}"],len_ihs) for val in globals()[f"sorted_ihs_{country}"]),total=len(globals()[f"sorted_ihs_{country}"]))]

# %% P-values are currently in the order of sorted iHS, need to reorder them back 
# based on the position of the original ihs value in ihs_vals_withoutnan array
for country in np.unique(df_samples.Country):
    
    globals()[f"reordered_pvals_{country}"] = np.empty(len(globals()[f"pvals_{country}"]))
    globals()[f"reordered_pvals_{country}"][globals()[f"sorted_indices_{country}"]] = globals()[f"pvals_{country}"] # reordered_pvals will now have the p-values in the original order of the ihs values


# %% Take log of p-values
for country in np.unique(df_samples.Country):
    globals()[f"neg_log_pvals_{country}"] = (-1*(np.log10(globals()[f"reordered_pvals_{country}"])))
    print("Computed log10 of p-values")

# %% Plotting adjusted p-values
################################### PLOT P VALS ############################################
all_data = []

for country in np.unique(df_samples.Country):
    print("Plotting neg log10 of p-values")

    # Define colors for each chromosome (for illustration)
    chromosome_colours = {
        '1': '#3d348b', '2': '#f18701', '3': '#f7b801', 'MT': '#7678ed'
    }
    chromosome_colours = {
    '1': '#A6C48A', '2': '#40476D', '3': '#51A3A3', 'MT': '#bc80bd'}
    # Plotting
    threshold = 4
    # Set up the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    # Create a list to hold the legend patches
    legend_patches=[]

    # Filter out 'Y_unplaced' or any other chromosomes not in your chromosome_lengths dictionary
    filtered_chroms = ['1', '2', '3', 'MT']

    # Iterate through each chromosome to plot its variants
    for chrom in filtered_chroms:
        # Create mask for the current chromosome
        mask = globals()[f"chrom_withoutnan_{country}"] == chrom
    
        # Apply the chromosome mask
        chrom_positions = globals()[f"pos_withoutnan_{country}"][mask]
        chrom_p_values = globals()[f"neg_log_pvals_{country}"][mask]
    
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
    plt.title(f'{country}')
    filename = f'iHS_pval_200_iterations_{country}.png'
    #plt.savefig(filename)
    plt.tight_layout()
    plt.show()

    print("iHS p-values plotted with all chromosomes")

    ## %% list all positions with iHS value over certain threshold

    mask_pvalues_above_threshold = globals()[f"neg_log_pvals_{country}"] >= threshold
    significant_ihs_positions = globals()[f"pos_withoutnan_{country}"][mask_pvalues_above_threshold]
    significant_ihs_chromosomes = globals()[f"chrom_withoutnan_{country}"][mask_pvalues_above_threshold]
    significant_ihs_values = globals()[f"neg_log_pvals_{country}"][mask_pvalues_above_threshold]
    country_array = np.full(len(significant_ihs_values), f"{country}")
    
    country_df = np.column_stack((country_array, significant_ihs_chromosomes, significant_ihs_positions, significant_ihs_values))
    country_df = pd.DataFrame(country_df, columns= ['country', 'chrom', 'pos', 'pval'])
    all_data.append(country_df)

    print("iHS p-values above threshold identified")


    ## Save positions and corresponding iHS values above the threshold to a text file

    with open(f"significant_{country}_iHS_threshold_{threshold}.txt", "w") as file:
        for chrom, position, ihs_value in zip(significant_ihs_chromosomes, significant_ihs_positions, significant_ihs_values):
            file.write(f"{chrom}\t{position}\t{ihs_value}\n")

#df = pd.DataFrame(all_data, columns= ['country', 'chrom', 'pos', 'pval'])
df = pd.concat(all_data, ignore_index = True)
df.to_csv("all_ihs_significant_4.csv", index= None)

# %% bring in the gff file to understand where each of these variants is

print("Using GFF file to bring in annotations for these positions")

# Parameters
#input_file_name = f"significant_{country}_iHS_threshold_{threshold}.txt"
input_file_name = "all_ihs_significant_4.csv"
#output_file_name = f"significant_{country}_iHS_threshold_{threshold}_GFF_annotated.txt"
output_file_name = "all_ihs_significant_4_annotated.csv"
gff_file = '/mnt/storage12/emma/ihs/genes_numeric_chrom.gff'

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
        heade =next(infile)
        for line in infile:
            parts = line.strip().split(",")
            country, chromosome, position, ihs_value = parts[0], parts[1], int(parts[2]), parts[3]
            
            # Find overlapping GFF lines for the position
            overlapping_gff_lines = find_overlapping_gff_lines(chromosome, position, gff_file)
            
            # Join all overlapping GFF lines into a single string
            gff_annotation = "; ".join(overlapping_gff_lines)
            
            # Write to the output file
            outfile.write(f"{country}\t{chromosome}\t{position}\t{ihs_value}\t{gff_annotation}\n")

print(f"iHS significant values identified and GFF annotations written here: {output_file_name}")


















