# ########### Cross-population extended haplotype homozygosity (XPEHH) ###########

## Compute the unstandardized cross-population extended haplotype homozygosity score (XPEHH) for each variant.
## allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
# create h1 and h2, selecting all variants instead of segregating variants only, which is what we did in iHS

# %%
import os
import numpy as np
import allel
import zarr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gffutils

# %% set wd
os.chdir('/mnt/storage12/emma/selection/xpehh/')
os.getcwd()

# %%
# convert phased, filtered, VCF file to zarr file
# already converted to zarr
#allel.vcf_to_zarr('all_aedes_norm_filt_miss0.5_mac3_minQ30.vcf.gz.recode.chrom.lmiss.filt.vcf.gz.recode.main_chrom_only_renamed.phased.vcf.gz', 'all_aedes_lmiss_phased_rename_chrom.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('all_aedes_lmiss_phased.zarr', mode='r')
#callset.tree(expand=True)

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
## VCF is phased so we can convert genotype arrays made earlier to haplotype array

sample_ids = callset['samples'][:]


#%% make array for each country

sample_indices = {}  # Dictionary to store sample indices

for country in np.unique(df_samples.Country):
    country_samples = df_samples[df_samples['Country'] == country]['ID'].values
    sample_indices[f"{country}_samples"] = country_samples
    #globals()[f"gt_{country}"] = filtered_gt.take(country_samples, axis=1)
    #globals()[f"h_{country}_seg"] = globals()[f"gt_{country}"].to_haplotypes().compute()
    globals()[f"{country}_indicies"] = np.array([np.where(sample_ids == id)[0][0] for id in country_samples if id in sample_ids]) 
    print("Max index:", globals()[f"{country}_indicies"].max(), "Sample array size:", len(sample_ids))
    globals()[f"gt_{country}_samples"] = gt.take(globals()[f"{country}_indicies"], axis=1)

    globals()[f"h_array_{country}"] = globals()[f"gt_{country}_samples"].to_haplotypes().compute()




# %% we need variant positions
pos = callset['variants/POS'][:]
chrom = callset['variants/CHROM'][:]


#%%
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import itertools

# %% compute xpehh
for i in range(0,56):
    country1 = list(itertools.permutations(np.unique(df_samples.Country),2))[i][0]
    country2 = list(itertools.permutations(np.unique(df_samples.Country),2))[i][1]
    globals()[f"xpehh_raw_{country1}_{country2}"] = allel.xpehh(globals()[f"h_array_{country1}"], globals()[f"h_array_{country2}"], pos, use_threads=True)



# %% look for where the biggest signal is
for i in range(0,56):
    country1 = list(itertools.permutations(np.unique(df_samples.Country),2))[i][0]
    country2 = list(itertools.permutations(np.unique(df_samples.Country),2))[i][1]

    #fig, ax = plt.subplots()
    #ax.hist(globals()[f"xpehh_raw_{country1}_{country2}"][~np.isnan(globals()[f"xpehh_raw_{country1}_{country2}"])], bins=20)
    #ax.set_xlabel('Raw XP-EHH')
    #ax.set_ylabel('Frequency (no. variants)');

    allele_counts_array = gt.count_alleles(max_allele=3).compute()
    globals()[f"xpehh_std_{country1}_{country2}"] = allel.standardize_by_allele_count(globals()[f"xpehh_raw_{country1}_{country2}"], allele_counts_array[:, 1])

    #fig, ax = plt.subplots()
    #ax.hist(globals()[f"xpehh_std_{country1}_{country2}"][0][~np.isnan(globals()[f"xpehh_std_{country1}_{country2}"][0])], bins=20)
    #ax.set_xlabel('Raw XP-EHH')
    #ax.set_ylabel('Frequency (no. variants)');

  # define chromosome lengths and colours 
    chromosome_lengths = {
    '035159': 16790,
    '035107': 310827022,
    '035108': 474425716,
    '035109': 409777670 }

    #Calculate cumulative offsets for each chromosome
    cumulative_lengths = {}
    cumulative_length = 0
    for chrom, length in chromosome_lengths.items():
        cumulative_lengths[chrom] = cumulative_length
        cumulative_length += length


    # Set threshold
    upper_threshold = 5
    lower_threshold = -5

    # Set up the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Ensure that pos, chrom, and xpehh_std are all numpy arrays to support advanced indexing
    pos = np.array(callset['variants/POS'][:])
    chrom = np.array(callset['variants/CHROM'][:])
    xpehh_standardised_values = np.array(globals()[f"xpehh_std_{country1}_{country2}"][0])

 
    # Define colors for each chromosome (for illustration)
    #chromosome_colours = {
    #'035107': '#3d348b', '035108': '#f18701', '035109': '#f7b801', '035159': '#f35b04'}
    chromosome_colours = {
    '035107': '#ffaa5c', '035108': '#F45739', '035109': '#40476D', '035159': '#51A3A3'}
    #DB7376
   #F45739 
    # Set up the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Create a list to hold the legend patches
    legend_patches = []

    # Filtered chromosomes list, assuming cumulative_lengths are defined for these
    #filtered_chroms = ['2L', '2R', '3L', '3R', 'anop_X', 'anop_mito']
    filtered_chroms = ['035107', '035108', '035109', '035159']

    # Iterate through each chromosome to plot its variants
    for unique_chrom in filtered_chroms:
        chrom_mask = chrom == unique_chrom
    
        chrom_positions = pos[chrom_mask]
        chrom_xpehh_values = xpehh_standardised_values[chrom_mask]
    
        non_nan_mask = ~np.isnan(chrom_xpehh_values)
        chrom_positions_no_nan = chrom_positions[non_nan_mask]
        chrom_xpehh_values_no_nan = chrom_xpehh_values[non_nan_mask]
    
        adjusted_positions = chrom_positions_no_nan + cumulative_lengths.get(unique_chrom, 0)

        # Conditions for plotting
        solid_mask = (chrom_xpehh_values_no_nan >= upper_threshold) | (chrom_xpehh_values_no_nan <= lower_threshold)
        faded_mask = ~solid_mask
    
        # Plot solid points for values above 5 or below -5
        ax.scatter(adjusted_positions[solid_mask], 
               chrom_xpehh_values_no_nan[solid_mask], 
               color=chromosome_colours[unique_chrom], alpha=1.0, s=10)
    
        # Plot faded points for other values
        ax.scatter(adjusted_positions[faded_mask], 
               chrom_xpehh_values_no_nan[faded_mask], 
               color=chromosome_colours[unique_chrom], alpha=0.1, s=10)
    
        # Add patch for the legend
        patch = mpatches.Patch(color=chromosome_colours[unique_chrom], label=unique_chrom)
        legend_patches.append(patch)

    # Add significance threshold lines and legend
    ax.axhline(y=upper_threshold, color='black', linestyle='--', label='{country1} Threshold')
    ax.axhline(y=lower_threshold, color='black', linestyle='--', label='{country2} Threshold')
    ax.legend(handles=legend_patches, title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Set labels
    ax.set_xlabel('Genomic Position (bp)')
    ax.set_ylabel('XP-EHH')
    plt.tight_layout()
    plt.title(f'Plot {country1} against {country2}')
    filename = f'plot_{country1}_{country2}.png'
    plt.savefig(filename)
    plt.show()
    plt.clf()

    globals()[f"{country1}_threshold_mask"] = xpehh_standardised_values >= upper_threshold
    globals()[f"{country2}_threshold_mask"] = xpehh_standardised_values >= lower_threshold

    country1_significant_chrom = chrom[globals()[f"{country1}_threshold_mask"]]
    country1_significant_pos = pos[globals()[f"{country1}_threshold_mask"]]
    country1_significant_xpehh = xpehh_standardised_values[globals()[f"{country1}_threshold_mask"]]
    country1_array = np.full(len(country1_significant_xpehh), country1)
    comparison1_array = np.full(len(country1_significant_xpehh), country2) 

    country2_significant_chrom = chrom[globals()[f"{country2}_threshold_mask"]]
    country2_significant_pos = pos[globals()[f"{country2}_threshold_mask"]]
    country2_significant_xpehh = xpehh_standardised_values[globals()[f"{country2}_threshold_mask"]]
    country2_array = np.full(len(country1_significant_xpehh), country2)
    comparison2_array = np.full(len(country2_significant_xpehh), country1)  

    # Combine the filtered data into a structured array
    country1_significant_xpehh_data = np.column_stack((country1_significant_chrom, country1_significant_pos, country1_significant_xpehh, country1_array, comparison1_array))
    country2_significant_xpehh_data = np.column_stack((country2_significant_chrom, country2_significant_pos, country2_significant_xpehh, country2_array, comparison2_array))

    filename1 = f'data_{country1}_sig_against_{country2}'
    filename2 = f'data_{country2}_sig_against_{country1}' 

    # Convert the structured array to a pandas DataFrame for easier handling
    df_significant_1_xpehh = pd.DataFrame(country1_significant_xpehh_data, columns=['Chromosome', 'Position', 'XPEHH', 'Country1', 'Country2'])
    df_significant_2_xpehh = pd.DataFrame(country2_significant_xpehh_data, columns = ['Chromosome', 'Position', 'XPEHH', 'Country1', 'Country2'])


    # Save to csv
    df_significant_1_xpehh.to_csv(f'df_significant_xpehh_threshold_{country1}_against{country2}.csv', index=False)
    df_significant_2_xpehh.to_csv(f'df_significant_xpehh_threshold_{country2}_against{country1}.csv', index=False)

    #all_data_array = []

    #all_data_array.append(df_significant_1_xpehh)

################################################
################################################
#%%
#bind dataframes




################################################
#%%
xpehh_hit_max = np.nanargmax(xpehh_raw)
xpehh_hit_max

# %% genomic position of top hit
pos[xpehh_hit_max]



# %% list all positions with xpehh value over or below a certain threshold
for i in range(0,56):
    country1 = list(itertools.permutations(np.unique(df_samples.Country),2))[i][0]
    country2 = list(itertools.permutations(np.unique(df_samples.Country),2))[i][1]

    globals()[f"{countrbijagos_threshold_mask = xpehh_standardised_values >= bijagos_threshold
    cameroon_threshold_mask = xpehh_standardised_values <= cameroon_threshold

# %% Apply the mask to filter the data
bij_significant_chrom = chrom[bijagos_threshold_mask]
bij_significant_pos = pos[bijagos_threshold_mask]
bij_significant_xpehh = xpehh_standardised_values[bijagos_threshold_mask]

cam_significant_chrom = chrom[cameroon_threshold_mask]
cam_significant_pos = pos[cameroon_threshold_mask]
cam_significant_xpehh = xpehh_standardised_values[cameroon_threshold_mask]

# %% Combine the filtered data into a structured array
bij_significant_xpehh_data = np.column_stack((bij_significant_chrom, bij_significant_pos, bij_significant_xpehh))
cam_significant_xpehh_data = np.column_stack((cam_significant_chrom, cam_significant_pos, cam_significant_xpehh))

# %% Convert the structured array to a pandas DataFrame for easier handling
df_significant_bij_xpehh = pd.DataFrame(bij_significant_xpehh_data, columns=['Chromosome', 'Position', 'XPEHH'])
df_significant_cam_xpehh = pd.DataFrame(cam_significant_xpehh_data, columns = ['Chromosome', 'Position', 'XPEHH'])

# %% Save to csv
df_significant_bij_xpehh.to_csv(f'df_significant_PR_xpehh_bijagos_threshold_{bijagos_threshold}.csv', index=False)
df_significant_cam_xpehh.to_csv(f'df_significant_Mexico_xpehh_cameroon_threshold_{cameroon_threshold}.csv', index=False)

# %% bring in the gff file to understand where each of these variants is

print("Using GFF file to bring in annotations for these positions")

## Annotate the Bijagos XPEHH file
# Parameters
input_file_name = f"df_significant_PR_xpehh_bijagos_threshold_{bijagos_threshold}.csv"
output_file_name = f"df_significant_PR_xpehh_bijagos_threshold_{bijagos_threshold}_annotated.csv"
gff_file = '/mnt/storage12/emma/Reference_files/genes.gff'

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
    outfile.write("Chromosome,Position,XPEHH,Gff_Annotation\n")

    # Open the file containing significant iHS positions to read
    with open(input_file_name, "r") as infile:
        next(infile) #skip header line
        for line in infile:
            parts = line.strip().split(",")
            chromosome, position, XPEHH = parts[0], int(parts[1]), parts[2]
            
            # Find overlapping GFF lines for the position
            overlapping_gff_lines = find_overlapping_gff_lines(chromosome, position, gff_file)
            
            # Join all overlapping GFF lines into a single string
            gff_annotation = "; ".join(overlapping_gff_lines)
            
            # Write to the output file
            outfile.write(f"{chromosome},{position},{XPEHH},{gff_annotation}\n")

print(f"XP-EHH significant values identified and GFF annotations written here: {output_file_name}")
# %%
## Annotate the Cameroon XPEHH file
# Parameters
input_file_name = f"df_significant_cam_xpehh_cameroon_threshold_{cameroon_threshold}.csv"
output_file_name = f"df_significant_cam_xpehh_cameroon_threshold_{cameroon_threshold}_annotated.csv"
gff_file = '/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3'

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
    outfile.write("Chromosome,Position,XPEHH,Gff_Annotation\n")

    # Open the file containing significant iHS positions to read
    with open(input_file_name, "r") as infile:
        next(infile) #skip header line
        for line in infile:
            parts = line.strip().split(",")
            chromosome, position, XPEHH = parts[0], int(parts[1]), parts[2]
            
            # Find overlapping GFF lines for the position
            overlapping_gff_lines = find_overlapping_gff_lines(chromosome, position, gff_file)
            
            # Join all overlapping GFF lines into a single string
            gff_annotation = "; ".join(overlapping_gff_lines)
            
            # Write to the output file
            outfile.write(f"{chromosome},{position},{XPEHH},{gff_annotation}\n")

print(f"XP-EHH significant values identified and GFF annotations written here: {output_file_name}")
# %%