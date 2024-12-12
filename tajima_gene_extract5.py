import csv
import subprocess
import re

output_file = 'tabix_tajima_output_genes.txt'  # Output file where results will be stored

# Open the output file in write mode
with open(output_file, 'w') as out_f:
    # Open the input CSV file and read each line
    with open('/mnt/storage12/emma/nuc_div/all_sig_2_tajimasD_100kb.csv', 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            chrom = row['CHROM']
            bin_start = int(row['BIN_START'])
            bin_end = bin_start + 100000
            tajima_d = row['TajimaD']  # Get the TajimaD value
            country = row['FileFrom']  # Get the country value

            # Format CHROM to match the GTF file format
            formatted_chrom = f"NC_0{chrom}.1"
            
            # Construct the tabix command
            region = f"{formatted_chrom}:{bin_start}-{bin_end}"
            command = f"tabix /mnt/storage12/emma/Reference_files/genes_sorted.gff.gz {region}"
            
            # Run the tabix command and capture the output
            try:
                output = subprocess.check_output(command, shell=True, text=True)
                if output:
                    # Write the region header to the file
                    out_f.write(f"Genes in region {region} (TajimaD: {tajima_d}, FileFrom: {country}):\n")
                    
                    # Add the extracted gene name, TajimaD, and country values to each gene line
                    for line in output.strip().split('\n'):
                        fields = line.split('\t')
                        if len(fields) > 8:  # Ensure the 9th column exists
                            attributes = fields[8]
                            match = re.search(r'gene=([^;]+)', attributes)
                            if match:
                                gene_name = match.group(1)
                            else:
                                gene_name = "N/A"
                            out_f.write(f"{line}\t{gene_name}\t{tajima_d}\t{country}\n")
                    
                    out_f.write('\n')  # Ensure outputs are separated by newlines
            except subprocess.CalledProcessError as e:
                print(f"Error processing region {region}: {e}")
                out_f.write(f"Error processing region {region}: {e}\n")

print(f"Output written to {output_file}")
