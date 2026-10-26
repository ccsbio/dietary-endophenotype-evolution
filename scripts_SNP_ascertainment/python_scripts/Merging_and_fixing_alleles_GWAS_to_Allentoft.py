# Load libaries
import sys
import os
sys.path.append('/home/csic/gpo/olg/ccs/packages/vegeta_last_version')
import pandas as pd
import numpy as np
import merge_GWAS_phenos as vg
import ancient_vegeta as avg
import work_vegeta as wvg
import argparse
import time


# ------------------ PARSE ARGUMENTS ------------------ 
# Create a parser object
parser = argparse.ArgumentParser(description="This function checks out if variants from a GWAS dataframe are the same ones than within a filtered VCF from Allentodt aDNA dataset. Also, it updates de GWAS file by filtering and adding the IDs found in the VCF")
# Add arguments
parser.add_argument("--gwas_data_file", type=str, required=True, help="The path and name of a tsv file with gwas data for a set of snps and phenotypes (rows and columns respectively)")
parser.add_argument("--in_dir", type=str, required=True, help="The input path with VCF files to use in the analysis")
parser.add_argument("--vcf", type=str, required=True, help="The filename of the vcf holding aDNA info for a set of SNPs of interest which alleles will be compared with gwas data")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add when saving outputs without file extension (ex: neolitic_diet_90_indpMarkers)")
parser.add_argument("--outfile_clump_with_plink", default=False, action="store_true", help="If clump with plink is going to be performed, an .assoc file will be generated using SNP ID that must be found in VCF data rather than GWAS")


# Parse the arguments
args = parser.parse_args()

# Store arguments
file_gwas = args.gwas_data_file
vcf_path = args.in_dir
vcf_filename = args.vcf
out_path = args.out_dir
filename_out = args.info_outfile
bool_plink = args.outfile_clump_with_plink



# ------------------ RUN THE ANALYSIS ------------------ 
# get the start time
st = time.time()



# Read VCF
VCF_clust = vg.read_vcf(vcf_path, vcf_filename)

# Read merged GWAS data with vegeta package ;)
df_gwas = pd.read_csv(file_gwas, sep='\t')
df_gwas = wvg.fix_data_genomic_coordinates_to_int('chr', 'pos', df_gwas)
print(df_gwas.shape)

# Merge GWAS with a smaller version of the VCF
VCF_smaller = VCF_clust.iloc[:,0:5].copy()
VCF_smaller = VCF_smaller.rename(columns={'CHROM':'chr', 'POS':'pos'})

# Now merge
common_gwas_clust_allentoft = wvg.merge_and_fix_SNPs_between_two_datasets(df_gwas, VCF_smaller, chr_cols='chr', pos_cols='pos', ref_x_col='ref', alt_x_col='alt', ref_y_col='REF', alt_y_col='ALT')

# If any common genetic positions depicts different alleles, remove them
all_ok = all(common_gwas_clust_allentoft['Same_alleles'])
print(all_ok)
print("Number of same alleles:")
print(common_gwas_clust_allentoft['Same_alleles'].value_counts())

if all_ok == False:
    # Remove SNPs with unmatching alleles
    SNPs_to_remove = common_gwas_clust_allentoft.loc[common_gwas_clust_allentoft['Same_alleles'] == False,'ID'].tolist()
    common_gwas_clust_allentoft = common_gwas_clust_allentoft.loc[~common_gwas_clust_allentoft['ID'].isin(SNPs_to_remove),:]


# Finally, save the merged allele data that will match Allentoft VCFs
df_gwas_fixed = common_gwas_clust_allentoft.loc[:,list(df_gwas.columns)+['ID']].copy()
df_gwas_fixed = df_gwas_fixed.rename(columns={"ID": "ID_in_merged_dataset"})


# Check if .assoc file for plink --clump is needed
if bool_plink == True:
    assoc_plink = common_gwas_clust_allentoft.loc[:,['ID', 'p.value']].copy() # Use ID expected in the VCF, not in the GWAS
    assoc_plink = assoc_plink.rename(columns={"ID": "SNP", "p.value": "P"})
    print("Option enabled to create an input text file for --clump with plink. The file will be saved as %s" % (out_path+filename_out+".assoc"))
    assoc_plink.to_csv(out_path + filename_out + '.assoc', sep=' ', header=True, index=False)


# save
df_gwas_fixed.to_csv(out_path + filename_out + '.tsv.gz', sep='\t', header=True, index=False, compression='gzip')



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')