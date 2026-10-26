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
parser = argparse.ArgumentParser(description="This function checks out if variants from a GWAS dataframe are the same ones than within a filtered VCF from 1000G. The same ones are kept and a fixed VCF for 1000G is created. Also, it updates de GWAS file by filtering and adding the IDs found in the VCF.")
# Add arguments
parser.add_argument("--gwas_data_file", type=str, required=True, help="The path and name of a tsv file with gwas data for a set of snps and phenotypes (rows and columns respectively)")
parser.add_argument("--in_dir", type=str, required=True, help="The input path with VCF files to use in the analysis")
parser.add_argument("--vcf", type=str, required=True, help="The filename of the vcf holding 1000G info for a set of SNPs of interest which alleles will be compared with gwas data")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add when saving outputs without file extension (ex: neolitic_diet_90_indpMarkers)")
#parser.add_argument("--outpath_updated_GWAS", type=str, required=True, help="Outpath to store updated GWAS matching previous data. No allele fixing, just filtering and adding 1000G ID")
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
#gwas_update = args.outpath_updated_GWAS



# ------------------ RUN THE ANALYSIS ------------------ 
# get the start time
st = time.time()



# Read VCF
VCF_clust = vg.read_vcf(vcf_path, vcf_filename)

# Read merged GWAS data with vegeta package ;)
df_gwas = pd.read_csv(file_gwas, sep='\t')
df_gwas = wvg.fix_data_genomic_coordinates_to_int('chr', 'pos', df_gwas)
df_gwas['Original_index_gwas'] = df_gwas.index
print(df_gwas.shape)

# Merge GWAS with a smaller version of the VCF
VCF_smaller = VCF_clust.iloc[:,0:9].copy()
VCF_smaller = VCF_smaller.rename(columns={'CHROM':'chr', 'POS':'pos'})

# Preprocess before doing merge
VCF_smaller['Original_index'] = VCF_clust.index
VCF_smaller, counter_indl = wvg.filter_INDELS(VCF_smaller, ['REF', 'ALT'])

# Now merge
common_gwas_1000g = wvg.merge_and_fix_SNPs_between_two_datasets(df_gwas, VCF_smaller, chr_cols='chr', pos_cols='pos', ref_x_col='ref', alt_x_col='alt', ref_y_col='REF', alt_y_col='ALT')

# If any common genetic positions depicts different alleles, remove them
all_ok = all(common_gwas_1000g['Same_alleles'])
print(all_ok)
print("Number of same alleles:")
print(common_gwas_1000g['Same_alleles'].value_counts())
lost_variants = df_gwas.shape[0] - common_gwas_1000g.shape[0]
perct = (lost_variants*100) / df_gwas.shape[0]
#print("Number of variants that did not overlap between aDNA and 1000G is %d which represents ~%f%" % (lost_variants, perct))
print(f"Number of variants that did not overlap between aDNA and 1000G is {lost_variants} which represents ~{perct:.2f}%")


if all_ok == False:
    # Remove SNPs with unmatching alleles
    SNPs_to_remove = common_gwas_1000g.loc[common_gwas_1000g['Same_alleles'] == False,'ID'].tolist()
    common_gwas_1000g = common_gwas_1000g.loc[~common_gwas_1000g['ID'].isin(SNPs_to_remove),:]


# Recover variants that are OK from the original VCF
VCF_filt = VCF_clust.loc[common_gwas_1000g['Original_index'].tolist(),]

# Save this version of the VCF
VCF_filt.to_csv(out_path + filename_out + '.vcf.gz', sep='\t', header=True, index=False, compression='gzip')


# Finally, save the merged allele data to 1000G (including the ID), but matching the allele data with aDNA
df_gwas_fixed = df_gwas.loc[common_gwas_1000g['Original_index_gwas'].tolist(),].copy()
df_gwas_fixed = df_gwas_fixed.assign(ID_1000G = common_gwas_1000g['ID'].values)
AF_EUR = wvg.parse_metadata_1000G_from_VCF_INFO(common_gwas_1000g, 'INFO', AA=False, AF=True, POP="EUR_AF")
df_gwas_fixed = df_gwas_fixed.assign(AF_EUR_1000G = AF_EUR.astype(float))
#df_gwas_fixed = df_gwas_fixed.rename(columns={"ID": "ID_in_merged_dataset"})


# Check if .assoc file for plink --clump is needed
if bool_plink == True:
    assoc_plink = common_gwas_1000g.loc[:,['ID', 'p.value']].copy() # Use ID expected in the VCF, not in the GWAS
    assoc_plink = assoc_plink.rename(columns={"ID": "SNP", "p.value": "P"})
    print("Option enabled to create an input text file for --clump with plink. The file will be saved as %s" % (out_path+filename_out+".assoc"))
    assoc_plink.to_csv(out_path + filename_out + '.assoc', sep=' ', header=True, index=False)


# save
gwas_filename = file_gwas.removesuffix(".tsv.gz")
df_gwas_fixed.to_csv(gwas_filename + '_update1000G.tsv.gz', sep='\t', header=True, index=False, compression='gzip')
print("GWAS data is updated by filtering SNPs present in 1000G and adding the ID. No allele fix is done. The file will be saved as %s" % (gwas_filename+'_update1000G.tsv.gz'))



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')