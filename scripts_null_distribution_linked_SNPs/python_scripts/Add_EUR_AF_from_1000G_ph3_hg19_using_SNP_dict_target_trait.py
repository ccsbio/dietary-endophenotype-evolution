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
parser = argparse.ArgumentParser(description="This function calculates proportion of risk allele in VCF with aDNA samples by using a SNP dict of independent prototypes and linked markers")
# Add arguments
parser.add_argument("--SNP_dict", type=str, required=True, help="A .pkl file for a python dictionary where keys are prototype SNPs and values are linked markers")
parser.add_argument("--vcf_file", type=str, required=True, help="A VCF file with data for a set of samples. This test is thought for dilpoid data info. Also, it expected to found SNPs from the null set")
parser.add_argument("--tabix_file", type=str, required=False, default=None, help="An index file for the sorted vcf. It is optional, but recommended")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs (it must end with a slash)")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add to the filename when saving outputs without file extension (ex: ASD_prototypes)")


# Parse the arguments
args = parser.parse_args()

# Store arguments
SNP_dict_f = args.SNP_dict
file_vcf_1000G = args.vcf_file
file_1000G_tabix = args.tabix_file
out_path = args.out_dir
filename_out = args.info_outfile



# ------------------ RUN THE ANALYSIS ------------------ 
# get the start time
st = time.time()



# Read SNP dict
SNP_dict = wvg.load_python_dict(SNP_dict_f)

# Fetch from 1000G VCF target SNPs
vcf_1000G = avg.get_pandas_df_from_indexed_and_bigger_VCF_using_python_dict(SNP_dict, file_vcf_1000G, file_1000G_tabix, samples=False)
vcf_1000G, indels_counter = wvg.filter_INDELS(vcf_1000G, ['REF', 'ALT'])
vcf_1000G = vcf_1000G.rename(columns={'CHROM':'chr', 'POS':'pos'})
print('Number of fetched positions from 1000G after INDEL filtering is: %d' % (vcf_1000G.shape[0]))
print('VCF columns recovered include:')
print(vcf_1000G.columns)
#print(vcf_1000G.iloc[0:4,0:6])

# Prepare prototype data for target trait
prototypes = pd.Series(SNP_dict.keys())
print("Dict length:")
print(len(prototypes) + len(avg.extract_values_from_dict_to_list(SNP_dict)))
#print(prototypes[0])
print("\n")
#print(avg.extract_values_from_dict_to_list(SNP_dict)[0:3])
df_prototypes = pd.DataFrame({'chr': prototypes.str.split('_').str[0], 'pos': prototypes.str.split('_').str[1], 'ref': prototypes.str.split('_').str[2], 'alt':prototypes.str.split('_').str[3], 'beta_sign':prototypes.str.split('_').str[4]})
df_prototypes = wvg.fix_data_genomic_coordinates_to_int('chr', 'pos', df_prototypes)
# You will need this columns to have allele annotation as in your data set (gwas or 1000G can report the complementary bases of the variant)
df_prototypes = df_prototypes.assign(original_ref = df_prototypes['ref'], original_alt = df_prototypes['alt'])
print(df_prototypes.iloc[0:4,])
print('Number of independent SNPs for %s is: %d' % (filename_out ,df_prototypes.shape[0]))


# Merge 1000G data with prototypes so you can add AF regarding the ALT
merged = wvg.merge_and_fix_SNPs_between_two_datasets(df_prototypes, vcf_1000G, 'chr', 'pos', 'ref', 'alt', 'REF', 'ALT')
print("\n")
print("\n")
#print("Same alleles:")
print(merged['Same_alleles'].value_counts())
print("\n")
print("\n")

if not all(merged['Same_alleles']):
    merged = merged.loc[merged["Same_alleles"] == True,:]
    print('Number of mismatching alleles filtered is %d' % (merged['Same_alleles'].value_counts()[False]))


# Now get AF for EUR populations for the ALT in 1000G
AF_EUR = wvg.parse_metadata_1000G_from_VCF_INFO(merged, 'INFO', AA=False, AF=True, POP="EUR_AF")
filt = merged.loc[:,['chr', 'pos', 'ref', 'alt', 'REF', 'ALT', 'original_ref', 'original_alt', 'beta_sign']]
filt['ALT_AF'] = AF_EUR.astype(float)


# Complementary bases were fixed with the previous function to match 1000G annotation, so you can just check if there are switched alleles to assign correct AF
switched = (filt['ref'] == filt['ALT']) & (filt['ALT'] == filt['ref'])
print('Switched alleles in merging:')
print(switched.value_counts())


# We assume biallelic positions, so just 1 derived allele
# We put AF regarding the allele in alt column
if any(switched):
    filt.loc[switched, 'ALT_AF'] = 1 - filt.loc[switched, 'ALT_AF'] 

# Now we have to be sure of which is the AF that matches the risk allele of the target trait
other_al_risk = filt['beta_sign'] == '-'
filt['Risk_allele_AF'] = filt['ALT_AF']
filt.loc[other_al_risk,'Risk_allele_AF'] = 1 - filt.loc[other_al_risk, 'ALT_AF']

# Finally remove columns from 1000G and use allele annotation as in your original dataset
filt = filt.drop(['REF', 'ALT', 'ref', 'alt'], axis=1)
filt = filt.rename(columns={'original_ref':'ref', 'original_alt':'alt'})

# Save prototypes with AF of european population from 1000G
print('AF successfully added from 1000G data. File is saved at %s' % (out_path))
filt.to_csv(out_path + 'SNPs_' + filename_out + '_' + str(len(SNP_dict.keys())) + '_indpMarkers_with_EUR_AF_risk_alleles.tsv', sep="\t", header=True, index=False)



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')