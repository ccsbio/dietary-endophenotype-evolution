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
parser = argparse.ArgumentParser(description="This function works with a cluster of genetically correlated traits with GWAS data for 1000G markers that were merged with vegeta package. A given trait info is subset to create a file with .assoc format of uniquely associated SNPs (non-pleiotropic) to perform clumping with plink")
# Add arguments
parser.add_argument("--gwas_data_file", type=str, required=True, help="The path and name of a tsv file with gwas data for a set of snps and phenotypes (rows and columns respectively). It is expected to contain an ID for 1000G data.")
parser.add_argument("--phenotype_name", type=str, required=True, help="The phenotype code that can be found in columns of the GWAS dataframe")
parser.add_argument("--pval_threshold", type=float, required=True, help="P.value threshold to consider SNPs as associated to the target trait")
parser.add_argument("--log_filename", type=str, required=True, help="Path to log file where the log info can be appended")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add when saving outputs without file extension (ex: neolitic_diet_90_indpMarkers)")

# Parse the arguments
args = parser.parse_args()

# Store arguments
#in_dir = args.in_dir
file_gwas = args.gwas_data_file
trait = args.phenotype_name
pval_threshold = args.pval_threshold
logfile = args.log_filename
out_path = args.out_dir
filename_out = args.info_outfile



# ------------------ RUN THE ANALYSIS ------------------ 
# get the start time
st = time.time()



# Read GWAS merged and fixed data
print("Warning! It is important that input GWAS file depicts a column called 'ID_1000G' that matches an ID found in VCFs of 1000G processed data")
df_gwas = pd.read_csv(file_gwas, sep='\t')
print(df_gwas.shape)

# 1. Extract SNPs uniquely associated to a trait
pvals = wvg.extract_field_merged_GWAS(df_gwas, 'p.value')
#print(pvals)
unique_SNPs, n_unique_by_trait = wvg.number_of_unique_SNPs_associated_by_trait_and_pval_threshold(pvals, pval_threshold)
associated_SNPs = unique_SNPs.loc[unique_SNPs['p.value_' + trait] == True,:].index
n_snps_associated = len(associated_SNPs)

# 2. Now subset GWAS data
df_trait = wvg.extract_GWAS_from_merged(df_gwas, trait, 'ID_1000G')
df_trait_associated = df_trait.loc[associated_SNPs,:]

# Before ascertainment be sure that data is ok
df_trait_associated = wvg.remove_pheno_from_colnames(df_trait_associated, trait)
df_trait_associated = wvg.fix_data_genomic_coordinates_to_int('chr', 'pos', df_trait_associated)
print(df_trait_associated.shape)
print(n_snps_associated)

# 3. Create .assoc file for plink
if df_trait_associated.shape[0] > 0:
    assoc_plink = df_trait_associated.loc[:,['ID_1000G', 'p.value']].copy() # Use ID expected in the VCF, not in the GWAS
    assoc_plink = assoc_plink.rename(columns={"ID_1000G": "SNP", "p.value": "P"})
    print("The number of associated SNPs to trait %s is: %d. The .assoc file will be saved as %s" % (trait, n_snps_associated, out_path+filename_out+".assoc"))
    assoc_plink.to_csv(out_path + filename_out + '.assoc', sep=' ', header=True, index=False)
else:
    snps=0
    command =  f"printf '{trait}\\t{snps}\\t{pval_threshold}\\n' >> {logfile}"
    os.system(command)
    raise ValueError("Warning! No SNPs were found uniquely associated to trait %s. Please check this out, process interrupted..." % (trait))


# logs
command =  f"printf '{trait}\\t{n_snps_associated}\\t{pval_threshold}\\n' >> {logfile}"
os.system(command)



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')