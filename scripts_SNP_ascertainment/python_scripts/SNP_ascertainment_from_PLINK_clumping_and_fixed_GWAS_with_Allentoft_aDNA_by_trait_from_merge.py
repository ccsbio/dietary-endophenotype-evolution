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
parser = argparse.ArgumentParser(description="This function creates a python dict where keys are index SNPs and values clumped ones by using a GWAS file and the clumped output from plink.")
# Add arguments
parser.add_argument("--gwas_data_file", type=str, required=True, help="The path and name of a tsv file with gwas data for a set of snps and phenotypes (rows and columns respectively). This data contains SNPs clumped with plink, including a column SNP ID matching the .clumped file")
parser.add_argument("--plink_clumped_file", type=str, required=True, help="The path and name of a tsv file with .clumped output from --clump by plink 1.9")
parser.add_argument("--phenotype_name", type=str, required=True, help="The phenotype code that can be found in columns of the GWAS dataframe")
parser.add_argument("--pval_threshold", default=None, required=False, help="P.value threshold to consider SNPs as associated to the target trait")
parser.add_argument("--log_filename", type=str, required=True, help="Path to log file where the log info can be appended")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add when saving outputs without file extension (ex: neolitic_diet_90_indpMarkers)")
parser.add_argument("--regions_qc", type=str, default=None, required=False, help="If a full path for .tsv file holding regions to filter out is provided, just the most significant SNP of that region is kept")


# Parse the arguments
args = parser.parse_args()

# Store arguments
file_gwas = args.gwas_data_file
plink_res = args.plink_clumped_file
trait = args.phenotype_name
pval_threshold = args.pval_threshold
logfile = args.log_filename
out_path = args.out_dir
filename_out = args.info_outfile
qc_file = args.regions_qc



# ------------------ RUN THE ANALYSIS ------------------ 
# get the start time
st = time.time()



# Read GWAS merged and fixed data
df_gwas = pd.read_csv(file_gwas, sep='\t')
df_clump = pd.read_csv(plink_res, sep='\t')
original_number_SNPs = df_clump.shape[0]
print(df_gwas.shape)
if pval_threshold != None:
    pval_threshold = float(pval_threshold)
else:
    pval_threshold = df_clump['P'].max()


# First check out if QC is needed
if qc_file != None:
    print("Warning! QC filtering enabled. A .tsv file is expected to have at least the following columns: chr, start, end)")
    df_qc = pd.read_csv(qc_file, sep='\t')
    df_clump = avg.filt_out_plink_clump_qc(df_qc, df_clump)


# Extract GWAS data from merged df given the input trait
# 1. Extract SNPs uniquely associated to a trait
pvals = wvg.extract_field_merged_GWAS(df_gwas, 'p.value')
unique_SNPs, n_unique_by_trait = wvg.number_of_unique_SNPs_associated_by_trait_and_pval_threshold(pvals, pval_threshold)
associated_SNPs = unique_SNPs.loc[unique_SNPs['p.value_' + trait] == True,:].index
# 2. Now subset columns with trait data
df_trait = wvg.extract_GWAS_from_merged(df_gwas, trait, 'ID_1000G')
df_trait_associated = df_trait.loc[associated_SNPs,:]
# 3. Before ascertainment be sure that data is ok
df_trait_associated = wvg.remove_pheno_from_colnames(df_trait_associated, trait)
df_trait_associated = wvg.fix_data_genomic_coordinates_to_int('chr', 'pos', df_trait_associated)
print("Shape of processed GWAS and number of uniquely associated SNPs (rows have to match)")
print(df_trait_associated.shape)
print(len(associated_SNPs))


# If there are no SNPs associated to the current trait for a given pval threshold, the creation of the SNP dictionary is interrupted
if df_trait_associated.shape[0] == 0:
    snps=0
    command =  f"printf '{trait}\\t{snps}\\t{snps}\\t{snps}\\t{pval_threshold}\\n' >> {logfile}"
    os.system(command)
    raise ValueError("Error! The p.value threshold of %f does not return SNPs associated to the trait: %s. Please check this out, process interrupted..." %(pval_threshold, trait))


# Now perform SNP ascertainment if clump also has independent markers for pval threshold
try:
    SNP_dict = avg.create_SNP_dict_from_PLINK_clumping_and_GWAS_with_merged_info_from_clumping_dataset(df_clump, df_trait_associated, pval_threshold=pval_threshold)
except ValueError:
    snps = 0
    command = f"printf '{trait}\\t{snps}\\t{snps}\\t{snps}\\t{pval_threshold}\\n' >> {logfile}"
    os.system(command)
    raise


# Stats
n_items_ascertainment = avg.count_total_items_dict(SNP_dict)
n_prot = len(SNP_dict.keys())
bool_snps_associated = df_trait_associated['p.value'] < pval_threshold
if any(bool_snps_associated):
    n_snps_associated = bool_snps_associated.value_counts()[True]
else:
    n_snps_associated = 0

# logs
command =  f"printf '{trait}\\t{n_snps_associated}\\t{n_items_ascertainment}\\t{n_prot}\\t{pval_threshold}\\n' >> {logfile}"
os.system(command)

# Save dict
wvg.save_python_dict(out_path + "dict_Allentoft_SNP_IDs_"+ str(pval_threshold) + "pval_" + trait + "_" + str(n_prot)+"indpMarkers" + "_" + filename_out + ".pkl", SNP_dict)

# PRINTS
print("A total of %d independent SNPs have been ascertained for trait %s from %d indep SNPs clumped with plink" % (n_prot, trait, original_number_SNPs))



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')