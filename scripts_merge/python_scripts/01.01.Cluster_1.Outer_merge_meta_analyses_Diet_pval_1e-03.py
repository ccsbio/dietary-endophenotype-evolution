import sys
sys.path.append("/home/csic/gpo/olg/ccs/packages/vegeta_last_version")
import merge_GWAS_phenos as vg
import work_vegeta as wvg
import pandas as pd
import time
import argparse
import time


# ------------------ PARSE ARGUMENTS ------------------ 
# Create a parser object
parser = argparse.ArgumentParser(description="This function performs outer merging of GWAS data between SNPs significantly associated by using a pval threshold")
# Add arguments
parser.add_argument("--array_file", type=str, required=True, help=".tsv file with GWAS data that has to be processed. Check array exaxmples for vegeta")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs. Please finish with a '/'")
parser.add_argument("--pval_threshold", type=float, required=True, help="P.value threshold for considering SNPs as significantly associated")

# Parse the arguments
args = parser.parse_args()

# Store arguments
phenos_file = args.array_file
dir_path = args.out_dir
pval_threshold = args.pval_threshold

# Run the program for merging
# get the start time
st = time.time()

merged_df = vg.outer_merge_GWAS_summary_statistics_by_trait_and_pval(dir_path, phenos_file, pval_threshold, filter_INDELS=True, only_autosomal=True, only_biallelic=True, save_tsv_liftOver=True)
merged_df = wvg.fix_data_genomic_coordinates_to_int('chr', 'pos', merged_df)

# Check for dup genomic region after merge
dups = merged_df.duplicated(subset=['chr', 'pos'])
if any(dups):
    n_dups = dups.value_counts()[True]
    if n_dups <10:
        print("Warning! additional merged file will be saved in the same out directory, because automatic filtering is performed.")
        print("There are %d genomic regions duplicated after filtering. The first row is kept." % (n_dups))
        merged_df_filt = merged_df.loc[dups == False,:].copy()
        merged_df_filt.to_csv(path_or_buf=dir_path + 'Merged_GWAS_meta_analyses_Diet_Cluster_1_logTG_prot_outer_pval_1e-03_dup_filt.tsv.gz', header=True, index=False, decimal='.', sep='\t', compression={'method': 'gzip', 'compresslevel': 4})
    else:
        print("Warning! more than 10 genomic regions are duplicated after the outer merge. Please, check it out.")

# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')


# Write the results
merged_df.to_csv(path_or_buf=dir_path + 'Merged_GWAS_meta_analyses_Diet_Cluster_1_logTG_prot_outer_pval_1e-03.tsv.gz', header=True, index=False, decimal='.', sep='\t', compression={'method': 'gzip', 'compresslevel': 4})