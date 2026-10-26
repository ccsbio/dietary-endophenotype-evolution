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
parser = argparse.ArgumentParser(description="This function creates a bedfile from GWAS data.")
# Add arguments
parser.add_argument("--merged_gwas", type=str, required=True, help="Full path to the .tsv file with merged GWAS data that has to be processed (expected a cluster of GWAS). Check array exaxmples for vegeta")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs. Please finish with a '/'")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add when saving outputs without file extension (ex: cluster1_logTG_prot)")

# Parse the arguments
args = parser.parse_args()

# Store arguments
merged_file = args.merged_gwas
dir_path = args.out_dir
file_out = args.info_outfile

# Run the program for merging
# get the start time
st = time.time()

# Read merged GWAS
df_clust =  pd.read_csv(merged_file, sep='\t')

# Create bedfile and add and ID 
bedfile_clust = wvg.create_bedfile_from_df(df_clust, 'chr', 'pos')
bedfile_clust['ID'] =  df_clust['chr'].astype(str) + "_" + df_clust['pos'].astype(str) + "_" + df_clust['ref'] + "_" + df_clust['alt']

# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')


# Write the results
bedfile_clust.to_csv(path_or_buf=dir_path + file_out + '.bed', header=False, index=False, sep='\t')