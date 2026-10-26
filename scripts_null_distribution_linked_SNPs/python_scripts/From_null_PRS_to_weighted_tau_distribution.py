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
import pysam



# ------------------ PARSE ARGUMENTS ------------------ 
# Create a parser object
parser = argparse.ArgumentParser(description="This function calculates weighted tau by y-variable bins from null genome-wide PRS for traits in a given genetic cluster.")
# Add arguments
parser.add_argument("--null_PRS_cluster_results", type=str, required=True, help="A tsv file where the first column is called 'Traits' and each row represents a replicate for a given trait with null proportion of risk alleles genome-wide. Therefore, the rest of columns are aDNA sample IDs that will match metadata.")
parser.add_argument("--n_cpus", type=int, required=True, default=None, help="This test uses parallel threding")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs (it must end with a slash)")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add to the filename when saving outputs without file extension (ex: ASD_prototypes)")
parser.add_argument("--metadata_samples", type=str, required=True, help="File in tsv format that include sample IDs to be included and the y variable to be tested against risk allele proportion (ex: Years BP)")
parser.add_argument("--y_variable", type=str, required=True, help="Name of the column in the metadata file for the y variable to be tested against proportion of risk alleles")
parser.add_argument("--log_filename", type=str, required=True, help="Path to log file where the log info can be appended")
parser.add_argument("--tau_log_observed", type=str, required=True, help="Path to log file where the results of the tau permutation test were recorded, included the observed tau. Two columns are required: 'trait' and 'tau_obs_res'")
parser.add_argument("--bin_start", type=int, required=True, help="If you want to calculate a weighted Tau, bins are needed for y-variable. This is the start value (ex: 0 years BP)")
parser.add_argument("--bin_end", type=int, required=True, help="If you want to calculate a weighted Tau, bins are needed for y-variable. This is the end value (ex: 15.000 years BP)")
parser.add_argument("--bin_size", type=int, required=True, help="If you want to calculate a weighted Tau, bins are needed for y-variable. This is the bin size to group observations (ex: within 500 years BP each)")
parser.add_argument("--filt_outlier_bins_by_y", type=int, required=False, help="If you want, some Bin categories can be removed from the weighted analysis because they give too much weight to some individuals")



# Parse the arguments
args = parser.parse_args()

# Store arguments
null_PRS_res_file = args.null_PRS_cluster_results
out_path = args.out_dir
filename_out = args.info_outfile
cpus = args.n_cpus
metadata_samples = args.metadata_samples
y_variable = args.y_variable
logfile = args.log_filename
tau_observed_file = args.tau_log_observed

# Weighting arguments
print("Tau calculation using weights is performed. Metadata of y-variable will be split into bins to calculate weights.")
print("Please, be sure to provide the following variables to create bins: start, end, size")
start = args.bin_start
end = args.bin_end
size = args.bin_size
print("*")
print("*")
print("*")

# Additional optional arguments for bin and weights management
filt_bins_o = args.filt_outlier_bins_by_y



# ------------------ RUN THE ANALYSIS ------------------ 
# get the start time
st = time.time()



# READ
# Read results of null PRS calculation for a set of traits from the same genetic cluster
null_PRS_cluster = pd.read_csv(null_PRS_res_file, sep="\t")
# Read aDNA metadata
metadata_trend_analysis = pd.read_csv(metadata_samples, sep="\t", encoding="utf-8", encoding_errors="ignore")
keep_samples = list(metadata_trend_analysis['Sample_ID'].tolist())
print("Warning! Tau calculation will be performed with samples that are depicted in the metadata, if you want to include or remove samples change this file or dataframe.")
print("Number of samples with null PRS %d. Number of samples in the metadata %d" % ((null_PRS_cluster.shape[1]-1), len(keep_samples)))
print("*")
print("*")
print("*")

# Read observed tau results
tau_observed = pd.read_csv(tau_observed_file, sep="\t")


# RUN 
# weight samples
w_metadata = avg.split_metadata_bins_add_weights(metadata_trend_analysis, y_variable, size, start, end)

# Prepare bin weights with or without filtering strategy
filter_strategy="None"
if filt_bins_o is not None:
    print("Bins from y-value cutoff of %d will be filtered out before weighted tau calculation" % (filt_bins_o))
    keep_samples = avg.filt_bin_outliers_by_y_variable(keep_samples, y_variable, filt_bins_o, w_metadata)
    print("Final number of samples is: %d" % (len(keep_samples)))
    filter_strategy="filtered bin outliers with "+str(filt_bins_o)+" cutoff"

    # Filter samples from weighted data
    w_metadata = w_metadata.loc[w_metadata['Sample_ID'].isin(keep_samples),:]


# Run weighted tau calculation for each one
res_null_PRS_wTau = test = avg.from_null_PRS_to_weighted_kendall_list_of_traits(null_PRS_cluster, cpus, w_metadata, y_variable)

# Calculate p-values for all tested traits using their computed nulls with observed tau
list_of_tested_traits = res_null_PRS_wTau['Trait'].unique().tolist()
test_type = "Weighted tau with null PRS genome-wide"
for t in list_of_tested_traits:
    res_t = res_null_PRS_wTau.loc[res_null_PRS_wTau['Trait'] == t,:]
    tau_t = tau_observed.loc[tau_observed['trait'].str.startswith(t, na=False),'tau_obs_res'].iloc[0]
    n_rep = res_t.shape[0]
    pval_greater = (res_t['Tau_weighted'] >= tau_t)
    if any(pval_greater):
        pval_greater = pval_greater.value_counts()[True] / n_rep
    else:
        pval_greater = 1 / (n_rep+1)
    
    pval_lower = (res_t['Tau_weighted'] <= tau_t)
    if any(pval_lower):
        pval_lower = pval_lower.value_counts()[True] / n_rep
    else:
        pval_lower = 1 / (n_rep+1)
    
    # Add tau results to logfile
    command =  f"printf '{t}\\t{tau_t}\\t{pval_greater}\\t{pval_lower}\\t{n_rep}\\t{test_type}\\t{filter_strategy}\\n' >> {logfile}"
    os.system(command)

# Prints
print("Weighted tau calculation between null PRS genome-wide and %s performed. Summary statistics were saved at %s" % (y_variable, logfile))


# Save output df
res_null_PRS_wTau.to_csv(out_path + 'Null_distribution_weighted_tau_using_genome_wide_null_PRS_SNPs_50Kb_dist_to_genes_matching_AF_' + filename_out + '.tsv', sep='\t', header=True, index=False)
print("Null distributions of weighted tau using null PRS for all tested traits were saved at %s" % (out_path))



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')