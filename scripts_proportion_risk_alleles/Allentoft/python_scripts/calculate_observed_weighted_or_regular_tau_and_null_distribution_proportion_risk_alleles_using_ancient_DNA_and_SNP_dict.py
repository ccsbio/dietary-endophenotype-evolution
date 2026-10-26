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
parser = argparse.ArgumentParser(description="This function calculates proportion of risk allele in VCF with aDNA samples by using a SNP dict of independent prototypes and linked markers. Also, Kendall tau is calculated between this and y-variable (usually time)")
# Add arguments
parser.add_argument("--in_dir", type=str, required=True, help="The input path with VCF files to use in the analysis")
parser.add_argument("--file", type=str, required=True, help="A VCF file with data for a set of samples. This vcf must be sorted and indexed")
parser.add_argument("--tabix", type=str, required=True, help="A tabix file for the input VCF.")
parser.add_argument("--SNP_dict", type=str, required=True, help="A .pkl file for a python dictionary where keys are prototype SNPs and values are linked markers")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add when saving outputs without file extension (ex: neolitic_diet_90_indpMarkers)")
parser.add_argument("--metadata_samples", type=str, required=True, help="File in tsv format that include sample IDs to be included and the y variable to be tested against risk allele proportion (ex: Years BP)")
parser.add_argument("--y_variable", type=str, required=True, help="Name of the column in the metadata file for the y variable to be tested against proportion of risk alleles")
parser.add_argument("--log_filename", type=str, required=True, help="Path to log file where the log info can be appended")
parser.add_argument("--number_replicates", type=int, required=True, help="Number of replicates for null distribution")
parser.add_argument("--cpus", type=int, required=True, help="Number of cpus to compute the null distribution")
parser.add_argument("--seed", type=int, required=True, default=None, help="A seed is required to create as many BitGenerators as replicates using spawn from numpy. This ensures independent random sampling")
parser.add_argument("--weight_y_var", default=False, action="store_true", help="If set to True, then weights based on number of obervations per bin of y-variable are used to calculate Tau")
parser.add_argument("--bin_start", type=int, required=False, help="If you want to calculate a weighted Tau, bins are needed for y-variable. This is the start value (ex: 0 years BP)")
parser.add_argument("--bin_end", type=int, required=False, help="If you want to calculate a weighted Tau, bins are needed for y-variable. This is the end value (ex: 15.000 years BP)")
parser.add_argument("--bin_size", type=int, required=False, help="If you want to calculate a weighted Tau, bins are needed for y-variable. This is the bin size to group observations (ex: within 500 years BP each)")
parser.add_argument("--filt_outlier_bins_by_y", type=int, required=False, help="If you want, some Bin categories can be removed from the weighted analysis because they give too much weight to some individuals")
parser.add_argument("--merge_Bins_in_tails", default=False, action="store_true", help="If set to True, then weights are recalculated by merging bins with one obeservation detected as a gap in the end of the distribution")
parser.add_argument("--bimodal_bins_by_y", type=int, required=False, help="If your samples are distributed as a bimodal and you want to decide the cutoff of the weights of both distributions use this option")


# Parse the arguments
args = parser.parse_args()

# Store arguments
in_dir = args.in_dir
vcf = args.file
dict_SNP = args.SNP_dict
out_path = args.out_dir
filename_out = args.info_outfile
w_bool = args.weight_y_var
metadata_samples = args.metadata_samples
y_variable = args.y_variable
logfile = args.log_filename
tabix = args.tabix
rep = args.number_replicates
cpus = args.cpus
seed = args.seed

# Additional optional arguments for bin and weights management
filt_bins_o = args.filt_outlier_bins_by_y
merge_bins = args.merge_Bins_in_tails
bimodal_b = args.bimodal_bins_by_y


print("Warning! this script performs Kendall Tau calculation of proportion of risk alleles against an y-variable. Please be sure to add a .tsv file in the --metadata_samples having the values to correlate with")
print("Warning! the optional arguments --filt_outlier_bins_by_y, --merge_Bins_in_tails OR --bimodal_bins_by_y just work for weighted tau")
print("The y-variable name is %s" % (y_variable))
print("Path/name of logfile is required to add Tau permutation results")
print("*")
print("*")
print("*")


if w_bool == True:
    print("Tau calculation with weights has been settled. Metadata of y-variable will be split into bins to calculate weights.")
    print("Please, be sure to provide the following variables to create bins: start, end, size")
    start = args.bin_start
    end = args.bin_end
    size = args.bin_size
    print("*")
    print("*")
    print("*")
    
    # If there are optional arguments provided for bin managament, just one can be provided at a time
    modes = [
    filt_bins_o is not None,
    merge_bins,
    bimodal_b is not None,]
    
    if sum(modes) > 1:
        raise ValueError(
            "Choose ONLY ONE of:\n"
            "--filt_outlier_bins_by_y\n"
            "--merge_Bins_in_tails\n"
            "--bimodal_bins_by_y"
            )



# ------------------ RUN THE ANALYSIS ------------------ 
# get the start time
st = time.time()



# SNP dict
SNP_selected = wvg.load_python_dict(dict_SNP)
prototypes = list(SNP_selected.keys())

# Fetch VCF. It is important to have the tabix file and a sorted VCF
print("Warning! Remember that the VCF must be sorted and indexed. Then, a tabix file is also needed. For memmory reasons, SNPs are fetched instead of loading the whole file.")
print("*")
print("*")
print("*")
df_vcf = avg.get_pandas_df_from_indexed_and_bigger_VCF_using_python_dict(SNP_selected, in_dir+vcf, in_dir+tabix)
df_vcf = wvg.fix_data_genomic_coordinates_to_int('CHROM','POS', df_vcf)
print(df_vcf.shape)

# Prepare VCF including coded info for risk SNP
VCF_with_risk_info = avg.subset_VCF_with_SNP_dict_data_add_risk_info(SNP_selected, df_vcf)
# Store sample IDs in a list
samples = VCF_with_risk_info.columns[12:]


# Read metadata for Tau calculation
metadata_trend_analysis = pd.read_csv(metadata_samples, sep="\t", encoding="utf-8", encoding_errors="ignore")
keep_samples = list(VCF_with_risk_info.columns[VCF_with_risk_info.columns.isin(metadata_trend_analysis['Sample_ID'].tolist())])
print("Warning! Tau calculation will be performed with samples that are depicted in the metadata, if you want to include or remove samples change this file.")
print("Number of samples in the VCF %d. Number of samples in the metadata %d" % (len(samples), len(keep_samples)))
print("*")
print("*")
print("*")


filter_strategy="None"
# Now get observed tau and null distribution
print(("Calculating proportion of risk alleles for %d samples with %d independent markers") % (len(keep_samples), len(SNP_selected.keys())))
if w_bool == False:
    # Observed tau
    tau_allentoft, pval, props, res_prototypes = avg.observed_Tau_corr_proportion_risk_alleles(VCF_with_risk_info, metadata_trend_analysis, y_variable, keep_samples, prototypes, SNP_selected)
    perc_t = "NaN"
    # Null distribution
    test_null = avg.null_distribution_Tau_corr_proportion_risk_alleles_parallel(rep, cpus, seed, VCF_with_risk_info, metadata_trend_analysis, y_variable, keep_samples, prototypes, SNP_selected)
else:
    # First bin y-variable and get weights
    # If you have a bimodal distribution of your samples weights can be done with 2 bins regarding one cutoff
    # Otherwise you can do it by uniform bins of x size
    if bimodal_b is not None:
        print("Warning! bimodal weight calculation enabled with an y-variable cutoff of %d" %(bimodal_b))
        w_metadata = avg.split_metadata_bins_bimodal_distribution_and_weights(metadata_trend_analysis, y_variable, bimodal_b)
        print(w_metadata['Bin'].value_counts())
        filter_strategy="bimodal bins with "+str(bimodal_b)+" cutoff"
    else:
        w_metadata = avg.split_metadata_bins_add_weights(metadata_trend_analysis, y_variable, size, start, end)
        print(w_metadata['Bin'].value_counts())

    # In this step we can either filter outlier bins or merge weights of a wider unrepresented area (in the tail of the bin distribution)
    if filt_bins_o is not None:
        print("Bins from y-value cutoff of %d will be filtered out before weighted tau calculation" % (filt_bins_o))
        keep_samples = avg.filt_bin_outliers_by_y_variable(keep_samples, y_variable, filt_bins_o, w_metadata)
        print("Final number of samples is: %d" % (len(keep_samples)))
        filter_strategy="filtered bin outliers with "+str(filt_bins_o)+" cutoff"
    
    if merge_bins == True:
        w_metadata = avg.check_tails_outliers_and_merge_bins(w_metadata, y_variable)
        print(w_metadata['Bin'].value_counts())
        filter_strategy="merge tail of bin outliers"

    # Observed tau
    tau_allentoft, props, res_prototypes, perc_t = avg.observed_Tau_weighted_corr_proportion_risk_alleles(VCF_with_risk_info, w_metadata, y_variable, keep_samples, prototypes, SNP_selected)
    pval="NaN"
    # Null distribution
    test_null = avg.null_distribution_weighted_Tau_corr_proportion_risk_alleles_parallel(rep, cpus, seed, VCF_with_risk_info, w_metadata, y_variable, keep_samples, prototypes, SNP_selected)


# Save results
test_null.to_csv(out_path + 'Null_distribution_Kendall_tau_prop_risk_' + filename_out + '.tsv', sep='\t', index=False, header=True)
print("*")
print("*")
print("*")
print("For file %s..." % (vcf))
if w_bool == False:
    pval_smaller = (test_null.loc[test_null['Tau'] < tau_allentoft,:].shape[0] / rep) 
    pval_greater = (test_null.loc[test_null['Tau'] > tau_allentoft,:].shape[0] / rep)
    print("Observed tau is: %e" % (tau_allentoft))
    print("Null distribution for Kendall tau correlation between proportion of risk alleles and %s computed. Saved at %s" % (y_variable, out_path))
else:
    pval_smaller = (test_null.loc[test_null['Tau_weighted'] <= tau_allentoft,:].shape[0] / rep) 
    pval_greater = (test_null.loc[test_null['Tau_weighted'] >= tau_allentoft,:].shape[0] / rep)
    print("Observed weighted tau is: %e" % (tau_allentoft))
    print("Null distribution for weighted Kendall tau correlation between proportion of risk alleles and %s computed. Saved at %s" % (y_variable, out_path))
print("One-tailed p.value when equal or smaller is: %e" % (pval_smaller))
print("One-tailed p.value when equal or greater is: %e" % (pval_greater))


# Add summarized results to logfile
command =  f"printf '{filename_out}\\t{tau_allentoft}\\t{pval_smaller}\\t{pval_greater}\\t{rep}\\t{w_bool}\\t{filter_strategy}\\n' >> {logfile}"
os.system(command)



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')