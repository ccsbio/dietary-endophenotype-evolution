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
parser = argparse.ArgumentParser(description="This function calculates proportion of risk allele in VCF with aDNA samples by using a SNP dict of independent prototypes and linked markers")
# Add arguments
parser.add_argument("--null_set_file", type=str, required=True, help="A tsv file with unique genetic coordinates of biallelic SNPs that will match a VCF (columns needed 'chr', 'pos', 'ALT_AF'). These snps are used to sample proxy alleles (controls) regarding prorotypes used in PRS calculation. The null set is expected to contain only those markers in close distance to genes (smaller of 50Kb). It also must contain at least the AF of the alternative allele for the studied population")
parser.add_argument("--prototype_SNPs_file", type=str, required=True, help="A tsv file with unique genetic coordinates of a set of independet markers associated to a target trait. It is expected only one prototype per position (columns needed 'chr', 'pos', 'Risk_allele_AF'). A column called 'Risk_allele_AF' holding AF of the risk allele used during observed PRS calculation in the target population is expected.")
parser.add_argument("--vcf_file", type=str, required=True, help="A VCF file with data for a set of samples. This test is thought for dilpoid data info. Also, it expected to found SNPs from the null set")
parser.add_argument("--tabix_file", type=str, required=False, default=None, help="An index file for the sorted vcf. It is optional, but recommended")
parser.add_argument("--n_replicates", type=int, required=False, default=1000, help="Number of replicates for the null distribution (by default 1000). The null PRS will be calculated with proxy alleles N times across all samples")
parser.add_argument("--n_cpus", type=int, required=True, default=None, help="This test uses parallel threding")
parser.add_argument("--seed", type=int, required=True, default=None, help="A seed is required to create as many BitGenerators as replicates using spawn from numpy. This ensures independent random sampling")
parser.add_argument("--AF_bin_span", type=float, required=True, help="Range to consider or deviation from the AF of the alleles used in the observed PRS. This test creates bins of AF to accept markers within this range as proxies when computing null PRS")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs (it must end with a slash)")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add to the filename when saving outputs without file extension (ex: ASD_prototypes)")
parser.add_argument("--samples_file", type=str, required=False, default=None, help="Optional sample file can be provided (.tsv without header and one column were rows are IDs present in the null VCF). Therefore, only sample IDs in the file will be considered for the test. By default, all samples in the VCF header are used")
parser.add_argument("--remove_snps_from_gwas", default=False, action="store_true", help="Optional step for checking all SNPs associated to the target trait using GWAS data for filtering them from the null set of SNPs if needed")
parser.add_argument("--pval_threshold", type=float, required=False, default=1e-05, help="p.value threshold used when filtering associated SNPs to a trait from the null set")
parser.add_argument("--gwas", type=str, required=False, help="GWAS file with a column called 'p.value' and the genomic coordinates in the 'chr' and 'pos' columns used when filtering associated SNPs to a trait from the null set")
parser.add_argument("--trait_label_merged_GWAS", type=str, required=False, default=None, help="If a trait label code is provided, then this program asumes gwas data provided has more than one trait merged with vegeta package")



# Parse the arguments
args = parser.parse_args()

# Store arguments
null_SNPs_file = args.null_set_file
prototype_SNPs_metadata = args.prototype_SNPs_file
vcf = args.vcf_file
tabix = args.tabix_file
out_path = args.out_dir
filename_out = args.info_outfile
filt_snps = args.remove_snps_from_gwas
span = args.AF_bin_span
rep = args.n_replicates
cpus = args.n_cpus
seed = args.seed
samples_file = args.samples_file


# If you have to check if SNPs associated to target trait are present in the dataset, get the rest of required arguments
if filt_snps == True:
    pval_threshold = args.pval_threshold
    gwas_file = args.gwas
    trait_vg = args.trait_label_merged_GWAS
    print("Warning! Filtering option of target trait SNPs from null set has been enabled, so gwas data is needed and a pval threshold. Therefore, associated SNPs will be removed from the null dataset if needed")
    gwas_data = pd.read_csv(gwas_file, sep="\t")
else:
    print("Warning! Filtering option of target trait SNPs from null set by default is disabled. Please be sure that SNPs statistically associated to the target trait are not present in the null dataset")
    print("Otherwise, you can enable the argument --remove_snps_from_gwas and provide a gwas file and a pval threshold, so the filtering is performed before running the test")


# If a set of input sample IDs is requiered, then read the sample file
if samples_file != None:
    print("Warning! Sample IDs list is provided. Please, be sure that the format of the sample file is tab delimited, where each row represents and ID with just one column and no header")
    samples = pd.read_csv(samples_file, sep="\t", header=None)



# ------------------ RUN THE ANALYSIS ------------------ 
# get the start time
st = time.time()



# Read SNP dictorionary of target trait and null SNPs set
target_prototypes_AF = pd.read_csv(prototype_SNPs_metadata, sep="\t")
linked_SNPs_AF = pd.read_csv(null_SNPs_file, sep="\t")

#print(linked_SNPs_AF.iloc[0:3,])
print("\n")
print("\n")
print("\n")

# Since it is expected to show unique genetic positions, create an 'ID' so you can filter target prototypes from the null dataset if necessary
target_prototypes_AF = target_prototypes_AF.assign(snp_id=target_prototypes_AF['chr'].astype(str) + "_" + target_prototypes_AF['pos'].astype(str))
n_prototypes = target_prototypes_AF.shape[0] # check out number of independent positions
list_targets = target_prototypes_AF['snp_id'].tolist()
# since you have the data of 1000G of the AF for the alternative, it is important get also the reference (the test looks for any allele matching the AF)
linked_SNPs_AF['REF_AF'] = 1 - linked_SNPs_AF['ALT_AF'] 
linked_SNPs_AF = linked_SNPs_AF.assign(snp_id=linked_SNPs_AF['chr'].astype(str) + "_" + linked_SNPs_AF['pos'].astype(str))
#print(linked_SNPs_AF.iloc[0:3,])


# If filtering option enabled do the following:
if filt_snps == True:
    # Check out associated by pval threshols to the target trait
    if trait_vg != None:
        # Check out SNPs associated to your target trait from merged GWAS data with vegeta
        df_trait = wvg.extract_GWAS_from_merged(gwas_data, trait_vg)
        print(df_trait.iloc[0:3,:])
        df_trait = wvg.remove_pheno_from_colnames(df_trait, trait_vg)
        associated_snps_bool = df_trait['p.value'] < pval_threshold
        associated_snps = df_trait.loc[associated_snps_bool,:] # filtered df
    else:
        associated_snps_bool = gwas_data['p.value'] < pval_threshold
        associated_snps = gwas_data.loc[associated_snps_bool,:] # filtered df
    
    associated_snps = (associated_snps['chr'].astype(str) + "_" + associated_snps['pos'].astype(str)).tolist() # Create an ID
    n_significant = len(associated_snps)
    print("Number of associated SNPs to the target trait when using a pval threshold of %f is: %d" % (pval_threshold, n_significant))

    # Check out if target SNPs are present
    associated_snps = pd.Series(associated_snps)
    bool_shared = associated_snps.isin(linked_SNPs_AF['snp_id'].tolist())
    shared = associated_snps.loc[bool_shared == True].tolist()
 
    n_shared = len(shared)
    if n_shared > 0:
        linked_SNPs_AF = linked_SNPs_AF.loc[~linked_SNPs_AF['snp_id'].isin(shared),:]
        print("A total of %d SNPs were found in the null dataset, they will be removed for downstream analyses" % (n_shared))


print("Number of independent markers for PRS calculation is %d." % (n_prototypes))
print("\n")
print("\n")
print("\n")

# Check out samples to be used. Default are all from the VCFs header
if samples_file != None:
    samples = samples.iloc[:,0].tolist()
else:
    print("By default, test is performed using all the samples available in the VCF header")
    vcf_fetch = pysam.VariantFile(vcf, index_filename=tabix)
    samples = list(vcf_fetch.header.samples)


print("Starting test using a null set of SNPs with size %d. PRS is computed for %d samples and the number of SNPs used each time is %d." % (linked_SNPs_AF.shape[0], len(samples), n_prototypes))
print("Number of replicates and cpus is %d and %d respectively" % (rep, cpus))
# Perform parallelized version of this test
null_PRS = avg.parallel_null_distribution_PRS_using_snps_increasing_trait_value_and_with_linkage_to_genes_by_matching_AF(rep, cpus, seed, linked_SNPs_AF, 
                                                                                                                         target_prototypes_AF, 'Risk_allele_AF', samples, 
                                                                                                                         span, vcf, tabix, n_with_backup=2)


try:
    # Save output df
    null_PRS.to_csv(out_path + 'Null_PRS_' + str(rep) + 'replicates_allentoft_' + str(n_prototypes) + '_independent_SNPs_50Kb_dist_to_genes_matching_AF_' + filename_out + '.tsv', sep='\t', header=True, index=False)
    print("Null distribution of PRS computed. Saved at %s" % (out_path))
    print("Number of replicates that succeded is: %d" % (null_PRS.shape[0]))
except Exception as e:
    # If any test failed, save the object with pickle
    print(".tsv file conversion failed. Instead python object as .pkl is stored at %s" % (out_path))
    print("This happens when any replicate fails, so the dataframe cannot be created when there are null rows")
    name_object = out_path + 'Null_PRS_' + str(rep) + 'replicates_allentoft_' + str(n_prototypes) + '_independent_SNPs_50Kb_dist_to_genes_matching_AF_' + filename_out + '.pkl'
    wvg.save_python_object(out_path, null_PRS)



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')