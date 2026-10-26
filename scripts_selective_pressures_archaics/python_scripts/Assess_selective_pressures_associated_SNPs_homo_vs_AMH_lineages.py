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
parser.add_argument("--vcf_archaic", type=str, required=True, help="A VCF file with pseudohaploid genotype data for one archaic homo individual")
parser.add_argument("--vcf_1000G", type=str, required=True, help="A VCF file with 1000G data for a set of variants. It can be used to check AF or AA from INFO.")
parser.add_argument("--tabix_file_archaic", type=str, required=False, default=None, help="An index file for the sorted vcf. It is optional, but recommended")
parser.add_argument("--tabix_file_1000G", type=str, required=False, default=None, help="An index file for the sorted vcf. It is optional, but recommended")
parser.add_argument("--out_dir", type=str, required=True, help="Output path (without output filename) to save outputs (it must end with a slash)")
parser.add_argument("--info_outfile", type=str, required=True, help="Info to add to the filename when saving outputs without file extension (ex: ASD_prototypes)")
parser.add_argument("--pval_threshold", type=float, required=True, default=1e-05, help="p.value threshold used when filtering associated SNPs to a trait from the null set")
parser.add_argument("--gwas", type=str, required=True, help="GWAS file with a column called 'p.value' and the genomic coordinates in the 'chr' and 'pos' columns used when filtering associated SNPs to a trait from the null set")
parser.add_argument("--trait_label_merged_GWAS", type=str, required=False, default=None, help="If a trait label code is provided, then this program asumes gwas data provided has more than one trait merged with vegeta package")
parser.add_argument("--log_filename", type=str, required=True, help="Path to log file where the log info can be appended")
#parser.add_argument("--alternative", type=str, required=False, default="two-sided", help="Type of alternative hypothesis testing with Fisher exact test")


# Parse the arguments
args = parser.parse_args()

# Store arguments
SNP_dict_f = args.SNP_dict
file_vcf_archaic = args.vcf_archaic
file_archaic_tabix = args.tabix_file_archaic
file_vcf_1000G = args.vcf_1000G
file_1000G_tabix = args.tabix_file_1000G
logfile = args.log_filename

out_path = args.out_dir
filename_out = args.info_outfile
gwas_file = args.gwas
trait_vg = args.trait_label_merged_GWAS
pval_threshold = args.pval_threshold
#alternative = args.alternative


# ------------------ RUN THE ANALYSIS ------------------ 
# get the start time
st = time.time()



# Read SNP dict
SNP_dict = wvg.load_python_dict(SNP_dict_f)

# Fetch from archaic homo VCF target SNPs
vcf_archaic = avg.get_pandas_df_from_indexed_and_bigger_VCF_using_python_dict(SNP_dict, file_vcf_archaic, file_archaic_tabix)
#vcf_archaic, indels_counter = wvg.filter_INDELS(vcf_archaic, ['REF', 'ALT']) # You cannot filter INDELS in archaics VCF. ALT columns has '.'
vcf_archaic = vcf_archaic.rename(columns={'CHROM':'chr', 'POS':'pos'})
homo = vcf_archaic.columns[9]
print('Number of fetched positions from %s is: %d' % (homo, vcf_archaic.shape[0]))
#print('VCF columns recovered include:')
#print(vcf_archaic.iloc[0:4,0:6])


# Fetch from 1000G VCF target SNP
vcf_1000G = avg.get_pandas_df_from_indexed_and_bigger_VCF_using_python_dict(SNP_dict, file_vcf_1000G, file_1000G_tabix, samples=False)
vcf_1000G, indels_counter = wvg.filter_INDELS(vcf_1000G, ['REF', 'ALT'])
vcf_1000G = vcf_1000G.rename(columns={'CHROM':'chr', 'POS':'pos'})
print('Number of fetched positions from 1000G after INDEL filtering is: %d' % (vcf_1000G.shape[0]))
#print('VCF columns recovered include:')
#print(vcf_1000G.columns)
#print(vcf_1000G.iloc[0:4,0:6])


# Prepare prototype data for target trait
list_prototypes = list(SNP_dict.keys())
print("Dict length:")
print(len(list_prototypes) + len(avg.extract_values_from_dict_to_list(SNP_dict)))
#print(list_prototypes[0])
print("\n")
#print(avg.extract_values_from_dict_to_list(SNP_dict)[0:3])
#df_prototypes = pd.DataFrame({'chr': prototypes.str.split('_').str[0], 'pos': prototypes.str.split('_').str[1], 'ref': prototypes.str.split('_').str[2], 'alt':prototypes.str.split('_').str[3], 'beta_sign':prototypes.str.split('_').str[4]})
#df_prototypes = wvg.fix_data_genomic_coordinates_to_int('chr', 'pos', df_prototypes)
print('Number of independent SNPs for %s is: %d' % (filename_out, len(list_prototypes)))


# Read GWAS data
gwas_data = pd.read_csv(gwas_file, sep="\t")
gwas_data = wvg.fix_data_genomic_coordinates_to_int('chr', 'pos', gwas_data)
# Check out SNPs associated to your target trait from merged GWAS data with vegeta
if trait_vg != None:
    df_trait = wvg.extract_GWAS_from_merged(gwas_data, trait_vg)
    #print(df_trait.iloc[0:3,:])
    df_trait = wvg.remove_pheno_from_colnames(df_trait, trait_vg)
    associated_snps_bool = df_trait['p.value'] < pval_threshold
    associated_snps = df_trait.loc[associated_snps_bool,:] # filtered df
else:
    associated_snps_bool = gwas_data['p.value'] < pval_threshold
    associated_snps = gwas_data.loc[associated_snps_bool,:] # filtered df
# add a custom ID is there are not IDs
associated_snps = associated_snps.assign(snp_id=associated_snps['chr'].astype(str) + "_" + associated_snps['pos'].astype(str))


# 1. MERGE 1000G DATA: extract AA
print("\n")
print("\n")
print("Merging 1000G with GWAS data")
AA_EUR_freq = wvg.parse_metadata_1000G_from_VCF_INFO(vcf_1000G, 'INFO', AA=True, AF=True, POP='EUR_AF')
AA_AF_df = vcf_1000G.loc[:,['chr', 'pos', 'REF', 'ALT']]
AA_AF_df = pd.concat([AA_AF_df, AA_EUR_freq], axis=1)
AA_AF_df['AF'] = pd.to_numeric(AA_AF_df['AF'])
# The order is important, so GWAS will follow allele annotation as reference genome from 1000G phase3
GWAS_AA = wvg.merge_and_fix_SNPs_between_two_datasets(associated_snps, AA_AF_df, 'chr', 'pos', 'ref', 'alt', 'REF', 'ALT')
print("\n")
print("\n")
#print("Same alleles:")
print(GWAS_AA['Same_alleles'].value_counts())
if not all(GWAS_AA['Same_alleles']):
    GWAS_AA = GWAS_AA.loc[GWAS_AA["Same_alleles"] == True,:]
    print('Number of mismatching alleles filtered when merging 1000G is %d' % (GWAS_AA['Same_alleles'].value_counts()[False]))
print("\n")
print("\n")

# Remove columns prior to next merging
GWAS_AA = GWAS_AA.drop(['REF', 'ALT', 'REF_complementary', 'ALT_complementary', 'Same_alleles'], axis=1)
#print(GWAS_AA.iloc[0:4,0:6])


# 2. MERGE ARCHAIC WITH GWAS AND 1000G FOR ASSOCIATED SNPS
print("Merging archaic data with 1000G and GWAS by genomic coordinates...")
merged = wvg.merge_and_fix_SNPs_between_two_datasets(vcf_archaic, GWAS_AA, 'chr', 'pos', 'REF', 'ALT', 'ref', 'alt')


# 3. GET INDEPENDENT POSITIONS
print("\n")
print("\n")
print("Recovering independent positions")
res = map(lambda x: avg.get_SNP_row_from_physical_position(x, merged, check_alleles=False), list_prototypes)
res = list(res)
ids = []
failed_prot = []
for idx, el in enumerate(res):
    if len(el) == 0:
        failed_prot += [list_prototypes[idx]]
        continue
    else:
        ids += [el['snp_id']]

# The dictionary have alternative candidates that can be used in case that some prototypes failed:
if len(failed_prot) > 0:
    recovered_SNP = []
    for i in failed_prot:
        candidates = SNP_dict[i]
        #print(i)
        #print(candidates)
        if candidates != '':
            for c in range(0,len(candidates)):
                #print(candidates[c])
                row = avg.get_SNP_row_from_physical_position(candidates[c], merged, check_alleles=False)
                if len(row) != 0:
                    #print(row['SNP'])
                    recovered_SNP.append(row['snp_id'])
                    break

if len(recovered_SNP) > 0:
    # Add recovered SNPs to ids
    ids += recovered_SNP


# Get the final merged data
merged_indep = merged.loc[merged['snp_id'].isin(ids),:]
merged_indep = merged_indep.drop(['REF_complementary', 'ALT_complementary', 'Same_alleles'], axis=1)
merged_indep = merged_indep.rename(columns={'chr':'CHROM', 'pos':'POS'})
print("Number of independent positions with GWAS, archaic and 1000G data is: %d" %(merged_indep.shape[0]))

# Parse the genotype of the archaic
#geno_archaic = wvg.parse_genotype_vcf_one_sample(homo, merged_indep, 0)
#merged_indep.loc[:,homo] = geno_archaic


# 4. COMPARE IF SAME ALLELES
print("\n")
print("\n")
print("Checking that archaic alleles match 1000G and GWAS...")
same_as_ref_A = merged_indep['ALT'] == '.'
print('Alleles same as reference for %s:' % (homo))
print(same_as_ref_A.value_counts())

# Check different allele
df_polymorphism = merged_indep.loc[same_as_ref_A == False, :].copy()
# Complamentary annotation was already fixed if needed. You just need to compare two combinations
df_polymorphism, n_switched = wvg.check_if_switched_alleles_with_dataset(df_polymorphism, 'REF', 'ALT', 'ref', 'alt')
df_polymorphism, n_same = wvg.check_if_same_alleles_with_dataset(df_polymorphism, 'REF', 'ALT', 'ref', 'alt')
if (n_same + n_switched) == df_polymorphism.shape[0]:
    print('All archaic variants have the same alleles as GWAS')
else:
    print("Warning! not all the alleles match archaic data")
    filter = df_polymorphism.loc[df_polymorphism['Same_alleles'] == False,'snp_id'].tolist()
    merged_indep = merged_indep.loc[~merged_indep['snp_id'].isin(filter),:]
    print('Number of mismatching alleles filtered when merging %s is %d' % (homo, len(filter)))


# 5. DECODE VCF GENOTYE
merged_indep = avg.from_int_to_allele_pseudohaploid_VCF(merged_indep, homo)
print("\n")
print("\n")
#print(merged_indep.iloc[0:4,:])


# 6. FINALLY, CHECK FOR DIFFERENTIAL SELECTIVE PRESSURES BETWEEN HOMO AND AMH LINEAGES
#oddsratio, pval = avg.fisher_exact_test_proportion_risk_alleles_homo_lineage_vs_AMH(merged_indep, snp_id_col="snp_id", alternative=alternative)
print("All available options of the fisher exact test will be tested")
oddsratio, pval_two_sided = avg.fisher_exact_test_proportion_risk_alleles_homo_lineage_vs_AMH(
    merged_indep,
    snp_id_col="snp_id",
    alternative="two-sided"
)

oddsratio, pval_greater = avg.fisher_exact_test_proportion_risk_alleles_homo_lineage_vs_AMH(
    merged_indep,
    snp_id_col="snp_id",
    alternative="greater"
)

oddsratio, pval_less = avg.fisher_exact_test_proportion_risk_alleles_homo_lineage_vs_AMH(
    merged_indep,
    snp_id_col="snp_id",
    alternative="less"
)


print('Test of selective pressures between homo and AMH lineages performed using archaic %s. Processed file used is saved at %s' % (homo, out_path))
merged_indep.to_csv(out_path + 'SNPs_' + filename_out + '_' + str(len(SNP_dict.keys())) + '_indpMarkers_with_' + homo + '_archaic_and_1000G_data_risk_alleles.tsv', sep="\t", header=True, index=False)


# Add summarized results to logfile
command =  f"printf '{filename_out}\\t{oddsratio}\\t{pval_two_sided}\\t{pval_greater}\\t{pval_less}\\t{homo}\\t{merged_indep.shape[0]}\\t{len(list_prototypes)}\\n' >> {logfile}"
os.system(command)



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')