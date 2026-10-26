import sys
sys.path.append("/home/csic/gpo/olg/ccs/packages/vegeta_last_version")
import merge_GWAS_phenos as vg
import pandas as pd
import time
import argparse
import time


# ------------------ PARSE ARGUMENTS ------------------ 
# Create a parser object
parser = argparse.ArgumentParser(description="This function performs build conversion of GWAS summary statistics considering both hg19 and hg38 genomes")
# Add arguments
parser.add_argument("--file", type=str, required=True, help="File with GWAS data that has to be processed")
parser.add_argument("--dir", type=str, required=True, help="Path to the original file that will be used also to save outputs. Please finish it with a '/'")
parser.add_argument("--chr_col", type=str, required=True, help="Chromosome column name")
parser.add_argument("--pos_col", type=str, required=True, help="Position column name")
parser.add_argument("--original_build", type=str, required=True, help="The original build of the GWAS. If hg38, it will be converted to hg19 and the other way around.")
parser.add_argument("--directory_liftOver", type=str, required=True, help="The path to a directory that has liftOver exe and chainfiles for build conversion (hg19 to hg38 and the other way around).")


# Parse the arguments
args = parser.parse_args()

# Store arguments
file = args.file
dir_path = args.dir
chr_col = args.chr_col
pos_col = args.pos_col
build = args.original_build
path_liftOver = args.directory_liftOver



# Run the program for merging
# get the start time
st = time.time()



# Read and preprocess file. Format is automatically detected and if more than one delimiter is found, the file is preprocessed to fix the problem:
(file_format, list_nrows_1_field_when_whitespace) = vg.check_file_format(dir_path, file)
if file_format == "vcf" or file_format == "tsv":
    file_to_use = vg.check_more_than_one_delimiter_tsv(dir_path, file, file_format, list_nrows_1_field_when_whitespace)
else:
    file_to_use = file
        
# Get a pandas df
df = vg.file_to_pd_df(dir_path, file_to_use, file_format)

# Checkpoint for missing genomic coordinates
missing_pos = df[pos_col].isna()
if any(missing_pos):
    df = df.loc[missing_pos == False,:]
    print("In total %d variants have been removed due to missing genomic position" % (missing_pos.value_counts()[True]))
missing_chr = df[chr_col].isna()
if any(missing_chr):
    df = df.loc[missing_chr == False,:]
    print("In total %d variants have been removed due to missing chromosome" % (missing_chr.value_counts()[True]))

# Perform conversion: 
print("Starting build conversion...")
(bedfile, filename) = vg.create_bedfile_from_chr_pos_liftOver_input(df, file_to_use, dir_path, chr_col=chr_col, pos_col=pos_col)
(exit_status, map_file, unmap_file) = vg.run_liftOver_hg38_to_hg19_and_viceversa(bedfile, filename, dir_path, original_build=build, path_liftOver=path_liftOver)
if exit_status == 1:
    bedfile_path = dir_path + bedfile
    print("An error ocurred during the conversion of genomic coordinates from build %s" % (build))
    print("Try to run the liftOver command line tool with the bedfile found at %s" % (bedfile_path))
    print('File %s failed to be processed. Aborting' % (file_to_use))
    print("***********************************************************************************")
    print("\n\n")
else:
    (df_conversion, missing_coordinates_counter) = vg.add_coordinates_from_liftOver_conversion(df, map_file, unmap_file, file_to_use, dir_path, original_build=build, chr_col=chr_col, pos_col=pos_col, save_tsv_liftOver=True)



# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
final_res = elapsed_time / 60
print('Execution time:', final_res, 'minutes')
print('Execution time:', elapsed_time, 'seconds')


# Write the results
#merged_df.to_csv(path_or_buf=dir_path + 'Merged_GWAS_meta_analyses_Diet_Cluster_1_Omega3_prot_outer_pval_1e-05.tsv.gz', header=True, index=False, decimal='.', sep='\t', compression={'method': 'gzip', 'compresslevel': 4})