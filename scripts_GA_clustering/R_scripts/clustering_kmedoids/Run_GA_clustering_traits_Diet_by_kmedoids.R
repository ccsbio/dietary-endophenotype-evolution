# Libraries
library(doParallel)
library(doRNG)
library(dplyr)
library(tidyr)
library(stringr)

# Import GA functions
source("/home/csic/gpo/olg/ccs/Diet_project/scripts_GA_clustering/R_scripts/clustering_kmedoids/Functions_GA_clustering_traits_by_kmedoids.R", keep.source = FALSE)

# Get command line arguments with the parameters to run the GA  algorithm
# Mandatory parameters
args <- commandArgs(trailingOnly = TRUE)
population_size <- as.numeric(args[1])
mutation_rate <- as.numeric(args[2])
iterations <- as.numeric(args[3])
rep_type <- as.character(args[4])
k_clusters <- as.numeric(args[5])
file_MDS_from_LDSC <- as.character(args[6])
out_dir <- as.character(args[7])
no_cores <- as.numeric(args[8])



# Optional argument `tournament_size`
if (length(args) == 9) {
  tournament_size <- as.numeric(args[9])
} else {
  # Default value if not provided
  tournament_size <- 0
}


# Prints
print("The following parameters are considered when running the GA:")
print(paste("Number of cores:", no_cores, sep = " "))
print(paste("Pop size:", population_size, sep = " "))
print(paste("Mutation rate:", mutation_rate, sep = " "))
print(paste("Generations:", iterations, sep = " "))
print(paste("Reproduction type:", rep_type, sep = " "))
print(paste("MDS input file:", file_MDS_from_LDSC, sep = " "))
print(paste("Tournament size:", tournament_size, sep = " "))


# Check output dir is ok
if (substring(out_dir, nchar(out_dir)) != "/"){
    stop("ERROR: out dir not valid. Please be sure to provide a directory with '/' in the end of the string")
  }


# ---- READ MATRICES ----
# Read fitted MDS for the pairwise genetic correlation matrix from LDSC for a set of GWAS
ldsc_MDS <- read.delim(file_MDS_from_LDSC, sep = "\t", header = TRUE, check.names = FALSE)


# Get Euclidean distances using the first two dimensions and store GWAS labels
coords <- ldsc_MDS[, c("D1", "D2")]
rownames(coords) <- ldsc_MDS$Codes  # assuming Codes column holds the names
dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
all_GWAS <- colnames(dist_matrix)
# Prints
print("Warning! input file must contain a 'Codes' column with the GWAS names. So, each row with MDS data (columns as D1, D2, ...) can be tracked to a particular datapoint")
print(paste("In total", as.character(length(all_GWAS)), "number of traits have been found", sep = " "))



# ---- PARALLEL SETTINGS ----
# Setup cluster: good practice not to use all the cores
clust <- makeCluster(no_cores, type = "PSOCK") #This line will take time

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = clust)

# (optional)
#check if it is registered 
foreach::getDoParRegistered()
#how many workers are available? 
foreach::getDoParWorkers()



# ---- RUN PARALLELIZED GA ----
results <- evolve.or.die.clusters.parallel.kmedoids(all_GWAS, dist_matrix, population_size, mutation_rate, iterations, k_clusters, optimization_by="Maximization", type_of_reproduction = rep_type, n_tournament = tournament_size)



# Stop the cluster settings used during parallelization
stopCluster(cl = clust)


# ---- EXPORT GA BEST RESULT ----
print(paste("Saving at:", out_dir, sep = " "))
save(results, file = paste(out_dir, "Best_clustering_GA_Diet_", as.character(length(all_GWAS)), "traits_", as.character(k_clusters), "clust_", rep_type, "rep_", iterations, "itera_", population_size, "pop_", mutation_rate, "MR.rda", sep = ""))