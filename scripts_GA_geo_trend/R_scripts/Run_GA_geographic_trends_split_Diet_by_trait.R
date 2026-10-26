# Libraries
library(doParallel)
library(doRNG)
library(dplyr)
library(tidyr)

# Import GA functions
source("/home/csic/gpo/olg/ccs/Diet_project/scripts_GA_geo_trend/R_scripts/Functions_GA_geographic_trend.R", keep.source = FALSE)

# Get command line arguments with the parameters to run the GA  algorithm
# Mandatory parameters
args <- commandArgs(trailingOnly = TRUE)
population_size <- as.numeric(args[1])
mutation_rate <- as.numeric(args[2])
iterations <- as.numeric(args[3])
rep_type <- as.character(args[4])
df_samples_metadata_file <- as.character(args[5])
x_variable_name <- as.character(args[6])
y_variable_name <- as.character(args[7])
longitud_col <- as.character(args[8])
latitud_col <- as.character(args[9])
out_dir <- as.character(args[10])
trait <- as.character(args[11])
no_cores <- as.numeric(args[12])

# Optional argument `tournament_size`
if (length(args) == 13) {
  tournament_size <- as.numeric(args[13])
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
print(paste("Sample data file:", df_samples_metadata_file, sep = " "))
print(paste("Independent variable:", x_variable_name, sep = " "))
print(paste("Dependent variable:", y_variable_name, sep = " "))
print(paste("Latitud column:", latitud_col, sep = " "))
print(paste("Longitud column:", longitud_col, sep = " "))
print(paste("Trait:", trait, sep = " "))
print(paste("Tournament size:", tournament_size, sep = " "))


# Check output dir is ok
if (substring(out_dir, nchar(out_dir)) != "/"){
    stop("ERROR: out dir not valid. Please be sure to provide a directory with '/' in the end of the string")
  }


# ---- READ MATRICES ----
# Read data with GWAS summary statistics
df_samples <- read.delim(df_samples_metadata_file, sep = "\t", header = TRUE)
df_trait <- df_samples %>% filter(Trait == trait)
# Store unique values of geographic coordinates in the whole population
long <- unique(df_trait[, longitud_col])
lat <- unique(df_trait[, latitud_col])
print(length(long))
print(length(lat))


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

# Export packages to use
clusterEvalQ(clust, {
  library(dplyr)
  library(tidyr)
})


# ---- RUN PARALLELIZED GA ----
results <- evolve_or_die_geo_trends_across_time(df_trait, x_variable_name, y_variable_name, long, lat, population_size, mutation_rate, iterations, rep_type, tournament_size)

# Stop the cluster settings used during parallelization
stopCluster(cl = clust)


# ---- EXPORT GA BEST RESULT ----
print(paste("Saving at:", out_dir, sep = " "))
save(results, file = paste(out_dir, "Best_geographic_trend_differences_GA_", trait, "_", rep_type, "rep_", iterations, "itera_", population_size, "pop_", mutation_rate, "MR_", x_variable_name, "Xvar.rda", sep = ""))