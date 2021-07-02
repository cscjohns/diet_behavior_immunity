#!/usr/bin/env Rscript

# Specify path to R package library
.libPaths(c("/gscratch/csde/smacklab/corbin/Rlibs",.libPaths()))

# Load packages needed
library(parallel)
library(doParallel)

# Store command line arguments
in_args <- commandArgs(trailingOnly = TRUE)

# Will complete 10 permutations of the data, getting a p-value distribution
# of the effect sizes of diet to compare against the true p-value distribution
# of effect sizes of diet

# Load in the R data
load("/gscratch/csde/smacklab/corbin/pval_distribution/pval_perms_in.Rdata")

# Create Z matrix
Z_matrix <- diag(1, nrow = 35)

# Set the parameters
nReps <- 10  # set the number of permutations
perm_pvals <- data.frame(matrix(rep(0, times = nReps * nrow(r_matrix)),
                                nrow = nrow(r_matrix)))
ncores <- detectCores(logical = TRUE)

# Loop through the permutations, resample sample IDs each time
for (i in 1:nReps) {
  # permute the original sample information data frame
  set.seed(10 * as.numeric(in_args[1]) + i)
  perm_ids <- sample(1:35, 35, replace = FALSE)
  perm_meta <- sample_info[perm_ids, ]
  clus <- makeCluster(ncores)
  registerDoParallel(cores = ncores)
  clusterExport(clus,
                varlist = c("r_matrix", "Z_matrix", "kinship",
                            "perm_meta"),
                envir = environment())
  lmm_perm <- parApply(clus, r_matrix, 1, function(y){
    library(EMMREML)
    emma <- emmreml(y = y,
                    X = as.matrix(model.matrix(~ perm_meta[, "DAB"])),
                    Z = as.matrix(Z_matrix),
                    K = as.matrix(kinship),
                    varbetahat = T, varuhat = T, PEVuhat = T, test = T)
    p <- emma$pvalbeta[2, "none"]
    return(p)
  })
  perm_pvals[, i] <- lmm_perm
  print(i)
}


save(perm_pvals, file = paste0("/gscratch/csde/smacklab/corbin/pval_distribution/pvals_DAB_", in_args[1], "_out.Rdata"))

q(save = "no")
