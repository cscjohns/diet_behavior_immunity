#!/usr/bin/env Rscript

# Specify path to R package library
.libPaths(c("/gscratch/csde/smacklab/corbin/Rlibs",.libPaths()))

load("/gscratch/csde/smacklab/corbin/pval_distribution/pvals_1_out.Rdata")

diet_pval_distribution <- as.vector(perm_pvals)

load("/gscratch/csde/smacklab/corbin/pval_distribution/pvals_DAB_1_out.Rdta")

DAB_pval_distribution <- as.vector(perm_pvals)

# Load in the R data of permuted data
for(i in 2:100){
  load(paste0("/gscratch/csde/smacklab/corbin/pval_distribution/pvals_", i, "_out.Rdata"))
  diet_pval_distribution <- c(diet_pval_distribution, as.vector(perm_pvals))

  load(paste0("/gscratch/csde/smacklab/corbin/pval_distribution/pvals_DAB_", i, "_out.Rdata"))
  DAB_pval_distribution <- c(DAB_pval_distribution, as.vector(perm_pvals))
}

# Load in original data
load("pval_perms_in.RData")

emma_diet$diet_FDR <- rep(1, times = nrow(emma_diet))
emma_DAB$DAB_FDR <- rep(1, times = nrow(emma_DAB))

for(i in 1:nrow(emma_diet){
  emma_diet[i, "diet_FDR"] <- sum(emma_diet[i, "diet_pval"] > diet_pval_distribution)/length(diet_pval_distribution)
  emma_DAB[i, "DAB_FDR"] <- sum(emma_DAB[i, "DAB_pval"] > DAB_pval_distribution)/length(DAB_pval_distribution)
}

save(emma_diet, emma_DAB, file = "emma_diet_DAB_out.Rdata")

q(save = "no")
