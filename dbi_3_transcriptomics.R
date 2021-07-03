################################################################################
#
# Transcriptomic analysis
#
################################################################################

# Format R session -------------------------------------------------------------
# Load packages
library(biomaRt)
library(cobs)
library(doParallel)
library(edgeR)
library(EMMREML)
library(limma)
library(tidyverse)
select <- dplyr::select

# Set theme and colors for graphing
theme_set(theme_classic(base_size = 20, base_family = ""))  # graphing defaults
wes_col <- "#D55E00"  # specify color for Western diet
med_col <- "#56B4E9"  # specify color for Mediterranean diet

# Import raw data
load("data_files/processed_data_2.RData")

# PCA on Gene Expression -------------------------------------------------------
# Conduct principal component analysis on the correlation matrix of normalized
# residual gene expression.

# PCA
GE_PCA <- prcomp(cor(r_matrix))
sum_GE_PCA <- summary(GE_PCA)  # store PC importance information
sample_info$ge_pc1 <- GE_PCA$x[, 1]  # add PC projection to Sample Info matrix

# T test comparing PC1 projection between diet groups
t.test(sample_info[sample_info$diet == "MD", "ge_pc1"],
       sample_info[sample_info$diet == "WD", "ge_pc1"])

# Generate Figure 1A
figure_1a <- sample_info %>%
  ggplot(aes(x = diet, y = ge_pc1,
             fill = diet, color = diet)) +
  geom_hline(yintercept = 0, size = 1, linetype = "dotted") +
  geom_boxplot(outlier.alpha = 0, col = "black", alpha = 0.3,
               coef = 0) +
  geom_jitter(size = 2, width = 0.1) +
  labs(y = "Projection onto Gene Expression PC1") +
  scale_fill_manual(name = "Diet",
                    breaks = c("WD", "MD"),
                    values = c(wes_col, med_col),
                    labels = c("WD", "MD")) +
  scale_color_manual(name = "Diet",
                     breaks = c("WD", "MD"),
                     values = c(wes_col, med_col),
                     labels = c("WD", "MD")) +
  scale_x_discrete(labels = c("WD" = "Western", "MD" = "Mediterranean")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")
saveRDS(figure_1a, file = "plots/figure_1a.RDS")

# Gene expression (GE) ~ diet --------------------------------------------------
# Model normalized residual gene expression as a function of diet using linear
# mixed effects model.

# Create Z matrix
Z_matrix <- diag(1, nrow = 35)  # identity matrix

# Set the parameters for parallel computing for lmm using EMMREML
ncores <- detectCores(logical = TRUE)  # find number of available cores
clus <- makeCluster(ncores, setup_strategy = "sequential")  # create cluster
registerDoParallel(cores = ncores)
clusterExport(clus,
              varlist = c("r_matrix", "sample_info", "Z_matrix", "kinship"),
              envir = environment())  # initiate local cluster

# Model residual gene expression as a function of diet in genewise models
emma_GE_diet <- data.frame(t(parApply(clus, r_matrix, 1, function(y){
  library(EMMREML)
  emma <- emmreml(y = y,  # model each gene
                  X = model.matrix(~ sample_info$diet_west),  # design matrix
                  Z = as.matrix(Z_matrix),  # identity matrix
                  K = as.matrix(kinship),  # relatedness matrix
                  varbetahat = T, varuhat = T, PEVuhat = T, test = T)
  return(c(emma$betahat[2], emma$pvalbeta[2, "none"], emma$varbetahat[2]))  
})))
colnames(emma_GE_diet) <- c("es_diet", "diet_pval", "diet_var")  # rename columns

# Calculate standardized betas
emma_GE_diet$std_beta_diet <- emma_GE_diet$es_diet / sqrt(emma_GE_diet$diet_var)

# Empirical False Discovery Rate (FDR) -----------------------------------------
# Model the effect of diet (or DAB) on gene expression in 1000 permutations,
# permuting the meta data in order to determine a null distribtuion of p-values.
# Then, use the permuted null distribution to calculate an empirical FDR,
# (Snyder-Mackler et al., 2016). Due to computational demands, this is best
# completed on a high performance computing cluster.


# Create a biomart object with ensembl mfas genes
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mfascicularis_gene_ensembl")

# Create a biomart object of GO terms associated with genes from study
mfas_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = unlist(emma_GE_diet$stable_ID),
                    mart = mart)

# Create data frame combining diet and DAB models
emma_GE_diet$stable_ID <- row.names(emma_GE_diet)
emma_GE_diet <- as_tibble(emma_GE_diet) %>%
  arrange(stable_ID)

# Add gene names to DEG_info table
DEG_info <- emma_GE_diet %>%
  left_join(mfas_genes, by = c("stable_ID" = "ensembl_gene_id"))

# Save output data -------------------------------------------------------------
# Save data for downstream analyses
save(gene_counts, sample_info, r_matrix, kinship, DEG_info, emma_GE_diet,
     file = "data_files/processed_data_3.RData")
