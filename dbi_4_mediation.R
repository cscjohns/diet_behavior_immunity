################################################################################
#
# Mediation analyses
#
################################################################################

# Format R session -------------------------------------------------------------
# Load packages
library(data.table)
library(doParallel)
library(EMMREML)
library(limma)
library(parallel)
library(tidyverse)
select <- dplyr::select

# Set theme and colors for graphing
theme_set(theme_classic(base_size = 20, base_family = ""))  # graphing defaults
wes_col <- "#D55E00"  # specify color for Western diet
med_col <- "#56B4E9"  # specify color for Mediterranean diet

# Import raw data
load("data_files/processed_data_3.RData")

# Gene expression (GE) ~ DAB score ---------------------------------------------
# Model gene expression as a function of DAB score using linear mixed effects
# models.

# Set the parameters for parallel computing for lmm using EMMREML
ncores <- detectCores(logical = TRUE)  # find number of available cores
clus <- makeCluster(ncores, setup_strategy = "sequential")  # create cluster
registerDoParallel(cores = ncores)
clusterExport(clus,
              varlist = c("r_matrix", "sample_info", "Z_matrix", "kinship"),
              envir = environment())  # initiate local cluster

# Model residual gene expression as a function of DAB in genewise models
emma_GE_dab <- data.frame(t(parApply(clus, r_matrix, 1, function(y){
  library(EMMREML)
  emma <- emmreml(y = y,  # model each gene
                  X = model.matrix(~ sample_info$DAB),  # design matrix
                  Z = as.matrix(Z_matrix),  # identity matrix
                  K = as.matrix(kinship),  # relatedness matrix
                  varbetahat = T, varuhat = T, PEVuhat = T, test = T)
  return(c(emma$betahat[2], emma$pvalbeta[2, "none"], emma$varbetahat[2]))
})))
colnames(emma_GE_dab) <- c("es_DAB", "DAB_pval", "DAB_var")  # rename

# Calculate standardized betas
emma_GE_dab$std_beta_DAB <- emma_GE_dab$es_DAB / sqrt(emma_GE_dab$DAB_var)

# Calculate FDR using q-value function (Storey & Tibshirani, 2003)
emma_GE_dab$DAB_FDR <- qvalue(emma_GE_dab$DAB_pval)$qvalues


# Gene expression (GE) ~ diet + DAB score --------------------------------------
# Model gene expression as a function of diet and DAB score using linear mixed
# effects models.

# Set the parameters for parallel computing for lmm using EMMREML
registerDoParallel(cores = ncores)
clusterExport(clus,
              varlist = c("r_matrix", "sample_info", "Z_matrix", "kinship"),
              envir = environment())  # initiate local cluster

# Model residual gene expression as a function of diet and DAB score in genewise
# models
emma_GE_diet_dab <- data.frame(t(parApply(clus, r_matrix, 1, function(y){
  library(EMMREML)
  emma <- emmreml(y = y,  # model each gene
                  X = model.matrix(~ sample_info$diet_west +
                                     sample_info$DAB),  # design matrix
                  Z = as.matrix(Z_matrix),  # identity matrix
                  K = as.matrix(kinship),  # relatedness matrix
                  varbetahat = T, varuhat = T, PEVuhat = T, test = T)
  return(c(emma$betahat[2], emma$pvalbeta[2, "none"], emma$varbetahat[2]))
})))
colnames(emma_GE_diet_dab) <- c("es_diet", "diet_pval", "diet_var")  # rename

# Calculate standardized betas
emma_GE_diet_dab$std_beta_diet <- emma_GE_diet_dab$es_diet /
  sqrt(emma_GE_diet_dab$diet_var)

# Calculate FDR using q-value function (Storey & Tibshirani, 2003)
emma_GE_diet_dab$diet_FDR <- qvalue(emma_GE_diet_dab$diet_pval)$qvalues

# DAB ~ GE ---------------------------------------------------------------------
# Model DAB score as a function of gene expression using linear mixed effects
# models.

# Set the parameters for parallel computing for lmm using EMMREML
registerDoParallel(cores = ncores)
clusterExport(clus,
              varlist = c("r_matrix", "sample_info", "Z_matrix", "kinship"),
              envir = environment())  # initiate local cluster

# Model DAB score as a funciton of residual gene expression in genewise models
emma_DAB_ge <- data.frame(t(parApply(clus, r_matrix, 1, function(x){
  library(EMMREML)
  emma <- emmreml(y = sample_info$DAB,  # model each gene
                  X = model.matrix(~ x),  # design matrix
                  Z = as.matrix(Z_matrix),  # identity matrix
                  K = as.matrix(kinship),  # relatedness matrix
                  varbetahat = T, varuhat = T, PEVuhat = T, test = T)
  return(c(emma$betahat[2], emma$pvalbeta[2, "none"], emma$varbetahat[2]))
})))
colnames(emma_DAB_ge) <- c("es_GE", "GE_pval", "GE_var")  # rename

# Calculate standardized betas
emma_DAB_ge$std_beta_GE <- emma_DAB_ge$es_GE / sqrt(emma_DAB_ge$GE_var)

# Calculate FDR using q-value function (Storey & Tibshirani, 2003)
emma_DAB_ge$GE_FDR <- qvalue(emma_DAB_ge$GE_pval)$qvalues

# DAB ~ diet -------------------------------------------------------------------
# Model DAB score as a function of diet using a linear mixed effects model.

library(EMMREML)
emma_DAB_diet <- emmreml(y = sample_info$DAB,
                         X = model.matrix(~ sample_info$diet_west),  # design
                         Z = as.matrix(Z_matrix),  # identity matrix
                         K = as.matrix(kinship),  # relatedness matrix
                         varbetahat = T, varuhat = T, PEVuhat = T, test = T)

# Calculate standardized beta
emma_DAB_diet$std_beta_diet <- emma_DAB_diet$betahat[2] /
  sqrt(emma_DAB_diet$varbetahat[2])


# DAB ~ diet + GE --------------------------------------------------------------
# Model DAB score as a function of diet and gene expression using linear mixed
# effects models.

# Set the parameters for parallel computing for lmm using EMMREML
registerDoParallel(cores = ncores)
clusterExport(clus,
              varlist = c("r_matrix", "sample_info", "Z_matrix", "kinship"),
              envir = environment())  # initiate local cluster

# Model DAB as a function of diet and residual gene expression
emma_DAB_diet_GE <- data.frame(t(parApply(clus, r_matrix, 1, function(x){
  library(EMMREML)
  emma <- emmreml(y = sample_info$DAB,  # model each gene
                  X = model.matrix(~ sample_info$diet_west +
                                     r_matrix[1,]),  # design matrix
                  Z = as.matrix(Z_matrix),  # identity matrix
                  K = as.matrix(kinship),  # relatedness matrix
                  varbetahat = T, varuhat = T, PEVuhat = T, test = T)
  return(c(emma$betahat[2], emma$pvalbeta[2, "none"], emma$varbetahat[2]))
})))
colnames(emma_DAB_diet_GE) <- c("es_diet", "diet_pval", "diet_var")  # rename

# Calculate standardized betas
emma_DAB_diet_GE$std_beta_diet <- emma_DAB_diet_GE$es_diet /
  sqrt(emma_DAB_diet_GE$diet_var)

# Calculate FDR using q-value function (Storey & Tibshirani, 2003)
emma_DAB_diet_GE$diet_FDR <- qvalue(emma_DAB_diet_GE$diet_pval)$qvalues

# Bootsrap effect sizes and errors----------------------------------------------
# This section of code was run as an array, calling the script below 1000 times
# with the following call in the submission job:

# taskID=${SLURM_ARRAY_TASK_ID}
# 
# sbatch -p ckpt -A csde-ckpt --mem=120G --time=20:00 mediation.R $taskID 10

# so that the job number is transmitted to each job and can be used to name the
# output. This will bootstrap the data 10000 times (1000 jobs * 10 reps per job)
# to calculate a confidence interval of the effect sizes of diet and DAB.

# Store command line arguments
in_args <- commandArgs(trailingOnly = TRUE)

# Create Z matrix
Z_matrix <- diag(1, nrow = 35)

# Set the parameters
nReps <- as.numeric(in_args[2])  # set the number of bootstrap permutations
ncores <- detectCores(logical = TRUE)

# Will store the betas and standardized betas for each of the following models:
# m_GE_diet: GE ~ diet
# m_GE_diet_dab: GE ~ diet + DAB
# m_DAB_diet: DAB ~ diet
# m_DAB_diet_ge: DAB ~ diet + GE

# Create empty data frames to hold betas and se from each model
m_GE_diet_beta <- data.frame(matrix(rep(0, times = nReps * nrow(r_matrix)),
                                 nrow = nrow(r_matrix)))
m_GE_diet_dab_beta <- m_GE_diet_beta
m_DAB_diet_beta <- m_GE_diet_beta
m_DAB_diet_ge_beta <- m_DAB_diet_beta
m_GE_diet_se <- m_GE_diet_beta
m_GE_diet_dab_se <- m_GE_diet_beta
m_DAB_diet_se <- m_DAB_diet_beta
m_DAB_diet_ge_se <- m_DAB_diet_beta

# Loop through the permutations, resample sample IDs each time
for (i in 1:nReps) {
  # set a unique seed value
  set.seed(nReps * (as.numeric(in_args[1]) - 1) + i)
  
  # permute the original sample information data frames
  bs_ids <- sample(1:35, 35, replace = TRUE)
  bs_kin <- kinship[bs_ids, bs_ids]
  bs_meta <- sample_info[bs_ids, ]
  bs_expr <- r_matrix[, bs_ids]
  bs_Zmat <- Z_matrix[bs_ids, bs_ids]
  
  # run models on a local cluster of nodes
  clus <- makeCluster(ncores)
  registerDoParallel(cores = ncores)
  clusterExport(clus,
                varlist = c("bs_expr", "bs_meta", "bs_Zmat", "bs_kin"),
                envir = environment())
  
  # m_GE_diet: GE ~ diet
  m_GE_diet <- data.frame(t(parApply(clus, bs_expr, 1, function(y){
    library(EMMREML)
    emma <- emmreml(y = y,
                    X = as.matrix(model.matrix(~ bs_meta$diet)),
                    Z = as.matrix(bs_Zmat),
                    K = as.matrix(bs_kin),
                    varbetahat = T)
    return(c(emma$betahat[2], emma$varbetahat[2]))})))
  m_GE_diet_beta[, i] <- m_GE_diet[, 1]
  m_GE_diet_se[, i] <- m_GE_diet[, 2]

  # m_GE_diet_dab: GE ~ diet + DAB
  m_GE_diet_dab <- data.frame(t(parApply(clus, bs_expr, 1, function(y){
    library(EMMREML)
    emma <- emmreml(y = y,
                    X = as.matrix(model.matrix(~ bs_meta$diet + bs_meta$DAB)),
                    Z = as.matrix(bs_Zmat),
                    K = as.matrix(bs_kin),
                    varbetahat = T)
    return(c(emma$betahat[2], emma$varbetahat[2]))})))
  m_GE_diet_dab_beta[, i] <- m_GE_diet_dab[, 1]
  m_GE_diet_dab_se[, i] <- m_GE_diet_dab[, 2]
  
  # m_DAB_diet: DAB ~ diet
  library(EMMREML)
  m_DAB_diet <- emmreml(y = bs_meta$DAB,
                    X = as.matrix(model.matrix(~ bs_meta$diet)),
                    Z = as.matrix(bs_Zmat),
                    K = as.matrix(bs_kin),
                    varbetahat = T)
  m_DAB_diet_beta[, i] <-
    rep(m_DAB_diet$betahat[2], times = nrow(m_DAB_diet_beta))
  m_DAB_diet_se[, i] <-
    rep(m_DAB_diet$varbetahat[2], times = nrow(m_DAB_diet_se))
  
  # m_DAB_diet_ge: DAB ~ diet + GE
  m_DAB_diet_ge <- data.frame(t(parApply(clus, bs_expr, 1, function(x){
    library(EMMREML)
    emma <- emmreml(y = bs_meta$DAB,
                    X = as.matrix(model.matrix(~ bs_meta$diet + x)),
                    Z = as.matrix(bs_Zmat),
                    K = as.matrix(bs_kin),
                    varbetahat = T)
    return(c(emma$betahat[2], emma$varbetahat[2]))})))
  m_DAB_diet_ge_beta[, i] <- m_DAB_diet_ge[, 1]
  m_DAB_diet_ge_se[, i] <- m_DAB_diet_ge[, 2]
  
  print(i)  # print iteration # for troubleshooting
}

# Create data frames of standardized betas
m_GE_diet_std_beta <- m_GE_diet_beta / sqrt(m_GE_diet_se)
m_GE_diet_dab_std_beta <- m_GE_diet_dab_beta / sqrt(m_GE_diet_dab_se)
m_DAB_diet_std_beta <- m_DAB_diet_beta / sqrt(m_DAB_diet_se)
m_DAB_diet_ge_std_beta <- m_DAB_diet_ge_beta / sqrt(m_DAB_diet_dab_se)
print("std_beta dfs done")

# Create data frames of differences in standardized betas
diff_GE_base_med_std_beta <- m_GE_diet_std_beta - m_GE_diet_dab_std_beta
diff_DAB_base_med_std_beta <- m_DAB_diet_std_beta - m_DAB_diet_ge_std_beta
print("beta differences calculated")

# Write data frames to tsv files
write_tsv(as.data.frame(t(diff_GE_base_med_std_beta)),
          paste0("/path_to_directory/diff_GE_base_med_std_beta_",
                 formatC(as.numeric(in_args[1]), width = 4, flag = 0), ".tsv"),
          col_names = FALSE)

write_tsv(as.data.frame(t(diff_DAB_base_med_std_beta)),
          paste0("/path_to_directory/diff_DAB_base_med_std_beta_",
                 formatC(as.numeric(in_args[1]), width = 4, flag = 0), ".tsv"),
          col_names = FALSE)

# Collate outputs---------------------------------------------------------------
# The outputs from the 1000 jobs were collated through running the following
# lines of code through a submitted script in command line.

# cat diff_GE_base_med_std_beta*.tsv > diff_GE_base_med_std_beta_out.tsv
# cat diff_DAB_base_med_std_beta*.tsv > diff_DAB_base_med_std_beta_out.tsv

# Calculate confidence intervals------------------------------------------------
# The following lines were run on the HPC due to the large data frames from the
# boostrapped iterations (12240 x 10000).

# Load data
diff_GE_base_med_std_beta <- t(fread(
  "/path_to_directory/diff_GE_base_med_std_beta_out.tsv"))
diff_DAB_base_med_std_beta <- t(fread(
  "/path_to_directory/diff_DAB_base_med_std_beta_out.tsv"))

print("files loaded")

# Create data frames of confidence intervals
diff_GE_base_med_std_beta_CI <- tibble(
  "stable_ID" = DEG_info$stable_ID,
  "lowCI" = apply(diff_GE_base_med_std_beta, 1, function(x){
    quantile(x, 0.05)}),
  "highCI" = apply(diff_GE_base_med_std_beta, 1, function(x){
    quantile(x, 0.95)}))
diff_DAB_base_med_std_beta_CI <- tibble(
  "stable_ID" = DEG_info$stable_ID,
  "lowCI" = apply(diff_DAB_base_med_std_beta, 1, function(x){
    quantile(x, 0.05)}),
  "highCI" = apply(diff_DAB_base_med_std_beta, 1, function(x){
    quantile(x, 0.95)}))

# Save data
save(diff_GE_base_med_std_beta_CI, diff_DAB_base_med_std_beta_CI,
     file = "/data_files/mediation_CI.RData")

# Mediation plot-------------
load("data_files/mediation_CI.RData")

mediation_table <- DEG_info %>%
  select(std_beta_diet, diet_FDR, DAB_FDR, stable_ID, external_gene_name) %>%
  left_join(diff_GE_base_med_std_beta_CI,
            by = "stable_ID")

mediation_slim <- mediation_table %>%
  filter(diet_FDR < 0.05,
         DAB_FDR < 0.05,
         lowCI * highCI > 0)

mediation_slim <- mediation_slim %>%
  mutate(percent_mediated = (std_beta_diet - mediated_std_beta_diet)/std_beta_diet)

mediation_slim %>%
  ggplot(aes(x = percent_mediated)) +
  geom_histogram(bins = 100)

mediation_slim %>%
  ggplot(aes(x = 100*percent_mediated,
             fill = factor(ifelse(std_beta_diet > 0, 1, 0)))) +
  geom_histogram(alpha = 0.5,
                 position = "identity",
                 bins = 60) +
  scale_fill_manual(breaks = c(0, 1),
                    values = c(wes_col, med_col),
                    labels = c("Western (n = 712 genes)","Mediterranean (n = 487 genes)"),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Mediation of Diet Effect on GE by DAB (%)",
       y = "Number of Genes") +
  theme(legend.position = "top")


# GE Mediation of DAB-----------------------------------------------------------
load("../2_data/2021_mediation_rev_quantiles.Rdata")

mediation_table2 <- DEG_info %>%
  select(std_beta_diet, diet_FDR, stable_ID, external_gene_name) %>%
  mutate(es_GE_diet = std_beta_diet)

mediation_table2$std_beta_diet <-
  rep(emma_DAB_diet$betahat[2]/sqrt(emma_DAB_diet$varbetahat[2]), times = 12240)
mediation_table2$GE_pval <- emma_DAB_ge$GE_pval
mediation_table2$GE_FDR <- emma_DAB_ge$GE_FDR
mediation_table2$mediated_std_beta_diet <- emma_DAB_diet_GE$std_beta_diet
mediation_table2$low95 <- diff_DAB_base_std_beta_quantiles$X2.5.
mediation_table2$hi95 <- diff_DAB_base_std_beta_quantiles$X97.5.
mediation_table2$low90 <- diff_DAB_base_std_beta_quantiles$X5.
mediation_table2$hi90 <- diff_DAB_base_std_beta_quantiles$X95.

mediation_slim2 <- mediation_table2 %>%
  filter(diet_FDR < 0.05,
         GE_FDR < 0.05,
         low90 * hi90 > 0)

mediation_slim2 <- mediation_slim2 %>%
  mutate(percent_mediated = (std_beta_diet - mediated_std_beta_diet)/std_beta_diet)

mediation_slim2 %>%
  ggplot(aes(x = percent_mediated)) +
  geom_histogram(bins = 100)

mediation_slim2 %>%
  ggplot(aes(x = 100*percent_mediated,
             fill = factor(ifelse(es_GE_diet > 0, 1, 0)))) +
  geom_histogram(alpha = 0.5,
                 position = "identity",
                 bins = 60) +
  scale_fill_manual(breaks = c(0, 1),
                    values = c(wes_col, med_col),
                    labels = c("Western (n = 523 genes)","Mediterranean (n = 375 genes)"),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Mediation of Diet Effect on DAB by GE (%)",
       y = "Number of Genes") +
  theme(legend.position = "top")

mediation_slim %>%
  mutate(p_med_1 = 100 * percent_mediated) %>%
  select(stable_ID, p_med_1) %>%
  inner_join(mediation_slim2 %>%
               mutate(p_med_2 = 100 * percent_mediated) %>%
               select(stable_ID, p_med_2),
             by = "stable_ID") %>%
  ggplot(aes(x = p_med_1, y = p_med_2)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "P1: DAB mediates GE ~ diet",
       y = "P2: GE mediates DAB ~ diet")
