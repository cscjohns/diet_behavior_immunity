################################################################################
#
# Prepare transcriptomic data for downstream analyses
#
################################################################################

# Format R session -------------------------------------------------------------
# Load packages
library(edgeR)
library(limma)
library(tidyverse)

# Set theme and colors for graphing
theme_set(theme_classic(base_size = 20, base_family = ""))  # graphing defaults
wes_col <- "#D55E00"  # specify color for Western diet
med_col <- "#56B4E9"  # specify color for Mediterranean diet

# Import raw data
load("data_files/raw_data.RData")

# Examine CD14 purity ----------------------------------------------------------
# Use CD14 (monocyte) and CD3 (T Cell) expression as markers of CD14 purity
# following bead-based selection. Generate plot of marker gene expression.

# Calculate marker gene expression, add to sample info matrix
gene_rpkm <- rpkm(gene_counts,
                  gene.length = gene_lengths,
                  normalized.lib.sizes = TRUE)
CD3 <- c('ENSMFAG00000040063', 'ENSMFAG00000045644', 'ENSMFAG00000040556',
         'ENSMFAG00000000076')  # genes associated w/ CD3 (Tcell marker)
CD14 <- 'ENSMFAG00000025943'  # gene associated w/ CD14 (monocyte marker)
sample_info$cd3r <- as.numeric(apply(gene_rpkm[CD3, ], 2, sum))
sample_info$cd14r <- as.numeric(gene_rpkm[CD14, ])

# Generate Appendix Figure 3
outliers <- c('m8361', 'm8365', 'm8385')  # outliers determined by CD3, CD14
appendix_figure_3 <- sample_info %>%
  mutate(outlier = ifelse(id %in% outliers, "outliers", "")) %>%
  ggplot(aes(x = cd3r, y = cd14r)) +
  geom_point(aes(col = diet), size = 3) +
  scale_color_manual(breaks = c("WD", "MD"),
                     values = c(wes_col, med_col),
                     labels = c("Western", "Mediterranean"),
                     name = NULL) +
  labs(x = "T Cell Marker (CD3)",
       y = "PBMC Marker (CD14)") +
  annotate("segment", x = 214.2751, xend = 400, y = 154.99798, yend = 300) +
  annotate("segment", x = 314.8900, xend = 400, y = 149.53622, yend = 300) +
  annotate("segment", x = 570.4983, xend = 400, y = 6.81035, yend = 300) +
  annotate("label", x = 400, y = 300, label = "Outliers") +
  theme(legend.position = "top")
saveRDS(appendix_figure_3, file = "plots/appendix_figure_3.RDS")

# Remove outliers from data frames
gene_counts <- gene_counts[, !colnames(gene_counts) %in% outliers]
sample_info <- sample_info[!sample_info$id %in% outliers, ]

# Filter lowly expressed genes -------------------------------------------------
# Filter genes that do not pass a threshold of a median expression of 1 RPKM in
# either diet group.

# Generate RPKM matrix
gene_rpkm <- rpkm(gene_counts,
                  gene.length = gene_lengths,
                  normalized.lib.sizes = TRUE)

# Remove genes with median normalized expression < 1
gene_rpkm <- gene_rpkm[apply(gene_rpkm, 1, function(x){
  median(x[sample_info$diet_west == 0]) >= 1 |
    median(x[sample_info$diet_west == 1]) >= 1}), ]

# Remove the same genes from the gene count matrix
gene_counts <- gene_counts[row.names(gene_rpkm), ]

# Regress batch effects out of model -------------------------------------------
# Model gene expression as a function of the CD14 and CD3 expression, RNA
# concentration, and RNA integrity (RIN). Then use the residuals from this model
# in transcriptomic analyses going forward (r_matrix).

# Create design matrix for the base hypothesis
d_diet <- model.matrix(~ sample_info$diet)

# Voom normalize gene counts
v_matrix <- voom(counts = gene_counts, d_diet)

# Model gene expression ~ batch effects to get a residual matrix
d_batch <- model.matrix(~ sample_info$cd14r +
                          sample_info$cd3r +
                          sample_info$rin +
                          sample_info$bioan_conc_ngul)
lm_batch <- eBayes(lmFit(v_matrix, d_batch))  # lin. mod. w/ batch effects
i_batch <- lm_batch$coefficients[, 1]  # intercepts of linear model
r_matrix <- apply(residuals.MArrayLM(object = lm_batch, v_matrix), 2,
                  function(x){x + i_batch})  # matrix of residuals from batch fx

# Kinship data -----------------------------------------------------------------
# Relatedness is estimated using the ngsRelate program (Hanghoj et al., 2019)
# with SNP genotypes inferred from the RNA-seq reads using bcftools mpileup (Li
# et al., 2009). A relatedness matrix is stored in the raw data file.

# Remove subjects not included in the study
kinship <- data.frame(cyno_kin[row.names(cyno_kin) %in% sample_info$id,
                               colnames(cyno_kin) %in% sample_info$id])

# Save output data -------------------------------------------------------------
# Save data for downstream analyses
save(gene_counts, sample_info, r_matrix, kinship,
     file = "data_files/processed_data_1.RData")
