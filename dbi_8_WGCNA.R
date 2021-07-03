################################################################################
#
# Weighted Gene Co-Expression Network Analysis (WGCNA)
#
################################################################################

# Display the current working directory
getwd()

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#Read in the WGCNA data
load("diet_pbmc_data/clean_wgcna_data.RData")

#===============================================================================
#
#  Code chunk 2
#
#===============================================================================

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(t_matrix, powerVector = powers, verbose = 5,
                         blockSize = 12240)
# Plot the results:
pdf("soft_thresholding.pdf")
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,cex = cex1, col = "red")
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.75,col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers,
     cex = cex1, col = "red")

dev.off()


#===============================================================================
#
#  Code chunk 3
#
#===============================================================================


softPower <- 6
adjacency <- adjacency(t_matrix, power = softPower)


#===============================================================================
#
#  Code chunk 4
#
#===============================================================================


# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM


#===============================================================================
#
#  Code chunk 5
#
#===============================================================================


# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
pdf("dendrogram_1.pdf")
sizeGrWindow(12, 9)
plot(geneTree, xlab = "", sub = "",
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


#===============================================================================
#
#  Code chunk 6
#
#===============================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize <- 30
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)


#===============================================================================
#
#  Code chunk 7
#
#===============================================================================


# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf("dendrogram_2.pdf")
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#  Code chunk 8
#
#===============================================================================


# Calculate eigengenes
MEList <- moduleEigengenes(t_matrix, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# Plot the result
pdf("ME_cluster.pdf")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#===============================================================================
#
#  Code chunk 9
#
#===============================================================================


MEDissThres <- 0.25
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
dev.off()
# Call an automatic merging function
merged <- mergeCloseModules(t_matrix, dynamicColors, cutHeight = MEDissThres,
                            verbose = 3)
# The merged module colors
mergedColors <- merged$colors
# Eigengenes of the new merged modules:
mergedMEs <- merged$newMEs


#===============================================================================
#
#  Code chunk 10
#
#===============================================================================

pdf(file = "dendrogram_3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#===============================================================================
#
#  Code chunk 11
#
#===============================================================================


# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, mergedMEs, dynamicColors, mergedColors, geneTree,
     file = "cyno_wgcna_output.RData")


#===============================================================================
#
#  Correlation between diet and module eigengenes
#
#===============================================================================

eigengene_diet_cor <- tibble(
  "module" = str_c("Module ", 1:16),
  "module_color" = str_remove(names(mergedMEs), "^ME"),
  "cor_diet" = apply(mergedMEs, 2, function(x){
    cor.test(sample_info$diet_med, x)$estimate}),
  "cor_pval" = apply(mergedMEs, 2, function(x){
    cor.test(sample_info$diet_med, x)$p.value}),
  "lowci" = apply(mergedMEs, 2, function(x){
    cor.test(sample_info$diet_med, x)$conf.int[1]}),
  "hici" = apply(mergedMEs, 2, function(x){
    cor.test(sample_info$diet_med, x)$conf.int[2]}),
  "ttest_t" = apply(mergedMEs, 2, function(x){
    t.test(x[sample_info$diet_med == 0],
           x[sample_info$diet_med == 1])$statistic}),
  "ttest_df" = apply(mergedMEs, 2, function(x){
    t.test(x[sample_info$diet_med == 0],
           x[sample_info$diet_med == 1])$parameter}),
  "ttest_pval" = apply(mergedMEs, 2, function(x){
    t.test(x[sample_info$diet_med == 0],
           x[sample_info$diet_med == 1])$p.value})
)
eigengene_diet_cor$padj <- p.adjust(eigengene_diet_cor$ttest_pval)

#===============================================================================
#
#  Enrichment for diet genes in modules
#
#===============================================================================

module_deg_info <- DEG_info %>%
  add_column("module_color" = mergedColors) %>%
  left_join(eigengene_diet_cor %>%
              select(module, module_color),
            by = "module_color")

module_fet <- module_deg_info %>%
  mutate(diet_dir = ifelse(diet_sig == "TRUE", diet_dir, "ns")) %>%
  group_by(module, diet_dir) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = "diet_dir",
              values_from = "n") %>%
  replace_na(list(MED = 0, WEST = 0, ns = 0)) %>%
  mutate(n_genes = MED + ns + WEST,
         WD_lor = log2(fisher.test(matrix(c(WEST,
                                            MED + ns,
                                            2236 - WEST,
                                            10004 - MED - ns),
                                          nrow = 2))$estimate),
         WD_lowci = log2(fisher.test(matrix(c(WEST,
                                              MED + ns,
                                              2236 - WEST,
                                              10004 - MED - ns),
                                            nrow = 2))$conf.int[1]),
         WD_hici = log2(fisher.test(matrix(c(WEST,
                                             MED + ns,
                                             2236 - WEST,
                                             10004 - MED - ns),
                                           nrow = 2))$conf.int[2]),
         WD_pval = fisher.test(matrix(c(WEST,
                                        MED + ns,
                                        2236 - WEST,
                                        10004 - MED - ns),
                                      nrow = 2))$p.value,
         MD_lor = log2(fisher.test(matrix(c(MED,
                                            WEST + ns,
                                            2664 - MED,
                                            9576 - WEST - ns),
                                          nrow = 2))$estimate),
         MD_lowci = log2(fisher.test(matrix(c(MED,
                                              WEST + ns,
                                              2664 - MED,
                                              9576 - WEST - ns),
                                            nrow = 2))$conf.int[1]),
         MD_hici = log2(fisher.test(matrix(c(MED,
                                             WEST + ns,
                                             2664 - MED,
                                             9576 - WEST - ns),
                                           nrow = 2))$conf.int[2]),
         MD_pval = fisher.test(matrix(c(MED,
                                        WEST + ns,
                                        2664 - MED,
                                        9576 - WEST - ns),
                                      nrow = 2))$p.value) %>%
  pivot_longer(WD_lor:MD_pval,
               names_to = "measure",
               values_to = "value") %>%
  mutate(value = ifelse(value == -Inf, -12, value)) %>%
  separate(measure,
           into = c("diet", "measure"),
           sep = "_") %>%
  pivot_wider(names_from = "measure",
              values_from = "value") %>%
  left_join(eigengene_diet_cor %>%
              select(module, cor_diet),
            by = "module") %>%
  mutate(module_label = as_factor(str_c(module, " (", n_genes, ")")),
         diet_genes = ifelse(diet == "MD", MED, WEST),
         diet = ifelse(diet == "MD", "Mediterranean", "Western"),
         order = as.numeric(str_remove(module, "Module ")))

module_fet$padj <- p.adjust(module_fet$pval)

module_fet %>%
  filter(n_genes > 30) %>%
  ggplot(aes(x = fct_reorder(module_label, desc(order)), y = lor)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
  geom_errorbar(aes(ymin = lowci, ymax = hici, col = diet),
                size = 3, alpha = 0.6, width = 0) +
  geom_label(aes(label = diet_genes, fill = diet), col = "white") +
  scale_color_manual(breaks = c("Mediterranean", "Western"),
                     values = c(med_col, wes_col),
                     guide = NULL) +
  scale_fill_manual(breaks = c("Mediterranean", "Western"),
                    values = c(med_col, wes_col),
                    guide = NULL) +
  scale_y_continuous(limits = c(-12, 4),
                     breaks = seq(-12, 4, 2),
                     labels = c("-Inf", seq(-10, 4, 2))) +
  labs(x = NULL,
       y = expression(Log[2]*"(Odds Ratio)"),
       title = "Enrichment of DE Genes in WGCNA Modules") +
  coord_flip()

module_fet %>%
  filter(n_genes > 30) %>%
  ggplot(aes(x = fct_reorder(module_label, desc(order)), y = lor)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
  geom_errorbar(aes(ymin = lowci, ymax = hici, col = ifelse(lowci * hici > 0, diet, "gray")),
                size = 3, alpha = 0.6, width = 0) +
  geom_label(aes(label = diet_genes, fill = ifelse(lowci * hici > 0, diet, "gray")), col = "white") +
  scale_color_manual(breaks = c("Mediterranean", "Western", "gray"),
                     values = c(med_col, wes_col, "gray"),
                     guide = NULL) +
  scale_fill_manual(breaks = c("Mediterranean", "Western", "gray"),
                    values = c(med_col, wes_col, "gray"),
                    guide = NULL) +
  scale_y_continuous(limits = c(-12, 4),
                     breaks = seq(-12, 4, 2),
                     labels = c("-Inf", seq(-10, 4, 2))) +
  labs(x = NULL,
       y = expression(Log[2]*"(Odds Ratio)"),
       title = "Enrichment of DE Genes in WGCNA Modules") +
  coord_flip() +
  facet_wrap(~diet)

#===============================================================================
#
#  GO enrichment
#
#===============================================================================

library(topGO)
library(biomaRt)

# Create a biomart object with ensembl mfas genes
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mfascicularis_gene_ensembl")

# Create a vector of gene names from expression matrix
ensembl_gene_names <- DEG_info$stable_ID

# Create a biomart object of GO terms associated with genes from study
mfas_GO <- getBM(attributes = c("ensembl_gene_id", "go_id"),
                 filters = "ensembl_gene_id",
                 values = ensembl_gene_names,
                 mart = mart)

# Create a list of 
gene_ID2GO <- lapply(unique(mfas_GO$ensembl_gene_id),
                     function(x){sort(mfas_GO[mfas_GO$ensembl_gene_id == x,
                                              "go_id"])})
names(gene_ID2GO) <- unique(mfas_GO$ensembl_gene_id)

# Create a list of module gene sets
module_genesets <- lapply(unique(mdi$module),
                          function(x){sort(mdi[mdi$module == x,
                                               "stable_ID"])})
names(module_genesets) <- unique(mdi$module)

all_genes <- DEG_info$stable_ID
names(all_genes) <- DEG_info$stable_ID

# Create topGOdata objects for each module geneset
module_topGOdata <- list()
module_go_FET <- list()
module_go_results <- list()
module_res2 <- list()

for(i in 1:length(module_genesets)){
  geneList <- factor(as.integer(all_genes %in% module_genesets[[i]]))
  names(geneList) <- names(all_genes)
  module_topGOdata[i] <- new("topGOdata",
                             description = "Simple session",
                             ontology = "BP",
                             allGenes = geneList,
                             nodeSize = 10,
                             annot = annFUN.gene2GO,
                             gene2GO = gene_ID2GO)
  this_topGOdata <- module_topGOdata[[i]]
  module_go_FET[i] <- runTest(this_topGOdata,
                              algorithm = "weight01",
                              statistic = "fisher")
  this_go_FET <- module_go_FET[[i]]
  module_go_results[[i]] <- GenTable(this_topGOdata,
                                     FET.weight01 = this_go_FET,
                                     orderBy = "FET.weight01",
                                     ranksOf = "FET.weight01",
                                     topNodes = this_go_FET@geneData[4],
                                     numChar = 1000)
  this_go_results <- module_go_results[[i]]
  module_res2[[i]] <- this_go_results %>%
    mutate(pval = as.numeric(FET.weight01),
           qval = qvalue(pval)$qvalues,
           padj = p.adjust(pval, method = "BH"),
           lor = log2(Significant / (Annotated - Significant)) -
             log2((this_go_FET@geneData[2] - Significant) /
                    (this_go_FET@geneData[1] - Annotated -
                       this_go_FET@geneData[2] - Significant))) %>%
    filter(pval < 0.01) %>%
    select(GO.ID, Term, Annotated, Significant, Expected, lor, pval, qval, padj) %>%
    arrange(desc(lor))
}

names(module_topGOdata) <- names(module_genesets)
names(module_go_FET) <- names(module_genesets)
names(module_go_results) <- names(module_genesets)
names(module_res2) <- names(module_genesets)


wgcna_topGO_results <- module_res2[[1]]
wgcna_topGO_results$module <- rep(names(module_res2)[1],
                                  times = nrow(module_res2[[1]]))

for(i in 2:length(module_genesets)){
  tempdf <- module_res2[[i]]
  tempdf$module <- rep(names(module_res2)[i],
                       times = nrow(module_res2[[i]]))
  wgcna_topGO_results <- rbind(wgcna_topGO_results,
                               tempdf)
}

wgcna_topGO_results %>%
  left_join(eigengene_diet_cor %>%
              select(module, cor_diet),
            by = "module") %>%
  arrange(desc(cor_diet), padj) %>%
  select(module, cor_diet, GO.ID, Term, Annotated, Significant, Expected, lor,
         pval, padj) %>%
  write_csv("wgcna_topGO_results.csv")

#===============================================================================
#
#  Enrichment for polarization genes in modules
#
#===============================================================================

module_pol <- module_deg_info %>%
  mutate(polarization = ifelse(polarization == "Pro-Inflammatory", "PRO",
                               ifelse(polarization == "Regulatory", "REG", "ns"))) %>%
  group_by(module, polarization) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = "polarization",
              values_from = "n") %>%
  replace_na(list(REG = 0, PRO = 0, ns = 0)) %>%
  mutate(n_genes = REG + ns + PRO,
         PR_lor = log2(fisher.test(matrix(c(PRO,
                                            REG + ns,
                                            698 - PRO,
                                            11542 - REG - ns),
                                          nrow = 2))$estimate),
         PR_lowci = log2(fisher.test(matrix(c(PRO,
                                              REG + ns,
                                              698 - PRO,
                                              11542 - REG - ns),
                                            nrow = 2))$conf.int[1]),
         PR_hici = log2(fisher.test(matrix(c(PRO,
                                             REG + ns,
                                             698 - PRO,
                                             11542 - REG - ns),
                                           nrow = 2))$conf.int[2]),
         PR_pval = fisher.test(matrix(c(PRO,
                                        REG + ns,
                                        698 - PRO,
                                        11542 - REG - ns),
                                      nrow = 2))$p.value,
         RG_lor = log2(fisher.test(matrix(c(REG,
                                            PRO + ns,
                                            138 - REG,
                                            12102 - PRO - ns),
                                          nrow = 2))$estimate),
         RG_lowci = log2(fisher.test(matrix(c(REG,
                                              PRO + ns,
                                              138 - REG,
                                              12102 - PRO - ns),
                                            nrow = 2))$conf.int[1]),
         RG_hici = log2(fisher.test(matrix(c(REG,
                                             PRO + ns,
                                             138 - REG,
                                             12102 - PRO - ns),
                                           nrow = 2))$conf.int[2]),
         RG_pval = fisher.test(matrix(c(REG,
                                        PRO + ns,
                                        138 - REG,
                                        12102 - PRO - ns),
                                      nrow = 2))$p.value) %>%
  pivot_longer(PR_lor:RG_pval,
               names_to = "measure",
               values_to = "value") %>%
  mutate(value = ifelse(value == -Inf, -10, value)) %>%
  separate(measure,
           into = c("polarization", "measure"),
           sep = "_") %>%
  pivot_wider(names_from = "measure",
              values_from = "value") %>%
  left_join(eigengene_diet_cor %>%
              select(module, cor_diet),
            by = "module") %>%
  mutate(module_label = as_factor(str_c(module, " (", n_genes, ")")),
         polarization_genes = ifelse(polarization == "RG", REG, PRO),
         polarization = ifelse(polarization == "RG", "Regulatory", "Pro-Inflammatory"),
         order = as.numeric(str_remove(module, "Module ")))

module_pol$padj <- p.adjust(module_pol$pval)

module_pol %>%
  filter(n_genes > 30) %>%
  ggplot(aes(x = fct_reorder(module_label, desc(order)), y = lor)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
  geom_errorbar(aes(ymin = lowci, ymax = hici, col = ifelse(lowci * hici > 0, polarization, "gray")),
                size = 3, alpha = 0.6, width = 0) +
  geom_label(aes(label = polarization_genes, fill = ifelse(lowci * hici > 0, polarization, "gray")), col = "white") +
  scale_color_manual(breaks = c("Regulatory", "Pro-Inflammatory", "gray"),
                     values = c(m2_col, m1_col, "gray"),
                     guide = NULL) +
  scale_fill_manual(breaks = c("Regulatory", "Pro-Inflammatory", "gray"),
                    values = c(m2_col, m1_col, "gray"),
                    guide = NULL) +
  scale_y_continuous(limits = c(-10, 8),
                     breaks = seq(-10, 8, 2),
                     labels = c("-Inf", seq(-8, 8, 2))) +
  labs(x = NULL,
       y = expression(Log[2]*"(Odds Ratio)"),
       title = "Enrichment of Polarization Genes in WGCNA Modules") +
  coord_flip() +
  facet_wrap(~polarization)
