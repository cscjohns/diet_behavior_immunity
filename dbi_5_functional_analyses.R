################################################################################
#
# Functional genomic analyses
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


# Create a biomart object with ensembl mfas genes
mart <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                dataset = paste0('mfascicularis', '_gene_ensembl'))

# Create a vector of gene names from expression matrix
ensembl.gene.names <- DEG_info$stable_ID

# Create a biomart object of GO terms associated with genes from study
mfas.GO <- getBM(attributes = c('ensembl_gene_id', 'go_id'),
                 filters = 'ensembl_gene_id',
                 values = ensembl.gene.names,
                 mart = mart)

# Create a list of 
gene.ID2GO <- lapply(unique(mfas.GO$ensembl_gene_id),
                     function(x){sort(mfas.GO[mfas.GO$ensembl_gene_id == x,
                                              'go_id'])})
names(gene.ID2GO) <- unique(mfas.GO$ensembl_gene_id)

# #try opposite list
# go_groups <- lapply(unique(mfas.GO$go_id),
#                      function(x){sort(mfas.GO[mfas.GO$go_id == x,
#                                               'ensembl_gene_id'])})
# names(go_groups) <- unique(mfas.GO$go_id)

## diet.tvals is a vector of diet standardized diet betas (beta/se(beta))
diet_tvals <-
  as.vector(unlist(DEG_info %>% arrange(stable_ID) %>% dplyr::select(std_beta_diet)))
names(diet_tvals) <- as.vector(unlist(DEG_info %>% arrange(stable_ID) %>% dplyr::select(stable_ID)))

# Calculate t-value threshold for FDR cutoff
med_tval_th <- (min(diet_tvals[unlist(DEG_info %>%
                                        filter(diet_FDR < 0.05,
                                               es_diet > 0) %>%
                                        dplyr::select(stable_ID))]) - 0.00001)
west_tval_th <- (max(diet_tvals[unlist(DEG_info %>%
                                         filter(diet_FDR < 0.05,
                                                es_diet < 0) %>%
                                         dplyr::select(stable_ID))]) + 0.00001)

# Create topGOdata objects for mediterranean and western diets
GOdata_med_de <- new('topGOdata',
                     description = 'Simple session',
                     ontology = 'BP',
                     allGenes = diet_tvals,
                     geneSel = function(diet_tvals) diet_tvals > med_tval_th,
                     nodeSize = 10,
                     annot = annFUN.gene2GO,
                     gene2GO = gene.ID2GO)

GOdata_west_de <- new('topGOdata',
                      description = 'Simple session',
                      ontology = 'BP',
                      allGenes = diet_tvals,
                      geneSel = function(diet_tvals) diet_tvals < west_tval_th,
                      nodeSize = 10,
                      annot = annFUN.gene2GO,
                      gene2GO = gene.ID2GO)

# Run Fisher Exact tests to look at enrichment
go_FET_med_de <- runTest(GOdata_med_de,
                         algorithm = 'weight01',
                         statistic = 'fisher')

go_FET_west_de <- runTest(GOdata_west_de,
                          algorithm = 'weight01',
                          statistic = 'fisher')

# Create results tables
go_med_de_results <- GenTable(GOdata_med_de,
                              FET.weight01 = go_FET_med_de,
                              orderBy = 'FET.weight01',
                              ranksOf = 'FET.weight01',
                              topNodes = go_FET_med_de@geneData[4],
                              numChar = 1000)
med_de_go_results <- go_med_de_results %>%
  mutate(pval = as.numeric(FET.weight01),
         qval = qvalue(pval)$qvalues,
         padj = p.adjust(pval, method = "BH"),
         lor = log2(Significant / (Annotated - Significant)) -
           log2((go_FET_med_de@geneData[2] - Significant) /
                  (go_FET_med_de@geneData[1] - Annotated -
                     go_FET_med_de@geneData[2] - Significant))) %>%
  filter(pval < 0.01) %>%
  select(GO.ID, Term, Annotated, Significant, Expected, lor, pval, qval, padj) %>%
  arrange(desc(lor))
write_csv(med_de_go_results, path = "GO_med.csv")

go_west_de_results <- GenTable(GOdata_west_de,
                               FET.weight01 = go_FET_west_de,
                               orderBy = 'FET.weight01',
                               ranksOf = 'FET.weight01',
                               topNodes = go_FET_west_de@geneData[4],
                               numChar = 1000)
west_de_go_results <- go_west_de_results %>%
  mutate(pval = as.numeric(FET.weight01),
         qval = qvalue(as.numeric(pval))$qvalues,
         padj = p.adjust(pval, method = "BH"),
         lor = log2(Significant / (Annotated - Significant)) -
           log2((go_FET_west_de@geneData[2] - Significant) /
                  (go_FET_west_de@geneData[1] - Annotated -
                     go_FET_west_de@geneData[2] - Significant))) %>%
  filter(pval < 0.01) %>%
  select(-c(FET.weight01)) %>%
  select(GO.ID, Term, Annotated, Significant, Expected, lor, pval, qval, padj) %>%
  arrange(desc(lor))
write_csv(west_de_go_results, path = "GO_west.csv")

# Mediation--------------------------------------------------------

# Create gene lists for background sets and mediated genes

west_mediated_genes <- as.vector(unlist(
  mediation_slim %>%
    filter(std_beta_diet < 0) %>%
    arrange(stable_ID) %>%
    select(std_beta_diet)))
names(west_mediated_genes) <- as.vector(unlist(
  mediation_slim %>%
    filter(std_beta_diet < 0) %>%
    arrange(stable_ID) %>%
    select(stable_ID)))

med_mediated_genes <- as.vector(unlist(
  mediation_slim %>%
    filter(std_beta_diet > 0) %>%
    arrange(stable_ID) %>%
    select(std_beta_diet)))

names(med_mediated_genes) <- as.vector(unlist(
  mediation_slim %>%
    filter(std_beta_diet > 0) %>%
    arrange(stable_ID) %>%
    select(stable_ID)))

all_mediated_genes <- c(west_mediated_genes, med_mediated_genes)

west_de_genes <- as.vector(unlist(DEG_info %>%
                                    filter(diet_FDR < 0.05,
                                           es_diet < 0) %>%
                                    arrange(stable_ID) %>%
                                    select(std_beta_diet)))
names(west_de_genes) <- as.vector(unlist(DEG_info %>%
                                           filter(diet_FDR < 0.05,
                                                  es_diet < 0) %>%
                                           arrange(stable_ID) %>%
                                           select(stable_ID)))

med_de_genes <- as.vector(unlist(DEG_info %>%
                                   filter(diet_FDR < 0.05,
                                          es_diet > 0) %>%
                                   arrange(stable_ID) %>%
                                   select(std_beta_diet)))
names(med_de_genes) <- as.vector(unlist(DEG_info %>%
                                          filter(diet_FDR < 0.05,
                                                 es_diet > 0) %>%
                                          arrange(stable_ID) %>%
                                          select(stable_ID)))

all_de_genes <- c(west_de_genes, med_de_genes)

# Create topGOdata objects for mediation
GOdata_med_mediation <- new('topGOdata',
                            description = 'Simple session',
                            ontology = 'BP',
                            allGenes = med_de_genes,
                            geneSel = function(med_de_genes) med_de_genes %in%
                              med_mediated_genes,
                            nodeSize = 5,
                            annot = annFUN.gene2GO,
                            gene2GO = gene.ID2GO)

GOdata_west_mediation <- new('topGOdata',
                             description = 'Simple session',
                             ontology = 'BP',
                             allGenes = west_de_genes,
                             geneSel = function(west_de_genes) west_de_genes %in%
                               west_mediated_genes,
                             nodeSize = 5,
                             annot = annFUN.gene2GO,
                             gene2GO = gene.ID2GO)

GOdata_all_mediation <- new('topGOdata',
                            description = 'Simple session',
                            ontology = 'BP',
                            allGenes = all_de_genes,
                            geneSel = function(all_de_genes) all_de_genes %in%
                              all_mediated_genes,
                            nodeSize = 5,
                            annot = annFUN.gene2GO,
                            gene2GO = gene.ID2GO)


# Run Fisher Exact Tests to look at enrichment
go_FET_med <- runTest(GOdata_med_mediation,
                      algorithm = 'weight01',
                      statistic = 'fisher')

go_FET_west <- runTest(GOdata_west_mediation,
                       algorithm = 'weight01',
                       statistic = 'fisher')

go_FET_all <- runTest(GOdata_all_mediation,
                      algorithm = 'weight01',
                      statistic = 'fisher')

# Create results tables
go_west_results <- GenTable(GOdata_west_mediation,
                            FET.weight01 = go_FET_west,
                            orderBy = 'FET.weight01',
                            ranksOf = 'FET.weight01',
                            topNodes = go_FET_west@geneData[4],
                            numChar = 1000)
west_go_results <- go_west_results %>%
  mutate(pval = as.numeric(FET.weight01),
         qval = qvalue(as.numeric(pval))$qvalues,
         padj = p.adjust(pval, method = "BH"),
         lor = log2(Significant / (Annotated - Significant)) -
           log2((go_FET_west@geneData[2] - Significant) /
                  (go_FET_west@geneData[1] - Annotated -
                     go_FET_west@geneData[2] - Significant))) %>%
  filter(pval < 0.01) %>%
  select(GO.ID, Term, Annotated, Significant, Expected, lor, pval, qval, padj) %>%
  arrange(desc(lor))
write_csv(west_go_results, path = "GO_west_mediated.csv")

go_med_results <- GenTable(GOdata_med_mediation,
                           FET.weight01 = go_FET_med,
                           orderBy = 'FET.weight01',
                           ranksOf = 'FET.weight01',
                           topNodes = go_FET_med@geneData[4],
                           numChar = 1000)
med_go_results <- go_med_results %>%
  mutate(pval = as.numeric(FET.weight01),
         qval = qvalue(as.numeric(pval))$qvalues,
         padj = p.adjust(pval, method = "BH"),
         lor = log2(Significant / (Annotated - Significant)) -
           log2((go_FET_med@geneData[2] - Significant) /
                  (go_FET_med@geneData[1] - Annotated -
                     go_FET_med@geneData[2] - Significant))) %>%
  filter(pval < 0.01) %>%
  select(GO.ID, Term, Annotated, Significant, Expected, lor, pval, qval, padj) %>%
  arrange(desc(lor))
write_csv(med_go_results, path = "GO_med_mediated.csv")

go_all_results <- GenTable(GOdata_all_mediation,
                           FET.weight01 = go_FET_all,
                           orderBy = 'FET.weight01',
                           ranksOf = 'FET.weight01',
                           topNodes = go_FET_all@geneData[4],
                           numChar = 1000)
all_go_results <- go_all_results %>%
  mutate(pval = as.numeric(FET.weight01),
         qval = qvalue(as.numeric(pval))$qvalues,
         padj = p.adjust(pval, method = "BH"),
         lor = log2(Significant / (Annotated - Significant)) -
           log2((go_FET_all@geneData[2] - Significant) /
                  (go_FET_all@geneData[1] - Annotated -
                     go_FET_all@geneData[2] - Significant))) %>%
  filter(pval < 0.01) %>%
  select(GO.ID, Term, Annotated, Significant, Expected, lor, pval, qval, padj) %>%
  arrange(desc(lor))
write_csv(all_go_results, path = "GO_all_mediated.csv")

# Old stuff ----------------------------------------
# vector of significance levels to test
sig.levels <- c(0.10, 0.05, 0.01)
for(i in 1:length(sig.levels)){
  # Set significance level
  sig.level <- sig.levels[i]
  
  # Calculate t-value threshold for FDR cutoff
  med.tval.th <-
    (min(diet.tvals[row.names(DEG_info[DEG_info$diet_FDR < sig.level &
                                         DEG_info$es_diet > 0, ])])
     - 0.00001)
  wes.tval.th <-
    (max(diet.tvals[row.names(DEG_info[DEG_info$diet_FDR < sig.level &
                                         DEG_info$es_diet < 0, ])])
     + 0.00001)
  
  # Create topGOdata objects for mediterranean and western diets
  GOdata.upreg.med <- new('topGOdata',
                          description = 'Simple session',
                          ontology = 'BP',
                          allGenes = diet.tvals,
                          geneSel = function(diet.tvals) diet.tvals > med.tval.th,
                          nodeSize = 10,
                          annot = annFUN.gene2GO,
                          gene2GO = gene.ID2GO)
  
  
  # Run four different tests to look at enrichment for MD
  FET.weight01 <- runTest(GOdata.upreg.med, algorithm = 'weight01',
                          statistic = 'fisher')
  FET.parentchild <- runTest(GOdata.upreg.med, algorithm = 'parentchild',
                             statistic = 'fisher')
  KS.weight01 <- runTest(GOdata.upreg.med, algorithm = 'weight01',
                         statistic = 'ks')
  KS.elim <- runTest(GOdata.upreg.med, algorithm = 'elim', statistic = 'ks')
  
  # Gather the top nodes from tests
  top.nodes <- max(c(FET.weight01@geneData[4], FET.parentchild@geneData[4],
                     KS.weight01@geneData[4], KS.elim@geneData[4]))
  
  # Create results object
  GOdata.upreg.med.results <- GenTable(GOdata.upreg.med,
                                       FET.weight01 = FET.weight01,
                                       FET.parentchild = FET.parentchild,
                                       KS.weight01 = KS.weight01,
                                       KS.elim = KS.elim,
                                       orderBy = 'FET.weight01',
                                       ranksOf = 'FET.weight01',
                                       topNodes = top.nodes)
  
  # rename column headings
  for(j in 1:length(colnames(GOdata.upreg.med.results))){
    colnames(GOdata.upreg.med.results)[j] <-
      paste0(as.character(sig.level), '_', colnames(GOdata.upreg.med.results)[j])
  }
  
  # Store results in data frame
  ifelse(i == 1, all.topGO.med.results <- GOdata.upreg.med.results[order(
    GOdata.upreg.med.results[,1]),],
    all.topGO.med.results <-
      cbind(all.topGO.med.results,
            GOdata.upreg.med.results[order(
              GOdata.upreg.med.results[,1]),]))
  
  # Graph the results
  printGraph(GOdata.upreg.med, KS.weight01, firstSigNodes = 5,
             fn.prefix = paste0(out_dir, 'topGO_med_', sig.level, '_KSweight01'),
             useInfo = 'all', pdfSW = TRUE)
  printGraph(GOdata.upreg.med, FET.weight01, firstSigNodes = 5,
             fn.prefix = paste0(out_dir, 'topGO_med_', sig.level, '_FETweight01'),
             useInfo = 'all', pdfSW = TRUE)
  printGraph(GOdata.upreg.med, KS.elim, firstSigNodes = 5,
             fn.prefix = paste0(out_dir, 'topGO_med_', sig.level, '_KSelim'),
             useInfo = 'all', pdfSW = TRUE)
  printGraph(GOdata.upreg.med, FET.parentchild, firstSigNodes = 5,
             fn.prefix = paste0(out_dir, 'topGO_med_', sig.level, '_FETparentchild'),
             useInfo = 'all', pdfSW = TRUE)
  
  # Repeat for WD
  GOdata.upreg.wes <- new('topGOdata',
                          description = 'Simple session',
                          ontology = 'BP',
                          allGenes = diet.tvals,
                          geneSel = function(diet.tvals) diet.tvals < wes.tval.th,
                          nodeSize = 10,
                          annot = annFUN.gene2GO,
                          gene2GO = gene.ID2GO)
  FET.weight01 <- runTest(GOdata.upreg.wes, algorithm = 'weight01',
                          statistic = 'fisher')
  FET.parentchild <- runTest(GOdata.upreg.wes, algorithm = 'parentchild',
                             statistic = 'fisher')
  KS.weight01 <- runTest(GOdata.upreg.wes, algorithm = 'weight01',
                         statistic = 'ks')
  KS.elim <- runTest(GOdata.upreg.wes, algorithm = 'elim', statistic = 'ks')
  top.nodes <- max(c(FET.weight01@geneData[4], FET.parentchild@geneData[4],
                     KS.weight01@geneData[4], KS.elim@geneData[4]))
  GOdata.upreg.wes.results <- GenTable(GOdata.upreg.wes,
                                       FET.weight01 = FET.weight01,
                                       FET.parentchild = FET.parentchild,
                                       KS.weight01 = KS.weight01,
                                       KS.elim = KS.elim,
                                       orderBy = 'FET.weight01',
                                       ranksOf = 'FET.weight01',
                                       topNodes = top.nodes)
  for(j in 1:length(colnames(GOdata.upreg.wes.results))){
    colnames(GOdata.upreg.wes.results)[j] <-
      paste0(as.character(sig.level), '_', colnames(GOdata.upreg.wes.results)[j])
  }
  
  ifelse(i == 1, all.topGO.wes.results <- GOdata.upreg.wes.results[order(
    GOdata.upreg.wes.results[,1]),],
    all.topGO.wes.results <-
      cbind(all.topGO.wes.results,
            GOdata.upreg.wes.results[order(
              GOdata.upreg.wes.results[,1]),]))
  
  
  # Graph the results
  printGraph(GOdata.upreg.wes, KS.weight01, firstSigNodes = 5,
             fn.prefix = paste0(out_dir, 'topGO_wes_', sig.level, '_KSweight01'),
             useInfo = 'all', pdfSW = TRUE)
  printGraph(GOdata.upreg.wes, FET.weight01, firstSigNodes = 5,
             fn.prefix = paste0(out_dir, 'topGO_wes_', sig.level, '_FETweight01'),
             useInfo = 'all', pdfSW = TRUE)
  printGraph(GOdata.upreg.wes, KS.elim, firstSigNodes = 5,
             fn.prefix = paste0(out_dir, 'topGO_wes_', sig.level, '_KSelim'),
             useInfo = 'all', pdfSW = TRUE)
  printGraph(GOdata.upreg.wes, FET.parentchild, firstSigNodes = 5,
             fn.prefix = paste0(out_dir, 'topGO_wes_', sig.level, '_FETparentchild'),
             useInfo = 'all', pdfSW = TRUE)
  
  # Write data to table
  write.table(GOdata.upreg.wes.results,
              file = paste0(out_dir, 'topGO_wes_', sig.level, '_upregulated.csv'),
              append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE,
              sep = ',')
  write.table(GOdata.upreg.med.results,
              file = paste0(out_dir, 'topGO_med_', sig.level, '_upregulated.csv'),
              append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE,
              sep = ',')
}



old.cols <- c('0.1_FET.weight01', '0.1_FET.parentchild',
              '0.1_KS.weight01', '0.1_KS.elim',
              '0.05_FET.weight01', '0.05_FET.parentchild',
              '0.05_KS.weight01', '0.05_KS.elim',
              '0.01_FET.weight01', '0.01_FET.parentchild',
              '0.01_KS.weight01', '0.01_KS.elim')
new.cols <- c('padj_0.1_FET.weight01', 'padj_0.1_FET.parentchild',
              'padj_0.1_KS.weight01', 'padj_0.1_KS.elim',
              'padj_0.05_FET.weight01', 'padj_0.05_FET.parentchild',
              'padj_0.05_KS.weight01', 'padj_0.05_KS.elim',
              'padj_0.01_FET.weight01', 'padj_0.01_FET.parentchild',
              'padj_0.01_KS.weight01', 'padj_0.01_KS.elim')

for(i in 1:12){
  all.topGO.med.results[, new.cols[i]] <-
    p.adjust(all.topGO.med.results[, old.cols[i]], method = 'BH')
  all.topGO.wes.results[, new.cols[i]] <-
    p.adjust(all.topGO.wes.results[, old.cols[i]], method = 'BH')
}

# Create a list of significant GO terms for each test and threshold
sig.GO.terms <- list()
sig.GO.terms[1:12] <- lapply(old.cols,
                             function(x){all.topGO.med.results[
                               (!is.na(all.topGO.med.results[,x]) &
                                  all.topGO.med.results[, x] < 0.05),
                               '0.1_Term']})
names(sig.GO.terms)[1:12] <-
  lapply(old.cols,function(x){as.character(paste0('med_upregulated_', x))})

sig.GO.terms[13:24] <- lapply(old.cols,
                              function(x){all.topGO.wes.results[
                                (!is.na(all.topGO.wes.results[,x]) &
                                   all.topGO.wes.results[, x] < 0.05),
                                '0.1_Term']})
names(sig.GO.terms)[13:24] <-
  lapply(old.cols,function(x){as.character(paste0('wes_upregulated_', x))})

# Create a list of significnat GO IDs for each test and threshold
sig.GO.IDs <- list()
sig.GO.IDs[1:12] <- lapply(old.cols,
                           function(x){all.topGO.med.results[
                             (!is.na(all.topGO.med.results[,x]) &
                                all.topGO.med.results[, x] < 0.05),
                             '0.1_GO.ID']})
names(sig.GO.IDs)[1:12] <-
  lapply(old.cols,function(x){as.character(paste0('med_upregulated_', x))})

sig.GO.IDs[13:24] <- lapply(old.cols,
                            function(x){all.topGO.wes.results[
                              (!is.na(all.topGO.wes.results[,x]) &
                                 all.topGO.wes.results[, x] < 0.05),
                              '0.1_GO.ID']})
names(sig.GO.IDs)[13:24] <-
  lapply(old.cols,function(x){as.character(paste0('wes_upregulated_', x))})

number.GO.sig <- data.frame(matrix(lapply(sig.GO.terms, length),
                                   byrow = FALSE, nrow = 12))
colnames(number.GO.sig) <- c('med_upregulated', 'wes_upregulated')
row.names(number.GO.sig) <- new.cols

sig.GO <- list(number.GO.sig, sig.GO.IDs, sig.GO.terms)

write.table(all.topGO.med.results[
  all.topGO.med.results$padj_0.05_FET.parentchild < 0.05,
  c('0.1_GO.ID', 'padj_0.05_FET.parentchild')],
  file = paste0(out_dir, 'GO_IDs_med.txt'),
  append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE,
  sep = '\t')

write.table(all.topGO.wes.results[
  all.topGO.wes.results$padj_0.05_FET.parentchild < 0.05,
  c('0.1_GO.ID', 'padj_0.05_FET.parentchild')],
  file = paste0(out_dir, 'GO_IDs_wes.txt'),
  append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE,
  sep = '\t')

write.table(all.topGO.med.results,
            file = paste0(out_dir, 'GO_med_upregulated_final.tsv'),
            append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE,
            sep = '\t')

write.table(all.topGO.wes.results,
            file = paste0(out_dir, 'GO_wes_upregulated_final.tsv'),
            append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE,
            sep = '\t')

# Testing enrichment for disease-associated genes ---------------

# Load in gene sets from ptwas paper
ptwas_in <- read_delim("~/Documents/2_data/ptwas_data.txt",
                       delim = " ") %>%
  clean_names()

# Load mart for converting human ensembl gene ids to M. fascicularis
c2h_table <- read_delim("~/Documents/2_data/cyno_to_human.txt",
                        delim = "\t") %>%
  clean_names()

# Rename the columns
colnames(c2h_table) <- c("human_stable_id", "cyno_stable_id")

# Identify MED/WEST/non_sig genes in the PTWAS gene set
ptwas_data <- ptwas_in %>%
  separate(col = "gene",
           into = c("human_stable_id", NA)) %>%  # remove gene variant id
  left_join(c2h_table,
            by = "human_stable_id") %>%  # use c2h table to ID orthologues
  left_join(DEG_info %>%
              mutate(diet_genelist = ifelse(diet_sig == "TRUE", diet_dir, "not_sig")) %>%
              select(stable_ID, diet_genelist),
            by = c("cyno_stable_id" = "stable_ID")) %>%
  replace_na(list(rep("missing", 6)))

ptwas_data$diet_genelist <- replace_na(ptwas_data$diet_genelist, "missing")

west_enrich <- ptwas_data %>%
  group_by(trait, diet_genelist) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = diet_genelist,
              values_from = n,
              values_fill = 0) %>%
  mutate(lor = log2(fisher.test(matrix(c(WEST,
                                         MED + not_sig,
                                         2236 - WEST,
                                         10004 - MED - not_sig),
                                       nrow = 2),
                                alternative = "g")$estimate),
         loglow_ci = log2(fisher.test(matrix(c(WEST,
                                               MED + not_sig,
                                               2236 - WEST,
                                               10004 - MED - not_sig),
                                             nrow = 2),
                                      alternative = "g")$conf.int[1]),
         loghi_ci = log2(fisher.test(matrix(c(WEST,
                                              MED + not_sig,
                                              2236 - WEST,
                                              10004 - MED - not_sig),
                                            nrow = 2),
                                     alternative = "g")$conf.int[2]),
         pval = fisher.test(matrix(c(WEST,
                                     MED + not_sig,
                                     2236 - WEST,
                                     10004 - MED - not_sig),
                                   nrow = 2),
                            alternative = "g")$p.value) %>%
  arrange(desc(lor))

west_enrich$qval <- qvalue(west_enrich$pval)$qvalues
west_enrich$padj <- p.adjust(west_enrich$pval)

med_enrich <- ptwas_data %>%
  group_by(trait, diet_genelist) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = diet_genelist,
              values_from = n,
              values_fill = 0) %>%
  mutate(lor = log2(fisher.test(matrix(c(MED,
                                         WEST + not_sig,
                                         2664 - MED,
                                         9576 - WEST - not_sig),
                                       nrow = 2),
                                alternative = "g")$estimate),
         pval = fisher.test(matrix(c(MED,
                                     WEST + not_sig,
                                     2664 - MED,
                                     9576 - WEST - not_sig),
                                   nrow = 2),
                            alternative = "g")$p.value) %>%
  arrange(desc(lor))

med_enrich$qval <- qvalue(med_enrich$pval)$qvalues

deg_enrich <- ptwas_data %>%
  group_by(trait, diet_genelist) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = diet_genelist,
              values_from = n,
              values_fill = 0) %>%
  mutate(or = fisher.test(matrix(c(MED + WEST,
                                   not_sig,
                                   4900 - MED - WEST,
                                   7340 - not_sig),
                                 nrow = 2),
                          alternative = "g")$estimate,
         pval = fisher.test(matrix(c(MED + WEST,
                                     not_sig,
                                     4900 - MED,
                                     7340 - not_sig),
                                   nrow = 2),
                            alternative = "g")$p.value) %>%
  arrange(desc(or))

deg_enrich$qval <- qvalue(deg_enrich$pval)$qvalues

west_report <- west_enrich %>%
  filter(qval < 0.2) %>%
  ungroup() %>%
  select(-trait) %>%
  add_column(
    source = c(
      "IMMUNOBASE, hg19",
      "ADIPOGen",
      "GLGC, MC",
      "GLGC, MC",
      "SSGAC",
      "UKB_6152_9",
      "Astle* et al., *2016",
      "UKB_23099",
      "UKB_21001"),
    trait = c(
      "Celiac Disease",
      "Adiponectin",
      "LDL-C",
      "HDL-C",
      "Education Years Pooled",
      "Hayfever, Allergic Rhinitis, or Eczema",
      "Eosinophil Counts",
      "Body Fat Percentage",
      "Body Mass Index (BMI)")) %>%
  select(trait, WEST, MED, not_sig, missing, lor, pval, qval, everything()) %>%
  mutate(trait = fct_reorder(trait, lor, .desc = FALSE)) %>%
  filter(trait != "Education Years Pooled")

west_report %>%
  mutate(sig = ifelse(qval < 0.2, "sig", "not_sig"),
         xtext = str_c("**", trait, "**,<br>*(", source, ")*"),
         xtext = fct_reorder(xtext, lor, .desc = FALSE)) %>%
  filter(qval < 0.2) %>%
  ggplot(aes(x = xtext, y = lor)) +
  geom_errorbar(aes(x = xtext,
                    ymin = loglow_ci,
                    ymax = loghi_ci),
                col = "black", size = 1, width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 1) +
  geom_point(aes(col = sig), size = 4) +
  scale_color_brewer(palette = "Dark2",
                     breaks = c("sig", "not_sig"),
                     labels = c("< 0.2", "> 0.2"),
                     name = "q Value",
                     guide = NULL) +
  labs(x = NULL,
       y = expression(Log[2]*"(Odds Ratio)"),
       title = "Traits Enriched in Western Genes") +
  coord_flip() +
  theme(legend.position = "top",
        axis.text.y = ggtext::element_markdown())

ggsave("~/Desktop/ptwas_plot.pdf",
       device = "pdf",
       width = 10,
       height = 5,
       units = "in")
