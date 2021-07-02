################################################################################
#
# Behavioral analysis
#
################################################################################

# Format R session -------------------------------------------------------------
# Load packages
library(tidyverse)

# Set theme and colors for graphing
theme_set(theme_classic(base_size = 20, base_family = ""))  # graphing defaults
wes_col <- "#D55E00"  # specify color for Western diet
med_col <- "#56B4E9"  # specify color for Mediterranean diet

# Import raw data
load("data_files/processed_data_1.RData")
load("data_files/raw_behavior.RData")

# Behavior vs. diet-------------------------------------------------------------
# Test for differences in diets in the baseline and treatment phases of the
# experiment.

# Test for differences in behavior b/w diets in during baseline phase
BL_pvals <- apply(BL_behav, 2, function(x){
  wilcox.test(x[sample_info$diet == "MD"],
              x[sample_info$diet == "WD"])$p.value})
sum(p.adjust(BL_pvals, method = "holm") < 0.05)

# Test for differences in behavior b/w diets in during treatment phase
TX_pvals <- apply(TX_behav, 2, function(x){
  wilcox.test(x[sample_info$diet == "MD"],
              x[sample_info$diet == "WD"])$p.value})
sum(p.adjust(TX_pvals, method = "holm") < 0.05)
TX_pvals[p.adjust(TX_pvals, method = "holm") < 0.05]

# Diet Altered Behavior (DAB) --------------------------------------------------
# Conduct a principal component analysis (PCA) of the behavioral data.

# PCA on weighted TX phase raw behaviors using FactoMineR's PCA function
behav_PCA <- PCA(TX_behav, ncp = 35, graph = FALSE,
                 scale.unit = TRUE)

# Add diet coding to the individual PC loadings to create correlation matrix
pca_df <- cbind(behav_PCA$ind$coord,
                "diet" = sample_info$diet,
                "rank" = sample_info$rel_rank) %>%
  as_tibble()

# Correlation matrix between PC loadings and diet
pca_cor <- cor(pca_df)

# Create a data frame to graph variance epxlained
pca_pov <- data.frame(cbind(behav_PCA$eig,
                            pca_cor[1:(ncol(pca_cor) - 2), "diet"],
                            c(1:nrow(behav_PCA$eig)))) %>%
  as_tibble()
colnames(pca_pov) <-
  c("eigenvalue", "pct_var", "cum_pct_var", "cor_diet", "PC")

# T tests of projections between diet groups
pc_diet_test <- apply(pca_df[, 1:20], 2, function(x){
  t.test(x[sample_info$diet == "WD"],
         x[sample_info$diet == "MD"])$p.value})
pca_pov$diet_pval <- pc_diet_test  # add p-values to data frame

# Add DAB score to sample_info matrix
sample_info$DAB <- -behav_PCA$ind$coord[, 2]

# Figures 4A and 4--figure supplement 1 (plots of behavior observed)------------
# The following funciton will generate boxplots with the underlying data
# represented as points for a given data frame (df) and set of behaviors
# (behav_set) with the ability to specifiy axis labels and limits. It is used
# to plot the behavioral observations during the baseline (BL) and treatment
# (TX) phases of the experiment.

# Define function for plotting
behavior.plots <- function(df, behav_set, axis_label, xmin, xmax){
  # Alter df input to prepare for plotting
  plotter <- df %>%
    as.data.frame() %>%
    add_column(diet = sample_info$diet) %>%
    gather(key = "code", value = "value", PCT_TM_B:RT_XSUBrec) %>%
    left_join(behav_codes %>%
                select(code, description),
              by = "code") %>%
    mutate(description = str_to_title(description))
  
  # Summarize the median values of each behavior in each group
  diet_meds <- plotter %>%
    group_by(diet, code) %>%
    summarize(median = median(value)) %>%
    spread(key = diet, value = median) %>%
    mutate(dif_meds = MD - WD, mag_dif_meds = abs(dif_meds))
  
  behavior_medians <- data.frame(
    "median" = apply(df, 2, median),
    "code" = colnames(BL_behav))
  
  behavior_pca <- data.frame(
    "loading" = behav_PCA$var$cor[, 2],
    "code" = names(behav_PCA$var$cor[, 2]))
  
  behav_plot <- plotter %>%
    left_join(behavior_medians, by = "code") %>%
    left_join(behavior_pca, by = "code") %>%
    left_join(diet_meds, by = "code") %>%
    mutate(description = str_remove(description, "Percent Of ") %>%
             str_remove(" Eyes Closed") %>%
             str_replace("Alone Nonstereotypic", "In") %>%
             str_remove("Time ") %>%
             str_remove("Spent ") %>%
             str_remove("In ") %>%
             str_remove("Rate Of ")) %>%
    filter(code %in% behav_set) %>%
    mutate(description = fct_reorder(description, loading, .desc = TRUE)) %>%
    ggplot(aes(x = description, y = value, color = diet)) +
    geom_boxplot(aes(fill = diet), color = "black", alpha = 0.3,
                 outlier.size = 0, coef = 0) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                               jitter.height = 0,
                                               dodge.width = 0.75)) +
    scale_color_manual(values = c(wes_col, med_col),
                       labels = c("WD", "MD"),
                       name = "Diet") +
    scale_fill_manual(values = c(wes_col, med_col),
                      labels = c("WD", "MD"),
                      name = "Diet") +
    scale_y_continuous(limits = c(xmin, xmax)) +
    labs(x = "",
         y = axis_label) +
    coord_flip() +
    theme(legend.position = "top")
  return(behav_plot)
}

# Generate plot for Figure 4A
figure_4a <- behavior.plots(TX_behav, c("PCT_TM_B", "PCT_TM_H", "PCT_TM_R"),
                            "Duration of Behavior (Percent of Time)", NA, NA)
saveRDS(figure_4a, file = "plots/figure_4a.RDS")

# Generate plot for Figure 4--figure supplement 1A
figure_4_figure_supplement_1a <-
  behavior.plots(BL_behav,
                 colnames(BL_behav)[!str_detect(colnames(BL_behav), "PCT_TM")],
                 "Frequency of Behavior (Events per Hour)", 0, 60)
saveRDS(figure_4_figure_supplement_1a,
        file = "plots/figure_4_figure_supplement_1a.RDS")

# Generate plot for Figure 4--figure supplement 1B
figure_4_figure_supplement_1b <-
  behavior.plots(TX_behav,
                 colnames(TX_behav)[!str_detect(colnames(TX_behav), "PCT_TM")],
                 "Frequency of Behavior (Events per Hour)", 0, 60)
saveRDS(figure_4_figure_supplement_1b,
        file = "plots/figure_4_figure_supplement_1b.RDS")

# Generate plot for Figure 4--figure supplement 1C
figure_4_figure_supplement_1c <-
  behavior.plots(BL_behav,
                 colnames(BL_behav)[!str_detect(colnames(BL_behav), "PCT_TM")],
                 "Duration of Behavior (Percent of Time)", 0, 100)
saveRDS(figure_4_figure_supplement_1c,
        file = "plots/figure_4_figure_supplement_1c.RDS")

# Generate plot for Figure 4--figure supplement 1D
figure_4_figure_supplement_1d <-
  behavior.plots(TX_behav,
                 colnames(TX_behav)[!str_detect(colnames(TX_behav), "PCT_TM")],
                 "Duration of Behavior (Percent of Time)", 0, 100)
saveRDS(figure_4_figure_supplement_1d,
        file = "plots/figure_4_figure_supplement_1d.RDS")

# Figure 4B (scree plot of behavioral PCA)--------------------------------------
# Plot the variance explained by each PC of the behavioral PCA and the p-value
# of the difference in projections between diet groups.

# Generate plot 4B
figure_4b <- pca_pov %>%
  ggplot(aes(x = PC, y = pct_var, fill = -log(pc_diet_test))) +
  geom_bar(stat = "identity") +
  labs(x = "Principal Component (PC)",
       y = "Variance Explained (%)") +
  scale_fill_gradient(low = "gray",
                      high = "DarkGreen",
                      name = "Diet Correlation\n-Log10(p)")
saveRDS(figure_4b, file = "plots/figure_4b.RDS")


# Figure 4C (plot of PC2 vs diet)-----------------------------------------------
# Plot individual PC2 projection (DAB score) by diet.

# Generate plot 4C
figure_4c <- sample_info %>%
  ggplot(aes(x = diet, y = DAB, fill = diet, color = diet)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.3, color = "black", coef = 0) +
  geom_jitter(size = 3, width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
  scale_fill_manual(name = "Diet",
                    values = c(wes_col, med_col),
                    labels = c("WD", "MD")) +
  scale_color_manual(name = "Diet",
                     values = c(wes_col, med_col),
                     labels = c("WD", "MD")) +
  scale_x_discrete(labels = c("WD" = "Western", "MD" = "Mediterranean")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  labs(y = "DAB Score")
saveRDS(figure_4c, file = "plots/figure_4c.RDS")

# Figure 4D (correlation between DAB score and behaviors)-----------------------
# Plot the correlation between individual PC2 projection (sign reversed; DAB
# score) and behaviors during the treatment phase. Only show those that are
# significantly correlated (Holm-Bonferroni-adjusted p-value < 0.05).

# Create data frame of correlations
pc2_behav_cor <- data.frame(
  "cor" = apply(TX_behav, 2, function(x){
    cor(sample_info$DAB, x)}),
  "cor_pval" = p.adjust(apply(TX_behav, 2, function(x){
    cor.test(sample_info$DAB, x)$p.value}), method = "holm"),
  "code" = colnames(TX_behav),
  "diet" = apply(TX_behav, 2, function(x){
    cor(sample_info$diet_med, x)})) %>%
  left_join(behav_codes %>%
              select(code, description),
            by = "code") %>%
  mutate(description = str_to_title(description) %>%
           str_remove("Percent Of ") %>%
           str_remove(" Eyes Closed") %>%
           str_replace("Alone Nonstereotypic", "In") %>%
           str_remove("Time ") %>%
           str_remove("Spent ") %>%
           str_remove("In ") %>%
           str_remove("Rate Of "))

# Generate plot 4D
figure_4d <- pc2_behav_cor %>%
  mutate(description = fct_reorder(description, cor)) %>%
  filter(cor_pval < 0.05) %>%
  ggplot(aes(x = description, y = cor)) +
  geom_point(aes(col = factor(ifelse(diet > 0, 1, 0))), size = 5) +
  geom_segment(aes(y = 0,
                   x = description, 
                   yend = cor, 
                   xend = description)) +
  scale_color_manual(breaks = c(0, 1),
                     values = c(wes_col, med_col),
                     guide = FALSE) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(y = "Correlation with Diet-Altered Behavior (DAB)",
       x = "") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "black")
saveRDS(figure_4d, file = "plots/figure_4d.RDS")

# Figure 4--figure supplement 2 (plot of PC1 vs rank)---------------------------
# Plot individual PC1 projections (sign reversed) vs relative rank.

# Test correlation between rank and PC1
cor.test(pca_df$rank, -pca_df$Dim.1)

# Generate figure 4--figure supplement 2
figure_4_figure_supplement_2 <- pca_df %>%
  ggplot(aes(x = rank, y = -Dim.1)) +
  geom_smooth(method = "lm", se = TRUE, col = "black") +
  geom_point(aes(col = as.factor(diet)), size = 4) +
  labs(x = "Relative Rank",
       y = "Projection onto Behavior PC1",
       caption = "r = 0.84") +
  scale_color_manual(name = NULL,
                     breaks = c(1, 2),
                     values = c(wes_col, med_col),
                     labels = c("Western", "Mediterranean")) +
  theme(legend.position = "top",
        plot.caption = element_text(face = "italic"))
saveRDS(figure_4_figure_supplement_2,
        file = "plots/figure_4_figure_supplement_2.RDS")

# Save output data--------------------------------------------------------------
# Save data for downstream analyses
save(gene_counts, sample_info, r_matrix, kinship,
     file = "data_files/processed_data_2.RData")
