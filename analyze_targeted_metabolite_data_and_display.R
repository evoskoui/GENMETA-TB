################################################################################
### Step 0: Import not log2 norm and median intensity normalized data
################################################################################
# rows are samples (non avg) and metabolites are columns
# Load libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggfortify)
library(plotly)
library(forcats)

# Set font
library(showtext)
font_add(family = "Times New Roman", 
         regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf")
showtext_auto()
showtext_opts(dpi = 300)

# Load pareto log 10 norm input data
lp_df_file <- "./Scaling/Targeted/log10_pareto_data_normalized.csv"
metab_df <- read_csv(lp_df_file)
names(metab_df)[1] <- "sample"

metadata_file <- "./Scaling/Targeted/norm_metadata.csv"
group_df <- read_csv(metadata_file) %>% 
  mutate(lineage = if_else(strain == "H37Rv", "Ref", lineage))


################################################################################
### Step 1: Get LFC values compared to reference H37Rv (L4) for KEGG ID met
################################################################################
# Calculate LFCs (optional, for later filtering)
# Reshape metabolite data
metab_long <- metab_df %>%
  pivot_longer(cols = -sample, names_to = "metabolite", values_to = "value")

# Join with group info to add lineage
annotated_df <- metab_long %>%
  left_join(group_df, by = "sample")

# Average over strain, reduce noise from technical repeats, include lineage
strain_avg_df <- annotated_df %>%
  group_by(metabolite, strain, lineage) %>%
  summarise(smean_value = mean(value, na.rm = TRUE), .groups = "drop")

strain_avg_df

# Find log10 fold changes using normalized values
lfc_df <- strain_avg_df %>%
  group_by(metabolite, lineage) %>%
  summarise(mean = mean(smean_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = lineage, values_from = mean) %>%
  mutate(
    L10FC_L1_Ref = L1 - Ref,
    L10FC_L2_Ref = L2 - Ref,
    L10FC_L4_Ref = L4 - Ref
  )

# Perform Welch's t-tests for each metabolite comparing each lineage to Ref
ttest_results <- annotated_df %>%
  filter(lineage %in% c("L1", "L2", "L4", "Ref")) %>%
  group_by(metabolite) %>%
  summarise(
    p_L1 = t.test(value[lineage == "L1"], value[lineage == "Ref"], var.equal = FALSE)$p.value,
    p_L2 = t.test(value[lineage == "L2"], value[lineage == "Ref"], var.equal = FALSE)$p.value,
    p_L4 = t.test(value[lineage == "L4"], value[lineage == "Ref"], var.equal = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    fdr_L1 = p.adjust(p_L1, method = "fdr"),
    fdr_L2 = p.adjust(p_L2, method = "fdr"),
    fdr_L4 = p.adjust(p_L4, method = "fdr")
  )

# Join LFC and p-value dfs
results <- lfc_df %>%
  left_join(ttest_results, by = "metabolite")

# Filter for significance (p < 0.05 and |L10FC| >= 0.3) in at least 1 group comp
sig_results <- results %>%
  filter(
    (fdr_L1 < 0.05 & abs(L10FC_L1_Ref) >= 0.3) |
      (fdr_L2 < 0.05 & abs(L10FC_L2_Ref) >= 0.3) |
      (fdr_L4 < 0.05 & abs(L10FC_L4_Ref) >= 0.3)
  )

################################################################################
### Step 2: Save a df of the metabolites with sig LFC (in at least 1 lineage)
################################################################################
lfc_filename <- './docu_data/Targeted/L10FC_sig_fdr_targeted_metabolites.csv'
write.csv(sig_results, lfc_filename, row.names = TRUE)

################################################################################
### Step 3: Display sig LFC values in a volcano plot for each lineage
################################################################################
# Filter to L1 df
L1_sig_results <- sig_results %>%
  mutate(neg_log10_fdr_L1 = -log10(fdr_L1)) %>% 
  select(metabolite, L10FC_L1_Ref, fdr_L1, neg_log10_fdr_L1) %>% 
  drop_na(fdr_L1)
# Create the L1 and Ref volcano plot with labels
library(ggplot2)
library(ggrepel)

# Define hex codes for each lineage
#lineage_colors <- c(
  #"L1" = "#DC143C",
  #"L2" = "#228B22",
  #"L4" = "#1E90FF"
#)

# Create the L1 and Ref volcano plot with labels only for significant points
L1_volcano_plot <- ggplot(L1_sig_results, aes(x = L10FC_L1_Ref, 
                                              y = neg_log10_fdr_L1)) +
  geom_point(aes(color = fdr_L1 < 0.05), alpha = 0.6) +
  geom_text_repel(
    data = subset(L1_sig_results, fdr_L1 < 0.05),
    aes(label = metabolite),
    size = 5,
    color = "#DC143C",
    max.overlaps = 25,
    family = "Times New Roman"
  ) +
  scale_color_manual(values = c("gray", "#DC143C")) +
  theme_minimal() +
  labs(title = "Lineage 1 Metabolites Compared to Reference H37Rv",
       x = "Log10 Fold Change",
       y = "-Log10(FDR-adj-p-value)",
       color = "Significant") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray60") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted", color = "gray60") +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  theme(
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
  )


# Print out the plot to examine
L1_volcano_plot

# Filter to L2 df
L2_sig_results <- sig_results %>%
  mutate(neg_log10_fdr_L2 = -log10(fdr_L2)) %>% 
  select(metabolite, L10FC_L2_Ref, fdr_L2, neg_log10_fdr_L2) %>% 
  drop_na(fdr_L2)
# Create the L2 and Ref volcano plot with labels
L2_volcano_plot <- ggplot(L2_sig_results, aes(x = L10FC_L2_Ref, 
                                              y = neg_log10_fdr_L2)) +
  geom_point(aes(color = fdr_L2 < 0.05), alpha = 0.6) +
  geom_text_repel(
    data = subset(L2_sig_results, fdr_L2 < 0.05),
    aes(label = metabolite),
    size = 5,
    color = "#228B22",
    max.overlaps = 25,
    family = "Times New Roman"
  ) +
  scale_color_manual(values = c("gray", "#228B22")) +
  theme_minimal() +
  labs(title = "Lineage 2 Metabolites Compared to Reference H37Rv",
       x = "Log10 Fold Change",
       y = "-Log10(FDR-adj-p-value)",
       color = "Significant") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray60") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted", color = "gray60") +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  theme(
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
  )

# Print out the plot to examine
L2_volcano_plot

# Filter to L4 df
L4_sig_results <- sig_results %>%
  mutate(neg_log10_fdr_L4 = -log10(fdr_L4)) %>% 
  select(metabolite, L10FC_L4_Ref, fdr_L4, neg_log10_fdr_L4) %>% 
  drop_na(fdr_L4)
# Create the L2 and Ref volcano plot with labels
L4_volcano_plot <- ggplot(L4_sig_results, aes(x = L10FC_L4_Ref, 
                                              y = neg_log10_fdr_L4)) +
  geom_point(aes(color = fdr_L4 < 0.05), alpha = 0.6) +
  geom_text_repel(
    data = subset(L4_sig_results, fdr_L4 < 0.05),
    aes(label = metabolite),
    size = 5,
    color = "#1E90FF",
    max.overlaps = 25,
    family = "Times New Roman"
  ) +
  scale_color_manual(values = c("gray", "#1E90FF")) +
  theme_minimal() +
  labs(title = "Lineage 4 Metabolites Compared to Reference H37Rv",
       x = "Log10 Fold Change",
       y = "-Log10(FDR-adj-p-value)",
       color = "Significant") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray60") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted", color = "gray60") +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  theme(
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
  )

# Print out the plot to examine
L4_volcano_plot

# Save all the volcano plots
ggsave("./docu_data/Targeted/L1_colored_metabolite_volcano.png", plot = L1_volcano_plot, 
       width = 12, height = 6, dpi = 300)
ggsave("./docu_data/Targeted/L2_colored_metabolite_volcano.png", plot = L2_volcano_plot, 
       width = 12, height = 6, dpi = 300)
ggsave("./docu_data/Targeted/L4_colored_metabolite_volcano.png", plot = L4_volcano_plot, 
       width = 12, height = 6, dpi = 300)

################################################################################
### Step 4: Cluster strains in a 2D PCA based on the metabolites
################################################################################
# Pivot wider so each metabolite becomes a column
strain_avg_wide_df <- strain_avg_df %>%
  pivot_wider(names_from = metabolite, values_from = smean_value) %>% 
  mutate(lineage = if_else(strain == "H37Rv", "L4", lineage))

# Remove label columns for PCA
strain_met_val_df <- strain_avg_wide_df %>% 
  select(-strain,-lineage)

# Perform PCA
pca_result <- prcomp(strain_met_val_df, scale. = TRUE)
# Extract loadings
loadings <- pca_result$rotation
# View loadings
print(loadings)

# Step 5: Create a data frame for plotting
pca_df <- data.frame(pca_result$x)

# Step 6: Add labels back to the data frame by their names
pca_df$Strain <- strain_avg_wide_df$strain
pca_df$Lineage <- strain_avg_wide_df$lineage

# Step 7: Calculate the percentage of variance explained
variance_explained <- summary(pca_result)$importance[3, ]
pc1_label <- paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)")
pc2_label <- paste0("PC2 (", round((variance_explained[2] - variance_explained[1]) * 100, 2), "%)")
pc3_label <- paste0("PC3 (", round((variance_explained[3] - variance_explained[2]) * 100, 2), "%)")

# Define hex codes for each lineage
lineage_colors <- c(
  "L1" = "#DC143C",
  "L2" = "#228B22",
  "L4" = "#1E90FF"
)

# Decide confidence intervals
# 95 good for statistical inference, 68 shows 1 sd
#stat_ellipse(type = "norm", level = 0.95)  # 95% CI (default)
#stat_ellipse(type = "norm", level = 0.68)  # 68% CI (1 SD)

# Make 2D PCA
pca_2d <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Lineage)) +
  geom_point() +
  geom_text(aes(label = Strain), vjust = -0.5, size = 7,
            family = "Times New Roman") +
  stat_ellipse(aes(group = Lineage), type = "norm", linetype = 2, level = 0.68) +
  labs(title = "PCA of Strain Metabolite Data by Lineage", x = pc1_label, y = pc2_label) +
  scale_color_manual(values = lineage_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    plot.title = element_text(size = 28),  # Change title font size
    axis.title = element_text(size = 28),  # Change axis titles font size
    axis.text = element_text(size = 28),   # Change axis text font size
    legend.text = element_text(size = 20), # Change legend text font size
    legend.title = element_text(size = 28) # Change legend title font size
  ) +
  guides(color = guide_legend(title = "Strain Group"))  # Customize legend title

pca_2d

#Save the 2D PCA
ggsave("./docu_data/Targeted/1sd_log10_mean_pareto_targeted_PCA.png", plot = pca_2d, 
       width = 10, height = 8, dpi = 300)

################################################################################
### Step 5: Pull out PC loadings and plot heatmap of PC metabolites
################################################################################
# Extract loadings
loadings <- pca_result$rotation

# View loadings
print(loadings)

# Calculate the absolute weights
abs_loadings <- abs(loadings)

# Initialize a list to store top 10 loadings dataframes for each PC
top_10_loadings_dfs <- list()

# Loop through each principal component to select top 10 loadings
for (pc in colnames(abs_loadings)) {
  # Rank the loadings for the current PC
  ranked_indices <- order(abs_loadings[, pc], decreasing = TRUE)
  
  # Select the top 10 loadings
  top_10_indices <- ranked_indices[1:10]
  
  # Create a dataframe for the top 10 loadings of the current PC
  top_10_loadings_df <- data.frame(Variable = rownames(loadings)[top_10_indices], 
                                   Loading = loadings[top_10_indices, pc])
  
  # Store the dataframe in the list
  top_10_loadings_dfs[[pc]] <- top_10_loadings_df
}

# Function to create bar plot for a given PC
create_bar_plot <- function(pc_loadings, pc_name) {
  # Melt the data frame for ggplot2
  loadings_melted <- melt(pc_loadings, id.vars = "Variable")
  names(loadings_melted)[c(1, 2, 3)] <- c("Metabolite", "PC", "Weight")
  
  # Filter to just one PC
  loadings_pc <- loadings_melted %>%
    filter(PC == "Loading")
  
  # Calculate and order by absolute weights
  loadings_pc$AbsWeight <- abs(loadings_pc$Weight)
  loadings_pc <- loadings_pc[order(-loadings_pc$AbsWeight), ]
  
  # Convert Metabolite to a factor with levels ordered by AbsWeight
  loadings_pc$Metabolite <- factor(loadings_pc$Metabolite, 
                                   levels = loadings_pc$Metabolite)
  
  # Plot the barplot
  ggplot(loadings_pc, aes(x = Metabolite, y = Weight, fill = Metabolite)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d() +
    labs(title = paste(pc_name, "Weights"), x = "Metabolites", y = "Weight") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    theme(axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.ticks.length = unit(0.1, "cm"))
}

# Create bar plots for PC1 and PC2
pc1_loadings <- top_10_loadings_dfs[["PC1"]]
pc2_loadings <- top_10_loadings_dfs[["PC2"]]
pc3_loadings <- top_10_loadings_dfs[["PC3"]]

pc1_plot <- create_bar_plot(pc1_loadings, "PC1")
pc2_plot <- create_bar_plot(pc2_loadings, "PC2")
pc3_plot <- create_bar_plot(pc3_loadings, "PC3")

# Print the plots
print(pc1_plot)
print(pc2_plot)
print(pc3_plot)

# Save the plots
ggsave("./docu_data/Targeted/pc1_top10_targeted_barplot.png", plot = pc1_plot, 
       width = 10, height = 5, dpi = 300)
ggsave("./docu_data/Targeted/pc2_top10_targeted_barplot.png", plot = pc2_plot, 
       width = 10, height = 5, dpi = 300)
ggsave("./docu_data/Targeted/pc3_top10_targeted_barplot.png", plot = pc3_plot, 
       width = 10, height = 5, dpi = 300)

################################################################################
### Step 6: Plot barplots for top LFC met from ref for L1, L2, and L4
################################################################################
# Filter for sig in L1 (p < 0.05 and |LFC| >= 1)
L1_sig_lfc <- results %>%
  filter((fdr_L1 < 0.05 & abs(L10FC_L1_Ref) >= 0.3)) %>% 
  select(metabolite, L10FC_L1_Ref, fdr_L1)

# Filter for sig in L2 (p < 0.05 and |LFC| >= 1)
L2_sig_lfc <- results %>%
  filter((fdr_L2 < 0.05 & abs(L10FC_L2_Ref) >= 0.3)) %>% 
  select(metabolite, L10FC_L2_Ref, fdr_L2)

# Filter for sig in L4 (p < 0.05 and |LFC| >= 1)
L4_sig_lfc <- results %>%
  filter((fdr_L4 < 0.05 & abs(L10FC_L4_Ref) >= 0.3)) %>% 
  select(metabolite, L10FC_L4_Ref, fdr_L4)

# L1 order met by abs(lfc)
L1_lfc_ordered <- L1_sig_lfc %>% 
  arrange(desc(abs(L10FC_L1_Ref)))
top_L1_lfc_ordered <- head(L1_lfc_ordered, 20)
top_L1_lfc_ordered$metabolite <- factor(
  top_L1_lfc_ordered$metabolite,
  levels = top_L1_lfc_ordered$metabolite[order(abs(top_L1_lfc_ordered$L10FC_L1_Ref))]
)

# L2 order met by abs(lfc)
L2_lfc_ordered <- L2_sig_lfc %>% 
  arrange(desc(abs(L10FC_L2_Ref)))
top_L2_lfc_ordered <- head(L2_lfc_ordered, 20)
top_L2_lfc_ordered$metabolite <- factor(
  top_L2_lfc_ordered$metabolite,
  levels = top_L2_lfc_ordered$metabolite[order(abs(top_L2_lfc_ordered$L10FC_L2_Ref))]
)

# L4 order met by abs(lfc)
L4_lfc_ordered <- L4_sig_lfc %>% 
  arrange(desc(abs(L10FC_L4_Ref)))
top_L4_lfc_ordered <- head(L4_lfc_ordered, 20)
top_L4_lfc_ordered$metabolite <- factor(
  top_L4_lfc_ordered$metabolite,
  levels = top_L4_lfc_ordered$metabolite[order(abs(top_L4_lfc_ordered$L10FC_L4_Ref))]
)

# L1 make the LFC barplot
L1_lfc_bar <- ggplot(top_L1_lfc_ordered, aes(x = metabolite, 
                                             y = L10FC_L1_Ref, 
                                             fill = fdr_L1)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#4B0082", high = "#9370DB", 
                      name = "FDR-adj p-value") + 
  labs(title = "Targeted Metabolites with Largest L1 L10FC from Ref", 
       x = "Metabolites", 
       y = "log10FC") +
  theme_minimal() +
  theme(text = element_text(family = "Times New Roman"),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.1, "cm")) +
  coord_flip()

# L2 make the LFC barplot
L2_lfc_bar <- ggplot(top_L2_lfc_ordered, aes(x = metabolite, 
                                             y = L10FC_L2_Ref, 
                                             fill = fdr_L2)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#4B0082", high = "#9370DB", 
                      name = "FDR-adj p-value") +
  labs(title = "Targeted Metabolites with Largest L2 L10FC from Ref", 
       x = "Metabolites", 
       y = "log10FC") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.1, "cm")) +
  coord_flip()

# L4 make the LFC barplot
L4_lfc_bar <- ggplot(top_L4_lfc_ordered, aes(x = metabolite, 
                                             y = L10FC_L4_Ref, 
                                             fill = fdr_L4)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#4B0082", high = "#9370DB", 
                      name = "FDR-adj p-value") +
  labs(title = "Targeted Metabolites with Largest L4 L10FC from Ref", 
       x = "Metabolites", 
       y = "log10FC") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.1, "cm")) +
  coord_flip()

# Print out the barplots to compare
L1_lfc_bar
L2_lfc_bar
L4_lfc_bar

# Save them too!
ggsave("./docu_data/Targeted/L1_L10FC_top10_targeted_barplot.png", plot = L1_lfc_bar, 
       width = 10, height = 5, dpi = 300)
ggsave("./docu_data/Targeted/L2_L10FC_top10_targeted_barplot.png", plot = L2_lfc_bar, 
       width = 10, height = 5, dpi = 300)
ggsave("./docu_data/Targeted/L4_L10FC_top10_targeted_barplot.png", plot = L4_lfc_bar, 
       width = 10, height = 5, dpi = 300)

################################################################################
### Step 7: Make all metabolite heatmap for lineage comparison
################################################################################
## Use any with sig LFC in one lineage
# Convert row names to a column (if needed)
h_table <- sig_results %>%
  select(metabolite, L10FC_L1_Ref, L10FC_L2_Ref, L10FC_L4_Ref)

# Pivot to long format
h_table_long <- h_table %>%
  pivot_longer(cols = -metabolite, names_to = "Category", values_to = "Value")

# Generate the heatmap
hmap <- ggplot(h_table_long, aes(x = Category, y = metabolite, fill = Value)) +
  geom_tile(color = "gray", size = 0.6) +
  geom_text(aes(label = sprintf("%.3f", Value)), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "#87CEEB", mid = "white", high = "#FF7F7F",
    midpoint = median(h_table_long$Value, na.rm = TRUE)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 18)
  ) +
  labs(title = "Heatmap with log2FC Values", x = "Lineage Group", 
       y = "Metabolite", fill = "")

print(hmap)

# Save the plot with large dimensions
ggsave("./docu_data/Targeted/tall_l10fc_sig_heatmap.png", plot = hmap, 
       width = 10, height = 40, dpi = 300)
