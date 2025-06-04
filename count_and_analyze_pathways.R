##############################################################################
# Step 1: Examine the pathways of all metabolites (w/ sig L10FC from ref)
##############################################################################
# Load libraries
library(tidyverse)

##############
# Set font
library(showtext)
font_add(family = "Times New Roman", 
         regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf")
showtext_auto()
showtext_opts(dpi = 300)
##############

# Load df
pathways_df_filename <- 'targeted_L10FC_sig_fdr_KEGG_metabolite_pathways.csv'
#pathways_df_filename <- 'untargeted_L10FC_sig_fdr_KEGG_metabolite_pathways.csv'
pathways_df <- read_csv(pathways_df_filename)

# Make small df, remove overview classes
pathways_only_df <- pathways_df %>% 
  select(Metabolite, KEGG_ID, KEGG_Class) %>% 
  distinct() %>% 
  filter(!str_detect(KEGG_Class, "Overview maps")) %>% 
  filter(!str_detect(KEGG_Class, "Global and overview maps"))

# Count metabolites per pathway
pathway_counts <- pathways_only_df %>%
  count(KEGG_Class, name = "metabolite_count")
  

# Plot
class_barplot <- ggplot(pathway_counts, aes(x = metabolite_count, 
                                            y = reorder(KEGG_Class, metabolite_count))) +
  geom_col(fill = "skyblue") +
  labs(
    title = "Untargeted Metabolite Count per Pathway Class",
    x = "Metabolite Count",
    y = "Pathway Class"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
    )

class_barplot

# Save the plot
plot_name <- "./docu_data/Targeted/metabolite_class_counts_targeted.png"
plot_name <- "./docu_data/Untargeted/metabolite_class_counts_untargeted.png"
ggsave(plot_name, 
       plot = class_barplot, 
       width = 12, height = 8, dpi = 300)

#'./docu_data/Targeted/metabolite_class_counts_targeted.png'

##############################################################################
# Step 2: Examine the pathways of the SNPS w/ matched metabolite(s)
##############################################################################
library(tidyverse)
# Load df
snp_matched_df_filename <- './docu_data/Targeted/targeted_KEGG_metabolite_matched_lineage_SNPs.csv'
#snp_matched_df_filename <- './docu_data/Untargeted/untargeted_KEGG_metabolite_matched_lineage_SNPs.csv'
snp_matched_df <- read_csv(snp_matched_df_filename)

# Subset and rename columns
snp_small_df <- snp_matched_df %>%
  unite(col = SNP, KEGG_Enzyme, gene_name, amino_acid_change, sep = "_", remove = FALSE) %>% 
  rename(Lineage_Group = "lineages_with_alt_allele") %>% 
  select(SNP, KEGG_Class, Lineage_Group) %>% 
  distinct() %>% 
  separate_rows(Lineage_Group, sep = ", ")

# Count SNPs per pathway class and group
snp_counts <- snp_small_df %>%
  count(KEGG_Class, Lineage_Group, name = "snp_count")

# Custom colors for the three groups
lineage_colors <- c(
  "L1" = "#DC143C",
  "L2" = "#228B22",
  "L4" = "#1E90FF"
)

# Ensure all combinations of KEGG_Class and Lineage_Group are present
snp_counts_complete <- snp_counts %>%
  complete(KEGG_Class, Lineage_Group, fill = list(snp_count = 0))

# Get unique KEGG classes in plotting order
class_levels <- levels(factor(snp_counts_complete$KEGG_Class))

# Plot with black lines between classes
snp_barplot <- ggplot(snp_counts_complete, aes(x = snp_count, 
                                               y = factor(KEGG_Class, levels = class_levels), 
                                               fill = Lineage_Group)) +
  geom_col(position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = lineage_colors) +
  labs(
    title = "Targeted Metabolite-Related SNPs per Pathway Class by Lineage",
    x = "Metabolite Related SNPs",
    y = "Pathway Class",
    fill = "Lineage"
  ) +
  theme_minimal() +
  geom_hline(yintercept = seq(1.5, length(class_levels) - 0.5, by = 1), color = "black", linewidth = 0.3)

snp_barplot


# Save the plot
snp_plot_path <- "./docu_data/Targeted/matched_snp_pathway_counts_targeted.png"
#snp_plot_path <- "./docu_data/Untargeted/untargeted_matched_snp_pathway_counts_targeted.png"
ggsave(snp_plot_path, 
       plot = snp_barplot, 
       width = 8, height = 10, dpi = 300)

##############################################################################
# Step 3: Make L10FC metabolite heatmaps for each pathway of interest
##############################################################################
# Load libraries
library(tidyverse)

##############
# Set font
library(showtext)
font_add(family = "Times New Roman", 
         regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf")
showtext_auto()
showtext_opts(dpi = 300)  # Adjust this value as needed otherwise font will be tiny
##############

# Load df
pathways_df_filename <- 'targeted_L10FC_sig_fdr_KEGG_metabolite_pathways.csv'
#pathways_df_filename <- 'untargeted_L10FC_sig_fdr_KEGG_metabolite_pathways.csv'

pathways_df <- read_csv(pathways_df_filename)

# Add FC values
fc_pathways_df <- pathways_df %>%
  mutate(FC_L1_Ref = sign(L10FC_L1_Ref) * 10^abs(L10FC_L1_Ref)) %>% 
  mutate(FC_L2_Ref = sign(L10FC_L2_Ref) * 10^abs(L10FC_L2_Ref)) %>% 
  mutate(FC_L4_Ref = sign(L10FC_L4_Ref) * 10^abs(L10FC_L4_Ref))

# Subset to include only data columns of interest
small_fc_pathways_df <- fc_pathways_df %>% 
  select(Metabolite, FC_L1_Ref, FC_L2_Ref, FC_L4_Ref, 
         KEGG_ID, KEGG_Class) %>%
  unite(col = Metabolite_ID, Metabolite, KEGG_ID, sep = "_", remove = TRUE) %>% 
  distinct()

# Make geat map with auto set colors - same across all plots and made for easy viz

# Set the KEGG class to analyze
selected_kegg_class <- "Amino acid metabolism" # change as needed

# Subset to include only pathway class related rows
subset_fc_pathways_df <- small_fc_pathways_df %>%
  filter(KEGG_Class == selected_kegg_class)

# Format data
# Convert row names to a column (if needed)
h_table <- subset_fc_pathways_df %>%
  select(Metabolite_ID, FC_L1_Ref, FC_L2_Ref, FC_L4_Ref)
# Pivot to long format
h_table_long <- h_table %>%
  pivot_longer(cols = -Metabolite_ID, names_to = "Category", values_to = "Value")

# Define bins and corresponding colors
h_table_long <- h_table_long %>%
  mutate(Bin = case_when(
    Value <= -30 ~ "< -30",
    Value <= -20 ~ "-30 to -20",
    Value <= -10 ~ "-20 to -10",
    Value <= -5 ~ "-10 to -5",
    Value <= -2 ~ "-5 to -2",
    Value <= 0 ~ "-2 to 0",
    Value <= 2 ~ "0 to 2",
    Value <= 5 ~ "2 to 5",
    Value <= 10 ~ "5 to 10",
    Value <= 20 ~ "10 to 20",
    Value <= 30 ~ "20 to 30",
    TRUE ~ "> 30"
  ))

# Define custom colors for each bin
bin_colors <- c(
  "< -30" = "#08306B",
  "-30 to -20" = "#08519C",
  "-20 to -10" = "#2171B5",
  "-10 to -5" = "#4292C6",
  "-5 to -2" = "#9ECAE1",
  "-2 to 0" = "#FFFFFF",
  "0 to 2" = "#FFFFFF",
  "2 to 5" = "#FCBBA1",
  "5 to 10" = "#FC9272",
  "10 to 20" = "#FB6A4A",
  "20 to 30" = "#EF3B2C",
  "> 30" = "#CB181D"
)

# Set text color for contrast
h_table_long <- h_table_long %>%
  mutate(TextColor = case_when(
    Bin %in% c("< -30", "-30 to -20", "-20 to -10", 
               "10 to 5", "20 to 30", "> 30") ~ "white",
    TRUE ~ "black"
  ))

# Generate the new heatmap
new_hmap <- ggplot(h_table_long, aes(x = Category, y = Metabolite_ID, 
                                     fill = Bin)) +
  geom_tile(color = "gray", size = 0.6) +
  geom_text(aes(label = sprintf("%.3f", Value), color = TextColor), size = 7,
            family = "Times New Roman") +
  scale_color_identity() +
  scale_fill_manual(values = bin_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 22),
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"), 
    axis.text.y = element_text(color = "black"),                      
    axis.title = element_text(size = 20, color = "black"),              
    axis.text = element_text(size = 20),                                 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20)
  ) +
  labs(title = paste(selected_kegg_class, "heatmap"), 
       x = "Lineage", y = "Metabolite", fill = "")

print(new_hmap)

# Save the plot with large dimensions
max_char_length <- max(nchar(h_table$Metabolite_ID))
nwidth <- max_char_length * 0.3 + 4.8  # Add padding for labels and margins
nheight = nrow(h_table)*0.5 + 1.5
# Create a safe file name
file_safe_class <- gsub(" ", "_", tolower(selected_kegg_class))
hmap_plotpath <- paste0("color_", file_safe_class, "_fc_heatmap_targeted.png")
ggsave(hmap_plotpath, plot = new_hmap, 
       width = nwidth, height = nheight, dpi = 300)

