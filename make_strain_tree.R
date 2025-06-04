# Step 0: Load libraries
library(ape)
library(ggtree)
library(dplyr)
library(ggplot2)

# Step 1: Read in the distance matrix for the 14 strains
output_filename <- 'strain_tree_distance_matrix.csv'
dist_matrix <- as.matrix(output_filename)

# Step 3: Create phylogenetic tree using Neighbor-Joining method
tree <- nj(dist_matrix)
phylo_tree <- as.phylo(tree)

# Step 4: Extract strain names
strain_names <- rownames(dist_matrix)

# Step 5: Extract group (prefix) and label (suffix)
group_prefix <- sub("_.*", "", strain_names)
label_suffix <- sub(".*_", "", strain_names)

# Step 6: Create a named vector for grouping
names(group_prefix) <- strain_names

#Step 7: Group the tree using groupOTU
grouped_tree <- groupOTU(phylo_tree, 
                         split(strain_names, group_prefix), 
                         group_name = "lineage")

# Step 8: Create a data frame for tip labels
group_data <- data.frame(
  label = label_suffix,
  original = strain_names
)
rownames(group_data) <- group_data$original

# Step 9: Define custom colors
lineage_colors <- c(
  "L1" = "#DC143C",  # Crimson
  "L2" = "#228B22",  # Forest Green
  "L4" = "#1E90FF"   # Dodger Blue
)

# Step 10: Plot the tree with suffix labels
c_tree <- ggtree(grouped_tree, layout = "circular", aes(color = lineage)) %<+% group_data +
  geom_tiplab(aes(label = label), size = 6, hjust = -0.1) +
  scale_color_manual(values = lineage_colors) +
  theme(legend.position = "none", 
        text = element_text(family = "Times New Roman"))

# Display the tree
print(c_tree)

# Save the plot
ggsave("circular_lineage_tree_colored.png", 
       plot = c_tree, width = 6, height = 6, units = "in", dpi = 300)
