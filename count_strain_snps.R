# Step 1: load library
library(tidyverse)

# Step 2: Read the data 
full_data_file <- "./my_input/genetics_input/annot_full_genetics_data.csv"
full_data <- read_csv(full_data_file)

full_data <- full_data %>%
  drop_na(gene_rv_code)
nrow(full_data)

# Make gcolumn to use as counting variable
# Make ref allele, alt allele, and gene accession (include nonsyn & syn mut)
full_data <- full_data %>%
  drop_na(gene_rv_code) %>% 
  mutate(variable = paste(reference_allele, alternate_allele, 
                          gene_rv_code, sep = "_"))

# Reformat data to subset to just gene name + mut, and strains with alt allele
data <- full_data %>%
  select(variable, strains_with_alt_allele) %>% 
  rename(group = strains_with_alt_allele)

# Step 3: Split the group column into individual rows and remove white spaces
data <- data %>%
  separate_rows(group, sep = ",") %>%
  mutate(group = str_trim(group))

# Step 4: Remove certain groups by name
groups_to_remove <- c("333", "355", "367", "411")  # Replace with actual group names to remove
data <- data %>%
  filter(!group %in% groups_to_remove)

# Step 5: Summarize the number of unique variables in each group
group_summary <- data %>%
  group_by(group) %>%
  summarise(count = n_distinct(variable))

print(group_summary)

# Step 6: Compare each pair of groups
compare_groups <- function(group1, group2, data) {
  vars_group1 <- data %>% filter(group == group1) %>% pull(variable)
  vars_group2 <- data %>% filter(group == group2) %>% pull(variable)
  
  shared_vars <- intersect(vars_group1, vars_group2)
  diff_vars_group1 <- setdiff(vars_group1, vars_group2)
  diff_vars_group2 <- setdiff(vars_group2, vars_group1)
  
  tibble(
    group1 = group1,
    group2 = group2,
    shared_count = length(shared_vars),
    diff_count_group1 = length(diff_vars_group1),
    diff_count_group2 = length(diff_vars_group2),
    total_diff_count = length(diff_vars_group1) + length(diff_vars_group2)
  )
}

# Get all unique groups
groups <- unique(data$group)

# Compare each pair of groups without duplication
comparison_results <- combn(groups, 2, function(pair) {
  compare_groups(pair[1], pair[2], data)
}, simplify = FALSE) %>%
  bind_rows()

print(comparison_results)

# Save as a csv to make a distance matrix
output_filename <- './14_strains_all_snps_counts.csv'
write.csv(comparison_results, output_filename, row.names = FALSE)
