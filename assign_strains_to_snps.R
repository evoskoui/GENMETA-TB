# This script takes the genetics data spreadsheet and assigns strains to SNPs
# essentially, it collapses the SNP by strain columns into one column that is a list of strains with the mutant SNP
# then it goes through and makes another column that lists the lineages
# for this, it needs a metadata file with the lineage - strain matches

# Step 0: Import necessary libraries
library(tidyverse)

# Step 1: Read in the input data file(s) (csv format)
metadata_file <- "strain_lineage_metadata.csv"
genetic_data_file <- "15_strains_genetics_data.csv"
metadata <- read_csv(metadata_file)
gendata <- read_csv(genetic_data_file)

# Step 2: Make function for searching alleles per row
find_columns_with_dynamic_value <- function(data, value_column, col_limit = 14) {
  data %>%
    rowwise() %>%
    mutate(columns_with_value = paste(names(data)[which(c_across(1:col_limit) == get(value_column))], collapse = ", ")) %>%
    ungroup() %>%
    select(columns_with_value)
}

# Step 3: Apply function and get mutant_allele lists per row
ref_result <- find_columns_with_dynamic_value(gendata, "reference_allele")

# Step 4: Apply function and get reference_allele lists per row
alt_result <- find_columns_with_dynamic_value(gendata, "alternate_allele")

# Step 5: Add the strain allele lists as new columns
annot_gendata <- bind_cols(gendata, ref_result) %>%
  rename("strains_with_ref_allele" = columns_with_value)
annot_gendata <- bind_cols(annot_gendata, alt_result) %>%
  rename("strains_with_alt_allele" = columns_with_value)

# Step 6: Make metadata table into a dict
meta_dict <- metadata %>%
  deframe()
print(meta_dict)

# Step 7: Add lineages as well
## Function to map semicolon-separated values
map_values <- function(val, dict) {
  keys <- unlist(strsplit(val, ",\\s*"))  # Split string by semicolon and remove spaces
  mapped <- unique(dict[keys])  # Map using the dictionary
  mapped[is.na(mapped)] <- "Unknown"  # Handle missing values
  return(paste(mapped, collapse = ", "))  # Rejoin into a single string
}
## Apply this map_values function to annotate using lapply
new_annot_gendata <- annot_gendata %>%
  mutate(lineages_with_ref_allele = sapply(annot_gendata$strains_with_ref_allele, 
                                      map_values, dict = meta_dict))
new_annot_gendata <- new_annot_gendata %>%
  mutate(lineages_with_alt_allele = sapply(annot_gendata$strains_with_alt_allele, 
                                      map_values, dict = meta_dict))
## View the result
print(new_annot_gendata)

# Step 8: Save the new csv
write.csv(new_annot_gendata, "annot_15_strains_genetics_data.csv", row.names = TRUE)
