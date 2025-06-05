# this script joins the metabolite data to the genetic data using KEGG information
###############################################################################
# Step 0: Load necessary libraries
library(tidyverse)
library(KEGGREST)
library(curl)

###############################################################################
# Step 1: Load input data files and reformat

# Load targeted analysis metabolites
df1_file <- './docu_data/Targeted/L10FC_sig_fdr_targeted_metabolites.csv'
df1 <- read_csv(df1_file)
names(df1)[1] <- "Rank"
names(df1)[2] <- "Metabolite"
df1 <- df1 %>%
  mutate(Metabolite = str_replace_all(Metabolite, "-", " "))
name_file <- "./my_input/reference/alternate_metabolite_names.csv"
name_df <- read_csv(name_file)
df1 <- df1 %>%
  left_join(name_df, by = c("Metabolite" = "Query")) %>%
  distinct()

# Load genetics data
df2_file <- "./my_input/genetics_input/annot_15_strains_genetics_data.csv"
df2 <- read_csv(df2_file)
selected_df2 <- df2 %>%
  select(gene_rv_code, gene_name, amino_acid_change, start_gene_pos, end_gene_pos,
         syn_or_non_mutation, protein_annotation, essential_or_not,
         strains_with_ref_allele, strains_with_alt_allele,
         lineages_with_ref_allele, lineages_with_alt_allele) %>%
  filter(syn_or_non_mutation == "nonsynonymous")

# Untargeted analysis (commented out on purpose)
# df1_file <- './docu_data/Untargeted/L10FC_sig_untargeted_metabolites.csv'
# df1 <- read_csv(df1_file)
# names(df1)[1] <- "Rank"
# names(df1)[2] <- "Metabolite_ID"
# df1 <- df1 %>% 
#   separate(Metabolite_ID, c("Metabolite", "KEGG_ID"), sep = "_")
# 
# df2_file <- "./my_input/genetics_input/annot_15_strains_genetics_data.csv"
# df2 <- read_csv(df2_file)
# selected_df2 <- df2 %>%
#   select(gene_rv_code, gene_name, amino_acid_change, start_gene_pos, end_gene_pos,
#          syn_or_non_mutation, protein_annotation, essential_or_not, 
#          strains_with_ref_allele, strains_with_alt_allele,
#          lineages_with_ref_allele, lineages_with_alt_allele) %>%
#   filter(syn_or_non_mutation == "nonsynonymous") 

###############################################################################
# Step 2: Pull pathways related to each KEGG compound
get_pathways_from_kegg <- function(kegg_compound) {
  pathways <- tryCatch({
    keggLink("pathway", kegg_compound)
  }, error = function(e) {
    print(paste("Error fetching pathways for:", kegg_compound, " - ", e$message))
    return(NULL)
  })
  
  if (is.null(pathways)) {
    return(NA)
  }
  
  print(paste("KEGG Compound:", kegg_compound))
  print(paste("Pathways returned:", pathways))
  
  pathway_values <- pathways
  extracted_codes <- sapply(pathway_values, function(value) {
    if (grepl("path:", value)) {
      code <- substr(value, 6, 13)  # Extract the 8 characters after "path:"
      mtu_code <- paste0("mtu", substr(code, 4, 8))  # Convert map to mtu
      # Check if the mtu pathway exists
      mtu_pathway <- tryCatch({
        keggGet(mtu_code)
      }, error = function(e) {
        print(paste("Error fetching mtu pathway for:", mtu_code, " - ", e$message))
        return(NULL)
      })
      if (!is.null(mtu_pathway)) {
        return(mtu_code)
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  })
  
  print(paste("Extracted codes:", extracted_codes))
  
  extracted_codes <- extracted_codes[!is.na(extracted_codes)]
  return(paste(extracted_codes, collapse = ","))
}
df1$KEGG_Pathways <- sapply(df1$KEGG_ID, get_pathways_from_kegg)

###############################################################################
# Step 3: Unnest pathways
df1_unnested_pathways <- df1 %>%
  separate_rows(KEGG_Pathways, sep = ",")

df1_unnested_pathways <- df1_unnested_pathways %>% 
  filter(KEGG_Pathways != "")

# Check unique pathways
unique_pathways <- df1_unnested_pathways %>%
  distinct(KEGG_Pathways)
print(nrow(unique_pathways)) # 79 or #62

# Join in the KEGG_pathway_names and classes
pathway_classes_df_file <- './docu_data/KEGG_pathway_names.csv'
pathway_classes_df <- read_csv(pathway_classes_df_file)
df1_labeled_pathways <- df1_unnested_pathways %>%
  left_join(pathway_classes_df, by = c("KEGG_Pathways" = "KEGG_Pathway_Codes")) %>% 
  drop_na()

checkpoint_filename <- './docu_data/Targeted/targeted_L10FC_sig_fdr_KEGG_metabolite_pathways.csv'
#checkpoint_filename <- './docu_data/Untargeted/untargeted_L10FC_sig_fdr_KEGG_metabolite_pathways.csv'

if (nrow(df1_labeled_pathways) >= 1) {
  write.csv(df1_labeled_pathways, checkpoint_filename, row.names = TRUE)
}

###############################################################################
# Reload df if necessary
reload_file <- './docu_data/Targeted/targeted_L10FC_sig_fdr_KEGG_metabolite_pathways.csv'
#reload_file <- './docu_data/Untargeted/untargeted_L10FC_sig_fdr_KEGG_metabolite_pathways.csv'
df1_unnested_pathways <- read_csv(reload_file)

# Step 4: Pull pathway enzymes
# Define function to get genes from KEGG (test with prints)
get_genes_from_kegg <- function(pathway_code) {
  print(paste("Processing pathway code:", pathway_code))
  
  # Retrieve pathway information from KEGG
  pathway_info <- keggGet(pathway_code)
  print("Retrieved pathway information")
  
  # Check if the pathway information is not empty and contains gene data
  if (length(pathway_info) > 0 && !is.null(pathway_info[[1]]$GENE)) {
    print("Pathway contains gene data")
    
    # Extract the gene information
    genes <- pathway_info[[1]]$GENE
    print(paste("Extracted genes:", genes))
    
    # Filter gene codes that start with "Rv" followed by 4 numbers and maybe a letter
    rv_codes <- grep("^Rv[0-9]{4}[A-Za-z]?$", genes, value = TRUE)
    print(paste("Filtered Rv codes:", rv_codes))
    
    # Return the Rv codes as a comma-separated list
    result <- paste(rv_codes, collapse = ",")
    print(paste("Final result:", result))
    return(result)
  } else {
    print("No gene data found")
    # Return NA if no gene information is found
    return(NA)
  }
}

#df1_small_test <- head(df1_unnested_pathways,3)
#df1_small_test$KEGG_Genes <- sapply(df1_small_test$KEGG_Pathways, 
                                    #get_genes_from_kegg)
#print(df1_small_test)

# Apply the function to the dataframe
df1_unnested_pathways$KEGG_Genes <- sapply(df1_unnested_pathways$KEGG_Pathways, 
                                           get_genes_from_kegg)

###############################################################################
# Step 5: Unnest enzymes
df1_genes_unnested <- df1_unnested_pathways %>%
  separate_rows(KEGG_Genes, sep = ",")
df1_genes_unnested <- df1_genes_unnested %>% 
  rename(KEGG_Enzyme = KEGG_Genes) %>% 
  rename(KEGG_Pathway = KEGG_Pathways)

# Save the metabolite - enzyme matched file?

# Load in biocyc info 
#filtered_df <- read_csv("./docu_data/all_snps_matched_biocyc.csv")
#df1_genes_unnested <- filtered_df

###############################################################################
# Step 6: Perform the join on the rv value column
joined_df <- df1_genes_unnested %>%
  left_join(selected_df2, by = c("KEGG_Enzyme" = "gene_rv_code")) %>%
  drop_na(amino_acid_change)

###############################################################################
# Step 7: Subset further to avoid strain-specific snps
df_final <- joined_df %>% 
  filter(nchar(strains_with_alt_allele) > 3)

###############################################################################
# Step 8: Save the path to SNP csv if the dataframe has more than 1 row
output_filename <- 'targeted_KEGG_metabolite_matched_lineage_SNPs.csv'
if (nrow(df_final) >= 1) {
  write.csv(df_final, output_filename, row.names = TRUE)
}
