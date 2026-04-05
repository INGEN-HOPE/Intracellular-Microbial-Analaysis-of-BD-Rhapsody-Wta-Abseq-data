# Load necessary libraries
library(dplyr)
library(tidyverse)
library(rstatix) 
anno = read.csv("anno_ra_bd.csv")
norm_species = read.csv("~/Downloads/50_NORM_AMONG_4_GROUP.csv", row.names = 1)

species_long <- norm_species %>%
  pivot_longer(cols = -Species, 
               names_to = "Sample", 
               values_to = "Abundance")

# Merge with the annotation data to add condition information
species_annotated <- species_long %>%
  inner_join(anno, by = "Sample") 

run_significance_tests <- function(data) {
  # Get unique species
  unique_species <- unique(data$Species)
  
  pairwise_results <- data.frame(
    Species = character(),
    group1 = character(),
    group2 = character(),
    p.value = numeric(),
    p.adj = numeric(),
    significant = logical(),
    stringsAsFactors = FALSE
  )
  
  # Run tests for each species
  for (sp in unique_species) {
    # Subset data for current species
    sp_data <- data %>% filter(Species == sp)
    
    
    # If Kruskal-Wallis is significant, run pairwise Wilcoxon tests
      pw_test <- sp_data %>%
        pairwise_wilcox_test(
          Abundance ~ Condition,
          p.adjust.method = "BH"
        )
      
      if (nrow(pw_test) > 0) {
        pw_test$Species <- sp
        pairwise_results <- rbind(
          pairwise_results,
          pw_test %>% select(Species, group1, group2, p, p.adj) %>%
            dplyr::rename(p.value = p) %>%
            mutate(significant = p.adj < 0.05)
        )
      }
  }
  
  return(list(
    pairwise = pairwise_results %>% arrange(p.value)
  ))
}

# Run the tests
results <- run_significance_tests(species_annotated)

write.csv(results$pairwise, "E:/Single cell Microbes/bd/4_group_overlap_pairwise_results_norm.csv", row.names = FALSE)
