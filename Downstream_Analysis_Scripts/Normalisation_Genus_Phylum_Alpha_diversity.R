#load library
library(pavian)
library(metagenomeSeq)
library(dplyr)
library(tidyr)
library(tibble) 
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(rstatix)

#load otu matrix from pavian and metadata
otu_matrix = read.csv("~/Downloads/raw_reads_counts.csv")
otu_matrix[, 8:ncol(otu_matrix)][is.na(otu_matrix[, 8:ncol(otu_matrix)])] <- 0
otu_matrix_reads <- otu_matrix[ ,c(8:ncol(otu_matrix))]

metadata_df <-read.csv("~/Downloads/Supplementary Table S1.xlsx - CLINICAL INFO & METADATA.csv", header = TRUE, skip = 2) %>% column_to_rownames("Sample.ID")

#Normalization
metaSeqObject = newMRexperiment(otu_matrix_reads)

metaSeqObject_CSS = cumNorm( metaSeqObject , p = cumNormStatFast(metaSeqObject) )

read_count_CSS = data.frame(MRcounts(metaSeqObject_CSS, norm = TRUE, log = TRUE))

read_count_CSS <- cbind(otu_matrix[, 1:7], read_count_CSS)

#write normalised reads
write.csv(read_count_CSS, "E:/Single cell Microbes/bd/NORMALIZED__raw_READS_sc_MICROBES.csv", row.names = FALSE)

#create phyloseq object
otu_mat <- read_count_CSS[ ,c(8:ncol(read_count_CSS))]
otu_mat = round(otu_mat)
rownames(otu_mat) <- read_count_CSS[, 7] 
tax_mat <- as.matrix(read_count_CSS[, 1:7])
rownames(tax_mat) <- read_count_CSS[, 7]
all(rownames(metadata_df) %in% colnames(otu_mat))

ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(metadata_df)
)

#write normalised otu matrix
write.csv(otu_table(otu_mat), "otu_matrix_bd.csv")

#calculate relative abundance
ps_relabund <- transform_sample_counts(ps, function(x) x / sum(x))
relabund_df <- cbind(as.data.frame(tax_table(ps_relabund)), as.data.frame(otu_table(ps_relabund)))

#write relative abundance sheet
write.csv(relabund_df, "Relative_Abundance_otu_matrix.csv")

# Summing by Genus per sample
genus_sum <- relabund_df %>%
  select(Genus, starts_with("b")) %>%
  pivot_longer(cols = starts_with("b"), names_to = "Sample", values_to = "Value") %>%
  group_by(Sample, Genus) %>%
  summarise(Total_Value = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_Value)

# Summing by phylum per sample
phylum_sum <- relabund_df %>%
  select(Phylum, starts_with("b")) %>%
  pivot_longer(cols = starts_with("b"), names_to = "Sample", values_to = "Value") %>%
  group_by(Sample, Phylum) %>%
  summarise(Total_Value = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_Value)

#write the summarized files of RA at Genus and Phylum level
write.csv(genus_sum, "GENUS_SUM_RA_BD.csv", row.names = FALSE)

write.csv(phylum_sum, "PHYLUM_SUM_RA_BD.csv", row.names = FALSE)

#calculate alpha diversity
alpha_diversity <- estimate_richness(ps)
alpha_diversity = merge(alpha_diversity, metadata_df[, 7, drop =FALSE], by = "row.names", drop=FALSE) %>% column_to_rownames("Row.names")

#write alpha diversity indices sheet
write.csv(data.alpha, "E:/Single cell Microbes/Alpha_diversity.csv")

# check statistical significance of alpha diversity indices among experimental groups 
stat_test <- alpha_diversity %>%
  wilcox_test(Shannon ~ Group) %>%
  add_significance(p.col = "p")

#select colours for visualization
pastel_colors <- c("Healthy" = "#6f917c", "Infected" = "#ff983e", "Recovered" = "#a793ac", "Long Recovered" = "#4EA6FF")

# Visualize Alpha diversity with their significance
p <- ggplot(data = alpha_diversity, ) +
  geom_violin(aes(x = Group, y = Shannon, color = Group), fill = NA, trim = FALSE, alpha = 0.7, linewidth = 1) +
  geom_boxplot(aes(x = Group, y = Shannon, color = Group), fill = NA, width = 0.2, alpha = 0.5, linewidth = 1, outliers = FALSE) +
  geom_jitter(aes(x = Group, y = Shannon, color = Group, , fill = Group), , shape = 21, width = 0.3, size = 1.8, alpha = 0.5) +
  scale_fill_manual(values = pastel_colors) +
  scale_color_manual(values = pastel_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"),  # Bold x-axis labels (Normal, Tumor)
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.title = element_text(size = 14, face = "bold", colour = "black"),
    axis.text = element_text(size = 14)
  ) +
  labs(x = "Group", y = "Shannon") +
  # Set order of Groups
  scale_x_discrete(limits = c("Healthy", "Infected", "Recovered", "Longitudinal Recovered")) +
  scale_y_continuous(breaks = seq(6, 8, by = 1),expand = expansion(mult = c(0, 0.05)), limits = c(6, 8))

# Set different heights for each comparison
y_pos <- c(7.3, 7.48, 7.57, 7.68, 7.79, 7.9)  

# Create a data frame for p-value annotations
p_values_alpha_diversity <- stat_test %>%
  add_xy_position(x = "Group", fun = "max") %>%
  mutate(y.position = y_pos)

# Add statistical annotation
p <- p + stat_pvalue_manual(
  p_values_alpha_diversity,
  label = "p.signif",
  size = 5,
  bracket.size = 0.8,
  hide.ns = TRUE,
  tip.length = 0.01,
  bracket.nudge.y = 0.07
) +
  theme(
    
  )

# Visualise
print(p)

# Save the plots
ggsave("E:/Single cell Microbes/FIGURES/FIG_NEW/Shannon_SC_MICROBES.svg", p, width = 10, height = 6, dpi = 600, bg = "white")
ggsave("E:/Single cell Microbes/FIGURES/FIG_NEW/Shannon_SC_MICROBES.png", p, width = 10, height = 6, dpi = 600, bg = "white")
ggsave("E:/Single cell Microbes/FIGURES/FIG_NEW/Shannon_SC_MICROBES.pdf", p, width = 10, height = 6, dpi = 600, bg = "white")