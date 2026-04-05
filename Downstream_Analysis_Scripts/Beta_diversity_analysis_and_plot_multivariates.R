library(vegan)
library(phyloseq)
library(ggplot2)
library(ggplot2)
library(tibble) 
library(patchwork)
library(svglite)

otu_matrix = read.csv("~/Downloads/otu_matrix_bd.csv", row.names = 1)
beta_dist <- vegdist(t(otu_matrix), method = "bray")

meta = read.csv("~/Downloads/Supplementary Table S1.xlsx - CLINICAL INFO & METADATA.csv", header = TRUE, skip = 2) %>% column_to_rownames("Sample.ID")

meta_ordered <- meta[colnames(otu_matrix), ]
all(rownames(meta_ordered) == colnames(otu_matrix))

# multivariate statistical test accounting for infection and age
permanova_multi <- adonis2(beta_dist ~ Group + Age, 
                           data = meta_ordered, 
                           permutations = 999, by = "margin")

pcoa <- cmdscale(beta_dist, eig = TRUE)

pcoa_df <- data.frame(
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
)


pcoa_df = merge(pcoa_df, meta[, 7, drop=FALSE], by = "row.names") %>% column_to_rownames("Row.names")

centroids <- aggregate(cbind(PC1, PC2) ~ Group, pcoa_df, mean) 

pcoa_df$Group <- factor(pcoa_df$Group, levels = c("Healthy", "Infected", "Recovered", "Long Recovered"))

pastel_colors <- c("Healthy" = "#6f917c", "Infected" = "#ff983e", "Recovered" = "#a793ac", "Long Recovered" = "#4EA6FF")


stats_label <- paste0("Infection: R2 = ", round(permanova_multi["Group", "R2"], 3) , ", p = ", permanova_multi["Group", "Pr(>F)"], "\n",
                      "Age: R2 = ", round(permanova_multi["Age", "R2"], 3), ", p = ", round(permanova_multi["Age", "Pr(>F)"], 2),"\n"
)

#visualisation of beta plot with statistical significance level
beta_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, level = 0.95, show.legend = FALSE) +
  theme_minimal() +
  scale_color_manual(values = pastel_colors) + 
  scale_fill_manual(values = pastel_colors) + 
  labs(
    title = "Beta Diversity (Adjusted for Age)",
    x = paste0("PC1 (", round(100 * pcoa$eig[1] / sum(pcoa$eig), 2), "%)"),
    y = paste0("PC2 (", round(100 * pcoa$eig[2] / sum(pcoa$eig), 2), "%)")
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    plot.title = element_blank(),
    legend.position = "none" 
  ) +
  annotate("text",
           x = min(pcoa_df$PC1), 
           y = max(pcoa_df$PC2), 
           label = stats_label,
           hjust = 1.3,
           vjust = 1,    
           size = 4, 
           color = "black", 
           fontface = "bold")

xdensity <- pcoa_df %>%
  
  ggplot(aes(x = PC1)) +
  
  geom_density(alpha = 0.5, aes(fill = Group, color = Group)) +
  
  scale_fill_manual(values = c("Healthy" = "#6f917c", "Infected" = "#ff983e", "Recovered" = "#a793ac", "Long Recovered" = "#4EA6FF")) + # Set manual fill colors
  
  scale_color_manual(values = c("Healthy" = "#6f917c", "Infected" = "#ff983e", "Recovered" = "#a793ac", "Long Recovered" = "#4EA6FF")) + # Set manual line colors
  
  theme_classic() + 
  
  theme(
    
    axis.title.x = element_blank(), 
    
    axis.text.x = element_blank(), 
    
    axis.ticks.x = element_blank(), 
    
    axis.line.x = element_blank(), 
    
    axis.text.y = element_blank(), 
    
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none" 
    
  )

ydensity <- pcoa_df %>%
  
  ggplot(aes(PC2)) +
  
  geom_density(alpha = 0.5, aes(fill = Group, color = Group)) + 
  
  scale_fill_manual(values = c("Healthy" = "#6f917c", "Infected" = "#ff983e", "Recovered" = "#a793ac", "Long Recovered" = "#4EA6FF")) + # Set manual fill colors
  
  scale_color_manual(values = c("Healthy" = "#6f917c", "Infected" = "#ff983e", "Recovered" = "#a793ac", "Long Recovered" = "#4EA6FF")) + # Set manual line colors
  
  theme_classic() + 
  
  theme(
    
    axis.text.x = element_blank(), 
    
    axis.ticks.x = element_blank(), 
    
    axis.title.y = element_blank(), 
    
    axis.text.y = element_blank(), 
    
    axis.ticks.y = element_blank(),
    
    axis.line.y = element_blank(), 
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    
    legend.position = "none" 
    
  ) +
  coord_flip()

# Combine plots
combined_plot <- xdensity + beta_plot + ydensity +  
  plot_layout(design = layout, guides = 'collect') & 
  labs(color = "Group", fill = "Group") & 
  theme(
    legend.position = "bottom", 
    plot.margin = margin(0,0,0,0)
  )

print(combined_plot)
pdf("beta_with_density.pdf")
print(combined_plot)
dev.off()
