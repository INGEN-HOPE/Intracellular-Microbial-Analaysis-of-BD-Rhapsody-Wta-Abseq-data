library(tidyverse)
library(ggplot2)

proportion_data <- read.csv("~/viral_celltype_analysis/microbe_proportions.csv", check.names = FALSE)
abundance_data <- read.csv("~/viral_celltype_analysis/virus_abundance.csv", check.names = FALSE)

proportion_long <- proportion_data %>%
  pivot_longer(
    cols = -Microbe,
    names_to = c("Cell_type", "Condition"),
    names_sep = "_",
    values_to = "Proportion"
  ) %>%
  mutate(
    Microbe_CellType = paste(Microbe, Cell_type, sep = "_")
  )

abundance_long <- abundance_data %>%
  pivot_longer(
    cols = -microbe_celltype_commensal,
    names_to = "Condition",
    values_to = "Abundance"
  ) %>%
  separate_wider_delim(
    microbe_celltype_commensal,
    delim = "_",
    names = c("Microbe", "Cell_type")
  ) %>%
  mutate(
    Microbe_CellType = paste(Microbe, Cell_type, sep = "_")
  )

combined_data <- proportion_long %>%
  left_join(
    abundance_long,
    by = c("Microbe_CellType", "Condition")
  ) %>%
  filter(!is.na(Abundance))  # Remove rows with missing abundance
combined_data$Condition <- factor(combined_data$Condition, 
                                  levels = c("Healthy", "Infected", "Recovered", "Longitudinal Recovered"))

p = ggplot(combined_data, aes(x = Condition, y = Microbe_CellType)) +
  geom_point(
    aes(size = Proportion, fill = Abundance), 
    shape = 21,      # Enables border + fill
    color = "black", # Border color
    alpha = 0.8, 
    stroke = 0.8
  ) +
  scale_size_continuous(
    range = c(2, 10), 
    name = "Proportion", 
    limits = c(0, 1)
  ) +
  scale_fill_gradient2(
    low = "#516A81", mid = "#8A949D", high = "#F4AF1A", 
    midpoint = median(combined_data$Abundance),
    name = "Normalized Abundance"
  ) +
  # Reorder and Style Legends
  guides(
    fill = guide_colorbar(order = 1), # Abundance on top
    size = guide_legend(
      order = 2,                      # Proportion below
      override.aes = list(
        fill = "black",               # Makes legend circles black
        color = "black",              # Ensures border is black
        alpha = 1                     # Makes them fully opaque
      )
    )
  ) +
  labs(
    x = "Condition",
    y = "Microbe + Cell Type",
    title = "virus Abundance and Cell Proportion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12, colour = "black")
  )
ggsave("~/viral_celltype_analysis/plot_viral.pdf", p, width =10, height=14)
ggsave("~/viral_celltype_analysis/plot_viral.png", p, width =10, height=14)
ggsave("~/viral_celltype_analysis/plot_viral.svg", p, width =10, height=14)
