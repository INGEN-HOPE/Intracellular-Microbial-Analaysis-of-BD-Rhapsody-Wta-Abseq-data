
library(ggpubr)
library(dplyr)
library(ggplot2)
library(rstatix)
library(tidyverse)
library(tibble)


df <- read.csv("~/Downloads/50_NORM_AMONG_4_GROUP.csv", row.names = 1)
anno <- read.csv("~/anno_ra_bd.csv") %>% column_to_rownames("Sample")

df= t(df)
df = df[-1,]
df = merge(df, anno, by = "row.names") %>% column_to_rownames("Row.names")
colnames(df) = gsub("-", "_", colnames(df))
colnames(df) = gsub(" ", "_", colnames(df))
# Get the names of numeric columns
cols_to_convert <- setdiff(colnames(df), "Condition")
# Convert those columns to numeric
df[cols_to_convert] <- lapply(df[cols_to_convert], function(x) as.numeric(as.character(x)))

#  colors
pastel_colors <- c("Healthy" = "#6f917c", "Infected" = "#ff983e", "Recovered" = "#a793ac", "Long Recovered" = "#4EA6FF")

# Get all microbe columns
microbe_cols <- setdiff(colnames(df), "Condition")

# Set the desired order of Condition
df$Condition <- factor(df$Condition, levels = c("Healthy", "Infected", "Recovered", "Long Recovered"))

# Loop through each microbe and save plot
for (microbe in microbe_cols) {
  # Perform Wilcoxon test
  stat_test <- df %>%
    wilcox_test(as.formula(paste(microbe, "~ Condition"))) %>%
    add_significance()
  
  # Adjust y-position for p-value labels
  y_max <- max(df[[microbe]], na.rm = TRUE)
  y_positions <- seq(1.1 * y_max, length.out = nrow(stat_test), by = 0.03 * y_max)
  
  # Create the plot
  p <- ggplot(data = df, aes_string(x = "Condition", y = microbe)) +
    geom_violin(aes(color = Condition, fill = Condition), trim = FALSE, alpha = 0.25, linewidth = 1) +
    geom_boxplot(aes(color = Condition), fill = NA, width = 0.2, alpha = 0.5, linewidth = 1) +
    geom_jitter(aes(color = Condition, fill = Condition), shape = 21, width = 0.3, size = 2, alpha = 0.7) +
    scale_fill_manual(values = pastel_colors) +
    scale_color_manual(values = pastel_colors) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14, colour = "black"),
      axis.text.y = element_text(size = 14, colour = "black"),
      axis.title = element_text(size = 14, colour = "black"),
      axis.text = element_text(size = 14)
    ) +
    stat_pvalue_manual(
      stat_test,
      label = "p",
      y.position = y_positions,
      tip.length = 0.01,
      size = 4
    ) +
    ylab(microbe)
  
  # Clean filename (remove special characters)
  safe_name <- gsub("[^A-Za-z0-9_]", "_", microbe)
  
  # Save the plot
  ggsave(
    filename = paste0("E:/Single cell Microbes/FIGURES/FIG_NEW/box/", microbe, "_sc_box_MICROBES.svg"),
    plot = p,
    width = 8,
    height = 8,
    dpi = 600,
    bg = "white"
  )
}
