# 1. Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# 2. Data Preparation
df_raw <- data.frame(
  Phylum = c('Pseudomonadota', 'Actinomycetota', 'Bacillota', 'Mycoplasmatota', 
             'Peploviricota', 'Artverviricota', 'Unclassified', 'Campylobacterota', 
             'Bacteroidota', 'Duplornaviricota', 'Nucleocytoviricota', 'Pisuviricota'),
  DHF = c(32.36, 20.11, 15.41, 9.62, 6.04, 3.53, 2.85, 2.17, 1.18, 1.14, 0, 0),
  PBMC = c(46.56, 8.78, 12.73, 9.33, 4.69, 3.89, 4.47, 1.49, 1.15, 1.01, 1.34, 1.03)
)

# 3. Reshape and Sort
df_long <- df_raw %>%
  pivot_longer(cols = c(DHF, PBMC),
               names_to = "Condition",
               values_to = "Abundance") %>%
  mutate(
    Phylum = reorder(Phylum, -Abundance),
    Condition = recode(Condition,
                       DHF = "Dengue Infected",
                       PBMC = "Healthy Control")
  )

# 4. Create the Plot
p = ggplot(df_long, aes(x = Phylum, y = Abundance, fill = Condition)) +
  # Decreased bar width (width = 0.6) and dodge (width = 0.7)
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.7), 
           width = 0.6) +
  
  # Add labels ONLY if Abundance is greater than 0
  geom_text(aes(label = ifelse(Abundance > 0, paste0(Abundance, "%"), "")), 
            position = position_dodge(width = 0.8), 
            vjust = -0.6, 
            size = 3.2, 
            fontface = "bold", 
            color = "#2C3E50") +
  
  # Colors
  scale_fill_manual(values = c("Dengue Infected" = "#FF6B6B", "Healthy Control" = "#4ECDC4")) +
  
  # Refined Styling
  theme_minimal(base_family = "sans") +
  labs(
    title = "Phylum Percent Relative Abundance",
    x = NULL,
    y = "Relative Abundance (%)",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    axis.text.x = element_text(angle = 40, hjust = 1, face = "italic", size = 16, color = "black"),
    axis.text.y = element_text( size = 16, color = "black"),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    legend.text = element_text( size = 16, color = "black"),
    legend.title = element_text( size = 16, color = "black")
  ) +
  # Added extra space at the top of the Y-axis for the labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))


svglite::svglite(
  paste0("~/public_dengue_public/phylum_barplot.svg"),
  width = 14, height = 8
)
print(p)
dev.off()






# 4. Create the Plot
p = ggplot(df_long, aes(x = Phylum, y = Abundance, fill = Condition)) +
  # Decreased bar width (width = 0.6) and dodge (width = 0.7)
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.8), 
           width = 0.6) +
  
  # Add labels ONLY if Abundance > 0 (Fixed dodge width to 0.7 to match bars)
  geom_text(aes(label = ifelse(Abundance > 0, paste0(Abundance, "%"), "")), 
            position = position_dodge(width = 0.9), 
            vjust = -0.6, 
            size = 4, # Slightly larger for the 14x8 SVG
            fontface = "bold", 
            color = "#2C3E50") +
  
  # Colors (Ensure these match the levels in your 'Condition' column)
  scale_fill_manual(values = c("Dengue Infected" = "#FF6B6B", "Healthy Control" = "#4ECDC4")) +
  
  # Refined Styling
  theme_minimal(base_family = "sans") +
  labs(
    title = "Phylum Percent Relative Abundance",
    x = NULL,
    y = "Relative Abundance (%)",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 18, color = "black"),
    axis.text.x = element_text(angle = 40, hjust = 1, face = "italic", size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 15)), # Space between title and axis
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black"),
    
    # --- ADJUST LEFT MARGIN HERE ---
    # margin(top, right, bottom, left, unit)
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") 
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Exporting
svglite::svglite(
  paste0("~/public_dengue_public/phylum_barplot.svg"),
  width = 15, height = 8
)
print(p)
dev.off()
