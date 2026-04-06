library(dplyr)
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

file_map <- c(
  "Healthy"                = "~/10reads/Microbe_neg_pos_10_reads_zero_Healthy.xlsx",
  "Infected"               = "~/10reads/Microbe_neg_pos_10_reads_zero_Infected.xlsx",
  "Longitudinal Recovered" = "~/10reads/Microbe_neg_pos_10_reads_zero_Longitudinal Recovered.xlsx",
  "Recovered"              = "~/10reads/Microbe_neg_pos_10_reads_zero_Recovered.xlsx"
)

P_ADJ_CUTOFF <- 0.05

# create list of lists of condition and their celltype within that data[[condition]][[cell_type]]

data <- imap(file_map, function(filepath, condition) {
  filepath <- path.expand(filepath)
  
  if (!file.exists(filepath)) {
    message("[SKIP] File not found: ", filepath)
    return(NULL)
  }
  
  message("[READ] ", condition, "  →  ", basename(filepath))
  sheet_names <- excel_sheets(filepath)
  
  sheets <- map(set_names(sheet_names), function(sheet) {
    read_excel(filepath, sheet = sheet)
  })
  
  message("       ", length(sheets), " cell types: ", 
          paste(names(sheets), collapse = ", "))
  
  sheets
})

# Remove any NULLs (missing files)
data <- compact(data)

# Filter p_val_adj < cutoff on every sheet  data_filtered[[condition]][[cell_type]]

data_filtered <- map(data, function(sheets) {
  map(sheets, function(df) {
    if (!"p_val_adj" %in% colnames(df)) {
      warning("'p_val_adj' column missing — skipping sheet")
      return(NULL)
    }
    df %>% filter(p_val_adj < P_ADJ_CUTOFF)
  }) %>% compact()
})

#Combine into one flat data frame where 'cell_type' and 'condition' columns 
combined_df <- imap(data_filtered, function(sheets, condition) {
  imap(sheets, function(df, cell_type) {
    if (nrow(df) == 0) return(NULL)
    df %>%
      mutate(
        cell_type = cell_type,
        condition  = condition,
        .before    = everything()
      )
  }) %>%
    compact() %>%
    bind_rows()
}) %>%
  bind_rows()


filtered_df <- combined_df %>%
  filter(
    (condition == "Healthy" & cell_type %in% c("Plasma Cell", "Regulatory T Cell")) |
      (condition == "Infected" & cell_type %in% c("Memory T Cell", "Plasma Cell", "Cytotoxic T Cell", "Regulatory T Cell", "Platelets", "Dendritic Cell")) |
      (condition == "Recovered" & cell_type %in% c("Plasma Cell", "B Cell", "Cytotoxic T Cell", "Monocyte", "Regulatory T Cell", "Platelets", "Dendritic Cell")) |
      (condition == "Longitudinal Recovered" & cell_type %in% c("Memory T Cell", "B Cell", "Monocyte"))
  )

plot_df <- plot_df %>%
  mutate(condition = factor(condition,
                            levels = c("Healthy", "Infected", "Recovered", "Longitudinal Recovered"))) %>%
  arrange(condition, cell_type) %>%   # 🔑 KEY STEP
  mutate(celltype_cond = factor(paste(condition, cell_type, sep = " | "),
                                levels = unique(paste(condition, cell_type, sep = " | "))))

x = ggplot(plot_df, aes(x = avg_log2FC, y = celltype_cond, color = regulation)) +
  geom_quasirandom(aes(group = interaction(celltype_cond, regulation)),
                   dodge.width = 0.6, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Down" = "#96CDCD", "Up" = "#8B668B")) +
  labs(x = "avg_log2FC", y = "Cell Type (by Condition)") +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))

bar_df <- plot_df %>%
  count(condition, cell_type, regulation) %>%
  mutate(
    condition = factor(condition,
                       levels = c("Healthy", "Infected", "Recovered", "Longitudinal Recovered"))
  ) %>%
  arrange(condition, cell_type) %>%
  mutate(celltype_cond = factor(paste(condition, cell_type, sep = " | "),
                                levels = unique(paste(condition, cell_type, sep = " | "))))

y = ggplot(bar_df, aes(x = celltype_cond, y = n, fill = regulation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Down" = "#96CDCD", "Up" = "#8B668B")) +
  theme_classic() +
  theme(axis.text.x = element_blank(),legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "", y = "Count")


z = (y / x) + plot_layout(heights = c(1, 3)) &
  theme(plot.margin = margin(0, 5, 0, 5))


ggsave("combined_plot.svg", z, height = 10, width = 7)
ggsave("combined_plot.png", z, height = 10, width = 7)
ggsave("combined_plot.pdf", z, height = 10, width = 7)
