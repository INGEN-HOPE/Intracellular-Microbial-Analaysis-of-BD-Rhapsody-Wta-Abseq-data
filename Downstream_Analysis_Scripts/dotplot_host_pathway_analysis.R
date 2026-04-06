library(ggplot2)
library(dplyr)
library(stringr)

files <- list.files("~/Downloads/pesudobulk_host_dotplot", pattern = "^dotplot_.*\\.csv$", full.names = TRUE)
conditions <- gsub("dotplot_|\\.csv", "", basename(files))

for (i in seq_along(files)) {
  cond <- conditions[i]
  
  df <- read.csv(files[i]) %>%
    na.omit() %>%
    mutate(
      GeneRatio = sapply(as.character(GeneRatio), function(x) eval(parse(text = x))),
      neg_log_padj = ifelse(is.infinite(-log10(p.adjust)), 
                            max(-log10(p.adjust[is.finite(-log10(p.adjust))])) + 1, 
                            -log10(p.adjust)),
      Description  = factor(str_wrap(Description, 45), levels = rev(unique(str_wrap(Description, 45))))
    )
  
  p <- ggplot(df, aes(x = GeneRatio, y = Description, size = Count, color = neg_log_padj)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(name = "Gene count", range = c(4, 12), breaks = c(5, 10, 16, 24, 30)) +
    scale_color_gradient(name = expression(-log[10](p.adjust)), low = "#96CDCD", high = "#8B668B") +
    facet_grid(regulation ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Gene Ratio", y = NULL, title = paste("Pathway Enrichment -", cond)) +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 14),
      axis.text.y        = element_text(size = 22, colour = "black"),
      axis.text.x = element_text(size = 22, colour = "black"),
      strip.text         = element_text(face = "bold", size = 22),
      strip.background   = element_rect(fill = "#f0f0f0", color = "gray70"),
      panel.grid.major.y = element_line(linetype = "dashed", color = "gray85"),
      legend.position    = "bottom"
    )
  
  for (ext in c("pdf", "png")) {
    ggsave(
      paste0("~/Downloads/pesudobulk_host_dotplot/pathway_dotplot_", cond, ".", ext),
      p, width = 13, height = 12, dpi = 600
    )
  }
  
  # SVG via svglite
  svglite::svglite(
    paste0("~/Downloads/pesudobulk_host_dotplot/pathway_dotplot_", cond, ".svg"),
    width = 13, height = 12
  )
  print(p)
  dev.off()
}

