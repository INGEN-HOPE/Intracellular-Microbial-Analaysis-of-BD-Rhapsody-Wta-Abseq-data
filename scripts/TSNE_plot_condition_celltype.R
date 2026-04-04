#load libraries
library(Seurat)
library(ggplot2)

bd = readRDS("bd.rds")

#Visualize Tsne coloured by condition
pdf("tsne_condition.pdf", width=12, height=8)
DimPlot(bd, reduction = "tsne", group.by = "condition", raster=FALSE) +
  scale_color_manual(values = c("Infected" = "#FF7F00",   
                                "Recovered" = "#9c87a1ff",  
                                "Healthy" = "#005f3c",    
                                "Longitudinal Recovered" = "#1E90FF")) +  
  theme(
    plot.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14, face = "bold")
  )

dev.off()


# Visualize Tsne coloured by celltype

pdf("tsne_plot_celltype.pdf", width = 10, height=8)

DimPlot(bd, reduction = "tsne", group.by = "celltype", raster=FALSE, label=FALSE) +
  scale_color_manual(values = c("Memory T Cell" = "#AFD3E2",
                                "Plasma Cell" = "#5DA7DF",
                                "Naive T Cell" = "#C7DCA7",
                                "B Cell" = "#7AA874",
                                "Activated T Cell" = "#FFC2B5",
                                "Cytotoxic T Cell" = "#D9534F",
                                "Monocyte" = "#A59D84",
                                "NK Cell" = "#F58634",
                                "Regulatory T Cell" = "#D9D7F1",
                                "Erythrocyte" = "#6F69AC",
                                "Platelets" = "#f2b72b",
                                "Dendritic Cell" = "#EAAC7F",
                                "NK T Cell" = "#AC7D88")) +
  theme(
    plot.title = element_blank(),
    legend.text = element_text(size = 14, colour="#808080"),
    axis.text = element_text(size = 14, face = "bold")
  )

dev.off()


