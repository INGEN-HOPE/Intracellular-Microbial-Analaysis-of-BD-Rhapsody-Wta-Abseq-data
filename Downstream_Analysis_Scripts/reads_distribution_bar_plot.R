# Load libraries
library(ggplot2)
library(dplyr)
library(ggsignif)
library(ggpubr)
library(rstatix)

reads = read.csv("E:/Single cell Microbes/bd/bd_reads_distribution.csv")
# Perform statistical test
stat_test <- reads %>%
  wilcox_test(Total_Reads ~ condition) %>%
  add_significance()

data <- aggregate(PER1000NONHUMAN ~ Condition, data = a, FUN = mean)
data <- rbind(data, data.frame(Condition = "Human", PER1000NONHUMAN = 1000))
data$group <- c(rep("Non Human", nrow(data) - 1), "Human")

#colors for condition
custom_colors <- c(
  "mild" = "#6f917cff", 
  "moderate" = "#ff983eff", 
  "Severe" = "#a793acff",
  "Mortality" = "#4EA6FF"
)

# Create the bar plot
p <- ggplot(reads, aes(x = condition, y = Total_Reads, fill = condition)) +
  geom_boxplot(stat = "identity", position = position_dodge(width = 0.75), width = 0.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = NULL,
    y = "Unmapped Read Count/\n1000 human reads",
    x = NULL,
    fill = NULL
  ) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1050)) +
  theme(
    legend.position = "right", 
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black")
  )

print(stat_test)

p = p + geom_signif(
  annotations = "0.00037", #add values accordingly from stat_test
  xmin = 1.7,
  xmax = 1.9, 
  y_position = 952,   
  tip_length = c(0.005, 0.005), 
  size = 0.7
)

p = p + geom_signif(
  annotations = "ns",  
  xmin = 1.7,    
  xmax = 2.1,  
  y_position = 959,  
  tip_length = c(0.005, 0.005), 
  size = 0.7
)

p = p + geom_signif(
  annotations = "0.01",  
  xmin = 1.7,    
  xmax = 2.27,  
  y_position = 993,  
  tip_length = c(0.01, 0.01), 
  size = 0.7
)

p = p + geom_signif(
  annotations = "0.00015",  
  xmin = 1.9,    
  xmax = 2.1,  
  y_position = 1021,  
  tip_length = c(0.005, 0.005), 
  size = 0.7
)

p = p + geom_signif(
  annotations = "ns",  
  xmin = 1.9,    
  xmax = 2.27,  
  y_position = 1050,  
  tip_length = c(0.005, 0.005), 
  size = 0.7
)

p = p + geom_signif(
  annotations = "0.01",  
  xmin = 2.1,    
  xmax = 2.27,  
  y_position = 900,  
  tip_length = c(0.01, 0.01), 
  size = 0.7
)


print(p)


# save plot
ggsave("E:/Single cell Microbes/FIGURES/FIG_NEW/read_count_non_human_per_1000_human_BD.svg", p, width = 10, height = 6, dpi = 600, bg="white")
ggsave("E:/Single cell Microbes/FIGURES/FIG_NEW/read_count_non_human_per_1000_human_BD.png", p, width = 10, height = 6, dpi = 600, bg="white")
