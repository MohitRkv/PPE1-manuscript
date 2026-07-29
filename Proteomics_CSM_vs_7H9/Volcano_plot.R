##############################################################
## Volcano Plot
##############################################################

library(readr)
library(ggplot2)
library(ggrepel)
library(grid)

##############################################################
# Read data
##############################################################

volcanof <- read_csv("CSM_vs_7H9.csv")

# Convert to numeric
volcanof$fold.change <- as.numeric(volcanof$fold.change)
volcanof$P.Value <- as.numeric(volcanof$P.Value)

# Remove missing values
volcanof <- na.omit(volcanof)

# Replace zero P-values
volcanof$P.Value[volcanof$P.Value == 0] <- 1e-300

##############################################################
# Thresholds
##############################################################

FCcutoff <- 0.2
pCutoff  <- 0.01

volcanof$minusLog10P <- -log10(volcanof$P.Value)

##############################################################
# Assign colours
##############################################################

volcanof$Group <- "NS"

volcanof$Group[
  volcanof$fold.change >= FCcutoff &
    volcanof$P.Value < pCutoff
] <- "Up"

volcanof$Group[
  volcanof$fold.change <= -FCcutoff &
    volcanof$P.Value < pCutoff
] <- "Down"

volcanof$Group <- factor(
  volcanof$Group,
  levels = c("Down", "NS", "Up")
)

##############################################################
# Genes to label
##############################################################

genes_to_label <- c("PPE1", "nrp", "fadD10")

label_df <- subset(
  volcanof,
  Gene.Symbol %in% genes_to_label
)

##############################################################
# Plot
##############################################################

p <- ggplot(
  volcanof,
  aes(
    x = fold.change,
    y = minusLog10P
  )
) +
  
  ############################################################
# Points
############################################################

geom_point(
  aes(colour = Group),
  size = 2.8,
  alpha = 0.75
) +
  
  ############################################################
# Colours and Legend
############################################################

scale_colour_manual(
  name = "Regulation",
  values = c(
    Down = "royalblue3",
    NS   = "grey75",
    Up   = "red2"
  ),
  labels = c(
    "Downregulated",
    "Not Significant",
    "Upregulated"
  )
) +
  
  ############################################################
# Threshold lines
############################################################

geom_vline(
  xintercept = c(-FCcutoff, FCcutoff),
  linetype = "dashed",
  linewidth = 0.6
) +
  
  geom_hline(
    yintercept = -log10(pCutoff),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  ############################################################
# Automatic labels
############################################################

geom_label_repel(
  data = label_df,
  aes(label = Gene.Symbol),
  fill = "white",
  colour = "black",
  label.size = 0.4,
  label.padding = unit(0.25, "lines"),
  box.padding = 0.5,
  point.padding = 0.3,
  segment.color = "black",
  segment.size = 0.5,
  min.segment.length = 0,
  max.overlaps = Inf,
  force = 2,
  seed = 123,
  size = 5,
  fontface = "plain"
) +
  
  ############################################################
# Axes
############################################################

coord_cartesian(
  xlim = c(-8, 8),
  ylim = c(0, 10)
) +
  
  labs(
    title = "Volcano Plot",
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~italic(P))
  ) +
  
  ############################################################
# Theme
############################################################

theme_classic(base_size = 18) +
  
  theme(
    plot.title = element_text(
      size = 24,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      size = 20,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 16,
      colour = "black"
    ),
    
    legend.position = "right",
    
    legend.title = element_text(size = 12, face = "bold"),
    
    legend.text = element_text(size = 12),
    
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.5
    ),
    
    legend.key = element_blank(),
    
    axis.line = element_line(
      linewidth = 0.8
    )
  )

##############################################################
# Display
##############################################################

print(p)

##############################################################
# Save
##############################################################

ggsave(
  "Volcano_plot7.pdf",
  plot = p,
  width = 11,
  height = 11,
  dpi = 600
)

ggsave(
  "Volcano_plot7.png",
  plot = p,
  width = 11,
  height = 11,
  dpi = 600
)
