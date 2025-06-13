#################################################
# CHAPTER 2. Investigating pangenome structure #
################################################

# Set the dir

main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

# Load/install required packages

if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(ggplot2, dplyr, tidyr, ggnewscale, ggrepel, patchwork, ggtext)

#######################
# Part 1: Donut chart #
#######################

gene_counts <- read.csv("pangenome/StAl_pangenome.tsv", sep = '\t')

# Compute percentages
gene_counts$fraction = gene_counts$Count / sum(gene_counts$Count)

# Compute the cumulative percentages (top of each rectangle)
gene_counts$ymax = cumsum(gene_counts$fraction)

# Compute the bottom of each rectangle
gene_counts$ymin = c(0, head(gene_counts$ymax, n = -1))

# Compute label position
gene_counts$labelPosition <- (gene_counts$ymax + gene_counts$ymin) / 2

# Compute a good label
gene_counts$label <- paste0(gene_counts$Category,
                            "\n value: ",
                            gene_counts$Count,
                            "\n",
                            round((gene_counts$fraction * 100), 2),
                            "%")

#3B7C70FF, #CE9642FF, #898E9FFF, #3B3A3EFF

colors <- c(
  "Core" = "#3B3A3EFF",
  "Soft-core" = "#898E9FFF",
  "Shell" = "#3B7C70FF",
  "Cloud" = "#CE9642FF"
)

# 4. Make the donut plot
pangenome_donut <- ggplot(gene_counts,
                          aes(
                            ymax = ymax,
                            ymin = ymin,
                            xmax = 4,
                            xmin = 3,
                            fill = Category
                          )) +
  geom_rect() +
  geom_text(x = 1.5,
            aes(y = labelPosition, label = label, color = Category),
            size = 5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  coord_polar(theta = "y") +
  xlim(c(-1, 4)) +
  labs(caption = "*Streptomyces albidoflavus* pangenome composition") +
  theme_void() +
  theme(legend.position = "none",
        plot.caption = element_markdown(size = 14, hjust = 0.5))


ggsave(
  "imgs/donut_chart.png",
  pangenome_donut,
  width = 8,
  height = 8,
  dpi = 600
)

#####################
### Part 2: SM254 ###
#####################


# Read the data
sm254_df <- read.delim("pangenome/sm254_pangenome.tsv",
                       sep = "\t",
                       header = TRUE)

# Set 'Category' as a factor with correct order and new labels
sm254_df$Category <- factor(
  sm254_df$Category,
  levels = c("Cloud", "Shell", "Soft-core", "Core"),
  # Desired order
  labels = c("Cloud genes", "Shell genes", "Soft-core genes", "Core genes")  # New names
)

# Plot with ggplot2
sm254_pangenome <- ggplot(sm254_df, aes(x = Count, y = Category, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), hjust = -0.2, size = 4) +
  theme_minimal() +
  labs(x = "", y = "") +
  scale_fill_manual(
    values = c(
      "Core genes" = "#3B3A3EFF",
      # Blue
      "Soft-core genes" = "#898E9FFF",
      # Green
      "Shell genes" = "#3B7C70FF",
      # Orange
      "Cloud genes" = "#CE9642FF"       # Red
    )
  ) +
  xlim(0, max(sm254_df$Count) * 1.05) +
  labs(caption = "Total genes (*Streptomyces albidoflavus* SM254)") +
  theme(
    legend.position = "none",
    plot.caption = element_markdown(size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  "imgs/sm254_pangenome.png",
  sm254_pangenome,
  width = 10,
  height = 6,
  dpi = 600
)

######################
###### COMBINED ######
######################

everything <- (pangenome_donut + sm254_pangenome) + plot_annotation(tag_levels = list(c("A", "B")))
ggsave(
  "imgs/pangenomes.png",
  plot = everything,
  width = 18,
  height = 8,
  dpi = 600
)
