# Set the dir

main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

# Load/install required packages

#install.packages('BiocManager')

if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  ggtree,
  ggimage,
  ggtext,
  phangorn,
  ggplot2,
  treeio,
  ggnewscale,
  viridis,
  phytools,
  patchwork
)

####################
# METADATA PARSING #
####################

# Load `metadata`

metadata <- read.table('metadata/welldone_metadata.tsv', sep = '\t', header = T)

metadata$strain[metadata$AN == "NZ_CP014485"] <- "NZ_CP014485 S. albidoflavus SM254"
metadata$strain[metadata$AN == "NC_020990"] <- "NC_020990 S. albidoflavus"
metadata$strain[metadata$AN == "NZ_CP133227"] <- "NZ_CP133227 S. albidoflavus RKJM-0023"
metadata$strain[metadata$AN == "NZ_CP113228"] <- "NZ_CP113228 S. albidoflavus J1074_D14"
metadata$strain[metadata$AN == "NZ_CP064783"] <- "NZ_CP064783 S. albidoflavus W68"
metadata$strain[metadata$AN == "NZ_CP170353"] <- "NZ_CP170353 S. albidoflavus S20"
metadata$strain[metadata$AN == "NZ_CP128384"] <- "NZ_CP128384 S. albidoflavus MGMM6"

metadata$strain[metadata$AN == "NZ_CP152102"] <- "NZ_CP152102 S. albidoflavus SC-3"
metadata$strain[metadata$AN == "NZ_CP109238"] <- "NZ_CP109238 S. albidoflavus NBC_01664"
metadata$strain[metadata$AN == "NZ_CP059254"] <- "NZ_CP059254 S. albidoflavus J1074/R2"
metadata$strain[metadata$AN == "NZ_CP108451"] <- "NZ_CP108451 S. albidoflavus NBC_01270"
metadata$strain[metadata$AN == "NZ_CP109173"] <- "NZ_CP109173 S. albidoflavus NBC_01722"
metadata$strain[metadata$AN == "NZ_CP108450"] <- "NZ_CP108450 S. albidoflavus NBC_01271"
metadata$strain[metadata$AN == "NZ_CP108379"] <- "NZ_CP108379 S. albidoflavus NBC_01331"

metadata$strain[metadata$AN == "NZ_CP140156"] <- "NZ_CP140156 S. albidoflavus ATCC 23899"
metadata$strain[metadata$AN == "NZ_CP108647"] <- "NZ_CP108647 S. albidoflavus NBC_01110"
metadata$strain[metadata$AN == "NZ_CP109235"] <- "NZ_CP109235 S. albidoflavus NBC_01665"
metadata$strain[metadata$AN == "NZ_CP109270"] <- "NZ_CP109270 S. albidoflavus NBC_01640"
metadata$strain[metadata$AN == "NZ_CP109243"] <- "NZ_CP109243 S. albidoflavus NBC_01661"
metadata$strain[metadata$AN == "NZ_CP109281"] <- "NZ_CP109281 S. albidoflavus NBC_01627"
metadata$strain[metadata$AN == "NZ_CP109226"] <- "NZ_CP109226 S. albidoflavus NBC_01671"

metadata$strain[metadata$AN == "NZ_CP109224"] <- "NZ_CP109224 S. albidoflavus NBC_01673"
metadata$strain[metadata$AN == "NZ_CP109594"] <- "NZ_CP109594 S. albidoflavus NBC_01675"
metadata$strain[metadata$AN == "NZ_CP109294"] <- "NZ_CP109294 S. albidoflavus NBC_01621"
metadata$strain[metadata$AN == "NZ_CP109570"] <- "NZ_CP109570 S. albidoflavus NBC_01374"
metadata$strain[metadata$AN == "NZ_CP108610"] <- "NZ_CP108610 S. albidoflavus NBC_01170"
metadata$strain[metadata$AN == "NZ_CP109088"] <- "NZ_CP109088 S. albidoflavus NBC_01790"
metadata$strain[metadata$AN == "NZ_CP109085"] <- "NZ_CP109085 S. albidoflavus NBC_01791"

metadata$strain[metadata$AN == "NZ_CP109203"] <- "NZ_CP109203 S. albidoflavus NBC_01692"
metadata$strain[metadata$AN == "NZ_CP109142"] <- "NZ_CP109142 S. albidoflavus NBC_01747"
metadata$strain[metadata$AN == "NZ_CP079112"] <- "NZ_CP079112 S. albidoflavus LGO-A16"
metadata$strain[metadata$AN == "NZ_CP079113"] <- "NZ_CP079113 S. albidoflavus LGO-A23"
metadata$strain[metadata$AN == "NZ_OX371412"] <- "NZ_OX371412 S. albidoflavus CCOS 2040 isolate Stup19_F108"
metadata$strain[metadata$AN == "NZ_CP109299"] <- "NZ_CP109299 S. albidoflavus NBC_01616"

# `Year` column

meta.year <- as.data.frame(metadata[, 'Year'])
colnames(meta.year) <- 'Year'
rownames(meta.year) <- metadata$AN
meta.year$Year[meta.year$Year == "ND"] <- NA

country_code_map <- c(
  "USA" = "us",
  "Germany" = "de",
  "Turkey" = "tr",
  "China" = "cn",
  "Russia" = "ru",
  "Netherlands" = "nl",
  "Brazil" = "br",
  "Switzerland" = "ch",
  "Denmark" = "dk"
)

metadata$flag_code <- country_code_map[metadata$Country]
metadata$flag_code[is.na(metadata$flag_code)] <- NA

if (!dir.exists("flags"))
  dir.create("flags")

  flag_urls <- paste0(
    "https://raw.githubusercontent.com/HatScripts/circle-flags/gh-pages/flags/",
    unique(metadata$flag_code[!is.na(metadata$flag_code)]),
    ".svg"
  )
  
  for (i in seq_along(flag_urls)) {
    country_code <- sub(".*/flags/(.*)\\.svg", "\\1", flag_urls[i])
    
    destfile <- paste0("flags/", country_code, ".svg")
    
    download.file(flag_urls[i], destfile, mode = "wb")
  }

metadata$flag_path <- paste0("flags/", metadata$flag_code, ".svg")
metadata$flag_path[metadata$flag_path == "flags/NA.svg"] <- NA

##############
# StAl TREE #
##############

# Read the tree file

StAl_tree <- read.tree("pangenome/tree/StAl_ufb.treefile")


# Draft tree

StAl_tree_fig <- ggtree(StAl_tree) %<+% metadata +
  xlim(0, 0.14) +
  #scale_color_identity() +
  geom_treescale(x = 0, y = -0.25, width = 0.01)

StAl_tree_fig <- StAl_tree_fig + geom_tiplab(
  aes(image = flag_path),
  geom = "image",
  offset = 0.0375,
  align = TRUE,
  size = 0.02,
  linesize = 0
)

StAl_boot <- StAl_tree_fig$data
StAl_boot <- StAl_boot[!StAl_boot$isTip, ]
StAl_boot$label <- as.numeric(StAl_boot$label)
StAl_boot$bootstrap <- '0'
StAl_boot$bootstrap[StAl_boot$label >= 70] <- '1'
StAl_boot$bootstrap[is.na(StAl_boot$label)] <- '1'

StAl_tree_fig <- StAl_tree_fig + new_scale_color() +
  geom_tree(data = StAl_boot, aes(color = bootstrap == '1')) +
  scale_color_manual(name = 'Bootstrap',
                     values = setNames(c("black", "grey"), c(T, F)),
                     guide = "none") +
  geom_hilight(
    mapping = aes(subset = node %in% c(40), fill = S),
    fill = "steelblue",
    alpha = .6,
    extend = 0.1
  ) +
  geom_tiplab(
    aes(
      label = strain,
      fontface = ifelse(grepl("^NZ_CP014485", strain), "bold", "plain"),
    ),
    align = TRUE,
    geom = "label",
    fill = "white",
    label.size = 0
  )

# Add `Year` heatmap

StAl_tree_fig <- gheatmap(
  StAl_tree_fig,
  meta.year,
  width = 0.025,
  offset = 0.04,
  color = "black",
  font.size = 4,
  colnames_offset_y = -0.05
) +
  scale_fill_viridis(
    option = "D",
    name = "Year",
    discrete = TRUE,
    na.translate = TRUE
  )

StAl_tree_fig <- StAl_tree_fig +
  annotate(
    "text",
    x = max(StAl_tree_fig$data$x) + 0.03825,
    y = -0.05,
    label = "Country",
    size = 4
  )

ggsave(
  'imgs/StAl.png',
  StAl_tree_fig,
  width = 25,
  height = 15,
  dpi = 600
)
