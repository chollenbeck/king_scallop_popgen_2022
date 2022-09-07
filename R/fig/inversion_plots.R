library(tidyverse)
library(patchwork)

##### Chromosome 2 Inversion Plots #####

chr2_plots <- read_rds(file = "data/derived/chr2_inversion_plots.rds")

chr2_plot <- chr2_plots$gt / ((chr2_plots$pca / chr2_plots$barplot) | chr2_plots$map) + plot_layout(heights = c(5, 20)) +
  plot_annotation(tag_levels = 'A')

ggsave(filename = "out/fig/chrom2_inversion.png", plot = chr2_plot, width = 28, height = 18,
       units = "cm")

ggsave(filename = "out/fig/chrom2_inversion.pdf", plot = chr2_plot, width = 28, height = 18,
       units = "cm")

##### Chromosome 12 Inversion Plots #####

chr12_plots <- read_rds(file = "data/derived/chr12_inversion_plots.rds")

chr12_plot <- chr12_plots$gt / ((chr12_plots$pca / chr12_plots$barplot) | chr12_plots$map) + plot_layout(heights = c(5, 20)) +
  plot_annotation(tag_levels = 'A')

ggsave(filename = "out/fig/chrom12_inversion.png", plot = chr12_plot, width = 28, height = 18,
       units = "cm")

ggsave(filename = "out/fig/chrom12_inversion.pdf", plot = chr12_plot, width = 28, height = 18,
       units = "cm")

##### Chromosome 8 Inversion Plots #####

chr8_pca <- read_rds(file = "data/derived/chr8_inversion_plots.rds")

chr8_pca

ggsave(filename = "out/fig/chrom8_inversion_pca.png", plot = chr8_pca, width = 12, height = 8,
       units = "cm")

ggsave(filename = "out/fig/chrom8_inversion_pca.pdf", plot = chr8_pca, width = 12, height = 8,
       units = "cm")