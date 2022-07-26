library(tidyverse)
library(patchwork)

plot_obj <- read_rds(here::here("data", "derived", "ldna_plot.rds"))

ldna_plot <- (plot_obj$chr2_network + ggtitle("Chrom 2") + plot_obj$chr2_chr) / (plot_obj$chr8_network + ggtitle("Chrom 8") + plot_obj$chr8_chr) / (plot_obj$chr12_network + ggtitle("Chrom 12") + plot_obj$chr12_chr) +
  plot_annotation(tag_levels = 'a')


ggsave(here::here("out", "fig", "ldna_plot.png"), plot = ldna_plot, width = 25, height = 32, units = "cm")