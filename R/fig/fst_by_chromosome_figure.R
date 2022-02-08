library(tidyverse)
library(eulerr)
library(vcfR)
library(adegenet)
library(hierfstat)
library(poppr)

source("bin/functions.R")

gs_tbl <- read_rds(here::here("data", "derived", "genome_scan_results.rds"))
env_tbl <- read_rds(here::here("data", "derived", "env_assoc_tbl.rds"))

# Join the tibbles
comb_tbl <- gs_tbl %>%
  left_join(env_tbl, by = c("locus", "chrom", "pos"))

gen <- read.vcfR(here::here("data", "raw", "out.17.recode.vcf")) %>%
  vcfR2genind()

pops <- tibble(ind = indNames(gen)) %>%
  extract(ind, "pop", "(\\w+)_", remove = FALSE) %>%
  pull(pop)

pop(gen) <- pops

# Calculate basic stats for each population
pop_stats <- gen %>%
  genind2hierfstat() %>%
  basic.stats()

global_stats <- pop_stats$perloc %>%
  rownames_to_column(var = "locus") %>%
  extract(locus, c("chrom", "pos"), "HiC_scaffold_(\\d+)_arrow_ctg1_(\\d+)", remove = FALSE) %>%
  as_tibble() %>%
  mutate(chrom = as.integer(chrom), pos = as.integer(pos)) %>%
  select(locus, chrom, pos, Ho, Hs, Fst, Fis)

comb_tbl <- comb_tbl %>%
  left_join(global_stats, by = c("locus", "chrom", "pos")) %>%
  mutate(outlier = if_else(bay_outlier == TRUE | pca_outlier == TRUE | lfmm_env_assoc == TRUE | rda_env_assoc == TRUE, TRUE, FALSE))


rects <- tibble(chrom = c(2, 8, 12),
                xstart = c(4.2E7, 0, 0),
                xend = c(5.6E7, 2.4E7, 1.4E7),
                color = c("#247ba0", "#247ba0", "#247ba0"))

chrom_names <- paste("Chr", 1:19, sep = "")
names(chrom_names) <- 1:19

comb_tbl %>%
  mutate(pos_mb = pos / 1000000) %>%
  ggplot() +
  geom_rect(data = rects, aes(x = NULL, y = NULL, xmin = xstart / 1000000, xmax = xend / 1000000, ymin = 0, ymax = 1), inherit.aes = FALSE, alpha = 0.3, fill = "grey") +
  geom_point(aes(x = pos_mb, y = Fst, col = outlier, shape = any_env_outlier)) +
  facet_wrap(~chrom, ncol = 5, labeller = labeller(chrom = chrom_names)) +
  theme_minimal() +
  scale_color_manual(name = "Selection Outlier", values = c(cmh_palette[5], cmh_palette[1]),
                     breaks = c(TRUE, FALSE)) +
  scale_shape_discrete(name = "Temperature Associated", breaks = c(TRUE, FALSE)) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
  labs(x = "Position (Mb)",
       y = expression(paste("Global ", F[ST], sep = " "))) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 5))


fst_plot <- comb_tbl %>%
  mutate(pos_mb = pos / 1000000) %>%
  ggplot() +
  geom_rect(data = rects, aes(x = NULL, y = NULL, xmin = xstart / 1000000, xmax = xend / 1000000, ymin = 0, ymax = 1), inherit.aes = FALSE, alpha = 0.3, fill = "grey") +
  geom_point(aes(x = pos_mb, y = Fst, col = outlier, shape = any_env_outlier)) +
  facet_wrap(~chrom, ncol = 19, labeller = labeller(chrom = chrom_names)) +
  theme_minimal() +
  scale_color_manual(name = "Selection outlier", values = c(cmh_palette[5], cmh_palette[1]),
                     breaks = c(TRUE, FALSE)) +
  scale_shape_discrete(name = "Temperature-associated", breaks = c(TRUE, FALSE)) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
  labs(x = "Chromosomal Position",
       y = expression(paste("Global ", F[ST], sep = " "))) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_blank(),
        panel.spacing.x = unit(0.1, "lines"))


ggsave(here::here("out", "fig", "fst_chromosome_plot.png"),
       fst_plot,
       width = 32,
       height = 5,
       units = "cm")
