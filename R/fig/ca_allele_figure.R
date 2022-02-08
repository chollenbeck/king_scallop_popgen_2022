library(tidyverse)
library(adegenet)
library(vcfR)
library(patchwork)

source("bin/functions.R")

# Data import

temp_tbl <- read_tsv(here::here("data", "derived", "mean_temps.tsv")) %>%
  mutate(temp_scaled = temp - mean(temp))


# Correspondence analysis

out <- read_rds(here::here("data", "derived", "outlier_genind.rds"))
# Fix the pop names

popNames(out) <- popNames(out) %>%
  str_replace("_\\d+", "")

# Read in the PCA colors
pca_col <- read_tsv(file = here::here("data", "raw", "colors.tsv"))



# Run a CA on the genind object
obj <- genind2genpop(out)
ca1 <- dudi.coa(tab(obj), scannf=FALSE, nf=10)
out_ca <- ca1$li[,1:2]
rownames(out_ca) <- popNames(out)

out_ca <- out_ca %>%
  rownames_to_column(var = "site") %>%
  mutate(site = factor(site, levels = c("ESP", "SW", "SE", "NW", "NE", "SLD", "SNO", "NNO"))) %>%
  select(site, Axis1, Axis2)

temps <- read_tsv(here::here("data", "raw", "mean_temps.tsv")) %>%
  left_join(out_ca, by = "site") %>%
  mutate(group = if_else(site == "SNO" | site == "NNO", "north", "south")) %>%
  mutate(site = factor(site, levels = c("ESP", "SW", "SE", "NW", "NE", "SLD", "SNO", "NNO")))


ca_plot <- temps %>%
  #filter(! site %in% c("SNO", "NNO")) %>%
  #mutate(site = as.character(site)) %>%
  ggplot(aes(x = Axis2, y = Axis1)) +
  geom_point(aes(color = site), size = 4) +
  scale_color_manual(values = pca_col$col[1:8]) +
  theme_minimal() +
  labs(x = "Outlier CA2",
       y = "Outlier CA1") +
  theme(legend.title = element_blank(),
        legend.position = "none")

temp_plot <- temps %>%
  #filter(! site %in% c("SNO", "NNO")) %>%
  #mutate(site = as.character(site)) %>%
  ggplot(aes(x = temp, y = Axis1)) +
  geom_smooth(method = "lm", col = "darkgrey", lty = 1) +
  geom_point(aes(color = site), size = 4) +
  scale_color_manual(values = pca_col$col[1:8]) +
  theme_minimal() +
  labs(x = "Mean Sea Surface Temperature (C)",
       y = "Outlier CA1") +
  theme(legend.title = element_blank(),
        legend.position = "none")


temp_plot2 <- temps %>%
  #filter(! site %in% c("SNO", "NNO")) %>%
  #mutate(site = as.character(site)) %>%
  ggplot(aes(x = -Axis1, y = temp)) +
  geom_smooth(method = "lm", col = "darkgrey", lty = 1) +
  geom_point(aes(color = site), size = 4) +
  scale_color_manual(values = pca_col$col[1:8]) +
  theme_minimal() +
  labs(y = "Mean Sea Surface Temperature (C)",
       x = "Outlier CA1") +
  theme(legend.title = element_blank(),
        legend.position = "none")

temps %>%
  #filter(! site %in% c("SNO", "NNO")) %>%
  #mutate(site = as.character(site)) %>%
  ggplot(aes(x = temp, y = Axis1)) +
  geom_smooth(method = "lm", col = "darkgrey", lty = 1) +
  geom_point(aes(color = site), size = 4) +
  scale_color_manual(values = pca_col$col[1:8]) +
  theme_minimal()

ca_loadings <- ca1$co %>%
  as_tibble(rownames = "allele") %>%
  extract(allele, c("locus", "chrom", "pos", "allele"), "(HiC_scaffold_(\\d+)_arrow_ctg1_(\\d+)).(\\d)") %>%
  mutate(pos = as.integer(pos)) %>%
  select(locus, chrom, pos, allele, Comp1, Comp2) %>%
  mutate(Comp1 = Comp1^2, Comp2 = Comp2^2) %>%
  arrange(desc(Comp1))

top_loadings <- ca_loadings %>%
  top_n(16, Comp1)

facet_order <- top_loadings %>%
  arrange(chrom, pos) %>%
  mutate(facet_title = glue::glue("Chr{chrom}:{pos}")) %>%
  pull(facet_title)

out %>%
  pantomime::get_allele_freqs()

af_plot <- out %>%
  pantomime::get_allele_freqs() %>%
  filter(allele == 0, locus %in% top_loadings$locus) %>%
  extract(locus, c("chrom", "pos"), "HiC_scaffold_(\\d+)_arrow_ctg1_(\\d+)", remove = FALSE) %>%
  mutate(pos = as.integer(pos)) %>%
  mutate(facet_title = glue::glue("Chr{chrom}:{pos}")) %>%
  mutate(facet_title = factor(facet_title, levels = facet_order)) %>%
  mutate(pop = factor(pop, levels = c("ESP", "SW", "SE", "NW", "NE", "SLD", "SNO", "NNO"))) %>%
  ggplot(aes(x = pop, y = freq, fill = pop)) +
  geom_col() +
  scale_fill_manual(values = pca_col$col) +
  facet_wrap(~facet_title) +
  theme_minimal() +
  labs(x = "Locality", y = "Allele frequency") +
  theme(axis.text.x = element_text(angle = 90),
        legend.title = element_blank(),
        strip.text = element_text(face = "bold"))

ca_freq_plot <- (ca_plot + temp_plot) / af_plot + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')

ggsave(filename = "out/fig/ca_frequency_figure.png", plot = ca_freq_plot, width = 22, height = 20,
       units = "cm")

