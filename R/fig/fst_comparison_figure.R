library(tidyverse)

source("bin/functions.R")

scot_nor_pw <- read_rds(here::here("data", "derived", "scot_nor_pw_fst.rds"))
scot_esp_pw <- read_rds(here::here("data", "derived", "scot_esp_pw_fst.rds"))
esp_nor_pw <- read_rds(here::here("data", "derived", "esp_nor_pw_fst.rds"))

all_reg_tbl <- bind_rows(scot_nor_pw, scot_esp_pw, esp_nor_pw) %>%
  filter(Fst > -0.05) %>%
  mutate(position_mb = pos / 1000000)

rects <- tibble(chrom = c(2, 2, 3, 8, 10, 10, 12, 13),
                xstart = c(0.8E7, 4E7, 3.2E7, 0, 1.5E7, 2.6E7, 0, 1.4E7),
                xend = c(1.2E7, 5.5E7, 3.8E7, 2.5E7, 1.8E7, 3.0E7, 1.5E7, 2E7),
                color = c("#247ba0", "#247ba0", "#69af6d", "#247ba0", "#69af6d", "#247ba0", "#247ba0","#69af6d" ))

chr_labs <- c("Chr2", "Chr3", "Chr8", "Chr10", "Chr12", "Chr13")
names(chr_labs) <- c(2, 3, 8, 10, 12, 13)

comp_labs <- c("Spain :: Norway", "Spain :: Scotland", "Scotland :: Norway")
names(comp_labs) <- c("esp_nor", "scot_esp", "scot_nor")

fst_comp_plot <- all_reg_tbl %>%
  filter(chrom %in% c(2, 3, 8, 10, 12, 13)) %>%
  ggplot(aes(x = position_mb, y = Fst, col = outlier)) +
  geom_rect(data = rects, aes(x = NULL, y = NULL, xmin = xstart / 1000000, xmax = xend / 1000000, ymin = 0, ymax = 1), inherit.aes = FALSE, alpha = 0.3, fill = rep(rects$color, 3)) +
  geom_point(size = 1) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  facet_grid(comp~chrom, labeller = labeller(chrom = chr_labs, comp = comp_labs)) +
  theme_minimal() +
  scale_color_manual(values = c(cmh_palette[1], cmh_palette[3]), labels = c("Neutral", "Outlier")) +
  labs(x = "Position (Mb)",
       y = expression(F[ST])) +
  theme(strip.text.y = element_text(angle = 0, face = "italic"),
        strip.text.x = element_text(face = "bold"),
        legend.position = "left",
        legend.title = element_blank())

ggsave(here::here("out", "fig", "fst_comparison_plot.png"),
       fst_comp_plot,
       width = 24,
       height = 12,
       units = "cm")