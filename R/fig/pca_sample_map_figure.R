library(tidyverse)
library(adegenet)
library(hierfstat)
library(ggmap)
library(PBSmapping)
library(RColorBrewer)
library(hrbrthemes)
library(patchwork)
library(gridExtra)

# Read in the neutral and outlier haplotype data
out <- read_rds(here::here("data", "derived", "outlier_genind.rds"))
neut <- read_rds(here::here("data", "derived", "neutral_genind.rds"))

# Fix the pop names
popNames(neut) <- popNames(neut) %>%
  str_replace("_\\d+", "")

popNames(out) <- popNames(out) %>%
  str_replace("_\\d+", "")

##### PCA #####


# Read in the PCA colors
pca_col <- read_tsv(file = here::here("data", "raw", "colors.tsv"))

# Plot a PCA for the neutral loci
x_neut <- scaleGen(neut, NA.method = "mean")
neut_pca <- dudi.pca(x_neut, center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)
eig_percent_neut <- round((neut_pca$eig / (sum(neut_pca$eig))) * 100, 2)

pca_dat_neut <- as.data.frame(neut_pca$li) %>%
  rownames_to_column(var = "ind") %>%
  cbind(pop = pop(neut)) %>%
  mutate(pop = factor(pop, levels = c("ESP", "SW", "SE", "NW", "NE", "SLD", "SNO", "NNO"))) %>%
  mutate(type = "neut", loc.type = "haps", location = "all")


neut_pca <- ggplot(data = pca_dat_neut) +
  geom_point(aes(x = Axis1, y = Axis2, colour = pop), size = 3, alpha = 0.8) +
  scale_colour_manual(values = pca_col$col) +
  theme_minimal() +
  labs(x = paste("", "Axis 1: ", eig_percent_neut[1], "%"), y = paste("", "Axis 2: ", eig_percent_neut[2], "%")) +
  coord_fixed(ratio = 1) +
  theme(legend.title = element_blank())

# Plot a PCA for the outlier loci
x_out <- scaleGen(out, NA.method = "mean")
out_pca <- dudi.pca(x_out, center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

eig_percent_out <- round((out_pca$eig / (sum(out_pca$eig))) * 100, 2)

pca_dat_out <- as.data.frame(out_pca$li) %>%
  rownames_to_column(var = "ind") %>%
  cbind(pop = pop(out)) %>%
  mutate(pop = factor(pop, levels = c("ESP", "SW", "SE", "NW", "NE", "SLD", "SNO", "NNO"))) %>%
  mutate(type = "out", loc.type = "haps", location = "all")

out_pca <- ggplot(data = pca_dat_out) +
  geom_point(aes(x = Axis1, y = Axis2, colour = pop), size = 3, alpha = 0.8) +
  scale_colour_manual(values = pca_col$col) +
  theme_minimal() +
  labs(x = paste("", "Axis 1: ", eig_percent_out[1], "%"), y = paste("", "Axis 2: ", eig_percent_out[2], "%")) +
  coord_fixed(ratio = 1) +
  theme(legend.title = element_blank())

##### Sample Map #####
coords <- read.csv(here::here("data", "raw", "coords.csv"))

coords$group <- as.factor(c(1, 2, 3, 4, 5, 6, 7, 8))

world <- map_data("world") %>%
  dplyr::rename(X = long, Y = lat, PID = group, POS = order) 

# Clip the polygons in the region of interest to avoid problems with reconnecting lines
world <- clipPolys(world, xlim = c(-12, 13), ylim = c(37, 66), keepExtra = TRUE)

#svg(filename = "out/fig/raw/sample_map.svg", width = 9, height = 6, pointsize = 12, bg = "white")

sample_map <- ggplot() +
  geom_polygon(data = world, aes(x = X, y = Y, group = PID), alpha = 0.75, colour = "white", size = 0.05) +
  geom_point(data = coords, aes(x = long, y = lat, colour = group), size = 4) +
  geom_text(data = coords, aes(label = site, x = long, y = lat), fontface = "bold",
            size = 3,
            nudge_x = c(-1.5, -1.5, -2, 1.5, 1.5, -1.9, 0, -2),
            nudge_y = c(0, 0, 0, -0.4, 0.4, 0.5, 1, 0)) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        plot.margin = margin(0.5, 0.1, 0, 0.5, "cm")) +
  scale_colour_manual(values = pca_col$col) +
  theme(legend.position = "none") +
  labs(x = "",
       y = "") +
  geom_segment(aes(x = 2.5, xend = 3.3, y = 63, yend = 57.5), size = 2, color = "red", linetype = 3) +
  geom_text(aes(x = 5.4, y = 67, label = "Norwegian", fontface = "bold")) +
  geom_text(aes(x = -7.5, y = 64, label = "Atlantic", fontface = "bold"))

map_no_lines <- ggplot() +
  geom_polygon(data = world, aes(x = X, y = Y, group = PID), alpha = 0.75, colour = "white", size = 0.05) +
  geom_point(data = coords, aes(x = long, y = lat, colour = group), size = 4) +
  geom_text(data = coords, aes(label = site, x = long, y = lat), fontface = "bold",
            size = 3,
            nudge_x = c(-1.5, -1.5, -2, 1.5, 1.5, -1.9, -1.8, -2),
            nudge_y = c(0, 0, 0, -0.4, 0.4, 0.5, 0.3, 0)) + 
  annotate("text", x = -4, y = 40, label = "SPAIN", color = "lightgrey", family = "Open Sans") +
  annotate("text", x = -2, y = 52.5, label = "U.K.", color = "lightgrey", family = "Open Sans") +
  annotate("text", x = 10, y = 62, label = "NORWAY", color = "lightgrey", family = "Open Sans") +
  theme_void() +
  theme(plot.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(0.5, 0.1, 0, 0.5, "cm")) +
  scale_colour_manual(values = pca_col$col) +
  theme(legend.position = "none") +
  labs(x = "",
       y = "")

map_no_points <- ggplot() +
  geom_polygon(data = world, aes(x = X, y = Y, group = PID), alpha = 0.75, colour = "white", size = 0.05) +
  annotate("text", x = -4, y = 40, label = "SPAIN", color = "lightgrey", family = "Open Sans") +
  annotate("text", x = -2, y = 52.5, label = "U.K.", color = "lightgrey", family = "Open Sans") +
  annotate("text", x = 10, y = 62, label = "NORWAY", color = "lightgrey", family = "Open Sans") +
  theme_void() +
  theme(plot.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(0.5, 0.1, 0, 0.5, "cm")) +
  theme(legend.position = "none") +
  labs(x = "",
       y = "")


#comb_plot <- sample_map + neut_pca / out_pca + plot_layout(ncol = 2, widths = c(4, 3))

comb_plot <- (neut_pca / out_pca) - sample_map +
  plot_layout(ncol = 2, widths = c(3, 4), guides = "collect") +
  plot_annotation(tag_levels = 'A')

ggsave("out/fig/sample_map.png", comb_plot, width = 20, height = 15, units = "cm")


