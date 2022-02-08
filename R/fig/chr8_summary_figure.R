library(tidyverse)
library(patchwork)
library(vcfR)
library(adegenet)
library(gaston)
library(vcfR)
library(gridExtra)

source("bin/functions.R")

# Read in the relevant data

outlier_tbl <- read_rds(here::here("data", "derived", "outlier_summary.rds"))
haplo_div_tbl <- read_rds(here::here("data", "derived", "haplotype_div.rds"))
gen_div_tbl <- read_rds(here::here("data", "derived", "genetic_div.rds"))
hapsel_tbl <- read_rds(here::here("data", "derived", "hap_selection_tbl.rds"))
rsb_tbl <- read_rds(here::here("data", "derived", "rsb_tbl.rds"))
scot_nor_pw <- read_rds(here::here("data", "derived", "scot_nor_pw_fst.rds"))
scot_esp_pw <- read_rds(here::here("data", "derived", "scot_esp_pw_fst.rds"))
esp_nor_pw <- read_rds(here::here("data", "derived", "esp_nor_pw_fst.rds"))

pca_col <- read_tsv(file = here::here("data", "raw", "colors.tsv"))

contig_bed <- read_tsv(here::here("data", "raw", "genome", "reference.fasta.fai"), col_names = FALSE)

chrom_vec <- contig_bed %>%
  tidyr::extract(X1, c("chrom"), "HiC_scaffold_(\\d+)_arrow_ctg1", remove = FALSE) %>%
  mutate(chrom = as.integer(chrom)) %>%
  filter(chrom <= 19) %>%
  pull(X1) %>%
  unique()

het_smooth <- chrom_vec %>%
  map(function(chrom) {
    #browser()
    gen_div_tbl %>%
      rename(chrom_num = chrom) %>%
      extract(locus, "CHR", "(HiC_scaffold_\\d+_arrow_ctg1)", remove = FALSE) %>%
      filter(CHR == chrom) %>%
      split(.$pop) %>%
      map(function(df) {
        
        df <- arrange(df, pos)
        he <- runner::runner(df$he, f = function(x){mean(x, na.rm = TRUE)}, k = 5)
        df$he_avg <- he
        
        ho <- runner::runner(df$ho, f = function(x){mean(x, na.rm = TRUE)}, k = 5)
        df$ho_avg <- ho
        
        return(df)
      }) %>%
      bind_rows()
  }) %>%
  bind_rows() %>%
  rename(chrom = CHR)

# Create the left side of the panel

env_tbl <- outlier_tbl %>%
  select(locus, chrom, pos, any_env_outlier)


scot_esp_pw <- scot_esp_pw %>%
  left_join(env_tbl, by = c("locus", "chrom", "pos"))

chr8_het <- het_smooth %>%
  filter(chrom_num == 8) %>%
  mutate(pop = factor(pop, levels = c("ESP", "SW", "SE", "NW", "NE", "SLD", "SNO", "NNO"))) %>%
  ggplot(aes(x = pos, y = he_avg, color = pop)) +
  geom_line() +
  scale_colour_manual(values = pca_col$col) +
  labs(y = "heterozygosity") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

chr8_fst <- scot_esp_pw %>%
  filter(chrom == 8) %>%
  ggplot(aes(x = pos, y = Fst, col = outlier)) +
  geom_point() +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_color_manual(values = c(cmh_palette[1], cmh_palette[3])) +
  labs(y = expression(F[ST]))

chr8_env <- scot_esp_pw %>%
  filter(chrom == 8) %>%
  ggplot(aes(x = pos, y = Fst, col = any_env_outlier)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c(cmh_palette[1], cmh_palette[4]))

chr8_ies <- hapsel_tbl %>%
  filter(chrom == 8) %>%
  rename(pop = site) %>%
  mutate(pop = factor(pop, levels = c("ESP", "SW", "SE", "NW", "NE", "SLD", "SNO", "NNO"))) %>%
  ggplot(aes(x = position, y = ies, color = pop)) +
  geom_line() +
  scale_colour_manual(values = pca_col$col) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  labs(y = "iES")

comp_labs <- c("Spain :: Norway", "Spain :: Scotland", "Scotland :: Norway")

chr8_rsb <- rsb_tbl %>%
  filter(chrom_num == 8) %>%
  ggplot(aes(x = pos / 1000000, y = logpvalue, color = comp)) +
  geom_line() +
  scale_color_manual(name = "Comparison",
                     labels = comp_labs, 
                     values = c(cmh_palette[3], cmh_palette[2], cmh_palette[5])) +
  labs(y = expression(log(P)[Rsb]),
       x = "Position (Mb)") +
  #labs(y = expression(paste("log(P) ",  [Rsb], sep = " "))) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.margin=margin()) +
  guides(color = guide_legend(nrow=2,byrow=TRUE))





# Calculate LD for the chromosome

# Read in the genind object with all of the SNPs
#gen_snps <- read_rds(here::here("data", "derived", "out_snps_gen.rds"))
gen <- read.vcfR(here::here("data", "derived", "out.17.phased.vcf")) %>%
  vcfR2genind()

pops <- tibble(ind = indNames(gen)) %>%
  extract(ind, "pop", "(\\w+)_", remove = FALSE) %>%
  pull(pop)

pop(gen) <- pops

# Combine Scotland pops and check LD
new_pop <- tibble(pop = as.character(pop(gen))) %>%
  mutate(new_pop = case_when(
    pop %in% c("SW", "NW", "SE", "NE", "SLD") ~ "SCOT",
    pop %in% c("SNO", "NNO") ~ "NOR",
    pop %in% c("ESP") ~ "ESP")) %>%
  pull(new_pop)

pop(gen) <- new_pop


snp_tbl <- tibble(locus = locNames(gen)) %>%
  extract(locus, c("scaffold", "pos"), "(HiC_scaffold_\\d+_arrow_ctg1)_(\\d+)", remove = FALSE) %>%
  extract(locus, "scaffold_num", "scaffold_(\\d+)_arrow_ctg1", remove = FALSE) %>%
  mutate(scaffold_num = as.character(scaffold_num)) %>%
  mutate(pos = as.integer(pos)) %>%
  select(snp_id = locus, chr = scaffold_num, pos)

bed_lst <- gen %>%
  seppop() %>%
  map(function(x) {
    
    pop_name <- popNames(x)
    df <- genind2df(x) %>%
      mutate_all(~str_replace_all(., "00", "0")) %>%
      mutate_all(~str_replace_all(., "01", "1")) %>%
      mutate_all(~str_replace_all(., "10", "1")) %>%
      mutate_all(~str_replace_all(., "11", "2")) %>%
      select(-pop) %>%
      as.matrix()
    gl <- as.genlight(df)
    gafftop::genlight_to_bed(gl, snp_tbl) 
    
  })

ld_all <- 1:length(bed_lst) %>%
  map(function(i) {
    
    pop <- names(bed_lst[i])
    bed <- bed_lst[[i]]
    
    ld <- LD(bed, c(1, 1920))
    #ld[is.nan(ld)] <- 0
    
    loc_order <- rownames(ld)
    
    
    ld_tbl <- ld %>% 
      as.data.frame() %>%
      rownames_to_column("locus1") %>%
      pivot_longer(-c(locus1), names_to = "locus2", values_to = "r2") %>%
      mutate(locus1 = fct_relevel(locus1, loc_order)) %>%
      mutate(locus2 = fct_relevel(locus2, loc_order)) %>%
      mutate(pop = pop) %>%
      extract(locus1, "chrom1", "HiC_scaffold_(\\d+)_arrow_ctg1", remove = FALSE) %>%
      extract(locus2, "chrom2", "HiC_scaffold_(\\d+)_arrow_ctg1", remove = FALSE)
    
    
    #ldplot <- do.call(grid.arrange, gglist)
    
    #ggsave(filename = here::here("out", glue::glue("{pop}_ld.pdf")), plot = ldplot )
  }) %>%
  bind_rows()

ld_8 <- ld_all %>%
  filter(chrom1 == 8, chrom2 == 8)

gg_8 <- ld_8 %>%
  split(.$pop) %>%
  map(function(x) {
    
    pop_name <- x$pop[1]
    
    plot_title <- case_when(pop_name == "NOR" ~ "Norway",
                            pop_name == "SCOT" ~ "Scotland",
                            pop_name == "ESP" ~ "Spain")
    
    gg <- ggplot(x, aes(x=locus1, y=locus2, fill=r2)) + 
      geom_raster() + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank()) +
      coord_fixed() +
      scale_fill_viridis_c(limits = c(0, 1)) +
      labs(title = plot_title,
           x = "Chromosome 8") +
      guides(fill = guide_legend(title = expression(r^2)))
    
    return(gg)
  })

plot1 <- (chr8_fst + chr8_het + chr8_ies + chr8_rsb + plot_layout(ncol = 1) + 
            plot_annotation(tag_levels = 'A'))  / 
  (gg_8$NOR + gg_8$SCOT + gg_8$ESP + plot_layout(ncol = 3, guides = "collect"))


ggsave(filename = "out/fig/chrom8_summary.png", plot = plot1, width = 15, height = 22,
       units = "cm")
