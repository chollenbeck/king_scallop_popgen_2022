library(tidyverse)
library(stringr)

# Filter genotypes with depth < 10 and genotype quality < 20
system("vcftools --gzvcf raw_snps.vcf.gz --out out.1 --minDP 10 --minGQ 20 --recode --recode-INFO-all")

# Filter out sites that were made monomorphic by the previous filter
system("vcftools --vcf out.1.recode.vcf --maf 0.001 --out out.2 --recode --recode-INFO-all")

# Remove sites with more than 50% missing data
system("vcftools --vcf out.2.recode.vcf --out out.3 --max-missing 0.5 --recode --recode-INFO-all")

# Produce a file with missingness per individual
system("vcftools --vcf out.3.recode.vcf --out out.3 --missing-indv")

# Load the data for the missingness file
out_3_imiss <- read_tsv("out.3.imiss")

# First get a list of all of the duplicate individuals so that these can be saved from filtering in this step
# All of the duplicates have a 'b' at the end of the individual name

dups_a <- str_subset(out_3_imiss$INDV, "b$")
dups_b <- gsub(pattern = "b", replacement = "", x = dups_a)
dups <- c(dups_a, dups_b)

# Plot a quick histogram of the data
qplot(out_3_imiss$F_MISS)

# Select individuals with more than 50% missing data
miss_70 <- filter(out_3_imiss, F_MISS > 0.7, ! INDV %in% dups) %>%
            select(INDV)

# Write the individuals to remove to a file
write_delim(miss_70, "remove.3.inds", col_names = FALSE)

# Also write the duplicate individuals to a file
write_delim(as.data.frame(dups_a), "dups_a.txt", col_names = FALSE)
write_delim(as.data.frame(dups_b), "dups_b.txt", col_names = FALSE)

# Remove individuals with >70% missing data
system("vcftools --vcf out.3.recode.vcf --out out.4 --remove remove.3.inds --recode --recode-INFO-all")

# Identify site discordance among duplicate individuals
system("vcftools --vcf out.4.recode.vcf --out dups.a --keep dups_a.txt --recode --recode-INFO-all")
system("vcftools --vcf out.4.recode.vcf --out dups.b --keep dups_b.txt --recode --recode-INFO-all")

# Remove the 'b' the names of the duplicate samples so that they can be compared
system("perl -p -i -e 's/_(\\d{3})b/_$1/g' dups.a.recode.vcf")

system("vcftools --vcf dups.a.recode.vcf --diff dups.b.recode.vcf --diff-site-discordance --out out.4")

# Read in the site discordance file
discord_raw <- read_tsv("out.4.diff.sites")

disc_sites <- filter(discord_raw, DISCORDANCE > 0) %>%
  select(CHROM, POS)

# Write the sites to remove to a file
write_delim(disc_sites, "remove.4.sites", col_names = FALSE)

# Remove discordant sites
system("vcftools --vcf out.4.recode.vcf --out out.5 --exclude-positions remove.4.sites --recode --recode-INFO-all")

# Remove duplicate individuals
system("vcftools --vcf out.5.recode.vcf --out out.6 --remove dups_a.txt --recode --recode-INFO-all")

# Remove loci with overall quality < 20
system("vcffilter -s -f 'QUAL > 20' out.6.recode.vcf > out.7.vcf")


# Calculate site depth
system("vcftools --vcf out.7.vcf --site-depth --out out.7")

# Read in the site depth file
site_depth_7 <- read_tsv("out.7.ldepth") %>%
                  mutate(MEAN_DEPTH = SUM_DEPTH / 217)

# Plot a histogram of the mean site depth per individual
qplot(site_depth_7$MEAN_DEPTH)

# Filter out loci with a mean site depth > 300
mean_site_depth_7 <- mean(site_depth_7$MEAN_DEPTH)
to_keep_7 <- filter(site_depth_7, MEAN_DEPTH < 300)
mean_site_depth_7_filt <- mean(to_keep_7$MEAN_DEPTH)

# Plot the distribution again
qplot(to_keep_7$MEAN_DEPTH)

# Make a list of the sites to filter
to_filter_7 <- filter(site_depth_7, MEAN_DEPTH >= 300) %>%
                  select(CHROM, POS)

# Write the sites to remove to a file
write_delim(to_filter_7, "remove.7.sites", col_names = FALSE)

# Remove the sites with VCFtools
system("vcftools --vcf out.7.vcf --out out.8 --exclude-positions remove.7.sites --recode --recode-INFO-all")

# Remove sites with more than 75% missing data
system("vcftools --vcf out.8.recode.vcf --out out.9 --max-missing 0.75 --recode --recode-INFO-all")

# Calculate individual missingness
system("vcftools --vcf out.9.recode.vcf --out out.9 --missing-indv")
       
# Load the data for the out.9.recode.vcf file
out_9_imiss <- read_tsv("out.9.imiss") %>%
  extract(INDV, "pop", "(\\w+)_", remove = FALSE)

out_9_imiss %>%
  count(pop)

# Plot a quick histogram of the data
ggplot(out_9_imiss, aes(x = F_MISS, fill = pop)) +
  geom_histogram()

# Select individuals with more than 40% missing data
miss_40 <- filter(out_9_imiss, F_MISS > 0.40) %>%
  select(INDV)

# Write the individuals to remove to a file
write_delim(miss_40, "remove.9.inds", col_names = FALSE)

# Remove the individuals with more than 40% missing genotypes
system("vcftools --vcf out.9.recode.vcf --out out.10 --remove remove.9.inds --recode --recode-INFO-all")

# Check for duplicate individuals using relatedness
system("vcftools --vcf out.10.recode.vcf --relatedness --out out.10")

# Load the data for the out.10.recode.vcf file
out_10_relate <- read_tsv("out.10.relatedness")

# Check for abnormally high relatedness
relatives <- filter(out_10_relate, INDV1 != INDV2) %>% 
              filter(RELATEDNESS_AJK > 0.7)

# Make a list of individuals that may be contaminated or clones
inds_to_remove <- rbind(relatives$INDV1, relatives$INDV2) %>% as.data.frame()

# Write the individuals to remove to a file
write_delim(inds_to_remove, "remove.10.inds", col_names = FALSE)

# Remove the individuals
system("vcftools --vcf out.10.recode.vcf --out out.11 --remove remove.10.inds --recode --recode-INFO-all")

# Remove sites with more than 10% missing data in any population
#system("sh pop_missing_filter.sh out.11.recode.vcf popmap 0.1 1 out.12")

# Decompose sites to allelic primitives
system("vcfallelicprimitives -k -g out.11.recode.vcf > out.12.vcf")

# Remove indels
system("vcftools --vcf out.12.vcf --out out.13 --remove-indels --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all")

# Filter for minor allele frequency
system("vcftools --vcf out.13.recode.vcf --out out.14 --maf 0.05 --recode --recode-INFO-all")

# Unphase phased sites
system("sed '/^##/! s/|/\//g' out.14.recode.vcf > out.14.unphased.vcf")

# Convert missing genotypes from '.' to './.'
system("cat out.14.unphased.vcf | perl -pe 's/\s\.:/\t.\/.:/g' > out.14.missrecode.vcf")

# Filter based on HWE in each population
system("perl filter_hwe_by_pop.pl -v out.14.recode.vcf --out out.15 -p popmap -h 0.001 -c 0.1")

# Thin the SNPs based on genomic position
system("bcftools +prune -w 1000bp -n 1 out.15.recode.vcf -o out.16.vcf")

# Get a list of the scaffolds and positions
system("grep -v '#' out.16.vcf | cut -f1,2 > all_sites_16.txt")

read_tsv("all_sites_16.txt", col_names = FALSE) %>%
  extract(X1, "chrom", "HiC_scaffold_(\\d+)_", remove = FALSE) %>%
  mutate(chrom = as.integer(chrom)) %>%
  filter(chrom <= 19) %>%
  select(X1, X2) %>%
  write_tsv("sites_to_keep_16.txt", col_names = FALSE)

system("vcftools --vcf out.16.vcf --out out.17 --positions sites_to_keep_16.txt --recode --recode-INFO-all")

