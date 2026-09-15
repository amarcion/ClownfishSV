library(tidyverse)

setwd("/Users/amarcion/Documents/PacBio_Sequencing_Assembly/10_Additional/15_Long_Reads_Mapping/")

PATH_SELF = "Statistics_windowsHet/"
OUTPATH="Heterozygosity_ByWindows/"


df_self <- read_csv(paste0(PATH_SELF,"all_het_windows_selfmapped.csv"))

# Filter low-quality windows
df_filt <- df_self %>% filter(total_snps >= 10, mean_depth >= 10, mean_depth <= 500)

# Mean heterozygosity per species
mean_het <- df_filt %>%
  group_by(species) %>%
  summarise(
    mean_het = mean(het_per_bp, na.rm = TRUE),
    sd_het   = sd(het_per_bp, na.rm = TRUE),
    median_het = median(het_per_bp, na.rm = TRUE),
    n_windows = n()
  ) %>%
  arrange(mean_het)

print(mean_het)

# Barplot of mean heterozygosity
ggplot(mean_het, aes(x = reorder(species, mean_het), y = mean_het)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  # geom_errorbar(aes(ymin = mean_het - sd_het,
  #                   ymax = mean_het + sd_het), width = 0.3) +
  labs(x = "Species", y = "Mean heterozygosity (het SNPs / bp)",
       title = "True individual heterozygosity per species") +
  theme_bw() +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(OUTPATH,"mean_heterozygosity_ByRef.pdf"), width = 7, height = 5)

# Genome-wide plot per species
for (SPECIES in unique(df_filt$species)) {
df_filt %>%
  mutate(chrom = factor(chrom)) %>%
  filter(species==SPECIES)%>%
  ggplot(aes(x = start / 1e6, y = het_per_bp)) +
  geom_line(linewidth = 0.3, color = "steelblue") +
  facet_wrap(~ chrom,  ncol = 5) +
  labs(x = "Position (Mb)", y = "Heterozygosity (het SNPs / bp)") +
  theme_bw(base_size = 7) +
  theme(strip.text.x = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(OUTPATH,"ByRef_eachSpecies/",SPECIES, ".het_per_bp_selfmapped.pdf"),
        width = 10, height = 8)
  }

# Genome-wide plot per species - Only 24 main scaffolds
for (SPECIES in unique(df_filt$species)) {
  df_filt %>%
    filter(grepl("^scf[0-9]+\\.1$", chrom) | grepl("chr", chrom))%>%
    mutate(chrom = factor(chrom)) %>%
    filter(species==SPECIES)%>%
    ggplot(aes(x = start / 1e6, y = het_per_bp)) +
    geom_line(linewidth = 0.3, color = "steelblue") +
    facet_wrap(~ chrom,  ncol = 5) +
    labs(x = "Position (Mb)", y = "Heterozygosity (het SNPs / bp)") +
    theme_bw(base_size = 7) +
    theme(strip.text.x = element_text(size = 5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(OUTPATH,"ByRef_eachSpecies/",SPECIES, ".MainScaff.het_per_bp_selfmapped.pdf"),
         width = 10, height = 8)
}

# This only with first 24
df_filt_renamed <- df_filt %>%
  filter(
    str_detect(chrom, "^scf[0-9]+\\.1$") |
      str_detect(chrom, "chr")
  ) %>%
  mutate(
    chrom = case_when(
      # scf01.1 -> chr1, scf02.1 -> chr2, ...
      str_detect(chrom, "^scf[0-9]+\\.1$") ~
        paste0("chr", as.integer(str_extract(chrom, "[0-9]+"))),
      
      # Aocel_chr7 -> chr7, Aocel_chr12 -> chr12, ...
      str_detect(chrom, "chr") ~
        str_extract(chrom, "chr[0-9A-Za-z]+"),
      
      TRUE ~ chrom))

pdf(paste0(OUTPATH,"het_per_bp_all_ByRef.pdf"), width = 20, height = 15)
df_filt_renamed %>%
  filter(grepl("^scf[0-9]+\\.1$", chrom) | grepl("chr", chrom))%>%
  mutate(chrom = factor(chrom)) %>%
  ggplot(aes(x = start / 1e6, y = het_per_bp)) +
  geom_line(linewidth = 0.3, color = "steelblue") +
  facet_grid(species ~ chrom, scales = "free_x", space = "free_x") +
  labs(x = "Position (Mb)", y = "Heterozygosity (het SNPs / bp)") +
  theme_bw(base_size = 7) +
  theme(strip.text.x = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Check mean coverage per species
df_filt_renamed %>%
  filter(grepl("^scf[0-9]+\\.1$", chrom) | grepl("chr", chrom))%>%
  mutate(chrom = factor(chrom)) %>%
  ggplot(aes(x = start / 1e6, y = mean_depth)) +
  geom_line(linewidth = 0.3, color = "steelblue") +
  facet_grid(species ~ chrom, scales = "free_x", space = "free_x") +
  labs(x = "Position (Mb)", y = "Heterozygosity (het SNPs / bp)") +
  theme_bw(base_size = 7) +
  theme(strip.text.x = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(OUTPATH,"meandepth_Max500_ByRef.pdf"), width = 20, height = 15)

