#!/usr/bin/env Rscript
# =============================================================================
# Plotting script for vcftools + mosdepth + bcftools stats outputs
# One VCF per species (self-mapped: AKA on AKA, AKY on AKY, etc.)
# =============================================================================

library(tidyverse)  
library(patchwork) 
library(ggplot2) 

# =============================================================================
# PARAMETERS - adjust to your setup
# =============================================================================

SPECIES <- c("AKA", "AKY", "ALL", "BIA", "CRP", "EPH", "FRE",
             "LAT", "LAZ", "MCC", "OMA", "POL", "PRC", "PRD", "SAN", "SEB")

# Get data on curnagl
STATSDIR <- "Stats_ByReference/"
OUTDIR   <- "SNP_clair3_BySpecies/plots"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# 
# Scaffolds to visualize in details
CHR18 <- c(
  AKA = "scf18.1", AKY = "scf18.1", ALL = "scf18.1", BIA = "scf18.1",
  CRP = "scf18.1", EPH = "scf18.1", FRE = "scf18.1", LAT = "scf18.1",
  LAZ = "scf18.1", MCC = "scf18.1", OMA = "scf18.1", POL = "scf18.1",
  PRC = "scf18.1", PRD = "scf18.1", SAN = "scf18.1", SEB = "scf18.1"
)

CHR01 <- c(
  AKA = "scf01.1", AKY = "scf01.1", ALL = "scf01.1", BIA = "scf01.1",
  CRP = "scf01.1", EPH = "scf01.1", FRE = "scf01.1", LAT = "scf01.1",
  LAZ = "scf01.1", MCC = "scf01.1", OMA = "scf01.1", POL = "scf01.1",
  PRC = "scf01.1", PRD = "scf01.1", SAN = "scf01.1", SEB = "scf01.1"
)

CHR03 <- c(
  AKA = "scf03.1", AKY = "scf03.1", ALL = "scf03.1", BIA = "scf03.1",
  CRP = "scf03.1", EPH = "scf03.1", FRE = "scf03.1", LAT = "scf03.1",
  LAZ = "scf03.1", MCC = "scf03.1", OMA = "scf03.1", POL = "scf03.1",
  PRC = "scf03.1", PRD = "scf03.1", SAN = "scf03.1", SEB = "scf03.1"
)


# Color palette for species (one color per species)
SP_COLORS <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(SPECIES)),
  SPECIES
)

# =============================================================================
# READ DATA
# =============================================================================

message("[INFO] Reading vcftools and mosdepth outputs...")

# --- SNP density (vcftools --SNPdensity) ---
# Columns: CHROM, BIN_START, SNP_COUNT, VARIANTS/KB
density_all <- map_dfr(SPECIES, function(sp) {
  f <- file.path(STATSDIR, sp, paste0(sp, ".snpden"))
  if (!file.exists(f)) { message("[WARNING] Missing: ", f); return(NULL) }
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(species = sp)
})

# --- Heterozygosity (vcftools --het) ---
# Columns: INDV, O(HOM), E(HOM), N_SITES, F
het_summary <- map_dfr(SPECIES, function(sp) {
  f <- file.path(STATSDIR, sp, paste0(sp, ".het"))
  if (!file.exists(f)) { message("[WARNING] Missing: ", f); return(NULL) }
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(species = sp)
})

# --- Allele frequency (vcftools --freq2) ---
# Columns: CHROM, POS, N_ALLELES, N_CHR, {FREQ}
af_all <- map_dfr(SPECIES, function(sp) {
  f <- file.path(STATSDIR, sp, paste0(sp, ".frq"))
  if (!file.exists(f)) { message("[WARNING] Missing: ", f); return(NULL) }
  read_tsv(f, show_col_types = FALSE,
           col_names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "REF_FREQ", "ALT_FREQ"),
           skip = 1) %>%
    mutate(species = sp,
           ALT_FREQ = as.numeric(ALT_FREQ))
})

# --- Per-site depth (vcftools --site-depth) ---
# Columns: CHROM, POS, SUM_DEPTH, SUMSQ_DEPTH
depth_site <- map_dfr(SPECIES, function(sp) {
  f <- file.path(STATSDIR, sp, paste0(sp, ".ldepth"))
  if (!file.exists(f)) { message("[WARNING] Missing: ", f); return(NULL) }
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(species = sp)
})

# --- Mosdepth windowed depth ---
# Columns: scaffold, start, end, depth
depth_wins <- map_dfr(SPECIES, function(sp) {
  f <- file.path(STATSDIR, sp, paste0(sp, "_mosdepth.regions.bed.gz"))
  if (!file.exists(f)) { message("[WARNING] Missing: ", f); return(NULL) }
  read_tsv(f, col_names = c("scaffold", "start", "end", "depth"),
           show_col_types = FALSE) %>%
    mutate(species = sp,
           mid = (start + end) / 2)
})

# --- bcftools stats ---
# Parse the SN (summary numbers) section for Ti/Tv and other stats
parse_bcftools_stats <- function(sp) {
  f <- file.path(STATSDIR, sp, paste0(sp, "_bcftools_stats.txt"))
  if (!file.exists(f)) { message("[WARNING] Missing: ", f); return(NULL) }
  lines <- readLines(f)
  sn    <- lines[grepl("^SN", lines)]
  vals  <- strsplit(sn, "\t")
  df_1    <- data.frame(
    key   = sapply(vals, `[`, 3),
    value = as.numeric(sapply(vals, `[`, 4)),
    stringsAsFactors = FALSE
  )
  df_1$species <- sp
  # tstv
  tstv <- lines[grepl("^TSTV", lines)]
  vals_tstv <- strsplit(tstv, "\t")
  tstv_ratio <- as.numeric(vals_tstv[[1]][5])
  tstv_ratio_alt <- as.numeric(vals_tstv[[1]][8])
  df_2 <- data.frame(
    key = c("tstv"),
    value=c(tstv_ratio),
    species=c(sp))
  df <- rbind(df_1, df_2)
  df
}

bcftools_stats <- map_dfr(SPECIES, parse_bcftools_stats)

message("[INFO] All data loaded.")

# =============================================================================
# PLOT 1: Summary heterozygosity per species (from vcftools --het)
# F statistic: positive = excess homozygosity, negative = excess heterozygosity
# Het rate = 1 - O(HOM)/N_SITES
# =============================================================================

message("[INFO] Plot 1: Heterozygosity summary per species...")

het_summary <- het_summary %>%
  rename(O_HOM = `O(HOM)`, E_HOM = `E(HOM)`) %>%
  mutate(het_rate = 1 - O_HOM / N_SITES)

p1 <- ggplot(het_summary, aes(x = reorder(species, het_rate), 
                              y = het_rate, fill = species)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(het_rate, 4)), hjust = -0.1, size = 3) +
  scale_fill_manual(values = SP_COLORS) +
  coord_flip() +
  expand_limits(y = max(het_summary$het_rate) * 1.15) +
  labs(
    title    = "Heterozygosity rate per species",
    subtitle = "Computed from vcftools --het; rate = 1 - O(HOM)/N_SITES",
    x        = NULL,
    y        = "Heterozygosity rate"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUTDIR, "01_heterozygosity_per_species.pdf"),
       p1, width = 8, height = 6)

# =============================================================================
# PLOT 2: Inbreeding coefficient F per species
# =============================================================================

message("[INFO] Plot 2: Inbreeding coefficient F...")

p2 <- ggplot(het_summary, aes(x = reorder(species, F), 
                              y = F, fill = F > 0)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4393c3")) +
  coord_flip() +
  labs(
    title    = "Inbreeding coefficient (F) per species",
    subtitle = "Positive F = excess homozygosity; negative F = excess heterozygosity",
    x        = NULL,
    y        = "F statistic"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUTDIR, "02_inbreeding_F_per_species.pdf"),
       p2, width = 8, height = 6)

# =============================================================================
# PLOT 3: SNP density along scaffolds - one panel per scaffold per species
# =============================================================================

message("[INFO] Plot 3: SNP density along scaffolds...")

# Get top scaffolds (by total SNP count across all species)
top_scaffolds <- density_all %>%
  group_by(CHROM) %>%
  summarise(total = sum(SNP_COUNT, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 24) %>%
  pull(CHROM)

for (sp in SPECIES) {
  
  sub <- density_all %>%
    filter(species == sp, CHROM %in% top_scaffolds)
  
  if (nrow(sub) == 0) next
  
  # Order scaffolds by name
  sub$CHROM <- factor(sub$CHROM, levels = top_scaffolds)
  
  p <- ggplot(sub, aes(x = BIN_START / 1e6, y = SNP_COUNT)) +
    geom_col(fill = SP_COLORS[sp], width = 0.1) +
    facet_wrap(~ CHROM, scales = "free_x", ncol = 4) +
    labs(
      title    = paste("SNP density along scaffolds —", sp),
      subtitle = "100 kb windows (vcftools --SNPdensity)",
      x        = "Position (Mb)",
      y        = "SNP count per 100 kb"
    ) +
    theme_bw(base_size = 9) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(size = 7),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(OUTDIR, paste0("03_snp_density_", sp, ".pdf")),
         p, width = 16, height = 12)
}

# =============================================================================
# PLOT 4: SNP density heatmap - all species x scaffolds
# =============================================================================

message("[INFO] Plot 4: SNP density heatmap all species...")

p4 <- density_all %>%
  filter(CHROM %in% top_scaffolds) %>%
  group_by(species, CHROM) %>%
  summarise(mean_density = mean(SNP_COUNT, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = CHROM, y = species, fill = mean_density)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", name = "Mean SNPs\nper 100kb") +
  labs(
    title    = "Mean SNP density across scaffolds and species",
    subtitle = "Top 24 scaffolds by total SNP count",
    x        = "Scaffold",
    y        = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid   = element_blank()
  )

ggsave(file.path(OUTDIR, "04_snp_density_heatmap.pdf"),
       p4, width = 14, height = 7)

# =============================================================================
# PLOT 5: Allele frequency distribution per species
# =============================================================================

message("[INFO] Plot 5: Allele frequency distributions...")

p5 <- ggplot(af_all[!is.na(af_all$ALT_FREQ), ],
             aes(x = ALT_FREQ)) +
  geom_histogram(bins = 50, fill = "#4393c3", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.5, color = "#d73027", 
             linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ species, scales = "free_y", ncol = 4) +
  labs(
    title    = "Allele frequency distribution per species",
    subtitle = "Dashed line = 0.5 (expected for true heterozygous sites)",
    x        = "Alternative allele frequency",
    y        = "Count"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUTDIR, "05_allele_frequency_distribution.pdf"),
       p5, width = 16, height = 12)

# =============================================================================
# PLOT 6: Depth distribution per species (from vcftools --site-depth)
# =============================================================================

message("[INFO] Plot 6: Depth distributions...")

p6 <- ggplot(depth_site, aes(x = SUM_DEPTH)) +
  geom_histogram(bins = 80, fill = "#74add1", color = "white", linewidth = 0.2) +
  facet_wrap(~ species, scales = "free", ncol = 4) +
  labs(
    title = "Per-site depth distribution per species",
    x     = "Depth",
    y     = "Count"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUTDIR, "06_depth_distribution.pdf"),
       p6, width = 16, height = 12)

# =============================================================================
# PLOT 7: Mosdepth windowed depth along scaffolds per species
# =============================================================================

message("[INFO] Plot 7: Mosdepth windowed depth along scaffolds...")

for (sp in SPECIES) {
  
  sub <- depth_wins %>%
    filter(species == sp, scaffold %in% top_scaffolds)
  
  if (nrow(sub) == 0) next
  
  sub$scaffold <- factor(sub$scaffold, levels = top_scaffolds)
  
  # Compute median depth for reference line
  med_depth <- median(sub$depth, na.rm = TRUE)
  
  p <- ggplot(sub, aes(x = mid / 1e6, y = depth)) +
    geom_line(color = SP_COLORS[sp], linewidth = 0.3, alpha = 0.8) +
    geom_hline(yintercept = med_depth, linetype = "dashed", 
               color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = med_depth * 1.5, linetype = "dotted",
               color = "#d73027", linewidth = 0.4) +
    facet_wrap(~ scaffold, scales = "free_x", ncol = 4) +
    labs(
      title    = paste("Windowed sequencing depth —", sp),
      subtitle = paste0("Dashed = median (", round(med_depth, 1), 
                        "x); dotted red = 1.5x median (potential collapsed repeats)"),
      x        = "Position (Mb)",
      y        = "Depth (100 kb windows)"
    ) +
    theme_bw(base_size = 9) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(size = 7),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(OUTDIR, paste0("07_depth_", sp, ".pdf")),
         p, width = 16, height = 12)
}

# =============================================================================
# PLOT 8: Ti/Tv ratio per species (from bcftools stats)
# =============================================================================

message("[INFO] Plot 8: Ti/Tv ratio per species...")

titv <- bcftools_stats %>%
  filter(grepl("tstv", key)) %>%
  mutate(key = trimws(key))

if (nrow(titv) > 0) {
  p8 <- ggplot(titv, aes(x = reorder(species, value), y = value, fill = species)) +
    geom_col(show.legend = FALSE) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = round(value, 2)), hjust = -0.1, size = 3) +
    scale_fill_manual(values = SP_COLORS) +
    coord_flip() +
    expand_limits(y = max(titv$value) * 1.15) +
    labs(
      title    = "Transition/Transversion ratio per species",
      subtitle = "Dashed line = expected Ti/Tv ~2 for vertebrates",
      x        = NULL,
      y        = "Ti/Tv ratio"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(file.path(OUTDIR, "08_titv_ratio.pdf"),
         p8, width = 8, height = 6)
}

# =============================================================================
# PLOT 9: Chromosome 18 — SNP density + depth across all species
# =============================================================================

message("[INFO] Plot 9: Chromosome 18 SNP density and depth...")

# SNP density on chr18
chr18_density <- map_dfr(SPECIES, function(sp) {
  scf <- CHR18[sp]
  density_all %>%
    filter(species == sp, CHROM == scf)
})

# Depth on chr18
chr18_depth <- map_dfr(SPECIES, function(sp) {
  scf <- CHR18[sp]
  depth_wins %>%
    filter(species == sp, scaffold == scf)
})

if (nrow(chr18_density) > 0) {
  
  p9a <- ggplot(chr18_density, 
                aes(x = BIN_START / 1e6, y = SNP_COUNT, color = species)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = SP_COLORS) +
    facet_wrap(~ species, ncol = 4) +
    labs(
      title    = "SNP density along chromosome 18",
      subtitle = "100 kb windows",
      x        = "Position (Mb)",
      y        = "SNP count per 100 kb"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  p9b <- ggplot(chr18_depth,
                aes(x = mid / 1e6, y = depth, color = species)) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    scale_color_manual(values = SP_COLORS) +
    facet_wrap(~ species, ncol = 4) +
    labs(
      title    = "Sequencing depth along chromosome 18",
      subtitle = "100 kb windows (mosdepth)",
      x        = "Position (Mb)",
      y        = "Depth"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  p9 <- p9a / p9b +
    plot_annotation(
      title    = "Chromosome 18 — SNP density and depth across species",
      subtitle = "Relevant for the large inversion described in the manuscript"
    )
  
  ggsave(file.path(OUTDIR, "09_chr18_snpdensity_depth.pdf"),
         p9, width = 16, height = 16)
}

# =============================================================================
# PLOT 9b: Chromosome 01 — SNP density + depth across all species for comparison
# =============================================================================

message("[INFO] Plot 9b: Chromosome 1 SNP density and depth...")

# SNP density on chr1
chr1_density <- map_dfr(SPECIES, function(sp) {
  scf <- CHR01[sp]
  density_all %>%
    filter(species == sp, CHROM == scf)
})

# Depth on chr1
chr1_depth <- map_dfr(SPECIES, function(sp) {
  scf <- CHR01[sp]
  depth_wins %>%
    filter(species == sp, scaffold == scf)
})

if (nrow(chr1_density) > 0) {
  
  p9a2 <- ggplot(chr1_density, 
                aes(x = BIN_START / 1e6, y = SNP_COUNT, color = species)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = SP_COLORS) +
    facet_wrap(~ species, ncol = 4) +
    labs(
      title    = "SNP density along chromosome 1",
      subtitle = "100 kb windows",
      x        = "Position (Mb)",
      y        = "SNP count per 100 kb"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  p9b2 <- ggplot(chr1_depth,
                aes(x = mid / 1e6, y = depth, color = species)) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    scale_color_manual(values = SP_COLORS) +
    facet_wrap(~ species, ncol = 4) +
    labs(
      title    = "Sequencing depth along chromosome 1",
      subtitle = "100 kb windows (mosdepth)",
      x        = "Position (Mb)",
      y        = "Depth"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  p92 <- p9a2 / p9b2 +
    plot_annotation(
      title    = "Chromosome 1 — SNP density and depth across species",
      subtitle = "Relevant for the large inversion described in the manuscript"
    )
  
  ggsave(file.path(OUTDIR, "09b_chr1_snpdensity_depth.pdf"),
         p92, width = 16, height = 16)
}

# =============================================================================
# PLOT 9c: Chromosome 03 — SNP density + depth across all species for comparison
# =============================================================================

message("[INFO] Plot 9c: Chromosome 3 SNP density and depth...")

# SNP density on chr3
chr3_density <- map_dfr(SPECIES, function(sp) {
  scf <- CHR03[sp]
  density_all %>%
    filter(species == sp, CHROM == scf)
})

# Depth on chr1
chr3_depth <- map_dfr(SPECIES, function(sp) {
  scf <- CHR03[sp]
  depth_wins %>%
    filter(species == sp, scaffold == scf)
})

if (nrow(chr3_density) > 0) {
  
  p9a3 <- ggplot(chr3_density, 
                 aes(x = BIN_START / 1e6, y = SNP_COUNT, color = species)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = SP_COLORS) +
    facet_wrap(~ species, ncol = 4) +
    labs(
      title    = "SNP density along chromosome 3",
      subtitle = "100 kb windows",
      x        = "Position (Mb)",
      y        = "SNP count per 100 kb"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  p9b3 <- ggplot(chr3_depth,
                 aes(x = mid / 1e6, y = depth, color = species)) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    scale_color_manual(values = SP_COLORS) +
    facet_wrap(~ species, ncol = 4) +
    labs(
      title    = "Sequencing depth along chromosome 3",
      subtitle = "100 kb windows (mosdepth)",
      x        = "Position (Mb)",
      y        = "Depth"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  p93 <- p9a3 / p9b3 +
    plot_annotation(
      title    = "Chromosome 3 — SNP density and depth across species",
      subtitle = "Relevant for the large inversion described in the manuscript"
    )
  
  ggsave(file.path(OUTDIR, "09c_chr3_snpdensity_depth.pdf"),
         p93, width = 16, height = 16)
}

# =============================================================================
# PLOT 10: Summary table — key stats per species
# =============================================================================

message("[INFO] Generating summary table...")

# Total SNPs per species
n_snps <- bcftools_stats %>%
  filter(grepl("number of SNPs", key)) %>%
  select(species, n_snps = value)

# Mean depth per species
mean_depth <- depth_site %>%
  group_by(species) %>%
  summarise(mean_depth = mean(SUM_DEPTH, na.rm = TRUE), .groups = "drop")

# Combine summary
summary_table <- het_summary %>%
  select(species, het_rate, F, N_SITES) %>%
  left_join(n_snps,    by = "species") %>%
  left_join(mean_depth, by = "species") %>%
  left_join(titv %>% select(species, titv = value), by = "species") %>%
  arrange(species)

write.csv(summary_table, 
          file.path(OUTDIR, "summary_stats_all_species.csv"), 
          row.names = FALSE)

message("\n[INFO] All plots saved to: ", OUTDIR)
message("[INFO] Output files:")
message("  01_heterozygosity_per_species.pdf")
message("  02_inbreeding_F_per_species.pdf")
message("  03_snp_density_<SP>.pdf          (one per species)")
message("  04_snp_density_heatmap.pdf")
message("  05_allele_frequency_distribution.pdf")
message("  06_depth_distribution.pdf")
message("  07_depth_<SP>.pdf                (one per species)")
message("  08_titv_ratio.pdf")
message("  09_chr18_snpdensity_depth.pdf")
message("  summary_stats_all_species.csv")
