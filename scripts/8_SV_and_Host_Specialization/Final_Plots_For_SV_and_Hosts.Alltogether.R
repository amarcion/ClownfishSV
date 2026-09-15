library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(stringr)
library(ape)
library(phytools)
library(vegan)
library(ggrepel)
library(colorspace)
library(patchwork)
library(purrr)

source("../6_WGA_and_SV/Final_Plots_For_SV.Functions.R")  # adjust path as needed
setwd(".")


CLA_Minimap_Prefix = "ClownSV.CLA.minimap.AllChrom."
BIA_Minimap_Prefix = "ClownSV.BIARef.minimap.AllChrom."

Ref_chromosome_Link_file = read.table("CLA_vs_BIARef.mapids.txt", h=F, 
                                      col.names= c("CLA_chr", "BIA_chr"))


############
# 1. Settings
############

Reference_Species = "BIA"
Mean_Length = TRUE
# We use the threshold of 1kb for merging

if (Reference_Species == "BIA") {CompleteSpeciesName="A. biaculeatus"}else{
  CompleteSpeciesName="A. clarkii"
}


############
# 2. Get data for SV proportion across clownfishes
###########

if (Mean_Length) {suffix = ".MeanLength.txt"}else { suffix = ".txt"}

sv_types <- c("INV", "DUP", "INVDP", "INVTR", "TRANS", "INS", "DEL")

CLA_mnp_all <- read_merged_sv(prefix = CLA_Minimap_Prefix,
                              label = "CLA", sv_types = sv_types,
                              suffix = suffix, add_meta_fn = add_meta)

BIA_mnp_all <- read_merged_sv(prefix = BIA_Minimap_Prefix,
                              label = "BIA", sv_types = sv_types,
                              suffix = suffix,add_meta_fn = add_meta)

SV_all <- dplyr::bind_rows(CLA_mnp_all, BIA_mnp_all)

############
# 2b. Clean and Check dataset - Overall Proportion
###########

# Check the presence of the correct species (no CLA with CLA as reference,
# No BIA when BIA is used as reference, no outgroups)
# Check Clarkii
species_CLA <- check_right_species(SV_all, "CLA")
length(species_CLA)  # should be 17
# Check BIA
species_BIARef <- check_right_species(SV_all, "BIA")
SV_all <- SV_all %>%
  mutate(Species = if_else(
    Reference == "BIA",
    clean_bia_species(Species),Species)) %>%
  filter(!is.na(Species), Species != "")
species_BIARef <- check_right_species(SV_all, "BIA") 
length(species_BIARef) # OK not

# Map CLA chromosomes to BIA chromosomes
SV_all <- SV_all %>% left_join(
  Ref_chromosome_Link_file,
  by = c("Chrom" = "BIA_chr")) %>%
  mutate(Chrom_original = Chrom,
         Chrom = if_else(
           Reference == "BIA",CLA_chr,Chrom)) %>%
  select(-CLA_chr)

# Standardize chromosome order
chrom_levels <- SV_all %>%
  distinct(Chrom) %>%
  mutate(num = as.numeric(gsub("chr", "", Chrom))) %>%
  arrange(num) %>%
  pull(Chrom)

SV_all <- SV_all %>%
  mutate(Chrom = factor(Chrom, levels = chrom_levels))

# Fix SV type
SV_all <- SV_all %>%
  mutate(
    SVType = case_when(
      str_starts(ID, "DEL") ~ "DEL",
      str_starts(ID, "INS") ~ "INS",
      TRUE ~ SV))

# "Merge" TRANS with INVTR et DUP with INVDP, as in the orignal SyRI output. 
# Merging in the sens of considering them for statistics of count and length

SV_all <- SV_all %>%
  mutate(
    SVType_Original = SVType,
    SVType = case_when(
      SVType %in% c("DUP", "INVDP") ~ "DUP",
      SVType %in% c("TRANS", "INVTR") ~ "TRANS",
      TRUE ~ SVType
    )
  )

# Keep consistent order in SV
sv_order <- c("INV", "DUP", "TRANS", "INS", "DEL")
SV_all$SVType <- factor(SV_all$SVType, levels=sv_order)

# Modify the column name if SV_Mean_Length
if (Mean_Length) {
  SV_all <- SV_all %>%
    dplyr::rename(SVLength = SV_Mean_Length) }

# Clean the Species and get number of species
SV_all <- SV_all %>%
  mutate(Species = str_split(Species, ",")) %>%
  mutate(Species = lapply(Species, function(x) unique(trimws(x)))) %>%
  mutate(Species = sapply(Species, paste, collapse = ","))  %>%
  mutate(n_species = sapply(str_split(Species, ","), length)) %>%
  mutate(Type = ifelse(n_species == 1, "Monomorphic", "Polymorphic"))


# Overall statistics 
SV_summary_global <- SV_all  %>%
  group_by(Reference, SVType) %>%
  summarise(
    Count = n(),
    Total_Length = sum(SVLength, na.rm = TRUE),
    .groups = "drop")

SV_summary_chrom <-  SV_all  %>%
  group_by(Chrom, Reference, SVType) %>% 
  summarise(
    Count = n(),
    Total_Length = sum(SVLength, na.rm = TRUE),
    .groups = "drop"
  )

# get the proportions
SV_summary_global_proportion <- SV_summary_global %>%
  pivot_longer(cols = c(Count, Total_Length),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric, Count = "Count", Total_Length = "Total Length"),
         Group = paste(Reference)) %>%
  group_by(Group, Metric) %>%
  mutate(Proportion = Value / sum(Value))


############
# 3. Host specialisation analyses
############

# ------------------------------------------------------------------
# 3a. Settings: host category per species abbreviation
#     ENT = Entacmaea, RAD = Radianthus, STI = Stichodactyla, GEN = Generalist
#     BIA and CLA are reference genomes; their host category is used to
#     correct species counts per host (see comment below).
# ------------------------------------------------------------------

Species_categories <- c(
  "AKA" = "RAD", "AKY" = "GEN", "ALL" = "GEN", "BIA" = "ENT",
  "CRP" = "GEN", "EPH" = "ENT", "FRE" = "ENT", "LAT" = "GEN",
  "LAZ" = "GEN", "MCC" = "ENT", "OMA" = "ENT", "POL" = "STI",
  "PRC" = "RAD", "PRD" = "RAD", "SAN" = "STI", "SEB" = "GEN",
  "OCE" = "RAD", "CLA" = "GEN"
)

# Colour palette shared across all host plots
host_colours <- c(
  "GEN" = "#000000",   # black
  "ENT" = "#DF0C0C",   # red
  "RAD" = "#F3C337",   # yellow
  "STI" = "#FFFFFF"    # white (use black border on bars)
)

sv_type_colours <-  c(
  INV   = "#F09837",
  DUP   = "#56BBF9",
  TRANS = "#B7D86E",
  INS   = "#D87B7B",
  DEL   = "#7B5EA7"
)

# Number of species per host category (used as threshold for "all specialists share the SV")
nb_species_per_host <- table(Species_categories)

# ------------------------------------------------------------------
# 3b. Filter dataset
#     We remove the reference species itself (BIA or CLA), since by
#     definition it can never carry an SV relative to itself.
#     We work on whatever Reference is currently active (Reference_Species).
# ------------------------------------------------------------------

SV_filt <- SV_all %>%
  filter(Reference == Reference_Species)

# Explode the species column so each row = one species carrying the SV
SV_filt_long <- SV_filt %>%
  mutate(Species_list = str_split(Species, ",")) %>%
  unnest(Species_list) %>%
  mutate(
    Species_list = str_trim(Species_list),
    Host         = Species_categories[Species_list]
  )

# ------------------------------------------------------------------
# 3c. Identify host-specific SVs
#     A SV is "host-specific" when:
#       (i)  all species carrying it belong to exactly one host category, AND
#       (ii) every species of that host category carries it.
#
#     Correction for the reference:
#       When BIA is the reference, a SV present in all 5 ENT species
#       does NOT appear in BIA (reference), but BIA is itself ENT.
#       Therefore the real ENT count to match against is 4 (not 5).
#       We handle this by subtracting 1 from nb_species_per_host for the
#       host of the reference species when filtering.
# ------------------------------------------------------------------

ref_host <- Species_categories[Reference_Species]   # host of the current reference

# Effective count of species per host that can carry the SV
#   (reference species cannot carry an SV relative to itself)
nb_species_effective <- nb_species_per_host
nb_species_effective[ref_host] <- nb_species_effective[ref_host] - 1

# Build SV-level summary: for each SV, which hosts carry it?
SV_host_summary <- SV_filt_long %>%
  group_by(ID, SVType, Chrom, Position, SVLength, Chrom_original) %>%
  summarise(
    species   = list(unique(Species_list)),
    hosts     = list(unique(Host)),
    n_species = n_distinct(Species_list),
    # host_type is set only when ALL carrying species share exactly one host
    host_type = ifelse(length(unique(Host)) == 1, unique(Host), NA_character_),
    .groups   = "drop"
  )

# Keep only SVs present in ALL species of a single host category
SV_host_specific <- SV_host_summary %>%
  filter(
    !is.na(host_type),
    n_species == nb_species_effective[host_type]
  )

# Get original position on chrom and  position for 

SV_host_summary_export <- SV_host_specific %>%
  mutate(
    species = sapply(species, paste, collapse = ","),
    hosts   = sapply(hosts, paste, collapse = ",")
  )

table_name = paste0("8_Final_Plots/", Reference_Species, ".HostSpecificSV.csv")
write.table(SV_host_summary_export[SV_host_summary_export$host_type != "ENT",],
            table_name, col.names = TRUE, row.names = F,
            quote = F, sep="\t")


# ------------------------------------------------------------------
# 3d. Count host-specific SVs and plot
# ------------------------------------------------------------------

SV_counts_host <- SV_host_specific %>%
  count(host_type, name = "SV_count")

SV_counts_host_bytype <- SV_host_specific %>%
  count(SVType, host_type, name = "SV_count")

# Here just a correction as no actually ENT
# This si because, no present with CLA, and not present with BIA because if we
# foud the 4 other species, still no present in BIA! 

SV_counts_host <- SV_counts_host %>% filter_out(host_type == "ENT")
SV_counts_host_bytype <- SV_counts_host_bytype %>% filter_out(host_type == "ENT")

plot_SV_by_HostType <- ggplot(
  SV_counts_host,
  aes(x = host_type, y = SV_count, fill = host_type)
) +
  geom_col(color = "black") +
  scale_fill_manual(values = host_colours) +
  theme_minimal() +
  labs(
    title = paste0("SV type composition per host"),
    subtitle = paste0("reference: ", Reference_Species),
    x     = "Host type",
    y     = "Number of SVs",
    fill  = "Host"
  ) +
  theme(legend.position = "bottom")

 
plot_SV_by_HostType_BySVType <- ggplot(
  SV_counts_host_bytype,
  aes(x = host_type, y = SV_count, fill = SVType)
) +
  geom_col(position = "stack", color = "black") +
  scale_fill_manual(values = sv_type_colours) +
  theme_minimal() +
  labs(
    title = paste0("SV type composition per host"),
    subtitle = paste0("reference: ", Reference_Species),
    x     = "Host",
    y     = "SV count",
    fill  = "SV Type"
  )

# ------------------------------------------------------------------
# 3e–3g. Random expectation — STI, RAD, ENT
#
# For each host category, data preparation runs unconditionally so the
# filtered tables are always available. Each block then checks whether
# the FOCAL group contains at least one SV before building plots or
# writing output files. If the focal group is empty the block is
# skipped entirely and a message is printed to the console.
# The shared y-axis maximum is computed only from the active blocks.
# ------------------------------------------------------------------

# ---- STI: focal pair POL-SAN ----

df_pairs <- SV_filt %>%
  mutate(Species_list = str_split(Species, ",")) %>%
  filter(map_int(Species_list, length) == 2) %>%
  mutate(
    sp1  = map_chr(Species_list, 1),
    sp2  = map_chr(Species_list, 2),
    pair = map2_chr(sp1, sp2, ~ paste(sort(c(.x, .y)), collapse = "-"))
  )

pair_counts        <- df_pairs %>% count(pair,         name = "count")
pair_counts_bytype <- df_pairs %>% count(SVType, pair, name = "count")

sti_target_pairs <- c(
  "POL-SAN",   # STI–STI  (focal)
  "AKA-POL",   # RAD–STI  (control)
  "AKA-SEB",   # RAD–GEN  (control)
  "SAN-SEB"    # STI–GEN  (control)
)

pair_counts_sti        <- pair_counts        %>%
  filter(pair %in% sti_target_pairs) %>%
  mutate(pair = factor(pair, levels = sti_target_pairs))
pair_counts_sti_bytype <- pair_counts_bytype %>%
  filter(pair %in% sti_target_pairs) %>%
  mutate(pair = factor(pair, levels = sti_target_pairs))

run_STI <- nrow(pair_counts_sti %>% filter(pair == "POL-SAN")) > 0 &&
  pair_counts_sti %>% filter(pair == "POL-SAN") %>% pull(count) > 0

if (!run_STI) message("STI: no SVs found in focal pair POL-SAN — skipping random expectation.")

# ---- RAD: focal cluster AKA-OCE-PRC-PRD ----

df_fourSpecies <- SV_filt %>%
  mutate(Species_list = str_split(Species, ",")) %>%
  filter(map_int(Species_list, length) == 4) %>%
  mutate(cluster = map_chr(Species_list, ~ paste(sort(.x), collapse = "-")))

fourspecies_counts        <- df_fourSpecies %>% count(cluster,         name = "count")
fourspecies_counts_bytype <- df_fourSpecies %>% count(SVType, cluster, name = "count")

rad_target_clusters <- c(
  "AKA-OCE-PRC-PRD",   # RAD–RAD–RAD–RAD  (focal)
  "OCE-PRC-PRD-SAN"    # RAD–RAD–RAD–STI  (control)
)

fourspecies_counts_rad        <- fourspecies_counts %>%
  filter(cluster %in% rad_target_clusters) %>%
  mutate(cluster = factor(cluster, levels = rad_target_clusters))
fourspecies_counts_rad_bytype <- fourspecies_counts_bytype %>%
  filter(cluster %in% rad_target_clusters) %>%
  mutate(cluster = factor(cluster, levels = rad_target_clusters))

run_RAD <- nrow(fourspecies_counts_rad %>% filter(cluster == "AKA-OCE-PRC-PRD")) > 0 &&
  fourspecies_counts_rad %>% filter(cluster == "AKA-OCE-PRC-PRD") %>% pull(count) > 0

if (!run_RAD) message("RAD: no SVs found in focal cluster AKA-OCE-PRC-PRD — skipping random expectation.")

# ---- ENT: focal cluster built dynamically ----

ent_n <- nb_species_effective["ENT"]   # 4 with BIA ref, 5 with CLA ref

df_entSpecies <- SV_filt %>%
  mutate(Species_list = str_split(Species, ",")) %>%
  filter(map_int(Species_list, length) == ent_n) %>%
  mutate(cluster = map_chr(Species_list, ~ paste(sort(.x), collapse = "-")))

entSpecies_counts        <- df_entSpecies %>% count(cluster,         name = "count")
entSpecies_counts_bytype <- df_entSpecies %>% count(SVType, cluster, name = "count")

ent_focal_species   <- sort(names(Species_categories[Species_categories == "ENT" &
                                                       names(Species_categories) != Reference_Species]))
ent_control_species <- sort(c(setdiff(ent_focal_species, "MCC"), "AKY"))
ent_focal_cluster   <- paste(ent_focal_species,   collapse = "-")
ent_control_cluster <- paste(ent_control_species, collapse = "-")
ent_target_clusters <- c(ent_focal_cluster, ent_control_cluster)

entSpecies_counts_ent <- entSpecies_counts %>%
  filter(cluster %in% ent_target_clusters) %>%
  mutate(cluster = factor(cluster, levels = ent_target_clusters))
entSpecies_counts_ent_bytype <- entSpecies_counts_bytype %>%
  filter(cluster %in% ent_target_clusters) %>%
  mutate(cluster = factor(cluster, levels = ent_target_clusters))

run_ENT <- nrow(entSpecies_counts_ent %>% filter(cluster == ent_focal_cluster)) > 0 &&
  entSpecies_counts_ent %>% filter(cluster == ent_focal_cluster) %>% pull(count) > 0
# Actually correct as never ENT: 
run_ENT <- FALSE
if (!run_ENT) message("ENT: no SVs found in focal cluster ", ent_focal_cluster, " — skipping random expectation.")

# ------------------------------------------------------------------
# Shared y-axis: maximum across all ACTIVE blocks only
# ------------------------------------------------------------------

ymax_inputs <- c()
if (run_STI) ymax_inputs <- c(ymax_inputs,
                              pair_counts_sti$count,
                              pair_counts_sti_bytype %>% group_by(pair) %>% summarise(t = sum(count)) %>% pull(t))
if (run_RAD) ymax_inputs <- c(ymax_inputs,
                              fourspecies_counts_rad$count,
                              fourspecies_counts_rad_bytype %>% group_by(cluster) %>% summarise(t = sum(count)) %>% pull(t))
if (run_ENT) ymax_inputs <- c(ymax_inputs,
                              entSpecies_counts_ent$count,
                              entSpecies_counts_ent_bytype %>% group_by(cluster) %>% summarise(t = sum(count)) %>% pull(t))

expectation_y_max <- if (length(ymax_inputs) > 0) max(ymax_inputs, na.rm = TRUE) else 1







# ------------------------------------------------------------------
# Plots and output — one guarded block per host category
# ------------------------------------------------------------------

############
# 4. Save host-specialisation plots
############

# ---- STI plots ----
if (run_STI) {
  
  plot_hostExpectation_STI <- ggplot(
    pair_counts_sti,
    aes(x = pair, y = count, fill = pair)
  ) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_manual(values = c(
      "POL-SAN" = "#FFFFFF", "AKA-POL" = "gray80",
      "AKA-SEB" = "gray80",  "SAN-SEB" = "gray80"
    )) +
    coord_cartesian(ylim = c(0, expectation_y_max)) +
    theme_minimal() +
    labs(
      title = paste0("SVs shared by STI vs control"),
      subtitle = paste0("reference: ", Reference_Species),
      x = "Species pair", y = "Number of SVs", fill = "Pair"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  plot_hostExpectation_STI_bySVType <- ggplot(
    pair_counts_sti_bytype,
    aes(x = pair, y = count, fill = SVType)
  ) +
    geom_col(position = "stack", color = "black", width = 0.8) +
    scale_fill_manual(values = sv_type_colours) +
    coord_cartesian(ylim = c(0, expectation_y_max)) +
    theme_minimal() +
    labs(
      title = paste0("SVs shared by STI vs control"),
      subtitle = paste0("reference: ", Reference_Species),
      x = "Species pair", y = "SV count", fill = "SV Type"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
}

# ---- RAD plots ----
if (run_RAD) {
  
  plot_hostExpectation_RAD <- ggplot(
    fourspecies_counts_rad,
    aes(x = cluster, y = count, fill = cluster)
  ) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_manual(values = c(
      "AKA-OCE-PRC-PRD" = "#F3C337",   # RAD colour (focal)
      "OCE-PRC-PRD-SAN" = "gray80"     # control
    )) +
    coord_cartesian(ylim = c(0, expectation_y_max)) +
    theme_minimal() +
    labs(
      title = paste0("SVs shared by RAD vs control"),
      subtitle = paste0("reference: ", Reference_Species),
      x = "Species group", y = "Number of SVs", fill = "Group"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  plot_hostExpectation_RAD_bySVType <- ggplot(
    fourspecies_counts_rad_bytype,
    aes(x = cluster, y = count, fill = SVType)
  ) +
    geom_col(position = "stack", color = "black", width = 0.8) +
    scale_fill_manual(values = sv_type_colours) +
    coord_cartesian(ylim = c(0, expectation_y_max)) +
    theme_minimal() +
    labs(
      title = paste0("SVs shared by RAD vs control"),
      subtitle = paste0("reference: ", Reference_Species),
      x = "Species group", y = "SV count", fill = "SV Type"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ---- ENT plots ----
if (run_ENT) {
  
  plot_hostExpectation_ENT <- ggplot(
    entSpecies_counts_ent,
    aes(x = cluster, y = count, fill = cluster)
  ) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_manual(values = setNames(
      c("#DF0C0C", "gray80"),
      c(ent_focal_cluster, ent_control_cluster)
    )) +
    coord_cartesian(ylim = c(0, expectation_y_max)) +
    theme_minimal() +
    labs(
      title = paste0("SVs shared by ENT vs control"),
      subtitle = paste0("reference: ", Reference_Species),
      x = "Species group", y = "Number of SVs", fill = "Group"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  plot_hostExpectation_ENT_bySVType <- ggplot(
    entSpecies_counts_ent_bytype,
    aes(x = cluster, y = count, fill = SVType)
  ) +
    geom_col(position = "stack", color = "black", width = 0.8) +
    scale_fill_manual(values = sv_type_colours) +
    coord_cartesian(ylim = c(0, expectation_y_max)) +
    theme_minimal() +
    labs(
      title = paste0("SVs shared by ENT vs control"),
      subtitle = paste0("reference: ", Reference_Species),
      x = "Species group", y = "SV count", fill = "SV Type"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  
}

# Save the plots

out_prefix <- paste0("8_Final_Plots/HostSV_", Reference_Species, "_")
SV_by_Host_name_op1 <- paste0(out_prefix,"SuppFig.Option1.pdf" )

pdf(SV_by_Host_name_op1, height = 5, width = 12 )

if (run_ENT) {   
  
  
  plot_SV_by_HostType + theme(legend.position = "none", 
                              axis.title.x = element_blank())  +
    plot_hostExpectation_ENT + theme(legend.position = "none", 
                                     axis.title.x = element_blank())  +
    plot_hostExpectation_RAD + theme(legend.position = "none", 
                                     axis.title.x = element_blank())  +
    plot_hostExpectation_STI + theme(legend.position = "none", 
                                     axis.title.x = element_blank()) +
    plot_layout(design ="1345")
  
}else { 
  
  
  plot_SV_by_HostType + theme(legend.position = "none", 
                              axis.title.x = element_blank())  +
    plot_hostExpectation_RAD + theme(legend.position = "none", 
                                     axis.title.x = element_blank())  +
    plot_hostExpectation_STI + theme(legend.position = "none", 
                                     axis.title.x = element_blank()) 
  
}
dev.off()




# Option 2

out_prefix <- paste0("8_Final_Plots/HostSV_", Reference_Species, "_")
SV_by_Host_name_op2 <- paste0(out_prefix,"SuppFig.Option2.pdf" )

pdf(SV_by_Host_name_op2, height = 5, width = 12 )

if (run_ENT) {
  plot_SV_by_HostType_BySVType + theme(legend.position = "none", 
                                       axis.title.x = element_blank())   +
    plot_hostExpectation_ENT_bySVType + theme(legend.position = "none", 
                                              axis.title.x = element_blank())  +
    plot_hostExpectation_RAD_bySVType + theme(legend.position = "none", 
                                              axis.title.x = element_blank())  +
    plot_hostExpectation_STI_bySVType + theme(axis.title.x = element_blank()) +
    plot_layout(design ="1345")
  
  
}else { 
  
  plot_SV_by_HostType_BySVType + theme(legend.position = "none", 
                                       axis.title.x = element_blank())   +
    plot_hostExpectation_RAD_bySVType + theme(legend.position = "none", 
                                              axis.title.x = element_blank())  +
    plot_hostExpectation_STI_bySVType + theme(axis.title.x = element_blank()) 
  
}

dev.off()



