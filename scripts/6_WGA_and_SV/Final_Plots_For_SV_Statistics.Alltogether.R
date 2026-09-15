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

source("Final_Plots_For_SV.Functions.R")  # sibling script in scripts/6_WGA_and_SV/ ; adjust path as needed
setwd(".")

# See Summarize_Stats_MergedSV.By_Methods.1kb.R and  Summarize_Stats_MergedSV.By_Methods.R 
# for results with both minimap2 and MUMmer4
# Also, here we only plot figures to be in the paper as main or supplementary


CLA_Minimap_Prefix = "ClownSV.CLA.minimap.AllChrom."
BIA_Minimap_Prefix = "ClownSV.BIARef.minimap.AllChrom."

# CLA_vs_BIARef.mapids.txt available in metadata
Ref_chromosome_Link_file = read.table("CLA_vs_BIARef.mapids.txt", h=F, 
                                      col.names= c("CLA_chr", "BIA_chr"))


############
# 1. Settings
############

Reference_Species = "CLA"
Mean_Length = TRUE
# We use the threshold of 1kb for merging

if (Reference_Species == "BIA") {CompleteSpeciesName="A. biaculeatus"}else{
  CompleteSpeciesName="A. clarkii"
}

############
# 1a. Data for SV by species
############

# Getting data
dataSV <- read.table(
  "SV_Statstics_DifferentSpecies_DifferentAligner.txt",  # Obtained with Extract_stats_from_SyRI.py used on the different reference outputs, and merged
  header = TRUE) %>%
  mutate(Reference = gsub("BIARef", "BIA", Reference))
  
if (Reference_Species=="BIA") {
  dataSV <- dataSV  %>%
    filter(
    !Species %in% c("PRCRef", "DTR"),
    Reference == Reference_Species,
    Aligner == "Minimap",
    Species != "BIA")} else {
    
    dataSV <- dataSV  %>%
        filter(
          !Species %in% c("PRCRef", "DTR"),
          Reference == Reference_Species,
          Aligner == "Minimap",
          Species != "CLA")
}  
      

# Get "easier" data to plot
data_full <- make_full_data(dataSV)
data_sv <- make_sv_data(dataSV)

# Correcte data_full to get to full proportion as not_aligned

data_full <- data_full %>%
  group_by(Species, Reference) %>%
  mutate(
    value = ifelse(
      category == "not_aligned",
      1 - sum(value[category != "not_aligned"]),
      value
    )
  ) %>%
  ungroup()


# Get more complicated data to plot (two plots together)
# Build plotting table
plot_data <- bind_rows(
  make_comp_data(dataSV),
  make_sv_data(dataSV)
)

plot_data <- make_positions(plot_data)

# Labels
label_df <- plot_data %>%
  group_by(Species) %>%
  summarise(
    x = mean(x),
    label = first(Species),
    .groups = "drop"
  ) %>%
  arrange(x)

# Category oder
plot_data$category <- factor(plot_data$category,
  levels = c(
    "not_aligned", "syn", "inv",
    "dup", "trans", "del", "ins",
    "prop_inv", "prop_dup",
    "prop_transloc", "prop_del",
    "prop_ins"))

plot_data <- plot_data %>%
  mutate(
    plot_type = factor(plot_type, levels = c("comp_withsv", "sv"))
  )


############
# 1b. PLOTS for SV by species
############

possible_title <- bquote(
  paste("Genome-wide Comparisons to ",
        italic(.(CompleteSpeciesName)), sep=""))

plot_SV_Species_Separated <- plot_sv_composition_perSpecies(
  plot_data, label_df, title_text = possible_title)

plot_SV_Species_Separated_V2 <- plot_genome_breakdown(data_full, reference = Reference_Species,
                                                      title_text = possible_title)

plot_SV_Species_Separated_V3 <- plot_sv_composition(data_sv, reference = Reference_Species, 
                                  title_text = possible_title)



############
# 2a. Get data for SV proportion across clownfishes
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

# Modify the column name if SV_Mean_Lengt
# Modify SV_Length if SV_Mean_Length
if (Mean_Length) {
  SV_all <- SV_all %>%
    rename(SVLength = SV_Mean_Length) }

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
# 2c. Performe PCoA (MDS)
###########

# Run PCoA
pcoa_SV <- run_pcoa_sv(SV_all, Reference_Species)

# Get color informations
# Information for coloration based on clades
clade_info <- data.frame(
  Species = c("FRE","EPH","CLA", "AKA","SAN","PRD", "OMA","ALL","LAT",
              "SEB","POL", "CRP","MCC","AKY", "LAZ", "PRC","OCE", "BIA"),
  Clade = c("Clade1","Clade1","Clade1", "Clade2","Clade2","Clade2",
    "Clade3","Clade3","Clade3", "Clade4","Clade4", "Clade5","Clade5",
    "Clade5", "Clade6","Clade7","Clade7","Clade8"))

clade_colors <- c("Clade1" = "#F28E2B", "Clade2" = "#EDC948", "Clade3" = "#9C4DCC",
  "Clade4" = "#3283C9", "Clade5" = "#0025DB", "Clade6" = "#09850E",
  "Clade7" = "#97B898", "Clade8" = "#66DB00")

# Generate color shades
# Generate the colors
species_colors <- c()
for (cl in unique(clade_info$Clade)) {
  sp <- clade_info$Species[clade_info$Clade == cl]
  base_col <- clade_colors[cl]

  species_colors <- c(species_colors, generate_clade_colors(sp, cl, base_col))
}

# Generate the gray colors
species_colors_gray <- c()
clades <- unique(clade_info$Clade)

for (i in seq_along(clades)) {
  
  cl <- clades[i]
  sp <- clade_info$Species[clade_info$Clade == cl]
  
  species_colors_gray <- c(
    species_colors_gray,
    generate_gray_clade_colors(sp, i, length(clades))
  )
}


############
# 2d. Stats and plot the data
###########

# Total Count and Total Length depending on Reference
totLength <- sum(SV_summary_global$Total_Length[SV_summary_global$Reference==Reference_Species])
totCount <- sum(SV_summary_global$Count[SV_summary_global$Reference==Reference_Species])
print(paste("Total SV length for", Reference_Species, "as reference: ",totLength))
print(paste("Total SV count for", Reference_Species, "as reference: ",totCount))

Possible_plot_Tiles = paste("SVs across clownfishes: ", totCount, " SVs , ~ ", 
                            round(totLength/1000000, 0), " Mb", sep="")

# Overall proportion in count and length
plot_SV_merged_Proportion <- plot_Merged_SV_proportion(SV_summary_global_proportion, "BIA", 
                                plottitle = Possible_plot_Tiles) 

# Species distribution
plot_Species_Distribution <- plot_SV_species_distribution(SV_all, Reference_Species, 
                                title = "Private and shared SVs across species")

# PCoA overall SV (presence/absence, all SVs)
plot_PCoA <- plot_pcoa_sv(pcoa_SV, Reference_Species, species_colors, 
                          title = "PCoA of Species Based on SVs")

plot_PCoA_gray <- plot_pcoa_sv(pcoa_SV, Reference_Species, species_colors_gray, 
                               title = "PCoA of Species Based on SVs")

############
# 3. ORGANIZE PLOTS For Figure
###########

# Option 1 
layout_op1 <- "
122
134
"
pdf_name_op1 = paste("8_Final_Plots/Figure_SV.", Reference_Species, 
                 ".Stats_Option_1.pdf", sep = "")
pdf(pdf_name_op1, width = 15, height = 8)
plot_SV_Species_Separated + theme(legend.position = "none") +
  plot_SV_merged_Proportion + theme(legend.position = "right") +
  plot_Species_Distribution + plot_PCoA +
  plot_layout(design = layout_op1,
              widths = c(1, 1, 1),
              heights = c(0.7, 1)) +
  plot_annotation(tag_levels = "A")
dev.off()

# Option 2
layout_op2 <- "
12
13
14
"
pdf_name_op2 = paste("8_Final_Plots/Figure_SV.", Reference_Species, 
                 ".Stats_Option_2.pdf", sep = "")
pdf(pdf_name_op2, width = 12, height = 11)
plot_SV_Species_Separated + theme(legend.position = "none") +
  plot_SV_merged_Proportion + theme(legend.position = "right", 
    axis.text.y = element_text(angle = 45)) +
  plot_Species_Distribution + plot_PCoA +
  plot_layout(design = layout_op2,
              widths = c(1.2, 1),
              heights = c(0.7, 0.8, 1.1)) +
  plot_annotation(tag_levels = "A")
dev.off()

theme(
  axis.title.y = element_text(angle = 0)
)

# Option 1b 
pdf_name_op1b = paste("8_Final_Plots/Figure_SV.", Reference_Species, 
                     ".Stats_Option_1b.pdf", sep = "")
pdf(pdf_name_op1b, width = 15, height = 8)
plot_SV_Species_Separated_V2 + theme(legend.position = "none") +
  plot_SV_merged_Proportion + theme(legend.position = "right") +
  plot_Species_Distribution + plot_PCoA +
  plot_layout(design = layout_op1,
              widths = c(1.2, 0.8, 0.8),
              heights = c(0.7, 1)) +
  plot_annotation(tag_levels = "A")
dev.off()

# Option 2b
pdf_name_op2b = paste("8_Final_Plots/Figure_SV.", Reference_Species, 
                     ".Stats_Option_2b.pdf", sep = "")
pdf(pdf_name_op2b, width = 12, height = 11)
plot_SV_Species_Separated_V2 + theme(legend.position = "none") +
  plot_SV_merged_Proportion + theme(legend.position = "right", 
                                    axis.text.y = element_text(angle = 45)) +
  plot_Species_Distribution + plot_PCoA +
  plot_layout(design = layout_op2,
              widths = c(1.2, 1),
              heights = c(0.7, 0.8, 1.1)) +
  plot_annotation(tag_levels = "A")
dev.off()

# Option 1c
pdf_name_op1c = paste("8_Final_Plots/Figure_SV.", Reference_Species, 
                      ".Stats_Option_1c.pdf", sep = "")
pdf(pdf_name_op1c, width = 15, height = 8)
plot_SV_Species_Separated_V3 + theme(legend.position = "none") +
  plot_SV_merged_Proportion + theme(legend.position = "right") +
  plot_Species_Distribution + plot_PCoA +
  plot_layout(design = layout_op1,
              widths = c(1.2, 0.8, 0.8),
              heights = c(0.7, 1)) +
  plot_annotation(tag_levels = "A")
dev.off()

# Option 2c
pdf_name_op2c = paste("8_Final_Plots/Figure_SV.", Reference_Species, 
                      ".Stats_Option_2c.pdf", sep = "")
pdf(pdf_name_op2c, width = 12, height = 11)
plot_SV_Species_Separated_V3 + theme(legend.position = "none") +
  plot_SV_merged_Proportion + theme(legend.position = "right", 
                                    axis.text.y = element_text(angle = 45)) +
  plot_Species_Distribution + plot_PCoA +
  plot_layout(design = layout_op2,
              widths = c(1.2, 1),
              heights = c(0.7, 0.8, 1.1)) +
  plot_annotation(tag_levels = "A")
dev.off()

# Option 3
# With legend, possibly final figure

pdf_name_op3 = paste("8_Final_Plots/Figure_SV.", Reference_Species, 
                      ".Stats_Option_3.pdf", sep = "")
pdf(pdf_name_op3, width = 12, height = 11)
plot_SV_Species_Separated_V2 + sv_legend_scale() + 
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) +
  plot_SV_merged_Proportion + theme(legend.position = "right", 
                                    axis.text.y = element_text(angle = 45)) +
  plot_Species_Distribution + plot_PCoA +
  plot_layout(design = layout_op2,
              widths = c(1.2, 1),
              heights = c(0.7, 0.8, 1.1)) +
  plot_annotation(tag_levels = "A")
dev.off()


# Option 3b
# With legend and gray PCA, possibly final figure
pdf_name_op3b = paste("8_Final_Plots/Figure_SV.", Reference_Species, 
                     ".Stats_Option_3b.pdf", sep = "")
pdf(pdf_name_op3b, width = 12, height = 11)
plot_SV_Species_Separated_V2 + sv_legend_scale() + 
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) +
  plot_SV_merged_Proportion + theme(legend.position = "right", 
                                    axis.text.y = element_text(angle = 45)) +
  plot_Species_Distribution + plot_PCoA_gray +
  plot_layout(design = layout_op2,
              widths = c(1.2, 1),
              heights = c(0.7, 0.8, 1.1)) +
  plot_annotation(tag_levels = "A")
dev.off()



############
# 4. SUPPLEMENTARY STATISTICS and PLOTS - SV LENGTH
###########

############
# 4a SUPPLEMENTARY STATISTICS 
###########

sv_colours <- c(
  INV   = "#F09837",
  DUP   = "#56BBF9",
  TRANS = "#B7D86E",
  INS   = "#D87B7B",
  DEL   = "#7B5EA7"
)

sv_summary_bySVType <- SV_all %>%
  group_by(Reference, SVType) %>%
  summarise(
    Count       = n(),
    TotalLength = sum(SVLength, na.rm = TRUE),
    .groups = "drop"
  )

# Change order of SVType for plotting
SV_all$SVType <- factor(SV_all$SVType, levels=rev(levels(SV_all$SVType)))

# Plot SV Count and Total Length by SVType
p_MergeSV_ByType <- plot_sv_by_type(SV_all, Reference_Species)

# Plot SV Length Distribution by SVType
p_SVLengthDist <- plot_sv_length_distribution(SV_all, Reference_Species)

# Plot SV Number and length by chromosomes
p_SV_per_Chrom <- plot_sv_chrom(SV_all, Reference_Species)

# Arrange plots

layout_supp1 <- "
123
444
555
"

pdf_name_suppOp1 = paste("8_Final_Plots/Figure_SV.SupplementaryStats.", Reference_Species, 
                      ".Stats.Option1.pdf", sep = "")
pdf(pdf_name_suppOp1, width = 12, height = 11)
p_MergeSV_ByType$count + labs(x = NULL) +
  p_MergeSV_ByType$length  + labs(x = NULL) +
  p_SVLengthDist$violin + labs(x = NULL) +
  p_SV_per_Chrom$count + theme(legend.position = "none")+
  p_SV_per_Chrom$length +
  plot_layout(design = layout_supp1, 
              widths = c(1,1,1), 
              heights = c(1.3, 0.8, 0.8)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.subtitle = element_blank(),
    plot.title = element_text(size = 10)
  )
dev.off()


############
# 5c PCoA per SV type
###########

pcoa_bytype <- run_pcoa_sv_bytype(SV_all, Reference_Species)

# Plot with clade colours
plots_color <- plot_pcoa_sv_bytype(pcoa_bytype, Reference_Species, species_colors)
plots_gray  <- plot_pcoa_sv_bytype(pcoa_bytype, Reference_Species, species_colors_gray)

# Combine all into one figure with patchwork

pdf_name_suppPCoA_op1 <- paste("8_Final_Plots/Figure_SV.SupplementaryStats.", Reference_Species, 
                               ".PCoA_ByType.Option1.pdf", sep = "")
pdf(pdf_name_suppPCoA_op1, width = 8, height = 9)
wrap_plots(plots_color, ncol = 2)
dev.off()

pdf_name_suppPCoA_op2 <- paste("8_Final_Plots/Figure_SV.SupplementaryStats.", Reference_Species, 
                               ".PCoA_ByType.Option2.pdf", sep = "")
pdf(pdf_name_suppPCoA_op2, width = 8, height = 9)
wrap_plots(plots_gray, ncol = 2)
dev.off()


############
# 6. STATISTICS ON SV DISTRIBUTION ACROSS THE GENOME
###########

# A) Full SV dataset
# Summarise

sv_summary_all <- SV_all %>% 
  filter(Reference == Reference_Species) %>%
  group_by(SVType, Chrom) %>%
  summarise(
    count = n(),
    total_length = sum(SVLength, na.rm = TRUE),
    .groups = "drop"
  )

# Normalize within SVType
sv_long_all <- sv_summary_all %>%
  group_by(SVType) %>%
  mutate(
    count_frac = count / sum(count),
    length_frac = total_length / sum(total_length)
  ) %>%
  ungroup() %>%
  pivot_longer(
    cols = c(count_frac, length_frac),
    names_to = "metric",
    values_to = "fraction"
  ) %>%
  mutate(metric = recode(metric,
                         count_frac = "Count fraction",
                         length_frac = "Length fraction"))

# Plot
plot_all_sv <- ggplot(sv_long_all, aes(x = Chrom, y = fraction, fill = fraction)) +
  geom_col() +
  facet_grid(metric ~ SVType, scales = "free_y") +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = paste("Polymorphic SVs (≥2 species): normalized distribution | ", Reference_Species ),
    x = "Chromosome",
    y = "Fraction within SVType"
  )

# B) Polymorphic SVs (≥ 2 species)

SV_poly <- SV_all %>% 
  filter(Reference == Reference_Species)  %>% 
  filter(n_species >= 2)

sv_summary_poly <- SV_poly %>%
  group_by(SVType, Chrom) %>%
  summarise(
    count = n(),
    total_length = sum(SVLength, na.rm = TRUE),
    .groups = "drop"
  )

sv_long_poly <- sv_summary_poly %>%
  group_by(SVType) %>%
  mutate(
    count_frac = count / sum(count),
    length_frac = total_length / sum(total_length)
  ) %>%
  ungroup() %>%
  pivot_longer(
    cols = c(count_frac, length_frac),
    names_to = "metric",
    values_to = "fraction"
  ) %>%
  mutate(metric = recode(metric,
                         count_frac = "Count fraction",
                         length_frac = "Length fraction"))

plot_poly_sv <- ggplot(sv_long_poly, aes(x = Chrom, y = fraction, fill = fraction)) +
  geom_col() +
  facet_grid(metric ~ SVType, scales = "free_y") +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = paste("Polymorphic SVs (≥2 species): normalized distribution | ", Reference_Species ),
    x = "Chromosome",
    y = "Fraction within SVType"
  )

# C) PLOT

pdf_name_supp3 = paste("Figure_SV.SupplementaryStats.", Reference_Species, 
                         ".Normalized_ByChrom.Option1.pdf", sep = "")
pdf(pdf_name_supp3, width = 12, height = 11)
plot_all_sv + plot_poly_sv +
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = "A") 
  
dev.off()


##########
# Save also the table
##########
# write.table(SV_all, file="Final_Dataset_SV.Minimap.BIA_and_CLA.csv", 
#             row.names = F, col.names = T, quote = F, sep="\t")



