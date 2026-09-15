library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

setwd(".")

dataSV <- read.table("SV_Statstics_DifferentSpecies_DifferentAligner.txt", h=T)



########################
## SV Length for all species
## For the different aligner and different references
########################

bp_labels <- function(x) {
  ifelse(x >= 1e6, paste0(x/1e6, " Mb"),
         ifelse(x >= 1e3, paste0(x/1e3, " kb"),
                paste0(x, " bp")))}

#####
####### With all species
#####

data_long <- dataSV %>%
  pivot_longer(
    cols = c(Inversion_length, Transloc_length, Duplication_length),
    names_to = "SV_type",
    values_to = "Length")

pdf(file="StructuralVariantsLenght_BySpecies_AllSettings.AllSpecies.pdf", 
    width = 12, height = 10)

ggplot(data_long, aes(x = Species, y = Length, fill = SV_type)) +
  geom_bar(stat = "identity") +
  facet_grid(Aligner ~ Reference ) +   #  key line
  labs(title = "Structural Variant Lengths by Species",
       y = "SV Length",
       x = "",
       fill = "SV Type") +
  scale_fill_manual(values = c(
    "Inversion_length" = "#F09837",
    "Transloc_length" = "#B7D86E",
    "Duplication_length" = "#56BBF9"
  )) +
  scale_y_continuous(labels = bp_labels) +
  theme_minimal() +
  coord_flip()
dev.off()

#####
####### Only Clownfishes
#####

data_long_filtered <- data_long %>%
  filter(!Species %in% c("PRCRef", "DTR"))

pdf(file="StructuralVariantsLenght_BySpecies_AllSettings.pdf", 
    width = 12, height = 10)

ggplot(data_long_filtered, aes(x = Species, y = Length, fill = SV_type)) +
  geom_bar(stat = "identity") +
  facet_grid(Reference ~ Aligner) +
  labs(title = "Structural Variant Lengths by Species",
       y = "SV Length",
       x = "",
       fill = "SV Type") +
  scale_fill_manual(values = c(
    "Inversion_length" = "#F09837",
    "Transloc_length" = "#B7D86E",
    "Duplication_length" = "#56BBF9"
  )) +
  scale_y_continuous(labels = bp_labels) +
  theme_minimal() +
  coord_flip()

dev.off()

########################
## Proportions of Alignment
## For the different aligner and different references
########################

# Remove DTR
dataSV <- dataSV %>%
  filter(!Species %in% c("DTR"))

# Remove also  BIA and CLA, depending on the reference
dataSV <- dataSV %>%
  filter(!(Reference == "BIARef" & Species == "BIA") &
           !(Reference == "CLA" & Species == "CLA"))


# Get proprotions
dataSV <- dataSV %>%
  mutate(
    SV_length = Inversion_length + Transloc_length + Duplication_length + long_DEL_length + long_INS_length,
    
    # 1️⃣ analyzed genome proportion
    prop_analyzed = Analzed_Genome_length / Original_Assembly_length,
    
    # 2️⃣ genome composition
    prop_synteny = (Synteny_length - long_INS_length - long_DEL_length) / Analzed_Genome_length,
    prop_not_aligned = NotAligned_length / Analzed_Genome_length,
    prop_sv = SV_length / Analzed_Genome_length,
    
    # 3️⃣ SV composition
    prop_inv = Inversion_length / SV_length,
    prop_transloc = Transloc_length / SV_length,
    prop_dup = Duplication_length / SV_length,
    prop_ins = long_INS_length / SV_length,
    prop_del = long_DEL_length / SV_length
  )

#####
# Analyzed proportion
#####

data_p1 <- dataSV %>%
  mutate(
    analyzed = Analzed_Genome_length / Original_Assembly_length,
    not_analyzed = 1 - analyzed
  ) %>%
  select(Species, Reference, Aligner, analyzed, not_analyzed) %>%
  pivot_longer(cols = c(analyzed, not_analyzed),
               names_to = "category",
               values_to = "value")
data_p1$category <- factor(data_p1$category,
                           levels = c("not_analyzed", "analyzed"))

pdf(file="1_Analyzed_Proportion_Genome.pdf", 
    width = 10, height = 10)
ggplot(data_p1, aes(x = Species, y = value, fill = category)) +
  geom_bar(stat = "identity") +
  facet_grid(Aligner ~ Reference ) +
  labs(title = "Analyzed vs Not Analyzed Genome",
       y = "Proportion",
       x = "",
       fill = "") +
  scale_fill_manual(values = c(
    "analyzed" = "#878787",
    "not_analyzed" = "#EBEBEB"
  )) +
  theme_minimal() +
  coord_flip() 
  #+ scale_y_continuous(labels = scales::percent)
dev.off()

#######
# Proportion of Aligned and not aligned
######

data_comp <- dataSV %>%
  mutate(
    SV_length = Inversion_length + Transloc_length + Duplication_length + long_DEL_length + long_INS_length,
    prop_synteny = (Synteny_length - long_INS_length - long_DEL_length) / Analzed_Genome_length,
    prop_sv = SV_length / Analzed_Genome_length,
    prop_not_aligned = NotAligned_length / Analzed_Genome_length
  ) %>%
  select(Species, Reference, Aligner,
         prop_synteny, prop_sv, prop_not_aligned) %>%
  pivot_longer(cols = starts_with("prop_"),
               names_to = "category",
               values_to = "value")

pdf(file="2_Proportion_Synteny_SV_NotAligned.pdf", 
    width = 10, height = 10)
ggplot(data_comp, aes(x = Species, y = value, fill = category)) +
  geom_bar(stat = "identity") +
  facet_grid(Aligner ~ Reference ) +
  labs(title = "Genome Composition",
       y = "Proportion",
       x = "",
       fill = "") +
  scale_fill_manual(values = c(
    "prop_synteny" = "#D87B7B",
    "prop_sv" = "#6AA3A7",
    "prop_not_aligned" = "#EBEBEB"
  )) +
  theme_minimal() +
  coord_flip()
dev.off()

#######
# Proportion of SV
######

data_sv <- dataSV %>%
  mutate(SV_length = Inversion_length + Transloc_length + Duplication_length + long_DEL_length + long_INS_length) %>%
  mutate(
    prop_inv = Inversion_length / SV_length,
    prop_transloc = Transloc_length / SV_length,
    prop_dup = Duplication_length / SV_length,
    prop_ins = long_INS_length  / SV_length,
    prop_del = long_DEL_length  / SV_length,
  ) %>%
  select(Species, Reference, Aligner,
         prop_inv, prop_transloc, prop_dup, prop_ins, prop_del) %>%
  pivot_longer(cols = starts_with("prop_"),
               names_to = "category",
               values_to = "value")

data_sv$category <- factor(data_sv$category,
                           levels = c("prop_inv", "prop_dup", "prop_transloc",
                                      "prop_del", "prop_ins"))

pdf(file="3_Proportion_SV.pdf",width = 10, height = 10)
ggplot(data_sv, aes(x = Species, y = value, fill = category)) +
  geom_bar(stat = "identity") +
  facet_grid(Aligner ~ Reference ) +
  labs(title = "SV Composition",
       y = "Proportion",
       x = "",
       fill = "SV Type") +
  scale_fill_manual(values = c(
    "prop_inv" = "#F09837",
    "prop_transloc" = "#B7D86E",
    "prop_dup" = "#56BBF9",
    "prop_ins" = "#949191",
    "prop_del" = "#CFCACA"
  )) +
  theme_minimal() +
  coord_flip()

dev.off()

#######
# All information together
######

data_full <- dataSV %>%
  mutate(
    syn = (Synteny_length - long_INS_length - long_DEL_length) / Analzed_Genome_length,
    inv = Inversion_length / Analzed_Genome_length,
    trans = Transloc_length / Analzed_Genome_length,
    dup = Duplication_length / Analzed_Genome_length,
    ins = long_INS_length / Analzed_Genome_length,
    del = long_DEL_length / Analzed_Genome_length,
    not_aligned = NotAligned_length / Analzed_Genome_length
  ) %>%
  select(Species, Reference, Aligner,
         syn, inv, trans, dup, not_aligned, ins, del) %>%
  pivot_longer(cols = -c(Species, Reference, Aligner),
               names_to = "category",
               values_to = "value")

data_full$category <- factor(data_full$category,
                           levels = c("not_aligned", "inv", "dup", "trans",
                                      "del", "ins", "syn"))

pdf(file="4_Proportion_AllTogether.pdf", 
    width = 10, height = 10)

ggplot(data_full, aes(x = Species, y = value, fill = category)) +
  geom_bar(stat = "identity") +
  facet_grid(Aligner ~ Reference ) +
  labs(title = "Full Genome Breakdown",
       y = "Proportion",
       x = "",
       fill = "") +
  scale_fill_manual(values = c(
    "not_aligned" = "#EBEBEB",
    "syn" = "#D87B7B",
    "inv" = "#F09837",
    "dup" = "#56BBF9",
    "trans" = "#B7D86E",
    "ins" = "#949191",
    "del" = "#CFCACA"
    
  )) +
  theme_minimal() +
  coord_flip()

dev.off()

#######
# SV Proportion over Aligned Genome
######

data_sv_grouped <- dataSV %>%
  filter(!(Reference == "BIARef" & Species == "BIA"),
         !(Reference == "CLA" & Species == "CLA")) %>%
  select(Species, Reference, Aligner, Percent_SV_over_Synteny) %>%
  mutate(
    Ref_Aligner = paste(Reference, Aligner, sep = "-"),
    Species = factor(Species, levels = unique(Species)),
    Ref_Aligner = factor(Ref_Aligner, levels = c("CLA-Mummer", "CLA-Minimap",
                                                 "BIARef-Mummer", "BIARef-Minimap"))
  )

# Ensure all species × Reference × Aligner exist
all_combinations <- expand_grid(
  Species = unique(data_sv_grouped$Species),
  Ref_Aligner = c("CLA-Mummer", "CLA-Minimap", "BIARef-Mummer", "BIARef-Minimap")
)

species_order <- c("AKA","AKY","ALL","BIA","CLA","CRP","EPH","FRE","LAT","LAZ","MCC", "OCE",
                   "OMA", "POL", "PRC", "PRCRef", "PRD", "SAN", "SEB")

data_sv_complete <- all_combinations %>%
  left_join(data_sv_grouped, by = c("Species", "Ref_Aligner")) %>%
  # Replace missing values with 0 (or leave NA if you prefer gaps)
  mutate(Percent_SV_over_Synteny = replace_na(Percent_SV_over_Synteny, 0))

data_sv_complete <- data_sv_complete %>%
  mutate(Species = factor(Species, levels = species_order))

pdf("5_Proportion_SV_Over_Aligned.pdf", width = 8, height = 12)
ggplot(data_sv_complete, aes(x = Species, y = Percent_SV_over_Synteny, fill = Ref_Aligner)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Percent SV (over Aligned Genome)",
       y = "Percent SV",
       x = "Species",
       fill = "Reference - Aligner") +
  scale_fill_manual(values = c(
    "CLA-Mummer" = "#F27979",
    "CLA-Minimap" = "#8A1717",
    "BIARef-Mummer" = "#6C9DE6",
    "BIARef-Minimap" = "#234B87"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  coord_flip()

dev.off()


#######
# SV Proportion over Analized (Aligned + Not Aligned) Genome
######

data_sv_grouped <- dataSV %>%
  filter(!(Reference == "BIARef" & Species == "BIA"),
         !(Reference == "CLA" & Species == "CLA")) %>%
  select(Species, Reference, Aligner, SV_Percentage_overAnalyzed) %>%
  mutate(
    Ref_Aligner = paste(Reference, Aligner, sep = "-"),
    Species = factor(Species, levels = unique(Species)),
    Ref_Aligner = factor(Ref_Aligner, levels = c("CLA-Mummer", "CLA-Minimap",
                                                 "BIARef-Mummer", "BIARef-Minimap"))
  )

# Ensure all species × Reference × Aligner exist
all_combinations <- expand_grid(
  Species = unique(data_sv_grouped$Species),
  Ref_Aligner = c("CLA-Mummer", "CLA-Minimap", "BIARef-Mummer", "BIARef-Minimap")
)

species_order <- c("AKA","AKY","ALL","BIA","CLA","CRP","EPH","FRE","LAT","LAZ","MCC", "OCE",
                   "OMA", "POL", "PRC", "PRCRef", "PRD", "SAN", "SEB")

data_sv_complete <- all_combinations %>%
  left_join(data_sv_grouped, by = c("Species", "Ref_Aligner")) %>%
  # Replace missing values with 0 (or leave NA if you prefer gaps)
  mutate(SV_Percentage_overAnalyzed = replace_na(SV_Percentage_overAnalyzed, 0))

data_sv_complete <- data_sv_complete %>%
  mutate(Species = factor(Species, levels = species_order))

pdf("6_Proportion_SV_Over_AnalyzedGenome.pdf", width = 8, height = 12)
ggplot(data_sv_complete, aes(x = Species, y = SV_Percentage_overAnalyzed, fill = Ref_Aligner)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Percent SV (over Analyzed Genome)",
       y = "Percent SV",
       x = "Species",
       fill = "Reference - Aligner") +
  scale_fill_manual(values = c(
    "CLA-Mummer" = "#F27979",
    "CLA-Minimap" = "#8A1717",
    "BIARef-Mummer" = "#6C9DE6",
    "BIARef-Minimap" = "#234B87"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  coord_flip()

dev.off()


#######
# SNPs Proportion over Synteny regions
#######

# Remove DTR
dataSV <- dataSV %>%
  filter(!Species %in% c("DTR"))

# Remove also  BIA and CLA, depending on the reference
dataSV <- dataSV %>%
  filter(!(Reference == "BIARef" & Species == "BIA") &
           !(Reference == "CLA" & Species == "CLA"))

data_sv_grouped <- dataSV %>%
  filter(!(Reference == "BIARef" & Species == "BIA"),
         !(Reference == "CLA" & Species == "CLA")) %>%
  select(Species, Reference, Aligner, SNP_divergence) %>%
  mutate(
    Ref_Aligner = paste(Reference, Aligner, sep = "-"),
    Species = factor(Species, levels = unique(Species)),
    Ref_Aligner = factor(Ref_Aligner, levels = c("CLA-Mummer", "CLA-Minimap",
                                                 "BIARef-Mummer", "BIARef-Minimap"))
  )

# Ensure all species × Reference × Aligner exist
all_combinations <- expand_grid(
  Species = unique(data_sv_grouped$Species),
  Ref_Aligner = c("CLA-Mummer", "CLA-Minimap", "BIARef-Mummer", "BIARef-Minimap")
)

species_order <- c("AKA","AKY","ALL","BIA","CLA","CRP","EPH","FRE","LAT","LAZ","MCC", "OCE",
                   "OMA", "POL", "PRC", "PRCRef", "PRD", "SAN", "SEB")

data_sv_complete <- all_combinations %>%
  left_join(data_sv_grouped, by = c("Species", "Ref_Aligner")) %>%
  # Replace missing values with 0 (or leave NA if you prefer gaps)
  mutate(SNP_divergence = replace_na(SNP_divergence, 0))

data_sv_complete <- data_sv_complete %>%
  mutate(Species = factor(Species, levels = species_order))

pdf("7_SNPs_Proportion_Over_Analyzed.pdf", width = 8, height = 12)
ggplot(data_sv_complete, aes(x = Species, y = SNP_divergence, fill = Ref_Aligner)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Percent SNPs (over Analyzed Genome)",
       y = "Percent SNPs",
       x = "Species",
       fill = "Reference - Aligner") +
  scale_fill_manual(values = c(
    "CLA-Mummer" = "#F27979",
    "CLA-Minimap" = "#8A1717",
    "BIARef-Mummer" = "#6C9DE6",
    "BIARef-Minimap" = "#234B87"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  coord_flip()

dev.off()


########################
## Plot only BIA vs BIARef
########################

dataSV_BIA <- dataSV %>%
  filter(Reference=="BIARef" & Species=="BIA")

BIA_vs_BIARef <- dataSV_BIA %>%
  select(Aligner, Synteny_length, Inversion_length, Transloc_length, 
         Duplication_length, NotAligned_length) %>%
  pivot_longer(
    cols = -Aligner,
    names_to = "Variant",
    values_to = "Length"
  )

BIA_vs_BIARef$category <- factor(BIA_vs_BIARef$Variant,
                           levels = c("Synteny_length", "Inversion_length", "Transloc_length",
                                      "Duplication_length", "NotAligned_length"))

BIARef_vs_BIA <- ggplot(BIA_vs_BIARef[BIA_vs_BIARef$Aligner=="Mummer",], aes(x = Aligner, y = Length, fill = Variant)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion of genome") +
  scale_fill_manual(values = c(
    "Synteny_length" = "#D87B7B",
    "Inversion_length" = "#F09837",
    "Transloc_length" = "#B7D86E",
    "Duplication_length" = "#56BBF9", 
    "NotAligned_length" = "#EBEBEB")) +
  xlab("Aligner") +
  theme_minimal() + 
  coord_flip() +
  labs(title = "BIA Reference vs BIA")

pdf("BIARef_vs_BIA.SV.pdf", width = 10, height =5)
BIARef_vs_BIA
dev.off()  

# ONLY COMPARISON BIARef - CLA and BIARef-BIA Mummer


dataSV_BIA_CLA_only <- dataSV %>%
  filter((Reference=="BIARef" & Species=="CLA") |
           (Reference=="BIARef" & Species=="BIA"))  %>%
  select(Aligner, Species, Synteny_length, Inversion_length, Transloc_length, 
         Duplication_length, NotAligned_length)


dataSV_BIA_CLA_only_long <- dataSV_BIA_CLA_only %>%
  pivot_longer(
    cols = c(Synteny_length, Inversion_length, Transloc_length, 
             Duplication_length, NotAligned_length),
    names_to = "Variant",
    values_to = "Length"
  )  %>%
  mutate(Group = paste(Aligner, Species, sep = "_"))

BIARef_vs_BIA_CLA <- ggplot(dataSV_BIA_CLA_only_long, aes(y = Species, x = Length, fill = Variant)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Aligner, nrow=2) +
  xlab("Cumulative length") +
  theme_minimal() +
  scale_fill_manual(values = c(
    "Synteny_length" = "#D87B7B",
    "Inversion_length" = "#F09837",
    "Transloc_length" = "#B7D86E",
    "Duplication_length" = "#56BBF9", 
    "NotAligned_length" = "#EBEBEB")) +
  labs(title = "BIA Reference with Mummer")


pdf("BIARef_vs_BIA_And_CLA.SV.pdf", width = 10, height =5)
BIARef_vs_BIA_CLA
dev.off()  

##### Only BIA and CLA for BIARef 

dataSV_BIA_CLA <- dataSV %>%
  filter(Reference=="BIARef" & (Species=="BIA" | Species=="CLA") &
           Aligner=="Minimap")

BIA_CLA_vs_BIARef <- dataSV_BIA_CLA %>%
  select(Species, Synteny_length, Inversion_length, Transloc_length, 
         Duplication_length, NotAligned_length) %>%
  pivot_longer(
    cols = -Species,
    names_to = "Variant",
    values_to = "Length"
  )

BIA_CLA_vs_BIARef$category <- factor(BIA_CLA_vs_BIARef$Variant,
                                 levels = rev(c("Synteny_length", "Inversion_length", "Transloc_length",
                                            "Duplication_length", "NotAligned_length")))
BIA_CLA_vs_BIARef$Species <- factor(BIA_CLA_vs_BIARef$Species,
                                    levels=c("CLA", "BIA"))


BIARef_vs_BIA_CLA <- ggplot(BIA_CLA_vs_BIARef, aes(x = Species, y = Length, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion of genome") +
  scale_fill_manual(values = c(
    "Synteny_length" = "#D87B7B",
    "Inversion_length" = "#F09837",
    "Transloc_length" = "#B7D86E",
    "Duplication_length" = "#56BBF9", 
    "NotAligned_length" = "#EBEBEB")) +
  xlab("Species") +
  theme_minimal() + 
  coord_flip() +
  labs(title = "BIA Reference vs BIA")

pdf("BIARef_vs_BIA_CLA.SV.pdf", width = 10, height =5)
BIARef_vs_BIA_CLA
dev.off()  




