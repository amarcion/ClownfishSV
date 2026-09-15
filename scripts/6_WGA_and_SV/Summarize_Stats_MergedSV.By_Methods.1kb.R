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

setwd(".")

CLA_Mummer_Prefix = "ClownSV.CLA.AllChrom."
CLA_Minimap_Prefix = "ClownSV.CLA.minimap.AllChrom."
BIA_Mummer_Prefix = "ClownSV.BIARef.Mummer.AllChrom."
BIA_Minimap_Prefix = "ClownSV.BIARef.minimap.AllChrom."

Ref_chromosome_Link_file = read.table("CLA_vs_BIARef.mapids.txt", h=F, 
                                      col.names= c("CLA_chr", "BIA_chr"))


##########################################
# GET THE DATA 
###########################################

Mean_Length = FALSE

if (Mean_Length) {
  suffix = ".MeanLength.txt"
}else {
  suffix = ".txt"
}


# Get the "normal" length"

# CLARKII MUMMER (ORIGINAL)
CLA_mum_INV <- read.table(paste(CLA_Mummer_Prefix,"INV.1kb_trsh.AllSVs",suffix, sep=""), h=T)
CLA_mum_DUP <- read.table(paste(CLA_Mummer_Prefix,"DUP.1kb_trsh.AllSVs", suffix, sep=""), h=T)
CLA_mum_INVDP <- read.table(paste(CLA_Mummer_Prefix,"INVDP.1kb_trsh.AllSVs", suffix , sep=""), h=T)
CLA_mum_INVTR <- read.table(paste(CLA_Mummer_Prefix,"INVTR.1kb_trsh.AllSVs", suffix, sep=""), h=T)
CLA_mum_TRANS <- read.table(paste(CLA_Mummer_Prefix,"TRANS.1kb_trsh.AllSVs", suffix , sep=""), h=T)
CLA_mum_INS <- read.table(paste(CLA_Mummer_Prefix,"INS.1kb_trsh.AllSVs", suffix, sep=""), h=T)
CLA_mum_DEL <- read.table(paste(CLA_Mummer_Prefix,"DEL.1kb_trsh.AllSVs", suffix, sep=""), h=T)

# Original Nb INS and DEL (before merging based on the localisation)
CLA_mum_original_Nb_INS = 64421
CLA_mum_original_Nb_DEL = 87771
CLA_mum_original_sum_INS = 19053203
CLA_mum_original_sum_DEL = 18547690

# CLARKII MINIMAP
CLA_mnp_INV <- read.table(paste(CLA_Minimap_Prefix,"INV.1kb_trsh.AllSVs",suffix, sep=""), h=T)
CLA_mnp_DUP <- read.table(paste(CLA_Minimap_Prefix,"DUP.1kb_trsh.AllSVs",suffix, sep=""), h=T)
CLA_mnp_INVDP <- read.table(paste(CLA_Minimap_Prefix,"INVDP.1kb_trsh.AllSVs",suffix, sep=""), h=T)
CLA_mnp_INVTR <- read.table(paste(CLA_Minimap_Prefix,"INVTR.1kb_trsh.AllSVs",suffix, sep=""), h=T)
CLA_mnp_TRANS <- read.table(paste(CLA_Minimap_Prefix,"TRANS.1kb_trsh.AllSVs",suffix , sep=""), h=T)
CLA_mnp_INS <- read.table(paste(CLA_Minimap_Prefix,"INS.1kb_trsh.AllSVs",suffix, sep=""), h=T)
CLA_mnp_DEL <- read.table(paste(CLA_Minimap_Prefix,"DEL.1kb_trsh.AllSVs", suffix,sep=""), h=T)

# Original Nb INS and DEL (before merging based on the localisation)
CLA_mnp_original_Nb_INS = 132344
CLA_mnp_original_Nb_DEL = 159850
CLA_mnp_original_sum_INS = 34250341
CLA_mnp_original_sum_DEL = 32038807


# BIARef MUMMER 
BIA_mum_INV <- read.table(paste(BIA_Mummer_Prefix,"INV.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mum_DUP <- read.table(paste(BIA_Mummer_Prefix,"DUP.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mum_INVDP <- read.table(paste(BIA_Mummer_Prefix,"INVDP.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mum_INVTR <- read.table(paste(BIA_Mummer_Prefix,"INVTR.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mum_TRANS <- read.table(paste(BIA_Mummer_Prefix,"TRANS.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mum_INS <- read.table(paste(BIA_Mummer_Prefix,"INS.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mum_DEL <- read.table(paste(BIA_Mummer_Prefix,"DEL.1kb_trsh.AllSVs",suffix, sep=""), h=T)

# Original Nb INS and DEL (before merging based on the localisation)
BIA_mum_original_Nb_INS = 21554
BIA_mum_original_Nb_DEL = 27344
BIA_mum_original_sum_INS = 5067895
BIA_mum_original_sum_DEL = 5435675

# BIA MINIMAP
BIA_mnp_INV <- read.table(paste(BIA_Minimap_Prefix,"INV.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mnp_DUP <- read.table(paste(BIA_Minimap_Prefix,"DUP.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mnp_INVDP <- read.table(paste(BIA_Minimap_Prefix,"INVDP.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mnp_INVTR <- read.table(paste(BIA_Minimap_Prefix,"INVTR.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mnp_TRANS <- read.table(paste(BIA_Minimap_Prefix,"TRANS.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mnp_INS <- read.table(paste(BIA_Minimap_Prefix,"INS.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mnp_DEL <- read.table(paste(BIA_Minimap_Prefix,"DEL.1kb_trsh.AllSVs",suffix, sep=""), h=T)
BIA_mnp_INS$Position <- as.integer(BIA_mnp_INS$Position)

# Original Nb INS and DEL (before merging based on the localisation)
BIA_mnp_original_Nb_INS = 74331
BIA_mnp_original_Nb_DEL = 109135
BIA_mnp_original_sum_INS = 17139392
BIA_mnp_original_sum_DEL = 23257913


##########################################
# ORGANIZE DATASET
###########################################

add_meta <- function(df, ref, aligner, svtype) {
  df %>%
    mutate(
      Reference = ref,
      Aligner = aligner,
      SV = svtype
    )
}

# CLA - MUMmer
CLA_mum_all <- bind_rows(
  add_meta(CLA_mum_INV, "CLA", "MUMmer", "INV"),
  add_meta(CLA_mum_DUP, "CLA", "MUMmer", "DUP"),
  add_meta(CLA_mum_INVDP, "CLA", "MUMmer", "INVDP"),
  add_meta(CLA_mum_INVTR, "CLA", "MUMmer", "INVTR"),
  add_meta(CLA_mum_TRANS, "CLA", "MUMmer", "TRANS"),
  add_meta(CLA_mum_INS, "CLA", "MUMmer", "INS"),
  add_meta(CLA_mum_DEL, "CLA", "MUMmer", "DEL")
)

# CLA - minimap
CLA_mnp_all <- bind_rows(
  add_meta(CLA_mnp_INV, "CLA", "Minimap", "INV"),
  add_meta(CLA_mnp_DUP, "CLA", "Minimap", "DUP"),
  add_meta(CLA_mnp_INVDP, "CLA", "Minimap", "INVDP"),
  add_meta(CLA_mnp_INVTR, "CLA", "Minimap", "INVTR"),
  add_meta(CLA_mnp_TRANS, "CLA", "Minimap", "TRANS"),
  add_meta(CLA_mnp_INS, "CLA", "Minimap", "INS"),
  add_meta(CLA_mnp_DEL, "CLA", "Minimap", "DEL")
)

# BIA - MUMmer
BIA_mum_all <- bind_rows(
  add_meta(BIA_mum_INV, "BIA", "MUMmer", "INV"),
  add_meta(BIA_mum_DUP, "BIA", "MUMmer", "DUP"),
  add_meta(BIA_mum_INVDP, "BIA", "MUMmer", "INVDP"),
  add_meta(BIA_mum_INVTR, "BIA", "MUMmer", "INVTR"),
  add_meta(BIA_mum_TRANS, "BIA", "MUMmer", "TRANS"),
  add_meta(BIA_mum_INS, "BIA", "MUMmer", "INS"),
  add_meta(BIA_mum_DEL, "BIA", "MUMmer", "DEL")
)

# BIA - minimap
BIA_mnp_all <- bind_rows(
  add_meta(BIA_mnp_INV, "BIA", "Minimap", "INV"),
  add_meta(BIA_mnp_DUP, "BIA", "Minimap", "DUP"),
  add_meta(BIA_mnp_INVDP, "BIA", "Minimap", "INVDP"),
  add_meta(BIA_mnp_INVTR, "BIA", "Minimap", "INVTR"),
  add_meta(BIA_mnp_TRANS, "BIA", "Minimap", "TRANS"),
  add_meta(BIA_mnp_INS, "BIA", "Minimap", "INS"),
  add_meta(BIA_mnp_DEL, "BIA", "Minimap", "DEL")
)

# Combine all
SV_all <- bind_rows(CLA_mum_all, CLA_mnp_all, BIA_mum_all, BIA_mnp_all)

# Apply mapping of chromsomes to BIA dataset
SV_all <- SV_all %>%
  left_join(Ref_chromosome_Link_file, by = c("Chrom" = "BIA_chr")) %>%
  mutate(
    Chrom = ifelse(Reference == "BIA", CLA_chr, Chrom)
  ) %>%
  select(-CLA_chr)

# Fix chrom order
SV_all$Chrom <- factor(SV_all$Chrom,
                       levels = paste0("chr", sort(as.numeric(gsub("chr", "", unique(SV_all$Chrom))))))

# Fix SV_type
SV_all <- SV_all %>%
  mutate(
    SVType = ifelse(grepl("^DEL", ID), "DEL", SV)
  )  %>%
  mutate(
    SVType = ifelse(grepl("^INS", ID), "INS", SV)
  ) 

# Keep only major SV
SV_all_filt <- SV_all  %>%
  filter(SVType %in% c("INS", "DEL", "DUP", "INV", "TRANS"))


SV_order <- c("INV", "DUP", "TRANS", "INS", "DEL")
SV_all_filt$SV <- factor(SV_all_filt$SV, levels = SV_order)

# Modify SV_Length if SV_Mean_Length

if (Mean_Length) {
  SV_all_filt <- SV_all_filt %>%
    rename(SVLength = SV_Mean_Length) 
}

if (Mean_Length) {
  SV_all <- SV_all %>%
    rename(SVLength = SV_Mean_Length) 
}



##########################################
# OVERALL STATISTICS DEPENDING ON REFERENCE AND ALIGNER
###########################################

summary_global <- SV_all_filt %>%
  group_by(Reference, Aligner, SV) %>%
  summarise(
    Count = n(),
    Total_Length = sum(SVLength, na.rm = TRUE),
    .groups = "drop"
  )
summary_global$SV <- factor(summary_global$SV, levels = SV_order)

summary_chr <- SV_all_filt %>%
  group_by(Chrom, Reference, Aligner, SV) %>%
  summarise(
    Count = n(),
    Total_Length = sum(SVLength, na.rm = TRUE),
    .groups = "drop"
  )

summary_chr$SV <- factor(summary_chr$SV, levels = SV_order)

summary_global_allSV <- SV_all %>%
  group_by(Reference, Aligner, SV) %>%
  summarise(
    Count = n(),
    Total_Length = sum(SVLength, na.rm = TRUE),
    .groups = "drop"
  )

summary_chr_allSV <- SV_all %>%
  group_by(Chrom, Reference, Aligner, SV) %>%
  summarise(
    Count = n(),
    Total_Length = sum(SVLength, na.rm = TRUE),
    .groups = "drop"
  )
summary_chr_allSV$SV <- factor(summary_chr_allSV$SV, levels = c("INV", "DUP", "TRANS", "INVDP",
                                                                "INVTR", "INS", "DEL"))

# Total count per Reference and Aligner
Total_SV <- summary_global_allSV %>%
  group_by(Reference, Aligner) %>%
  summarise(
    Total_Count = sum(Count),
    Total_Length = sum(Total_Length),
    .groups = "drop"
  )

Total_SV_NoComplexSV <- summary_global %>%
  group_by(Reference, Aligner) %>%
  summarise(
    Total_Count = sum(Count),
    Total_Length = sum(Total_Length),
    .groups = "drop"
  )


Total_SV_NoComplexSV
Total_SV

#######
# "Clean" Dataset without DTR
######

# Remove duplicated species in the "Species" columns
SV_clean_all <- SV_all %>%
  mutate(Species = str_split(Species, ",")) %>%
  mutate(Species = lapply(Species, function(x) setdiff(x, "DTR"))) %>%
  mutate(Species = sapply(Species, paste, collapse = ",")) %>%
  filter(Species != "" & !is.na(Species))

SV_clean_filt  <- SV_all_filt %>%
  mutate(Species = str_split(Species, ",")) %>%
  mutate(Species = lapply(Species, function(x) setdiff(x, "DTR"))) %>%
  mutate(Species = sapply(Species, paste, collapse = ",")) %>%
  filter(Species != "" & !is.na(Species))

Total_SV_NoDTR <- SV_clean_all %>%
  group_by(Reference, Aligner, SV) %>%
  summarise(
    Count = n(),
    Total_Length = sum(SVLength, na.rm = TRUE),
    .groups = "drop") %>%
  group_by(Reference, Aligner) %>%
  summarise(
    Total_Count = sum(Count),
    Total_Length = sum(Total_Length),
    .groups = "drop")

Total_SV_NoComplexSV_NoDTR <- SV_clean_filt %>%
  group_by(Reference, Aligner, SV) %>%
  summarise(
    Count = n(),
    Total_Length = sum(SVLength, na.rm = TRUE),
    .groups = "drop") %>%
  group_by(Reference, Aligner) %>%
  summarise(
    Total_Count = sum(Count),
    Total_Length = sum(Total_Length),
    .groups = "drop")

Total_SV_NoComplexSV_NoDTR
Total_SV_NoDTR



#########
## PLOTS OVERALL STATISTICS
#########

#####
# Overall number and Length of SVs
#####

p1 <- ggplot(summary_global, aes(x = SV, y = Count, fill = Reference)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Aligner) +
  theme_bw() +
  labs(title = "Number of SVs", y = "Count")

# Total length of SVs
p2 <- ggplot(summary_global, aes(x = SV, y = Total_Length, fill = Reference)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Aligner) +
  theme_bw() +
  labs(title = "Total SV Length", y = "Total Length (bp)")

#pdf("7_Stats_and_Plots_1kb/8b_MergedStats_NumberSV_MeanLength.pdf", width = 7, height = 5)
pdf("7_Stats_and_Plots_1kb/8_MergedStats_NumberSV.pdf", width = 7, height = 5)
p1
dev.off()
#browseURL("7_Stats_and_Plots_1kb/8_MergedStats_NumberSV.pdf")

#pdf("7_Stats_and_Plots_1kb/9b_MergedStats_LengthSV_MeanLength.pdf", width = 7, height = 5)
pdf("7_Stats_and_Plots_1kb/9_MergedStats_LengthSV.pdf", width = 7, height = 5)
p2
dev.off()
#browseURL("7_Stats_and_Plots_1kb/9_MergedStats_LengthSV.pdf")


#pdf("7_Stats_and_Plots_1kb/10b_MergedStats_Nb_and_Len_SV_MeanLength.pdf", width = 7, height = 7)
pdf("7_Stats_and_Plots_1kb/10_MergedStats_Nb_and_Len_SV.pdf", width = 7, height = 7)
p1 / p2  #
dev.off()
browseURL("7_Stats_and_Plots_1kb/10_MergedStats_Nb_and_Len_SV.pdf")

#####
# Number and Length By chromosomes
#####
SV_colors <- c(
  INV   = "#F09837",
  DUP   = "#56BBF9",
  TRANS = "#B7D86E",
  INS   = "grey",
  DEL   = "grey60"
)


# Total number of SV by chromosomes
p3 <- ggplot(summary_chr, aes(x = Chrom, y = Count, fill = SV)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = SV_colors, drop = FALSE) +
  facet_grid(Aligner ~ Reference) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "SV Count per Chromosome (CLA-aligned coordinates)")


p4 <- ggplot(summary_chr, aes(x = Chrom, y = Total_Length, fill = SV)) +
  geom_bar(stat = "identity")  +
  scale_fill_manual(values = SV_colors, drop = FALSE) +
  facet_grid(Reference ~ Aligner) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "SV Length per Chromosome (CLA-aligned coordinates)")


#pdf("7_Stats_and_Plots_1kb/11b_MergedStats_NumberSV_ByChrom_MeanLength.pdf", width = 7, height = 5)
pdf("7_Stats_and_Plots_1kb/11_MergedStats_NumberSV_ByChrom.pdf", width = 7, height = 5)
p3
dev.off()
#browseURL("7_Stats_and_Plots_1kb/11_MergedStats_NumberSV_ByChrom.pdf")

#pdf("7_Stats_and_Plots_1kb/12b_MergedStats_LengthSV_ByChrom_MeanLength.pdf", width = 7, height = 5)
pdf("7_Stats_and_Plots_1kb/12_MergedStats_LengthSV_ByChrom.pdf", width = 7, height = 5)
p4
dev.off()
#browseURL("7_Stats_and_Plots_1kb/12_MergedStats_LengthSV_ByChrom.pdf")

#pdf("7_Stats_and_Plots_1kb/13b_MergedStats_Nb_and_Len_SV_ByChrom_MeanLength.pdf", width = 12, height =12)
pdf("7_Stats_and_Plots_1kb/13_MergedStats_Nb_and_Len_SV_ByChrom.pdf", width = 12, height =12)
p3 / p4  #
dev.off()
#browseURL("7_Stats_and_Plots_1kb/13_MergedStats_Nb_and_Len_SV_ByChrom.pdf")


#####
# Overall number of SV  per chromsomes, with separated statistics
#####

df_total <- summary_chr_allSV %>%
  group_by(Chrom, Reference, Aligner) %>%
  summarise(
    Count = sum(Count),
    .groups = "drop"
  ) %>%
  mutate(SV = "TOTAL")

summary_chr_allSV_withTotal <- bind_rows(summary_chr_allSV, df_total) %>%
  mutate(Group = paste(Reference, Aligner, sep = " - "))

summary_chr_allSV_withTotal$SV <- factor(summary_chr_allSV_withTotal$SV,
       levels = c("TOTAL", "DEL", "DUP", "INS", "INV", "INVDP", "INVTR","TRANS"))

p5 <- ggplot(summary_chr_allSV_withTotal, aes(x = Chrom, y = Count, fill = Group)) +
  geom_col(position = "dodge") +
  labs(
    title = "SV Counts per Chromosome",
    x = "Chromosome",
    y = "Count",
    fill = "Reference + Aligner"
  ) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("7_Stats_and_Plots_1kb/14_MergedStats_Nb_SV_ByChrom_BySV.pdf", width = 12, height =12)
p5 + facet_wrap(~ SV, scales = "free_y", ncol = 2) #
dev.off()

pdf("7_Stats_and_Plots_1kb/14b_MergedStats_Nb_SV_ByChrom_BySV.pdf", width = 12, height =12)
p5 + facet_wrap(~ SV, ncol = 2) #
dev.off()

#####
# Check difference of INS and DEL depending on merging at the threshold
# or overall count
#####

INDELS_noMerging <- tribble(
  ~Reference, ~Aligner, ~SV,  ~Count,   ~Total_Length,
  "CLA", "Minimap", "INS", 132344, 34250341,
  "CLA", "Minimap", "DEL", 159850, 32038807,
  "CLA", "MUMmer",  "INS", 64421,  19053203,
  "CLA", "MUMmer",  "DEL", 87771,  18547690,
  "BIA", "Minimap", "INS", 74331,  17139392,
  "BIA", "Minimap", "DEL", 109135, 23257913,
  "BIA", "MUMmer",  "INS", 21554,  5067895,
  "BIA", "MUMmer",  "DEL", 27344,  5435675
)


INDELS_noMerging <- INDELS_noMerging %>%
  mutate(Dataset = "Original")

summary_global_INDELS <- summary_global %>%
  filter(SV %in% c("INS", "DEL")) %>%
  rename(Count = Count,
         Total_Length = Total_Length) %>%
  mutate(Dataset = "Merged 10kb")

INDELS_Combined <-  bind_rows(INDELS_noMerging, 
                             summary_global_INDELS)

INDELS_long <- INDELS_Combined %>%
  pivot_longer(cols = c(Count, Total_Length),
               names_to = "Metric",
               values_to = "Value")

p_indel <- ggplot(INDELS_long,
         aes(x = interaction(Reference, Aligner),
             y = Value,
             fill = Dataset)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.8)) +
    facet_grid(Metric ~ SV, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = "Original vs Merged (10kb) INS and DEL",
         x = "Reference - Aligner",
         y = NULL) +
    aes(x = paste(Reference, Aligner, sep = " - "))

pdf("7_Stats_and_Plots_1kb/15_MergedStats_INDELS.pdf", width = 10, height =8)
p_indel
dev.off()


##########################################
# PROPORTION OF SV (Number and length)
# Different Reference and Aligner
###########################################

SV_prop <- summary_global %>%
  group_by(Aligner, Reference) %>%
  mutate(Proportion = Count / sum(Count),
         Metric = paste(Aligner, Reference, sep = " - "))

SV_prop_len <- summary_global %>%
  group_by(Aligner, Reference) %>%
  mutate(Proportion = Total_Length / sum(Total_Length),
         Metric = paste(Aligner, Reference, sep = " - "))


# Proportion of SV
p_prop <- ggplot(SV_prop,aes(x = Metric, y = Proportion,
                             fill = SV)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "INV" = "#F09837", "DUP" = "#56BBF9", "TRANS" = "#B7D86E",
    "INS" = "grey","DEL" = "grey60")) +
  labs(
    title = "Proportions of Structural Variant Types",
    x = "", y = "Proportion", fill = "SV Type" ) +
  theme_minimal() + coord_flip() +
  theme(axis.text.y = element_text(size = 12))

p_prop_len <- ggplot(SV_prop_len, aes(x = Metric, y = Proportion,
                                      fill = SV)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "INV" = "#F09837", "DUP" = "#56BBF9", "TRANS" = "#B7D86E",
    "INS" = "grey","DEL" = "grey60")) +
  labs(
    title = "Proportions of Structural Variant Types",
    x = "", y = "Proportion", fill = "SV Type") +
  theme_minimal() + coord_flip() +
  theme(axis.text.y = element_text(size = 12))

pdf("7_Stats_and_Plots_1kb/16_Proportion_SV.pdf", width = 10, height =8)
p_prop / p_prop_len
dev.off()


#### Count and Length, side by side

SV_prop_V2 <- summary_global %>%
  pivot_longer(cols = c(Count, Total_Length),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                    Count = "Count",
                    Total_Length = "Total Length"),
    Group = paste(Reference, Aligner, sep = " - ")) %>%
  group_by(Group, Metric) %>%
  mutate(Proportion = Value / sum(Value))

p_prop2 <- ggplot(SV_prop_V2,
       aes(x = Metric, y = Proportion, fill = SV)) +
  geom_bar(stat = "identity", width = 0.6) +
  facet_wrap(~Group, nrow = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "INV" = "#F09837", "DUP" = "#56BBF9", "TRANS" = "#B7D86E",
    "INS" = "grey", "DEL" = "grey60" )) +
  labs(
    title = "Proportion of SV Types: Count vs Total Length",
    x = "", y = "Proportion", fill = "SV Type" ) +
  theme_minimal() + coord_flip() +
  theme(axis.text.x = element_text(size = 12))

pdf("7_Stats_and_Plots_1kb/16b_Proportion_SV.pdf", width = 9, height =7)
p_prop2
dev.off()

##########################################
# PROPORTION OF SV BY SPECIES
# MONOPHYLETIC / PARAPHYLETIC / SINGLESPECIES
###########################################

tree <- read.tree("PacBioSpecies.Tree.tree")

is_monophyletic_set <- function(species, tree) {
  species <- unlist(strsplit(species, ","))
  species <- intersect(species, tree$tip.label)
  
  if (length(species) <= 1) return(TRUE)
  
  mrca_node <- getMRCA(tree, species)
  clade_tips <- extract.clade(tree, mrca_node)$tip.label
  
  return(setequal(species, clade_tips))
}

# No INVDP, INVTR, ETC: SV_all_filt
# Remove duplicated species in the "Species" columns
SV_clean <- SV_all %>%
  select(-SVCat)  %>%
  mutate(Species = str_split(Species, ",")) %>%
  mutate(Species = lapply(Species, function(x) unique(trimws(x)))) %>%
  mutate(Species = sapply(Species, paste, collapse = ","))

# Remove DTR as not considered in all analyses and not clownfish (and not in the tree)
SV_clean <- SV_clean %>%
  mutate(Species = str_split(Species, ",")) %>%
  mutate(Species = lapply(Species, function(x) setdiff(x, "DTR"))) %>%
  mutate(Species = sapply(Species, paste, collapse = ",")) %>%
  filter(Species != "" & !is.na(Species))

SV_clean <- SV_clean %>%
  mutate(n_species = sapply(str_split(Species, ","), length)) %>%
  mutate(Type = ifelse(n_species == 1, "Monomorphic", "Polymorphic"))
 
# Get New category for Monophyletic/Paraphyletic
SV_clean$Monophyletic <- sapply(
  SV_clean$Species,
  is_monophyletic_set,
  tree = tree
)

########## 
##### Plot
########## 

# Prepare the data
SV_plot_type_phylo <- SV_clean %>%
  mutate(
    Type = ifelse(n_species == 1, "Monomorphic", "Polymorphic"),
    Phylo = ifelse(Monophyletic, "Monophyletic", "Paraphyletic")
  )

# Create a stacking variable
SV_plot_type_phylo <- SV_plot_type_phylo %>%
  mutate(
    Stack = case_when(
      Type == "Monomorphic" ~ "Monomorphic",
      Type == "Polymorphic" & Phylo == "Monophyletic" ~ "Monophyletic",
      Type == "Polymorphic" & Phylo == "Paraphyletic" ~ "Paraphyletic"
    )
  )

# Order Factor
SV_plot_type_phylo$Type <- factor(SV_plot_type_phylo$Type, levels = c("Monomorphic", "Polymorphic"))

SV_plot_type_phylo$Stack <- factor(
  SV_plot_type_phylo$Stack,
  levels = c("Monomorphic", "Monophyletic", "Paraphyletic")
)

plot_type_phylo <-  ggplot(SV_plot_type_phylo, aes(x = Type, fill = Stack)) +
  geom_bar() +
  facet_grid(Reference ~ Aligner) +
  theme_bw() +
  labs(
    title = "Monomorphic vs Polymorphic variants",
    x = "",
    y = "Count",
    fill = "Category"
  )

plot_type_phylo_bySV <- ggplot(SV_plot_type_phylo, aes(x = Type, fill = Stack)) +
  geom_bar() +
  facet_grid(Reference ~ Aligner + SVType) +
  theme_bw() +
  labs(
    title = "Variant classification by SV type",
    x = "",
    y = "Count",
    fill = "Category"
  )


pdf("7_Stats_and_Plots_1kb/17_Proportion_Monomorphic_Polymorphic.pdf", width = 7, height =8)
plot_type_phylo
dev.off()

pdf("7_Stats_and_Plots_1kb/17b_Proportion_Monomorphic_Polymorphic_BySV.pdf", width = 14, height =8)
plot_type_phylo_bySV + theme(
  axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



##########################################
### Number of private and shared SV
###########################################

SV_plot_private_shared <- SV_clean %>%
  mutate(
    SpeciesGroup = ifelse(n_species == 1, "1", ">1")
  )

pdf("7_Stats_and_Plots_1kb/18_NuberSpecies_perSV.pdf", width = 9, height =8)
ggplot(SV_plot_private_shared, aes(x = factor(n_species), fill = SpeciesGroup)) +
  geom_bar() +
  facet_grid(Reference ~ Aligner, scales = "free_y") +
  scale_fill_manual(
    values = c("1" = "grey30", ">1" = "grey60"),
    guide = "none"
  ) +
  theme_bw() +
  labs(
    title = "Number of species per SV",
    x = "Number of species",
    y = "Count"
  )
dev.off()


pdf("7_Stats_and_Plots_1kb/18b_NuberSpecies_perSV_BySVType.pdf", width = 15, height =11)
ggplot(SV_plot_private_shared, aes(x = factor(n_species), fill = SpeciesGroup)) +
  geom_bar(color = "black") +
  facet_grid(SVType ~ Reference + Aligner, scales = "free_y") +
  scale_fill_manual(
    values = c("1" = "grey30", ">1" = "grey80"),
    guide = "none"
  ) +
  theme_bw() +
  labs(
    title = "Number of species per SV by SV type",
    x = "Number of species",
    y = "Count"
  )
dev.off()


##########################################
# SV Length
###########################################

plot_SVlength <- ggplot(SV_all_filt, aes(x = "", y = SVLength)) +
  geom_boxplot() +
  facet_grid(Reference ~ Aligner) +
  theme_bw() +
  labs( title = "SV length distribution", x = "",y = "SV length")

plot_SVlength_Log <- ggplot(SV_all_filt, aes(x = "", y = SVLength)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_grid(Reference ~ Aligner) +
  theme_bw() +
  labs(title = "SV length distribution (log scale)", x = "", y = "SV length (log10)")

plot_SVlength_Log_bySV <- ggplot(SV_all_filt, aes(x = SVType, y = SVLength, fill = SVType)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_grid(Reference ~ Aligner) +
  theme_bw() +
  labs(title = "SV length by type", x = "SV type", y = "SV length (log10)")

# All SV
plot_allSVlength <- ggplot(SV_clean, aes(x = "", y = SVLength)) +
  geom_boxplot() +
  facet_grid(Reference ~ Aligner) +
  theme_bw() +
  labs( title = "SV length distribution", x = "",y = "SV length")

plot_allSVlength_Log <- ggplot(SV_clean, aes(x = "", y = SVLength)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_grid(Reference ~ Aligner) +
  theme_bw() +
  labs(title = "SV length distribution (log scale)", x = "", y = "SV length (log10)")

plot_allSVlength_Log_bySV <- ggplot(SV_clean, aes(x = SVType, y = SVLength, fill = SVType)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_grid(Reference ~ Aligner) +
  theme_bw() +
  labs(title = "SV length by type", x = "SV type", y = "SV length (log10)")

pdf("7_Stats_and_Plots_1kb/19_SVLength_BySVType.pdf", width = 15, height =11)
plot_SVlength_Log_bySV
dev.off()

pdf("7_Stats_and_Plots_1kb/19b_SVLength_BySVType_allSV.pdf", width = 15, height =11)
plot_allSVlength_Log_bySV
dev.off()

#########################
#### CLUSTERINC (PCA and PCoA)
#########################

##############
#### Prepare data
##############

# Information for coloration
clade_info <- data.frame(
  Species = c("FRE","EPH","CLA",
              "AKA","SAN","PRD",
              "OMA","ALL","LAT",
              "SEB","POL",
              "CRP","MCC","AKY",
              "LAZ",
              "PRC","OCE",
              "BIA"),
  Clade = c(
    "Clade1","Clade1","Clade1",
    "Clade2","Clade2","Clade2",
    "Clade3","Clade3","Clade3",
    "Clade4","Clade4",
    "Clade5","Clade5","Clade5",
    "Clade6",
    "Clade7","Clade7",
    "Clade8"
  )
)

clade_colors <- c(
  "Clade1" = "#F28E2B",  
  "Clade2" = "#EDC948",  
  "Clade3" = "#9C4DCC",  
  "Clade4" = "#3283C9", 
  "Clade5" = "#0025DB",
  "Clade6" = "#09850E",
  "Clade7" = "#97B898", 
  "Clade8" = "#66DB00"
)

# Function to generate shades
generate_clade_colors <- function(species, clade, base_color) {
  n <- length(species)
  shades <- lighten(base_color, seq(0, 0.4, length.out = n))
  names(shades) <- species
  return(shades)
}

# Generate the colors
species_colors <- c()
for (cl in unique(clade_info$Clade)) {
  sp <- clade_info$Species[clade_info$Clade == cl]
  base_col <- clade_colors[cl]
  
  species_colors <- c(
    species_colors,
    generate_clade_colors(sp, cl, base_col)
  )
}

# Remove BIA when BIARef
SV_clean_NoBIA <- SV_clean %>%
  mutate(Species = str_split(Species, ",")) %>%
  mutate(Species = mapply(function(species, ref) {
    if (ref == "BIA") {
      species <- setdiff(species, "BIA")
    }
    return(species)
  }, Species, Reference, SIMPLIFY = FALSE)) %>%
  mutate(Species = sapply(Species, paste, collapse = ",")) %>%
  filter(Species != "")

SV_clean_NoBIA_filt <- SV_clean_NoBIA  %>%
  filter(SVType %in% c("INS", "DEL", "DUP", "INV", "TRANS"))

##############
#### Do PCA
##############

pca_list <- list()
for (ref in unique(SV_clean_NoBIA$Reference)) {
  for (alg in unique(SV_clean_NoBIA$Aligner)) {
    
    # Subset
    df_sub <- SV_clean_NoBIA %>%
      filter(Reference == ref, Aligner == alg)
    
    # Skip if too few SVs
    if (nrow(df_sub) < 2) next
    
    # Expand species
    df_long <- df_sub %>%
      separate_rows(Species, sep = ",")
    
    # Presence/absence matrix
    SV_pa <- df_long %>%
      mutate(presence = 1) %>%
      select(ID, Species, presence) %>%
      distinct() %>%
      pivot_wider(
        names_from = Species,
        values_from = presence,
        values_fill = 0)
    
    # Convert to matrix
    mat <- as.data.frame(SV_pa)
    rownames(mat) <- mat$ID
    mat$ID <- NULL
    
    # Transpose (species × SVs)
    mat_t <- t(mat)
    mat_t <- mat_t[, apply(mat_t, 2, var) != 0]
    
    # Skip if too few species
    if (nrow(mat_t) < 2) next
    
    # PCA
    pca_res <- prcomp(mat_t, scale. = TRUE)
    
    # Store results
    pca_df <- as.data.frame(pca_res$x)
    pca_df$Species <- rownames(pca_df)
    pca_df$Reference <- ref
    pca_df$Aligner <- alg
    
    pca_list[[paste(ref, alg, sep = "_")]] <- pca_df
  }
}
# Combine all PCA results
pca_all <- bind_rows(pca_list)

pca_all <- pca_all %>%
  left_join(clade_info, by = "Species")

plot_pca <- ggplot(pca_all, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Species), max.overlaps = Inf) +
  facet_grid(Reference ~ Aligner) +
  scale_color_manual(values = species_colors) +
  theme_bw()

plot_pca

#########################
#### Do PCoA (or MDS)
#########################

pcoa_list <- list()
for (ref in unique(SV_clean_NoBIA$Reference)) {
  for (alg in unique(SV_clean_NoBIA$Aligner)) {
    
    # Subset
    df_sub <- SV_clean_NoBIA %>%
      filter(Reference == ref, Aligner == alg)
    
    if (nrow(df_sub) < 2) next
    
    # Expand species
    df_long <- df_sub %>%
      separate_rows(Species, sep = ",")
    
    # Presence/absence matrix
    SV_pa <- df_long %>%
      mutate(presence = 1) %>%
      select(ID, Species, presence) %>%
      distinct() %>%
      pivot_wider(
        names_from = Species,
        values_from = presence,
        values_fill = 0
      )
    
    # Convert to matrix
    mat <- as.data.frame(SV_pa)
    rownames(mat) <- mat$ID
    mat$ID <- NULL
    
    # Transpose (species × SVs)
    mat_t <- t(mat)
    
    # Skip if not enough species
    if (nrow(mat_t) < 2) next
    
    # Jaccard distance
    dist_mat <- vegdist(mat_t, method = "jaccard")
    
    # PCoA
    pcoa_res <- cmdscale(dist_mat, k = 4, eig = TRUE)
    
    # Store results
    pcoa_df <- as.data.frame(pcoa_res$points)
    colnames(pcoa_df) <- c("PCoA1", "PCoA2","PCoA3", "PCoA4")
    pcoa_df$Species <- rownames(pcoa_df)
    pcoa_df$Reference <- ref
    pcoa_df$Aligner <- alg
    
    # Variance explained (important!)
    eig_vals <- pcoa_res$eig
    var_explained <- eig_vals / sum(eig_vals)
    
    pcoa_df$Var1 <- var_explained[1]
    pcoa_df$Var2 <- var_explained[2]
    
    pcoa_list[[paste(ref, alg, sep = "_")]] <- pcoa_df
  }
}

# Combine all
pcoa_all <- bind_rows(pcoa_list)
pcoa_all <- pcoa_all %>%
  left_join(clade_info, by = "Species")

# Plot all together
plot_MDS <- ggplot(pcoa_all, aes(x = PCoA1, y = PCoA2, color = Species)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Species), size = 3, max.overlaps = Inf) +
  facet_grid(Reference ~ Aligner) +
  scale_color_manual(values = species_colors) +
  theme_bw() + theme(legend.position = "none")

plot_MDS_MDS3_4 <- ggplot(pcoa_all, aes(x = PCoA3, y = PCoA4, color = Species)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Species), size = 3, max.overlaps = Inf) +
  facet_grid(Reference ~ Aligner) +
  scale_color_manual(values = species_colors) +
  theme_bw() + theme(legend.position = "none")

pdf("7_Stats_and_Plots_1kb/20_MDS_AllSV.pdf", width = 8, height =8)
plot_MDS
dev.off()

pdf("7_Stats_and_Plots_1kb/20b_MDS_MDS3_4_AllSV.pdf", width = 8, height =8)
plot_MDS_MDS3_4
dev.off()

# Plot separately
plot_pcoa <- function(df) {
  
  ref <- unique(df$Reference)
  alg <- unique(df$Aligner)
  
  ggplot(df, aes(x = PCoA1, y = PCoA2, color = Species)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = Species), size = 3, max.overlaps = Inf) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      title = paste("PCoA -", ref, "-", alg),
      x = paste0("PCoA1 (", round(df$Var1[1] * 100, 1), "%)"),
      y = paste0("PCoA2 (", round(df$Var2[1] * 100, 1), "%)")
    )
}
plots <- split(pcoa_all, list(pcoa_all$Reference, pcoa_all$Aligner))
plot_list <- lapply(plots, plot_pcoa)

pdf("7_Stats_and_Plots_1kb/20c_MDS_AllSV_MethodsSeparately.pdf", width = 8, height =8)
for (p in plot_list) {
  print(p)
}
dev.off()

#########################
#### Do PCoA (or MDS) By SV Types
#########################

groups <- SV_clean_NoBIA_filt %>%
  group_by(Reference, Aligner, SVType) %>%
  group_split()

pcoa_list_SVtype <- list()
for (df_sub in groups) {
  
  if (nrow(df_sub) < 2) next
  
  ref <- unique(df_sub$Reference)
  alg <- unique(df_sub$Aligner)
  svt <- unique(df_sub$SVType)
  
  df_long <- df_sub %>%
    separate_rows(Species, sep = ",")
  
  SV_pa <- df_long %>%
    mutate(presence = 1) %>%
    select(ID, Species, presence) %>%
    distinct() %>%
    pivot_wider(
      names_from = Species,
      values_from = presence,
      values_fill = 0
    )
  
  mat <- as.data.frame(SV_pa)
  rownames(mat) <- mat$ID
  mat$ID <- NULL
  
  mat_t <- t(mat)
  
  if (nrow(mat_t) < 2) next
  
  dist_mat <- vegdist(mat_t, method = "jaccard")
  pcoa_res <- cmdscale(dist_mat, k = 4, eig = TRUE)
  
  pcoa_df <- as.data.frame(pcoa_res$points)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2","PCoA3", "PCoA4")
  pcoa_df$Species <- rownames(pcoa_df)
  pcoa_df$Reference <- ref
  pcoa_df$Aligner <- alg
  pcoa_df$SVType <- svt
  
  eig_vals <- pcoa_res$eig
  var_explained <- eig_vals / sum(eig_vals)
  
  pcoa_df$Var1 <- var_explained[1]
  pcoa_df$Var2 <- var_explained[2]
  
  pcoa_list_SVtype[[paste(ref, alg, svt, sep = "_")]] <- pcoa_df
}

pcoa_all_SVtype <- do.call(rbind, pcoa_list_SVtype)

plot_MDS_SVtype <- ggplot(pcoa_all_SVtype, aes(x = PCoA1, y = PCoA2, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = species_colors) +
  facet_grid(SVType ~ Reference + Aligner) +
  theme_bw() +
  labs(
    title = "PCoA (Jaccard distance)",
    x = "PCoA1",
    y = "PCoA2"
  )

pdf("7_Stats_and_Plots_1kb/20d_MDS_BySVType.pdf", width = 8, height =8)
plot_MDS_SVtype
dev.off()

# plots separately per SVType
plot_list_ByType <- list()
for (svt in unique(pcoa_all_SVtype$SVType)) {
  
  df_plot <- pcoa_all_SVtype %>%
    filter(SVType == svt)
  
  p <- ggplot(df_plot, aes(x = PCoA1, y = PCoA2, color = Species)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = Species), size = 3, max.overlaps = Inf) +
    
    facet_grid(Reference ~ Aligner) +
    
    scale_color_manual(values = species_colors) +
    
    labs(
      title = paste("SVType:", svt),
      x = paste0("PCoA1"),
      y = paste0("PCoA2")
    ) +
    
    theme_bw() +
    theme(legend.position = "none")
  
  plot_list_ByType[[svt]] <- p
}

pdf("7_Stats_and_Plots_1kb/20e_MDS_BySVType_Separately.pdf", width = 8, height =8)
for (p in plot_list_ByType) {
  print(p)}
dev.off()
#plot_list_ByType[["INV"]]   # or "DUP", "INS", etc.


#########################
#### SV and Host Specialization
########################

################
#### Analyze the data
################

Species_categories <- c(
  "AKA"="RAD","AKY"="GEN","ALL"="GEN","BIA"="ENT",
  "CRP"="GEN","EPH"="ENT","FRE"="ENT","LAT"="GEN",
  "LAZ"="GEN","MCC"="ENT","OMA"="ENT","POL"="STI",
  "PRC"="RAD","PRD"="RAD","SAN"="STI","SEB"="GEN",
  "OCE"="RAD", "CLA"="GEN")


# Number of species per host
nb_species_per_host <- table(Species_categories)
nb_species_per_host

# Refilter
SV_clean_NoBIA_SpcList <- SV_clean_NoBIA_filt %>%
  mutate(Species_list = str_split(Species, ",")) %>% 
  unnest(Species_list) %>%
  mutate(Species_list = str_trim(Species_list))%>%
  # Assign SV_clean_NoBIA_SpcList
  mutate(Host = Species_categories[Species_list])

host_counts <- SV_clean_NoBIA_SpcList %>%
  distinct(Reference, Aligner, Species_list, Host) %>%
  count(Reference, Aligner, Host, name = "n_species_host") %>%
  # We need to add +1 for BIA and CLA in the reference to
  # avoid having, for instance ENT-specific SV with BIARef (which is not correct)
  mutate(
    n_species_host = case_when(
      Reference == "BIA" & Host == "ENT" ~ n_species_host + 1,
      Reference == "CLA"    & Host == "GEN" ~ n_species_host + 1,
      TRUE ~ n_species_host
    )
  )

# Reconstruct SV-level summaries
SV_clean_NoBIA_SpcList_Host_Summary <- SV_clean_NoBIA_SpcList %>%
  group_by(Reference, Aligner, ID, SVType) %>%
  summarise(
    species = list(unique(Species_list)),
    hosts = list(unique(Host)),
    n_species = n_distinct(Species_list),
    host_type = ifelse(length(unique(Host)) == 1, unique(Host), NA),
    .groups = "drop"
  )
# Even if the number of SV per Host category depends on the Reference used (BIA vs CLA), 
# we should still be able to compare the overall number of species per Host Cat.
# This is also because we added +1 in the count (depending on the ref)
# Example for BIARef: 
# On CLA reference: we need 5 ENT species to keep the SV. 
# On BIA ref: we will need 4 ENT more, but if we have the SV it means that BIA 
# does not have the SV, so we shuld not keep the SV. So we still can compare to "5" 
# The same works for generalists and CLA reference

SV_host_specific <- SV_clean_NoBIA_SpcList_Host_Summary %>%
  filter(
    !is.na(host_type),  # only one host
    n_species == nb_species_per_host[host_type]  # present in ALL species of that host
  )

# Count SV
SV_counts <- SV_host_specific %>%
  count(Reference, Aligner, host_type, name = "SV_count")

SV_counts_bytype <- SV_host_specific %>%
  count(Reference, Aligner, SVType, host_type, name = "SV_count")

plot_SV_by_HostType <- ggplot(SV_counts, aes(x = host_type, y = SV_count, fill = host_type)) +
  geom_col(color = "black") +  # black border helps especially for white bars
  facet_grid(Reference ~ Aligner) +
  scale_fill_manual(values = c(
    "GEN" = "#000000",   # black
    "ENT" = "#DF0C0C",   # red
    "RAD" = "#F3C337",   # yellow
    "STI" = "#FFFFFF"    # white
  )) +
  theme_minimal() +
  labs(title = "Host-specific SV per Reference × Aligner",
       x = "Host type",
       y = "Number of SV",
       fill = "Host") +
  theme(
    legend.position = "bottom")


plot_SV_by_HostType_BySVType <- ggplot(SV_counts_bytype,
       aes(x = host_type, y = SV_count, fill = SVType)) +
  geom_col(position = "stack") +
  facet_grid(Reference ~ Aligner) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "INV" = "#F09837", "DUP" = "#56BBF9", "TRANS" = "#B7D86E",
    "INS" = "grey", "DEL" = "grey60" ))+
  labs(title = "SV type composition per host",
       x = "Host",
       y = "SV count",
       fill = "SV Type")


pdf("7_Stats_and_Plots_1kb/21_SV_ByHost.pdf", width = 8, height =8)
plot_SV_by_HostType
dev.off()

pdf("7_Stats_and_Plots_1kb/21b_SV_ByHost_BySVType.pdf", width = 8, height =8)
plot_SV_by_HostType_BySVType
dev.off()


################
#### SV and HOSTS: RANDOM EXPECTATION
################

########
# STICHODACTYLA
########

df_pairs <- SV_clean_NoBIA_filt %>%
  mutate(Species_list = str_split(Species, ",")) %>%
  filter(map_int(Species_list, length) == 2)

df_pairs <- df_pairs %>%
  mutate(
    sp1 = map_chr(Species_list, 1),
    sp2 = map_chr(Species_list, 2))

df_pairs <- df_pairs %>%
  mutate(
    pair = map2_chr(sp1, sp2, ~ paste(sort(c(.x, .y)), collapse = "-")))

pair_counts <- df_pairs %>%
  count(Reference, Aligner,  pair, name = "count")

pair_counts_byType <- df_pairs %>%
  count(Reference, Aligner, SVType, pair, name = "count")

target_pairs <- c("POL-SAN","AKA-POL", "AKA-SEB", "SAN-SEB")

pair_counts_filtered <- pair_counts %>%
  filter(pair %in% target_pairs)  %>%
  mutate(pair = factor(pair, levels = target_pairs))

pair_counts_byType_filtered <- pair_counts_byType %>%
  filter(pair %in% target_pairs) %>%
  mutate(pair = factor(pair, levels = target_pairs))

plot_hostExpectation_STI <- 
  ggplot(pair_counts_filtered, aes(x = pair, y = count, fill = pair)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  facet_grid(Reference ~ Aligner)  +
  theme_minimal()+
  scale_fill_manual(values = c(
    "AKA-POL" = "gray80",   # black
    "AKA-SEB" = "gray80",   # red
    "POL-SAN" = "#FFFFFF",   # yellow
    "SAN-SEB" = "gray80"    # white
  ))


plot_hostExpectation_STI_bySVType <- 
  ggplot(pair_counts_byType_filtered, aes(x = pair, y = count, fill = SVType)) +
  geom_col(position = "stack", width = 0.8, color = "black") +
  facet_grid(Reference ~ Aligner)  +
  theme_minimal()+
  scale_fill_manual(values = c(
    "INV" = "#F09837", "DUP" = "#56BBF9", "TRANS" = "#B7D86E",
    "INS" = "grey", "DEL" = "grey60" ))


pdf("7_Stats_and_Plots_1kb/22_SV_Expectation_STI.pdf", width = 8, height =8)
plot_hostExpectation_STI
dev.off()

pdf("7_Stats_and_Plots_1kb/22b_SV_Expectation_STI_BySVType.pdf", width = 8, height =8)
plot_hostExpectation_STI_bySVType
dev.off()

########
# RADIANTHUS
########

df_fourSpecies <- SV_clean_NoBIA_filt %>%
  mutate(Species_list = str_split(Species, ",")) %>%
  filter(map_int(Species_list, length) == 4)

df_fourSpecies <- df_fourSpecies %>%
  mutate(
    clusters = map_chr(Species_list, ~ paste(sort(.x), collapse = "-"))
  )

fourspecies_counts <- df_fourSpecies  %>%
  count(Reference, Aligner, clusters, name = "count")

fourspecies_counts_byType <- df_fourSpecies %>%
  count(Reference, Aligner, SVType, clusters, name = "count")

# Get the Number with no host sharing but same phylogenetic distance

target_fourSpecies <- c("AKA-OCE-PRC-PRD", "OCE-PRC-PRD-SAN")

fourspecies_counts_filtered <- fourspecies_counts %>%
  filter(clusters %in% target_fourSpecies)  %>%
  mutate(clusters = factor(clusters, levels = target_fourSpecies))

fourspecies_counts_byType_filtered <- fourspecies_counts_byType %>%
  filter(clusters %in% target_fourSpecies) %>%
  mutate(clusters = factor(clusters, levels = target_fourSpecies))


plot_hostExpectation_RAD <- ggplot(fourspecies_counts_filtered, aes(x = clusters, y = count, fill = clusters)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  facet_grid(Reference ~ Aligner)  +
  theme_minimal()+
  scale_fill_manual(values = c(
    "OCE-PRC-PRD-SAN" = "gray80",  
    "AKA-OCE-PRC-PRD" = "#FFFFFF")) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) 

plot_hostExpectation_RAD_bySVType <- 
  ggplot(fourspecies_counts_byType_filtered, aes(x = clusters, y = count, fill = SVType)) +
  geom_col(position = "stack", width = 0.8, color = "black") +
  facet_grid(Reference ~ Aligner)  +
  theme_minimal()+
  scale_fill_manual(values = c(
    "INV" = "#F09837", "DUP" = "#56BBF9", "TRANS" = "#B7D86E",
    "INS" = "grey", "DEL" = "grey60" )) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) 


pdf("7_Stats_and_Plots_1kb/23_SV_Expectation_RAD.pdf", width = 8, height =8)
plot_hostExpectation_RAD
dev.off()

pdf("7_Stats_and_Plots_1kb/23b_SV_Expectation_RAD_BySVType.pdf", width = 8, height =8)
plot_hostExpectation_RAD_bySVType
dev.off()



################################################
#### INSIGHTS IN LONGEST SV
################################################

# Get 10 longest SV for each Reference and Aligner
top10_sv <- SV_clean %>%
  group_by(Reference, Aligner) %>%
  arrange(desc(SVLength), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write.table(top10_sv, "7_Stats_and_Plots_1kb/LongestSV_ByRefandMethod_MeanLength.txt",
            sep="\t", quote = F, row.names = F)

# SV longer than 1 Mb
long_sv_counts <- SV_clean %>%
  filter(SVLength > 1e6) %>%
  group_by(Reference, Aligner) %>%
  summarise(n_SV = n(), .groups = "drop")

# Around 100 per Aligner-Reference

sv_long_1mb <- SV_clean %>%
  filter(SVLength > 1e6) %>%
  group_by(Reference, Aligner) %>%
  arrange(desc(SVLength), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Chr_Pos = paste0(Chrom, ":", Position))

write.table(sv_long_1mb, "7_Stats_and_Plots_1kb/SV_1Mb_ByRefandMethod_MeanLength.txt",
            sep="\t", quote = F, row.names = F)
