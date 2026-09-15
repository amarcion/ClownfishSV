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
library(GenomicRanges)
library(rtracklayer)
library(cowplot)

setwd(".")

Ref_chromosome_Link_file = read.table("CLA_vs_BIARef.mapids.txt", h=F, 
                                      col.names= c("CLA_chr", "BIA_chr"))

BIA_Ref_chromosome_Link_file = read.table("BIARef_vs_BIARef.Chrom.txt", h=F, 
                                          col.names= c("BIA_Chr", "BIA_Scaff"))


sv_type_colours <-  c(
  INV   = "#F09837",
  DUP   = "#56BBF9",
  TRANS = "#B7D86E",
  INS   = "#D87B7B",
  DEL   = "#7B5EA7"
)

###############
### SETTINGS
###############

Reference_Species = "BIA"
out_prefix = "."

gff <- import("Premnas_biaculeatus.genome.gff")
genes <- gff[gff$type == "gene"]
cds   <- gff[gff$type == "CDS"]
exons <- gff[gff$type == "exon"]

##################
### DATA
##################

SV_all <- read.table("Final_Dataset_SV.Minimap.BIA_and_CLA.csv",
                     header = T)


SV_all <- SV_all %>% 
  filter(Reference==Reference_Species)

SV_private <- SV_all %>% filter(n_species == 1)
# 74109 for BIARef

##################################################
# PRIVATE SV
##################################################

# 74109 SV for BIARef
SV_private %>% count(SVType, sort = TRUE)

SV_private <- SV_private %>%
  left_join(
    BIA_Ref_chromosome_Link_file,
    by = c("Chrom_original" = "BIA_Scaff")
  )


# Check species distribution
sv_counts <- SV_private %>%
  count(SVType, Species)

species_total <- SV_private %>%
  count(Species, name = "n") %>%
  mutate(SVType = "TOTAL")

sv_counts_all <- bind_rows(sv_counts, species_total)

sv_counts_all$SVType <- factor(sv_counts_all$SVType, 
                               labels=c("TOTAL", "INS", "DEL",
                                        "TRANS", "DUP", "INV"))


p_privateSV_inSpecies <- 
  ggplot(sv_counts_all, aes(x = Species, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ SVType, scales = "free_x") +
  #facet_wrap(~ SVType) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "SV counts per species by SV type - Private SV",
    x = "Species",
    y = "Count"
  )


# Plot in file
name_plot1 <- paste0(out_prefix, Reference_Species, ".PrivateSV.BySpeciesandType.pdf")
pdf(name_plot1, height = 6, width = 7)
p_privateSV_inSpecies
dev.off()


# Plot with stacked
sv_counts_plot <- sv_counts %>%
  group_by(Species) %>%
  mutate(total = sum(n)) %>%
  ungroup()

sv_counts_plot$SVType <- factor(
  sv_counts_plot$SVType,
  levels = c("DEL", "INS", "DUP", "INV", "TRANS")
)

# Stacked barplot

p_privateSV_inSpecies_stacked <- 
  ggplot(sv_counts_plot, aes(x = Species, y = n, fill = SVType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sv_type_colours) +
  theme_bw() +
  labs(
    x = "Species",
    y = "Number of SVs",
    fill = "SV type",
    title = "Private SV | SV counts per species by SV type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_flip()


# Plot in file
name_plot1b <- paste0(out_prefix, Reference_Species, ".PrivateSV.BySpeciesandType.Stacked.pdf")
pdf(name_plot1b, height = 6, width = 7)
p_privateSV_inSpecies_stacked
dev.off()

# Save y-max for plotting later
ylim_for_SVNumer <- max(sv_counts_plot$total)


# For BIA: Most of the private SV: OCE, LAZ, PRC, and CLA

##############
# For Private SV: see if the affect some genes.
# See how many in regions were we have genes
#############

sv_gr <- GRanges(
  seqnames = SV_private$BIA_Chr,
  ranges = IRanges(
    start = SV_private$Position,
    end   = SV_private$Position + SV_private$SVLength
  ),
  SVType = SV_private$SVType,
  Species = SV_private$Species
)

# Overlapping to coding regions
hits_cds <- findOverlaps(sv_gr, cds, ignore.strand=T, type="any")
SV_private$in_coding <- FALSE
SV_private$in_coding[unique(queryHits(hits_cds))] <- TRUE
nrow((SV_private[SV_private$in_coding == "TRUE",]))
# 6452 in Coding regions for BIA

# Overlapping to Genes
hits_gene <- findOverlaps(sv_gr, genes, ignore.strand=T, type="any")
SV_private$in_gene <- FALSE
SV_private$in_gene[unique(queryHits(hits_gene))] <- TRUE
nrow((SV_private[SV_private$in_gene == "TRUE",]))
# 45263 in Genes for BIA  --> We only consider in Coding Regions

# Focus on the one that are affecting coding regions
affected_cds <- split(
  mcols(cds)$Parent[subjectHits(hits_cds)],
  queryHits(hits_cds)
)

SV_private$Affected_CDS <- NA
SV_private$Affected_CDS[as.integer(names(affected_cds))] <- sapply(
  affected_cds,
  function(genes) {
    paste(sort(unique(genes)), collapse = ",")
  }
)

##########
### SV Affecting only CDS - PLOTS and STATS
##########

SV_private_AffectingCDS <- SV_private %>%
  filter(in_coding == TRUE)
# 6452 SV in BIA

# Filter if CDS rported more than once
SV_private_AffectingCDS$Affected_CDS <- lapply(
  SV_private_AffectingCDS$Affected_CDS,
  function(x) unique(x)
)

SV_private_AffectingCDS$Affected_CDS <- sapply(
  SV_private_AffectingCDS$Affected_CDS,
  function(x) {
    if (length(x) == 0 || all(is.na(x))) return(NA)
    paste(unique(x), collapse = ",")
  },
  USE.NAMES = FALSE
)


sv_counts_inCDS <- SV_private_AffectingCDS %>%
  count(SVType, Species)

species_total_inCDS <- SV_private_AffectingCDS %>%
  count(Species, name = "n") %>%
  mutate(SVType = "TOTAL")

sv_counts_inCDS <- bind_rows(sv_counts_inCDS, species_total_inCDS)

sv_counts_inCDS$SVType <- factor(sv_counts_inCDS$SVType, 
                                 levels = c("TOTAL", "INS", "DEL",
                                            "TRANS", "DUP", "INV"))

p_privateSV_inCDS_inSpecies <- ggplot(sv_counts_inCDS, aes(x = Species, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ SVType, scales = "free_x") +
  #facet_wrap(~ SVType) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "SV counts per species by SV type - Private SV, affecting CDS",
    x = "Species",
    y = "Count"
  )

# PRIVATE SV AFFECTING CDS - STACKED PLOT 

sv_counts_inCDS <- SV_private_AffectingCDS %>%
  count(SVType, Species)

sv_counts_inCDS_plot <- sv_counts_inCDS %>%
  group_by(Species) %>%
  mutate(total = sum(n)) %>%
  ungroup()

sv_counts_inCDS_plot$SVType <- factor(
  sv_counts_inCDS_plot$SVType,
  levels = c("DEL", "INS", "DUP", "INV", "TRANS")
)

# Stacked barplot

p_privateSV_inSpecies_inCDS_stacked <- 
  ggplot(sv_counts_inCDS_plot, aes(x = Species, y = n, fill = SVType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sv_type_colours) +
  #coord_flip(ylim = c(0, ylim_for_SVNumer)) +
  theme_bw() +
  labs(
    x = "Species",
    y = "Number of SVs",
    fill = "SV type",
    title = "Private SV - Within CDS | SV counts per species by SV type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_flip()

p_privateSV_inSpecies_inCDS_stacked



##########
### Save Results in Table
##########

table1_name <- paste0(out_prefix,Reference_Species, ".PrivateSV.AffectingCDS.csv")

write.table(SV_private_AffectingCDS, table1_name,
            sep="\t", quote = F, row.names = F, col.names = T)


# Also save everything if needed
# Filter if CDS rported more than once
SV_private$Affected_CDS <- lapply(
  SV_private$Affected_CDS,
  function(x) unique(x)
)

SV_private$Affected_CDS <- sapply(
  SV_private$Affected_CDS,
  function(x) {
    if (length(x) == 0 || all(is.na(x))) return(NA)
    paste(unique(x), collapse = ",")
  },
  USE.NAMES = FALSE
)


table1b_name <- paste0(out_prefix,Reference_Species, ".PrivateSV.AllSV_WithCDSInfo.csv")

write.table(SV_private, table1b_name,
            sep="\t", quote = F, row.names = F, col.names = T)





########################################################
# SHARED SV
########################################################

SV_shared <- SV_all %>% filter(n_species >= 2)
# 54572 in BIA

SV_shared <- SV_shared %>%
  left_join(
    BIA_Ref_chromosome_Link_file,
    by = c("Chrom_original" = "BIA_Scaff")
  )

# Overall distribution: does it reflect the overall?

sv_counts_shared <- SV_shared %>%
  count(SVType)

### Check if affecting Coding regions

sv_gr_shared <- GRanges(
  seqnames = SV_shared$BIA_Chr,
  ranges = IRanges(
    start = SV_shared$Position,
    end   = SV_shared$Position + SV_shared$SVLength
  ),
  SVType = SV_shared$SVType,
  Species = SV_shared$Species
)

# Overlapping to coding regions
hits_cds_shared <- findOverlaps(sv_gr_shared, cds, ignore.strand=T, type="any")
SV_shared$in_coding <- FALSE
SV_shared$in_coding[unique(queryHits(hits_cds_shared))] <- TRUE
nrow((SV_shared[SV_shared$in_coding == "TRUE",]))
# 3'625 in Coding regions for BIA

# Overlapping to Genes
hits_gene_shared <- findOverlaps(sv_gr_shared, genes, ignore.strand=T, type="any")
SV_shared$in_gene <- FALSE
SV_shared$in_gene[unique(queryHits(hits_gene_shared))] <- TRUE
nrow((SV_shared[SV_shared$in_gene == "TRUE",]))
# 33'441 in Genes in BIA --> We only consider in Coding Regions

# Focus on the one that are affecting coding regions
affected_cds_shared <- split(
  mcols(cds)$Parent[subjectHits(hits_cds_shared)],
  queryHits(hits_cds_shared)
)

SV_shared$Affected_CDS <- NA
SV_shared$Affected_CDS[as.integer(names(affected_cds_shared))] <- sapply(
  affected_cds_shared,
  function(genes) {
    paste(sort(unique(genes)), collapse = ",")
  }
)

SV_shared_AffectingCDS <- SV_shared %>%
  filter(in_coding == TRUE)
# 3625 SV in BIA

# Filter if CDS rported more than once
SV_shared_AffectingCDS$Affected_CDS <- lapply(
  SV_shared_AffectingCDS$Affected_CDS,
  function(x) unique(x))

SV_shared_AffectingCDS$Affected_CDS <- sapply(
  SV_shared_AffectingCDS$Affected_CDS,
  function(x) {
    if (length(x) == 0 || all(is.na(x))) return(NA)
    paste(unique(x), collapse = ",")
  },
  USE.NAMES = FALSE
)



##########
### Save Results in Table
##########

table2_name <- paste0(out_prefix,Reference_Species, ".SharedSV.AffectingCDS.csv")

write.table(SV_shared_AffectingCDS, table2_name,
            sep="\t", quote = F, row.names = F, col.names = T)


# Also save everthin if needed
# Filter if CDS rported more than once

SV_shared$Affected_CDS <- lapply(
  SV_shared$Affected_CDS,
  function(x) unique(x)
)

SV_shared$Affected_CDS <- sapply(
  SV_shared$Affected_CDS,
  function(x) {
    if (length(x) == 0 || all(is.na(x))) return(NA)
    paste(unique(x), collapse = ",")
  },
  USE.NAMES = FALSE
)


table2b_name <- paste0(out_prefix,Reference_Species, ".SharedSV.AllSV_WithCDSInfo.csv")

write.table(SV_shared, table2b_name,
            sep="\t", quote = F, row.names = F, col.names = T)





##############################
#### Check if we can do some plots of the species
##############################

SV_shared$Species <- sapply(strsplit(SV_shared$Species, ","), function(x) {
  paste(sort(x), collapse = ",")
})

sv_counts_shared <- SV_shared %>%
  count(SVType, Species)


# Plot with stacked
sv_counts_shared_plot <- sv_counts_shared %>%
  group_by(Species) %>%
  mutate(total = sum(n)) %>%
  ungroup()

sv_counts_plot$SVType <- factor(
  sv_counts_plot$SVType,
  levels = c("DEL", "INS", "DUP", "INV", "TRANS")
)

# Stacked barplot

p_privateSV_inSpecies_stacked <- 
  ggplot(sv_counts_plot, aes(x = Species, y = n, fill = SVType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sv_type_colours) +
  theme_bw() +
  labs(
    x = "Species",
    y = "Number of SVs",
    fill = "SV type",
    title = "Private SV | SV counts per species by SV type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_flip()





###########
#### Check length of SV not affecting SV
###########
SV_shared_NoCDS <- SV_shared %>%
  filter_out(in_coding == TRUE)


SV_private_noCDS <- SV_private %>%
  filter_out(in_coding == TRUE)


all_lengths <- c(
  SV_shared_AffectingCDS$SVLength,
  SV_shared_NoCDS$SVLength,
  SV_private_AffectingCDS$SVLength,
  SV_private_noCDS$SVLength
)

ylim <- max(log10(all_lengths), na.rm = TRUE)

p_sharedCDS <- ggplot(SV_shared_AffectingCDS, aes(x = SVType, y = SVLength, fill = SVType)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.6) +
  scale_y_log10() +
  scale_fill_manual(values = sv_type_colours) +
  theme_bw(base_size = 12) +
  labs(
    title = "SV Length - Shared SV, Affecting CDS",
    x = "SV Type",
    y = "SV Length (log10)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) + coord_cartesian(ylim = c(1, 10^(ylim+1)))


p_shared_noCDS <- ggplot(SV_shared_NoCDS, aes(x = SVType, y = SVLength, fill = SVType)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.6) +
  scale_y_log10() +
  scale_fill_manual(values = sv_type_colours) +
  theme_bw(base_size = 12) +
  labs(
    title = "SV Length - Shared SV, Not affecting CDS",
    x = "SV Type",
    y = "SV Length (log10)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) + coord_cartesian(ylim = c(1, 10^(ylim+1)))


p_privateCDS <- ggplot(SV_private_AffectingCDS, aes(x = SVType, y = SVLength, fill = SVType)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.6) +
  scale_y_log10() +
  scale_fill_manual(values = sv_type_colours) +
  theme_bw(base_size = 12) +
  labs(
    title = "SV Length - Private SV, Affecting CDS",
    x = "SV Type",
    y = "SV Length (log10)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) + coord_cartesian(ylim = c(1, 10^(ylim+1)))



p_private_noCDS <- ggplot(SV_private_noCDS, aes(x = SVType, y = SVLength, fill = SVType)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.6) +
  scale_y_log10() +
  scale_fill_manual(values = sv_type_colours) +
  theme_bw(base_size = 12) +
  labs(
    title = "SV Length - Private SV, Not affecting CDS",
    x = "SV Type",
    y = "SV Length (log10)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) + coord_cartesian(ylim = c(1, 10^(ylim+1)))


title_plot3 <- paste0(out_prefix, Reference_Species, ".Length_by_Category.pdf")
pdf(title_plot3, width = 10, height = 10)
(p_privateCDS + p_private_noCDS) / (p_sharedCDS + p_shared_noCDS) 
dev.off()


##########
### PLOT ALL INFO TOGETHER
##########


##########
### PRIVATE SV
##########

all_lengths_private <- c(
  SV_private_AffectingCDS$SVLength,
  SV_private$SVLength
)

ylim_private <- max(log10(all_lengths_private), na.rm = TRUE)

SV_private_AffectingCDS$SVType <- factor(SV_private_AffectingCDS$SVType,
                                         levels=c("DEL", "INS", "DUP", "TRANS", "INV"))

p_privateCDS_Length <- ggplot(SV_private_AffectingCDS, aes(x = SVType, y = SVLength, fill = SVType)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.6) +
  scale_y_log10() +
  scale_fill_manual(values = sv_type_colours) +
  theme_bw(base_size = 12) +
  labs(
    title = "Private SV - Within CDS | SV Length Distribution by SV Type",
    x = "SV Type",
    y = "SV Length (log10)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) + coord_cartesian(ylim = c(1, 10^(ylim_private+1)))

SV_private$SVType <- factor(SV_private$SVType,
                            levels=c("DEL", "INS", "DUP", "TRANS", "INV"))

p_private_allSV_Length <- ggplot(SV_private, aes(x = SVType, y = SVLength, fill = SVType)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.6) +
  scale_y_log10() +
  scale_fill_manual(values = sv_type_colours) +
  theme_bw(base_size = 12) +
  labs(
    title = "Private SV | SV Length Distribution by SV Type",
    x = "SV Type",
    y = "SV Length (log10)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) + coord_cartesian(ylim = c(1, 10^(ylim_private+1)))


# Save the plots

A <- p_privateSV_inSpecies_stacked + labs(title = NULL)+ 
  p_private_allSV_Length + labs(title = NULL) +
  plot_annotation(title = "A. Private SVs | SV counts per species and SV Length distribution") 

B <- p_privateSV_inSpecies_inCDS_stacked + labs(title = NULL)+
  p_privateCDS_Length + labs(title = NULL) +
  plot_annotation(title = "B. Private SVs in CDS | SV counts per species and SV Length distribution")

p_name_PrivateSV_allInfo <- paste0(out_prefix, Reference_Species, ".PrivateSV.CountsAndLength.pdf")
pdf(p_name_PrivateSV_allInfo, width = 10, height = 9)
plot_grid(A, B, ncol=1)
dev.off()  


#########

######  SHARED SV IN CDS
# Check if something specific in species 

head(SV_shared_AffectingCDS)

SV_shared_AffectingCDS$Species <- sapply(
  SV_shared_AffectingCDS$Species,
  function(x) {
    paste(sort(trimws(strsplit(x, ",")[[1]])), collapse = ",")
  }
)

normalize_complex <- function(x) {
  paste(sort(trimws(strsplit(x, ",")[[1]])), collapse = ",")
}

SV_shared_AffectingCDS$Species_complex <- sapply(
  SV_shared_AffectingCDS$Species,
  normalize_complex
)

complex_counts <- sort(table(SV_shared_AffectingCDS$Species_complex), decreasing = TRUE)

complex_counts





###### ALL SHARED SV
# Check if something specific in species 

SV_shared$Species <- sapply(
  SV_shared$Species,
  function(x) {
    paste(sort(trimws(strsplit(x, ",")[[1]])), collapse = ",")
  }
)

normalize_complex <- function(x) {
  paste(sort(trimws(strsplit(x, ",")[[1]])), collapse = ",")
}

SV_shared$Species_complex <- sapply(
  SV_shared$Species,
  normalize_complex
)

complex_counts <- sort(table(SV_shared$Species_complex), decreasing = TRUE)

complex_counts