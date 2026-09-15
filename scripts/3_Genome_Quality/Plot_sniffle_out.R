library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/Users/amarcion/Documents/PacBio_Sequencing_Assembly/10_Additional/15_Long_Reads_Mapping/sniffle_out")

# Read summaries
summary <- read.csv("SV_summary_by_genotype.csv")
raw     <- read.csv("all_SVs_raw.csv")

# Define colors for genotypes
gt_colors <- c(
  "0/0" = "#2196F3",  # blue  = hom ref
  "0/1" = "#FF9800",  # orange = het
  "1/1" = "#F44336"   # red   = hom alt
)

# SV type order
sv_order <- c("DEL", "INS", "DUP", "INV", "TRANS")

# Plot 1: Count per SV type and genotype per species
p1 <- ggplot(summary, 
             aes(x = factor(svtype, sv_order), 
                 y = count, 
                 fill = genotype)) +
  geom_bar(stat = "identity", 
           position = "dodge") +
  scale_fill_manual(values = gt_colors) +
  facet_wrap(~species, ncol = 4) +
  theme_bw() +
  labs(title = "SV counts by type and genotype",
       subtitle = "Self-mapped reads",
       x = "SV Type", y = "Count",
       fill = "Genotype") +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1))

ggsave("SV_counts_by_genotype.pdf", p1, 
       width = 16, height = 12)

# Plot 2: Total length per SV type and genotype
p2 <- ggplot(summary,
             aes(x = factor(svtype, sv_order),
                 y = total_len_mb,
                 fill = genotype)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values = gt_colors) +
  facet_wrap(~species, ncol = 4) +
  theme_bw() +
  labs(title = "Total SV length by type and genotype",
       subtitle = "Self-mapped reads",
       x = "SV Type", y = "Total length (Mb)",
       fill = "Genotype") +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1))

ggsave("SV_length_by_genotype.pdf", p2,
       width = 16, height = 12)

# Plot 3: Focus on heterozygous SVs
# These are the most interesting
het <- summary %>% filter(genotype == "0/1")

p3 <- ggplot(het,
             aes(x = species,
                 y = count,
                 fill = factor(svtype, sv_order))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "Heterozygous SVs per species",
       subtitle = "Potentially missed by purge_dups",
       x = "Species", y = "Count",
       fill = "SV Type") +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1))

ggsave("heterozygous_SVs.pdf", p3,
       width = 10, height = 6)

# Plot 4: Proportion of genotypes per SV type
prop <- summary %>%
  group_by(svtype, genotype) %>%
  summarise(total_count = sum(count),
            total_len   = sum(total_len_mb)) %>%
  group_by(svtype) %>%
  mutate(prop_count = total_count/sum(total_count),
         prop_len   = total_len/sum(total_len))

p4 <- ggplot(prop,
             aes(x = factor(svtype, sv_order),
                 y = prop_count,
                 fill = genotype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = gt_colors) +
  theme_bw() +
  labs(title = "Genotype proportions per SV type",
       x = "SV Type",
       y = "Proportion",
       fill = "Genotype")

ggsave("genotype_proportions.pdf", p4,
       width = 8, height = 5)
