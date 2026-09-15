# ============================================================
# Pairwise Robinson-Foulds distances + MDS + heatmap
# across chromosome-level AND inversion-genotype tree topologies
#
# Handles trees with unequal species sets (e.g. inversion-genotype
# trees restricted to a subset of species) via pairwise dynamic
# tip-pruning + normalized RF, rather than forcing every tree down
# to one global common tip set.
# ============================================================

library(ape)
library(phytools)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)
library(viridis)

# ---------------------------------------------------------------
# Parameters -- edit these
# ---------------------------------------------------------------
tree_dir         <- "Chrom_and_INV_Topologies/treefiles/"
tree_pattern     <- "\\.treefile$"
wholegenome_file <- "Partitions.Concat_1to1OG.ClownDams.txt.treefile"
chrom_map_file   <- "CLA_vs_BIARef.mapids.txt"   # set to NULL to skip accession->chr renaming
exclude_patterns <- c()                    # trees dropped entirely (partial name match), e.g. poor SyRI alignment
min_shared_tips  <- 4                             # minimum shared tips required to compare a pair (need >=4 for >=1 internal edge)
normalize_rf     <- TRUE                          # normalize RF by max possible RF for the shared tip count
out_dir          <- "Chrom_and_INV_Topologies/Plot_topo_Inconsistency_Chrom/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# Helper: map CM accession -> chr name, preserving any suffix
# (e.g. "CM132566.1_INV" -> "chr18_INV"). Falls back to the raw
# file basename if no accession match / no chrom_map supplied.
# ---------------------------------------------------------------
load_chrom_map <- function(mapids_file) {
  if (is.null(mapids_file) || !file.exists(mapids_file)) return(NULL)
  read.table(mapids_file, header = FALSE, sep = "\t",
             col.names = c("chr_name", "accession"))
}
chrom_map <- load_chrom_map(chrom_map_file)

extract_name <- function(filename, chrom_map = NULL) {
  base <- sub(tree_pattern, "", basename(filename))
  if (!is.null(chrom_map)) {
    acc <- regmatches(base, regexpr("CM[0-9]+\\.[0-9]+", base))
    if (length(acc) == 1 && nzchar(acc)) {
      hit <- chrom_map[chrom_map$accession == acc, ]
      if (nrow(hit) > 0) {
        suffix <- sub(acc, "", base, fixed = TRUE)  # keeps _INV / _STD / etc.
        return(paste0(hit$chr_name[1], suffix))
      }
    }
  }
  base
}

# ---------------------------------------------------------------
# Read all chromosome / inversion-genotype trees
# ---------------------------------------------------------------
files <- list.files(tree_dir, pattern = tree_pattern, full.names = TRUE, recursive = TRUE)
files <- files[normalizePath(files) != normalizePath(wholegenome_file, mustWork = FALSE)]
stopifnot(length(files) > 0)

tree_names <- vapply(files, extract_name, character(1), chrom_map = chrom_map)
if (any(duplicated(tree_names))) {
  stop("Duplicate tree names after extraction -- check extract_name() logic:\n  ",
       paste(tree_names[duplicated(tree_names)], collapse = ", "))
}
trees_chr <- setNames(lapply(files, read.tree), tree_names)

# Drop explicitly excluded trees (e.g. chr17, poor SyRI alignment)
if (length(exclude_patterns) > 0) {
  keep <- !Reduce(`|`, lapply(exclude_patterns, function(p) grepl(p, names(trees_chr))))
  message(sprintf("Excluding %d tree(s) matching: %s",
                  sum(!keep), paste(exclude_patterns, collapse = ", ")))
  trees_chr <- trees_chr[keep]
}

# Add whole-genome tree
tree_wg <- read.tree(wholegenome_file)
trees   <- c(trees_chr, list(WholeGenome = tree_wg))
trees   <- lapply(trees, unroot)

# ---------------------------------------------------------------
# Species coverage metadata (kept in lockstep with `trees` throughout)
# ---------------------------------------------------------------
tree_meta <- data.frame(
  ShortName = names(trees),
  Ntips     = sapply(trees, function(t) length(t$tip.label)),
  stringsAsFactors = FALSE
)
tree_meta$Coverage <- ifelse(tree_meta$Ntips == max(tree_meta$Ntips), "Full", "Reduced")
tree_meta$Type     <- ifelse(tree_meta$ShortName == "WholeGenome", "WholeGenome", "Chromosome/Inversion")
tree_meta$Label    <- paste0(tree_meta$ShortName, " (n=", tree_meta$Ntips, ")")
rownames(tree_meta) <- tree_meta$Label
names(trees)        <- tree_meta$Label

cat("\n--- Tree species coverage ---\n")
print(tree_meta[, c("ShortName", "Ntips", "Coverage", "Type")], row.names = FALSE)

# ---------------------------------------------------------------
# Keep only full-species-coverage trees (drop reduced-taxa windows).
# "Full" is defined relative to the other chromosome/window trees,
# not the whole-genome tree, since the whole-genome tree can
# legitimately have a slightly different (usually smaller) tip set
# -- it is always kept regardless, as the comparison anchor.
# ---------------------------------------------------------------
only_full_species <- TRUE  # set FALSE to restore the dynamic-pruning behaviour for all trees

if (only_full_species) {
  full_n <- max(tree_meta$Ntips[tree_meta$Type != "WholeGenome"])
  keep_full <- (tree_meta$Ntips == full_n) | (tree_meta$Type == "WholeGenome")
  n_dropped <- sum(!keep_full)
  if (n_dropped > 0) {
    message(sprintf("\nDropping %d reduced-coverage tree(s) (< %d taxa): %s",
                    n_dropped, full_n, paste(tree_meta$ShortName[!keep_full], collapse = ", ")))
  }
  trees     <- trees[keep_full]
  tree_meta <- tree_meta[keep_full, ]
  cat(sprintf("\nKept %d tree(s) with full species coverage (%d taxa), plus the whole-genome tree (%d taxa).\n",
              sum(tree_meta$Type != "WholeGenome"), full_n,
              tree_meta$Ntips[tree_meta$Type == "WholeGenome"]))
}

# ---------------------------------------------------------------
# Pairwise RF distance with DYNAMIC per-pair tip pruning
# (each pair is pruned only to the tips *those two* trees share,
#  not to one global common set -- preserves resolution for
#  full-species comparisons while still allowing reduced trees in)
# ---------------------------------------------------------------
compute_rf_matrix <- function(trees, min_shared_tips = 4, normalize = TRUE) {
  nm <- names(trees); n <- length(trees)
  rf_mat     <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
  shared_mat <- matrix(NA_integer_, n, n, dimnames = list(nm, nm))
  diag(rf_mat) <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      ti <- trees[[i]]; tj <- trees[[j]]
      shared <- intersect(ti$tip.label, tj$tip.label)
      shared_mat[i, j] <- shared_mat[j, i] <- length(shared)
      
      if (length(shared) < min_shared_tips) {
        message(sprintf("  [skip] %s vs %s: only %d shared tips (< %d)",
                        nm[i], nm[j], length(shared), min_shared_tips))
        next
      }
      
      pi <- drop.tip(ti, setdiff(ti$tip.label, shared))
      pj <- drop.tip(tj, setdiff(tj$tip.label, shared))
      # drop.tip() can leave a pruned subtree with a bifurcating root
      # (ape then treats it as "rooted" even though it isn't) --
      # re-unroot before comparing to avoid the spurious dist.topo() warning
      pi <- unroot(pi)
      pj <- unroot(pj)
      rf_raw <- as.numeric(dist.topo(pi, pj))
      
      if (normalize) {
        max_rf <- 2 * (length(shared) - 3)  # max RF for an unrooted binary tree
        rf_raw <- if (max_rf > 0) rf_raw / max_rf else 0
      }
      rf_mat[i, j] <- rf_mat[j, i] <- rf_raw
    }
  }
  list(rf = rf_mat, shared = shared_mat)
}

cat("\n--- Computing pairwise RF distances (dynamic pruning) ---\n")
rf_out     <- compute_rf_matrix(trees, min_shared_tips = min_shared_tips, normalize = normalize_rf)
rf_matrix  <- rf_out$rf
shared_mat <- rf_out$shared

# cmdscale() requires a complete matrix -- drop any tree left with
# unresolved (NA) comparisons to *every* other tree it couldn't be compared to
if (anyNA(rf_matrix)) {
  still_na <- rownames(rf_matrix)[apply(rf_matrix, 1, function(r) anyNA(r))]
  if (length(still_na) > 0) {
    message("\nDropping from heatmap/MDS (insufficient shared tips with other trees): ",
            paste(still_na, collapse = ", "))
    keep <- !(rownames(rf_matrix) %in% still_na)
    rf_matrix <- rf_matrix[keep, keep]
    trees     <- trees[keep]
    tree_meta <- tree_meta[keep, ]
  }
}

# ---------------------------------------------------------------
# Count distinct topologies among the compared trees
# (trees with RF = 0 relative to each other -- i.e. identical on
#  their shared taxa -- are grouped as the same topology).
# Complete linkage + cutree(h = tol) only merges a group if EVERY
# pairwise distance within it is 0, avoiding chained false-merges.
# ---------------------------------------------------------------
count_topologies <- function(rf_matrix, tol = 1e-8) {
  d  <- as.dist(rf_matrix)
  hc <- hclust(d, method = "complete")
  groups <- cutree(hc, h = tol)
  data.frame(Label = rownames(rf_matrix), TopologyGroup = groups, row.names = NULL)
}

topo_groups <- count_topologies(rf_matrix)
n_unique_topologies <- length(unique(topo_groups$TopologyGroup))

cat(sprintf("\n--- Distinct topologies ---\nNumber of distinct topologies among %d compared trees: %d\n",
            nrow(topo_groups), n_unique_topologies))

topo_summary <- topo_groups %>%
  group_by(TopologyGroup) %>%
  summarise(N_trees = n(), Trees = paste(Label, collapse = ", "), .groups = "drop") %>%
  arrange(desc(N_trees))
print(topo_summary, n = Inf)

write.csv(topo_summary, file.path(out_dir, "Topology_groups.ChromAndINV.csv"), row.names = FALSE)

# ---------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------
hc <- hclust(as.dist(rf_matrix))
order <- hc$order
rf_matrix_ord <- rf_matrix[order, order]

rf_df_ord <- melt(rf_matrix_ord, varnames = c("Tree1", "Tree2"), value.name = "RF")
rf_df_ord <- rf_df_ord %>%
  mutate(
    Tree1 = factor(Tree1, levels = rownames(rf_matrix_ord)),
    Tree2 = factor(Tree2, levels = colnames(rf_matrix_ord)),
    IsWholeGenome = grepl("^WholeGenome ", Tree1) | grepl("^WholeGenome ", Tree2)
  )

p_heatmap <- ggplot(rf_df_ord, aes(Tree1, Tree2, fill = RF)) +
  geom_tile(color = "white", size = 0.2) +
  geom_tile(data = rf_df_ord %>% filter(IsWholeGenome),
            aes(Tree1, Tree2), color = "black", size = 0.5, fill = NA) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1,
                       name = if (normalize_rf) "Normalized RF" else "RF distance") +
  coord_fixed() +
  labs(
    title = "Pairwise Robinson-Foulds Distances",
    subtitle = "Chromosome + inversion-genotype topologies vs. whole genome",
    x = "", y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "grey20"),
    axis.text.y = element_text(size = 9, color = "grey20"),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey40")
  )

#pdf(file.path(out_dir, "Heatmap_Trees.ChromAndINV.pdf"), width = 9, height = 9)
pdf(file.path(out_dir, "Heatmap_Trees.ChromAndINV.AllSpecies.pdf"), width = 9, height = 9)
print(p_heatmap)
dev.off()

# ---------------------------------------------------------------
# MDS
# ---------------------------------------------------------------
mds <- cmdscale(as.dist(rf_matrix), k = 2)
mds_df <- data.frame(MDS1 = mds[, 1], MDS2 = mds[, 2], Label = rownames(mds))
mds_df <- merge(mds_df, tree_meta, by = "Label")

p_MDS <- ggplot(mds_df, aes(MDS1, MDS2, color = Type, shape = Coverage, label = ShortName)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(size = 4, fontface = "bold", max.overlaps = 30, show.legend = FALSE) +
  scale_color_manual(values = c("Chromosome/Inversion" = "steelblue", "WholeGenome" = "firebrick")) +
  scale_shape_manual(values = c("Full" = 16, "Reduced" = 17)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_line(color = "grey80"),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    title = "MDS of Phylogenetic Tree Topologies",
    subtitle = if (normalize_rf) "Distances = normalized Robinson-Foulds (RF)" else "Distances = Robinson-Foulds (RF)",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2"
  )

pdf(file.path(out_dir, "MDS_Trees.ChromAndINV.pdf"), width = 8, height = 8)
print(p_MDS)
dev.off()

p_heatmap
p_MDS
