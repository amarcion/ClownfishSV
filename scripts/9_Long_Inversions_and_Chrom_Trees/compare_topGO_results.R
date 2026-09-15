base    <- "GO_Enrichment"
res_dir <- file.path(base, "Results_perSpecies")

chr18_species <- c("AKA","AKY","ALL","CRP","LAT","MCC","OMA","POL","SAN","SEB")
chr9_species  <- c("POL","CLA")

args <- commandArgs(trailingOnly = TRUE)
PVAL_CUT <- if (length(args) >= 1) as.numeric(args[1]) else 0.05
pval_tag <- sub("^0\\.", "p", format(PVAL_CUT, scientific = FALSE))

compare_inversion <- function(species_vec, inv, node_size) {
  tabs <- lapply(species_vec, function(sp) {
    f <- file.path(res_dir, paste0(sp, "_", inv, "_nodeSize", node_size, ".txt"))
    d <- read.delim(f, stringsAsFactors = FALSE)
    d[d$classicFisher_num < PVAL_CUT & !is.na(d$classicFisher_num), c("GO.ID","Term","classicFisher_num")]
  })
  names(tabs) <- species_vec

  all_ids <- unique(unlist(lapply(tabs, function(d) d$GO.ID)))
  term_lookup <- c()
  for (d in tabs) {
    for (i in seq_len(nrow(d))) term_lookup[d$GO.ID[i]] <- d$Term[i]
  }

  mat <- matrix(NA_real_, nrow = length(all_ids), ncol = length(species_vec),
                dimnames = list(all_ids, species_vec))
  for (sp in species_vec) {
    d <- tabs[[sp]]
    mat[d$GO.ID, sp] <- d$classicFisher_num
  }

  nSig <- rowSums(!is.na(mat))
  out <- data.frame(
    GO.ID = all_ids,
    Term = term_lookup[all_ids],
    N_species_significant = nSig,
    N_species_tested = length(species_vec),
    stringsAsFactors = FALSE
  )
  out <- cbind(out, mat)
  out <- out[order(-out$N_species_significant, out$GO.ID), ]
  rownames(out) <- NULL
  out
}

for (ns in c(5, 10)) {
  chr18_cmp <- compare_inversion(chr18_species, "chr18", ns)
  write.table(chr18_cmp, file = file.path(res_dir, paste0("Comparison_chr18_nodeSize", ns, "_", pval_tag, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  chr9_cmp <- compare_inversion(chr9_species, "chr9", ns)
  write.table(chr9_cmp, file = file.path(res_dir, paste0("Comparison_chr9_nodeSize", ns, "_", pval_tag, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  cat(sprintf("\n=== nodeSize=%d (p<%s) ===\n", ns, PVAL_CUT))
  cat(sprintf("chr18: %d unique significant GO terms across %d species\n", nrow(chr18_cmp), length(chr18_species)))
  cat("  Terms significant in >=5/10 species:\n")
  top18 <- chr18_cmp[chr18_cmp$N_species_significant >= 5, c("GO.ID","Term","N_species_significant")]
  print(top18, row.names = FALSE)

  cat(sprintf("\nchr9: %d unique significant GO terms across %d species\n", nrow(chr9_cmp), length(chr9_species)))
  cat("  Terms significant in BOTH species:\n")
  top9 <- chr9_cmp[chr9_cmp$N_species_significant == 2, c("GO.ID","Term","N_species_significant")]
  print(top9, row.names = FALSE)
}
