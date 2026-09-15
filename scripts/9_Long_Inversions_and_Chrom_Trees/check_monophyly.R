library(ape)

inv_sp  <- read.table("INV_and_Species.txt", header=TRUE, sep="\t",
                      stringsAsFactors=FALSE, quote="\"")

tree_dir <- "INV_SNPs_Tree/"

results <- data.frame(
  Accession      = character(),
  Chunk          = integer(),
  N_carriers     = integer(),
  Carriers       = character(),
  N_tips_tree    = integer(),
  Tips_in_tree   = character(),
  Carriers_present = integer(),
  Is_monophyletic  = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(inv_sp))) {
  acc   <- inv_sp$Accession[i]
  chunk <- inv_sp$Chunck_Nb[i]

  # Parse carriers — fix "ALL.LAT" typo and deduplicate
  raw_sp  <- gsub("\\s", "", inv_sp$Species_INV[i])
  raw_sp  <- gsub("\\.", ",", raw_sp)          # fix dot separator
  carriers <- unique(trimws(strsplit(raw_sp, ",")[[1]]))
  carriers <- carriers[carriers != ""]

  # Find matching treefile
  pattern  <- paste0(acc, ".Chunk_", chunk, ".INV")
  files    <- list.files(tree_dir, pattern=pattern, full.names=TRUE)
  files    <- files[grepl("\\.treefile$", files)]

  if (length(files) == 0) {
    cat(sprintf("  NO TREE: %s Chunk %d\n", acc, chunk))
    results <- rbind(results, data.frame(
      Accession=acc, Chunk=chunk,
      N_carriers=length(carriers),
      Carriers=paste(sort(carriers), collapse=";"),
      N_tips_tree=NA, Tips_in_tree=NA,
      Carriers_present=NA, Is_monophyletic="NO_TREE",
      stringsAsFactors=FALSE))
    next
  }

  tree <- read.tree(files[1])
  tips <- tree$tip.label

  # Carriers present in this tree
  present <- intersect(carriers, tips)

  if (length(present) < 2) {
    mono_result <- "TOO_FEW"
  } else if (length(present) == length(tips)) {
    mono_result <- "ALL_TIPS"        # trivially monophyletic
  } else {
    mono_result <- ifelse(is.monophyletic(tree, present), "TRUE", "FALSE")
  }

  cat(sprintf("  %s Ck%d  carriers=%d  in_tree=%d/%d  monophyletic=%s\n",
              acc, chunk, length(carriers), length(present),
              length(tips), mono_result))

  results <- rbind(results, data.frame(
    Accession        = acc,
    Chunk            = chunk,
    N_carriers       = length(carriers),
    Carriers         = paste(sort(carriers), collapse=";"),
    N_tips_tree      = length(tips),
    Tips_in_tree     = paste(sort(tips), collapse=";"),
    Carriers_present = length(present),
    Is_monophyletic  = mono_result,
    stringsAsFactors = FALSE
  ))
}

out <- "Plotted_Trees/Monophyly_INV_Carriers.csv"
write.csv(results, out, row.names=FALSE, quote=TRUE)
cat(sprintf("\nResults written to %s\n", out))

# Summary
cat("\n=== Summary ===\n")
tab <- table(results$Is_monophyletic)
print(tab)

cat("\nMonophyletic (TRUE):\n")
sub <- results[results$Is_monophyletic == "TRUE", ]
for (j in seq_len(nrow(sub)))
  cat(sprintf("  %s  Chunk %d  (%d carriers)\n",
              sub$Accession[j], sub$Chunk[j], sub$N_carriers[j]))

cat("\nNOT monophyletic (FALSE):\n")
sub <- results[results$Is_monophyletic == "FALSE", ]
for (j in seq_len(nrow(sub)))
  cat(sprintf("  %s  Chunk %d  (%d carriers)\n",
              sub$Accession[j], sub$Chunk[j], sub$N_carriers[j]))
