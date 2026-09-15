suppressMessages(library(topGO))

base    <- "GO_Enrichment"
go_dir  <- file.path(base, "GeneOntologies_perSpecies")
inv_dir <- file.path(base, "GenesInInversion_perSpecies")
out_dir <- file.path(base, "Results_perSpecies")
dir.create(out_dir, showWarnings = FALSE)

chr18_species <- c("AKA","AKY","ALL","CRP","LAT","MCC","OMA","POL","SAN","SEB")
chr9_species  <- c("POL","CLA")

runs <- c(
  lapply(chr18_species, function(sp) list(species = sp, inv = "chr18")),
  lapply(chr9_species,  function(sp) list(species = sp, inv = "chr9"))
)

node_sizes <- c(5, 10)

log_con <- file(file.path(out_dir, "run_log.txt"), open = "wt")

for (r in runs) {
  sp  <- r$species
  inv <- r$inv

  go_file  <- file.path(go_dir,  paste0(sp, "_GeneOntologies.txt"))
  inv_file <- file.path(inv_dir, paste0(sp, "_", inv, "_InversionGenes.txt"))

  geneID2GO <- readMappings(file = go_file)
  geneUniverse <- names(geneID2GO)

  genesOfInterest <- readLines(inv_file)
  genesOfInterest <- genesOfInterest[nzchar(genesOfInterest)]

  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse

  nInUniverse <- sum(geneUniverse %in% genesOfInterest)
  writeLines(sprintf("%s %s: universe=%d genesOfInterest_total=%d genesOfInterest_inUniverse=%d",
                      sp, inv, length(geneUniverse), length(genesOfInterest), nInUniverse), log_con)

  for (ns in node_sizes) {
    tag <- paste0(sp, "_", inv, "_nodeSize", ns)

    myGOdata <- tryCatch(
      suppressMessages(new("topGOdata", description = paste(sp, inv, "nodeSize", ns),
                            ontology = "BP", allGenes = geneList,
                            annot = annFUN.gene2GO, gene2GO = geneID2GO,
                            nodeSize = ns)),
      error = function(e) { writeLines(paste("ERROR building topGOdata for", tag, ":", conditionMessage(e)), log_con); NULL }
    )
    if (is.null(myGOdata)) next

    resultFisher <- tryCatch(
      runTest(myGOdata, algorithm = "weight01", statistic = "fisher"),
      error = function(e) { writeLines(paste("ERROR running test for", tag, ":", conditionMessage(e)), log_con); NULL }
    )
    if (is.null(resultFisher)) next

    nTested <- length(score(resultFisher))
    allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "classicFisher",
                        ranksOf = "classicFisher", topNodes = nTested)

    allRes$classicFisher_num <- suppressWarnings(as.numeric(sub("^< ", "", allRes$classicFisher)))
    allRes$Species   <- sp
    allRes$Inversion <- inv
    allRes$NodeSize  <- ns

    outfile <- file.path(out_dir, paste0(tag, ".txt"))
    write.table(allRes, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)

    writeLines(sprintf("  -> %s: %d GO terms tested, %d with p<0.05, written to %s",
                        tag, nTested, sum(allRes$classicFisher_num < 0.05, na.rm = TRUE), basename(outfile)),
               log_con)
  }
}

close(log_con)
cat("DONE\n")
