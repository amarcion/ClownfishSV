setwd("/Users/amarcion/Documents/PacBio_Sequencing_Assembly/10_Additional/15_Long_Reads_Mapping/mapping_Stats_BySpecies")

species=c("AKA", "AKY", "ALL", "BIA", "CLA", "CRP", "EPH", "FRE", "LAT",
          "LAZ", "MCC","OCE", "OMA", "POL", "PRC", "PRD", "SAN", "SEB")

COV_SUFFIX=".COV.stats"
FILT_COV_SUFFIX=".filt.COV.stats"
MAP_SUFFIX=".MAPQ.stats"
FILT_MAP_SUFFIX=".filt.MAPQ.stats"


#COV_SUFFIX
for (sp in species) {
  data <- read.table(paste0(sp, COV_SUFFIX))
  pdf(paste0(sp, COV_SUFFIX, ".pdf"), width = 8)
  par(mfrow=c(1,2))
  plot(data$V2, data$V3, type="l", main=sp, xlab="COV", ylab="Counts")
  plot(data$V2, data$V3, type="l", main=sp, xlab="COV", ylab="Counts", xlim=c(0,100))
  dev.off()
}

#COV_SUFFIX
for (sp in species) {
  data <- read.table(paste0(sp, COV_SUFFIX))
  pdf(paste0(sp, COV_SUFFIX, ".pdf"), width = 8)
  par(mfrow=c(1,2))
  plot(data$V2, data$V3, type="l", main=sp, xlab="COV", ylab="Counts")
  plot(data$V2, data$V3, type="l", main=sp, xlab="COV", ylab="Counts", xlim=c(0,100))
  dev.off()
}

#MAP_SUFFIX
for (sp in species) {
  data <- read.table(paste0(sp, MAP_SUFFIX))
  pdf(paste0(sp, MAP_SUFFIX, ".pdf"), width = 8)
  par(mfrow=c(1,2))
  plot(data$V2, data$V3, type="l", main=sp, xlab="MAPQ", ylab="Counts")
  plot(data$V2, data$V3, type="l", main=sp, xlab="MAPQ", ylab="Counts", xlim=c(0,100))
  dev.off()
}

#FILT_MAP_SUFFIX
for (sp in species) {
  data <- read.table(paste0(sp, FILT_MAP_SUFFIX))
  pdf(paste0(sp, FILT_MAP_SUFFIX, ".pdf"), width = 8)
  par(mfrow=c(1,2))
  plot(data$V2, data$V3, type="l", main=sp, xlab="MAPQ", ylab="Counts")
  plot(data$V2, data$V3, type="l", main=sp, xlab="MAPQ", ylab="Counts", xlim=c(0,100))
  dev.off()
}

