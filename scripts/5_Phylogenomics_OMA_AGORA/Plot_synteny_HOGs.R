library(syntenyPlotteR)
library(patchwork)


draw.ideogram_V2 <- function(file_data, sizefile, output, directory = NULL, fileformat = "png", colours = colours.default, w = 8.5, h = 10, ps = 5, species="NA") {
  
  if (is.null(directory)) {
    directory <- tempdir()
  }
  
  colours.default <- c(
    "1" = "#BFD73B", "2" = "#39ACE2", "3" = "#F16E8A",
    "4" = "#2DB995", "5" = "#855823", "6" = "#A085BD",
    "7" = "#2EB560", "8" = "#D79128", "9" = "#FDBB63",
    "10" = "#AFDFE5", "11" = "#BF1E2D", "12" = "purple4",
    "13" = "#B59F31", "14" = "#F68B1F", "15" = "#EF374B",
    "16" = "#D376FF", "17" = "#009445", "18" = "#CE4699",
    "19" = "#7C9ACD", "20" = "#84C441", "21" = "#404F23",
    "22" = "#607F4B", "23" = "#EBB4A9", "24" = "#F6EB83",
    "25" = "#915F6D", "26" = "#602F92", "27" = "#81CEC6",
    "28" = "#F8DA04", "29" = "peachpuff2", "30" = "gray85", "33" = "peachpuff3",
    "W" = "#9590FF", "Z" = "#666666", "Y" = "#9590FF", "X" = "#666666",
    "LGE22" = "grey", "LGE64" = "gray64",
    "1A" = "pink", "1B" = "dark blue", "4A" = "light green",
    "Gap" = "white"
  )
  
    size <- tarstart <- tarend <- refchr <- ystart <- yend <- NULL
  data <- utils::read.delim(file_data, header = FALSE)
  
  #colnames(data) <- c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "orien", "tar", "ref")
  colnames(data) <- c("refchr", "refstart", "refend", "tarchr", "tarstart", "tarend", "orien", "ref", "tar")
  
  data$tarstart <- as.numeric(gsub(",", "", data$tarstart))
  data$tarend <- as.numeric(gsub(",", "", data$tarend))
  data$refstart <- as.numeric(gsub(",", "", data$refstart))
  data$refend <- as.numeric(gsub(",", "", data$refend))
  
  sizes <- utils::read.delim(sizefile, header = FALSE) # to be consistent with naming in EH
  names(sizes) <- c("chromosome", "size", "species")
  sizes$size <- as.numeric(gsub(",", "", sizes$size))
  
  ref <- unique(data$ref)
  tar <- unique(data$tar)
  
  tar_sizes <- sizes[sizes$species == tar, ]
  ref_sizes <- sizes[sizes$species == ref, ]
  
  colnames(tar_sizes) <- c("tarchr", "size")
  colnames(ref_sizes) <- c("refchr", "size")
  
  #data$tarchr <- factor(data$tarchr, levels = tar_sizes$tarchr)
  data$tarchr <- factor(
    data$tarchr,
    levels = unique(data$tarchr[order(as.numeric(gsub("\\D", "", data$tarchr)))])
  )
  data$refchr <- factor(
    data$refchr,
    levels = unique(data$refchr[order(as.numeric(gsub("\\D", "", data$refchr)))])
  )
  
  
  for (i in c(1:nrow(data))) {
    dir <- data[i, "orien"]
    chr <- data[i, "refchr"]
    full_len <- ref_sizes[ref_sizes$refchr == chr, 2]
    y1 <- round(data[i, "refstart"] / full_len, digits = 4)
    y2 <- round(data[i, "refend"] / full_len, digits = 4)
    
    inverted <- grepl("-", dir, fixed = TRUE)
    if (inverted == TRUE) {
      data[i, "ystart"] <- y2
      data[i, "yend"] <- y1
    } else {
      data[i, "ystart"] <- y1
      data[i, "yend"] <- y2
    }
  }

  plots <- ggplot2::ggplot( data = data) +
    ggplot2::geom_rect(
      data = tar_sizes, mapping = ggplot2::aes(xmin = 1, xmax = size, ymin = -0.1, ymax = 1.1),
      fill = "white", color = "black", alpha = 0.85, linewidth = 0.2
    ) +
    ggplot2::geom_rect(
      data = data, mapping = ggplot2::aes(xmin = tarstart, xmax = tarend, ymin = -0.1, ymax = 1.1, fill = refchr),
      color = NA, alpha = 0.85, linewidth = 0.2
    ) +
    ggplot2::geom_segment(data = data, mapping = ggplot2::aes(x = tarstart, y = ystart, xend = tarend, yend = yend), linewidth = 0.2) +
    ggplot2::facet_grid(as.factor(tarchr) ~ .) +
    ggplot2::labs(fill = "Reference", x = "Chomosome length (Genes)", title=species) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(size = 0.2),
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold", size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.title.align = 0.5
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::scale_fill_manual(values = colours, drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = c(
        0, 2.5e+07, 5e+07, 7.5e+07, 1e+08, 1.25e+08, 1.5e+08, 1.75e+08,
        2e+08, 2.25e+08, 2.5e+08, 2.75e+08, 3e+08, 3.25e+08, 3.5e+08
      ),
      labels = c(
        "0", "25", "50", "75", "100", "125", "150", "175",
        "200", "225", "250", "275", "300", "325", "350"
      )
    )
  #message(paste0("Saving ideogram image to ", directory))
  #print(plots)
  #ggplot2::ggsave(paste0(directory,"/",output, ".", fileformat), plots, device = fileformat, width = w, height = h, pointsize = ps)
  #return(plots)
}

colours_ancestralKaryo <- c(
  "CAR_1" = "#BFD73B", "CAR_2" = "#39ACE2", "CAR_3" = "#F16E8A",
  "CAR_4" = "#2DB995", "CAR_5" = "#855823", "CAR_6" = "#A085BD",
  "CAR_7" = "#2EB560", "CAR_8" = "#D79128", "CAR_9" = "#FDBB63",
  "CAR_10" = "#AFDFE5", "CAR_11" = "#BF1E2D", "CAR_12" = "purple4",
  "CAR_13" = "#B59F31", "CAR_14" = "#F68B1F", "CAR_15" = "#EF374B",
  "CAR_16" = "#D376FF", "CAR_17" = "#009445", "CAR_18" = "#CE4699",
  "CAR_19" = "#7C9ACD", "CAR_20" = "#84C441", "CAR_21" = "#404F23",
  "CAR_22" = "#607F4B", "CAR_23" = "#EBB4A9", "CAR_24" = "#F6EB83",
  "CAR_25" = "#915F6D", "CAR_26" = "#602F92", "CAR_27" = "#81CEC6",
  "CAR_28" = "#F8DA04")


# Try for AKA
file_data = "A4_vs_AKA.AlignmentFile.txt"
sizefile = "A4_vs_AKA.ChromLength.txt"
output="A4_vs_AKA.Synteny"
draw.ideogram_V2(file_data, sizefile, output, directory=".",fileformat="pdf", 
                 colours=colours_ancestralKaryo, species="AKA")

# All Species
species_list <- c("AKA", "AKY", "ALL", "BIA", "CLA", 
                  "CRP", "EPH", "FRE", "LAT", "LAZ", "MCC", 
                  "OCE", "OMA", "POL", "PRC", "PRD", "SAN", "SEB")
#Only "missing" "AMPPE" from plot

# For the following, need to uncomment last lines in function draw.ideogram_V2
for (sp in species_list) {
  
  file_data <- paste0("A4_vs_", sp, ".AlignmentFile.txt")
  sizefile  <- paste0("A4_vs_", sp, ".ChromLength.txt")
  output    <- paste0("A4_vs_", sp, ".Synteny")
  
  draw.ideogram_V2(
    file_data,
    sizefile,
    output,
    directory = ".",
    fileformat = "pdf",
    colours = colours_ancestralKaryo, 
    species = sp
  )
}

# Plot all together
plot_list <- list()

for (sp in species_list) {
  
  file_data <- paste0("A4_vs_", sp, ".AlignmentFile.txt")
  sizefile  <- paste0("A4_vs_", sp, ".ChromLength.txt")
  output    <- paste0("A4_vs_", sp, ".Synteny")
  
  plot_list[[sp]] <- draw.ideogram_V2(
    file_data,
    sizefile,
    output,
    directory = ".",
    fileformat = "pdf",
    colours = colours_ancestralKaryo, 
    species = sp
  )
}


combined_plot <- wrap_plots(plot_list, ncol = 6) +
  plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "right")

pdf("Synteny_Alltogether.pdf", width = 40, height = 60)
print(combined_plot)
dev.off()


