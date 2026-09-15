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


# FUNCTION for plotting SV statistics

### For data organization and plots of SV per species

make_comp_data <- function(df) {
  
  df %>%
    transmute(
      Species,
      plot_type = "comp_withsv",
      
      not_aligned = NotAligned_length / Analzed_Genome_length,
      
      syn = (Synteny_length - long_INS_length - long_DEL_length) /
        Analzed_Genome_length,
      
      inv = Inversion_length / Analzed_Genome_length,
      dup = Duplication_length / Analzed_Genome_length,
      trans = Transloc_length / Analzed_Genome_length,
      del = long_DEL_length / Analzed_Genome_length,
      ins = long_INS_length / Analzed_Genome_length
    ) %>%
    
    pivot_longer(
      cols = -c(Species, plot_type),
      names_to = "category",
      values_to = "value"
    )
}

make_sv_data <- function(df) {
  
  df %>%
    
    mutate(
      SV_length =
        Inversion_length +
        Transloc_length +
        Duplication_length +
        long_DEL_length +
        long_INS_length
    ) %>%
    
    transmute(
      Species,
      Reference,   # 👈 KEEP THIS
      plot_type = "sv",
      
      prop_inv = Inversion_length / SV_length,
      prop_dup = Duplication_length / SV_length,
      prop_transloc = Transloc_length / SV_length,
      prop_del = long_DEL_length / SV_length,
      prop_ins = long_INS_length / SV_length
    ) %>%
    
    pivot_longer(
      cols = -c(Species, Reference, plot_type),
      names_to = "category",
      values_to = "value"
    ) %>%
    
    mutate(
      category = factor(
        category,
        levels = c(
          "prop_inv",
          "prop_dup",
          "prop_transloc",
          "prop_del",
          "prop_ins"
        )
      )
    )
}

make_positions <- function(df,
                           within_gap = 1.5,
                           between_gap = 1.5) {
  
  species_order <- rev(unique(df$Species))
  plot_order <- rev(c("comp_withsv", "sv"))
  
  pos_df <- do.call(
    rbind,
    lapply(seq_along(species_order), function(i) {
      
      base <- (i - 1) * (2 * within_gap + between_gap)
      
      data.frame(
        Species = species_order[i],
        plot_type = plot_order,
        x = c(base, base + within_gap)
      )
    })
  )
  
  df %>%
    left_join(pos_df, by = c("Species", "plot_type"))
}

make_full_data <- function(df) {
  
  df %>%
    mutate(
      syn = (Synteny_length - long_INS_length - long_DEL_length) / Analzed_Genome_length,
      inv = Inversion_length / Analzed_Genome_length,
      trans = Transloc_length / Analzed_Genome_length,
      dup = Duplication_length / Analzed_Genome_length,
      ins = long_INS_length / Analzed_Genome_length,
      del = long_DEL_length / Analzed_Genome_length,
      not_aligned = NotAligned_length / Analzed_Genome_length
    ) %>%
    
    select(
      Species,
      Reference,
      syn, inv, trans, dup, not_aligned, ins, del
    ) %>%
    
    pivot_longer(
      cols = -c(Species, Reference),
      names_to = "category",
      values_to = "value"
    ) %>%
    
    mutate(
      category = factor(
        category,
        levels = c(
          "inv", "dup", "trans",
          "del", "ins", "not_aligned","syn"
        )
      )
    )
}

plot_sv_composition_perSpecies <- function(plot_data, label_df, title_text) {
  
  ggplot(plot_data, aes(x = x, y = value, fill = category)) +
    
    geom_bar(stat = "identity", width = 0.9) +
    
    scale_x_continuous(
      breaks = label_df$x,
      labels = label_df$label
    ) +
    
    coord_flip() +
    
    labs(
      title = title_text,
      x = "",
      y = "Proportion",
      fill = ""
    ) +
    
    # scale_fill_manual(values = c(
    #   "not_aligned" = "#EBEBEB",
    #   "syn" = "#D87B7B",
    #   "inv" = "#F09837",
    #   "dup" = "#56BBF9",
    #   "trans" = "#B7D86E",
    #   "del" = "#CFCACA",
    #   "ins" = "#949191",
    #   
    #   "prop_inv" = "#F09837",
    #   "prop_dup" = "#56BBF9",
    #   "prop_transloc" = "#B7D86E",
    #   "prop_del" = "#CFCACA",
    #   "prop_ins" = "#949191"
    # )) +
    # 
    
    scale_fill_manual(values = c(
      "not_aligned" = "#EBEBEB",
      "syn" = "#A9A9A9",
      "inv" = "#F09837",
      "dup" = "#56BBF9",
      "trans" = "#B7D86E",
      "del" = "#7B5EA7",
      "ins" = "#D87B7B",
      
      "prop_inv" = "#F09837",
      "prop_dup" = "#56BBF9",
      "prop_transloc" = "#B7D86E",
      "prop_del" = "#7B5EA7",
      "prop_ins" = "#D87B7B"
    )) +
    
    theme_minimal() +
    
    theme(
      axis.text.y = element_text(size = 9)
    )
}

plot_genome_breakdown <- function(df, reference, title_text = "Full Genome Breakdown") {
  
  ggplot(df,
         aes(x = Species, y = value, fill = category)) +
    
    geom_bar(stat = "identity") +
    
    coord_flip() +
    
    labs(
      title = title_text,
      x = "",
      y = "Proportion",
      fill = ""
    ) +
    
    # Old colors
    # scale_fill_manual(values = c(
    #   "not_aligned" = "#EBEBEB",
    #   "syn" = "#D87B7B",
    #   "inv" = "#F09837",
    #   "dup" = "#56BBF9",
    #   "trans" = "#B7D86E",
    #   "ins" = "#949191",
    #   "del" = "#CFCACA"
    # )) +
    
    scale_fill_manual(values = c(
      "not_aligned" = "#EBEBEB",
      "syn" = "#A9A9A9",
      "inv" = "#F09837",
      "dup" = "#56BBF9",
      "trans" = "#B7D86E",
      "ins" = "#D87B7B",
      "del" = "#7B5EA7"
    )) +
    
    theme_minimal()
}

sv_legend_scale <- function() {
  
  scale_fill_manual(
    
    breaks = c("syn", "not_aligned", "del", "ins", "trans", "dup", "inv"),
    
    values = c("not_aligned" = "#EBEBEB", "syn" = "#A9A9A9", "inv" = "#F09837",
      "dup" = "#56BBF9", "trans" = "#B7D86E", "ins" = "#D87B7B",
      "del" = "#7B5EA7"),
    
    labels = c(
      "syn" = "Synteny",
      "not_aligned" = "Unaligned",
      "del" = "DEL",
      "ins" = "INS",
      "trans" = "TRANS",
      "dup" = "DUP",
      "inv" = "INV"
    )
  )
}

plot_sv_composition <- function(df,reference, title_text = "SV Composition") {
  
  ggplot(df,
         aes(x = Species, y = value, fill = category)) +
    
    geom_bar(stat = "identity") +
    coord_flip() +
    
    labs(
      title = title_text,
      x = "",
      y = "Proportion",
      fill = "SV Type"
    ) +
    
    # scale_fill_manual(values = c(
    #   "prop_inv" = "#F09837",
    #   "prop_transloc" = "#B7D86E",
    #   "prop_dup" = "#56BBF9",
    #   "prop_ins" = "#949191",
    #   "prop_del" = "#CFCACA"
    
    scale_fill_manual(values = c(
      "prop_inv" = "#F09837",
      "prop_transloc" = "#B7D86E",
      "prop_dup" = "#56BBF9",
      "prop_ins" = "#D87B7B",
      "prop_del" = "#7B5EA7"
    )) +
    
    theme_minimal()
}

### For data organization and plots for merged SV

read_merged_sv <- function(prefix, label, sv_types, suffix, add_meta_fn) {
  
  sv_list <- lapply(sv_types, function(sv) {
    df <- read.table(
      paste0(prefix, sv, ".1kb_trsh.AllSVs", suffix),
      header = TRUE
    )
    
    add_meta_fn(df, label, "Minimap", sv)
  })
  
  dplyr::bind_rows(sv_list)
}

add_meta <- function(df, ref, aligner, svtype) {
  df %>%
    mutate(
      Reference = ref,
      Aligner = aligner,
      SV = svtype
    )
}

check_right_species <- function(data, ref) {
  data %>%
    filter(Reference == ref) %>%
    separate_rows(Species, sep = ",") %>%
    pull(Species) %>%
    unique() %>%
    sort()
}

clean_bia_species <- function(species_str) {
  species_str %>%
    str_split(",") %>%
    map(~ setdiff(.x, "BIA")) %>%
    map_chr(~ paste(.x, collapse = ",")) %>%
    unlist()
}

plot_Merged_SV_proportion <- function(data, reference_name, plottitle="Proportion of SV Types") {
  
  ggplot(
    dplyr::filter(data, Reference == reference_name),
    aes(x = Metric, y = Proportion, fill = SVType)
  ) +
    geom_bar(stat = "identity", width = 0.6) +
    
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    
    # scale_fill_manual(values = c(
    #   "INV" = "#F09837",
    #   "DUP" = "#56BBF9",
    #   "TRANS" = "#B7D86E",
    #   "INS" = "grey",
    #   "DEL" = "grey60"
    # )) +
    
    scale_fill_manual(values = c(
      "INV" = "#F09837",
      "DUP" = "#56BBF9",
      "TRANS" = "#B7D86E",
      "INS" = "#D87B7B",
      "DEL" = "#7B5EA7"
    )) +
    
    labs(
      title = plottitle,
      x = "",
      y = "Proportion",
      fill = "SV Type"
    ) +
    
    theme_minimal() +
    coord_flip() +
    theme(axis.text.x = element_text(size = 12))
}

plot_SV_species_distribution <- function(data, reference_name, title=paste0("Number of species per SV - ", reference_name)) {
  
  plot_data <- data %>%
    dplyr::filter(Reference == reference_name) %>%
    dplyr::mutate(
      SpeciesGroup = ifelse(n_species == 1, "1", ">1")
    )
  
  ggplot(plot_data, aes(x = factor(n_species), fill = SpeciesGroup)) +
    geom_bar(
      width = 0.8,
      color = "#777777",   # darker outline
      linewidth = 0.3
    )  +
    
    scale_fill_manual(
      # values = c("1" = "grey30", ">1" = "grey60"),
      values = c("1" = "#DCDCDC", ">1" = "#DCDCDC"),
      guide = "none"
    ) +
    
    theme_bw() +
    
    labs(
      title=title,
      x = "Number of species",
      y = "Count"
    ) +
    theme(
      axis.line = element_blank(),      # removes axis lines
      panel.border = element_blank()   # removes the box
    )
}

# PCoA functions

# Functions to generate color shades
generate_clade_colors <- function(species, clade, base_color) {
  n <- length(species)
  shades <- lighten(base_color, seq(0, 0.4, length.out = n))
  names(shades) <- species
  return(shades)
}

generate_gray_clade_colors <- function(species, clade_index, n_clades) {
  
  n <- length(species)
  
  # split grayscale into non-overlapping ranges
  band_width <- 0.8 / n_clades
  
  start <- 0.1 + (clade_index - 1) * band_width
  end   <- start + band_width * 0.8
  
  shades <- gray(seq(start, end, length.out = n))
  
  names(shades) <- species
  
  return(shades)
}

# Function to perform and plot PCoA
run_pcoa_sv <- function(data, reference, k = 4) {
  
  # -------------------------
  # 1. Subset by reference
  # -------------------------
  df_ref <- data %>%
    dplyr::filter(Reference == reference)
  
  if (nrow(df_ref) < 2) return(NULL)
  
  # -------------------------
  # 2. Expand species
  # -------------------------
  df_long <- df_ref %>%
    tidyr::separate_rows(Species, sep = ",")
  
  # -------------------------
  # 3. Presence/absence matrix
  # -------------------------
  SV_pa <- df_long %>%
    dplyr::mutate(presence = 1) %>%
    dplyr::select(ID, Species, presence) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = Species,
      values_from = presence,
      values_fill = 0
    )
  
  mat <- as.data.frame(SV_pa)
  rownames(mat) <- mat$ID
  mat$ID <- NULL
  
  mat_t <- t(mat)
  
  if (nrow(mat_t) < 2) return(NULL)
  
  # -------------------------
  # 4. Distance + PCoA
  # -------------------------
  dist_mat <- vegan::vegdist(mat_t, method = "jaccard")
  pcoa_res <- cmdscale(dist_mat, k = k, eig = TRUE)
  
  # -------------------------
  # 5. Format output
  # -------------------------
  pcoa_df <- as.data.frame(pcoa_res$points)
  colnames(pcoa_df) <- paste0("PCoA", seq_len(k))
  
  eig_vals <- pcoa_res$eig
  var_explained <- eig_vals / sum(eig_vals)
  
  pcoa_df$Var1 <- var_explained[1]
  pcoa_df$Var2 <- var_explained[2]
  
  pcoa_df$Species <- rownames(pcoa_df)
  pcoa_df$Reference <- reference
  
  return(pcoa_df)
}

plot_pcoa_sv <- function(pcoa_data, reference, species_colors = NULL,
                        title=paste0("PCoA - ", reference)) {
  
  df <- pcoa_data %>%
    dplyr::filter(Reference == reference)
  
  var1 <- round(df$Var1[1] * 100, 1)
  var2 <- round(df$Var2[1] * 100, 1)
  
  ggplot(df, aes(x = PCoA1, y = PCoA2, color = Species)) +
    
    geom_point(size = 3) +
    
    ggrepel::geom_text_repel(
      aes(label = Species),
      size = 3,
      max.overlaps = Inf
    ) +
    
    {
      if (!is.null(species_colors)) {
        scale_color_manual(values = species_colors)
      } else {
        scale_color_discrete()
      }
    } +
    
    labs(
      title = title,
      x = paste0("PCoA1 (", var1, "%)"),
      y = paste0("PCoA2 (", var2, "%)")
    ) +
    
    theme_bw() +
    theme(legend.position = "none")
}


# Functions for supplementary table

sv_by_type_theme <- theme_bw(base_size = 13) +
  theme(
    strip.background   = element_rect(fill = "#2C3E50", colour = NA),
    strip.text         = element_text(colour = "white", face = "bold", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(face = "bold"),
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 15),
    plot.subtitle      = element_text(colour = "grey40", size = 11)
  )


plot_sv_by_type <- function(data, Reference_species) {
  
  df <- data %>% filter(Reference == Reference_species)
  
  if (nrow(df) == 0) {
    stop(paste("No data found for Reference:", Reference_species))
  }
  
  if (Reference_Species == "BIA") {CompleteSpeciesName="A. biaculeatus"}else{
    CompleteSpeciesName="A. clarkii"
  }
  
  sv_summary <- df %>%
    group_by(SVType) %>%
    summarise(
      Count       = n(),
      TotalLength = sum(SVLength, na.rm = TRUE),
      .groups = "drop"
    )
  
  title_count <- bquote(paste("SV Count by Type — relative to ",
                               italic(.(CompleteSpeciesName)), sep=""))
  
  p_count <- ggplot(sv_summary, aes(x = SVType, y = Count, fill = SVType)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = Count), vjust = -0.4, fontface = "bold", size = 3.8) +
    scale_fill_manual(values = sv_colours) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = title_count,
      subtitle = "Total number of structural variants per SV type",
      x        = "SV Type",
      y        = "Count"
    ) +
    sv_by_type_theme
  
  title_length <- bquote(paste("Total SV Length by Type — relative to ",
          italic(.(CompleteSpeciesName)), sep=""))
  
  p_length <- ggplot(sv_summary, aes(x = SVType, y = TotalLength / 1e6, fill = SVType)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
    geom_text(
      aes(label = sprintf("%.1f", TotalLength / 1e6)),
      vjust = -0.4, fontface = "bold", size = 3.8
    ) +
    scale_fill_manual(values = sv_colours) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = title_length,
      subtitle = "Cumulative length of structural variants (Mb)",
      x        = "SV Type",
      y        = "Total Length (Mb)"
    ) +
    sv_by_type_theme
  
  list(count = p_count, length = p_length)
}

plot_sv_length_distribution <- function(data, Reference_species) {
  
  if (Reference_Species == "BIA") {CompleteSpeciesName="A. biaculeatus"}else{
    CompleteSpeciesName="A. clarkii"
  }

  df <- data %>% filter(Reference == Reference_species)
  
  if (nrow(df) == 0) {
    stop(paste("No data found for Reference:", Reference_species))
  }
  
  
  title_violin <- bquote(paste("SV Length Distribution by Type — relative to ",
                              italic(.(CompleteSpeciesName)), sep=""))
  
  p_violin <- ggplot(df, aes(x = SVType, y = SVLength, fill = SVType)) +
    geom_violin(trim = FALSE, alpha = 0.8, colour = "white", linewidth = 0.3) +
    geom_boxplot(width = 0.08, fill = "white", colour = "grey30",
                 outlier.shape = NA, linewidth = 0.4) +
    scale_fill_manual(values = sv_colours) +
    scale_y_log10(
      labels = scales::label_log(digits = 2),
      breaks = scales::breaks_log(n = 6)
    ) +
    labs(
      title    = title_violin,
      subtitle = "Log10-scaled SV length per type; inner boxplot shows median and IQR",
      x        = "SV Type",
      y        = "SV Length (bp, log10)"
    ) +
    sv_by_type_theme
  
  list(violin = p_violin)
}

plot_sv_chrom <- function(data, Reference_species) {
  
  df <- data %>% filter(Reference == Reference_species)
  
  if (Reference_Species == "BIA") {CompleteSpeciesName="A. biaculeatus"}else{
    CompleteSpeciesName="A. clarkii"
  }
  
  if (nrow(df) == 0) {
    stop(paste("No data found for Reference:", Reference_species))
  }
  
  # Order chromosomes naturally (chr1, chr2, ..., chrX, chrY)
  chrom_order <- unique(df$Chrom[order(as.integer(gsub("chr([0-9]+).*", "\\1", df$Chrom)),
                                       gsub("chr", "", df$Chrom))])
  
  sv_summary <- df %>%
    group_by(Chrom, SVType) %>%
    summarise(
      Count       = n(),
      TotalLength = sum(SVLength, na.rm = TRUE) / 1e6,
      .groups     = "drop"
    ) %>%
    mutate(Chrom = factor(Chrom, levels = chrom_order))
  
  chrom_theme <- theme_bw(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold", size = 10),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(colour = "grey40", size = 10)
    )
  
  p_count_tile <- bquote(paste("SV Count by Chromosome — relative to ",
                               italic(.(CompleteSpeciesName)), sep=""))
  
  p_count <- ggplot(sv_summary, aes(x = Chrom, y = Count, fill = SVType)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = sv_colours) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = p_count_tile,
      subtitle = "Stacked by SV type",
      x        = "Chromosome",
      y        = "Count",
      fill     = "SV type"
    ) +
    chrom_theme
  
  p_length_tile <- bquote(paste("SV Length by Chromosome — relative to ",
                               italic(.(CompleteSpeciesName)), sep=""))
  
  p_length <- ggplot(sv_summary, aes(x = Chrom, y = TotalLength, fill = SVType)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = sv_colours) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = p_length_tile,
      subtitle = "Stacked by SV type (Mb)",
      x        = "Chromosome",
      y        = "Total length (Mb)",
      fill     = "SV type"
    ) +
    chrom_theme
  
  list(count = p_count, length = p_length)
}

# PCoA by type

run_pcoa_sv_bytype <- function(data, reference, k = 4) {
  
  sv_types <- unique(data$SVType)
  
  results <- lapply(sv_types, function(svt) {
    
    df_ref <- data %>%
      dplyr::filter(Reference == reference, SVType == svt)
    
    if (nrow(df_ref) < 2) return(NULL)
    
    df_long <- df_ref %>%
      tidyr::separate_rows(Species, sep = ",")
    
    SV_pa <- df_long %>%
      dplyr::mutate(presence = 1) %>%
      dplyr::select(ID, Species, presence) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(
        names_from  = Species,
        values_from = presence,
        values_fill = 0
      )
    
    mat <- as.data.frame(SV_pa)
    rownames(mat) <- mat$ID
    mat$ID <- NULL
    
    mat_t <- t(mat)
    
    if (nrow(mat_t) < 2) return(NULL)
    
    dist_mat <- vegan::vegdist(mat_t, method = "jaccard")
    pcoa_res <- cmdscale(dist_mat, k = k, eig = TRUE)
    
    pcoa_df <- as.data.frame(pcoa_res$points)
    colnames(pcoa_df) <- paste0("PCoA", seq_len(k))
    
    eig_vals      <- pcoa_res$eig
    var_explained <- eig_vals / sum(eig_vals)
    
    pcoa_df$Var1      <- var_explained[1]
    pcoa_df$Var2      <- var_explained[2]
    pcoa_df$Species   <- rownames(pcoa_df)
    pcoa_df$Reference <- reference
    pcoa_df$SVType    <- svt
    
    pcoa_df
  })
  
  names(results) <- sv_types
  Filter(Negate(is.null), results)  # drop SVTypes with too few data
}

plot_pcoa_sv_bytype <- function(pcoa_bytype, reference, species_colors = NULL) {
  
  lapply(names(pcoa_bytype), function(svt) {
    
    df <- pcoa_bytype[[svt]] %>%
      dplyr::filter(Reference == reference)
    
    if (nrow(df) == 0) return(NULL)
    
    var1 <- round(df$Var1[1] * 100, 1)
    var2 <- round(df$Var2[1] * 100, 1)
    
    ggplot(df, aes(x = PCoA1, y = PCoA2, color = Species)) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(
        aes(label = Species),
        size        = 3,
        max.overlaps = Inf
      ) +
      {
        if (!is.null(species_colors)) {
          scale_color_manual(values = species_colors)
        } else {
          scale_color_discrete()
        }
      } +
      labs(
        title = paste0("PCoA - ", reference, " | ", svt),
        x     = paste0("PCoA1 (", var1, "%)"),
        y     = paste0("PCoA2 (", var2, "%)")
      ) +
      theme_bw() +
      theme(legend.position = "none")
    
  }) |> setNames(names(pcoa_bytype))
}



