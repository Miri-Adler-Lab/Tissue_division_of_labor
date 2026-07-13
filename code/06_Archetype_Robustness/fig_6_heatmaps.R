############################################################
# Jaccard Heatmap: Global vs Tissue-Specific Archetypes
# ----------------------------------------------------------
# This script creates Jaccard similarity heatmaps comparing
# global and tissue-specific archetypes
############################################################

## -----------------------------
## 1. Load required packages
## -----------------------------
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(stringr)

## -----------------------------
## 2. Setup paths (relative)
## -----------------------------

# Robust function to get script directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && ctx$path != "") {
      return(dirname(normalizePath(ctx$path)))
    }
  }
  if (sys.nframe() >= 1) {
    ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile)))
    }
  }
  return(getwd())
}

# Setup paths
input_root <- Sys.getenv("PUBLICATION_INPUT_ROOT")
output_root <- Sys.getenv("PUBLICATION_OUTPUT_ROOT")
if (!nzchar(input_root) || !nzchar(output_root)) {
  stop("PUBLICATION_INPUT_ROOT and PUBLICATION_OUTPUT_ROOT must both be set")
}

# Define directories
input_dir <- file.path(normalizePath(input_root), "06_Archetype_Robustness")
output_dir <- file.path(normalizePath(output_root), "Figure_6_Archetype_Robustness")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Cell types to process
cell_types <- c("Endothelial", "Fibroblast", "Macrophage")

## -----------------------------
## 3. User-defined settings
## -----------------------------

# Tissue colors for row annotation
tissue_colors <- c(
  liver = "#FF7F00",
  pancreas = "red",
  lung = "yellow",
  spleen = "yellow2",
  bladder = "red3",
  fat = "tan1",
  heart = "purple",
  small = "orange3",
  small_intestine = "orange3",
  thymus = "seagreen",
  large = "turquoise",
  large_intestine = "turquoise",
  trachea = "purple3",
  tongue = "pink"
)

# Legend settings
legend_max <- 0.4

## -----------------------------
## 4. Helper functions
## -----------------------------

format_sample_names <- function(names_vec) {
  names_vec %>%
    str_replace("(.*?)(\\d+)$", function(m) {
      tissue <- str_match(m, "^(.*?)(\\d+)$")[, 2]
      num <- str_match(m, "^(.*?)(\\d+)$")[, 3]
      formatted_tissue <- str_replace_all(tissue, "_", " ")
      paste0(formatted_tissue, "_", num)
    })
}

compute_jaccard_matrix <- function(global_archetypes, tissue_archetypes) {
  jmat <- matrix(
    0,
    nrow = length(global_archetypes),
    ncol = length(tissue_archetypes),
    dimnames = list(
      names(global_archetypes),
      names(tissue_archetypes)
    )
  )

  for (g_name in names(global_archetypes)) {
    g_genes <- global_archetypes[[g_name]]
    for (t_name in names(tissue_archetypes)) {
      t_genes <- tissue_archetypes[[t_name]]
      intersection <- length(intersect(g_genes, t_genes))
      union_size <- length(union(g_genes, t_genes))
      jaccard <- ifelse(union_size > 0, intersection / union_size, 0)
      jmat[g_name, t_name] <- jaccard
    }
  }
  jmat
}

############################################################
## Section 1: Global vs Tissue-Specific Archetypes
############################################################

for (cell_type in cell_types) {
  # Map cell type to filename pattern and global prefix
  if (cell_type == "Endothelial") {
    cell_type_file <- "endothelials"
    global_prefix <- "AE"
  } else if (cell_type == "Fibroblast") {
    cell_type_file <- "fibroblasts"
    global_prefix <- "AF"
  } else if (cell_type == "Macrophage") {
    cell_type_file <- "macrophages"
    global_prefix <- "AM"
  }

  # The merged files are in the root of 06_Archetype_Robustness
  input_csv <- file.path(input_dir, paste0("merged_archetypes_", cell_type_file, ".csv.gz"))
  output_prefix <- paste0(cell_type, "_jaccard")

  message("Processing: ", cell_type)
  if (!file.exists(input_csv)) {
    stop("File not found: ", input_csv)
  }

  # Read CSV
  df <- read_csv(input_csv, col_types = cols(.default = col_character()))

  # Convert each column to a character vector of genes (drop NAs)
  archetypes <- df %>% map(~ .x[!is.na(.x)])

  # Split into global vs tissue-specific archetypes
  global_archetypes <- archetypes[startsWith(names(archetypes), global_prefix)]
  tissue_archetypes <- archetypes[!startsWith(names(archetypes), global_prefix)]

  # Compute Jaccard similarity matrix (rows: global, cols: tissue)
  jaccard_mat <- compute_jaccard_matrix(global_archetypes, tissue_archetypes)
  heatmap_data <- t(jaccard_mat)

  # Format row names for plotting
  rownames(heatmap_data) <- format_sample_names(rownames(heatmap_data))

  # Extract tissue names for row annotation
  tissue_names <- str_extract(rownames(heatmap_data), "^[a-zA-Z]+")

  ## Build heatmap

  # Row annotation: tissue of origin
  anno_tissue <- rowAnnotation(
    Tissue = tissue_names,
    col = list(Tissue = tissue_colors),
    show_annotation_name = FALSE,
    width = unit(5, "mm")
  )

  # Continuous color map
  max_val <- max(heatmap_data, na.rm = TRUE)
  scale_max <- min(max_val, legend_max)
  col_fun <- colorRamp2(
    c(0, scale_max),
    c("white", "black")
  )

  # Apply col_fun manually per cell
  cell_fun <- function(j, i, x, y, width, height, fill) {
    grid.rect(
      x, y, width, height,
      gp = gpar(
        col  = NA,
        fill = col_fun(heatmap_data[i, j])
      )
    )
  }

  # Define legend ticks
  legend_breaks <- seq(0, legend_max, by = 0.1)
  legend_labels <- sprintf("%.1f", legend_breaks)

  ht_main <- Heatmap(
    heatmap_data,
    name = "Jaccard Similarity",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_rot = 0,
    column_names_gp = gpar(fontsize = 7),
    column_names_centered = TRUE,
    show_row_names = FALSE,
    cell_fun = cell_fun,
    heatmap_legend_param = list(
      at = legend_breaks,
      labels = legend_labels,
      direction = "horizontal",
      title_position = "topcenter"
    ),
    width = unit(ncol(heatmap_data) * 6, "mm"),
    height = unit(nrow(heatmap_data) * 6, "mm")
  )

  # Row names as separate annotation
  rowname_anno <- rowAnnotation(
    RowNames = anno_text(
      rownames(heatmap_data),
      gp = gpar(fontsize = 10),
      just = "left"
    ),
    show_annotation_name = FALSE,
    width = max_text_width(rownames(heatmap_data)) + unit(2, "mm")
  )

  # Combined heatmap object
  ht_combined <- ht_main + anno_tissue + rowname_anno

  ## Save heatmap

  # PNG
  png(
    file.path(output_dir, paste0(output_prefix, "_continuous_heatmap.png")),
    width  = 2000,
    height = 4000,
    res    = 300
  )
  draw(ht_combined,
    heatmap_legend_side = "top",
    annotation_legend_side = "right"
  )
  dev.off()

  # PDF
  pdf(
    file.path(output_dir, paste0(output_prefix, "_continuous_heatmap.pdf")),
    width = 10,
    height = 20
  )
  draw(ht_combined,
    heatmap_legend_side = "top",
    annotation_legend_side = "right"
  )
  dev.off()

  ## Standalone tissue legend

  tissue_legend <- Legend(
    labels    = names(tissue_colors),
    legend_gp = gpar(fill = tissue_colors),
    title     = "Tissue",
    ncol      = 1
  )

  # Tissue legend PNG
  png(
    file.path(output_dir, paste0(output_prefix, "_tissue_legend.png")),
    width  = 800,
    height = 600,
    res    = 150
  )
  draw(tissue_legend,
    x = unit(1, "cm"),
    y = unit(1, "cm"),
    just = c("left", "bottom")
  )
  dev.off()

  # Tissue legend PDF
  pdf(
    file.path(output_dir, paste0(output_prefix, "_tissue_legend.pdf")),
    width = 5.3,
    height = 4
  )
  draw(tissue_legend,
    x = unit(1, "cm"),
    y = unit(1, "cm"),
    just = c("left", "bottom")
  )
  dev.off()
} # End of Section 1

############################################################
## Section 2: Tissue-Specific vs Tissue-Specific
############################################################

for (cell_type in cell_types) {
  if (cell_type == "Endothelial") {
    cell_type_file <- "endothelials"
    global_prefix <- "AE"
  } else if (cell_type == "Fibroblast") {
    cell_type_file <- "fibroblasts"
    global_prefix <- "AF"
  } else if (cell_type == "Macrophage") {
    cell_type_file <- "macrophages"
    global_prefix <- "AM"
  }

  input_csv_tissue_only <- file.path(input_dir, paste0("merged_archetypes_", cell_type_file, ".csv.gz"))

  if (!file.exists(input_csv_tissue_only)) stop("File not found: ", input_csv_tissue_only)

  # Read data
  df_tissue <- read.csv(input_csv_tissue_only, stringsAsFactors = FALSE)

  # Convert each column to a cleaned gene set
  gene_sets <- lapply(df_tissue, function(col) {
    col <- as.character(col)
    col <- col[!is.na(col) & col != ""]
    unique(col)
  })

  # Filter out global archetypes
  gene_sets <- gene_sets[!startsWith(names(gene_sets), global_prefix)]

  # Archetype and tissue names
  archetype_names <- names(gene_sets)
  tissue_names <- gsub("[0-9]+$", "", archetype_names)

  # Pairwise Jaccard index
  n_arch <- length(gene_sets)
  jaccard_index <- matrix(
    0,
    nrow = n_arch,
    ncol = n_arch,
    dimnames = list(archetype_names, archetype_names)
  )

  for (i in seq_along(gene_sets)) {
    for (j in seq_along(gene_sets)) {
      set1 <- gene_sets[[i]]
      set2 <- gene_sets[[j]]
      intersection <- length(intersect(set1, set2))
      union_size <- length(union(set1, set2))
      jaccard_index[i, j] <- ifelse(union_size == 0, 0, intersection / union_size)
    }
  }

  # Build color scale from OFF-DIAGONAL values only
  off_diag <- jaccard_index
  diag(off_diag) <- NA

  min_val <- suppressWarnings(min(off_diag, na.rm = TRUE))
  max_val <- suppressWarnings(max(off_diag, na.rm = TRUE))

  if (!is.finite(min_val) || !is.finite(max_val)) {
    min_val <- 0
    max_val <- 1
  } else if (min_val == max_val) {
    min_val <- max(0, min_val - 1e-6)
    max_val <- min(1, max_val + 1e-6)
  }

  col_fun_self <- colorRamp2(c(min_val, max_val), c("white", "black"))

  # Mask diagonal
  plot_mat <- jaccard_index
  diag(plot_mat) <- NA

  # Annotations
  top_annot <- HeatmapAnnotation(
    Tissue = tissue_names,
    col = list(Tissue = tissue_colors),
    show_annotation_name = FALSE,
    height = unit(3, "mm")
  )

  left_annot <- rowAnnotation(
    Tissue = tissue_names,
    col = list(Tissue = tissue_colors),
    show_annotation_name = FALSE,
    width = unit(3, "mm")
  )

  # Draw heatmap
  ht_self <- Heatmap(
    plot_mat,
    name = "Jaccard",
    col = col_fun_self,
    na_col = "grey95",
    top_annotation = top_annot,
    left_annotation = left_annot,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 90,
    border = TRUE,
    column_names_gp = gpar(fontsize = 9),
    row_names_gp = gpar(fontsize = 9),
    width = unit(14, "cm"),
    height = unit(14, "cm"),
    row_dend_width = unit(0.5, "cm"),
    column_dend_height = unit(0.5, "cm"),
    heatmap_legend_param = list(title = "Jaccard")
  )

  # Save
  pdf(file.path(output_dir, paste0(cell_type, "_tissue_jaccard_heatmap.pdf")), width = 16, height = 8)
  draw(ht_self)
  dev.off()

  png(file.path(output_dir, paste0(cell_type, "_tissue_jaccard_heatmap.png")), width = 2000, height = 1000, res = 300)
  draw(ht_self)
  dev.off()
} # End of Section 2
