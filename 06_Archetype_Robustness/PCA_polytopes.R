############################################################
# PCA + Archetype Polytope Visualization Pipeline
# ----------------------------------------------------------
# This script:
#   1. Reads cell-by-gene expression matrices (CSV) per
#      tissue and cell type.
#   2. Reads archetype coordinates (arcOrig) from .mat files.
#   3. Runs a joint PCA on cells + archetypes.
#   4. Plots cells in 3D PCA space, with archetypes
#      forming a polytope and numbered vertices.
#   5. Exports HTML, PNG, and PDF per
#      (tissue, cell type).
#
# Expected input directory structure:
#   base_dir/
#     <tissue>/
#       <cell_type>/
#         <tissue>_scrna.csv       # cells × genes
#         <tissue>_arcOrig.mat     # MATLAB file with 'arcOrig'
#
# Outputs:
#   output_dir/
#     <tissue>_<cell_type>.html  # 3D plot
#     <tissue>_<cell_type>.png   # PNG snapshot
#     <tissue>_<cell_type>.pdf   # PDF snapshot
############################################################


# install.packages("R.matlab")

####################################################
### ---PCA plots + archetypes(export all plots)---###
####################################################


library(R.matlab)
library(data.table)
library(plotly)
library(webshot2)
library(magick)
library(htmlwidgets)

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

# Create export folder
output_dir <- file.path(normalizePath(output_root), "Figure_6_Archetype_Robustness")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define paths
base_dir <- file.path(normalizePath(input_root), "06_Archetype_Robustness")

# Cell types and expected tissues
cell_types <- c("Macrophage", "Fibroblast", "Endothelial")
tissues <- c("Fat", "Liver", "Lung", "Bladder", "Heart", "Large_Intestine", "Pancreas", "Small_Intestine", "Spleen", "Thymus", "Tongue", "Trachea")

panel_manifest <- expand.grid(
  tissue = tissues,
  cell_type = cell_types,
  stringsAsFactors = FALSE
)
excluded_panels <- c(
  "Liver/Fibroblast",
  "Heart/Endothelial",
  "Pancreas/Fibroblast",
  "Spleen/Fibroblast"
)
panel_manifest$key <- paste(panel_manifest$tissue, panel_manifest$cell_type, sep = "/")
panel_manifest <- panel_manifest[!panel_manifest$key %in% excluded_panels, ]

panel_file <- function(cell_type, tissue, kind) {
  panel_dir <- file.path(base_dir, tissue, cell_type)
  patterns <- c(
    expression = "_scrna(_endo)?\\.csv\\.gz$",
    arc = "_arc\\.mat$",
    arcOrig = "_arcOrig\\.mat$"
  )
  matches <- list.files(panel_dir, pattern = patterns[[kind]], full.names = TRUE)
  if (length(matches) != 1) {
    stop(sprintf(
      "Expected exactly one %s file for %s/%s; found %d",
      kind, tissue, cell_type, length(matches)
    ))
  }
  matches[[1]]
}

for (i in seq_len(nrow(panel_manifest))) {
  for (kind in c("expression", "arc", "arcOrig")) {
    panel_file(panel_manifest$cell_type[[i]], panel_manifest$tissue[[i]], kind)
  }
}

# Load expression data and run PCA
load_expression_and_pca <- function(cell_type, tissue) {
  file_path <- panel_file(cell_type, tissue, "expression")
  expr <- fread(file_path, header = FALSE)
  pca <- prcomp(expr, center = TRUE, scale. = TRUE)
  return(pca$x[, 1:3]) # only first 3 PCs
}

# Load archetype data from .mat file
load_archetypes <- function(cell_type, tissue) {
  mat_path <- panel_file(cell_type, tissue, "arc")
  mat_data <- readMat(mat_path)
  archetypes <- mat_data$arc # <-- adjusted to your file's key
  return(archetypes)
}

# Updated plot function with clean background
plot_pca_polytope <- function(pca_coords, archetypes) {
  if (is.null(pca_coords)) {
    return(plotly_empty())
  }

  p <- plot_ly() %>%
    add_markers(
      x = pca_coords[, 1], y = pca_coords[, 2], z = pca_coords[, 3],
      type = "scatter3d", mode = "markers",
      marker = list(size = 1, color = "black"),
      name = "Cells", showlegend = FALSE
    )

  if (!is.null(archetypes)) {
    archetypes <- as.data.frame(archetypes)
    colnames(archetypes) <- c("x", "y", "z")

    p <- p %>%
      add_markers(
        data = archetypes, x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        marker = list(size = 4, color = "red"),
        name = "Archetypes", showlegend = FALSE
      )

    for (i in seq_len(nrow(archetypes) - 1)) {
      for (j in seq.int(i + 1, nrow(archetypes))) {
        p <- p %>%
          add_trace(
            x = c(archetypes$x[i], archetypes$x[j]),
            y = c(archetypes$y[i], archetypes$y[j]),
            z = c(archetypes$z[i], archetypes$z[j]),
            type = "scatter3d", mode = "lines",
            line = list(color = "blue", width = 2),
            showlegend = FALSE
          )
      }
    }
  }

  # Minimal clean layout
  p <- layout(p,
    scene = list(
      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
      zaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
      bgcolor = "white"
    ),
    margin = list(l = 0, r = 0, b = 0, t = 0)
  )

  return(p)
}

# Main export loop

html_paths_basic <- character()
png_paths_basic <- character()
pdf_paths_basic <- character()
for (panel_i in seq_len(nrow(panel_manifest))) {
    tissue <- panel_manifest$tissue[[panel_i]]
    cell_type <- panel_manifest$cell_type[[panel_i]]
    message("Exporting: ", tissue, " / ", cell_type)

    pca_coords <- load_expression_and_pca(cell_type, tissue)
    archetypes <- load_archetypes(cell_type, tissue)

    # Build plot
    p <- plot_pca_polytope(pca_coords, archetypes)

    # Create tissue-specific output directory
    tissue_dir <- file.path(output_dir, tissue)
    dir.create(tissue_dir, showWarnings = FALSE)

    # Paths
    base_name <- paste0(tissue, "_", cell_type)
    html_path <- file.path(tissue_dir, paste0(base_name, ".html"))
    png_path <- file.path(tissue_dir, paste0(base_name, ".png"))
    pdf_path <- file.path(tissue_dir, paste0(base_name, ".pdf"))

    # Export
    saveWidget(p, html_path, selfcontained = FALSE)
    html_paths_basic <- c(html_paths_basic, html_path)
    png_paths_basic <- c(png_paths_basic, png_path)
    pdf_paths_basic <- c(pdf_paths_basic, pdf_path)
}
webshot2::webshot(html_paths_basic, png_paths_basic, vwidth = 800, vheight = 800)
for (i in seq_along(png_paths_basic)) {
  image_write(image_read(png_paths_basic[[i]]), path = pdf_paths_basic[[i]], format = "pdf")
}


####################################################
##### ---PCA plots + archetypes(using arcOrig)---####
####################################################

# Re-load libraries just in case (as per original script structure)
library(R.matlab)
library(data.table)
library(plotly)
library(webshot2)
library(magick)
library(htmlwidgets)

# Create export folder (already done, but keeping structure)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

load_archetypes_orig <- function(cell_type, tissue) {
  mat_path <- panel_file(cell_type, tissue, "arcOrig")
  mat_data <- readMat(mat_path)
  if (!"arcOrig" %in% names(mat_data)) stop("arcOrig not found in .mat file")
  archetypes <- mat_data$arcOrig # archetypes in gene expression space
  return(archetypes)
}

run_joint_pca <- function(expr_cells, expr_archetypes) {
  combined <- rbind(expr_cells, expr_archetypes)
  pca <- prcomp(combined, center = TRUE, scale. = TRUE)
  pcs <- pca$x[, 1:3]
  n_cells <- nrow(expr_cells)
  list(
    cell_pcs = pcs[1:n_cells, ],
    arc_pcs = pcs[(n_cells + 1):nrow(pcs), ]
  )
}


plot_pca_polytope_joint <- function(cell_pcs, arc_pcs) {
  if (is.null(cell_pcs)) {
    return(plotly_empty())
  }

  arc_df <- as.data.frame(arc_pcs)
  colnames(arc_df) <- c("x", "y", "z")

  p <- plot_ly() %>%
    add_markers(
      x = cell_pcs[, 1], y = cell_pcs[, 2], z = cell_pcs[, 3],
      type = "scatter3d", mode = "markers",
      marker = list(size = 3, color = "black"),
      name = "Cells", showlegend = FALSE
    ) %>%
    add_markers(
      data = arc_df, x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      marker = list(size = 4, color = "red"),
      name = "Archetypes", showlegend = FALSE
    )

  # Add lines between archetypes
  for (i in seq_len(nrow(arc_df) - 1)) {
    for (j in seq.int(i + 1, nrow(arc_df))) {
      p <- p %>%
        add_trace(
          x = c(arc_df$x[i], arc_df$x[j]),
          y = c(arc_df$y[i], arc_df$y[j]),
          z = c(arc_df$z[i], arc_df$z[j]),
          type = "scatter3d", mode = "lines",
          line = list(color = "blue", width = 2),
          showlegend = FALSE
        )
    }
  }


  p <- layout(p,
    scene = list(
      xaxis = list(
        showgrid = FALSE, zeroline = FALSE,
        showticklabels = FALSE, title = "", visible = FALSE
      ),
      yaxis = list(
        showgrid = FALSE, zeroline = FALSE,
        showticklabels = FALSE, title = "", visible = FALSE
      ),
      zaxis = list(
        showgrid = FALSE, zeroline = FALSE,
        showticklabels = FALSE, title = "", visible = FALSE
      ),
      bgcolor = "white"
    ),
    margin = list(l = 0, r = 0, b = 0, t = 0),
    paper_bgcolor = "white",
    plot_bgcolor = "white"
  )

  return(p)
}


html_paths_arc_orig <- character()
png_paths_arc_orig <- character()
pdf_paths_arc_orig <- character()
for (panel_i in seq_len(nrow(panel_manifest))) {
    tissue <- panel_manifest$tissue[[panel_i]]
    cell_type <- panel_manifest$cell_type[[panel_i]]
    message("Exporting (arcOrig): ", tissue, " / ", cell_type)

    csv_path <- panel_file(cell_type, tissue, "expression")
    expr_cells <- fread(csv_path, header = FALSE)
    expr_arcs <- load_archetypes_orig(cell_type, tissue)

    if (nrow(expr_cells) == 0) stop(sprintf("No cells for %s/%s", tissue, cell_type))

    if (ncol(expr_cells) != ncol(expr_arcs)) {
      stop(sprintf("Column mismatch for %s/%s: cells=%d, arcs=%d", tissue, cell_type, ncol(expr_cells), ncol(expr_arcs)))
    }

    pca_out <- run_joint_pca(expr_cells, expr_arcs)

    # Plot
    p <- plot_pca_polytope_joint(pca_out$cell_pcs, pca_out$arc_pcs)

    # Create tissue-specific output directory
    tissue_dir <- file.path(output_dir, tissue)
    dir.create(tissue_dir, showWarnings = FALSE)

    base_name <- paste0(tissue, "_", cell_type, "_arcOrig")
    html_path <- file.path(tissue_dir, paste0(base_name, ".html"))
    png_path <- file.path(tissue_dir, paste0(base_name, ".png"))
    pdf_path <- file.path(tissue_dir, paste0(base_name, ".pdf"))

    saveWidget(p, html_path, selfcontained = FALSE)
    html_paths_arc_orig <- c(html_paths_arc_orig, html_path)
    png_paths_arc_orig <- c(png_paths_arc_orig, png_path)
    pdf_paths_arc_orig <- c(pdf_paths_arc_orig, pdf_path)
}
webshot2::webshot(html_paths_arc_orig, png_paths_arc_orig, vwidth = 800, vheight = 800)
for (i in seq_along(png_paths_arc_orig)) {
  image_write(image_read(png_paths_arc_orig[[i]]), path = pdf_paths_arc_orig[[i]], format = "pdf")
}


##################################################################
### ---PCA plots + archetypes(using arcOrig)+color per tissue---###
##################################################################



library(R.matlab)
library(data.table)
library(plotly)
library(webshot2)
library(magick)
library(htmlwidgets)

# Define fixed color palette per tissue
tissue_colors <- c(
  Liver            = "#FF7F00", # orange
  Pancreas         = "#FF0000", # red
  Lung             = "#FFFF00", # yellow
  Spleen           = "#EEEE00", # yellow2
  Bladder          = "#CD0000", # red3
  Fat              = "#FFA54F", # tan1
  Heart            = "#A020F0", # purple
  Small_Intestine  = "#CD8500", # orange3
  Thymus           = "#2E8B57", # seagreen
  Large_Intestine  = "#40E0D0", # turquoise
  Trachea          = "#7D26CD", # purple3
  Tongue           = "#FFC0CB" # pink
)

get_tissue_color <- function(t) {
  if (!t %in% names(tissue_colors)) {
    stop(sprintf("No color defined for tissue '%s'. Check tissue_colors keys.", t))
  }
  unname(tissue_colors[[t]])
}

# Plot function
plot_pca_polytope_color <- function(cell_pcs, arc_pcs, tissue) {
  if (is.null(cell_pcs)) {
    return(plotly_empty())
  }

  arc_df <- as.data.frame(arc_pcs)
  colnames(arc_df) <- c("x", "y", "z")
  arc_df$arc_id <- seq_len(nrow(arc_df)) # 1,2,3,...

  p <- plot_ly() %>%
    # Cells
    add_markers(
      x = cell_pcs[, 1], y = cell_pcs[, 2], z = cell_pcs[, 3],
      type = "scatter3d", mode = "markers",
      marker = list(
        size = 4, color = get_tissue_color(tissue), opacity = 0.85,
        line = list(color = "white", width = 0.3)
      ),
      name = paste("Cells -", tissue), showlegend = FALSE
    ) %>%
    # Archetype points
    add_markers(
      data = arc_df, x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      marker = list(size = 5, color = "red"),
      name = "Archetypes", showlegend = FALSE
    ) %>%
    # Archetype labels
    # Archetype labels with offset in Z
    add_trace(
      data = arc_df,
      x = ~x, y = ~y, z = ~ I(z + 0.05), # small bump in z
      type = "scatter3d", mode = "text",
      text = ~arc_id,
      textfont = list(size = 20, color = "black"),
      textposition = "top center",
      showlegend = FALSE
    )


  # Lines between archetypes (unchanged)
  for (i in seq_len(nrow(arc_df) - 1)) {
    for (j in seq.int(i + 1, nrow(arc_df))) {
      p <- p %>%
        add_trace(
          x = c(arc_df$x[i], arc_df$x[j]),
          y = c(arc_df$y[i], arc_df$y[j]),
          z = c(arc_df$z[i], arc_df$z[j]),
          type = "scatter3d", mode = "lines",
          line = list(color = "black", width = 2),
          showlegend = FALSE
        )
    }
  }

  layout(
    p,
    scene = list(
      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "", visible = FALSE),
      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "", visible = FALSE),
      zaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "", visible = FALSE),
      bgcolor = "white"
    ),
    margin = list(l = 0, r = 0, b = 0, t = 0),
    paper_bgcolor = "white",
    plot_bgcolor = "white"
  )
}

# Main loop
html_paths_color <- character()
png_paths_color <- character()
pdf_paths_color <- character()
for (panel_i in seq_len(nrow(panel_manifest))) {
    tissue <- panel_manifest$tissue[[panel_i]]
    cell_type <- panel_manifest$cell_type[[panel_i]]
    message("Exporting (color): ", tissue, " / ", cell_type)

    csv_path <- panel_file(cell_type, tissue, "expression")
    expr_cells <- fread(csv_path, header = FALSE)
    expr_arcs <- load_archetypes_orig(cell_type, tissue)
    if (nrow(expr_cells) == 0) stop(sprintf("No cells for %s/%s", tissue, cell_type))

    if (ncol(expr_cells) != ncol(expr_arcs)) {
      stop(sprintf("Column mismatch for %s/%s: genes=%d, arcs=%d", tissue, cell_type, ncol(expr_cells), ncol(expr_arcs)))
    }

    pca_out <- run_joint_pca(expr_cells, expr_arcs)
    p <- plot_pca_polytope_color(pca_out$cell_pcs, pca_out$arc_pcs, tissue)

    # Create tissue-specific output directory
    tissue_dir <- file.path(output_dir, tissue)
    dir.create(tissue_dir, showWarnings = FALSE)

    base_name <- paste0(tissue, "_", cell_type, "_color")
    html_path <- file.path(tissue_dir, paste0(base_name, ".html"))
    png_path <- file.path(tissue_dir, paste0(base_name, ".png"))
    pdf_path <- file.path(tissue_dir, paste0(base_name, ".pdf"))

    saveWidget(p, html_path, selfcontained = FALSE)
    html_paths_color <- c(html_paths_color, html_path)
    png_paths_color <- c(png_paths_color, png_path)
    pdf_paths_color <- c(pdf_paths_color, pdf_path)
}
webshot2::webshot(html_paths_color, png_paths_color, vwidth = 800, vheight = 800)
for (i in seq_along(png_paths_color)) {
  image_write(image_read(png_paths_color[[i]]), path = pdf_paths_color[[i]], format = "pdf")
}
