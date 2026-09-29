# ============================================================================
# Cell Type Identification Functions
# ============================================================================
# 
# Purpose: Identify cell types by scoring clusters based on marker expression
# 
# Functions:
#   - assign_celltype_from_dotplot(): borrow dotplot info for typing
#   - score_clusters_for_markers(): Score clusters for a set of markers
#   - identify_celltype(): Wrapper for single cell type
#   - identify_all_celltypes(): Batch process multiple cell types
# ============================================================================

# TODO: pulled from elsewhere, makes more sense here, cleanup
# Helper function to extract markers
sctb_extract_celltype_markers <- function(celltype_settings, celltypes_of_interest) {
  markers_list <- list()

  for (celltype in celltypes_of_interest) {
    # Find the row for this celltype
    row <- celltype_settings[celltype_settings$celltype == celltype, ]

    if (nrow(row) == 0) {
      warning(sprintf("Cell type '%s' not found in celltype_settings", celltype))
      next
    }

    # Parse comma-separated markers
    markers_str <- row$markers_used[1]
    markers <- trimws(strsplit(markers_str, ",")[[1]])

    markers_list[[celltype]] <- markers
  }

  return(markers_list)
}

# faster for purposes of sensitivity testing
#TODO: insert function signature description here
sctb_assign_celltype_from_dotplot <- function(seurat_obj, marker_lists, min_cells = 50,
                                        min_score = 3, min_mean_pct = 2,
                                         cluster_col = "seurat_clusters") {
  assignments <- data.frame()

  cluster_counts <- table(Idents(seurat_obj))
  clusters_to_type <- names(cluster_counts)[
    cluster_counts >= min_cells &
    names(cluster_counts) != "singleton"
  ]

  # avoid subset() — select cells directly
  cells_to_use <- colnames(seurat_obj)[Idents(seurat_obj) %in% clusters_to_type]
  typing_obj <- seurat_obj[, cells_to_use]

  for (celltype in names(marker_lists)) {
    markers <- marker_lists[[celltype]]
    dp <- DotPlot(typing_obj, features = markers)
    dp_data <- dp$data %>%
      group_by(id) %>%
      summarize(
        mean_exp = mean(avg.exp.scaled),
        mean_pct = mean(pct.exp),
        score = mean(avg.exp.scaled * pct.exp)
      ) %>%
      mutate(celltype = celltype)
    assignments <- rbind(assignments, dp_data)
  }

  best <- assignments %>%
    group_by(id) %>%
    slice_max(score, n = 1) %>%
    rename(cluster = id) %>%
    mutate(
        celltype = ifelse(score < min_score | mean_pct < min_mean_pct,
            "Unassigned", celltype),
        n_cells  = as.integer(cluster_counts[as.character(cluster)])
    )

  return(list(
    assignments = best,
    all_scores = assignments
  ))
}

#' Score clusters based on marker gene expression
#' 
#' @param all_markers Data frame from FindAllMarkers() with columns:
#'   gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj
#' @param marker_genes Character vector of marker gene symbols
#' @param celltype_name Optional name for the cell type (for output)
#' 
#' @return Data frame with cluster scores, ranked by best match
#' 
#' @export
sctb_score_clusters_for_markers <- function(all_markers,
                                       marker_genes,
                                       celltype_name = NULL) {

  # Find markers in cluster markers
  markers_in_clusters <- all_markers %>%
    dplyr::filter(gene %in% marker_genes)  # Explicit namespace

  if (nrow(markers_in_clusters) == 0) {
    warning(sprintf("No markers found for %s",
                   ifelse(is.null(celltype_name), "this cell type", celltype_name)))
    return(data.frame())
  }

  # Score each cluster
  cluster_scores <- markers_in_clusters %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(
      n_markers = n(),
      n_markers_total = length(marker_genes),
      mean_log2FC = mean(avg_log2FC),
      mean_pct1 = mean(pct.1),
      mean_pct2 = mean(pct.2),
      specificity = mean(pct.1) / (mean(pct.2) + 0.01),
      markers_found = paste(gene, collapse = ", "),
      markers_missing = paste(setdiff(marker_genes, gene), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      marker_coverage = n_markers / n_markers_total,
      score = marker_coverage * mean_log2FC * specificity
    ) %>%
    dplyr::arrange(desc(score))

  if (!is.null(celltype_name)) {
    cluster_scores$celltype <- celltype_name
  }

  return(cluster_scores)
}

#' Identify best cluster(s) for a single cell type
#' Intent is quick check for cell clusters matching marker genes list
#' 
#' @param all_markers Data frame from FindAllMarkers()
#' @param marker_genes Character vector of marker genes
#' @param celltype_name Name of cell type
#' @param min_markers Minimum number of markers required for assignment (default: 2)
#' @param min_score Minimum score required for assignment (default: 5)
#' @param verbose Print detailed output (default: TRUE)
#' 
#' @return List with:
#'   - best_cluster: Top ranked cluster
#'   - all_scores: Full scoring table
#'   - marker_details: Per-marker statistics for best cluster
#'   - assignment: Suggested cluster assignment
#' 
#' @export
identify_celltype <- function(all_markers,
                              marker_genes,
                              celltype_name,
                              min_markers = 2,
                              min_score = 5,
                              verbose = TRUE) {
  
  if (verbose) {
    cat(sprintf("\n=== Identifying %s ===\n", celltype_name))
    cat(sprintf("Markers: %s\n", paste(marker_genes, collapse = ", ")))
  }
  
  # Score clusters
  scores <- score_clusters_for_markers(all_markers, marker_genes, celltype_name)
  
  if (nrow(scores) == 0) {
    if (verbose) {
      cat(sprintf("WARNING: No clusters found for %s\n", celltype_name))
    }
    return(list(
      best_cluster = NA,
      all_scores = scores,
      marker_details = data.frame(),
      assignment = data.frame(
        cluster = NA,
        celltype = celltype_name,
        confidence = "none",
        n_markers = 0,
        score = 0
      )
    ))
  }
  
  # Get best match
  best <- scores[1, ]
  
  # Get per-marker details for best cluster
  marker_details <- all_markers %>%
  dplyr::filter(cluster == best$cluster, gene %in% marker_genes) %>%
  dplyr::mutate(
    specificity = pct.1 / (pct.2 + 0.01),
    fold_change = 2^avg_log2FC
  ) %>%
  dplyr::select(gene, avg_log2FC, fold_change, pct.1, pct.2, specificity, p_val_adj) %>%
  dplyr::arrange(desc(avg_log2FC))
  
  # Identify missing markers
  missing_markers <- setdiff(marker_genes, marker_details$gene)
  
  # Determine confidence
  confidence <- case_when(
    # High: ALL markers found AND decent score
    best$n_markers == best$n_markers_total & 
      best$score >= min_score ~ "high",
    
    # Medium: Most markers found AND decent score
    best$n_markers >= min_markers & 
      best$score >= min_score ~ "medium",
    
    # Low: At least 1 marker
    best$n_markers >= 1 ~ "low",
    
    TRUE ~ "none"
  )
  
  if (verbose) {
    cat(sprintf("Found %d / %d markers across %d clusters\n",
               nrow(marker_details),
               length(marker_genes),
               nrow(scores)))
    
    cat(sprintf("\nBest match: Cluster %s (confidence: %s)\n",
               best$cluster, confidence))
    cat(sprintf("  Overall score: %.2f\n", best$score))
    
    # === PER-MARKER DETAILS ===
    cat("\n  Per-marker statistics:\n")
    cat("  ----------------------\n")
    
    for (i in 1:nrow(marker_details)) {
      m <- marker_details[i, ]
      cat(sprintf("  %s:\n", m$gene))
      cat(sprintf("    log2FC: %.2f (%.1f-fold upregulation)\n",
                 m$avg_log2FC, m$fold_change))
      cat(sprintf("    In cluster: %.1f%% of cells\n", m$pct.1 * 100))
      cat(sprintf("    Elsewhere:  %.1f%% of cells\n", m$pct.2 * 100))
      cat(sprintf("    Specificity: %.1fx enriched\n", m$specificity))
      
      # Flag if marker is weak
      if (m$specificity < 2.0) {
        cat("    [!] Low specificity - also common elsewhere\n")
      }
      if (m$pct.1 < 0.5) {
        cat("    [!] Low coverage - only in minority of cluster\n")
      }
      if (m$avg_log2FC < 1.0) {
        cat("    [!] Low fold-change - weak marker\n")
      }
      cat("\n")
    }
    
    if (length(missing_markers) > 0) {
      cat(sprintf("  Missing markers: %s\n", paste(missing_markers, collapse = ", ")))
      cat("  (Not in top differentially expressed genes for any cluster)\n\n")
    }
    
    # Summary statistics
    cat("  Summary:\n")
    cat(sprintf("    Mean log2FC: %.2f (%.1f-fold)\n",
               mean(marker_details$avg_log2FC),
               mean(marker_details$fold_change)))
    cat(sprintf("    Mean coverage: %.1f%% in cluster vs %.1f%% elsewhere\n",
               mean(marker_details$pct.1) * 100,
               mean(marker_details$pct.2) * 100))
    cat(sprintf("    Mean specificity: %.1fx enriched\n",
               mean(marker_details$specificity)))
    
    # Show other matches if they exist
    if (nrow(scores) > 1) {
      cat("\n  Other potential matches:\n")
      for (i in 2:min(3, nrow(scores))) {
        cat(sprintf("    Cluster %s: %s (score: %.2f)\n",
                   scores$cluster[i],
                   scores$markers_found[i],
                   scores$score[i]))
      }
    }
  }
  
  # Create assignment
  assignment <- data.frame(
    cluster = best$cluster,
    celltype = celltype_name,
    confidence = confidence,
    n_markers = best$n_markers,
    n_markers_total = best$n_markers_total,
    score = best$score,
    markers_found = best$markers_found,
    markers_missing = ifelse(length(missing_markers) > 0,
                            paste(missing_markers, collapse = ", "),
                            ""),
    mean_log2FC = best$mean_log2FC,
    mean_pct1 = best$mean_pct1,
    mean_pct2 = best$mean_pct2,
    mean_specificity = best$specificity
  )
  
  return(list(
    best_cluster = best$cluster,
    all_scores = scores,
    marker_details = marker_details,  # NEW: Per-marker details
    assignment = assignment
  ))
}

#' Identify multiple cell types at once
#' 
#' @param all_markers Data frame from FindAllMarkers()
#' @param celltype_markers Named list where names are cell types and values are marker vectors
#'   Example: list(Oligodendrocytes = c("MBP", "MOBP", "PLP1"), ...)
#' @param min_markers Minimum markers for assignment
#' @param min_score Minimum score for assignment
#' @param verbose Print detailed output
#' 
#' @return Data frame with all cell type assignments
#' 
#' @export
sctb_id_alli_celltypes <- function(all_markers,
                                   celltype_markers,
                                   min_markers = 2,
                                   min_score = 5,
                                   verbose = TRUE,
				   print_conflicts = TRUE) {
  
  cat("\n=== Identifying cell types ===\n")
  cat(sprintf("Processing %d cell types...\n", length(celltype_markers)))
  
  all_assignments <- list()
  all_scores_list <- list()
  
  for (celltype in names(celltype_markers)) {
    markers <- celltype_markers[[celltype]]
    
    result <- identify_celltype(
      all_markers,
      markers,
      celltype,
      min_markers = min_markers,
      min_score = min_score,
      verbose = verbose
    )
    
    all_assignments[[celltype]] <- result$assignment
    all_scores_list[[celltype]] <- result$all_scores
  }
  
  # Combine assignments
  assignments_df <- bind_rows(all_assignments)
  
  # Check for conflicts (same cluster assigned to multiple types)
  conflicts <- assignments_df %>%
    group_by(cluster) %>%
    filter(n() > 1, !is.na(cluster)) %>%
    arrange(cluster, desc(score))
  
  if (nrow(conflicts) > 0 && print_conflicts) {  # ← Add check
    cat("\n!!! WARNING: Conflicting assignments !!!\n")
    cat("The following clusters match multiple cell types:\n")

    # Convert to regular dataframe first
    conflicts_to_print <- conflicts %>%
      ungroup() %>%
      select(cluster, celltype, confidence, score) %>%
      #bug? 
      #as.data.frame()
      data.frame()

    print(conflicts_to_print)
    cat("\nReview these clusters manually.\n")
  }
  return(list(
    assignments = assignments_df,
    all_scores = bind_rows(all_scores_list),
    conflicts = conflicts
  ))
}

#' Apply cell type labels and optionally filter removed cell types
#'
#' @param seurat_obj Seurat object with clusters
#' @param celltype_settings Data frame with cluster, celltype, and action columns
#' @param remove_flagged If TRUE, removes cells with action == "remove"
#' @param cluster_col Name of cluster column in metadata (default: "seurat_clusters")
#' @param verbose Print summary
#'
#' @return Seurat object with cell types applied (and optionally filtered)
#'
#' @export
sctb_apply_celltype_labels <- function(seurat_obj,
                                  celltype_settings,
                                  remove_flagged = TRUE,
                                  cluster_col = "seurat_clusters",
                                  verbose = TRUE) {

  # Get cluster assignments
  cluster_ids <- as.character(seurat_obj@meta.data[[cluster_col]])

  # Match to cell type settings
  celltype_idx <- match(cluster_ids, as.character(celltype_settings$cluster))
  seurat_obj$celltype <- celltype_settings$celltype[celltype_idx]
  seurat_obj$celltype_action <- celltype_settings$action[celltype_idx]

  if (verbose) {
    cat("\n=== Cell Type Assignment Summary ===\n")
    celltype_counts <- table(seurat_obj$celltype, seurat_obj$celltype_action)
    print(celltype_counts)
  }

  # Filter if requested
  if (remove_flagged) {
    cells_before <- ncol(seurat_obj)

    # Identify cells to remove
    remove_celltypes <- celltype_settings$celltype[celltype_settings$action == "remove"]
    cells_to_keep <- !seurat_obj$celltype %in% remove_celltypes

    # Subset
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])

    cells_after <- ncol(seurat_obj)
    cells_removed <- cells_before - cells_after

    if (verbose) {
      cat(sprintf("\nRemoved %d cells (%.1f%%) with flagged cell types:\n",
                 cells_removed, 100 * cells_removed / cells_before))
      cat(sprintf("  - %s\n", paste(remove_celltypes, collapse = ", ")))
      cat(sprintf("Retained %d cells across %d cell types\n",
                 cells_after, length(unique(seurat_obj$celltype))))
    }
  }

  # Set Idents to celltype
  Idents(seurat_obj) <- seurat_obj$celltype

  return(seurat_obj)
}
