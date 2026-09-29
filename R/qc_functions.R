#' Check PC correlations with technical covariates
#' @param seurat_obj Merged Seurat object with PCA reduction
#' @param n_pcs Number of PCs to test (default 10)
#' @param reduction Name of reduction to use (default "pca")
#' @return ggplot object showing correlation heatmap
sctb_pc_technical_bias <- function(seurat_obj, n_pcs = 10, reduction = "pca") {
  # Technical covariates to test
  covariates <- c("percent.mt", "nCount_RNA", "nFeature_RNA", 
                  "TSS.enrichment", "nCount_ATAC", "nFeature_ATAC")
  
  # Check which covariates exist
  available_covariates <- covariates[covariates %in% colnames(seurat_obj@meta.data)]
  
  # Get PC embeddings
  pcs <- Embeddings(seurat_obj, reduction = reduction)[, 1:n_pcs]
  
  # Compute correlations
  cor_matrix <- matrix(NA, nrow = n_pcs, ncol = length(available_covariates))
  rownames(cor_matrix) <- paste0("PC", 1:n_pcs)
  colnames(cor_matrix) <- available_covariates
  
  for (i in 1:n_pcs) {
    for (j in seq_along(available_covariates)) {
      cov_name <- available_covariates[j]
      cor_matrix[i, j] <- cor(pcs[, i], seurat_obj@meta.data[[cov_name]], 
                              use = "complete.obs")
    }
  }
  
  # Melt for ggplot
  cor_df <- melt(cor_matrix)
  colnames(cor_df) <- c("PC", "Covariate", "Correlation")
  
  # Create heatmap
  p1 <- ggplot(cor_df, aes(x = Covariate, y = PC, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("PC Correlation with Technical Covariates -", 
                       seurat_obj@project.name),
         x = "Technical Covariate", y = "Principal Component")
  
  # Create line plot (absolute correlations)
  cor_df$AbsCorrelation <- abs(cor_df$Correlation)
  p2 <- ggplot(cor_df, aes(x = PC, y = AbsCorrelation, 
                           color = Covariate, group = Covariate)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(title = "Absolute Correlation by PC (Technical Bias Decay)",
         x = "Principal Component", y = "Absolute Correlation",
         color = "Covariate") +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
    annotate("text", x = n_pcs * 0.8, y = 0.12, 
             label = "Threshold (|r| = 0.1)", color = "gray50")
  
  # Print summary stats
  cat("\n=== PC1 Correlations (Technical Bias Check) ===\n")
  pc1_cors <- cor_matrix[1, ]
  print(sort(abs(pc1_cors), decreasing = TRUE))
  
  cat("\n=== PC2 Correlations (Should be lower) ===\n")
  pc2_cors <- cor_matrix[2, ]
  print(sort(abs(pc2_cors), decreasing = TRUE))
  
  cat("\n=== Maximum Absolute Correlation by PC ===\n")
  max_cors <- apply(abs(cor_matrix), 1, max)
  print(max_cors)
  
  # Return both plots
  list(heatmap = p1, lineplot = p2, cor_matrix = cor_matrix)
}

#' Check reduction-component correlations with technical covariates
#'
#' Reduction-agnostic replacement for check_pc_technical_bias(). Works on any
#' dimensional reduction (pca, lsi, harmony, ...). Two things make it modality-
#' neutral where the old function was RNA-biased:
#'
#'   1. Component labels are taken from the reduction's OWN embedding column
#'      names, so an LSI run is labeled LSI_* and a PCA run PC_*. The plot can no
#'      longer mislabel LSI components as "PC1, PC2, ..." the way the hardcoded
#'      paste0("PC", ...) did.
#'
#'   2. The summary drops the "PC2 should be lower" editorializing. That assumed
#'      every covariate behaves like RNA depth -- peak on component 1, decay
#'      after -- which is only true for a covariate matched to the reduction's
#'      dominant technical axis. A guest covariate (nucleosome_signal on a PCA
#'      reduction, say) can legitimately peak anywhere. It is replaced by a
#'      "peak component per covariate" table that reports where each technical
#'      axis actually lives instead of presuming component 1.
#'
#' Covariates are tested as a MIXED RNA+ATAC set on purpose: a technical depth
#' axis loads on one modality's depth and is flat elsewhere; a shared biological
#' axis loads on both modalities' complexity together; the near-zero cross-modal
#' cells are the negative control that distinguishes the two. Absent covariates
#' are dropped by the availability filter, so projects that never computed (e.g.)
#' nucleosome_signal still run.
#'
#' @param seurat_obj      Seurat object with the requested reduction computed.
#' @param n_comps         Number of components to test (clamped to what exists).
#' @param reduction       Name of reduction to use (default "pca").
#' @param flag_threshold  |r| at/above which a component is flagged to inspect.
#' @param line_threshold  Dashed reference line on the decay lineplot.
#' @return (invisibly) list(heatmap, lineplot, cor_matrix)
sctb_technical_bias_both <- function(seurat_obj, n_comps = 10,
                                           reduction = "pca",
                                           flag_threshold = 0.3,
                                           line_threshold = 0.1) {
  # Mixed RNA + ATAC covariates by design (see header). Add ATAC-specific
  # quality proxies here (pct_reads_in_peaks, blacklist_ratio) if/when computed;
  # the availability filter keeps them optional.
  covariates <- c("percent.mt", "nCount_RNA", "nFeature_RNA", "nucleosome_signal",
                  "TSS.enrichment", "nCount_ATAC", "nFeature_ATAC")
  available_covariates <- covariates[covariates %in% colnames(seurat_obj@meta.data)]
  if (length(available_covariates) == 0L) {
    stop("check_reduction_technical_bias: none of the expected covariates are ",
         "present in meta.data.")
  }

  # Pull embeddings; let component labels come from the DATA, not a hardcoded
  # prefix. Clamp n_comps to what the reduction actually has, and say so loudly
  # rather than letting the [, 1:n] index silently error or overrun.
  emb <- Embeddings(seurat_obj, reduction = reduction)
  n_avail <- ncol(emb)
  if (n_comps > n_avail) {
    warning(sprintf(
      "check_reduction_technical_bias: requested %d components but reduction '%s' ",
      n_comps, reduction),
      sprintf("has only %d; clamping to %d.", n_avail, n_avail))
    n_comps <- n_avail
  }
  emb <- emb[, seq_len(n_comps), drop = FALSE]
  comp_labels <- colnames(emb)                       # e.g. "PC_1" / "LSI_1"

  # Short label for titles: strip the trailing underscore off the reduction key,
  # falling back to the (upper-cased) reduction name if the key is unset.
  red_key <- tryCatch(Key(seurat_obj[[reduction]]), error = function(e) "")
  reduction_label <- if (nzchar(red_key)) sub("_$", "", red_key) else toupper(reduction)

  # Correlation matrix: components x covariates
  cor_matrix <- matrix(NA_real_, nrow = n_comps, ncol = length(available_covariates),
                       dimnames = list(comp_labels, available_covariates))
  for (i in seq_len(n_comps)) {
    for (j in seq_along(available_covariates)) {
      cor_matrix[i, j] <- cor(emb[, i], seurat_obj@meta.data[[available_covariates[j]]],
                              use = "complete.obs")
    }
  }

  # Melt for ggplot; fix factor levels to row order so component 10 doesn't sort
  # before component 2.
  cor_df <- melt(cor_matrix)
  colnames(cor_df) <- c("Component", "Covariate", "Correlation")
  cor_df$Component <- factor(cor_df$Component, levels = comp_labels)

  # Heatmap
  p1 <- ggplot(cor_df, aes(x = Covariate, y = Component, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1)) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(reduction_label,
                        " Correlation with Technical Covariates - ",
                        seurat_obj@project.name),
         x = "Technical Covariate", y = paste(reduction_label, "Component"))

  # Decay lineplot (absolute correlation). linewidth (not size) for current
  # ggplot2; if on <3.4 change linewidth -> size.
  cor_df$AbsCorrelation <- abs(cor_df$Correlation)
  p2 <- ggplot(cor_df, aes(x = Component, y = AbsCorrelation,
                           color = Covariate, group = Covariate)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(reduction_label, "- Absolute Correlation by Component"),
         x = paste(reduction_label, "Component"), y = "Absolute Correlation",
         color = "Covariate") +
    geom_hline(yintercept = line_threshold, linetype = "dashed", color = "gray50") +
    annotate("text", x = n_comps * 0.8, y = line_threshold + 0.02,
             label = sprintf("Threshold (|r| = %.2f)", line_threshold),
             color = "gray50")

  # ---- Modality-neutral summary ----
  abs_mat <- abs(cor_matrix)

  cat(sprintf("\n=== %s: Component 1 correlations (sorted |r|) ===\n", reduction_label))
  print(round(sort(abs_mat[1, ], decreasing = TRUE), 3))

  cat(sprintf("\n=== %s: Max |r| by component ===\n", reduction_label))
  print(round(apply(abs_mat, 1, max), 3))

  # Where each covariate PEAKS -- the modality-neutral replacement for
  # "component 2 should be lower". Tells you where each technical axis actually
  # lives instead of assuming it lives on component 1.
  cat(sprintf("\n=== %s: Peak component per covariate ===\n", reduction_label))
  peak_tbl <- data.frame(
    covariate      = available_covariates,
    peak_component = comp_labels[apply(abs_mat, 2, which.max)],
    peak_abs_r     = round(apply(abs_mat, 2, max), 3),
    row.names      = NULL,
    stringsAsFactors = FALSE
  )
  peak_tbl <- peak_tbl[order(-peak_tbl$peak_abs_r), ]
  print(peak_tbl, row.names = FALSE)

  # Components worth inspecting: any at/above flag_threshold on any covariate,
  # named with the responsible covariate and signed r.
  max_by_comp <- apply(abs_mat, 1, max)
  flagged_idx <- which(max_by_comp >= flag_threshold)
  cat(sprintf("\n=== %s: Components with |r| >= %.2f (inspect, do not auto-drop) ===\n",
              reduction_label, flag_threshold))
  if (length(flagged_idx) == 0L) {
    cat("  none\n")
  } else {
    for (i in flagged_idx) {
      jmax <- which.max(abs_mat[i, ])
      cat(sprintf("  %-8s  top: %-18s r = %+.2f\n",
                  comp_labels[i], available_covariates[jmax], cor_matrix[i, jmax]))
    }
  }

  invisible(list(heatmap = p1, lineplot = p2, cor_matrix = cor_matrix))
}

#' Run parameter sweep and display plots for visual inspection
#' @param seurat_obj Harmonized Seurat object
#' @param dims_range List of dimension ranges
#' @param knn_values Vector of k values
#' @param res_values Vector of resolution values
#' @param alg Clustering algorithm
#' @param plots boolean (default: TRUE) of whether to create dimplots during the sweep
#' @param typing boolean (default: FALSE) of whether to do a quick typing check during sweep
#' @param cluster_seed Random seed
#' @return Data frame with basic results (for reference only)
sctb_param_sweep <- function(seurat_obj, dims_range, knn_values, res_values, alg,
                                     plots = TRUE, typing = FALSE, cluster_seed) {

  n_combos <- length(dims_range) * length(knn_values) * length(res_values)
  cat(sprintf("\n=== Parameter Sweep: %d combinations ===\n", n_combos))
  if (plots) {
    cat("Plots will be displayed in browser\n\n")
  }

  total_cells <- ncol(seurat_obj)
  results <- list()
  counter <- 1

  # OUTER LOOP: dims + knn (run FMMN once per combination)
  for (dims_idx in seq_along(dims_range)) {
    dims <- dims_range[[dims_idx]]
    dims_min <- min(dims)
    dims_max <- max(dims)
    dims_str <- sprintf("%d-%d", dims_min, dims_max)

      for (knn in knn_values) {
        cat(sprintf("\n=== FMMN: dims=%s, knn=%d ===\n",
                   gsub("-", ":", dims_str), knn))

        # NOW passing dims to FMMN_task!
        obj_fmmn <- FMMN_task(seurat_obj, knn = knn, dims = dims)

      # INNER LOOP: resolutions (reuse FMMN result)
      for (res in res_values) {
        cat(sprintf("Clustering loop: [%d/%d] dims=%s, knn=%d, res=%.3f ... ",
                   counter, n_combos, dims_str, knn, res))

        # Copy FMMN result and cluster
        obj_clustered <- obj_fmmn
        obj_clustered <- cluster_data(
          obj_clustered,
          alg = alg,
          res = res,
          cluster_seed = cluster_seed,
          singleton_handling = "keep",
          run_umap = TRUE  # Need UMAP for plotting
        )

        # Get cluster info
        cluster_assignments <- obj_clustered@meta.data$seurat_clusters
        singleton_count <- sum(cluster_assignments == "singleton")
        cluster_count <- length(setdiff(unique(cluster_assignments), "singleton"))

        cat(sprintf("%d clusters (%d singletons)\n", cluster_count, singleton_count))

        # ---------------------------------------------------------------
        # TYPING BLOCK
        # ---------------------------------------------------------------
        if (typing) {
          typing_result <- assign_celltype_from_dotplot(
            seurat_obj   = obj_clustered,
            marker_lists = tissue_markers$marker_lists
          )

          # Merge parameter metadata into assignments
          sweep_row <- typing_result$assignments %>%
            dplyr::rename(cluster_id = cluster) %>%
            dplyr::mutate(
              sweep_id    = counter,
              dims_min    = dims_min,
              dims_max    = dims_max,
              knn         = knn,
              res         = res,
              n_clusters  = cluster_count,
              n_singletons = singleton_count,
              pct_of_total = n_cells / total_cells * 100
            ) %>%
            dplyr::select(
              sweep_id, dims_min, dims_max, knn, res,
              n_clusters, n_singletons,
              cluster_id, n_cells, pct_of_total,
              mean_exp, mean_pct, score, celltype
            )

          # Write: overwrite + header on first sweep, append thereafter
          write.table(
            sweep_row,
            file      = param_sweep_typing_file,
            append    = (counter > 1),
            sep       = ",",
            row.names = FALSE,
            col.names = (counter == 1),
            quote     = TRUE
          )

          #this should be percent of data assigned, not clusters, or should just not exist as
          #an inline report at all, as comparing number of clusters assigned isn't very meaningful.
          #e.g. could have 12 great clusters covering 98% of the data, with 50 tiny unassigned.
          cat(sprintf("  Typed: %d clusters assigned, %d unassigned\n",
                     sum(sweep_row$celltype != "Unassigned"),
                     sum(sweep_row$celltype == "Unassigned")))
        }

        # Store basic sweep results
        results[[counter]] <- list(
          dims_min = dims_min,
          dims_max = dims_max,
          dims_str = dims_str,
          knn = knn,
          resolution = res,
          n_clusters = cluster_count,
          n_singletons = singleton_count
        )

        if (plots) {
          plot_title <- sprintf("dims=%s, knn=%d, res=%.3f\n%d clusters (%d singletons)",
                             dims_str, knn, res, cluster_count, singleton_count)

          p <- DimPlot(obj_clustered, reduction = "wnn.umap",
                    group.by = "seurat_clusters", label = TRUE, raster = FALSE) +
            ggtitle(plot_title)

          print(p)
          rm(p)
          }

        # Clean up
        rm(obj_clustered)
        gc(verbose = FALSE)

        counter <- counter + 1
      }

      # Clean up FMMN result after all resolutions
      # R auto garbage collection isn't great, bleed happens, make explicit
      rm(obj_fmmn)
      gc(verbose = FALSE, full = TRUE)
    }
  }

  # Convert results to data frame for reference
  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      dims_min = x$dims_min,
      dims_max = x$dims_max,
      knn = x$knn,
      resolution = x$resolution,
      n_clusters = x$n_clusters,
      n_singletons = x$n_singletons,
      stringsAsFactors = FALSE
    )
  }))

  return(results_df)
}

#' Compare typing results across parameter sweep combinations
#'
#' Reads param_sweep_typing_file (written by run_parameter_sweep with
#' typing = TRUE), filters out small clusters, and summarizes at the
#' celltype level per sweep combination. Produces one row per
#' sweep_id x celltype with weighted score, total cells, and cluster
#' count — giving a single view of which parameter set best identified
#' each population.
#'
#' @param min_pct_of_total Minimum percent of total cells for a cluster
#'   to be included in the comparison (default 0.25, i.e. 250 cells in
#'   a 100k dataset). Clusters below this are excluded from the summary
#'   but were still typed — this is a reporting filter, not a typing
#'   filter.
#' @param typing_file Path to the sweep typing CSV. Defaults to
#'   param_sweep_typing_file from system_settings.R.
#' @return Data frame with one row per sweep_id x celltype, sorted by
#'   sweep_id then descending pct_of_total_assigned.
sctb_comp_typing_results <- function(min_pct_of_total = 0.25,
                                   typing_file = param_sweep_typing_file) {

  if (!file.exists(typing_file)) {
    stop(sprintf(
      "Typing results file not found: %s\n",
      "Run run_parameter_sweep(typing = TRUE) first.",
      typing_file
    ))
  }

  raw <- read.csv(typing_file, stringsAsFactors = FALSE)

  # Filter small clusters and Unassigned before aggregating
  filtered <- raw %>%
    dplyr::filter(
      pct_of_total >= min_pct_of_total,
      celltype != "Unassigned"
    )

  if (nrow(filtered) == 0) {
    warning(sprintf(
      "No clusters pass min_pct_of_total = %.2f%%. ",
      "Try lowering the threshold.",
      min_pct_of_total
    ))
    return(invisible(data.frame()))
  }

  # Aggregate to sweep_id x celltype
  summary <- filtered %>%
    dplyr::group_by(sweep_id, dims_min, dims_max, knn, res, celltype) %>%
    dplyr::summarize(
      n_clusters_assigned  = dplyr::n(),
      n_cells_total        = sum(n_cells),
      pct_of_total_assigned = sum(pct_of_total),
      weighted_score       = sum(score * n_cells) / sum(n_cells),
      mean_pct_expressed   = stats::weighted.mean(mean_pct, n_cells),
      .groups = "drop"
    ) %>%
    dplyr::arrange(sweep_id, dplyr::desc(pct_of_total_assigned))

  # Also compute pct_unassigned per sweep for context
  unassigned_summary <- raw %>%
    dplyr::group_by(sweep_id) %>%
    dplyr::summarize(
      total_clusters      = dplyr::n(),
      n_unassigned        = sum(celltype == "Unassigned"),
      pct_unassigned_clusters = mean(celltype == "Unassigned") * 100,
      .groups = "drop"
    )

  summary <- summary %>%
    dplyr::left_join(unassigned_summary, by = "sweep_id")

  cat(sprintf("\n=== Typing Comparison: %d sweep combinations ===\n",
             dplyr::n_distinct(summary$sweep_id)))
  cat(sprintf("Min cluster size filter: %.2f%% of total cells\n", min_pct_of_total))
  cat(sprintf("Cell types represented: %s\n\n",
             paste(sort(unique(summary$celltype)), collapse = ", ")))

  return(summary)
}

# TODO: description
stb_basic_elbow <- function(seurat_obj) {
  pc_sd <- seurat_obj[["pca"]]@stdev
  pc_table <- data.frame(PC = 1:length(pc_sd), StandardDeviation = pc_sd)

  # Custom ElbowPlot with grid
  ggplot(pc_table, aes(x = PC, y = StandardDeviation)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    scale_x_continuous(breaks = 1:length(pc_sd)) +
    scale_y_continuous(breaks = pretty(pc_sd)) +
    labs(title = "Elbow Plot of PCs", x = "Principal Component", y = "Standard Deviation") +
    theme(panel.grid.major = element_line(color = "gray", size = 0.5),
          panel.grid.minor = element_line(color = "lightgray", size = 0.2))
}

#TODO: description
#QC density plots against ATAC assay
sctb_QCDensity_ATAC <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "ATAC"
  nCA_TSS_density <- DensityScatter(seurat_obj,
                                    x = 'nCount_ATAC',
                                    y = 'TSS.enrichment',
                                    log_x = TRUE,
                                    quantiles = TRUE
  ) +
    labs(
      title = paste("ATAC counts + TSS density plot for ",
                    seurat_obj@project.name,
                    " for debris and doublets")
    ) +
    stat_density_2d(
      color = "red",
      linewidth = 1,   # updated from size
      bins = 6,
      contour = TRUE
    )
  nFA_TSS_density <- DensityScatter(seurat_obj,
                                    x = 'nFeature_ATAC',
                                    y = 'TSS.enrichment',
                                    log_x = TRUE,
                                    quantiles = TRUE
  ) +
    labs(
    title = paste("ATAC features + TSS density plot for ",
                  seurat_obj@project.name,
                  " for debris and doublets")) + 
    stat_density_2d(
      color = "red",
      linewidth = 1,
      bins = 6, # Increase bins for smoother contours
      contour = TRUE
    )
  nss_TSS_density <- DensityScatter(seurat_obj,
                                    x = 'nucleosome_signal',
                                    y = 'TSS.enrichment',
                                    quantiles = TRUE
  ) +
    labs(
    title = paste("nucleosome signal vs TSS enrichment density plot for ",
                  seurat_obj@project.name,
                  " to find nucleosome-rich cells")) + 
    stat_density_2d(
      color = "red",
      linewidth = 1,
      bins = 6, # Increase bins for smoother contours
      contour = TRUE
    )

  print(nCA_TSS_density)
  print(nFA_TSS_density)
  print(nss_TSS_density)
}

sctb_QCDensity_RNA <- function(seurat_obj) {
  #DefaultAssay(seurat_obj) <- "RNA"
  pMT_nCR_density <- DensityScatter(seurat_obj,
                                    x = 'percent.mt',
                                    y = 'nCount_RNA',
                                    log_y = TRUE,
                                    quantiles = c(5, 10, 50, 90, 95)
  ) +
    labs(
      title = paste("percent.mt + RNA counts density plot for ",
                  seurat_obj@project.name,
                  " for cells with high mitochondrial content")) + 
    stat_density_2d(
      color = "red",
      linewidth = 1,
      bins = 10, # Increase bins for smoother contours
      contour = TRUE
    )
  pMT_nFR_density <- DensityScatter(seurat_obj,
                                    x = 'percent.mt',
                                    y = 'nFeature_RNA',
                                    log_y = TRUE,
                                    quantiles = c(5, 10, 50, 90, 95)
  ) +
    labs(
    title = paste("percent.mt + RNA features density plot for ",
                  seurat_obj@project.name,
                  " for cells with high mitochondrial content")) + 
    stat_density_2d(
      color = "red",
      linewidth = 1,
      bins = 10, # Increase bins for smoother contours
      contour = TRUE
    )
  nCR_nFR_density <- DensityScatter(seurat_obj,
                                    x = 'nCount_RNA',
                                    y = 'nFeature_RNA',
                                    log_x = TRUE,
                                    log_y = TRUE,
                                    quantiles = TRUE
  ) +
    ggtitle(paste("RNA count vs RNA feature density plot for ",
                  seurat_obj@project.name,
                  " for doublets and general QC"))
  print(nCR_nFR_density)
  print(pMT_nCR_density)
  print(pMT_nFR_density)
}

#QC Vln plots against ATAC assay
sctb_QCVlnA <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "ATAC"
  VPsn2 <- VlnPlot(
    object = seurat_obj,
    features = c("nCount_ATAC", "nCount_RNA", "nFeature_ATAC", "nFeature_RNA",
                 "TSS.enrichment", "nucleosome_signal", "percent.mt"),
    ncol = 4,
    pt.size = 0
  ) 
  print(VPsn2)
}

#QC Violin plots against RNA assay
sctb_QCVlnR <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  VP2sn2 <- VlnPlot(seurat_obj, features = "percent.mt") +
    ggtitle(paste("MT Density plot for", seurat_obj@project.name))
  print(VP2sn2)
}
