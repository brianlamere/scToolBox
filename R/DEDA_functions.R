#' Run pseudo-bulk differential expression
#' 
#' Aggregates RNA counts by sample, then tests at sample level.
#' Statistically appropriate for multi-sample studies.
#' 
#' @param seurat_obj Seurat object with RNA assay
#' @param celltype_col Column name for cell types
#' @param celltype_value Cell type to test (e.g., "Oligodendrocytes")
#' @param sample_col Column name for sample IDs (e.g., "orig.ident")
#' @param group_col Column name for experimental groups
#' @param ident.1 First group
#' @param ident.2 Second group
#' @param min_cells_per_sample Minimum cells per sample to include
#' 
#' @return Data frame with DE results
sctb_pseudobulk_DE <- function(seurat_obj,
                              celltype_col = "celltype",
                              celltype_value,
                              sample_col = "orig.ident",
                              group_col = "group",
                              ident.1,
                              ident.2,
                              min_cells_per_sample = 50) {
  
  # Subset to cell type
  cell_subset <- subset(seurat_obj, subset = !!as.name(celltype_col) == celltype_value)
  
  # Further subset to groups of interest
  cells_to_keep <- cell_subset@meta.data[[group_col]] %in% c(ident.1, ident.2)
  cell_subset <- subset(cell_subset, cells = colnames(cell_subset)[cells_to_keep])
  
  cat(sprintf("  %s: %d cells across %d samples\n",
             celltype_value, ncol(cell_subset),
             length(unique(cell_subset@meta.data[[sample_col]]))))
  
  # Get RNA counts
  DefaultAssay(cell_subset) <- "RNA"
  
  # Get counts (handle Seurat v5)
  counts <- tryCatch({
    GetAssayData(cell_subset, assay = "RNA", layer = "counts")
  }, error = function(e) {
    cell_subset[["RNA"]]@counts
  })
  
  # Create sample metadata
  sample_meta <- cell_subset@meta.data %>%
    group_by(!!sym(sample_col)) %>%
    summarize(
      group = first(!!sym(group_col)),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    filter(n_cells >= min_cells_per_sample)
  
  # Check if we have enough samples
  n_group1 <- sum(sample_meta$group == ident.1)
  n_group2 <- sum(sample_meta$group == ident.2)
  
  if (n_group1 < 2 || n_group2 < 2) {
    warning(sprintf("  Insufficient samples for %s: %s (n=%d) vs %s (n=%d). Need ≥2 per group.",
                   celltype_value, ident.1, n_group1, ident.2, n_group2))
    return(data.frame())
  }
  
  cat(sprintf("  Aggregating to %d samples: %s (n=%d) vs %s (n=%d)\n",
             nrow(sample_meta), ident.1, n_group1, ident.2, n_group2))
  
  # Aggregate counts by sample (pseudo-bulk)
  pseudobulk_counts <- matrix(0, nrow = nrow(counts), ncol = nrow(sample_meta))
  rownames(pseudobulk_counts) <- rownames(counts)
  colnames(pseudobulk_counts) <- sample_meta[[sample_col]]
  
  for (sample in sample_meta[[sample_col]]) {
    sample_cells <- colnames(cell_subset)[cell_subset@meta.data[[sample_col]] == sample]
    if (length(sample_cells) > 0) {
      pseudobulk_counts[, sample] <- Matrix::rowSums(counts[, sample_cells, drop = FALSE])
    }
  }
  
  # Filter genes: must be expressed in at least 1 sample per group
  group1_samples <- sample_meta[[sample_col]][sample_meta$group == ident.1]
  group2_samples <- sample_meta[[sample_col]][sample_meta$group == ident.2]
  
  genes_group1 <- rowSums(pseudobulk_counts[, group1_samples, drop = FALSE] > 0) > 0
  genes_group2 <- rowSums(pseudobulk_counts[, group2_samples, drop = FALSE] > 0) > 0
  genes_to_test <- genes_group1 | genes_group2
  
  pseudobulk_counts <- pseudobulk_counts[genes_to_test, ]
  
  cat(sprintf("  Testing %d expressed genes\n", nrow(pseudobulk_counts)))
  
  # Create DESeq2 dataset
  coldata <- data.frame(
    sample = sample_meta[[sample_col]],
    group = factor(sample_meta$group, levels = c(ident.2, ident.1)),
    row.names = sample_meta[[sample_col]]
  )
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(pseudobulk_counts),
    colData = coldata,
    design = ~ group
  )
  
  # Run DESeq2
  cat("  Running DESeq2...\n")
  dds <- DESeq(dds, fitType = "local", quiet = TRUE)
  
  # Get results
  res <- results(dds, contrast = c("group", ident.1, ident.2))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Reorder columns
  res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE",
                       "stat", "pvalue", "padj")]
  
  # Sort by adjusted p-value
  res_df <- res_df[order(res_df$padj), ]
  
  n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  cat(sprintf("  → %d DE genes at FDR < 0.05 (%.1f%% of tested)\n",
             n_sig, 100 * n_sig / nrow(res_df)))
  
  return(res_df)
}

#' Wrapper for pseudo-bulk DE with export
run_pseudobulk_DE_and_export <- function(seurat_obj,
                                         celltype_col = "celltype",
                                         celltype_value,
                                         sample_col = "orig.ident",
                                         group_col = "group",
                                         ident.1,
                                         ident.2,
                                         output_prefix = "DiffExpress_results",
                                         min_cells_per_sample = 50) {
  
  # Run pseudo-bulk DE
  results <- run_pseudobulk_DE(
    seurat_obj = seurat_obj,
    celltype_col = celltype_col,
    celltype_value = celltype_value,  # ← Changed
    sample_col = sample_col,
    group_col = group_col,
    ident.1 = ident.1,
    ident.2 = ident.2,
    min_cells_per_sample = min_cells_per_sample
  )
  
  # Check if we got results
  if (nrow(results) == 0) {
    cat("  No results (insufficient samples)\n\n")
    return(invisible(NULL))
  }
  
  # Save results
  fname <- sprintf("%s_pseudobulk_%s_%s_vs_%s.csv",
                  output_prefix, celltype_value, ident.1, ident.2)  # ← Changed
  write.csv(results, fname, row.names = FALSE)
  
  cat(sprintf("  Saved: %s\n\n", basename(fname)))
  
  invisible(results)
}

#' Run pseudo-bulk differential accessibility
run_pseudobulk_DA <- function(seurat_obj,
                              celltype_col = "celltype",
                              celltype_value,
                              sample_col = "orig.ident",
                              group_col = "group",
                              ident.1,
                              ident.2,
                              min_cells_per_sample = 50) {
  
  # Subset to cell type
  cell_subset <- subset(seurat_obj, subset = !!as.name(celltype_col) == celltype_value)
  
  # Further subset to groups of interest
  cells_to_keep <- cell_subset@meta.data[[group_col]] %in% c(ident.1, ident.2)
  cell_subset <- subset(cell_subset, cells = colnames(cell_subset)[cells_to_keep])
  
  cat(sprintf("  %s: %d cells across %d samples\n",
             celltype_value, ncol(cell_subset),
             length(unique(cell_subset@meta.data[[sample_col]]))))
  
  # Get ATAC counts
  DefaultAssay(cell_subset) <- "ATAC"
  
  counts <- tryCatch({
    GetAssayData(cell_subset, assay = "ATAC", layer = "counts")
  }, error = function(e) {
    tryCatch({
      GetAssayData(cell_subset, assay = "ATAC", slot = "counts")
    }, error = function(e2) {
      cell_subset[["ATAC"]]@counts
    })
  })
  
  # Create sample metadata
  sample_meta <- cell_subset@meta.data %>%
    group_by(!!sym(sample_col)) %>%
    summarize(
      group = first(!!sym(group_col)),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    filter(n_cells >= min_cells_per_sample)
  
  # Check if we have enough samples
  n_group1 <- sum(sample_meta$group == ident.1)
  n_group2 <- sum(sample_meta$group == ident.2)
  
  if (n_group1 < 2 || n_group2 < 2) {
    warning(sprintf("  Insufficient samples for %s: %s (n=%d) vs %s (n=%d). Need ≥2 per group.",
                   celltype_value, ident.1, n_group1, ident.2, n_group2))
    return(data.frame())
  }
  
  cat(sprintf("  Aggregating to %d samples: %s (n=%d) vs %s (n=%d)\n",
             nrow(sample_meta), ident.1, n_group1, ident.2, n_group2))
  
  # Aggregate counts by sample (pseudo-bulk)
  pseudobulk_counts <- matrix(0, nrow = nrow(counts), ncol = nrow(sample_meta))
  rownames(pseudobulk_counts) <- rownames(counts)
  colnames(pseudobulk_counts) <- sample_meta[[sample_col]]
  
  for (sample in sample_meta[[sample_col]]) {
    sample_cells <- colnames(cell_subset)[cell_subset@meta.data[[sample_col]] == sample]
    if (length(sample_cells) > 0) {
      pseudobulk_counts[, sample] <- Matrix::rowSums(counts[, sample_cells, drop = FALSE])
    }
  }
  
  # Filter peaks: must be accessible in at least 1 sample per group
  group1_samples <- sample_meta[[sample_col]][sample_meta$group == ident.1]
  group2_samples <- sample_meta[[sample_col]][sample_meta$group == ident.2]
  
  peaks_group1 <- rowSums(pseudobulk_counts[, group1_samples, drop = FALSE] > 0) > 0
  peaks_group2 <- rowSums(pseudobulk_counts[, group2_samples, drop = FALSE] > 0) > 0
  peaks_to_test <- peaks_group1 | peaks_group2
  
  pseudobulk_counts <- pseudobulk_counts[peaks_to_test, ]
  
  cat(sprintf("  Testing %d accessible peaks\n", nrow(pseudobulk_counts)))
  
  # Create DESeq2 dataset
  coldata <- data.frame(
    sample = sample_meta[[sample_col]],
    group = factor(sample_meta$group, levels = c(ident.2, ident.1)),
    row.names = sample_meta[[sample_col]]
  )
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(pseudobulk_counts),
    colData = coldata,
    design = ~ group
  )
  
  # Run DESeq2
  cat("  Running DESeq2...\n")
  dds <- DESeq(dds, quiet = TRUE)
  
  # Get results
  res <- results(dds, contrast = c("group", ident.1, ident.2))
  res_df <- as.data.frame(res)
  res_df$peak <- rownames(res_df)
  
  # Reorder columns
  res_df <- res_df[, c("peak", "baseMean", "log2FoldChange", "lfcSE",
                       "stat", "pvalue", "padj")]
  
  # Sort by adjusted p-value
  res_df <- res_df[order(res_df$padj), ]
  
  n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  cat(sprintf("  → %d DA peaks at FDR < 0.05 (%.1f%% of tested)\n",
             n_sig, 100 * n_sig / nrow(res_df)))
  
  return(res_df)
}

#' Wrapper for pseudo-bulk DA with export
run_pseudobulk_DA_and_export <- function(seurat_obj,
                                         celltype_col = "celltype",
                                         celltype_value,
                                         sample_col = "orig.ident",
                                         group_col = "group",
                                         ident.1,
                                         ident.2,
                                         output_prefix = "DiffAccess_results",
                                         min_cells_per_sample = 50) {
  
  # Run pseudo-bulk DA
  results <- run_pseudobulk_DA(
    seurat_obj = seurat_obj,
    celltype_col = celltype_col,
    celltype_value = celltype_value,  # ← Changed
    sample_col = sample_col,
    group_col = group_col,
    ident.1 = ident.1,
    ident.2 = ident.2,
    min_cells_per_sample = min_cells_per_sample
  )
  
  # Check if we got results
  if (nrow(results) == 0) {
    cat("  No results (insufficient samples)\n")
    return(invisible(NULL))
  }
  
  # Save results
  fname <- sprintf("%s_pseudobulk_%s_%s_vs_%s.csv",
                  output_prefix, celltype_value, ident.1, ident.2)  # ← Changed
  write.csv(results, fname, row.names = FALSE)
  
  cat(sprintf("  Saved: %s\n\n", basename(fname)))
  
  invisible(results)
}
