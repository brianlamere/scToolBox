#' Read CellBender-corrected RNA counts from H5 file
#'
#' @param sample Sample identifier
#' @return dgCMatrix of CellBender RNA counts (genes × barcodes)
sctb_read_cellbender_rna_counts <- function(sample) {
  cb_h5_path <- file.path(cb_datadir, sample, cellbender_rna_h5filename)
  
  if (!file.exists(cb_h5_path)) {
    stop(sprintf("CellBender H5 file not found for sample '%s': %s", sample, cb_h5_path))
  }
  
  # Check if scCustomize is available
  if (!requireNamespace("scCustomize", quietly = TRUE)) {
    stop("Package 'scCustomize' is required for reading CellBender output. Install with: install.packages('scCustomize')")
  }
  
  cat(sprintf("Reading CellBender RNA counts from: %s\n", cb_h5_path))
  cb_rna <- scCustomize::Read_CellBender_h5_Mat(file_name = cb_h5_path)
  
  cat(sprintf("Loaded CellBender RNA: %d genes, %d cells\n", nrow(cb_rna), ncol(cb_rna)))
  
  return(cb_rna)
}

#' Splice CellBender RNA into multiome data by intersecting barcodes
#'
#' @param rna_orig CellRanger RNA counts matrix
#' @param atac_orig CellRanger ATAC counts matrix
#' @param cb_rna CellBender RNA counts matrix
#' @param sample Sample name (for error messages)
#' @return Named list with rna_counts, atac_counts, common_cells, in_mat, cb_mat2
sctb_splice_cellbender_rna_into_multiome <- function(rna_orig, atac_orig, cb_rna, sample) {
  # Get barcode sets
  atac_cells <- colnames(atac_orig)
  cb_cells <- colnames(cb_rna)
  
  # Intersect barcodes
  common_cells <- intersect(atac_cells, cb_cells)
  
  if (length(common_cells) == 0) {
    stop(sprintf("Sample '%s': No common barcodes between ATAC (%d cells) and CellBender RNA (%d cells).", 
                 sample, length(atac_cells), length(cb_cells)))
  }
  
  # Report coverage
  prop_atac_covered <- length(common_cells) / length(atac_cells)
  cat(sprintf("Barcode intersection: %d common cells (%.1f%% of ATAC barcodes)\n", 
              length(common_cells), prop_atac_covered * 100))
  
  if (prop_atac_covered < 0.80) {
    warning(sprintf("Sample '%s': Only %.1f%% of ATAC barcodes found in CellBender output (< 80%%). Consider re-running CellBender.", 
                    sample, prop_atac_covered * 100))
  }
  
  # Subset both modalities to common cells
  rna_counts <- cb_rna[, common_cells, drop = FALSE]
  atac_counts <- atac_orig[, common_cells, drop = FALSE]
  
  # Also subset the original CellRanger RNA to common cells (for metrics comparison)
  in_mat <- rna_orig[, common_cells, drop = FALSE]
  cb_mat2 <- cb_rna[, common_cells, drop = FALSE]
  
  return(list(
    rna_counts = rna_counts,
    atac_counts = atac_counts,
    common_cells = common_cells,
    in_mat = in_mat,
    cb_mat2 = cb_mat2
  ))
}

#' Compute CellBender merge metrics
#'
#' @param sample Sample name
#' @param in_mat CellRanger RNA (subsetted to common cells)
#' @param cb_mat2 CellBender RNA (subsetted to common cells)
#' @param atac_orig Original ATAC counts (full)
#' @param cb_rna_full CellBender RNA counts (full)
#' @param common_cells Vector of intersected barcodes
#' @return 1-row data.frame with merge quality metrics
sctb_compute_cb_merge_metrics <- function(sample, in_mat, cb_mat2, atac_orig, cb_rna_full, common_cells) {
  # Barcode overlap metrics
  n_atac <- ncol(atac_orig)
  n_cb <- ncol(cb_rna_full)
  n_intersection <- length(common_cells)
  prop_atac_covered <- n_intersection / n_atac
  
  # Per-cell UMI totals (for correlation)
  cellranger_umis <- Matrix::colSums(in_mat)
  cellbender_umis <- Matrix::colSums(cb_mat2)
  
  # Correlation metrics
  pearson <- cor(cellranger_umis, cellbender_umis, method = "pearson")
  spearman <- cor(cellranger_umis, cellbender_umis, method = "spearman")
  
  # Linear regression (CellBender ~ CellRanger)
  lm_fit <- lm(cellbender_umis ~ cellranger_umis)
  r2 <- summary(lm_fit)$r.squared
  slope <- coef(lm_fit)[2]
  intercept <- coef(lm_fit)[1]
  
  # Removal metrics (weighted by original counts)
  removed_counts <- cellranger_umis - cellbender_umis
  weighted_removed <- sum(removed_counts * cellranger_umis) / sum(cellranger_umis^2)
  
  # Ratio distribution (CB / CellRanger per cell)
  ratios <- cellbender_umis / cellranger_umis
  ratio_median <- median(ratios)
  ratio_q05 <- quantile(ratios, 0.05)
  ratio_q95 <- quantile(ratios, 0.95)
  
  # Non-zero ratio (sparsity comparison)
  nnz_cellranger <- Matrix::nnzero(in_mat)
  nnz_cellbender <- Matrix::nnzero(cb_mat2)
  nnz_ratio <- nnz_cellbender / nnz_cellranger
  
  # Assemble metrics into 1-row data.frame
  metrics <- data.frame(
    sample = sample,
    n_atac = n_atac,
    n_cb = n_cb,
    n_intersection = n_intersection,
    prop_atac_covered = prop_atac_covered,
    pearson = pearson,
    spearman = spearman,
    r2 = r2,
    slope = slope,
    intercept = intercept,
    weighted_removed = weighted_removed,
    ratio_median = ratio_median,
    ratio_q05 = ratio_q05,
    ratio_q95 = ratio_q95,
    nnz_ratio = nnz_ratio,
    stringsAsFactors = FALSE
  )
  
  return(metrics)
}

#' Emit CellBender merge report (display and/or write)
#'
#' @param metrics 1-row data.frame from compute_cb_merge_metrics()
#' @param mode "write", "display", or "none"
#' @param sample Sample name (for filename)
sctb_emit_cb_report <- function(metrics, mode = c("write", "display", "none"), sample) {
  mode <- match.arg(mode)
  
  if (mode == "none") {
    return(invisible(NULL))
  }
  
  # Determine output directory based on mode
  if (mode == "display") {
    # QC mode: write to tmp and display
    out_dir <- file.path(tmpfiledir, "cellbender_merge_reports")
    
    # Ensure directory exists
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Write CSV
    out_path <- file.path(out_dir, paste0(sample, ".csv"))
    write.csv(metrics, out_path, row.names = FALSE)
    
    # Display to console
    cat("\n=== CellBender Merge Report ===\n")
    cat(sprintf("Sample: %s\n", sample))
    cat(sprintf("ATAC cells: %d | CellBender cells: %d | Intersection: %d (%.1f%%)\n",
                metrics$n_atac, metrics$n_cb, metrics$n_intersection, 
                metrics$prop_atac_covered * 100))
    cat(sprintf("Correlation - Pearson: %.3f | Spearman: %.3f | R²: %.3f\n",
                metrics$pearson, metrics$spearman, metrics$r2))
    cat(sprintf("Regression - Slope: %.3f | Intercept: %.1f\n",
                metrics$slope, metrics$intercept))
    cat(sprintf("Removal - Weighted: %.3f | Ratio median: %.3f [Q05: %.3f, Q95: %.3f]\n",
                metrics$weighted_removed, metrics$ratio_median, 
                metrics$ratio_q05, metrics$ratio_q95))
    cat(sprintf("Sparsity - NNZ ratio: %.3f\n", metrics$nnz_ratio))
    cat(sprintf("Report saved to: %s\n", out_path))
    cat("================================\n\n")
    
    ## View in RStudio (safe in interactive sessions)
    #if (interactive()) {
    #  utils::View(metrics, title = paste("CellBender Merge:", sample))
    #}
    
  } else if (mode == "write") {
    # Pipeline mode: write to export directory (no display)
    out_dir <- cellbender_report_dir
    
    # Ensure directory exists
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    out_path <- file.path(out_dir, paste0(sample, ".csv"))
    write.csv(metrics, out_path, row.names = FALSE)
    cat(sprintf("CellBender merge report written to: %s\n", out_path))
  }
  
  return(invisible(NULL))
}

#' Create the initial base seurat object
#'
#' @param samplename Sample identifier
#' @param cb_report CellBender merge report mode: "write", "display", or "none"
#' @return base Seurat5 object with percent.mt added
base_object <- function(samplename, cb_report = c("write", "display", "none")) {
  cb_report <- match.arg(cb_report)
  
  # Always read the multiome H5 (needed for ATAC, and for RNA if not using CellBender)
  fullrna <- paste(cra_outdir, samplename, h5filename, sep = "/")
  fullatac <- paste(cra_outdir, samplename, atacfilename, sep = "/")
  
  counts <- Read10X_h5(filename = fullrna)
  rna_counts_orig <- counts$`Gene Expression`
  atac_counts_orig <- counts$Peaks
  
  # Determine final RNA and ATAC counts based on use_cellbender setting
  if (use_cellbender) {
    cat("CellBender mode enabled: loading corrected RNA counts...\n")
    
    # Read CellBender RNA
    cb_rna <- read_cellbender_rna_counts(samplename)
    
    # Splice CellBender RNA into multiome (intersect barcodes, subset both modalities)
    spliced <- splice_cellbender_rna_into_multiome(
      rna_orig = rna_counts_orig,
      atac_orig = atac_counts_orig,
      cb_rna = cb_rna,
      sample = samplename
    )
    
    # Compute merge metrics
    metrics <- compute_cb_merge_metrics(
      sample = samplename,
      in_mat = spliced$in_mat,
      cb_mat2 = spliced$cb_mat2,
      atac_orig = atac_counts_orig,
      cb_rna_full = cb_rna,
      common_cells = spliced$common_cells
    )
    
    # Emit report (display/write/none)
    emit_cb_report(metrics, mode = cb_report, sample = samplename)
    
    # Use spliced (intersected) counts
    rna_counts <- spliced$rna_counts
    atac_counts <- spliced$atac_counts
    
    cat(sprintf("Using CellBender RNA counts: %d cells, %d genes\n", 
                ncol(rna_counts), nrow(rna_counts)))
    cat(sprintf("ATAC counts subsetted to common cells: %d cells, %d peaks\n", 
                ncol(atac_counts), nrow(atac_counts)))
    
  } else {
    # Use original CellRanger counts
    rna_counts <- rna_counts_orig
    atac_counts <- atac_counts_orig
  }
  
  # Create Seurat object (same logic whether using CellBender or not)
  print("Creating the RNA assay for the Seurat object...")
  baseSeuratObj <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = samplename
  )
  
  print("Adding the ATAC assay for the Seurat object...")
  baseSeuratObj[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fullatac,
    annotation = EnsDbAnnos
  )
  
  print("Calculating a slot for percent.mt for downstream QC")
  DefaultAssay(baseSeuratObj) <- "RNA"
  baseSeuratObj[["percent.mt"]] <- PercentageFeatureSet(baseSeuratObj, pattern = "^MT-")
  
  return(baseSeuratObj)
}

base_qc_object <- function(sample, EnsDbAnnos, save = FALSE, cb_report = "display") {
  cat(sprintf("\nProcessing sample: %s\n", sample))
  base_path <- get_rds_path(sample, "base")
  
  # Create base object with CellBender reporting
  base_obj <- base_object(sample, cb_report = cb_report)
  
  cat("Adding chromosome mapping information to ATAC assay.\n")
  base_obj <- chromosome_mapping(base_obj, rna_annos = EnsDbAnnos)
  
  DefaultAssay(base_obj) <- "ATAC"
  
  cat("Calculating Nucleosome Signal...\n")
  base_obj <- NucleosomeSignal(base_obj)
  
  cat("Calculating TSS Enrichment...\n")
  base_obj <- TSSEnrichment(base_obj)
  
  if (save) {
    saveRDS(base_obj, base_path)
    cat(sprintf("Saved base object to: %s\n", base_path))
  }
  
  return(base_obj)
}

#' Annotate RNA features by EnsDb annotation and ATAC peaks by parsing peak names (hyphen format).
#'
#' @param seurat_obj Seurat object.
#' @param rna_annos GRanges object of gene annotation (from loadannotations/EnsDb).
#' @param warn_threshold Warn if fewer features mapped than this (RNA).
#' @return Seurat object with chromosome info in @misc$feature.info for both RNA and ATAC.
sctb_chromosome_mapping <- function(seurat_obj, rna_annos, warn_threshold = 16000) {
  # RNA: Map gene symbols to chromosomes using annotation
  if (!is.null(seurat_obj[["RNA"]])) {
    feature_names <- rownames(seurat_obj[["RNA"]])
    anno_symbols <- mcols(rna_annos)$gene_name
    #anno_ensembl <- mcols(rna_annos)$gene_id # or $gene_id depending on EnsDb version
    anno_chroms <- as.character(seqnames(rna_annos))
    
    chrom_vec <- anno_chroms[match(feature_names, anno_symbols)]
    names(chrom_vec) <- feature_names
    
    n_mapped <- sum(!is.na(chrom_vec))
    n_total <- length(chrom_vec)
    
    if (n_mapped < warn_threshold) {
      warning(sprintf("chromosome_mapping (RNA): Only %d/%d features mapped to a chromosome (expected 18,000-24,000).", n_mapped, n_total))
    } else {
      message(sprintf("chromosome_mapping (RNA): %d/%d features mapped to a chromosome.", n_mapped, n_total))
    }
    
    # Ensure feature.info is a data.frame with rownames = feature names
    feature_info <- seurat_obj[["RNA"]]@misc$feature.info
    if (is.null(feature_info) || !is.data.frame(feature_info) || nrow(feature_info) == 0) {
      feature_info <- data.frame(row.names = feature_names)
    }
    feature_info <- feature_info[feature_names, , drop=FALSE]
    feature_info$chromosome <- chrom_vec
    seurat_obj[["RNA"]]@misc$feature.info <- feature_info
  }
  
  # ATAC: Parse chromosome from peak names (chrN-...)
  if (!is.null(seurat_obj[["ATAC"]])) {
    peak_names <- rownames(seurat_obj[["ATAC"]])
    chrom_vec <- sub("-.*", "", peak_names)
    names(chrom_vec) <- peak_names
    
    feature_info <- seurat_obj[["ATAC"]]@misc$feature.info
    if (is.null(feature_info) || !is.data.frame(feature_info) || nrow(feature_info) == 0) {
      feature_info <- data.frame(row.names = peak_names)
    }
    feature_info <- feature_info[peak_names, , drop=FALSE]
    feature_info$chromosome <- chrom_vec
    seurat_obj[["ATAC"]]@misc$feature.info <- feature_info
    message(sprintf("chromosome_mapping (ATAC): Chromosome parsed for %d peaks.", length(peak_names)))
  }
  
  return(seurat_obj)
}

merge_sample_objects <- function(samplelist, suffix = "pipeline1", project_name = "opioid", path_fun = get_rds_path) {
  # Get file paths for all samples
  file_paths <- sapply(samplelist, function(sample) path_fun(sample, suffix))
  
  # Check file existence before proceeding
  missing_files <- file_paths[!file.exists(file_paths)]
  if (length(missing_files) > 0) {
    stop(sprintf("Missing RDS files for samples: %s", paste(missing_files, collapse=", ")))
  }
  
  # Read sample objects
  sample_objs <- lapply(file_paths, readRDS)
  
  # Merge sample objects
  merged_seurat <- merge(x = sample_objs[[1]], y = sample_objs[-1])
  
  invisible(merged_seurat)
}
                                                
#TODO: divorce from pipeline1_settings
#' Remove doublets using scDblFinder
#'  
#' Detects and removes doublets from a Seurat object using scDblFinder.
#' Looks up expected doublet rate from pipeline1_settings dataframe.
#'
#' @param seurat_obj A Seurat object (with @project.name set to sample name)
#' @param pipeline1_settings (optional) Data frame with sample-specific settings
#'        If not supplied, will look for `pipeline1_settings` in global environment
#' @param qc_report If TRUE, prints detailed output and plots (for run_qc.R)
#'        If FALSE, minimal output (for run_pipeline1.R)
#' @return List with:
#'   - obj: Filtered Seurat object containing only singlets
#'   - stats: List of doublet detection statistics for reporting
#' @export
doubletRemoveSample <- function(seurat_obj, pipeline1_settings = NULL, qc_report = FALSE) {

  # Get sample name from Seurat object
  sample_name <- seurat_obj@project.name

  # Fetch pipeline1_settings if not provided
  if (is.null(pipeline1_settings)) {
    if (!exists("pipeline1_settings", envir = .GlobalEnv)) {
      stop("pipeline1_settings not found in global environment, and not provided as argument.")
    }
    pipeline1_settings <- get("pipeline1_settings", envir = .GlobalEnv)
  }

  # Look up parameters for this sample
  params <- pipeline1_settings[pipeline1_settings$sample == sample_name, ]
  if (nrow(params) == 0) {
    stop(sprintf("No doublet settings found for sample '%s'.", sample_name))
  }

  # Extract expected doublet rate
  expected_dbr <- params$expected_dbr
  if (is.na(expected_dbr)) {
    stop(sprintf("expected_dbr not set for sample '%s'. Run QC first to calculate this value.",
                sample_name))
  }

  # Get global settings
  if (!exists("random_seed", envir = .GlobalEnv)) {
    stop("random_seed not found. Source project_settings.R first.")
  }
  if (!exists("doublet_rate_sd", envir = .GlobalEnv)) {
    stop("doublet_rate_sd not found. Source project_settings.R first.")
  }

  random_seed <- get("random_seed", envir = .GlobalEnv)
  doublet_rate_sd <- get("doublet_rate_sd", envir = .GlobalEnv)

  # Calculate parameters
  n_cells <- ncol(seurat_obj)
  dbr_rate <- expected_dbr
  est_doubs <- (n_cells * dbr_rate / 100)

  if (qc_report) {
    cat(sprintf("\n=== scDblFinder: %s ===\n", sample_name))
    cat(sprintf("Cells before doublet removal: %d\n", n_cells))
    cat(sprintf("Expected doublets: %.2f (%.2f%% of %.2f cells)\n",
    est_doubs, dbr_rate, n_cells))
    cat(sprintf("Running scDblFinder with dbr=%.3f, dbr.sd=%.3f\n",
        doublet_rate_per_1000, doublet_rate_sd))
  }

  # Convert to SingleCellExperiment (using Seurat's conversion method)
  sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")

  # Set seed for reproducibility
  set.seed(random_seed)

  # Run scDblFinder
  sce <- scDblFinder::scDblFinder(
    sce,
    clusters = FALSE,
    dbr = (dbr_rate / 100),
    dbr.sd = doublet_rate_sd,
    verbose = qc_report
  )

  # Extract results
  doublet_class <- sce$scDblFinder.class
  doublet_score <- sce$scDblFinder.score

  n_doublets <- sum(doublet_class == "doublet")
  n_singlets <- sum(doublet_class == "singlet")
  pct_doublets <- 100 * n_doublets / n_cells

  # Calculate statistics
  threshold <- min(doublet_score[doublet_class == "doublet"])
  singlet_score_median <- median(doublet_score[doublet_class == "singlet"])
  doublet_score_median <- median(doublet_score[doublet_class == "doublet"])

  # Report results
  if (qc_report) {
    cat(sprintf("\n=== Results ===\n"))
    cat(sprintf("Singlets: %d (%.1f%%)\n", n_singlets, 100*n_singlets/n_cells))
    cat(sprintf("Doublets: %d (%.1f%%)\n", n_doublets, pct_doublets))
    cat(sprintf("Expected: %.1f (%.2f%%)\n", est_doubs, dbr_rate))
    cat(sprintf("Difference: %+.1f doublets (%+.2f%%)\n\n",
               n_doublets - est_doubs, pct_doublets - dbr_rate))

    cat("Doublet score summary:\n")
    print(summary(doublet_score))
    cat("\nSinglet scores:\n")
    print(summary(doublet_score[doublet_class == "singlet"]))
    cat("Doublet scores:\n")
    print(summary(doublet_score[doublet_class == "doublet"]))

    # Two-panel visualization
    par(mfrow = c(1, 2))

    # Panel 1: Full distribution
    hist(doublet_score, breaks = 50,
         main = paste("All Scores -", sample_name),
         xlab = "scDblFinder Score",
         ylab = "Frequency",
         col = "lightblue")
    abline(v = threshold, col = "red", lwd = 2, lty = 2)
    legend("topright",
           legend = sprintf("Threshold: %.3f", threshold),
           col = "red", lty = 2, lwd = 2, cex = 0.8)

    # Panel 2: Doublet region zoomed (0.7-1.0)
    doublet_region <- doublet_score[doublet_score > 0.7]
    if (length(doublet_region) > 0) {
      hist(doublet_region, breaks = 30,
           main = "Doublet Region (0.7-1.0)",
           xlab = "scDblFinder Score",
           ylab = "Frequency",
           col = "salmon",
           xlim = c(0.7, 1.0))
      abline(v = threshold, col = "red", lwd = 2, lty = 2)
      abline(v = doublet_score_median, col = "darkred", lwd = 2, lty = 3)
      legend("topleft",
             legend = c(sprintf("Threshold: %.3f", threshold),
                       sprintf("Median: %.3f", doublet_score_median)),
             col = c("red", "darkred"),
             lty = c(2, 3),
             lwd = 2,
             cex = 0.7)
      text(0.85, max(hist(doublet_region, breaks = 30, plot = FALSE)$counts) * 0.9,
           sprintf("%d doublets\n(%.1f%%)", n_doublets, pct_doublets),
           cex = 0.9)
    }

    par(mfrow = c(1, 1))

    cat("\nFiltering to singlets only...\n")
  }

  # Filter Seurat object to singlets only
  singlets_idx <- which(doublet_class == "singlet")
  filtered_obj <- seurat_obj[, singlets_idx]

  if (qc_report) {
    cat(sprintf("Cells after doublet removal: %d\n", ncol(filtered_obj)))
  } else {
    cat(sprintf("Removed %d doublets from %s, %d cells remaining\n",
               n_doublets, sample_name, ncol(filtered_obj)))
  }

  # Compile statistics for reporting
  doublet_stats <- list(
    n_cells_before = n_cells,
    expected_dbr = expected_dbr,
    n_doublets = n_doublets,
    n_singlets = n_singlets,
    pct_doublets = pct_doublets,
    threshold = threshold,
    singlet_score_median = singlet_score_median,
    doublet_score_median = doublet_score_median
  )

  # Return both the filtered object and statistics
  return(list(
    obj = filtered_obj,
    stats = doublet_stats
  ))
}
# Contains rna and ATAC modality functions, and harmony/batch effects reductions

#' @param seurat_obj Seurat object.
#' @return Seurat object with FindVariableFeatures, ScaleData, and RunPCA
stb_RNA_FSR <- function(rna_obj) {
  DefaultAssay(rna_obj) <- "RNA"
  rna_obj <- FindVariableFeatures(rna_obj, assay = "RNA")
  rna_obj <- ScaleData(rna_obj, assay = "RNA",
                        features = VariableFeatures(rna_obj))
  rna_obj <- RunPCA(pm_rna_obj, assay = "RNA",
                      features = VariableFeatures(rna_obj))
  return(rna_obj)
}

#' @param seurat_obj Seurat object.
#' @param string with quality score value for min cutoff
#' @return Chromatin Assay with RunTFIDF, FindTopFeatures, and RunSVD
stb_ATAC_RFR <- function(atac_obj, qvalue = "q0") {
  DefaultAssay(atac_obj) <- "ATAC"
  pm_atac_obj <- RunTFIDF(atac_obj)
  pm_atac_obj <- FindTopFeatures(atac_obj, min.cutoff = qvalue)
  pm_atac_obj <- RunSVD(atac_obj)
  return(atac_obj)
}

stb_harmonize_both <- function(harmony_obj, harmony_max_iter = 50,
                         harmony_project.dim = FALSE,
                         harmony_dims = NULL, random_seed = NULL,
			 plot_convergence = FALSE) { 
  DefaultAssay(harmony_obj) <- "RNA"
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }  
  harmony_obj <- RunHarmony(
    harmony_obj,
    group.by.vars = "orig.ident",
    reduction.use = "pca",
    plot_convergence = plot_convergence,
    max_iter = harmony_max_iter,
    reduction.save = reduction.save.RNA,
    project.dim = harmony_project.dim,
    dims.use = harmony_dims
  )

  DefaultAssay(harmony_obj) <- "ATAC"
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }  
  harmony_obj <- RunHarmony(
    object = harmony_obj,
    group.by.vars = "orig.ident",
    reduction.use = "lsi",
    project.dim = harmony_project.dim,
    max_iter = harmony_max_iter,
    reduction.save = reduction.save.ATAC,
    dims.use = harmony_dims
  )
  return(harmony_obj)
}

stb_FMMN_task <- function(FMMN_obj, knn, dims) {
  FMMN_obj <- FindMultiModalNeighbors(
    object = FMMN_obj,
    reduction.list = list(reduction.save.RNA, reduction.save.ATAC),
    dims.list = list(dims, dims),  # Use all harmony dims
    k.nn = knn,
    knn.graph.name = "wknn",
    snn.graph.name = "wsnn",
    weighted.nn.name = "weighted.nn"
  )
  return(FMMN_obj)
}

stb_cluster_data <- function(harmony_obj, alg, res, run_umap = FALSE, cluster_seed,
                         singleton_handling = c("discard", "merge", "keep")) {
  singleton_handling <- match.arg(singleton_handling)

  # The DefaultAssay is being set for consistent behavior, not because we're doing
  # assay-specific actions; I don't want a minor unintended change to occur just
  # because the assay was ATAC or something else, instead of RNA.
  DefaultAssay(harmony_obj) <- "RNA"

  # Determine group.singletons argument for FindClusters
  group_singletons <- (singleton_handling == "merge")

  # Clustering step
  harmony_obj <- FindClusters(
    harmony_obj,
    graph.name = "wsnn",
    algorithm = alg,
    resolution = res,
    group.singletons = group_singletons,
    # I do not like needing the below, and will work on the data until it isn't
    # needed.  Stable data doesn't change with new random seeds.  Hardcoding the
    # seed is cheating.  All but one cluster/cell type are very very stable, so
    # setting this allows parts of the project to move forward without it
    random.seed = cluster_seed
  )

  # Run UMAP BEFORE discarding singletons
  # This ensures weighted.nn neighbor graph is still present
  if (run_umap) {
    harmony_obj <- RunUMAP(
      harmony_obj,
      nn.name = "weighted.nn",            # This matches weighted.nn.name from FMMN_task
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_"
    )
  }

  # NOW discard singletons (after UMAP is calculated)
  # The UMAP coordinates for retained cells will be preserved during subsetting
  if (singleton_handling == "discard") {
    if ("singleton" %in% levels(harmony_obj$seurat_clusters)) {
      singleton_cells <- WhichCells(harmony_obj, idents = "singleton")
      harmony_obj <- subset(harmony_obj, cells = setdiff(colnames(harmony_obj), singleton_cells))
      # Drop unused cluster level
      harmony_obj$seurat_clusters <- droplevels(harmony_obj$seurat_clusters)

      cat(sprintf("Discarded %d singleton cells\n", length(singleton_cells)))
    }
  }

  return(harmony_obj)
}
