# restrict an atac assay to peaks in a bedfile
sctb_atac_peak_subset <- function(dmr_peaks_bed = dmr_peaks_bed, seurat_obj = seurat_obj) {
    cat("Reading DMR-overlapping peak set...\n")
    dmr_peaks <- read.table(dmr_peaks_bed, header = FALSE, sep = "\t",
                        col.names = c("chr", "start", "end", "name", "score"),
                        stringsAsFactors = FALSE)
    cat(sprintf("  DMR peak set: %d peaks\n", nrow(dmr_peaks)))

    dmr_gr <- GRanges(
     seqnames = dmr_peaks$chr,
     ranges   = IRanges(start = dmr_peaks$start + 1,   # convert to 1-based
                     end   = dmr_peaks$end)
    )

    # Build GRanges for this object's ATAC peaks
    obj_peaks    <- granges(seurat_obj[["ATAC"]])
    obj_peak_ids <- rownames(seurat_obj[["ATAC"]])

    # findOverlaps: any overlap between object peaks and DMR peaks
    overlaps     <- findOverlaps(obj_peaks, dmr_gr)
    matched_idx  <- unique(queryHits(overlaps))
    matched_peaks <- obj_peak_ids[matched_idx]

    cat(sprintf("  Object ATAC peaks: %d\n", length(obj_peak_ids)))
    cat(sprintf("  Object peaks overlapping DMR set: %d\n", length(matched_peaks)))

    if (length(matched_peaks) < 500) {
        stop(sprintf(
        "Only %d peaks matched — check coordinate systems or peak file format.",
        length(matched_peaks)
    ))
    }

    cat("Subsetting ATAC to DMR-overlapping peaks...\n")
    seurat_sub <- seurat_obj
    seurat_sub[["ATAC"]] <- subset(seurat_obj[["ATAC"]], features = matched_peaks)
    cat(sprintf("  ATAC assay now: %d peaks x %d cells\n",
            nrow(seurat_sub[["ATAC"]]), ncol(seurat_sub)))
    return(seurat_sub)
}

#TODO: description
#can switch from v86 to v116
sctb_loadannotations <- function(ensdb = EnsDb.Hsapiens.v116) {
  annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  return(annotation)
}
