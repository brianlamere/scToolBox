# ============================================================================
# Consolidated_Markers.R
# ============================================================================
# Tissue-specific marker sets for cell type identification.
# Each function returns a list with two elements:
#   $reference_table : data.frame with columns gene, source, confidence,
#                      notes, celltype -- one row per marker gene
#   $marker_lists    : named list of character vectors, one entry per cell
#                      type, used directly by identify_all_celltypes()
#
# To add a new tissue:
#   1. Write a new function following the pattern below
#   2. Add an entry to the TISSUE_MARKER_FUNCTIONS list in Rmultiome-main.R
#
# ============================================================================

# ============================================================================
# prefrontal cortex / cortical nuclei
# ============================================================================
#
# Primary Sources:
# 1. Mathys et al. 2019 Nature (DOI: 10.1038/s41586-019-1195-2)
#    - Prefrontal cortex snRNA-seq, Alzheimer's disease
# 2. Sankowski et al. 2019 Nat Neurosci (DOI: 10.1038/s41593-019-0532-y)
#    - Microglia vs macrophage distinction in human brain
# 3. Zillich et al. 2025 Nat Commun (DOI: 10.1038/s41467-024-55467-0)
#    - Caudate nucleus multiome, validation reference
# ============================================================================

cortex_consolidated_markers <- function() {

  oligodendrocyte_markers <- data.frame(
    gene = c("MBP", "MOBP", "PLP1", "MOG", "MAG"),
    source = "Mathys_2019",
    confidence = c("high", "high", "high", "medium", "medium"),
    notes = "Mature myelinating oligodendrocytes",
    stringsAsFactors = FALSE
  )

  opc_markers <- data.frame(
    gene = c("VCAN", "PDGFRA", "PCDH15", "CSPG4", "SOX10"),
    source = "Mathys_2019",
    confidence = c("high", "high", "high", "medium", "medium"),
    notes = "Oligodendrocyte precursors, proliferative",
    stringsAsFactors = FALSE
  )

  astrocyte_markers <- data.frame(
    gene = c("GFAP", "ALDH1L1", "GLUL", "AQP4", "SLC1A2", "SLC4A4"),
    source = c(rep("Mathys_2019", 4), "Zillich_2025", "Zillich_2025"),
    confidence = rep("high", 6),
    notes = "Mature protoplasmic astrocytes; GFAP may be low in homeostatic state",
    stringsAsFactors = FALSE
  )

  microglia_markers <- data.frame(
    gene = c("P2RY12", "TMEM119", "CX3CR1", "GPR34", "OLFML3",
             "APBB1IP", "CSF1R", "SALL1"),
    source = c(rep("Sankowski_2019", 5), "Zillich_2025",
               "Mathys_2019", "Sankowski_2019"),
    confidence = c("high", "high", "high", "high", "high",
                  "medium", "low_specificity", "high"),
    notes = c(
      "THE definitive microglia marker (homeostatic)",
      "Microglia-specific transmembrane protein",
      "Fractalkine receptor, high in microglia",
      "GPCR, microglia-specific",
      "Olfactomedin-like, microglia-specific",
      "Microglia marker from Zillich",
      "Pan-myeloid, NOT specific (both microglia and macrophages)",
      "Transcription factor, microglia identity"
    ),
    stringsAsFactors = FALSE
  )

  macrophage_markers <- data.frame(
    gene = c("CD163", "MRC1", "LYVE1", "CD14", "FCGR3A", "MS4A7", "F13A1"),
    source = rep("Sankowski_2019", 7),
    confidence = c("high", "high", "high", "high", "medium", "medium", "low"),
    notes = c(
      "Scavenger receptor, perivascular macrophages (CRITICAL for distinguishing from microglia)",
      "Mannose receptor (CD206), macrophages",
      "Lymphatic marker, border-associated macrophages",
      "Monocyte marker, blood-derived macrophages",
      "CD16, non-classical monocytes",
      "Monocyte/macrophage marker",
      "Factor XIII, tissue macrophages"
    ),
    stringsAsFactors = FALSE
  )

  excitatory_markers <- data.frame(
    gene = c("SLC17A7", "CAMK2A", "SATB2", "TBR1", "NRGN"),
    source = c(rep("Mathys_2019", 4), "General"),
    confidence = rep("high", 5),
    notes = "Glutamatergic neurons; SATB2/TBR1 for layer identity",
    stringsAsFactors = FALSE
  )

  gabaergic_core <- data.frame(
    gene = c("GAD1", "GAD2", "SLC32A1"),
    source = rep("Mathys_2019", 3),
    confidence = rep("high", 3),
    notes = "Pan-GABAergic markers (all interneurons)",
    stringsAsFactors = FALSE
  )

  gabaergic_subtypes <- data.frame(
    gene = c("PVALB", "SST", "VIP", "LAMP5", "CALB2"),
    source = rep("Mathys_2019", 5),
    confidence = rep("high", 5),
    notes = c(
      "Parvalbumin+ fast-spiking interneurons",
      "Somatostatin+ Martinotti cells",
      "Vasoactive intestinal peptide+ interneurons",
      "LAMP5+ neurogliaform cells",
      "CALB2 (calretinin) in VIP/LAMP5 subsets"
    ),
    stringsAsFactors = FALSE
  )

  endothelial_markers <- data.frame(
    gene = c("FLT1", "CLDN5", "KDR"),
    source = rep("Mathys_2019", 4),
    confidence = c("high", "high", "high", "medium"),
    notes = "Vascular endothelial cells; CLDN5 = tight junctions (BBB)",
    stringsAsFactors = FALSE
  )

  pericyte_markers <- data.frame(
    gene = c("PDGFRB", "RGS5", "NOTCH3", "CARMN", "ACTA2"),
    source = rep("Mathys_2019", 5),
    confidence = c("high", "high", "high", "medium", "low"),
    notes = c(
      "PDGF receptor beta, definitive pericyte marker",
      "Regulator of G-protein signaling 5",
      "NOTCH3 mutations cause CADASIL (pericyte dysfunction)",
      "Cardiac mesoderm marker",
      "Smooth muscle actin (shared with smooth muscle)"
    ),
    stringsAsFactors = FALSE
  )

  lymphocyte_markers <- data.frame(
    gene = c("CD96", "CD3D", "IL2RG"),
    source = rep("Zillich_2025", 3),
    confidence = c("medium", "high", "high"),
    notes = c(
      "NK cell and T-cell marker",
      "T-cell receptor component (definitive)",
      "Common gamma chain, lymphocyte marker"
    ),
    stringsAsFactors = FALSE
  )

  stress_qc_markers <- data.frame(
    gene = c("VIM", "GFAP", "NEFM", "NEFL", "NRGN"),
    source = "Your_pipeline",
    confidence = "QC_flag",
    notes = "Combination suggests stressed/dying cells with mixed identity (remove from analysis)",
    stringsAsFactors = FALSE
  )

  reference_table <- dplyr::bind_rows(
    oligodendrocyte_markers %>% dplyr::mutate(celltype = "Oligodendrocytes"),
    opc_markers             %>% dplyr::mutate(celltype = "OPCs"),
    astrocyte_markers       %>% dplyr::mutate(celltype = "Astrocytes"),
    microglia_markers       %>% dplyr::mutate(celltype = "Microglia"),
    macrophage_markers      %>% dplyr::mutate(celltype = "Macrophages"),
    excitatory_markers      %>% dplyr::mutate(celltype = "Excitatory_Neurons"),
    gabaergic_core          %>% dplyr::mutate(celltype = "GABAergic_Core"),
    gabaergic_subtypes      %>% dplyr::mutate(
                                  celltype = paste0("GABAergic_", gene)),
    endothelial_markers     %>% dplyr::mutate(celltype = "Endothelial"),
    pericyte_markers        %>% dplyr::mutate(celltype = "Pericytes"),
    lymphocyte_markers      %>% dplyr::mutate(celltype = "Lymphocytes"),
    stress_qc_markers       %>% dplyr::mutate(celltype = "QC_Flag")
  )

  marker_lists <- list(
    Oligodendrocytes    = c("MBP", "MOBP", "PLP1"),
    OPCs                = c("VCAN", "PDGFRA", "PCDH15"),
    Astrocytes          = c("GFAP", "ALDH1L1", "AQP4", "SLC1A2"),
    Microglia           = c("P2RY12", "TMEM119", "CX3CR1"),
    Macrophages         = c("CD163", "MRC1", "LYVE1"),
    Excitatory_Neurons  = c("SLC17A7", "CAMK2A", "SATB2", "TBR1"),
    GABAergic_PVALB     = c("GAD1", "GAD2", "PVALB"),
    GABAergic_SST       = c("GAD1", "GAD2", "SST"),
    GABAergic_VIP       = c("GAD1", "GAD2", "VIP"),
    GABAergic_LAMP5     = c("GAD1", "GAD2", "LAMP5"),
    Endothelial         = c("FLT1", "CLDN5"),
    Pericytes           = c("PDGFRB", "RGS5", "NOTCH3")
  )

  list(reference_table = reference_table, marker_lists = marker_lists)
}


# ============================================================================
# SPLEEN
# ============================================================================
#
# Primary Sources:
# 1. Madissoon et al. 2020 Genome Biology (DOI: 10.1186/s13059-019-1906-x)
#    - Human spleen scRNA-seq, healthy donors
# 2. Brown et al. 2019 Science Immunology (DOI: 10.1126/sciimmunol.aax4385)
#    - Human spleen macrophage subsets including red pulp macrophages
# 3. Sankowski et al. 2019 Nat Neurosci (DOI: 10.1038/s41593-019-0532-y)
#    - Monocyte/macrophage distinction (cross-tissue applicable)
# 4. General immunology consensus markers (textbook-level confidence)
#
# NOTE ON HIV/OPIOID COHORT:
# These markers identify cell IDENTITY, not activation state. HIV and opioid
# exposure alter proportions and activation states but do not abolish lineage
# markers. Exhaustion/activation markers are noted separately.
# ============================================================================

spleen_consolidated_markers <- function() {

  t_cell_markers <- data.frame(
    gene = c("CD3D", "CD3E", "CD3G", "TRAC"),
    source = rep("General_consensus", 4),
    confidence = rep("high", 4),
    notes = "Pan-T cell; TRAC = TCR alpha constant region",
    stringsAsFactors = FALSE
  )

  cd4_t_markers <- data.frame(
    gene = c("CD4", "IL7R", "CCR7", "TCF7"),
    source = c("General_consensus", rep("Madissoon_2020", 3)),
    confidence = c("high", "high", "medium", "medium"),
    notes = c("CD4 coreceptor", "IL-7 receptor, naive/memory",
              "CCR7, lymph node homing", "TCF1, naive T cell TF"),
    stringsAsFactors = FALSE
  )

  cd8_t_markers <- data.frame(
    gene = c("CD8A", "CD8B", "GZMK", "GZMB", "PRF1", "NKG7"),
    source = c(rep("General_consensus", 2), rep("Madissoon_2020", 4)),
    confidence = c("high", "high", "medium", "medium", "medium", "medium"),
    notes = c("CD8 alpha chain", "CD8 beta chain",
              "Granzyme K, effector memory CD8",
              "Granzyme B, cytotoxic; elevated in HIV",
              "Perforin, cytotoxic machinery",
              "NK/cytotoxic granule protein 7"),
    stringsAsFactors = FALSE
  )

  treg_markers <- data.frame(
    gene = c("FOXP3", "IL2RA", "IKZF2", "CTLA4"),
    source = c(rep("General_consensus", 2), rep("Madissoon_2020", 2)),
    confidence = c("high", "high", "medium", "medium"),
    notes = c("FOXP3 - definitive Treg TF", "CD25, activated Tregs",
              "Helios, thymic Tregs",
              "CTLA-4, immune checkpoint; elevated in HIV"),
    stringsAsFactors = FALSE
  )

  gammadelta_t_markers <- data.frame(
    gene = c("TRGC1", "TRGC2", "TRDC", "ZBTB16"),
    source = rep("Madissoon_2020", 4),
    confidence = c("high", "high", "high", "medium"),
    notes = c("TCR gamma constant 1", "TCR gamma constant 2",
              "TCR delta constant - definitive gamma-delta T marker",
              "PLZF TF, innate-like gamma-delta T"),
    stringsAsFactors = FALSE
  )

  b_cell_markers <- data.frame(
    gene = c("MS4A1", "CD19", "PAX5", "CD79A", "CD79B"),
    source = c(rep("General_consensus", 3), rep("Madissoon_2020", 2)),
    confidence = rep("high", 5),
    notes = c("CD20, pan-B", "CD19, pan-B co-receptor",
              "PAX5, B lineage TF", "CD79a, BCR signaling",
              "CD79b, BCR signaling"),
    stringsAsFactors = FALSE
  )

  gc_b_markers <- data.frame(
    gene = c("BCL6", "AICDA", "CXCR5", "MEF2B"),
    source = rep("Madissoon_2020", 4),
    confidence = c("high", "high", "medium", "medium"),
    notes = c("BCL6 - master GC TF",
              "AID (AICDA) - somatic hypermutation; definitive GC",
              "CXCR5 - follicular homing",
              "MEF2B - GC B cell TF"),
    stringsAsFactors = FALSE
  )

  mz_b_markers <- data.frame(
    gene = c("FCRL4", "FCRL5"),
    source = rep("Madissoon_2020", 2),
    confidence = rep("high", 2),
    notes = c("FcRL4 - marginal zone B (spleen-specific)",
              "FcRL5 - marginal zone/memory B"),
    stringsAsFactors = FALSE
  )

  plasma_markers <- data.frame(
    gene = c("MZB1", "SDC1", "JCHAIN", "PRDM1", "XBP1"),
    source = rep("Madissoon_2020", 5),
    confidence = c("high", "high", "high", "high", "medium"),
    notes = c("MZB1 - plasma cell marker",
              "CD138 (SDC1) - definitive plasma cell surface marker",
              "J chain - Ig joining chain",
              "BLIMP1 (PRDM1) - plasma cell master TF",
              "XBP1 - ER stress, plasma cells"),
    stringsAsFactors = FALSE
  )

  nk_markers <- data.frame(
    gene = c("NCAM1", "KLRB1", "KLRD1", "GNLY", "NKG7", "FCGR3A"),
    source = c("General_consensus", rep("Madissoon_2020", 5)),
    confidence = c("high", "high", "high", "high", "high", "medium"),
    notes = c("CD56 (NCAM1) - NK marker; also on NKT",
              "CD161 (KLRB1) - NK/NKT",
              "CD94 (KLRD1) - NK, forms heterodimer with NKG2",
              "Granulysin - NK cytotoxicity; reduced in HIV",
              "NKG7 - NK/cytotoxic granule protein",
              "CD16 (FCGR3A) - mature NK cells"),
    stringsAsFactors = FALSE
  )

  classical_mono_markers <- data.frame(
    gene = c("CD14", "LYZ", "CST3", "FCN1", "S100A8", "S100A9", "VCAN"),
    source = c("General_consensus", rep("Madissoon_2020", 6)),
    confidence = rep("high", 7),
    notes = c("CD14 - classical monocyte (CD14++CD16-)",
              "LYZ - lysozyme, myeloid",
              "CST3 - cystatin C, monocyte",
              "FCN1 - ficolin 1, classical monocyte",
              "S100A8 - calprotectin; elevated in HIV",
              "S100A9 - calprotectin; elevated in HIV",
              "VCAN - versican, classical monocyte"),
    stringsAsFactors = FALSE
  )

  nonclassical_mono_markers <- data.frame(
    gene = c("FCGR3A", "CX3CR1", "CDKN1C", "LST1"),
    source = rep("Madissoon_2020", 4),
    confidence = c("high", "high", "high", "medium"),
    notes = c("CD16 (FCGR3A) - non-classical (CD14+CD16++)",
              "CX3CR1 - fractalkine receptor, patrolling",
              "CDKN1C - p57, non-classical monocyte",
              "LST1 - non-classical monocyte"),
    stringsAsFactors = FALSE
  )

  rpm_markers <- data.frame(
    gene = c("VSIG4", "SLC40A1", "HMOX1", "TIMD4", "FOLR2",
             "CD163", "MRC1", "F13A1"),
    source = c(rep("Brown_2019", 5), rep("Sankowski_2019", 3)),
    confidence = c("high", "high", "high", "high", "high",
                  "high", "medium", "medium"),
    notes = c(
      "VSIG4 - definitive red pulp macrophage",
      "SLC40A1 - ferroportin, iron export; red pulp specialization",
      "HMOX1 - heme oxygenase; heme catabolism from RBC clearance",
      "TIMD4 - TIM4; tissue-resident macrophage retention",
      "FOLR2 - folate receptor beta; tissue macrophages",
      "CD163 - scavenger receptor; red pulp/perivascular macro",
      "MRC1 - CD206, mannose receptor; tissue macrophages",
      "F13A1 - factor XIII; tissue macrophages"
    ),
    stringsAsFactors = FALSE
  )

  cdc1_markers <- data.frame(
    gene = c("CLEC9A", "XCR1", "CADM1"),
    source = rep("Madissoon_2020", 3),
    confidence = rep("high", 3),
    notes = c("CLEC9A (DNGR-1) - cDC1 specific",
              "XCR1 - definitive cDC1 marker",
              "CADM1 - cell adhesion, cDC1"),
    stringsAsFactors = FALSE
  )

  cdc2_markers <- data.frame(
    gene = c("CLEC10A", "CD1C", "FCER1A"),
    source = rep("Madissoon_2020", 3),
    confidence = rep("high", 3),
    notes = c("CLEC10A (CD301) - cDC2 lectin",
              "CD1c - cDC2 antigen presentation",
              "FCER1A - Fc epsilon receptor; cDC2"),
    stringsAsFactors = FALSE
  )

  pdc_markers <- data.frame(
    gene = c("LILRA4", "CLEC4C", "IL3RA"),
    source = rep("Madissoon_2020", 3),
    confidence = rep("high", 3),
    notes = c("LILRA4 - ILT7; definitive pDC",
              "CLEC4C - BDCA-2; pDC surface",
              "IL3RA - CD123, pDC"),
    stringsAsFactors = FALSE
  )

  endothelial_markers <- data.frame(
    gene = c("PECAM1", "FLT1", "CDH5", "CLDN5"),
    source = c(rep("General_consensus", 2), rep("Madissoon_2020", 2)),
    confidence = rep("high", 4),
    notes = c("CD31 (PECAM1) - pan-endothelial",
              "VEGFR1 (FLT1) - vascular endothelial",
              "VE-cadherin (CDH5) - adherens junctions",
              "Claudin-5 - tight junctions"),
    stringsAsFactors = FALSE
  )

  stromal_markers <- data.frame(
    gene = c("PDPN", "CCL19", "CCL21", "CXCL13"),
    source = rep("Madissoon_2020", 4),
    confidence = c("high", "high", "high", "high"),
    notes = c("Podoplanin (PDPN) - fibroblastic reticular cells (FRC)",
              "CCL19 - T cell zone FRC chemokine",
              "CCL21 - T cell zone FRC chemokine",
              "CXCL13 - follicular DC / B cell zone organizer"),
    stringsAsFactors = FALSE
  )

  # --- Contamination markers ---
  # Present in marker_lists so identify_all_celltypes() can score them.
  # Clusters scoring highest on these get typed as Contamination_* and
  # removed via apply_celltype_labels(..., remove_flagged = TRUE).

  rbc_markers <- data.frame(
    gene = c("HBA1", "HBA2", "HBB", "ALAS2"),
    source = rep("General_consensus", 4),
    confidence = rep("high", 4),
    notes = c("Hemoglobin alpha 1 - definitive RBC",
              "Hemoglobin alpha 2 - RBC",
              "Hemoglobin beta - RBC",
              "ALAS2 - erythroid heme synthesis"),
    stringsAsFactors = FALSE
  )

  platelet_markers <- data.frame(
    gene = c("PPBP", "PF4", "GP9", "ITGA2B","SNAP23",
        "HSPD1","HSPH1","HSP90AB1"),
    source = rep("General_consensus", 8),
    confidence = rep("high", 8),
    notes = c("CXCL7 (PPBP) - definitive platelet",
              "CXCL4 (PF4) - definitive platelet",
              "GP9 - platelet surface glycoprotein",
              "CD41 (ITGA2B) - platelet/megakaryocyte",
              "SNAP23", "HSPD1","HSPH1","HSP90AB1"),
    stringsAsFactors = FALSE
  )

  neutrophil_markers <- data.frame(
    gene = c("FCGR3B", "CXCR2", "S100A12", "CSF3R"),
    source = rep("General_consensus", 4),
    confidence = rep("high", 4),
    notes = c("CD16b (FCGR3B) - neutrophil-specific CD16 isoform",
              "IL-8 receptor (CXCR2) - neutrophil chemotaxis",
              "S100A12 - neutrophil-specific calprotectin",
              "G-CSF receptor (CSF3R) - granulocyte"),
    stringsAsFactors = FALSE
  )

  reference_table <- dplyr::bind_rows(
    t_cell_markers          %>% dplyr::mutate(celltype = "T_cell_Pan"),
    cd4_t_markers           %>% dplyr::mutate(celltype = "T_CD4"),
    cd8_t_markers           %>% dplyr::mutate(celltype = "T_CD8"),
    treg_markers            %>% dplyr::mutate(celltype = "T_Regulatory"),
    gammadelta_t_markers    %>% dplyr::mutate(celltype = "T_GammaDelta"),
    b_cell_markers          %>% dplyr::mutate(celltype = "B_cell_Pan"),
    gc_b_markers            %>% dplyr::mutate(celltype = "B_GerminalCenter"),
    mz_b_markers            %>% dplyr::mutate(celltype = "B_MarginalZone"),
    plasma_markers          %>% dplyr::mutate(celltype = "B_Plasma"),
    nk_markers              %>% dplyr::mutate(celltype = "NK_cell"),
    classical_mono_markers  %>% dplyr::mutate(celltype = "Monocyte_Classical"),
    nonclassical_mono_markers %>% dplyr::mutate(celltype = "Monocyte_NonClassical"),
    rpm_markers             %>% dplyr::mutate(celltype = "Macrophage_RedPulp"),
    cdc1_markers            %>% dplyr::mutate(celltype = "DC_cDC1"),
    cdc2_markers            %>% dplyr::mutate(celltype = "DC_cDC2"),
    pdc_markers             %>% dplyr::mutate(celltype = "DC_Plasmacytoid"),
    endothelial_markers     %>% dplyr::mutate(celltype = "Endothelial"),
    stromal_markers         %>% dplyr::mutate(celltype = "Stromal_FRC"),
    rbc_markers             %>% dplyr::mutate(celltype = "Contamination_RBC"),
    platelet_markers        %>% dplyr::mutate(celltype = "Contamination_Platelet"),
    neutrophil_markers      %>% dplyr::mutate(celltype = "Contamination_Neutrophil")
  )

  marker_lists <- list(
    T_cell_Pan              = c("CD3D", "CD3E", "TRAC"),
    T_CD4                   = c("CD4", "IL7R", "CCR7"),
    T_CD8                   = c("CD8A", "CD8B", "GZMK"),
    T_Regulatory            = c("FOXP3", "IL2RA", "IKZF2"),
    T_GammaDelta            = c("TRGC2", "TRDC"),
    B_cell_Pan              = c("MS4A1", "CD19", "PAX5", "CD79A"),
    B_GerminalCenter        = c("BCL6", "AICDA", "CXCR5"),
    B_MarginalZone          = c("FCRL4", "FCRL5"),
    B_Plasma                = c("MZB1", "SDC1", "JCHAIN", "PRDM1"),
    NK_cell                 = c("NCAM1", "KLRB1", "GNLY", "NKG7", "KLRD1"),
    Monocyte_Classical      = c("CD14", "LYZ", "FCN1", "S100A8", "S100A9"),
    Monocyte_NonClassical   = c("FCGR3A", "CX3CR1", "CDKN1C"),
    Macrophage_RedPulp      = c("VSIG4", "SLC40A1", "HMOX1", "TIMD4", "FOLR2"),
    DC_cDC1                 = c("CLEC9A", "XCR1", "CADM1"),
    DC_cDC2                 = c("CLEC10A", "CD1C", "FCER1A"),
    DC_Plasmacytoid         = c("LILRA4", "CLEC4C", "IL3RA"),
    Endothelial             = c("PECAM1", "CDH5", "FLT1"),
    Stromal_FRC             = c("PDPN", "CCL19", "CCL21", "CXCL13"),
    Contamination_RBC       = c("HBA1", "HBA2", "HBB", "ALAS2"),
    Contamination_Platelet  = c("PPBP", "PF4", "GP9", "ITGA2B","SNAP23","HSPD1"),
    Contamination_Neutrophil = c("FCGR3B", "CXCR2", "S100A12", "CSF3R")
  )

  list(reference_table = reference_table, marker_lists = marker_lists)
}
