require(SingleCellExperiment)
require(cellassign)

AD_SCE <- readRDS("data/Processed_10x_AD_data/AD_SCE_normalized.rds")
brainMarkers <- readRDS("data/markers/Brain_markers.RDS")$Celltype

mgi <- marker_list_to_mat(brainMarkers, include_other = FALSE)

marker_in_sce <- match(rownames(mgi), rowData(AD_SCE)$hgnc_symbol)
marker_in_sce <- marker_in_sce[!is.na(marker_in_sce)]

marker_SCE <- AD_SCE[marker_in_sce[!is.na(marker_in_sce)]]
counts(marker_SCE) <- as.matrix(counts(marker_SCE))

mgi <- mgi[!rownames(mgi) %in% setdiff(rownames(mgi), rownames(marker_SCE)), ]

fit <- cellassign(
  exprs_obj = marker_SCE,
  marker_gene_info = mgi,
  s = sizeFactors(AD_SCE),
  shrinkage = TRUE,
  max_iter_adam = 50,
  min_delta = 2,
  n_batches = 10,
  verbose = TRUE
)

saveRDS(fit, "data/Processed_10x_AD_data/cellAssignFit.rds")
AD_SCE$cellassign.type <- fit$cell_type
saveRDS(AD_SCE, "data/Processed_10x_AD_data/AD_CellAssign_Types.rds")
