require(SingleCellExperiment)
require(cellassign)

AD_SCE <- readRDS("data/Processed_10x_AD_data/AD_SCE_normalized.rds")
brainMarkers <- readRDS("data/markers/Brain_markers.RDS")$Celltype

mgi <- marker_list_to_mat(brainMarkers, include_other = FALSE)
marker_in_sce <- match(rownames(mgi), rowData(AD_SCE)$hgnc_symbol)
marker_SCE <- AD_SCE[marker_in_sce[!is.na(marker_in_sce)]]
counts(marker_SCE) <- as.matrix(counts(marker_SCE))

marker_SCE <- AD_SCE[unlist(brainMarkers)]
