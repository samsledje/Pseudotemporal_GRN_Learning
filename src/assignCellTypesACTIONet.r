library(SingleCellExperiment)
library(ACTIONet)

AD_SCE <- readRDS("data/Processed_10x_AD_data/AD_SCE_normalized.rds")
AD_AN <- readRDS("data/Processed_10x_AD_data/AD_ACTIONet.rds")
brainMarkers <- readRDS("data/markers/Brain_markers.RDS")$Celltype

AD_Annotation_ACTIONet <- annotate.cells.using.markers(AD_AN, AD_SCE, brainMarkers)

saveRDS(AD_Annotation_ACTIONet, "data/Processed_10x_AD_data/AD_ACTIONet_Types.rds")