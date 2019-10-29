AD_SCE <- readRDS("data/Processed_10x_AD_data/AD_SCE_normalized.rds")
AD_AN <- readRDS("data/Processed_10x_AD_data/AD_ACTIONet.rds")
ACTIONetCellTypes <- readRDS("data/Processed_10x_AD_data/AD_ACTIONet_Types.rds")
CellAssignCellTypes <- readRDS("data/Processed_10x_AD_data/AD_CellAssign_Types.rds")
brainMarkers <- readRDS("data/markers/Brain_markers.RDS")$Celltype
