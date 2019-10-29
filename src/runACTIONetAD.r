require(ACTIONet)

reduced.dim <- 50
k.max <- 20

AD_SCE <- readRDS("data/Processed_10x_AD_data/AD_SCE.rds")
sce.reduced <- reduce.sce(AD_SCE, reduced_dim = reduced.dim)
ACTIONet.out <- run.ACTIONet(sce.reduced, k_max = k.max)

saveRDS(sce.reduced, "data/Processed_10x_AD_data/AD_SCE_reduced.rds")
saveRDS(ACTIONet.out, "data/Processed_10x_AD_data/AD_ACTIONet.rds")
