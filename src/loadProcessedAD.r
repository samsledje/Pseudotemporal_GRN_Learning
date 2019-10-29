require(Matrix)
require(SingleCellExperiment)
require(DropleUtils)
require(biomaRt)
require(scater)
require(scran)

mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Read in counts and metadata and clean metadata
filtered.counts <- readMM("data/Processed_10x_AD_data/filtered_count_matrix.mtx")
filtered.colMetadata <- read.delim("data/Processed_10x_AD_data/filtered_column_metadata.txt")
filtered.colMetadata$pre.cluster <- factor(filtered.colMetadata$pre.cluster)
idMapping = read.table("data/Processed_10x_AD_data/id_mapping.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
filtered.colMetadata$Subject <- lapply(filtered.colMetadata$projid,function(x) { idMapping$Subject[[match(x,idMapping$projid)]] })

# Assign row names (genes) and column names (cell tags)
rownames(filtered.counts) <- readLines("data/Processed_10x_AD_data/filtered_gene_row_names.txt")
colnames(filtered.counts) <- filtered.colMetadata$TAG

# Create tsne array
tsneDims <- filtered.colMetadata[c("tsne1", "tsne2")]

# Create Single Cell Experiment
AD_SCE = SingleCellExperiment(assays = list(counts = filtered.counts),colData=filtered.colMetadata[c("projid","Subject","pre.cluster","broad.cell.type","Subcluster")])
reducedDim(AD_SCE) <- tsneDims
reducedDimNames(AD_SCE) <- "tsne"

# Annotate Genes
rowData(AD_SCE)$hgnc_symbol <- rownames(AD_SCE)
AD_SCE <- getBMFeatureAnnos(AD_SCE, filters="hgnc_symbol",
                            attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "chromosome_name"),
                            dataset = "hsapiens_gene_ensembl")

# Compute Size Factors
AD_SCE <- computeSumFactors(AD_SCE, BPPARAM= MulticoreParam(10))

# Compute log counts
AD_SCE <- normalize(AD_SCE)

# Save File
saveRDS(AD_SCE,"data/Processed_10x_AD_data/AD_SCE_normalized.rds")

