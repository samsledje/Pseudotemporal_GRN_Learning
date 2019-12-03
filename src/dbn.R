setwd("/home/samsl/Pseudotemporal-GRNs")
#setwd("D:/Drive/PhD/GitHub/Pseudotemporal_GRN_Learning")

require(bnlearn)
require(dplyr)
require(Rgraphviz)

read.gmt <- function(file) {
  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  name = unlist(geneSetDB[1])
  desc = strsplit(geneSetDB[2][[1]],"> ")[[1]][2]
  geneSetDB = unlist(geneSetDB[-(1:2)])
  geneSetDB = list("name"=name,"description"=desc,"genes"=geneSetDB)
  return(geneSetDB)
}

# Select which data to load
cell.type = "Mic"
clusts <- readRDS(paste0("processed-data/",cell.type,"/clusters.rds"))
df <- as.data.frame(readRDS(paste0("processed-data/",cell.type,"/top_gene_counts.rds")))
K.clusters = 4

# Load gene lists for network generation
marker.root <- "processed-data/markers/"
kegg <- read.gmt(paste0(marker.root,"KEGG_AD.gmt"))
adLib <- list("name"="library","description"="Library in AD SCE exeriment","genes"=colnames(df))

blalock_incpt_dn <- read.gmt(paste0(marker.root,"BLALOCK_INCPT_AD_DN.gmt"))
blalock_incpt_up <- read.gmt(paste0(marker.root,"BLALOCK_INCPT_AD_UP.gmt"))
blalock_dn <- read.gmt(paste0(marker.root,"BLALOCK_AD_DN.gmt"))
blalock_up <- read.gmt(paste0(marker.root,"BLALOCK_AD_UP.gmt"))

gene.sets <- list("kegg"=kegg$genes,"incpt_dn"=blalock_incpt_dn$genes,"incpt_up"=blalock_incpt_up$genes,"dn"=blalock_dn$genes,"up"=blalock_up$genes)
mat <- make_comb_mat(gene.sets, mode = "distinct")
UpSet(mat,comb_order = order(-comb_size(mat)))

# Select subset of genes
gene.df <- df[,intersect(adLib$genes,kegg$genes)]

# Write output for GRNVBEM
pseudotime.sorted <- t(df)[,order(clusts$DiseaseScore)]
write.csv(pseudotime.sorted, file=paste0("processed-data/",cell.type,"/GRNVBEM.csv"))

# Divide cells into layers based on clusters
N.genes = 25
df <- df[,1:N.genes]
clusts$ScoreClusts = as.numeric(clusts$DiseaseRange)

layers = list()
for (i in seq(K.clusters)) {
  layers[[i]] <- df[clusts$ScoreClusts == i,]
  colnames(layers[[i]]) <- paste0(colnames(layers[[i]]),"_t",i)
}

# Rows = cells, Columns = genes

simulate.Observations <- function(matrix){
  rownames(matrix) <- c()
  rawdf <- as.data.frame(matrix)
  rawdf[] <- lapply(rawdf, factor)
  
  # network learning requires > 1 possible value for each var
  df <- select_if(rawdf,function(x) return(nlevels(x)>1))
  
  # learn from gene expression at each time step to generate data
  network <- hc(df)
  simdf <- rbn(network, data = df, 10000)
  return(simdf)
}

bl <- data.frame()
wl <- data.frame()
bndf <- NULL
prevcols <- c()

for (i in (1:K.clusters)){
  matrix <- layers[[i]]
  newdf <- simulate.Observations(matrix)
  
  # only allow edges to future time points
  bl <- rbind(bl,
              expand.grid(colnames(newdf),prevcols),
              expand.grid(colnames(newdf),colnames(newdf)))
#  wl <- rbind(wl,
#              expand.grid(prevcols, colnames(newdf)))
  prevcols <- c(prevcols,colnames(newdf))
  if (is.null(bndf)){
    bndf <- newdf
  }
  else{
    bndf <- cbind(bndf,newdf)
  }
}

colnames(wl) <- c("from","to")

# many different methods available for network learning
#network <- pc.stable(bndf,blacklist = bl)
#network <- gs(bndf, blacklist = bl)
network <- iamb(bndf, blacklist = bl)

# output edges of graph
network[["arcs"]]
plot(network)
groups <- list()
for (i in seq(K.clusters)) {
  groups[[i]] <- colnames(layers[[i]])
}
graphviz.plot(network,groups=groups)
