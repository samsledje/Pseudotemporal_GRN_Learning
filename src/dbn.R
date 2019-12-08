#setwd("/home/samsl/Pseudotemporal-GRNs")
#setwd("D:/Drive/PhD/GitHub/Pseudotemporal_GRN_Learning")
setwd("/Users/dkaur/Documents/Pseudotemporal_GRN_Learning")

require(bnlearn)
require(dplyr)
require(Rgraphviz)
require(data.table)
require(ComplexHeatmap)

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

# Load gene lists for network generation
marker.root <- "processed-data/markers/"
kegg <- read.gmt(paste0(marker.root,"KEGG_AD.gmt"))
adLib <- list("name"="library","description"="Library in AD SCE exeriment","genes"=colnames(df))
blalock_incpt_dn <- read.gmt(paste0(marker.root,"BLALOCK_INCPT_AD_DN.gmt"))
blalock_incpt_up <- read.gmt(paste0(marker.root,"BLALOCK_INCPT_AD_UP.gmt"))
blalock_dn <- read.gmt(paste0(marker.root,"BLALOCK_AD_DN.gmt"))
blalock_up <- read.gmt(paste0(marker.root,"BLALOCK_AD_UP.gmt"))
blalock <- list("name"="blalock","description"="All blalock","genes"=c(blalock_incpt_dn$genes, blalock_incpt_up$genes, blalock_dn$genes, blalock_up$genes))
gene.sets <- list("library"=adLib,"kegg"=kegg$genes,"incpt_dn"=blalock_incpt_dn$genes,"incpt_up"=blalock_incpt_up$genes,"dn"=blalock_dn$genes,"up"=blalock_up$genes,"blalock"=blalock)

# Select subset of genes
gene.df <- df[,intersect(adLib$genes,kegg$genes)]

# Divide cells into layers based on clusters
clusts$ScoreClusts = as.numeric(clusts$DiseaseRange)
K.clusters = length(unique(clusts$ScoreClusts))

layers = list()
for (i in seq(K.clusters)) {
  layers[[i]] <- gene.df[clusts$ScoreClusts == i,]
  colnames(layers[[i]]) <- paste0(colnames(layers[[i]]),"_t",i)
}
gene.names <- unlist(list(cbind(unlist(lapply(layers,colnames)))))

mat <- NULL
for (i in seq(K.clusters)) {
  layer.mean <- apply(layers[[i]],2,mean)
  if (is.null(mat)) {
    mat <- data.table(layer.mean)
  } else {
    mat <- cbind(mat, data.table(layer.mean))
  }
}

# Find genes which cluster together with high change in expression across layers
avg.expression <- data.frame(mat)
rownames(avg.expression) <- names(gene.df)
colnames(avg.expression) <- c(1,2,3,4)
ht = Heatmap(avg.expression, cluster_columns = FALSE)
top.difeq.genes <- rownames(avg.expression)[row_order(ht)][1:30]

chosen.genes <- top.difeq.genes
dbn.df <- gene.df[,chosen.genes]

# Recalculate layers with only chosen genes
reduced.layers = list()
for (i in seq(K.clusters)) {
  reduced.layers[[i]] <- dbn.df[clusts$ScoreClusts == i,]
  colnames(reduced.layers[[i]]) <- paste0(colnames(reduced.layers[[i]]),"_t",i)
}
gene.names <- unlist(list(cbind(unlist(lapply(reduced.layers,colnames)))))


# Rows = cells, Columns = genes
sim.df <- function(matrix){
  rownames(matrix) <- c()
  rawdf <- as.data.frame(matrix)
  rawdf[] <- lapply(rawdf, as.numeric)
  
  # network learning requires > 1 possible value for each var
  #df <- select_if(rawdf,function(x) return(nlevels(x)>1))
  
  # learn from gene expression at each time step to generate data
  network <- hc(rawdf)
  simdf <- rbn(network, data = rawdf, 1000)
  return(simdf)
}

sample.df <- function(matrix,samplesize){
  rownames(matrix) <- c()
  rawdf <- as.data.frame(matrix)
  rawdf[] <- lapply(rawdf, as.numeric)
  sampled <- sample_n(rawdf,samplesize,replace = TRUE)
  #sampled <- select_if(sampled,function(x) return(nlevels(x)>1))
  return(sampled)
}

BOOTSTRAP_SAMPLES = FALSE

create_df <- function(){
  bl <- data.frame()
  bndf <- NULL
  prevdf <- NULL
  prevcols <- c()
  
  for (i in (1:K.clusters)){
    matrix <- reduced.layers[[i]]
    
    if (BOOTSTRAP_SAMPLES) {
      newdf <- sample.df(matrix,2000)
    } else {
      newdf <- sim.df(matrix)
    }

    # only allow edges to future time points
    bl <- rbind(bl,
                expand.grid(colnames(newdf),prevcols),
                expand.grid(colnames(newdf),colnames(newdf)))
  
    prevcols <- c(prevcols,colnames(newdf))
    if (is.null(bndf)){
      bndf <- newdf
    }
    else{
      bndf <- cbind(bndf,newdf)
    }
    prevdf <- newdf
  }
  res <- list("df" = bndf, "bl" = bl)
  return(res)
}

# many different methods available for network learning
# network <- pc.stable(bndf,blacklist = bl)
# pc_network <- pc.stable(bndf, blacklist = bl, alpha = .05)[["arcs"]]
# gs_network <- gs(bndf, blacklist = bl, alpha = .05)[["arcs"]]

# -------------------------- #
set.seed(12345)
ITERS = 10
#ALPHA = 0.22
ALPHA = 0.05
interactions <- data.frame()

for (i in 1:ITERS){
  print(i)
  res <- create_df()
  bndf <- res$df
  bl <- res$bl
  fi_network <- fast.iamb(bndf, blacklist = bl)
  interactions <- rbind(interactions,fi_network[["arcs"]])
}


View(interactions[duplicated(interactions),])
nrow(interactions[duplicated(interactions),])
# -------------------------- #



# SIF
write.sif <- function(file,network,df) {
  arcs = arc.strength(network, bndf)
  arcs$relation = rep.int("->",dim(arcs)[1])
  arcs$parentTime = apply(as.matrix(arcs$from), MARGIN=1, FUN=function (x) {as.numeric(strsplit(x,"_t")[[1]][2])})
  arcs = arcs[,c("from","relation","to","strength")]
  names(arcs) = c("Parent","-","Child","Weight")
  write.table(arcs,file=file, row.names = FALSE, quote = FALSE,sep=" ")
  nodeTimes = list("names"=names(network$nodes),"times"=apply(as.matrix(unlist(names(network$nodes))), MARGIN=1, FUN=function (x) {as.numeric(strsplit(x,"_t")[[1]][2])}))
  write.table(nodeTimes,paste0(file,".timetable"), row.names = FALSE, quote = FALSE, sep=" ")
}

write.sif(paste0("processed-data/",cell.type,"/DBN.sif"), network, df)
