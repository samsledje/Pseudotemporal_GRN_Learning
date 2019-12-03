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
pseudotime.sorted <- t(gene.df)[,order(clusts$DiseaseScore)]
write.csv(pseudotime.sorted, file=paste0("processed-data/",cell.type,"/GRNVBEM.csv"))

# Divide cells into layers based on clusters
clusts$ScoreClusts = as.numeric(clusts$DiseaseRange)
K.clusters = length(unique(clusts$ScoreClusts))

layers = list()
for (i in seq(K.clusters)) {
  layers[[i]] <- gene.df[clusts$ScoreClusts == i,]
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
  prevcols <- c(prevcols,colnames(newdf))
  if (is.null(bndf)){
    bndf <- newdf
  }
  else{
    bndf <- cbind(bndf,newdf)
  }
}

# many different methods available for network learning
#network <- pc.stable(bndf,blacklist = bl)
#network <- gs(bndf, blacklist = bl)
network <- iamb(bndf, blacklist = bl)

# output edges of graph
network[["arcs"]]
plot(network)
groups <- list()
for (i in seq(K.clusters)) {
  groups[[i]] <- names(network$nodes)[grep(paste0(".*_t",i),names(network$nodes))]
  #groups[[i]] <- colnames(layers[[i]])
}
graphviz.plot(network,groups=groups)

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
