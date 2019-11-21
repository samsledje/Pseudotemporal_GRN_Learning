setwd("/home/samsl/Pseudotemporal-GRNs")

files <- c("data/Oli/noAD_Oli.rds",
           "data/Oli/possibleAD_Oli.rds",
           "data/Oli/probableAD_Oli.rds",
           "data/Oli/definiteAD_Oli.rds"
           )

library(bnlearn)
library(dplyr)

cell.type = "Mic"
clusts <- readRDS(paste0("data/",cell.type,"/clusters.rds"))
df <- readRDS(paste0("data/",cell.type,"/top_gene_counts.rds"))
K.clusters = 4
clusts$ScoreClusts = as.numeric(clusts$DiseaseRange)

layers = list()
names(layers) <- seq(K.clusters)
for (i in seq(K.clusters)) {
  layers[[i]] <- t(df[clusts$ScoreClusts == i,])
  rownames(layers[[i]]) <- paste0(rownames(layers[[i]]),"_t",i)
}

mat_to_df <- function(matrix,time){
  # matrix should have genes as rows, cells as columns
  transposed <- t(matrix)
  rownames(transposed) <- c()
  rawdf <- as.data.frame(transposed)
  rawdf[] <- lapply(rawdf, factor)
  
  # remove duplicate columns
  df <- rawdf[, !duplicated(colnames(rawdf), fromLast = TRUE)] 
  
  # network learning requires > 1 possible value for each var
  df <- select_if(df,function(x) return(nlevels(x)>1))
  
  # learn from gene expression at each time step to generate data
  network <- hc(df)
  simdf <- rbn(network, data = df, 10000)
  colnames(simdf) <- paste(colnames(simdf), time, sep = "_")
  return(simdf)
}

bl <- data.frame()
bndf <- NULL
prevcols <- c()

for (i in (1:length(files))){
  matrix <- readRDS(files[i])
  newdf <- mat_to_df(matrix,paste("t",i,sep=""))
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
graphviz.plot(network,groups=list(t1=rownames(t1),t2=rownames(t2),t3=rownames(t3)))
