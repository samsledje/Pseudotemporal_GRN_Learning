setwd("/home/samsl/Pseudotemporal-GRNs")

files <- c("data/Oli/noAD_Oli.rds",
           "data/Oli/possibleAD_Oli.rds",
           "data/Oli/probableAD_Oli.rds",
           "data/Oli/definiteAD_Oli.rds"
           )

library(bnlearn)
library(dplyr)

clusts <- readRDS(paste0("data/",cell.type,"/clusters.rds"))
df <- readRDS(paste0("data/",cell.type,"/top_gene_counts.rds"))
cell.type = "Mic"
K.clusters = 3

t1 = t(df[clusts$Cluster == 1,])
rownames(t1) <- paste0(rownames(t1),"_t1")
t2 = t(df[clusts$Cluster == 2,])
rownames(t2) <- paste0(rownames(t2),"_t2")
t3 = t(df[clusts$Cluster == 3,])
rownames(t3) <- paste0(rownames(t3),"_t3")


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
network <- pc.stable(bndf,blacklist = bl)
network <- gs(bndf, blacklist = bl)
network <- iamb(bndf, blacklist = bl)

# output edges of graph
network[["arcs"]]
graphviz.plot(network,groups=list(t1=rownames(t1),t2=rownames(t2),t3=rownames(t3)))
