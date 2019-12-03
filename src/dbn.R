setwd("/Users/dkaur/Documents/Pseudotemporal_GRN_Learning")

library(bnlearn)
library(dplyr)

cell.type = "DataMatrices"
clusts <- readRDS(paste0("data/",cell.type,"/clusters.rds"))
df <- readRDS(paste0("data/",cell.type,"/top_gene_counts.rds"))

K.clusters = 4
df <- df[,1:25]
clusts$ScoreClusts = as.numeric(clusts$DiseaseRange)

layers = list()
for (i in seq(K.clusters)) {
  layers[[i]] <- df[clusts$ScoreClusts == i,]
  colnames(layers[[i]]) <- paste0(colnames(layers[[i]]),"_t",i)
}

sim_timestep <- function(matrix){
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
bndf <- NULL
prevcols <- c()

for (i in (1:K.clusters)){
  matrix <- layers[[i]]
  newdf <- sim_timestep(matrix)
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

# whitelist edges between time steps if weighted 

# many different methods available for network learning
# network <- pc.stable(bndf,blacklist = bl)
# network <- gs(bndf, blacklist = bl)
network <- fast.iamb(bndf, blacklist = bl)

# output edges of graph
network[["arcs"]]
plot(network)
graphviz.plot(network,groups=list(t1=rownames(t1),t2=rownames(t2),t3=rownames(t3)))
