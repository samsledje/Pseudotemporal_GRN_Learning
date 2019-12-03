setwd("/Users/dkaur/Documents/Pseudotemporal_GRN_Learning")

require(bnlearn)
require(dplyr)

cell.type = "Mic"
clusts <- readRDS(paste0("processed-data/",cell.type,"/clusters.rds"))
df <- as.data.frame(readRDS(paste0("processed-data/",cell.type,"/top_gene_counts.rds")))

#pseudotime.sorted <- t(df)[,order(clusts$DiseaseScore)]
#write.table(pseudotime.sorted, file="loremIpsum.csv", sep=",")

K.clusters = 4
N.genes = 100

df <- df[,1:N.genes]
clusts$ScoreClusts = as.numeric(clusts$DiseaseRange)

layers = list()
for (i in seq(K.clusters)) {
  layers[[i]] <- df[clusts$ScoreClusts == i,]
  colnames(layers[[i]]) <- paste0(colnames(layers[[i]]),"_t",i)
}


#mincells = nrow(layers[[1]])
#for (i in 2:K.clusters){
#  mincells = min(mincells,nrow(layers[[i]]))
#}

sample_df <- function(matrix,samplesize){
  rownames(matrix) <- c()
  rawdf <- as.data.frame(matrix)
  rawdf[] <- lapply(rawdf, factor)
  
  # network learning requires > 1 possible value for each var
  sampled <- sample_n(rawdf,samplesize,replace = TRUE)
  sampled <- select_if(sampled,function(x) return(nlevels(x)>1))
  return(sampled)
}


bl <- data.frame()
bndf <- NULL
prevdf <- NULL
prevcols <- c()

for (i in (1:K.clusters)){
  matrix <- layers[[i]]
  newdf <- sample_df(matrix,10000)
  # only allow edges to future time points
  bl <- rbind(bl,
              expand.grid(colnames(newdf),prevcols),
              expand.grid(colnames(newdf),colnames(newdf)))
  #wl <- rbind(wl,
  #expand.grid(colnames(prevdf),colnames(newdf)))
  
  prevcols <- c(prevcols,colnames(newdf))
  if (is.null(bndf)){
    bndf <- newdf
  }
  else{
    bndf <- cbind(bndf,newdf)
  }
  prevdf <- newdf
}

network <- fast.iamb(bndf, blacklist = bl,alpha = .05)
network[["arcs"]]
arc.strength(network,bndf)

