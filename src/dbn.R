setwd("/Users/dkaur/Documents/Pseudotemporal_GRN_Learning")

#setwd("/home/samsledje/Pseudotemporal-GRNs")
#setwd("D:/Drive/PhD/GitHub/Pseudotemporal_GRN_Learning")

require(bnlearn)
require(dplyr)
#require(Rgraphviz)

cell.type = "Mic"
clusts <- readRDS(paste0("processed-data/",cell.type,"/clusters.rds"))
df <- as.data.frame(readRDS(paste0("processed-data/",cell.type,"/top_gene_counts.rds")))

pseudotime.sorted <- t(df)[,order(clusts$DiseaseScore)]
write.table(pseudotime.sorted, file="loremIpsum.csv", sep=",")


K.clusters = 4
N.genes = 100

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
  simdf <- rbn(network, data = df, 1000)
  return(simdf)
}



create_df <- function(){
  bl <- data.frame()
  bndf <- NULL
  prevdf <- NULL
  prevcols <- c()
  
  for (i in (1:K.clusters)){
    matrix <- layers[[i]]
    newdf <- simulate.Observations(matrix)
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
  res <- list("df" = bndf, "bl" = bl)
  return(res)
}

res <- create_df()
bndf <- res$df
bl <- res$bl

# many different methods available for network learning
# network <- pc.stable(bndf,blacklist = bl)
# pc_network <- pc.stable(bndf, blacklist = bl, alpha = .05)[["arcs"]]
# gs_network <- gs(bndf, blacklist = bl, alpha = .05)[["arcs"]]
fi_network <- fast.iamb(bndf, blacklist = bl, alpha = .05)[["arcs"]]

for (i in 1:10){
  fi_network <- fast.iamb(bndf, blacklist = bl, alpha = .05)[["arcs"]]
  
}

# output edges of graph
network[["arcs"]]
#plot(network)
groups <- list()
for (i in seq(K.clusters)) {
  groups[[i]] <- colnames(layers[[i]])
}

intersect()
#graphviz.plot(network,groups=groups)
arc.strength(network,bndf)

