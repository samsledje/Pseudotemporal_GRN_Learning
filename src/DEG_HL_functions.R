
substr.last <- function(string){
  substr(string, nchar(string), nchar(string))
}

filter.sce <- function(sce, min_cells_per_feat = NULL, min_feats_per_cell = NULL, min_umis_per_cell = NULL, max_umis_per_cell = NULL, max_iter = 5){

  i = 0
  repeat{
    prev_dim = dim(sce)
  
    if(is.numeric(min_umis_per_cell)){
      umi_counts = Matrix::colSums(sce@assays[["counts"]])
      sce <- sce[, umi_counts >= min_umis_per_cell]
    }
  
    if(is.numeric(max_umis_per_cell)){
      umi_counts = Matrix::colSums(sce@assays[["counts"]])
      sce <- sce[, umi_counts <= max_umis_per_cell]
    }
  
    if(is.numeric(min_cells_per_feat)){
      cell_counts = Matrix::rowSums(sce@assays[["counts"]] > 0)
      sce <- sce[cell_counts >= min_cells_per_feat, ]
    }
  
    if(is.numeric(min_feats_per_cell)){
      feature_counts = Matrix::colSums(sce@assays[["counts"]] > 0)
      sce <- sce[, feature_counts >= min_feats_per_cell]
    }

    i = i+1
    print(i)
    if( all( dim(sce) == prev_dim) ){break}
  
    if(i == max_iter){
      print(paste('Failed to filter SCE in', max_iter, 'iterations.', sep = " "))
      break
    }
  }
  
  return(sce)
}

filter.sce.by.group <- function(sce, group_list, by = 'group', group_attr = '', min_cells_per_feat = NULL, min_feats_per_cell = NULL, min_umis_per_cell = NULL, max_umis_per_cell = NULL, max_iter = 50){

  remaining_feats = list()
  remaining_cells = list()
  
  if(group_attr == ''){
    stop(print('No filtering criteria given.'))
  } else{
  
    for(group in group_list){
      if(by == 'group'){
        this_sce <- sce[ , sce[[group_attr]] == group ]
      } else if(by == 'batch'){
        this_sce <- sce[ , as.numeric(sce[[group_attr]]) %in% group ]
      } else{
        print('Not a valid subset')
        return(sce)
        }
      
      this_sce <- filter.sce(this_sce, min_cells_per_feat = min_cells_per_feat, 
                             min_umis_per_cell = min_umis_per_cell, 
                             max_umis_per_cell = max_umis_per_cell, 
                             min_feats_per_cell = min_feats_per_cell, 
                             max_iter = max_iter)
      remaining_feats[[length(remaining_feats) + 1]] <- rownames(this_sce)
      remaining_cells[[length(remaining_cells) + 1]] <- colnames(this_sce)
      
    }
    
    remaining_feats <- Reduce(intersect, remaining_feats)
    remaining_cells <- do.call(c, remaining_cells)
    sce <- sce[remaining_feats, remaining_cells]
    return(sce)
    
  }
  
}

import.sce.cellranger <- function(input_path, mtx_file = 'matrix.mtx.gz', gene_annotations = 'features.tsv.gz', sample_annotations = 'barcodes.tsv.gz', sep = '\t', aggr = FALSE, prefilter = FALSE, prefil_params = NULL) {
  require(Matrix)
  require(scran)

  counts = readMM(paste(input_path, mtx_file, sep='/'));
  gene.table = read.table(paste(input_path, gene_annotations, sep='/'), header = FALSE, sep=sep, as.is = TRUE)
  sample_annotations = read.table(paste(input_path, sample_annotations, sep='/'), header = FALSE, sep=sep, as.is = TRUE)

  rownames(counts) = gene.table[, 1]
  colnames(counts) = sample_annotations[, 1]

  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = sample_annotations,
    rowData = gene.table
  )

  if(aggr){
    sce$batch = factor( substr.last( colnames(sce) ) )
  }

  if( prefilter){
    sce <-  do.call(filter.sce, append( list(sce = sce), as.list(prefil_params) ) )
  }

  return(sce)
}


combine.tech.replicates <- function(sce_list){

  cell_list <- do.call(intersect, lapply(sce_list, function(x) colnames(x)))
  feat_list <- do.call(intersect, lapply(sce_list, function(x) rownames(x)))
  if( (length(feat_list) & length(cell_list)) ){

    counts_intersect <- Reduce('+', lapply(sce_list, function(x) x[feat_list, cell_list]@assays[['counts']]))

    common_cell <- do.call(intersect, lapply(sce_list, function(x) colData(x)$V1))
    cell_data <- subset(colData(sce_list[[1]]),  colData(sce_list[[1]])$V1 %in% common_cell)
    rownames(cell_data) <- cell_data$V1

    common_feat <- do.call(intersect, lapply(sce_list, function(x) rowData(x)$V1))
    feat_data <- as.data.frame( rowData(sce_list[[1]])[common_feat, ] )
    rownames(feat_data) <- feat_data$V1

    rownames(counts_intersect) = feat_data[, 1]
    colnames(counts_intersect) = cell_data[, 1]

    sce <- SingleCellExperiment(
      assays = list(counts = counts_intersect),
      colData = cell_data,
      rowData = feat_data
    )

    return(sce)

  }else{
    if( !length(feat_list) ) print('Feature intersection is empty Are row names in the same format?')
    if( !length(cell_list) ) print('No common cell barcodes')
  }

}


reduce.and.batch.with.harmony <- function(sce, reduced_dim = 50, batch_attr = 'batch', harmony_pca = FALSE, reduce_iters = 5, harmony_iters = 10){
  
  sce <- reduce.sce(sce, reduced_dim = reduced_dim, max.iter = reduce_iters)
  sce@reducedDims$S_r <- HarmonyMatrix(sce@reducedDims$S_r, sce[[batch_attr]], do_pca = harmony_pca, max.iter.harmony = harmony_iters)
  
  return(sce)
}

log.norm.sce <- function(sce){
  
  sce.norm = sce
  A = as(sce@assays[["counts"]], 'dgTMatrix')
  cs = Matrix::colSums(A)
  cs[cs == 0] = 1
  B = Matrix::sparseMatrix(i = A@i+1, j = A@j+1, x = log(1 + median(cs)*(A@x / cs[A@j + 1])), dims = dim(A))
  sce.norm@assays[["logcounts"]] = B

  return(sce.norm)
  
}


reduce.norm.sce <- function(sce, reduced_dim = 50, max.iter = 10){
  require(ACTIONet)
  reduction.out = reduceGeneExpression(as(sce@assays[["logcounts"]], 'sparseMatrix'), reduced_dim = reduced_dim, method = 1, iters = max.iter)

  SingleCellExperiment::reducedDim(sce, "S_r") <- t(reduction.out$S_r)
  
  V = reduction.out$V
  colnames(V) = sapply(1:dim(V)[2], function(i) sprintf('PC%d', i))
  
  X = rowData(sce)
  PC.idx = -grep('^PC', colnames(X))
  if(length(PC.idx) > 0)
    X = X[, PC.idx]
  rowData(sce) = cbind(V, X)
  
  return(sce)
}




change.ACTIONet.gene.name.format <- function(ACTIONet, sce = NULL, symType = 'GeneSym', gene.name.vec = NULL){

  if(symType == 'GeneSym'){
    cd = 'V2'
  } else if(symType == 'ENSMBL'){
    cd = 'V1'
  } else {
    print('Invalid gene name format.')
    return(ACTIONet.out)
  }
  # print(cd)
  if( !is.null(gene.name.vec) ){

    rownames(ACTIONet$reconstruct.out$archetype_profile) = gene.name.vec
    rownames(ACTIONet$signature.profile) = gene.name.vec
    rownames(ACTIONet$sce@assays[['counts']]) = gene.name.vec
    return(ACTIONet)

  } else if( !is.null(sce) ){

    rownames(ACTIONet$reconstruct.out$archetype_profile) = rowData(sce)[[cd]]
    rownames(ACTIONet$signature.profile) = rowData(sce)[[cd]]
    # rownames(ACTIONet$sce@assays[['counts']]) = rowData(ACTIONet$sce)[[cd]]
    
    return(ACTIONet)

  } else if( !is.null(ACTIONet$sce) ){
    
    rownames(ACTIONet$reconstruct.out$archetype_profile) = rowData(ACTIONet$sce)[[cd]]
    rownames(ACTIONet$signature.profile) = rowData(ACTIONet$sce)[[cd]]
    rownames(ACTIONet$sce) = rowData(ACTIONet$sce)[[cd]]
    rownames(ACTIONet$sce@assays[['counts']]) = rowData(ACTIONet$sce)[[cd]]
    
    if( !is.null(sce) ){
      cat('ACTIONet.out object contains SingleCellExperiment object. \'sce\' argument ignored.')
    }
    
    return(ACTIONet)
  } else {
    print('Invalid arguments')
    return(ACTIONet)
    }

}

# make.cell.attributes.table <- function(ACTIONet, sce = NULL, annotations, sce_attr = NULL, add_attr = NULL, filter_cells = FALSE, fil_params = NULL){
make.cell.attributes.table <- function(sce, annotations, sce_attr = NULL, add_attr = NULL, filter_cells = FALSE, fil_params = NULL){
  # if( !is.null(sce) ){
  #   sce = sce
  # } else {
  #   sce = ACTIONet$sce
  # }

  cell.attributes = data.frame("Barcode" = as.character(colnames(sce)), 
                            "sce_idx" = colnames(sce), 
                            "CellType" = annotations$Labels, 
                            "Confidence" = annotations$Labels.confidence, 
                            stringsAsFactors = FALSE)

  ########### THESE 2 ARE REDUNDANT
  if( !is.null(sce_attr) ){
    for(i in sce_attr){
        cell.attributes[[i]] <- sce[[i]]
    }
  }

  if( class(add_attr) == "list" ){
    for(this_attr in names(add_attr)){
      cell.attributes[[this_attr]] <- add_attr[[this_attr]]
    }
    
  }
########################################
  
  # if(filter_cells){
  #     sce <-  do.call(filter.annotated.cells, append( list(cell.attributes = cell.attributes), as.list(fil_params) ) )
  # }

  return(cell.attributes)
}

filter.annotated.cells <- function(cell.attributes, conf_th = NULL, min_cells = NULL, keep_only= NULL, excluded_types = NULL){

  if( !is.null(conf_th) ){
    cell.attributes <- subset( cell.attributes, (cell.attributes$Confidence >= conf_th) )
  }

  if( !is.null(min_cells) ){
    ct_above_th <- names(subset(table(cell.attributes$CellType), table(cell.attributes$CellType) > min_cells) )
    cell.attributes <- subset(cell.attributes, cell.attributes$CellType %in% ct_above_th)
  }

  if( !is.null(keep_only) ){
    cell.attributes <- subset(cell.attributes, cell.attributes$CellType %in% keep_only)
  } else if( !is.null(excluded_types) ){
    cell.attributes <- subset(cell.attributes, !(cell.attributes$CellType %in% excluded_types) )
  }

for(attr in names(cell.attributes)){
  if( is.factor(cell.attributes[[attr]]) ){
    cell.attributes[[attr]] <- droplevels(cell.attributes[[attr]])
  }
}

  return(cell.attributes)
}


remove.cells.from.annotations <- function(annotations, which_cells, keep_cells = FALSE){
  
  if(keep_cells == TRUE){
    annotations$Labels <- annotations$Labels[which_cells]
    annotations$Labels.confidence <- annotations$Labels.confidence[which_cells]
    annotations$cell2celltype.mat <- annotations$cell2celltype.mat[which_cells,]
    annotations$imputed.marker.expression <- annotations$imputed.marker.expression[which_cells,]
    if( !is.null( annotations$clusters ) ){
      annotations$clusters <- annotations$clusters[which_cells]
    }
  } else{
    annotations$Labels <- annotations$Labels[-which_cells]
    annotations$Labels.confidence <- annotations$Labels.confidence[-which_cells]
    annotations$cell2celltype.mat <- annotations$cell2celltype.mat[-which_cells,]
    annotations$imputed.marker.expression <- annotations$imputed.marker.expression[-which_cells,]
    if( !is.null( annotations$clusters ) ){
      annotations$clusters <- annotations$clusters[-which_cells]
    }
    
  }
  
  return(annotations)
  
}

# 
# remove.cells.from.annotations <- function(sce, annotations, hashes, keep_cells = FALSE){
#   
#   which_cells = which(!(sce$cell.hashtag %in% hashes))
#   
#   if(keep_cells == TRUE){
#     annotations$Labels <- annotations$Labels[which_cells]
#     annotations$Labels.confidence <- annotations$Labels.confidence[which_cells]
#     annotations$cell2celltype.mat <- annotations$cell2celltype.mat[which_cells,]
#     annotations$imputed.marker.expression <- annotations$imputed.marker.expression[which_cells,]
#     if( !is.null( annotations$clusters ) ){
#       annotations$clusters <- annotations$clusters[which_cells]
#     }
#   } else{
#     annotations$Labels <- annotations$Labels[-which_cells]
#     annotations$Labels.confidence <- annotations$Labels.confidence[-which_cells]
#     annotations$cell2celltype.mat <- annotations$cell2celltype.mat[-which_cells,]
#     annotations$imputed.marker.expression <- annotations$imputed.marker.expression[-which_cells,]
#     if( !is.null( annotations$clusters ) ){
#       annotations$clusters <- annotations$clusters[-which_cells]
#     }
#     
#   }
#   
#   return(annotations)
#   
# }
