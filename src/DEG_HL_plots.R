
plot.ACTIONplot <- function(ACTIONet.out, color.attr = NA, transparency.attr = NA, size.attr = NA, cex = 2, CPal = "d3", legend.pos = "bottomright", Legend = TRUE, Title = NULL) {
  require(ggplot2)
  require(ggpubr)

  if(is.igraph(ACTIONet.out))
    ACTIONet = ACTIONet.out
  else
    ACTIONet = ACTIONet.out$ACTIONet

  nV = length(V(ACTIONet))
  coor = cbind(V(ACTIONet)$x, V(ACTIONet)$y)

  # if( length(color.attr) == length(V(ACTIONet)) ){
  #   vCol = color.attr
  # }
  # else
  if(is.numeric(color.attr)) {
    v = sort(unique(color.attr))
    Annot = as.character(v)
    Pal = ggpubr::get_palette(CPal, length(Annot))
    vCol = Pal[color.attr]
  }
  else if(is.factor(color.attr)) {
    Annot = levels(color.attr)
    Pal = ggpubr::get_palette(CPal, length(Annot))
    vCol = Pal[color.attr]
  }
  else if(is.character(color.attr)) {
    color.attr = factor(color.attr, levels = sort(unique(color.attr)))
    Annot = levels(color.attr)
    Pal = ggpubr::get_palette(CPal, length(Annot))
    vCol = Pal[color.attr]
  }
  else {
    vCol = V(ACTIONet)$color
    Annot = NA
    Pal = NULL
  }


  HSV = rgb2hsv(col2rgb(vCol))
  HSV[3, ] = HSV[3, ]*0.7
  vCol.border = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))

  if(is.numeric(transparency.attr)) {
    trans.factor = 1 / ( 1 + exp(-scale(transparency.attr)) )
    trans.factor = trans.factor ^ 2

    trans.factor[transparency.attr == 0] = 0

    vCol = ggplot2::alpha(vCol, trans.factor)
    vCol.border = ggplot2::alpha(vCol.border, trans.factor)
  }

  # if(is.numeric(size.attr)) {
  #   size.factor = 1 / ( 1+ exp(-scale(transparency.attr)) )
  #   size.factor = trans.factor ^ 2

  if(is.numeric(size.attr)) {
    size.factor = 1 / ( 1 + exp(-scale(size.attr)) )
    size.factor = size.factor ^ 2

    cex = cex * size.factor

  }

  sketch.graph = ACTIONet
  sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
  V(sketch.graph)$size = cex
  V(sketch.graph)$color = vCol
  V(sketch.graph)$frame.color = vCol.border


  plot(sketch.graph, vertex.label=NA, layout=coor, main = Title)
  if(Legend)
    legend(legend.pos, legend = Annot, fill=Pal, cex = 0.5)

}

plot.ACTIONplot.interactive <- function(ACTIONet.out, sce, labels = NULL, top.arch.genes = 10, blacklist.pattern = '\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP', marker.per.cell = 5, node.size = 2, CPal = "Spectral", show.legend = TRUE, annotate.cells = FALSE, opacity = 1.0, title = 'ACTIONet', Alt_Text = NULL) {

  require(plotly)
  require(ACTIONet)

  signature.profile = ACTIONet.out$signature.profile[, ACTIONet.out$core.out$core.archs]

  filtered.row.mask = grepl(blacklist.pattern, toupper(rownames(sce)))
  signature.profile = signature.profile[!filtered.row.mask, ]

  sorted.top.genes = apply(signature.profile, 2, function(x) rownames(signature.profile)[order(x, decreasing = TRUE)[1:top.arch.genes]])

  selected.genes = sort(unique(as.character(sorted.top.genes)))

  if( !is.null(Alt_Text) ){
    node.annotations = Alt_Text
  } else if(annotate.cells == TRUE) {
    #imputed.markers = impute.genes.using.ACTIONet(ACTIONet.out, sce, selected.genes, prune = TRUE, rescale = TRUE)
    imputed.markers = t(ACTIONet.out$signature.profile[selected.genes, ACTIONet.out$core.out$core.archs] %*% ACTIONet.out$core.out$H)
    node.annotations = apply(imputed.markers, 1, function(x) {
      top.genes = colnames(imputed.markers)[order(x, decreasing = TRUE)[1:marker.per.cell]]
      label = paste(top.genes, collapse = ';')
      return(label)
    })
  } else {
    node.annotations = ''
  }
  if(is.null(labels)) {
    labels = as.factor(array(1, dim(sce)[2]))
  } else if(is.character(labels) | is.numeric(labels)) {
    labels = factor(labels, levels = sort(unique(labels)))
  }

  # Setup visualization parameters
  sketch.graph = ACTIONet.out$ACTIONet
  sketch.graph = delete.edges(sketch.graph, E(sketch.graph))

  node.data <- get.data.frame(sketch.graph, what="vertices")
  edge.data <- get.data.frame(sketch.graph, what="edges")

  Nv <- dim(node.data)[1]
  Ne <- dim(edge.data)[1]

  edge_shapes <- list()

  # Adjust parameters
  node.data$size = node.size

  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)

  Pal = ggpubr::get_palette(CPal, length(levels(labels)))
  names(Pal) = levels(labels)
  node.data$type = labels
  network <- plot_ly(node.data, x = ~x, y = ~y, opacity = opacity, color = ~type, colors = Pal, marker = list(size = ~size, opacity = 1.0, alpha = 1, line = list(width = 0.1*node.size, alpha = 0.5, color = 'rgb(0, 0, 0)')), text = node.annotations, mode = "markers", type = 'scatter', hoverinfo = "text")

  p <- plotly::layout(
    network,
    title = title,
    shapes = edge_shapes,
    xaxis = axis,
    yaxis = axis,
    showlegend=show.legend
  )
}


plot.ACTIONet.3D <- function(ACTIONet.out, color.attr = NA, transparency.attr = NA, size.attr = NA, node.size = 0.2, CPal = "d3") {
  require(ggplot2)
  require(ggpubr)
  require(threejs)

  if(is.igraph(ACTIONet.out))
    ACTIONet = ACTIONet.out
  else
    ACTIONet = ACTIONet.out$ACTIONet

  nV = length(V(ACTIONet))
  coor = cbind(V(ACTIONet)$x3D, V(ACTIONet)$y3D, V(ACTIONet)$z3D)

  # Annot = NA
  # if(is.numeric(color.attr)) {
  #   color.attr = as.character(color.attr)
  #   color.attr = factor(color.attr, levels = sort(unique(color.attr)))
  #   Annot = levels(color.attr)
  #
  #   if(!is.na(CPal)) {
  #     Pal = CPal[1:length(Annot)]
  #   } else {
  #     Pal = ggpubr::get_palette(CPal, length(Annot))
  #   }
  #   names(Pal) = Annot
  # }
  # else if(is.factor(color.attr)) {
  #   Annot = levels(color.attr)
  #   if(!is.na(CPal)) {
  #     Pal = CPal[1:length(Annot)]
  #   } else {
  #     Pal = ggpubr::get_palette(CPal, length(Annot))
  #   }
  #   names(Pal) = Annot
  # }
  # else if(is.character(color.attr)) {
  #   color.attr = factor(color.attr, levels = sort(unique(color.attr)))
  #   Annot = levels(color.attr)
  #   if(!is.na(CPal)) {
  #     Pal = CPal[1:length(Annot)]
  #   } else {
  #     Pal = ggpubr::get_palette(CPal, length(Annot))
  #   }
  #   names(Pal) = Annot
  # }
  # if(length(Annot) > 1) {
  #   vCol = as.character(Pal[color.attr])
  # } else {
  #   vCol = V(ACTIONet)$color
  # }

  if(is.numeric(color.attr)) {
    v = sort(unique(color.attr))
    Annot = as.character(v)
    Pal = ggpubr::get_palette(CPal, length(Annot))
    vCol = Pal[color.attr]
  }
  else if(is.factor(color.attr)) {
    Annot = levels(color.attr)
    Pal = ggpubr::get_palette(CPal, length(Annot))
    vCol = Pal[color.attr]
  }
  else if(is.character(color.attr)) {
    color.attr = factor(color.attr, levels = sort(unique(color.attr)))
    Annot = levels(color.attr)
    Pal = ggpubr::get_palette(CPal, length(Annot))
    vCol = Pal[color.attr]
  }
  else {
    vCol = V(ACTIONet)$color
    Annot = NA
  }

  HSV = rgb2hsv(col2rgb(vCol))
  HSV[3, ] = HSV[3, ]*0.7
  vCol.border = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))

  if(is.numeric(transparency.attr)) {
    vCol = ggplot2::alpha(vCol, 1 / (1 + exp(-scale(transparency.attr))))
    vCol.border = ggplot2::alpha(vCol.border, 1 / (1 + exp(-scale(transparency.attr))))
  }

  if(is.numeric(size.attr)) {
    node.size = node.size * (1 / (1 + exp(-scale(size.attr))))
  }


  scatterplot3js(x = coor[, 1], y = coor[, 2], z = coor[, 3], axis.scales = FALSE, size = node.size, axis = F, grid = F, color = vCol, stroke = vCol.border, bg='black')
}

plot.ACTIONet.3Dscatter <- function(ACTIONet.out, color.attr = NA, transparency.attr = NA, size.attr = NA, node.size = 0.2, CPal = "d3") {
  require(ggplot2)
  require(ggpubr)
  require(threejs)

  if(is.igraph(ACTIONet.out))
    ACTIONet = ACTIONet.out
  else
    ACTIONet = ACTIONet.out$ACTIONet

  nV = length(V(ACTIONet))
  coor = cbind(V(ACTIONet)$x3D, V(ACTIONet)$y3D, V(ACTIONet)$z3D)

  if(is.numeric(color.attr)) {
    v = sort(unique(color.attr))
    Annot = as.character(v)
    Pal = ggpubr::get_palette(CPal, length(Annot))
    vCol = Pal[color.attr]
  }
  else if(is.factor(color.attr)) {
    Annot = levels(color.attr)
    Pal = ggpubr::get_palette(CPal, length(Annot))
    vCol = Pal[color.attr]
  }
  else if(is.character(color.attr)) {
    color.attr = factor(color.attr, levels = sort(unique(color.attr)))
    Annot = levels(color.attr)
    Pal = ggpubr::get_palette(CPal, length(Annot))
    vCol = Pal[color.attr]
  }
  else {
    vCol = V(ACTIONet)$color
    Annot = NA
  }

  HSV = rgb2hsv(col2rgb(vCol))
  HSV[3, ] = HSV[3, ]*0.7
  vCol.border = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))

  if(is.numeric(transparency.attr)) {
    vCol = ggplot2::alpha(vCol, 1 / (1 + exp(-scale(transparency.attr))))
    vCol.border = ggplot2::alpha(vCol.border, 1 / (1 + exp(-scale(transparency.attr))))
  }

  if(is.numeric(size.attr)) {
    node.size = node.size * (1 / (1 + exp(-scale(size.attr))))
  }


  scatterplot3d(x = coor[, 1], y = coor[, 2], z = coor[, 3], axis = F, grid = T, color = vCol, bg='black', pch = 16)
}
