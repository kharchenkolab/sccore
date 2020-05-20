#' @import igraph
NULL

#' Collapse Graph PAGA
#' 
#' @description Collapse graph using PAGA 1.2 algorithm
#' @param graph graph to be collapsed
#' @param groups factor on vertices describing cluster assignment (can specify integer vertex ids, or character vertex names which will be matched)
#' @param linearize should normally be always `TRUE` (default=TRUE)
#' @param winsorize winsorize final connectivity statistics value. (default=FALSE) Note: Original PAGA has it always `TRUE`,
#'   but in this case there is no way to distinguish level of connectivity for highly connected groups. 
collapseGraphPaga <- function(graph, groups, linearize=TRUE, winsorize=FALSE) {

  if ((!(is(graph, "Matrix") || is(graph, "matrix")) || ncol(graph) != nrow(graph)) && !igraph::is.igraph(graph)) {
    stop("Unknown graph format. Only adjacency matrix or igraph are supported")
  }

  if (!igraph::is.igraph(graph)) {
    ones <- graph %>% as('dgTMatrix') %>% as('symmetricMatrix')
    ones@x <- rep(1, length(ones@x))
    graph <- igraph::graph_from_adjacency_matrix(ones, mode='directed', weighted=TRUE)
  } else {
    graph %<>% igraph::as.directed()
    igraph::edge.attributes(graph)$weight <- rep(1, igraph::gsize(graph))
  }

  if (!is.null(igraph::V(graph)$name)) {
    groups %<>% .[igraph::V(graph)$name]
    graph %<>% igraph::remove.vertex.attribute("name")
  }

  groups %<>% .[!is.na(.)] %>% as.factor() %>% droplevels()
  if (length(groups) != length(igraph::V(graph))){
    stop("Groups must be provided for all graph vertices")
  }

  cluster.names <- levels(groups)
  groups %<>% as.integer()

  vc <- igraph::make_clusters(graph, membership = groups, algorithm = 'conos', merges = NULL, modularity = FALSE)
  ns <- igraph::sizes(vc)
  n <- sum(ns)

  es.inner.cluster <- rep(0, length(vc))
  es.inner.cluster[as.integer(names(igraph::groups(vc)))] <- sapply(igraph::groups(vc), function(s)
    igraph::induced_subgraph(graph, s, impl = 'copy_and_delete') %>% igraph::gsize() %>% `*`(2))

  inter.es <- igraph::contract(graph, vc$membership) %>%
    igraph::simplify(remove.multiple=FALSE, remove.loops=TRUE, edge.attr.comb=igraph::igraph_opt("weight")) %>%
    igraph::as_adj(attr="weight") %>% as('dgTMatrix')

  es <- es.inner.cluster + Matrix::rowSums(inter.es)

  inter.es %<>% `+`(Matrix::t(.))
  expected.random.null.vals <- (es[inter.es@i + 1]*ns[inter.es@j + 1] + es[inter.es@j + 1]*ns[inter.es@i + 1])/(n - 1)
  if (!linearize) {
    sd.random.null.vals <- (es[inter.es@i + 1]*ns[inter.es@j + 1]*(n-ns[inter.es@j + 1]-1) +
                              es[inter.es@j + 1]*ns[inter.es@i + 1])*(n-ns[inter.es@i + 1]-1) / ((n - 1)^2)
    scaled.values <- ifelse(sd.random.null.vals < 1e-10, 0, (inter.es@x - expected.random.null.vals)/sd.random.null.vals) %>% as.numeric()
  } else {
    scaled.values <- ifelse(expected.random.null.vals < 1e-10, 1, inter.es@x / expected.random.null.vals) %>% as.numeric()
  }

  if (!is.logical(winsorize) || winsorize) {
    scaled.values %<>% pmin(winsorize)
  }
  inter.es@x <- scaled.values
  dimnames(inter.es) <- list(cluster.names, cluster.names)
  return(inter.es)

}

#' Collapse Graph By Sum
#'
#' @param normalize boolean whether to recalculate edge weight as observed/expected (default=TRUE)
#' @inheritParams collapseGraphPaga
collapseGraphSum <- function(graph, groups, normalize=TRUE) {

  gcon <- contract.vertices(graph, groups, vertex.attr.comb=list('num'='sum',"ignore")) %>%
    simplify(edge.attr.comb=list(weight="sum","ignore"))

  if(normalize) {
    ex <- outer(V(gcon)$num,V(gcon)$num)/(sum(V(gcon)$num)*(sum(V(gcon)$num)-1)/2)*sum(E(graph)$weight)
    gcon2 <- graph_from_adjacency_matrix(as(as_adjacency_matrix(gcon, attr = "weight", sparse = FALSE)/ex,'dgCMatrix'), mode = "undirected", weighted=TRUE)
    V(gcon2)$num <- V(gcon)$num
    gcon <- gcon2;
  }

  if(is.factor(groups)) {
    V(gcon)$name <- levels(groups)
  } else {
    # not sure when this was actually needed
    gcon <- induced.subgraph(gcon, unique(groups))
  }

  return(gcon)
}

#' Get Cluster Graph
#'
#' @description Collapse vertices belonging to each cluster in a graph
#'
#' @param plot whether to show collapsed graph plot
#' @param method either "sum" or "paga" (default="sum")
#' @inheritParams collapseGraphPaga
#' @return collapsed graph
#' @export
getClusterGraph <- function(graph, groups, method="sum", plot=FALSE, node.scale=50, edge.scale=50, edge.alpha=0.3, ...) {
  
  V(graph)$num <- 1

  if (is.integer(groups) && is.null(names(groups))) {
    nv <- vcount(graph)
    if (length(groups)!=nv) {
      stop('Length of groups should be equal to the number of vertices')
    }
    if (max(groups)>nv) {
      stop('Groups specifies ids that are larger than the number of vertices in the graph')
    }
    if (any(is.na(groups))) {
      # remove vertices that are not part of the groups
      vi <- which(!is.na(groups))
      g <- induced.subgraph(graph,vi)
      groups <- groups[vi]
    } else {
      g <- graph
    }
  } else {
    gn <- V(graph)$name
    groups <- na.omit(groups[names(groups) %in% gn])
    if (length(groups)<2) {
      stop('Valid names of groups elements include too few cells')
    }
    if (length(groups)<length(gn)) {
      g <- induced.subgraph(graph,names(groups))
    } else {
      g <- graph;
    }
    if(is.factor(groups)) {
      groups <- groups[V(g)$name]
    } else {
      groups <- as.factor(setNames(as.character(groups[V(g)$name]),V(g)$name))
    }
  }

  if (method == "sum") {
    gcon <- collapseGraphSum(g, groups, ...)
  } else if (method == "paga") {
    gcon <- collapseGraphPaga(g, groups, ...)
  } else {
    stop("Unknown method: ", method)
  }

  if (plot) {
    set.seed(1)
    par(mar = rep(0.1, 4))
    plot.igraph(gcon, layout=layout_with_fr(gcon), vertex.size=V(gcon)$num/(sum(V(gcon)$num)/node.scale), 
        edge.width=E(gcon)$weight/sum(E(gcon)$weight/edge.scale), adjustcolor('black', alpha.f=edge.alpha))
  }

  return(invisible(gcon))

}
