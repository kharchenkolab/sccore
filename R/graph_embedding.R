#' splitVectorByNodes
#' 
#' @param vec input vector to be divided
#' @param nodes nodes used to divide the vector 'vec' via split()
#' @param n.nodes numeric The number of nodes for splitting
#' @return list from vec with names given by the nodes
#' @examples
#' adjList = graphToAdjList(conosGraph)
#' print(names(adjList))
#' ## [1] "idx" "probabilities" "names" 
#' length(adjList$names)
#' ## [1] 12000
#'
#' @export
splitVectorByNodes <- function(vec, nodes, n.nodes) {
  res <- lapply(1:n.nodes, function(x) list())
  splitted <- split(vec, nodes)
  res[as.integer(names(splitted))] <- splitted
  return(res)
}

#' Convert igraph graph into an adjacency list
#' 
#' @param graph input igraph object
#' @return adjacency list, defined by list(idx=adj.list, probabilities=probs, names=edge.list.fact$levels
#' @examples
#' library(dplyr)
#' edge.list.fact <- igraph::as_edgelist(conosGraph) %>% as_factor()
#' edge.list <- matrix(edge.list.fact$values, ncol=2)
#' n.nodes <- length(igraph::V(conosGraph))
#' splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes)
#' 
#' @export
graphToAdjList <- function(graph) {
  edge.list.fact <- igraph::as_edgelist(graph) %>% as_factor()
  edge.list <- matrix(edge.list.fact$values, ncol=2)
  n.nodes <- length(igraph::V(graph))
  adj.list <- mapply(c, splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes),
                     splitVectorByNodes(edge.list[,2], edge.list[,1], n.nodes)) %>%
    lapply(unlist) %>% lapply(`-`, 1)

  probs <- mapply(c, splitVectorByNodes(igraph::E(graph)$weight, edge.list[,2], n.nodes),
                  splitVectorByNodes(igraph::E(graph)$weight, edge.list[,1], n.nodes)) %>%
    lapply(unlist) %>%
    lapply(function(x) x / sum(x))

  if (any(sapply(probs, function(x) sum(is.na(x))))){
    stop("NAs in transition probabilities")
  }

  return(list(idx=adj.list, probabilities=probs, names=edge.list.fact$levels))
}

#' Embed a k-nearest neighbor (kNN) graph within a UMAP. Used within embedGraphUmap(). Please see McInnes et al <doi:10.21105/joss.00861> for the UMAP description and implementation.
#' 
#' @param commute.times graph commute times from get_nearest_neighbors(). The definition of commute_time(u, v) is the expected time starting at u = to reach v and then return to u . 
#' @param n.neighbors numeric Number of neighbors
#' @param names vector of names for UMAP rownames (default=NULL)
#' @param n.cores numeric Number of cores to use (except during stochastic gradient descent) (default=1). See 'n_threads' in uwot::umap()
#' @param n.epochs numeric Number of epochs to use during the optimization of the embedded coordinates (default=1000). See 'n_epochs' in uwot::umap()
#' @param spread numeric The effective scale of embedded points (numeric default=15). See 'spread' in uwot::umap()
#' @param min.dist numeric The effective minimum distance between embedded points (default=0.001). See 'min.dist' in uwot::umap()
#' @param n.sgd.cores  numeric Number of cores to use during stochastic gradient descent. If set to > 1, then results will not be reproducible, even if 'set.seed' is called with a fixed seed before running (default=n.cores) See 'n_sgd_threads' in uwot::umap()
#' @param target.dims numeric Dimensions for 'n_components' in uwot::umap(n_components=target.dims) (default=2)
#' @param verbose boolean Verbose output (default=TRUE)
#' @param ... arguments passed to uwot::umap()
#' @return resulting kNN graph embedding within a UMAP
#' 
#' @export
embedKnnGraph <- function(commute.times, n.neighbors, names=NULL, n.cores=1, n.epochs=1000, spread=15, min.dist=0.001, n.sgd.cores=n.cores, target.dims=2, verbose=TRUE, ...) {
  min.n.neighbors <- sapply(commute.times$idx, length) %>% min()
  if (min.n.neighbors < n.neighbors) {
    n.neighbors <- min.n.neighbors
    warning("Maximal number of estimated neighbors is ", min.n.neighbors, ". Consider increasing min.visited.verts, min.prob or min.prob.lower.")
  }

  ct.top <- sapply(commute.times$dist, `[`, 1:n.neighbors) %>% t() + 1
  ct.top.ids <- sapply(commute.times$idx, `[`, 1:n.neighbors) %>% t() + 1

  ct.top.ids <- cbind(1:nrow(ct.top.ids), ct.top.ids)
  ct.top <- cbind(rep(0, nrow(ct.top)), ct.top)

  umap <- uwot::umap(data.frame(x=rep(0, nrow(ct.top))), nn_method=list(idx=ct.top.ids, dist=ct.top), n_components=target.dims, n_threads=n.cores, n_epochs=n.epochs, spread=spread, min_dist=min.dist, n_sgd_threads=n.sgd.cores, verbose=verbose, ...)
  rownames(umap) <- names
  return(umap)
}

#' Embed a graph into a UMAP, Uniform Manifold Approximation and Projection for Dimension Reduction, <https://github.com/lmcinnes/umap>, <doi:10.21105/joss.00861>
#' 
#' @param graph input igraph object
#' @param min.prob numeric Minimum probability for proximity when calculating hitting time per neighbors (default=1e-3)
#' @param min.visited.verts numeric Minimum number of vertices visted when calculating hitting time per neighbors (default=1000)
#' @param n.cores numeric Number of cores to use (default=1)
#' @param max.hitting.nn.num numeric Maximum adjacencies for calculating hitting time per neighbor, hitting_time_per_neighbors() (default=0)
#' @param max.commute.nn.num numeric Maximum adjacencies for calculating commute time per neighbor, commute_time_per_node() (default=0)
#' @param min.prob.lower numeric Probability threshold to continue iteration in depth first search hitting time, dfs_hitting_time() (default=1e-7)
#' @param n.neighbors numeric Number of neighbors (default=40)
#' @param n.epochs numeric Number of epochs to use during the optimization of the embedded coordinates (default=1000). See 'n_epochs' in uwot::umap()
#' @param spread numeric The effective scale of embedded points (numeric default=15). See 'spread' in uwot::umap()
#' @param min.dist numeric The effective minimum distance between embedded points (default=0.001). See 'min.dist' in uwot::umap()
#' @param return.all boolean If TRUE, return list(adj.info=adj.info, commute.times=commute.times, umap=umap). Otherwise, just return UMAP(default=FALSE)
#' @param n.sgd.cores numeric Number of cores to use during stochastic gradient descent. If set to > 1, then results will not be reproducible, even if 'set.seed' is called with a fixed seed before running (default=n_threads) See 'n_sgd_threads' in uwot::umap()
#' @param verbose boolean Verbose output (default=TRUE)
#' @param ... Additional arguments passed to embedKnnGraph()
#' @return resulting UMAP embedding
#' 
#' @export
embedGraphUmap <- function(graph, min.prob=1e-3, min.visited.verts=1000, n.cores=1,
                           max.hitting.nn.num=0, max.commute.nn.num=0, min.prob.lower=1e-7,
                           n.neighbors=40, n.epochs=1000, spread=15, min.dist=0.001, return.all=FALSE,
                           n.sgd.cores=n.cores, verbose=TRUE, ...) {
  conn.comps <- igraph::components(graph)
  if (conn.comps$no > 1) {
    warning("Conos graph is not connected. Embedding may behave unexpectedly. ",
            "Please, consider increasing 'k' and/or 'k.self' parameters of 'buildGraph'")
  }
  min.visited.verts = min(min.visited.verts, min(conn.comps$csize) - 1)
  if (max.hitting.nn.num == 0) {
    max.hitting.nn.num <- length(igraph::V(graph)) - 1
  }

  if (verbose) message("Convert graph to adjacency list...")
  adj.info <- graphToAdjList(graph)
  if (verbose) message("Done")

  if (verbose) message("Estimate nearest neighbors and commute times...")
  commute.times <- get_nearest_neighbors(adj.info$idx, adj.info$probabilities, min_prob=min.prob,
                                         min_visited_verts=min.visited.verts, n_cores=n.cores, max_hitting_nn_num=max.hitting.nn.num,
                                         max_commute_nn_num=max.commute.nn.num, min_prob_lower=min.prob.lower, verbose=verbose)
  if (verbose) message("Done")

  if (verbose) message("Estimate UMAP embedding...")
  umap <- embedKnnGraph(commute.times, n.neighbors=n.neighbors, names=adj.info$names, n.cores=n.cores, n.epochs=n.epochs, spread=spread, min.dist=min.dist, verbose=verbose, n.sgd.cores=n.sgd.cores, ...)
  if (verbose) message("Done")

  if (return.all){
    return(list(adj.info=adj.info, commute.times=commute.times, umap=umap))
  }

  return(umap)
}