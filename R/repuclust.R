


#' perform greedy clustering into reporting units to visualize relationships
#'
#' @inheritParams self_assign
#' @export
repuclust <- function(reference, gen_start_col) {

  # first compute and get all the scaled likelihoods
  sa <- self_assign(reference = reference, gen_start_col = gen_start_col)


  # for each individual and each pair of collection-and-inferred-collection
  # we compute the "normalized scaled likelihood."  That is, the posterior for
  # the individual belonging to the correct versus the other collection.  Oh yeah,
  # we also initialize everyone's `repuclust` to just be their true collection.
  # note that we try to ensure that if the collection and repunit are already factors, we
  # preserve those numbers with them.
  pair_posts <- sa %>%
    dplyr::mutate(repuclust = as.integer(factor(collection))) %>%
    dplyr::mutate(inferred_repuclust = as.integer(factor(inferred_collection, levels = levels(factor(collection))))) %>%
    dplyr::group_by(indiv) %>%
    dplyr::mutate(self_post = scaled_likelihood[collection == inferred_collection]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pair_post = scaled_likelihood / (scaled_likelihood + self_post)) %>%
    dplyr::filter(collection != inferred_collection) %>%
    dplyr::ungroup()


  # now, having done that, we want to also get a tibble with information about the
  # tips on this tree
  leaves <- pair_posts %>%
    dplyr::count(collection, repuclust) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n) %>%
    dplyr::arrange(repuclust)


  # now, put the tips on the list that will store the nodes in the tree
  nodes <- lapply(leaves$repuclust, function(x) list(d1 = NA, d2 = NA, members = x))
  max_node <- length(nodes)  # also, get the number of the final leaf node (will get incremented, later)
  ntips <- nrow(leaves)  # just store the number of leaves that we are dealing with


  # now initialize some things that will get successively overwritten:
  # toss the columns that we don't need so we can get through things faster...
  pp <- pair_posts %>%
    dplyr::select(indiv, repuclust, inferred_repuclust, pair_post)

  # now, cycle over the interior nodes and figure out who shall be in them
  for (i in (ntips + 1):(2 * ntips - 1)) {
    rmp <- top_repclust_mean_pair(pp)       # find pair to merge
    nodes[[i]] <- new_node(rmp, nodes, i)      # make the new node that has them merged
    pp <- merge_repuclusts(pp = pp,
                           new_repu = i,
                           d1 = nodes[[i]]$d1,
                           d2 = nodes[[i]]$d2)
  }


  # now, nodes is a list that has the info that we want
  # we will return that, along with some other things
  list(ntips = ntips,
       leaves = leaves,
       nodes = nodes)
}



#' find the pair that includes the highest incorrect normalized scaled likelihood
#'
#' @param pp a data frame like pair_posts
#' @keywords internal
top_repclust_mean_pair <- function(pp) {
  pp %>%
    dplyr::group_by(repuclust, inferred_repuclust) %>%
    dplyr::summarise(mean_pairwise = mean(pair_post)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(inferred_repuclust != repuclust) %>%
    dplyr::arrange(desc(mean_pairwise)) %>%
    dplyr::filter(repuclust %in% c(repuclust[1], inferred_repuclust[1]),
           inferred_repuclust %in% c(repuclust[1], inferred_repuclust[1])
    )
}


#' given the top reclust mean pair, make a new node out of them
#' @param rmp  the output of top_repclust_mean_pair
#' @param nodes a list of nodes that are already there (so the members can be obtained)
#' @param idx the index that this node should take
#' @keywords internal
new_node <- function(rmp, nodes, idx) {
  n1 <- rmp$repuclust[1]  # store the nodes of these
  n2 <- rmp$repuclust[2]
  ret <- list()
  ret$idx <- idx
  ret$d1 <- n1
  ret$p1 <- rmp$mean_pairwise[1]
  ret$d2 <- n2
  ret$p2 <- rmp$mean_pairwise[2]
  ret$members <- sort(base::union(nodes[[n1]]$members, nodes[[n2]]$members))
  ret
}



#' simple function to merge the repuclusts in the pair_posts data frame
#' @param pp the pair_posts data frame
#' @param new_repu the index of the new reporting unit node
#' @param d1 the repuclust index of the first daughter of new_repu
#' @param d2 the repuclust index of the second daughter of new_repu
#' @keywords internal
merge_repuclusts <- function(pp, new_repu, d1, d2) {
  pp %>%
    dplyr::mutate(repuclust = ifelse(repuclust %in% c(d1, d2), new_repu, repuclust),
           inferred_repuclust = ifelse(inferred_repuclust %in% c(d1, d2), new_repu, inferred_repuclust))
}



#' convert the output of repuclust to an ape phylo object
#' @param r the output from repuclust()
repuclust2phylo <- function(r) {
  edge_df <- lapply(r$nodes[-(1:boing$ntips)], function(x) {
    tibble::tibble(parent = x$idx,
                   daughter = c(x$d1, x$d2),
                   probs = c(x$p1, x$p2))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(edge_lengths = 2 * (0.5 - probs)) %>%
    dplyr::mutate(edge_lengths = ifelse(edge_lengths < 0.9, .005, edge_lengths))

  ret <- list()
  ret$edge <- cbind(edge_df$parent, edge_df$daughter)
  ret$Nnode <- length(r$nodes) - r$ntips
  ret$tip.label <- as.character(r$leaves$collection)
  ret$edge.length <- edge_df$edge_lengths

  class(ret) <- "phylo"

}
