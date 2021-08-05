#' Add ranks for the fusions.
#'
#' Add ranks for the fusions to results reported by the specified tool. The
#' total number of fusions called by that tool will also be reported. The
#' resulting data frame will be sorted by rank. This is an extension of the
#' equalize function.
#'
#' @param tool The tool that produced the results.
#' @param results The results data frame.
#' @return A data frame with rank information columns added to the results.
#' @export
#' @examples
#' add_ranks('ericScript', results)
add_ranks <- function(tool, results) {
  if (tool == "ericScript") {
    results <- add_ranks_eric_script(results)
  } else if (tool == "fusionCatcher") {
    results <- add_ranks_fusion_catcher(results)
  } else if (tool == "fusionMap") {
    results <- add_ranks_fusion_map(results)
  } else if (tool == "jaffa") {
    results <- add_ranks_jaffa(results)
  } else if (tool == "mapSplice") {
    results <- add_ranks_map_splice(results)
  } else if (tool == "soapFuse") {
    results <- add_ranks_soap_fuse(results)
  } else if (tool == "starFusion") {
    results <- add_ranks_star_fusion(results)
  } else if (tool == "tophatFusion") {
    results <- add_ranks_tophat_fusion(results)
  } else if (tool == "dragen") {
    results <- add_ranks_dragen(results)
  } else if (tool == "arriba") {
    results <- add_ranks_arriba(results)
  } else if (tool == "cicero") {
    results <- add_ranks_cicero(results)
  } else {
    print(paste0("Tool ", tool, " is not a valid tool (add_ranks)."))
    stop()
  }
  return(results)
}

add_ranks_eric_script <- function(results) {
  results <- dplyr::arrange(results, dplyr::desc(EricScore))
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_fusion_catcher <- function(results) {
  results <- dplyr::arrange(results, dplyr::desc(SupportingReads))
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_fusion_map <- function(results) {
  results <- dplyr::arrange(results, dplyr::desc(SupportingReads))
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_jaffa <- function(results) {
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_map_splice <- function(results) {
  results <- dplyr::arrange(results, dplyr::desc(coverage))
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_soap_fuse <- function(results) {
  results <- dplyr::arrange(results, dplyr::desc(SupportingReads))
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_star_fusion <- function(results) {
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_tophat_fusion <- function(results) {
  results <- dplyr::arrange(results, dplyr::desc(Score))
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_dragen <- function(results) {
  results <- dplyr::arrange(results, dplyr::desc(Score))
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_arriba <- function(results) {
  results <- dplyr::arrange(results, dplyr::desc(SupportingReads))
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}

add_ranks_cicero <- function(results) {
  results <- dplyr::arrange(results, dplyr::desc(SupportingReads))
  results <- dplyr::mutate(results, Rank = seq.int(nrow(results)), Total = nrow(results))
  return(results)
}
