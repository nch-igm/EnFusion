#' Get the overlapping fusions.
#'
#' Get an array of overlapping fusions that appear in at least two of the
#' supplied fusion arrays.
#'
#' @param fusions A list of character arrays with fusions.
#' @export
#' @return A character array with the overlapping fusions. The character vector
#'   will be empty if there is no overlap, character(0).
#' @examples
#' get_overlapping_from_fusions(fusions)
get_overlapping_from_fusions <- function(fusions) {
  overlapping <- character(0)
  for (i in seq_along(fusions)) {
    for (j in i:length(fusions)) {
      if (i != j) {
        overlapping <- c(overlapping, intersect(fusions[[i]], fusions[[j]]))
      }
    }
  }
  return(unique(overlapping))
}

#' Get the report from the data.
#'
#' Get all of the information from the data that will be included in the report,
#' given the list of overlapping fusions.
#'
#' @param overlapping A character array of overlapping fusions.
#' @param all_data The list that contains all of the info for all the results.
#' @return A data frame with information to write to the report.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' get_report(overlapping, all_data)
get_report <- function(overlapping, all_data) {
  report <- tibble::tibble(UnorderedFusion = character(), OrderedFusion = character(), 
    Tool = character(), Rank = integer(), Total = integer(), SupportingReads = integer(),
    HGNCSymbol1 = character(), HGNCSymbol2 = character(), multimatch1 = character(),
    multimatch2 = character(), alias_match_1 = character(), alias_match_2 = character(), 
    Gene1 = character(), Chr1 = character(), Break1 = character(), Strand1 = character(), 
    Gene2 = character(), Chr2 = character(), Break2 = character(), Strand2 = character(), 
    FMFusionID <- character(), FrameShiftClass <- character(), FMKnownTranscript1 <- character(), 
    FMKnownExonNumber1 <- character(), FMKnownTranscript2 <- character(), FMKnownExonNumber2 <- character(), 
    FusionJunctionSequence <- character())
  for (i in 1:length(all_data)) {
    to_add <- all_data[[i]]$df %>% dplyr::filter(UnorderedFusion %in% overlapping) %>% 
      dplyr::group_by(UnorderedFusion) %>% dplyr::filter(Rank == min(Rank)) %>% 
      dplyr::mutate(Tool = all_data[[i]]$tool, Break1 = as.character(Break1), 
        Break2 = as.character(Break2)) %>% dplyr::select(UnorderedFusion, 
      OrderedFusion, Tool, Rank, Total, SupportingReads, HGNCSymbol1, Gene1, multimatch1, alias_match_1, Chr1, Break1, 
      Strand1, HGNCSymbol2, Gene2, multimatch2, alias_match_2, Chr2, Break2, Strand2, FMFusionID, FrameShiftClass, FMKnownTranscript1, 
      FMKnownExonNumber1, FMKnownTranscript2, FMKnownExonNumber2, FusionJunctionSequence)
    report <- dplyr::bind_rows(report, to_add)
  }
  report <- report %>% dplyr::group_by(UnorderedFusion) %>% dplyr::mutate(NumTools = length(unique(Tool)))
  return(report)
}

#' Get the report for all breakpoints from the data.
#'
#' Get all of the information from the data that will be included in the report,
#' given the list of overlapping fusions, for all breakpoints between the gene
#' pairs.
#'
#' @param overlapping A character array of overlapping fusions.
#' @param all_data The list that contains all of the info for all the results.
#' @return A data frame with information to write to the report for all
#'   breakpoints.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' get_report_all_breakpoints(overlapping, all_data)
get_report_all_breakpoints <- function(overlapping, all_data) {
  report <- tibble::tibble(UnorderedFusion = character(), OrderedFusion = character(), 
    Tool = character(), Rank = integer(), Total = integer(), SupportingReads = integer(), 
    Gene1 = character(), Chr1 = character(), Break1 = character(), Strand1 = character(), 
    Gene2 = character(), Chr2 = character(), Break2 = character(), Strand2 = character(), 
    HGNCSymbol1 = character(), HGNCSymbol2 = character(), multimatch1 = character(),
    multimatch2 = character(), alias_match_1 = character(), alias_match_2 = character(), 
    FMFusionID <- character(), FrameShiftClass <- character(), FMKnownTranscript1 <- character(), 
    FMKnownExonNumber1 <- character(), FMKnownTranscript2 <- character(), FMKnownExonNumber2 <- character(), 
    FusionJunctionSequence <- character())
  for (i in 1:length(all_data)) {
    to_add <- all_data[[i]]$df %>% dplyr::filter(UnorderedFusion %in% overlapping) %>% 
      dplyr::mutate(Tool = all_data[[i]]$tool, Break1 = as.character(Break1), 
        Break2 = as.character(Break2)) %>% dplyr::select(UnorderedFusion, 
          OrderedFusion, Tool, Rank, Total, SupportingReads, HGNCSymbol1, Gene1, multimatch1, alias_match_1, Chr1, Break1, 
          Strand1, HGNCSymbol2, Gene2, multimatch2, alias_match_2, Chr2, Break2, Strand2, FMFusionID, FrameShiftClass, FMKnownTranscript1, 
          FMKnownExonNumber1, FMKnownTranscript2, FMKnownExonNumber2, FusionJunctionSequence)
    report <- dplyr::bind_rows(report, to_add)
  }
  report <- report %>% dplyr::group_by(UnorderedFusion) %>% dplyr::mutate(NumTools = length(unique(Tool)))
  return(report)
}

#' Order the report.
#'
#' Order the report so that it's in a format ready to be written to output.
#'
#' @param report A data frame with the report information.
#' @return A data frame with the ordered information to write to output.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' order_report(report)
order_report <- function(report) {
  report <- report %>% dplyr::group_by(UnorderedFusion) %>% dplyr::mutate(AveRank = mean(Rank)) %>% 
    dplyr::ungroup() %>% dplyr::arrange(dplyr::desc(NumTools), AveRank, UnorderedFusion, 
    Rank, Total, dplyr::desc(SupportingReads), Tool) %>% dplyr::select(UnorderedFusion, 
          OrderedFusion, NumTools, Tool, Rank, Total, SupportingReads, HGNCSymbol1, Gene1, multimatch1, alias_match_1, Chr1, Break1,                                                             
          Strand1, HGNCSymbol2, Gene2, multimatch2, alias_match_2, Chr2, Break2, Strand2, FMFusionID, FrameShiftClass, FMKnownTranscript1,
          FMKnownExonNumber1, FMKnownTranscript2, FMKnownExonNumber2, FusionJunctionSequence)
  return(report)
}

#' Order the report for all breakpoints.
#'
#' Order the report for all breakpoints so that it's in a format ready to be
#' written to output. The ordering is taken from the original report.
#'
#' @param report A data frame with the report information.
#' @return A data frame with the ordered information to write to output.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' order_report_all_breakpoints(report)
order_report_all_breakpoints <- function(report, original) {
  if (nrow(original) > 0) {
    FusionCall = seq.int(nrow(original))
  } else {
    FusionCall = integer(0)
  }
  original_ordering <- original %>% dplyr::mutate(FusionCall = FusionCall) %>% 
    dplyr::select(UnorderedFusion, Tool, FusionCall)
  report <- report %>% dplyr::left_join(original_ordering, by = c("UnorderedFusion", 
    "Tool")) %>% dplyr::arrange(FusionCall, Rank) %>% dplyr::select(UnorderedFusion, 
    OrderedFusion, NumTools, Tool, Rank, Total, SupportingReads, HGNCSymbol1, Gene1, multimatch1, alias_match_1, Chr1, Break1,  
    Strand1, HGNCSymbol2, Gene2, multimatch2, alias_match_2, Chr2, Break2, Strand2, FMFusionID, FrameShiftClass, FMKnownTranscript1,
    FMKnownExonNumber1, FMKnownTranscript2, FMKnownExonNumber2, FusionJunctionSequence)
  return(report)
}

#' @param fusions A list of character arrays with fusions.
#' @export
#' @return A character array with the union fusions. The character vector
#'   will be empty if there are no fusions to union, character(0).
#' @examples
#' get_union_from_fusions(fusions)
get_union_from_fusions <- function(fusions) {
  fusion_union <- character(0)
  for (i in seq_along(fusions)) {
    for (j in i:length(fusions)) {
      if (i != j) {
        fusion_union <- c(fusion_union, union(fusions[[i]], fusions[[j]]))
      }
    }
  }
  return(unique(fusion_union))
}

#' Get the report from union data.
#'
#' Get all of the information from the data that will be included in the report,
#' given the list of overlapping fusions.
#'
#' @param uniond A character array of uniond fusions.
#' @param all_data The list that contains all of the info for all the results.
#' @return A data frame with information to write to the report.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' get_report(overlapping, all_data)
get_union_report <- function(fusion_union, all_data) {
  union_report <- tibble::tibble(UnorderedFusion = character(), OrderedFusion = character(), 
                                 Tool = character(), Rank = integer(), Total = integer(), SupportingReads = integer(),
                                 HGNCSymbol1 = character(), HGNCSymbol2 = character(), multimatch1 = character(),
                                 multimatch2 = character(), alias_match_1 = character(), alias_match_2 = character(), 
                                 Gene1 = character(), Chr1 = character(), Break1 = character(), Strand1 = character(), 
                                 Gene2 = character(), Chr2 = character(), Break2 = character(), Strand2 = character(), 
                                 FMFusionID <- character(), FrameShiftClass <- character(), FMKnownTranscript1 <- character(), 
                                 FMKnownExonNumber1 <- character(), FMKnownTranscript2 <- character(), FMKnownExonNumber2 <- character(), 
                                 FusionJunctionSequence <- character())
  for (i in 1:length(all_data)) {
    to_add <- all_data[[i]]$df %>% dplyr::filter(UnorderedFusion %in% fusion_union) %>% 
      dplyr::group_by(UnorderedFusion) %>% dplyr::filter(Rank == min(Rank)) %>% 
      dplyr::mutate(Tool = all_data[[i]]$tool, Break1 = as.character(Break1), 
        Break2 = as.character(Break2)) %>% dplyr::select(UnorderedFusion, 
                                                         OrderedFusion, Tool, Rank, Total, SupportingReads, HGNCSymbol1, Gene1, multimatch1, alias_match_1, Chr1, Break1, 
                                                         Strand1, HGNCSymbol2, Gene2, multimatch2, alias_match_2, Chr2, Break2, Strand2, FMFusionID, FrameShiftClass, FMKnownTranscript1, 
                                                         FMKnownExonNumber1, FMKnownTranscript2, FMKnownExonNumber2, FusionJunctionSequence)
    report <- dplyr::bind_rows(report, to_add)
    report <- dplyr::bind_rows(report, to_add)
    union_report <- dplyr::bind_rows(union_report, to_add)
  }
  union_report <- union_report %>% dplyr::group_by(UnorderedFusion) %>% dplyr::mutate(NumTools = length(unique(Tool)))
  return(union_report)
}

#' Order the union_report.
#'
#' Order the union_report so that it's in a format ready to be written to output.
#'
#' @param union_report A data frame with the union_report information.
#' @return A data frame with the ordered information to write to output.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' order_union_report(union_report)
order_union_report <- function(union_report) {
  union_report <- union_report %>% dplyr::group_by(UnorderedFusion) %>% dplyr::mutate(AveRank = mean(Rank)) %>% 
    dplyr::ungroup() %>% dplyr::arrange(dplyr::desc(NumTools), AveRank, UnorderedFusion, 
    Rank, Total, dplyr::desc(SupportingReads), Tool) %>% dplyr::select(UnorderedFusion,
    OrderedFusion, NumTools, Tool, Rank, Total, SupportingReads, HGNCSymbol1, Gene1, multimatch1, alias_match_1,
    Chr1, Break1, Strand1, HGNCSymbol2, Gene2, multimatch2, alias_match_2, Chr2, Break2, Strand2)
  return(union_report)
}
