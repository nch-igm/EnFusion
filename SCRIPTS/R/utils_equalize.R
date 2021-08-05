#' Equalize the results.
#'
#' Equalize the results so the tools will have a set of common columns that can
#' be used in downstream analysis.
#'
#' The fusion detection tools write results that have different formats, but
#' there is a set of information that is common between the tools. This includes
#' the fusion genes and breakpoints that were called. In the equalize functions,
#' this common set of information is converted to a standard format so that the
#' tools can be treated equally in downstream analysis.
#'
#' @param tool The tool that produced the results.
#' @param results The results data frame.
#' @return A data frame with additional columns added to the results.
#' @export
#' @examples
#' equalize_results("ericScript", eric_script)
#' equalize_results("fusionCatcher", fusion_catcher)
#' equalize_results("fusionMap", fusion_map)
#' equalize_results("jaffa", jaffa)
#' equalize_results("mapSplice", map_splice)
#' equalize_results("soapFuse", soap_fuse)
#' equalize_results("starFusion", star_fusion)
#' equalize_results("tophatFusion", tophat_fusion)
equalize_results <- function(tool, results) {
  if (tool == "ericScript") {
    results <- equalize_eric_script(results)
  } else if (tool == "fusionCatcher") {
    results <- equalize_fusion_catcher(results)
  } else if (tool == "fusionMap") {
    results <- equalize_fusion_map(results)
  } else if (tool == "jaffa") {
    results <- equalize_jaffa(results)
  } else if (tool == "mapSplice") {
    results <- equalize_map_splice(results)
  } else if (tool == "soapFuse") {
    results <- equalize_soap_fuse(results)
  } else if (tool == "starFusion") {
    results <- equalize_star_fusion(results)
  } else if (tool == "tophatFusion") {
    results <- equalize_tophat_fusion(results)
  } else if (tool == "dragen") {
    results <- equalize_dragen(results)
  } else if (tool == "arriba") {
    results <- equalize_arriba(results)
  } else if (tool == "cicero") {
    results <- equalize_cicero(results)
  } else {
    print(paste0("Tool ", tool, " is not a valid tool (equalize_results)."))
    stop()
  }
  return(results)
}


equalize_fusion_catcher <- function(fusion_catcher) {
  fusion_catcher$Gene1 <- fusion_catcher$`Gene_1_symbol(5end_fusion_partner)`
  fusion_catcher$Gene2 <- fusion_catcher$`Gene_2_symbol(3end_fusion_partner)`
  
  fusion_catcher$SupportingReads <- fusion_catcher$Spanning_pairs
  # Spanning_pairs: Count of pairs of reads supporting the fusion (including
  # also the multimapping reads). Spanning_unique_reads: Count of unique reads
  # (i.e. unique mapping positions) mapping on the fusion junctions. Shortly,
  # here are counted all the reads which map on the fusion junction minus the
  # PCR duplicate reads.

  fusion_catcher$Chr1 <- vapply(strsplit(fusion_catcher$`Fusion_point_for_gene_1(5end_fusion_partner)`, ":",
                                         fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  fusion_catcher$Chr1 <- equalize_chrom(fusion_catcher$Chr1)
  fusion_catcher$Break1 <- vapply(strsplit(fusion_catcher$`Fusion_point_for_gene_1(5end_fusion_partner)`, ":",
                                           fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  fusion_catcher$Strand1 <- vapply(strsplit(fusion_catcher$`Fusion_point_for_gene_1(5end_fusion_partner)`, ":",
                                            fixed = TRUE), `[`, 3, FUN.VALUE = character(1))
  fusion_catcher$Chr2 <- vapply(strsplit(fusion_catcher$`Fusion_point_for_gene_2(3end_fusion_partner)`, ":",
                                         fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  fusion_catcher$Chr2 <- equalize_chrom(fusion_catcher$Chr2)
  fusion_catcher$Break2 <- vapply(strsplit(fusion_catcher$`Fusion_point_for_gene_2(3end_fusion_partner)`, ":",
                                           fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  fusion_catcher$Strand2 <- vapply(strsplit(fusion_catcher$`Fusion_point_for_gene_2(3end_fusion_partner)`, ":",
                                            fixed = TRUE), `[`, 3, FUN.VALUE = character(1))
  
  alias_list <- readRDS("/SCRIPTS/R/alias_list_071221.rds")
  fusion_catcher$aliasbreak1 <- paste0("chr",fusion_catcher$Chr1, ":", fusion_catcher$Break1)
  fusion_catcher$aliasbreak2 <- paste0("chr",fusion_catcher$Chr2, ":", fusion_catcher$Break2)
  
  fusion_catcher$HGNCSymbol1 <- unlist(lapply(fusion_catcher$Gene1, alias_matching, breakpoint = fusion_catcher$aliasbreak1, alias_list = alias_list))
  fusion_catcher$HGNCSymbol2 <- unlist(lapply(fusion_catcher$Gene2, alias_matching, breakpoint = fusion_catcher$aliasbreak2, alias_list = alias_list))
  
  fusion_catcher$multimatch1 <- ifelse(grepl('\\|', fusion_catcher$HGNCSymbol1), "Yes", "No")
  fusion_catcher$multimatch2 <- ifelse(grepl('\\|', fusion_catcher$HGNCSymbol2), "Yes", "No")
  fusion_catcher$alias_match_1 <- ifelse(fusion_catcher$Gene1 == fusion_catcher$HGNCSymbol1, "Yes", "No")
  fusion_catcher$alias_match_2 <- ifelse(fusion_catcher$Gene2 == fusion_catcher$HGNCSymbol2, "Yes", "No")
  
  fusion_catcher$OrderedFusion <- get_ordered(fusion_catcher$HGNCSymbol1, fusion_catcher$HGNCSymbol2)
  fusion_catcher$UnorderedFusion <- get_unordered(fusion_catcher$HGNCSymbol1, fusion_catcher$HGNCSymbol2)

  fusion_catcher$FMFusionID <- character(nrow(fusion_catcher))
  fusion_catcher$FrameShiftClass <- character(nrow(fusion_catcher))
  fusion_catcher$FMKnownTranscript1 <- character(nrow(fusion_catcher))
  fusion_catcher$FMKnownExonNumber1 <- character(nrow(fusion_catcher))
  fusion_catcher$FMKnownTranscript2 <- character(nrow(fusion_catcher))
  fusion_catcher$FMKnownExonNumber2 <- character(nrow(fusion_catcher))
  fusion_catcher$FusionJunctionSequence <- character(nrow(fusion_catcher))

  return(fusion_catcher)
}

equalize_fusion_map <- function(fusion_map) {
  fusion_map$Gene1 <- fusion_map$KnownGene1
  fusion_map$Gene2 <- fusion_map$KnownGene2
  
  fusion_map$SupportingReads <- unname(unlist(fusion_map[, 3]))
  # SeedCount: Specifies the count of reads, for a particular FusionID, that
  # were considered a Seed Read. (A Seed Read is a read that has at least the
  # MinimalFusionAlignmentLength for each gene. >= a nt mapped in both sources.)
  # RescuedCount: Specifies the count of reads, for a particular FusionID, that
  # were considered a Rescued Read. (A Rescued Read is a read that does not have
  # at least the MinimalFusionAlignmentLength for each gene. < a nt mapped in
  # either source.)

  fusion_map$Chr1 <- equalize_chrom(fusion_map$Chromosome1)
  fusion_map$Break1 <- fusion_map$Position1
  fusion_map$Strand1 <- substr(fusion_map$Strand, 1, 1)
  fusion_map$Chr2 <- equalize_chrom(fusion_map$Chromosome2)
  fusion_map$Break2 <- fusion_map$Position2
  fusion_map$Strand2 <- substr(fusion_map$Strand, 2, 2)
  
  alias_list <- readRDS("/SCRIPTS/R/alias_list_071221.rds")
  fusion_map$aliasbreak1 <- paste0("chr",fusion_map$Chr1, ":", fusion_map$Break1)
  fusion_map$aliasbreak2 <- paste0("chr",fusion_map$Chr2, ":", fusion_map$Break2)
  
  fusion_map$HGNCSymbol1 <- unlist(lapply(fusion_map$Gene1, alias_matching, breakpoint = fusion_map$aliasbreak1, alias_list = alias_list))
  fusion_map$HGNCSymbol2 <- unlist(lapply(fusion_map$Gene2, alias_matching, breakpoint = fusion_map$aliasbreak2, alias_list = alias_list))
  
  fusion_map$multimatch1 <- ifelse(grepl('\\|', fusion_map$HGNCSymbol1), "Yes", "No")
  fusion_map$multimatch2 <- ifelse(grepl('\\|', fusion_map$HGNCSymbol2), "Yes", "No")
  fusion_map$alias_match_1 <- ifelse(fusion_map$Gene1 == fusion_map$HGNCSymbol1, "Yes", "No")
  fusion_map$alias_match_2 <- ifelse(fusion_map$Gene2 == fusion_map$HGNCSymbol2, "Yes", "No")
  
  fusion_map$OrderedFusion <- get_ordered(fusion_map$HGNCSymbol1, fusion_map$HGNCSymbol2)
  fusion_map$UnorderedFusion <- get_unordered(fusion_map$HGNCSymbol1, fusion_map$HGNCSymbol2)

  fusion_map$FMFusionID <- fusion_map$FusionID
  fusion_map$FrameShiftClass <- fusion_map$FrameShiftClass
  fusion_map$FMKnownTranscript1 <- fusion_map$KnownTranscript1
  fusion_map$FMKnownExonNumber1 <- as.character(fusion_map$KnownExonNumber1)
  fusion_map$FMKnownTranscript2 <- fusion_map$KnownTranscript2
  fusion_map$FMKnownExonNumber2 <- as.character(fusion_map$KnownExonNumber2)
  fusion_map$FusionJunctionSequence <- fusion_map$FusionJunctionSequence

  return(fusion_map)
}

equalize_jaffa <- function(jaffa) {
  jaffa$Gene1 <- vapply(strsplit(jaffa$`fusion genes`, ":", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  jaffa$Gene2 <- vapply(strsplit(jaffa$`fusion genes`, ":", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  
  jaffa$SupportingReads <- jaffa$`spanning pairs` + jaffa$`spanning reads`
  # spanning pairs: The number of read-pairs, where each read in the pair aligns
  # entirely on either side of the breakpoint. You might see a '-' in some of
  # these. This indicates that no spanning pairs were found, but that the contig
  # had only a small amount of flanking sequence to align reads to. i.e. the
  # spanning pairs results may not be indicative of the true support for the
  # fusion event.  spanning reads: The number of reads aligning to the
  # breakpoint, with at least 15 bases of flanking sequence either side (by
  # default).

  jaffa$Chr1 <- equalize_chrom(jaffa$chrom1)
  jaffa$Break1 <- jaffa$base1
  jaffa$Strand1 <- jaffa$strand1
  jaffa$Chr2 <- equalize_chrom(jaffa$chrom2)
  jaffa$Break2 <- jaffa$base2
  jaffa$Strand2 <- jaffa$strand2

  alias_list <- readRDS("/SCRIPTS/R/alias_list_071221.rds")
  jaffa$aliasbreak1 <- paste0("chr",jaffa$Chr1, ":", jaffa$Break1)
  jaffa$aliasbreak2 <- paste0("chr",jaffa$Chr2, ":", jaffa$Break2)
  
  jaffa$HGNCSymbol1 <- unlist(lapply(jaffa$Gene1, alias_matching, breakpoint = jaffa$aliasbreak1, alias_list = alias_list))
  jaffa$HGNCSymbol2 <- unlist(lapply(jaffa$Gene2, alias_matching, breakpoint = jaffa$aliasbreak2, alias_list = alias_list))

  jaffa$multimatch1 <- ifelse(grepl('\\|', jaffa$HGNCSymbol1), "Yes", "No")
  jaffa$multimatch2 <- ifelse(grepl('\\|', jaffa$HGNCSymbol2), "Yes", "No")
  jaffa$alias_match_1 <- ifelse(jaffa$Gene1 == jaffa$HGNCSymbol1, "Yes", "No")
  jaffa$alias_match_2 <- ifelse(jaffa$Gene2 == jaffa$HGNCSymbol2, "Yes", "No")
  
  jaffa$OrderedFusion <- get_ordered(jaffa$HGNCSymbol1, jaffa$HGNCSymbol2)
  jaffa$UnorderedFusion <- get_unordered(jaffa$HGNCSymbol1, jaffa$HGNCSymbol2)
  
  jaffa$FMFusionID <- character(nrow(jaffa))
  jaffa$FrameShiftClass <- character(nrow(jaffa))
  jaffa$FMKnownTranscript1 <- character(nrow(jaffa))
  jaffa$FMKnownExonNumber1 <- character(nrow(jaffa))
  jaffa$FMKnownTranscript2 <- character(nrow(jaffa))
  jaffa$FMKnownExonNumber2 <- character(nrow(jaffa))
  jaffa$FusionJunctionSequence <- character(nrow(jaffa))

  return(jaffa)
}

equalize_map_splice <- function(map_splice) {
  map_splice$Gene1 <- gsub(",$", "", map_splice$annotated_gene_donor)
  map_splice$Gene2 <- gsub(",$", "", map_splice$annotated_gene_acceptor)
  
  map_splice$SupportingReads <- map_splice$coverage

  map_splice$Chr1 <- vapply(strsplit(map_splice$chrom, "~", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  map_splice$Chr1 <- equalize_chrom(map_splice$Chr1)
  map_splice$Break1 <- map_splice$doner_end
  map_splice$Strand1 <- substr(map_splice$strand, 1, 1)
  map_splice$Chr2 <- vapply(strsplit(map_splice$chrom, "~", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  map_splice$Chr2 <- equalize_chrom(map_splice$Chr2)
  map_splice$Break2 <- map_splice$acceptor_start
  map_splice$Strand2 <- substr(map_splice$strand, 2, 2)
  
  alias_list <- readRDS("/SCRIPTS/R/alias_list_071221.rds")
  map_splice$aliasbreak1 <- paste0("chr",map_splice$Chr1, ":", map_splice$Break1)
  map_splice$aliasbreak2 <- paste0("chr",map_splice$Chr2, ":", map_splice$Break2)
  
  map_splice$HGNCSymbol1 <- unlist(lapply(map_splice$Gene1, alias_matching, breakpoint = map_splice$aliasbreak1, alias_list = alias_list))
  map_splice$HGNCSymbol2 <- unlist(lapply(map_splice$Gene2, alias_matching, breakpoint = map_splice$aliasbreak2, alias_list = alias_list))
  
  map_splice$multimatch1 <- ifelse(grepl('\\|', map_splice$HGNCSymbol1), "Yes", "No")
  map_splice$multimatch2 <- ifelse(grepl('\\|', map_splice$HGNCSymbol2), "Yes", "No")
  map_splice$alias_match_1 <- ifelse(map_splice$Gene1 == map_splice$HGNCSymbol1, "Yes", "No")
  map_splice$alias_match_2 <- ifelse(map_splice$Gene2 == map_splice$HGNCSymbol2, "Yes", "No")
  
  map_splice$OrderedFusion <- get_ordered(map_splice$HGNCSymbol1, map_splice$HGNCSymbol2)
  map_splice$UnorderedFusion <- get_unordered(map_splice$HGNCSymbol1, map_splice$HGNCSymbol2)

  map_splice$FMFusionID <- character(nrow(map_splice))
  map_splice$FrameShiftClass <- character(nrow(map_splice))
  map_splice$FMKnownTranscript1 <- character(nrow(map_splice))
  map_splice$FMKnownExonNumber1 <- character(nrow(map_splice))
  map_splice$FMKnownTranscript2 <- character(nrow(map_splice))
  map_splice$FMKnownExonNumber2 <- character(nrow(map_splice))
  map_splice$FusionJunctionSequence <- character(nrow(map_splice))

  return(map_splice)
}

equalize_soap_fuse <- function(soap_fuse) {
  soap_fuse$Gene1 <- soap_fuse$up_gene
  soap_fuse$Gene2 <- soap_fuse$dw_gene
  soap_fuse$OrderedFusion <- get_ordered(soap_fuse$Gene1, soap_fuse$Gene2)
  soap_fuse$UnorderedFusion <- get_unordered(soap_fuse$Gene1, soap_fuse$Gene2)
  soap_fuse$SupportingReads <- soap_fuse$Span_reads_num + soap_fuse$Junc_reads_num

  soap_fuse$Chr1 <- equalize_chrom(soap_fuse$up_chr)
  soap_fuse$Break1 <- soap_fuse$up_Genome_pos
  soap_fuse$Strand1 <- soap_fuse$up_strand
  soap_fuse$Chr2 <- equalize_chrom(soap_fuse$dw_chr)
  soap_fuse$Break2 <- soap_fuse$dw_Genome_pos
  soap_fuse$Strand2 <- soap_fuse$dw_strand

  soap_fuse$FMFusionID <- character(nrow(soap_fuse))
  soap_fuse$FrameShiftClass <- character(nrow(soap_fuse))
  soap_fuse$FMKnownTranscript1 <- character(nrow(soap_fuse))
  soap_fuse$FMKnownExonNumber1 <- character(nrow(soap_fuse))
  soap_fuse$FMKnownTranscript2 <- character(nrow(soap_fuse))
  soap_fuse$FMKnownExonNumber2 <- character(nrow(soap_fuse))
  soap_fuse$FusionJunctionSequence <- character(nrow(soap_fuse))

  return(soap_fuse)
}

equalize_star_fusion <- function(star_fusion) {
  star_fusion$Gene1 <- vapply(strsplit(star_fusion$FusionName, "--", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  star_fusion$Gene2 <- vapply(strsplit(star_fusion$FusionName, "--", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  star_fusion$SupportingReads <- star_fusion$JunctionReadCount
  
  star_fusion$Chr1 <- vapply(strsplit(star_fusion$LeftBreakpoint, ":", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  star_fusion$Chr1 <- equalize_chrom(star_fusion$Chr1)
  star_fusion$Break1 <- vapply(strsplit(star_fusion$LeftBreakpoint, ":", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  star_fusion$Strand1 <- vapply(strsplit(star_fusion$LeftBreakpoint, ":", fixed = TRUE), `[`, 3, FUN.VALUE = character(1))
  star_fusion$Chr2 <- vapply(strsplit(star_fusion$RightBreakpoint, ":", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  star_fusion$Chr2 <- equalize_chrom(star_fusion$Chr2)
  star_fusion$Break2 <- vapply(strsplit(star_fusion$RightBreakpoint, ":", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  star_fusion$Strand2 <- vapply(strsplit(star_fusion$RightBreakpoint, ":", fixed = TRUE), `[`, 3, FUN.VALUE = character(1))

  alias_list <- readRDS("/SCRIPTS/R/alias_list_071221.rds")
  star_fusion$aliasbreak1 <- paste0("chr",star_fusion$Chr1, ":", star_fusion$Break1)
  star_fusion$aliasbreak2 <- paste0("chr",star_fusion$Chr2, ":", star_fusion$Break2)
  
  star_fusion$HGNCSymbol1 <- unlist(lapply(star_fusion$Gene1, alias_matching, breakpoint = star_fusion$aliasbreak1, alias_list = alias_list))
  star_fusion$HGNCSymbol2 <- unlist(lapply(star_fusion$Gene2, alias_matching, breakpoint = star_fusion$aliasbreak2, alias_list = alias_list))
  
  star_fusion$multimatch1 <- ifelse(grepl('\\|', star_fusion$HGNCSymbol1), "Yes", "No")
  star_fusion$multimatch2 <- ifelse(grepl('\\|', star_fusion$HGNCSymbol2), "Yes", "No")
  star_fusion$alias_match_1 <- ifelse(star_fusion$Gene1 == star_fusion$HGNCSymbol1, "Yes", "No")
  star_fusion$alias_match_2 <- ifelse(star_fusion$Gene2 == star_fusion$HGNCSymbol2, "Yes", "No")
  
  star_fusion$OrderedFusion <- get_ordered(star_fusion$HGNCSymbol1, star_fusion$HGNCSymbol2)
  star_fusion$UnorderedFusion <- get_unordered(star_fusion$HGNCSymbol1, star_fusion$HGNCSymbol2)
  
  star_fusion$FMFusionID <- character(nrow(star_fusion))
  star_fusion$FrameShiftClass <- character(nrow(star_fusion))
  star_fusion$FMKnownTranscript1 <- character(nrow(star_fusion))
  star_fusion$FMKnownExonNumber1 <- character(nrow(star_fusion))
  star_fusion$FMKnownTranscript2 <- character(nrow(star_fusion))
  star_fusion$FMKnownExonNumber2 <- character(nrow(star_fusion))
  star_fusion$FusionJunctionSequence <- character(nrow(star_fusion))

  return(star_fusion)
}

equalize_tophat_fusion <- function(tophat_fusion) {
  map_strand <- function(strand) {
    result <- purrr::map_chr(strand, function(x) {
      result <- list(f = "+", r = "-")[[x]]
      return(result)
    })
    return(result)
  }
  tophat_fusion <- dplyr::mutate(
    tophat_fusion,
    Gene1 = Gene1,
    Gene2 = Gene2,
    OrderedFusion = get_ordered(Gene1, Gene2),
    UnorderedFusion = get_unordered(Gene1, Gene2),
    SupportingReads = SpanningReads + SpanningPairs + SpanningPairsInFusion,
    Chr1 = equalize_chrom(Chr1),
    Break1 = Break1,
    Strand1 = map_strand(Strand1),
    Chr2 = equalize_chrom(Chr2),
    Break2 = Break2,
    Strand2 = map_strand(Strand2),

    FMFusionID = character(nrow(tophat_fusion)),
    FrameShiftClass = character(nrow(tophat_fusion)),
    FMKnownTranscript1 = character(nrow(tophat_fusion)),
    FMKnownExonNumber1 = character(nrow(tophat_fusion)),
    FMKnownTranscript2 = character(nrow(tophat_fusion)),
    FMKnownExonNumber2 = character(nrow(tophat_fusion)),
    FusionJunctionSequence = character(nrow(tophat_fusion))
  )
  return(tophat_fusion)
}

equalize_dragen <- function(dragen) {
  dragen$Gene1 <- vapply(strsplit(dragen$FusionName, "--", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  dragen$Gene2 <- vapply(strsplit(dragen$FusionName, "--", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  dragen$OrderedFusion <- get_ordered(dragen$Gene1, dragen$Gene2)
  dragen$UnorderedFusion <- get_unordered(dragen$Gene1, dragen$Gene2)
  dragen$SupportingReads <- lengths(lapply(dragen$ReadNames, grepRaw, pattern = ";", all = TRUE, fixed = TRUE))

  dragen$Chr1 <- vapply(strsplit(dragen$LeftBreakpoint, ":", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  dragen$Chr1 <- equalize_chrom(dragen$Chr1)
  dragen$Break1 <- vapply(strsplit(dragen$LeftBreakpoint, ":", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  dragen$Strand1 <- vapply(strsplit(dragen$LeftBreakpoint, ":", fixed = TRUE), `[`, 3, FUN.VALUE = character(1))
  dragen$Chr2 <- vapply(strsplit(dragen$RightBreakpoint, ":", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  dragen$Chr2 <- equalize_chrom(dragen$Chr2)
  dragen$Break2 <- vapply(strsplit(dragen$RightBreakpoint, ":", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  dragen$Strand2 <- vapply(strsplit(dragen$RightBreakpoint, ":", fixed = TRUE), `[`, 3, FUN.VALUE = character(1))

  dragen$FMFusionID <- character(nrow(dragen))
  dragen$FrameShiftClass <- character(nrow(dragen))
  dragen$FMKnownTranscript1 <- character(nrow(dragen))
  dragen$FMKnownExonNumber1 <- character(nrow(dragen))
  dragen$FMKnownTranscript2 <- character(nrow(dragen))
  dragen$FMKnownExonNumber2 <- character(nrow(dragen))
  dragen$FusionJunctionSequence <- character(nrow(dragen))

  return(dragen)
}

equalize_arriba <- function(arriba) {
  arriba$Gene1 <- arriba$Gene1
  arriba$Gene2 <- arriba$Gene2

  arriba$Chr1 <- vapply(strsplit(arriba$breakpoint1, ":", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  arriba$Chr1 <- equalize_chrom(arriba$Chr1)
  arriba$Break1 <- vapply(strsplit(arriba$breakpoint1, ":", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  arriba$Strand1 <- substr(arriba$strand1, 1, 1)
  arriba$Chr2 <- vapply(strsplit(arriba$breakpoint2, ":", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  arriba$Chr2 <- equalize_chrom(arriba$Chr2)
  arriba$Break2 <- vapply(strsplit(arriba$breakpoint2, ":", fixed = TRUE), `[`, 2, FUN.VALUE = character(1))
  arriba$Strand2 <- substr(arriba$strand2, 3, 3)

  alias_list <- readRDS("/SCRIPTS/R/alias_list_071221.rds")
  arriba$aliasbreak1 <- paste0("chr",arriba$Chr1, ":", arriba$Break1)
  arriba$aliasbreak2 <- paste0("chr",arriba$Chr2, ":", arriba$Break2)
  
  arriba$HGNCSymbol1 <- unlist(lapply(arriba$Gene1, alias_matching, breakpoint = arriba$aliasbreak1, alias_list = alias_list))
  arriba$HGNCSymbol2 <- unlist(lapply(arriba$Gene2, alias_matching, breakpoint = arriba$aliasbreak2, alias_list = alias_list))
  
  arriba$multimatch1 <- ifelse(grepl('\\|', arriba$HGNCSymbol1), "Yes", "No")
  arriba$multimatch2 <- ifelse(grepl('\\|', arriba$HGNCSymbol2), "Yes", "No")
  arriba$alias_match_1 <- ifelse(arriba$Gene1 == arriba$HGNCSymbol1, "Yes", "No")
  arriba$alias_match_2 <- ifelse(arriba$Gene2 == arriba$HGNCSymbol2, "Yes", "No")
  
  arriba$OrderedFusion <- get_ordered(arriba$HGNCSymbol1, arriba$HGNCSymbol2)
  arriba$UnorderedFusion <- get_unordered(arriba$HGNCSymbol1, arriba$HGNCSymbol2)  
  
  arriba$SupportingReads <- arriba$split_reads1 + arriba$split_reads2 + arriba$discordant_mates
  arriba$FMFusionID <- character(nrow(arriba))
  arriba$FrameShiftClass <- arriba$reading_frame
  arriba$FMKnownTranscript1 <- character(nrow(arriba))
  arriba$FMKnownExonNumber1 <- character(nrow(arriba))
  arriba$FMKnownTranscript2 <- character(nrow(arriba))
  arriba$FMKnownExonNumber2 <- character(nrow(arriba))
  arriba$FusionJunctionSequence <- arriba$fusion_transcript

  return(arriba)
}

equalize_cicero <- function(cicero) {
  cicero$Gene1 <- vapply(strsplit(cicero$geneA, ",", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))
  cicero$Gene2 <- vapply(strsplit(cicero$geneB, ",", fixed = TRUE), `[`, 1, FUN.VALUE = character(1))

  cicero$Chr1 <- equalize_chrom(cicero$chrA)
  cicero$Break1 <- cicero$posA
  cicero$Strand1 <- cicero$ortA
  cicero$Chr2 <- equalize_chrom(cicero$chrB)
  cicero$Break2 <- cicero$posB
  cicero$Strand2 <- cicero$ortB
  
  alias_list <- readRDS("/SCRIPTS/R/alias_list_071221.rds")
  cicero$aliasbreak1 <- paste0("chr",cicero$Chr1, ":", cicero$Break1)
  cicero$aliasbreak2 <- paste0("chr",cicero$Chr2, ":", cicero$Break2)
  
  cicero$HGNCSymbol1 <- unlist(lapply(cicero$Gene1, alias_matching, breakpoint = cicero$aliasbreak1, alias_list = alias_list))
  cicero$HGNCSymbol2 <- unlist(lapply(cicero$Gene2, alias_matching, breakpoint = cicero$aliasbreak2, alias_list = alias_list))
  
  cicero$multimatch1 <- ifelse(grepl('\\|', cicero$HGNCSymbol1), "Yes", "No")
  cicero$multimatch2 <- ifelse(grepl('\\|', cicero$HGNCSymbol2), "Yes", "No")
  cicero$alias_match_1 <- ifelse(cicero$Gene1 == cicero$HGNCSymbol1, "Yes", "No")
  cicero$alias_match_2 <- ifelse(cicero$Gene2 == cicero$HGNCSymbol2, "Yes", "No")
  
  cicero$OrderedFusion <- get_ordered(cicero$HGNCSymbol1, cicero$HGNCSymbol2)
  cicero$UnorderedFusion <- get_unordered(cicero$HGNCSymbol1, cicero$HGNCSymbol2)  
  
  cicero$SupportingReads <- cicero$readsA + cicero$readsB
  cicero$FMFusionID <- character(nrow(cicero))
  cicero$FrameShiftClass <- character(nrow(cicero))
  cicero$FMKnownTranscript1 <- character(nrow(cicero))
  cicero$FMKnownExonNumber1 <- character(nrow(cicero))
  cicero$FMKnownTranscript2 <- character(nrow(cicero))
  cicero$FMKnownExonNumber2 <- character(nrow(cicero))
  cicero$FusionJunctionSequence <- cicero$contig

  return(cicero)
}
