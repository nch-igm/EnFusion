#' Import fusion detection results.
#'
#' Import fusion detection results for the given tool from the specified
#' location. The import file must be in a valid format for the tool. Note: The
#' import file for the tophatFusion results is named result.txt. There must be
#' an additional file named potential_fusion.txt in the same folder in order for
#' the import to succeed.
#'
#' @param tool The tool that produced the results.
#' @param location The location of the results.
#' @return A tibble with the results. If the results contain no rows of data
#'   then a 0x0 tibble will be returned.
#' @export
#' @examples
#' import_results('ericScript', 'ericscript-sample.results.filtered.tsv')
#' import_results('fusionCatcher', 'final-list_candidate-fusion-genes.txt')
#' import_results('fusionMap', 'fusionmap-sample.FusionReport.txt')
#' import_results('jaffa', 'jaffa_results.csv')
#' import_results('mapSplice', 'fusions_well_annotated.txt')
#' import_results('soapFuse', 'soapfuse-sample.final.Fusion.specific.for.genes')
#' import_results('starFusion', 'star-fusion.fusion_predictions.abridged.tsv')
#' import_results('tophatFusion', 'result.txt')
#' import_results('dragen', 'DRAGEN.fusion_candidates.final')
import_results <- function(tool, location) {
  if (tool == "ericScript") {
    results <- import_eric_script(location)
  } else if (tool == "fusionCatcher") {
    results <- import_fusion_catcher(location)
  } else if (tool == "fusionMap") {
    results <- import_fusion_map(location)
  } else if (tool == "jaffa") {
    results <- import_jaffa(location)
  } else if (tool == "mapSplice") {
    results <- import_map_splice(location)
  } else if (tool == "soapFuse") {
    results <- import_soap_fuse(location)
  } else if (tool == "starFusion") {
    results <- import_star_fusion(location)
  } else if (tool == "tophatFusion") {
    results <- import_tophat_fusion(location)
  } else if (tool == "dragen") {
    results <- import_dragen(location)
  } else if (tool == "arriba") {
    results <- import_arriba(location)
  } else if (tool == "cicero") {
    results <- import_cicero(location)
  } else {
    print(paste0("Tool ", tool, " is not a valid tool (import_results)."))
    stop()
  }
  return(results)
}

import_eric_script <- function(location) {
  eric_script <- readr::read_delim(location, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_types = readr::cols(.default = readr::col_character(), crossingreads = readr::col_integer(), 
      spanningreads = readr::col_integer(), mean.insertsize = readr::col_double(), 
      GeneExpr1 = readr::col_double(), GeneExpr2 = readr::col_double(), GeneExpr_Fused = readr::col_double(), 
      ES = readr::col_double(), GJS = readr::col_double(), US = readr::col_double(), 
      EricScore = readr::col_double()))
  if (nrow(eric_script) == 0) {
    eric_script <- tibble::tibble()
  }
  return(eric_script)
}

import_fusion_catcher <- function(location) {
  fusion_catcher <- readr::read_delim(location, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_types = readr::cols(.default = readr::col_character(), Counts_of_common_mapping_reads = readr::col_integer(), 
      Spanning_pairs = readr::col_integer(), Spanning_unique_reads = readr::col_integer(), 
      Longest_anchor_found = readr::col_integer()))
  if (nrow(fusion_catcher) == 0) {
    fusion_catcher <- tibble::tibble()
  }
  return(fusion_catcher)
}

import_fusion_map <- function(location) {
  fusion_map <- readr::read_delim(location, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_types = readr::cols(.default = readr::col_integer(), FusionID = readr::col_character(), 
      Strand = readr::col_character(), Chromosome1 = readr::col_character(), 
      Chromosome2 = readr::col_character(), KnownGene1 = readr::col_character(), 
      KnownTranscript1 = readr::col_character(), KnownExonNumber1 = readr::col_number(), 
      KnownTranscriptStrand1 = readr::col_character(), KnownGene2 = readr::col_character(), 
      KnownTranscript2 = readr::col_character(), KnownExonNumber2 = readr::col_number(), 
      KnownTranscriptStrand2 = readr::col_character(), FusionJunctionSequence = readr::col_character(), 
      FusionGene = readr::col_character(), SplicePattern = readr::col_character(), 
      SplicePatternClass = readr::col_character(), FrameShift = readr::col_character(), 
      FrameShiftClass = readr::col_character(), OnExonBoundary = readr::col_character(), 
      Filter = readr::col_character()))
  if (nrow(fusion_map) == 0) {
    fusion_map <- tibble::tibble()
  }
  return(fusion_map)
}

import_jaffa <- function(location) {
  jaffa <- readr::read_delim(location, ",", escape_double = TRUE, trim_ws = TRUE, 
    col_types = readr::cols(.default = readr::col_character(), base1 = readr::col_integer(), 
      base2 = readr::col_integer(), `gap (kb)` = readr::col_double(), `spanning pairs` = readr::col_integer(), 
      `spanning reads` = readr::col_integer(), inframe = readr::col_logical(), 
      aligns = readr::col_logical(), rearrangement = readr::col_logical(), 
      `contig break` = readr::col_integer()))
  if (nrow(jaffa) == 0) {
    jaffa <- tibble::tibble()
  }
  return(jaffa)
}

import_map_splice <- function(location) {
  map_splice_colnames <- c("chrom", "doner_end", "acceptor_start", "id", "coverage", 
    "strand", "rgb", "block_count", "block_size", "block_distance", "entropy", 
    "flank_case", "flank_string", "min_mismatch", "max_mismatch", "ave_mismatch", 
    "max_min_suffix", "max_min_prefex", "min_anchor_difference", "unique_read_count", 
    "multi_read_count", "paired_read_count", "left_paired_read_count", "right_paired_read_count", 
    "multiple_paired_read_count", "unique_paired_read_count", "single_read_count", 
    "encompassing_read_pair_count", "donor_start", "acceptor_end", "doner_isoforms", 
    "acceptor_isoforms", "donor_uniformity_score_obsolete", "acceptor_uniformity_score_obsolete", 
    "donor_uniformity_KS-test_score_obsolete", "acceptor_uniformity_KS-test_score_obsolete", 
    "minimal_doner_isoform_length", "maximal_donor_isoform_length", "minimal_acceptor_isoform_length", 
    "maximal_acceptor_isoform_length", "paired_reads_entropy", "mismatch_per_bp", 
    "anchor_score", "max_doner_fragment", "max_acceptor_fragment", "max_cur_fragment", 
    "min_cur_fragment", "ave_cur_fragment", "donor_encompass_unique", "donor_encompass_multiple", 
    "acceptor_encompass_unique", "acceptor_encompass_multiple", "donor_match_to_normal", 
    "acceptor_match_to_normal", "donor_seq", "acceptor_seq", "match_gene_strand", 
    "annotated_type", "fusion_type", "gene_strand", "annotated_gene_donor", "annotated_gene_acceptor", 
    "blank")
  map_splice <- readr::read_delim(location, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_names = map_splice_colnames, col_types = readr::cols(.default = readr::col_integer(), 
      chrom = readr::col_character(), id = readr::col_character(), strand = readr::col_character(), 
      rgb = readr::col_character(), block_size = readr::col_character(), block_distance = readr::col_character(), 
      entropy = readr::col_double(), flank_string = readr::col_character(), 
      ave_mismatch = readr::col_double(), doner_isoforms = readr::col_character(), 
      acceptor_isoforms = readr::col_character(), `donor_uniformity_KS-test_score_obsolete` = readr::col_double(), 
      `acceptor_uniformity_KS-test_score_obsolete` = readr::col_double(), paired_reads_entropy = readr::col_double(), 
      mismatch_per_bp = readr::col_double(), anchor_score = readr::col_double(), 
      min_cur_fragment = readr::col_double(), ave_cur_fragment = readr::col_double(), 
      donor_match_to_normal = readr::col_character(), acceptor_match_to_normal = readr::col_character(), 
      donor_seq = readr::col_character(), acceptor_seq = readr::col_character(), 
      annotated_type = readr::col_character(), fusion_type = readr::col_character(), 
      gene_strand = readr::col_character(), annotated_gene_donor = readr::col_character(), 
      annotated_gene_acceptor = readr::col_character(), blank = readr::col_skip()))
  return(map_splice)
}

import_soap_fuse <- function(location) {
  soap_fuse <- readr::read_delim(location, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_types = readr::cols(.default = readr::col_character(), up_Genome_pos = readr::col_integer(), 
      dw_Genome_pos = readr::col_integer(), Span_reads_num = readr::col_integer(), 
      Junc_reads_num = readr::col_integer()))
  if (nrow(soap_fuse) == 0) {
    soap_fuse <- tibble::tibble()
  }
  return(soap_fuse)
}

import_star_fusion <- function(location) {
  star_fusion <- readr::read_delim(location, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_type = readr::cols(.default = readr::col_character(), JunctionReadCount = readr::col_integer(), 
      SpanningFragCount = readr::col_integer(), LeftBreakEntropy = readr::col_double(), 
      RightBreakEntropy = readr::col_double(), FFPM = readr::col_double()))
  if (nrow(star_fusion) == 0) {
    star_fusion <- tibble::tibble()
  } else {
    colnames(star_fusion)[1] <- "FusionName"
  }
  return(star_fusion)
}

import_tophat_fusion <- function(location) {
  tophat_fusion_colnames <- c("Sample", "Gene1", "Chr1", "Break1", "Gene2", "Chr2", 
    "Break2", "SpanningReads", "SpanningPairs", "SpanningPairsInFusion", "Score")
  tophat_fusion <- readr::read_delim(location, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_names = tophat_fusion_colnames, col_type = readr::cols(.default = readr::col_character(), 
      Break1 = readr::col_integer(), Break2 = readr::col_integer(), SpanningReads = readr::col_integer(), 
      SpanningPairs = readr::col_integer(), SpanningPairsInFusion = readr::col_integer(), 
      Score = readr::col_double()))
  # assume potential_fusion.txt is in the same directory as result.txt
  potential_loc <- file.path(dirname(location), "potential_fusion.txt")
  potential_lines <- readr::read_lines(potential_loc)
  potential_tbl <- tibble::tibble(Sample = character(), Gene1 = character(), Chr1 = character(), 
    Break1 = integer(), Strand1 = character(), Gene2 = character(), Chr2 = character(), 
    Break2 = integer(), Strand2 = character(), SpanningReads = integer(), SpanningPairs = integer(), 
    SpanningPairsInFusion = integer())
  if (length(potential_lines) > 0) {
    for (i in 1:length(potential_lines)) {
      if (i%%6 == 1) {
        p_split <- stringr::str_split(potential_lines[i], " ", simplify = TRUE)
        p_sample <- p_split[, 1]
        p_chr <- stringr::str_split(p_split[, 2], "-", simplify = TRUE)
        p_chr1 <- p_chr[, 1]
        p_chr2 <- p_chr[, 2]
        p_break1 <- as.integer(p_split[, 3])
        p_break2 <- as.integer(p_split[, 4])
        p_strand1 <- stringr::str_sub(p_split[, 5], 1, 1)
        p_strand2 <- stringr::str_sub(p_split[, 5], 2, 2)
        p_spanning_reads <- as.integer(p_split[, 6])
        p_spanning_pairs <- as.integer(p_split[, 7])
        p_spanning_pairs_in_fusion <- as.integer(p_split[, 8])
      } else if (i%%6 == 5) {
        p_split <- stringr::str_split(potential_lines[i], " ", simplify = TRUE)
        p_gene1 <- p_split[, 1]
        p_gene2 <- p_split[, 3]
      } else if (i%%6 == 0) {
        potential_tbl <- tibble::add_row(potential_tbl, Sample = p_sample, 
          Gene1 = p_gene1, Chr1 = p_chr1, Break1 = p_break1, Strand1 = p_strand1, 
          Gene2 = p_gene2, Chr2 = p_chr2, Break2 = p_break2, Strand2 = p_strand2, 
          SpanningReads = p_spanning_reads, SpanningPairs = p_spanning_pairs, 
          SpanningPairsInFusion = p_spanning_pairs_in_fusion)
      }
    }
  }
  tophat_fusion <- dplyr::left_join(tophat_fusion, potential_tbl, by = c("Sample", 
    "Gene1", "Chr1", "Break1", "Gene2", "Chr2", "Break2", "SpanningReads", "SpanningPairs", 
    "SpanningPairsInFusion"))
  stopifnot(!any(is.na(tophat_fusion$Strand1) | is.na(tophat_fusion$Strand2)))
  return(tophat_fusion)
}

import_dragen <- function(location) {
  dragen <- readr::read_delim(location, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_type = readr::cols(.default = readr::col_character(), Score = readr::col_integer()))
  if (nrow(dragen) == 0) {
    dragen <- tibble::tibble()
  } else {
    colnames(dragen)[1] <- "FusionName"
  }
  return(dragen)
}

import_arriba <- function(location) {
  arriba <- readr::read_delim(location, "\t", escape_double = FALSE, trim_ws = TRUE, 
    col_type = readr::cols(.default = readr::col_character(), split_reads1 = readr::col_integer(), 
      split_reads2 = readr::col_integer(), discordant_mates = readr::col_integer()))
  colnames(arriba) <- c("Gene1", "Gene2", "strand1", "strand2", "breakpoint1", 
    "breakpoint2", "site1", "site2", "type", "direction1", "direction2", "split_reads1", 
    "split_reads2", "discordant_mates", "coverage1", "coverage2", "confidence", 
    "closest_genomic_breakpoint1", "closest_genomic_breakpoint2", "filters", 
    "fusion_transcript", "reading_frame", "peptide_sequence", "read_identifiers")
  head(arriba)
  if (nrow(arriba) == 0) {
    arriba <- tibble::tibble()
  }
  return(arriba)
}

import_cicero <- function(location) {
  cicero <- readr::read_delim(
    location, "\t", escape_double = FALSE, trim_ws = TRUE, na = "",
    col_type = readr::cols(
      .default = readr::col_character(),
      readsA = readr::col_integer(), 
      readsB = readr::col_integer()))
  colnames(cicero) <- c("sample", "geneA", "chrA", "posA", "ortA", "featureA",
    "geneB", "chrB", "posB", "ortB", "featureB", "sv_ort", "readsA", "readsB",
    "matchA", "matchB", "repeatA", "repeatB", "coverageA", "coverageB",
    "ratioA", "ratioB", "qposA", "qposB", "total_readsA", "total_readsB",
    "contig", "type")
  head(cicero)
  if (nrow(cicero) == 0) {
    cicero <- tibble::tibble()
  }
  return(cicero)
}
