#!/usr/bin/env Rscript
print(getwd())

srcdir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = T)))
suppressWarnings({
  source(file.path(srcdir, "utils_import.R"))
  source(file.path(srcdir, "utils_equalize.R"))
  source(file.path(srcdir, "utils_equalize_helpers.R"))
  source(file.path(srcdir, "utils_add_ranks.R"))
  source(file.path(srcdir, "utils_results.R"))
  source(file.path(srcdir, "connect_to_db.R"))

  library(optparse, quietly = TRUE, warn.conflicts = FALSE)
  library(VennDiagram, quietly = TRUE, warn.conflicts = FALSE)
  library(gridExtra, quietly = TRUE, warn.conflicts = FALSE)
  library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  library(readr, quietly = TRUE, warn.conflicts = FALSE)
  library(DBI, quietly = TRUE, warn.conflicts = FALSE)
  library(dbplyr, quietly = TRUE, warn.conflicts = FALSE)
  library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
  library(purrr, quietly = TRUE, warn.conflicts = FALSE)
})

print(srcdir)

################################################################################

# Process Input

usage <- "assemble_results.R --sample s1"
description <- paste0("assemble_results.R reads a set of fusion detection result files and ",
  "produces an aggregated report with the overlapping fusion events predicted ",
  "in the input files. The program searches recursively from the current ", "directory to find known files. The input files can also be specified ",
  "manually.")
option_list <- list(make_option("--sample", type = "character", help = "Name of the sample. (required)"),
  make_option("--frequency", type = "double", help = "Fusion gene pair cohort frequency, used for cutoff filtering. Use percentage as decimal value. [default = 0.10]."),
  make_option("--baseDir", type = "character", help = "Base directory to search for input files. [default = ./]"),
  make_option("--outReport", type = "character", help = "Location of the output report. [default = overlap_$sample.tsv]"),
  make_option("--collapseoutReport", type = "character", help = "Location of the output collapse report. [default = collapse_overlap_$sample.tsv]"),
  make_option("--outSingleton", type = "character", help = "Location of the Singleton output report. [default = Singleton_KnownFusions_$sample.tsv]"),
  make_option("--outmutlimatch", type = "character", help = "Location of the multimatch output report. [default = Multimatch_$sample.tsv]"),
  make_option("--outNoWL", type = "character", help = "Location of the output no knownfusionlist report. [default = no_knownfusionlist_collapse_$sample.tsv]"),
  make_option("--outNoWL2", type = "character", help = "Location of the output no knownfusionlist report. [default = no_knownfusionlist_2callers_collapse_$sample.tsv]"),
  make_option("--foutReport", type = "character", help = "Location of the output filtered, two callers report. [default = filtered_overlap_knownfusionlist_2callers_$sample.tsv]"),
  make_option("--foutReport3", type = "character", help = "Location of the output filtered, three callers report. [default = filtered_overlap_knownfusionlist_3callers_$sample.tsv]"),
  make_option("--ericScript", type = "character", default = NULL, help = "Specific location of the ericScript results file ([sample].results.filtered.tsv)."),
  make_option("--fusionCatcher", type = "character", default = NULL, help = "Specific location of the fusionCatcher results file (final-list_candidate-fusion-genes.txt)."),
  make_option("--fusionMap", type = "character", default = NULL, help = "Specific location of the fusionMap results file (results/FusionDetection.FusionReport.Table.txt)."),
  make_option("--jaffa", type = "character", default = NULL, help = "Specific location of the jaffa results file (jaffa_results.csv)."),
  make_option("--mapSplice", type = "character", default = NULL, help = "Specific location fo the mapSplice results file (fusions_well_annotated.txt)."),
  make_option("--soapFuse", type = "character", default = NULL, help = "Specific location of the soapFuse results file (out/final_fusion_genes/[sample]/[sample].final.Fusion.specific.for.genes)."),
  make_option("--starFusion", type = "character", default = NULL, help = "Specific location of the starFusion results file (star-fusion.fusion_predictions.abridged.tsv)."),
  make_option("--tophatFusion", type = "character", default = NULL, help = "Specific location of the tophatFusion results file (fusion_candidates.final.txt), please note that if file is empty, overlap will fail."),
#  make_option("--dragen", type = "character", default = NULL, help = "Specific location of the dragen results file (DRAGEN.fusion_candidates.final.txt)."),
  make_option("--arriba", type = "character", default = NULL, help = "Specific location of the arriba results file (fusions.tsv)."),
  make_option("--cicero", type = "character", default = NULL, help = "Specific location of the cicero results file (annotated.fusion.txt)."),
  make_option("--frequencyFile", type = "character", default = NULL, help = "Location of the gene pair frequency file. [default = {srcdir}/GenePairCounts_2021-08-05.tsv]"),
  make_option("--frequencyDbConnect", action = "store_true", default = FALSE, help = "Flag to connect to the frequency database instead of use the frequency file."))
opt <- parse_args(OptionParser(usage = usage, description = description, option_list = option_list))

stopifnot(is.character(opt$sample))
sample <- opt$sample

if (is.null(opt$frequency)) {
  frequency <- 0.1
} else {
  frequency <- opt$frequency
}
if (is.null(opt$baseDir)) {
  baseDir <- "."
} else {
  baseDir <- opt$baseDir
}
if (is.null(opt$outReport)) {
  outReport <- paste0("overlap_", sample, ".tsv")
} else {
  outReport <- opt$outReport
}
if (is.null(opt$frequencyFile)) {
  frequencyFile <- file.path(srcdir, "GenePairCounts_2021-08-05.tsv")
} else {
  frequencyFile <- opt$frequencyFile
}
# TODO: add else statements?
if (is.null(opt$collapseoutReport)) {
  collapseoutReport <- paste0("collapse_overlap_", sample, ".tsv")
}
if (is.null(opt$foutReport)) {
  foutReport <- paste0("filtered_overlap_knownfusionlist_2callers_", sample, ".tsv")
}
if (is.null(opt$foutReport3)) {
  foutReport3 <- paste0("filtered_overlap_knownfusionlist_3callers_", sample, ".tsv")
}
if (is.null(opt$outSingleton)) {
  outSingleton <- paste0("Singleton_KnownFusions_", sample, ".tsv")
}
if (is.null(opt$outNoWL)) {
  outNoWL <- paste0("no_knownfusionlist_collapse_", sample, ".tsv")
}
if (is.null(opt$outNoWL2)) {
  outNoWL2 <- paste0("no_knownfusionlist_2callers_collapse_", sample, ".tsv")
}
if (is.null(opt$outmutlimatch)) {
  outmutlimatch <- paste0("Multimatch_", sample, ".tsv")
}

# Assemble a list of files to merge

cat(paste0("Assembing a list of files to merge...", "\n"))

find_results <- function(base) {
  results_patterns <- list(starFusion = "star-fusion\\.fusion_predictions\\.abridged\\.tsv",
    fusionMap = "FusionDetection\\.FusionReport\\.Table\\.txt", tophatFusion = "result\\.txt",
    fusionCatcher = "final-list_candidate-fusion-genes\\.txt", jaffa = "jaffa_results\\.csv",
    mapSplice = "fusions_well_annotated\\.txt", soapFuse = "\\.final\\.Fusion\\.specific\\.for\\.genes",
 #   dragen = "DRAGEN\\.fusion_candidates\\.final",
    arriba = "fusions\\.tsv",
    cicero = "annotated\\.fusion\\.txt")
  map(results_patterns, function(p) {
    list.files(base, pattern = p, recursive = TRUE, full.names = TRUE)
  })
}

list_results <- find_results(baseDir)
#uncomment to test local
#list_results <- find_results(data_dir)

tools <- c("ericScript", "fusionCatcher", "fusionMap", "jaffa", "mapSplice", "soapFuse",
  "starFusion", "tophatFusion", #"dragen",
  "arriba", "cicero")
custom_tools <- match(tools, names(opt))
custom_tools <- custom_tools[!is.na(custom_tools)]
#un comment to test
#custom_input <- tools
custom_input <- opt[custom_tools]

input_list <- list_results
for (tool in tools) {
  if (tool %in% names(custom_input)) {
    input_list[[tool]] <- custom_input[[tool]]
  }
  if (tool %in% names(input_list)) {
    if (length(input_list[[tool]]) > 0) {
      input_list[[tool]] <- input_list[[tool]][1]
    } else {
      input_list[[tool]] <- NULL
    }
  }
}

print(names(input_list))
cat(paste0("\n"))

# Put together all_data list

all_data <- vector("list", length(input_list))
size <- 0
for (tool in names(input_list)) {
  size <- size + 1
  all_data[[size]] <- list(tool = tool, file = input_list[[tool]])
}

if (size == 0) {
  stop("At least 1 tool must be provided as input.")
}

################################################################################

# Import datasets defined in all_data

calls <- vector("list", length(all_data))
for (i in length(all_data):1) {
  new_df <- import_results(all_data[[i]]$tool, all_data[[i]]$file)
  calls[[i]] <- list(tool = all_data[[i]]$tool, calls = nrow(new_df))
  if (all(dim(new_df) == c(0, 0))) {
    all_data[[i]] <- NULL
  } else {
    all_data[[i]]$df <- new_df
    all_data[[i]]$df <- equalize_results(all_data[[i]]$tool, all_data[[i]]$df)
    all_data[[i]]$df <- add_ranks(all_data[[i]]$tool, all_data[[i]]$df)
  }
}

if (length(all_data) == 0) {
  stop("At least 1 of the tools must have non-empty results.")
}

#write.csv(all_data, "test_output.csv")


# Get the fusions

fusions <- list()
for (i in 1:length(all_data)) {
  tool <- all_data[[i]]$tool
  fusions[[tool]] <- unique(all_data[[i]]$df$UnorderedFusion)
}


#print(fusions)

################################################################################

# Get the 'population frequency' for all previously reported fusion pairs

if (opt$frequencyDbConnect) {
  # pull frequency data from the DB if the user requested the DB connection
  con <- connect_to_db()
  
  fusion_db <- tbl(con, "fusion")
  analysis_db <- tbl(con, "analysis")
  
  ordered_db <- left_join(fusion_db, analysis_db, "analysis_id") %>% distinct(sample_id,
                                                                              gene1, gene2) %>% collect()
  
  unordered_db <- ordered_db %>% mutate(unordered = get_unordered(gene1, gene2)) %>%
    distinct(sample_id, unordered)
} else {
  # pull frequency data from the frequency file
  gene_pair_counts <- read_tsv(
    frequencyFile,
    col_types = cols(
      .default = col_character(),
      Count = col_integer(),
      Samples = col_integer(),
      Freq = col_double()))
}

################################################################################
# Get the overlapping genes and results




overlapping <- get_overlapping_from_fusions(fusions)
print(overlapping)
report <- get_report(overlapping, all_data)
report <- order_report(report)
all_breakpoints <- get_report_all_breakpoints(overlapping, all_data)
all_breakpoints <- order_report_all_breakpoints(all_breakpoints, report)

#### for singleton work
fusion_union <- get_union_from_fusions(fusions)
union_report <- get_union_report(fusion_union, all_data)
union_report <- order_union_report(union_report)

#all_union_breakpoints <- get_union_report_all_breakpoints(fusion_union, all_data)
#all_union_breakpoints <- order_union_report_all_breakpoints(all_breakpoints, union_report)
union_multimatch <- union_report

union_report <- filter(union_report, NumTools < 2)


################################################################################

# Add population frequencies to the report
pop_frequency_db <- function(unordered_fusions, unordered_db) {
  result <- purrr::map_chr(unordered_fusions, function(uf) {
    return(sum(unordered_db$unordered == uf)/length(unique(unordered_db$sample_id)))
  })
  return(result)
}
pop_count_db <- function(unordered_fusions, unordered_db) {
  result <- purrr::map_chr(unordered_fusions, function(uf) {
    return(sum(unordered_db$unordered == uf))
  })
}

pop_frequency_file <- function(unordered_fusions, gene_pair_counts) {
  result <- purrr::map_chr(unordered_fusions, function(uf) {
    freq <- gene_pair_counts[gene_pair_counts$Genes == uf, "Freq"]
    if (nrow(freq) == 0) {
      freq <- 0
    }
    return(as.character(freq))
  })
  return(result)
}
pop_count_file <- function(unordered_fusions, gene_pair_counts) {
  result <- purrr::map_chr(unordered_fusions, function(uf) {
    count <- gene_pair_counts[gene_pair_counts$Genes == uf, "Count"]
    if (nrow(count) == 0) {
      count <- 0
    }
    return(as.character(count))
  })
}

if (opt$frequencyDbConnect) {
  report <- report %>% mutate(GenePairFrequency = pop_frequency_db(UnorderedFusion, unordered_db),
    GenePairCount = pop_count_db(UnorderedFusion, unordered_db),
    SampleCount = length(unique(unordered_db$sample_id))) %>%
    select(UnorderedFusion, OrderedFusion, NumTools, GenePairFrequency,
           GenePairCount, SampleCount, HGNCSymbol1, Gene1, alias_match_1, multimatch1, Chr1, Break1,
           Strand1, HGNCSymbol2, Gene2, alias_match_2, multimatch2, Chr2, Break2, Strand2,
           everything())
  union_report <- union_report %>% mutate(GenePairFrequency = pop_frequency_db(UnorderedFusion, unordered_db),
    GenePairCount = pop_count_db(UnorderedFusion, unordered_db),
    SampleCount = length(unique(unordered_db$sample_id))) %>%
      select(UnorderedFusion, OrderedFusion, NumTools, GenePairFrequency,
             GenePairCount, SampleCount, HGNCSymbol1, Gene1, alias_match_1, multimatch1, Chr1, Break1,
             Strand1, HGNCSymbol2, Gene2, alias_match_2, multimatch2, Chr2, Break2, Strand2,
             everything())
} else {
  report <- report %>% mutate(GenePairFrequency = pop_frequency_file(UnorderedFusion, gene_pair_counts),
    GenePairCount = pop_count_file(UnorderedFusion, gene_pair_counts),
    SampleCount = as.character(gene_pair_counts[1,  "Samples"])) %>%
    select(UnorderedFusion, OrderedFusion, NumTools, GenePairFrequency,
           GenePairCount, SampleCount, HGNCSymbol1, Gene1, alias_match_1, multimatch1, Chr1, Break1,
           Strand1, HGNCSymbol2, Gene2, alias_match_2, multimatch2, Chr2, Break2, Strand2,
           everything())
  union_report <- union_report %>% mutate(GenePairFrequency = pop_frequency_file(UnorderedFusion, gene_pair_counts),
    GenePairCount = pop_count_db(UnorderedFusion, gene_pair_counts),
    SampleCount = as.character(gene_pair_counts[1,  "Samples"])) %>%
      select(UnorderedFusion, OrderedFusion, NumTools, GenePairFrequency,
             GenePairCount, SampleCount, HGNCSymbol1, Gene1, alias_match_1, multimatch1, Chr1, Break1,
             Strand1, HGNCSymbol2, Gene2, alias_match_2, multimatch2, Chr2, Break2, Strand2,
             everything())
}

################################################################################

# Create filtered report
knownfusionlist <- unlist(read.table(file = "/SCRIPTS/known_fusion_list.txt", sep = '\t', header = FALSE), use.names = FALSE)
knownsinglegene1 <- read.table(file = "/SCRIPTS/Singleton1List.txt", sep = '\t', header = TRUE)
knownsinglegene2 <- read.table(file = "/SCRIPTS/Singleton2List.txt", sep = '\t', header = TRUE)

# add "KnownFusion" column to annotate if fusion is known pathogenic
report$KnownFusion <- ifelse(report$UnorderedFusion %in% knownfusionlist, "yes", "no")  #create column identifying knownfusionlist genes as 'known fusions'
union_report$KnownFusion <- ifelse(union_report$UnorderedFusion %in% knownfusionlist, "yes", "no")

union_multimatch$Multimatch <- ifelse(grepl('|',union_multimatch$UnorderedFusion, fixed=TRUE), "yes", "no")
union_multimatch <- dplyr::filter(union_multimatch, union_multimatch$Multimatch == "yes")

report$KnownGene1 <- ifelse(report$Gene1 %in% knownsinglegene1$Gene1, "yes", "no")  #create column identifying knownfusionlist genes as 'known fusions'
union_report$KnownGene1 <- ifelse(union_report$Gene1 %in% knownsinglegene1$Gene1, "yes", "no")
report <- dplyr::left_join(report,knownsinglegene1, by="Gene1")
union_report <- dplyr::left_join(union_report,knownsinglegene1, by="Gene1")

report$KnownGene2 <- ifelse(report$Gene2 %in% knownsinglegene2$Gene2, "yes", "no")  #create column identifying knownfusionlist genes as 'known fusions'
union_report$KnownGene2 <- ifelse(union_report$Gene2 %in% knownsinglegene2$Gene2, "yes", "no")
report <- dplyr::left_join(report,knownsinglegene2, by="Gene2")
union_report <- dplyr::left_join(union_report,knownsinglegene2, by="Gene2")


report$Distance = ifelse(report$Chr1 == report$Chr2 & report$Strand1 == report$Strand2,
  abs((as.numeric(report$Break1)) - (as.numeric(report$Break2))), NA)
union_report$Distance = ifelse(union_report$Chr1 == union_report$Chr2 &
  union_report$Strand1 == union_report$Strand2,
  abs((as.numeric(union_report$Break1)) - (as.numeric(union_report$Break2))), NA)

report$ReadThrough = ifelse(report$Distance > 2e+05 | is.na(report$Distance),
  "no", "yes")  #create readthrough column
union_report$ReadThrough = ifelse(union_report$Distance > 2e+05 | is.na(union_report$Distance),
  "no", "yes")  #create readthrough column

report <- report[c("UnorderedFusion", "OrderedFusion", "KnownFusion", "NumTools", "GenePairFrequency",
  "GenePairCount", "SampleCount", "Tool", "Rank", "Total", "SupportingReads", "HGNCSymbol1", "Gene1", "alias_match_1", "multimatch1", "KnownGene1", "Gene1Score", "Gene1Type",
  "Chr1", "Break1", "Strand1", "HGNCSymbol2", "Gene2", "alias_match_2", "multimatch2", "KnownGene2", "Gene2Score", "Gene2Type",
  "Chr2", "Break2", "Strand2", "Distance", "ReadThrough",
  "FusionJunctionSequence")]
union_report <- union_report[c("UnorderedFusion", "OrderedFusion", "KnownFusion", "NumTools", "GenePairFrequency",
  "GenePairCount", "SampleCount", "Tool", "Rank", "Total", "SupportingReads", "HGNCSymbol1", "Gene1", "alias_match_1", "multimatch1", "KnownGene1", "Gene1Score", "Gene1Type",
  "Chr1", "Break1", "Strand1", "HGNCSymbol2", "Gene2", "alias_match_2", "multimatch2", "KnownGene2", "Gene2Score", "Gene2Type",
  "Chr2", "Break2", "Strand2", "Distance", "ReadThrough")]

#filter by frequency and knownfusionlist and create freport, which is the "filtered report"
freport <- filter(report, GenePairFrequency < frequency | !report$KnownFusion == "no")  # Filter based on database frequency for fusion < 10% and not in knownfusionlist
union_report <- filter(union_report, union_report$KnownFusion == "yes")

#filter by readthrough with knownfusionlist
freport_read_throughs <- filter(freport, freport$ReadThrough == "yes")
fusions_read_throughs <- unique(freport_read_throughs$UnorderedFusion)
freport <- filter(freport, !(UnorderedFusion %in% fusions_read_throughs) | !freport$KnownFusion ==
  "no")  # Filter out read through events unless in knownfusionlist

#filter by minimum number of reads with knownfusionlist
freport_read <- filter(freport, freport$SupportingReads >= 4 | !freport$KnownFusion ==
  "no")  # At least 1 caller must have 4 reads of evidence or be in knownfusionlist
fusions_with_enough_reads <- unique(freport_read$UnorderedFusion)
freport <- filter(freport, UnorderedFusion %in% fusions_with_enough_reads)  #filter fusions with out enough reads and not in knownfusionlist

#create comment lines for number of calls per tool
comment_line <- paste0("# Sample: ", sample)
comment_line <- paste0(comment_line, "\n", "# NumToolsAggregated: ", "\t", length(calls))
for (i in 1:length(calls)) {
  comment_line <- paste0(comment_line, "\n", "# - ", calls[[i]]$tool, "Calls = ", "\t",
    calls[[i]]$calls)
}


# Write unfiltered results
report_conn_un <- file(outReport)
writeLines(comment_line, report_conn_un)
close(report_conn_un)
write_tsv(report, outReport, append = TRUE, col_names = TRUE)

# Write filtered results (2 callers required)
report_conn_f <- file(foutReport)
writeLines(comment_line, report_conn_f)
close(report_conn_f)
write_tsv(freport, foutReport, append = TRUE, col_names = TRUE)

# Write filtered results (3 callers required)
freport3 <- filter(freport, freport$NumTools > 2 | !freport$KnownFusion == "no")  #filter out fusions with less than 3 callers unless in knownfusionlist

report_conn_f3 <- file(foutReport3)
writeLines(comment_line,report_conn_f3)
close(report_conn_f3)
write_tsv(freport3, foutReport3, append = TRUE, col_names = TRUE)

# Write union_report for singletons
singleton_con <- file(outSingleton)
writeLines(comment_line, singleton_con)
close(singleton_con)
write_tsv(union_report, outSingleton, append = TRUE, col_names = TRUE)

# Write union_report for multimatch
multimatch_con <- file(outmutlimatch)
writeLines(comment_line, multimatch_con)
close(multimatch_con)
write_tsv(union_multimatch, outmutlimatch, append = TRUE, col_names = TRUE)

### write results that are collapsed
#collapse function:
collapse_on_rank <- function(collapsed_data) {
  top_rank <- report %>%
    dplyr::group_by(UnorderedFusion) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup()
  top_rank_rows <- top_rank %>%
    dplyr::select(UnorderedFusion, OrderedFusion, KnownFusion, NumTools, GenePairFrequency,
      GenePairCount, SampleCount, HGNCSymbol1, Gene1, alias_match_1, multimatch1, KnownGene1, Gene1Score, Gene1Type, Chr1, Break1,
      Strand1, HGNCSymbol2, Gene2, alias_match_2, multimatch2, KnownGene2, Gene2Score, Gene2Type, Chr2, Break2, Strand2, Distance) %>%
    dplyr::group_by(UnorderedFusion, OrderedFusion, KnownFusion, NumTools, GenePairFrequency,
      GenePairCount, SampleCount, HGNCSymbol1, Gene1, alias_match_1, multimatch1, KnownGene1, Gene1Score, Gene1Type, Chr1, Break1,
      Strand1, HGNCSymbol2, Gene2, alias_match_2, multimatch2, KnownGene2, Gene2Score, Gene2Type, Chr2, Break2, Strand2, Distance) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup()
  top_rank_rows
}



#if there are results in the report, collapse and write a report tsv
if (nrow(report) > 0) {
  tools <- aggregate(Tool ~ UnorderedFusion, data = report, paste, collapse = ",")
  collapsed_df <- collapse_on_rank(report)
  merged_collapsed_df <- merge(collapsed_df, tools, "UnorderedFusion")
  collapsed_report <- merged_collapsed_df[order(merged_collapsed_df$NumTools, decreasing = TRUE),]
  # write collasped report on collapsed unfiltered overlap data
  report_conn_un_consensus <- file(collapseoutReport)
  writeLines(comment_line, report_conn_un_consensus)
  close(report_conn_un_consensus)
  write_tsv(collapsed_report, collapseoutReport, append = TRUE, col_names = TRUE)
  } else {
    print("no data to collapse for report")
}

if (nrow(freport) > 0) {
  ftools <- aggregate(Tool ~ UnorderedFusion, data = freport, paste, collapse = ",")
  fcollapsed_df <- collapse_on_rank(freport)
  fmerged_collapsed_df <- merge(fcollapsed_df, ftools, "UnorderedFusion")
  fcollapsed_report <- fmerged_collapsed_df[order(fmerged_collapsed_df$NumTools, decreasing = TRUE),]
  foutReport_col <- paste0("collapse_",foutReport)
  report_conn_f_col <- file(foutReport_col)
  writeLines(comment_line, report_conn_f_col)
  close(report_conn_f_col)
  write_tsv(fcollapsed_report, foutReport_col, append = TRUE, col_names = TRUE)
  fcollapsed_report3 <- filter(fcollapsed_report, fcollapsed_report$NumTools > 2 | !fcollapsed_report$KnownFusion == "no")
  unique_filtered_fusions <-  length(fcollapsed_report3$UnorderedFusion)
  comment_line <- paste0(comment_line, "\n", "# filtered_overlap = ", "\t", unique_filtered_fusions)
  print(comment_line)
  foutReport_col3 <- paste0("collapse_",foutReport3)
  report_conn_f_col3 <- file(foutReport_col3)
  writeLines(comment_line, report_conn_f_col3)
  close(report_conn_f_col3)
  write_tsv(fcollapsed_report3, foutReport_col3, append = TRUE, col_names = TRUE)
  } else {
    print("no data to collapse for freport")
}



##### RUN WITHOUT knownfusionlist, 3 caller minimum, AND COLLAPSE FOR PAPER #####
#filter by frequency no knownfusionlist
reportNW <- filter(report, GenePairFrequency < frequency)  # Database frequency for fusion < 10%.

reportNW_read_throughs <- filter(reportNW, reportNW$ReadThrough == "yes")
reportNW_read_throughs <- unique(reportNW_read_throughs$UnorderedFusion)
reportNW <- filter(reportNW, !(UnorderedFusion %in% reportNW_read_throughs))  # Filter out read through events unless in knownfusionlist

reportNW_read <- filter(reportNW, SupportingReads >= 4)  # At least 1 caller must have 4 reads of evidence.
fusions_with_enough_reads <- unique(reportNW_read$UnorderedFusion)
reportNW <- filter(reportNW, UnorderedFusion %in% fusions_with_enough_reads)

reportNW2 <- reportNW

NW2tools <- aggregate(Tool ~ UnorderedFusion, data = reportNW2, paste, collapse = ",")
NW2collapsed_df <- collapse_on_rank(NW2tools)
NW2collapsed_df <- merge(NW2collapsed_df, NW2tools, "UnorderedFusion")
NW2collapsed_df_report <- NW2collapsed_df[order(NW2collapsed_df$NumTools, decreasing = TRUE),]

if (nrow(NW2collapsed_df_report) > 0) {
    print(NW2collapsed_df_report)
  # Write filtered report
    write_tsv(NW2collapsed_df_report, outNoWL2)
    } else {
      print("no data to collapse for nowhitelist2 report")
}

reportNW <- filter(reportNW, NumTools >2)
#print(reportNW)

NWtools <- aggregate(Tool ~ UnorderedFusion, data = reportNW, paste, collapse = ",")
NWcollapsed_df <- collapse_on_rank(NWtools)
NWcollapsed_df <- merge(NWcollapsed_df, NWtools, "UnorderedFusion")
NWcollapsed_df_report <- NWcollapsed_df[order(NWcollapsed_df$NumTools, decreasing = TRUE),]
#print(NWcollapsed_df_report)



if (nrow(NWcollapsed_df_report) > 0) {
    print(NWcollapsed_df_report)
  # Write filtered report
  write_tsv(NWcollapsed_df_report, outNoWL)
    } else {
      print("no data to collapse for nowhitelist report")
}






################################################################################
