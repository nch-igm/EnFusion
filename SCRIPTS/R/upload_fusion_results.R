#!/usr/bin/env Rscript

srcdir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = T)))
suppressWarnings({
  source(file.path(srcdir, "utils_import.R"))
  source(file.path(srcdir, "utils_equalize.R"))
  source(file.path(srcdir, "utils_equalize_helpers.R"))
  source(file.path(srcdir, "utils_add_ranks.R"))
  source(file.path(srcdir, "utils_upload.R"))
  source(file.path(srcdir, "connect_to_db.R"))
})

library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(DBI, quietly = TRUE, warn.conflicts = FALSE)
library(purrr, quietly = TRUE, warn.conflicts = FALSE)

# Process input

description <- paste0("This script acts as an importer for fusion detection ", "results to the fusion_detection database.")
option_list <- list(make_option("--subject", type = "character", help = "Name of the subject, e.g. Tumor_pt40_RNA. (required)"), 
  make_option("--sample", type = "character", help = "Name of the sample, e.g. 18-0246-01_0029-A. (required)"), 
  make_option("--phenotype", type = "character", help = "Phenotype of the subject. (optional)"), 
  make_option("--storage", type = "character", help = paste0("Storage method for the sample, such as 'Frozen' ", 
    "or 'FFPE'. (optional)")), make_option("--reads", type = "integer", help = "Number of read pairs for the sample. (optional)"), 
  make_option("--tool", type = "character", help = paste0("Tool for which the results are being entered, ", 
    "possible values are 'StarFusion', 'FusionMap', ", "'TophatFusion', 'FusionCatcher', 'SoapFuse', 'arriba', 'cicero', ", 
    "'EricScript', and 'Jaffa'. (optional)")), make_option("--results", type = "character", 
    help = paste0("Location of the results file for the associated ", "tool, e.g. 'results.txt'. (optional)")))
opt_parser <- OptionParser(description = description, option_list = option_list)
opt <- parse_args(opt_parser)
opt$help <- NULL

stopifnot(is.character(opt$subject))
stopifnot(is.character(opt$sample))

cat(paste0("\n"))
cat(paste0("Running upload_fusion_results.R with the following parameters...", "\n"))
for (param in names(opt)) {
  cat(paste0(str_pad(param, max(nchar(names(opt))), "right"), " : ", opt[[param]], 
    "\n"))
}
cat("\n")

null_to_na <- function(x) {
  if (is.null(x)) {
    x <- NA
  }
  return(x)
}
opt$phenotype <- null_to_na(opt$phenotype)
opt$storage <- null_to_na(opt$storage)
opt$reads <- null_to_na(opt$reads)

# Make database connection

con <- connect_to_db()

# Checking subject and sample data in the database

cat(paste0("Checking the (subject, sample) combination...", "\n"))
cat(paste0("('", opt$subject, "', '", opt$sample, "')", "\n"))
subject_id <- get_subject_id(con, opt$subject)
sample_id <- get_sample_id(con, opt$sample)
cat(paste0("(", subject_id, ", ", sample_id, ")", "\n"))
cat(paste0("\n"))

if (subject_id > 0 & sample_id > 0) {
  query <- paste0("SELECT subject_id FROM sample ", "WHERE sample_id = ", sample_id, 
    ";")
  sample <- dbGetQuery(con, query)
  test_subject_id <- sample$subject_id
  stopifnot(subject_id == test_subject_id)
}

# Upload subject and sample data

if (subject_id < 0 | sample_id < 0) {
  
  cat(paste0("Uploading subject and/or sample data...", "\n"))
  if (subject_id < 0) {
    subject_id <- upload_subject(con, opt$subject, opt$phenotype)
    cat(paste0("Subject (subject_id = ", subject_id, ") has been uploaded", "\n"))
  }
  if (sample_id < 0) {
    sample_id <- upload_sample(con, opt$sample, opt$subject, opt$storage, opt$reads)
    cat(paste0("Sample (sample_id = ", sample_id, ") has been uploaded", "\n"))
  }
  cat(paste0("\n"))
  
}

# Assemble a list of files to upload

cat(paste0("Assembing a list of files to upload...", "\n"))

find_results <- function(base) {
  results_patterns <- list(StarFusion = "star-fusion\\.fusion_predictions\\.abridged\\.tsv", 
    FusionMap = "FusionDetection\\.FusionReport\\.Table\\.txt", TophatFusion = "result\\.txt", 
    FusionCatcher = "final-list_candidate-fusion-genes\\.txt", Jaffa = "jaffa_results\\.csv", 
    MapSplice = "fusions_well_annotated\\.txt", SoapFuse = "\\.final\\.Fusion\\.specific\\.for\\.genes", 
    dragen = "DRAGEN\\.fusion_candidates\\.final", arriba = "fusions\\.tsv",
    cicero = "annotated\\.fusion\\.txt")
  map(results_patterns, function(p) {
    list.files(base, pattern = p, recursive = TRUE, full.names = TRUE)
  })
}

if (is.character(opt$tool) & is.character(opt$results)) {
  upload_files <- list()
  upload_files[[opt$tool]] = opt$results
} else {
  upload_files <- find_results(".")
  upload_files <- upload_files[map_lgl(upload_files, function(x) length(x) > 0)]
  upload_files <- map(upload_files, function(x) x[1])
}

print(names(upload_files))
cat(paste0("\n"))

# Loop through list of files to upload

for (tool in names(upload_files)) {
  
  # Checking analysis data in the database
  
  cat(paste0("==> ", upload_files[[tool]], " <==", "\n"))
  cat(paste0("Checking for an analysis with '", tool, "' as tool...", "\n"))
  analysis_id <- get_analysis_id(con, opt$sample, tool)
  cat(paste0("analysis_id = ", analysis_id, "\n"))
  cat("\n")
  
  # Upload data to the db
  
  if (analysis_id > 0) {
    
    cat(paste0("The analysis is already in the database.", "\n"))
    cat(paste0("\n"))
    
  } else {
    
    lc_tool <- paste0(tolower(str_sub(tool, 1, 1)), str_sub(tool, 2, -1))  # lowercase first letter
    fusions <- import_results(lc_tool, upload_files[[tool]])
    if (nrow(fusions) > 0) {
      fusions <- equalize_results(lc_tool, fusions)
      fusions <- add_ranks(lc_tool, fusions)
      up_fusions <- select(fusions, Rank, Gene1, Gene2, Chr1, Break1, Strand1, 
        Chr2, Break2, Strand2)
      print(head(up_fusions))
      check_fusions(up_fusions)
      cat(paste0("\n"))
    } else {
      print(fusions)
    }
    
    analysis_id <- upload_analysis(con, opt$sample, tool)
    cat(paste0("Analysis (analysis_id = ", analysis_id, ") has been uploaded", 
      "\n"))
    
    if (nrow(fusions) > 0) {
      output <- upload_fusions(con, analysis_id, up_fusions)
      cat(paste0(output, " fusions have been uploaded", "\n"))
      cat(paste0("\n"))
    } else {
      cat(paste0("There were no fusions to upload", "\n"))
      cat(paste0("\n"))
    }
    
  }
}

# Close database connection

cat("Closing the database connection...\n")
result <- dbDisconnect(con)
if (result) {
  cat("Success\n")
} else {
  cat("Failure\n")
}
