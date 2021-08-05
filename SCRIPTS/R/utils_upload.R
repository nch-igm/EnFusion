library(DBI, quietly = TRUE, warn.conflicts = FALSE)
library(RPostgres, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)

#' Get the ID of the subject from the database.
#'
#' Get the ID of the subject from the PostgreSQL database if the subject
#' already exists, otherwise return -1.
#'
#' @param con The PostgreSQL database connection.
#' @param subject_name The name of the subject.
#' @return The subject_id of the subject, or -1 if it doesn't exist in the
#'   database.
#' @examples
#' get_subject_id(con, 'Tumor_pt40_RNA')
get_subject_id <- function(con, subject_name) {
  query <- paste0("SELECT subject_id FROM subject ", "WHERE subject_name = '", 
    subject_name, "';")
  res <- dbGetQuery(con, query)
  subject_id <- -1
  if (nrow(res) == 1) {
    subject_id <- res$subject_id
  }
  return(subject_id)
}

#' Get the ID of the sample from the database.
#'
#' Get the ID of the sample from the PostgreSQL database if the sample
#' already exists, otherwise return -1.
#'
#' @param con The PostgreSQL database connection.
#' @param sample_name The name of the sample.
#' @return The sample_id of the sample, or -1 if it doesn't exist in the
#'   database.
#' @examples
#' get_sample_id(con, '18-0246-01_0029-A')
get_sample_id <- function(con, sample_name) {
  query <- paste0("SELECT sample_id FROM sample ", "WHERE sample_name = '", sample_name, 
    "';")
  res <- dbGetQuery(con, query)
  sample_id <- -1
  if (nrow(res) == 1) {
    sample_id <- res$sample_id
  }
  return(sample_id)
}

#' Get the ID of the tool from the database.
#'
#' Get the ID of the tool from the PostgreSQL database.
#'
#' @param con The PostgreSQL database connection.
#' @param tool The name of the tool.
#' @return The tool_id of the tool, or -1 if it doesn't exist in the database.
#' @examples
#' get_tool_id(con, 'FusionMap')
get_tool_id <- function(con, tool) {
  query <- paste0("SELECT tool_id FROM tool ", "WHERE tool_name = '", tool, "';")
  res <- dbGetQuery(con, query)
  tool_id <- -1
  if (nrow(res) == 1) {
    tool_id <- res$tool_id
  }
  return(tool_id)
}

#' Get the ID of the analysis from the database.
#'
#' Get the ID of the analysis from the PostgreSQL database if the analysis
#' already exists, otherwise return -1.
#'
#' @param con The PostgreSQL database connection.
#' @param sample_name The name of the sample.
#' @param tool The name of the tool.
#' @return The analysis_id of the analysis, or -1 if it doesn't exist in the
#'   database.
#' @examples
#' get_analysis_id(con, '18-0246-01_0029-A', 'FusionMap')
get_analysis_id <- function(con, sample_name, tool) {
  sample_id <- get_sample_id(con, sample_name)
  tool_id <- get_tool_id(con, tool)
  query <- paste0("SELECT analysis_id FROM analysis ", "WHERE sample_id = ", sample_id, 
    " AND tool_id = ", tool_id, ";")
  res <- dbGetQuery(con, query)
  analysis_id <- -1
  if (nrow(res) == 1) {
    analysis_id <- res$analysis_id
  }
  return(analysis_id)
}

#' Upload the subject information.
#'
#' Upload the subject information to the PostgreSQL database, if it isn't
#' already in the database.
#'
#' @param con The PostgreSQL database connection.
#' @param subject_name The name of the subject to add to the database.
#' @param phenotype (optional) The phenotype of the subject.
#' @return The subject_id of the subject, whether it was uploaded itself or
#'   already present in the database.
#' @examples
#' upload_subject(con, 'Tumor_pt40_RNA')
upload_subject <- function(con, subject_name, phenotype = NA) {
  subject_id <- get_subject_id(con, subject_name)
  if (subject_id == -1) {
    subject <- tibble(subject_name, phenotype)
    dbWriteTable(con, "subject", subject, append = TRUE, row.names = FALSE)
    subject_id <- get_subject_id(con, subject_name)
  }
  return(subject_id)
}

#' Upload the sample information.
#'
#' Upload the sample information to the PostgreSQL database, if it isn't
#' already in the database.
#'
#' @param con The PostgreSQL database connection.
#' @param sample_name The name of the sample to add to the database.
#' @param subject_name The name of the subject associated with the sample.
#' @param storage_method (optional) The storage method of the sample.
#' @param num_reads (optional) The number of reads for the sample.
#' @return The sample_id of the subject, whether it was uploaded itself or
#'   already present in the database.
#' @examples
#' upload_sample(con, '18-0246-01_0029-A', 'Tumor_pt40_RNA')
upload_sample <- function(con, sample_name, subject_name, storage_method = NA, num_reads = NA) {
  sample_id <- get_sample_id(con, sample_name)
  if (sample_id == -1) {
    subject_id <- get_subject_id(con, subject_name)
    sample <- tibble(subject_id, sample_name, storage_method, num_reads)
    dbWriteTable(con, "sample", sample, append = TRUE, row.names = FALSE)
    sample_id <- get_sample_id(con, sample_name)
  }
  return(sample_id)
}

#' Upload the analysis information.
#'
#' Upload the analysis information to the PostgreSQL database, if it isn't
#' already in the database.
#'
#' @param con The PostgreSQL database connection.
#' @param sample_name The name of the sample.
#' @param tool The name of the tool.
#' @param tool_version (optional) The version of the tool used in the analysis.
#' @return The analysis_id of the analysis, whether it was uploaded itself or
#'   already present in the database.
#' @examples
#' upload_analysis(con, '18-0246-01_0029-A', 'FusionMap')
upload_analysis <- function(con, sample_name, tool, tool_version = NA) {
  analysis_id <- get_analysis_id(con, sample_name, tool)
  if (analysis_id == -1) {
    sample_id <- get_sample_id(con, sample_name)
    tool_id <- get_tool_id(con, tool)
    analysis <- tibble(sample_id, tool_id, tool_version)
    dbWriteTable(con, "analysis", analysis, append = TRUE, row.names = FALSE)
    analysis_id <- get_analysis_id(con, sample_name, tool)
  }
  return(analysis_id)
}

#' Check that the fusions are ready for upload.
#'
#' Check that the fusions are formatted correctly so they can be uploaded to the
#' database.
#'
#' @param fusions A data.frame of fusions to upload to the database.
#' @return TRUE or FALSE for whether it passes the test.
#' @examples
#' check_fusions(fusions)
check_fusions <- function(fusions) {
  chroms <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
    "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M")
  pattern <- paste0("(", str_c(chroms, collapse = "|"), ")(?:_.+)?")
  strands <- c("+", "-", "+/+", "+/-", "+/.", "-/+", "-/-", "-/.", "./+", "./-", "./.")
  
  stopifnot("Rank" %in% colnames(fusions))
  stopifnot("Gene1" %in% colnames(fusions))
  stopifnot("Gene2" %in% colnames(fusions))
  stopifnot("Chr1" %in% colnames(fusions))
  stopifnot("Break1" %in% colnames(fusions))
  stopifnot("Strand1" %in% colnames(fusions))
  stopifnot("Chr2" %in% colnames(fusions))
  stopifnot("Break2" %in% colnames(fusions))
  stopifnot("Strand2" %in% colnames(fusions))
  
  stopifnot(all(fusions$Rank > 0))
  stopifnot(all(fusions$Break1 > 0 | is.na(fusions$Break1)))
  stopifnot(all(fusions$Break2 > 0 | is.na(fusions$Break2)))
  
  stopifnot(all(str_detect(fusions$Chr1, pattern)))
  stopifnot(all(str_detect(fusions$Chr2, pattern)))
  
  stopifnot(all(fusions$Strand1 %in% strands))
  stopifnot(all(fusions$Strand2 %in% strands))
}

#' Upload the fusion data.
#'
#' Upload the fusion data to the PostgreSQL database.
#'
#' @param con The PostgreSQL database connection.
#' @param analysis_id The id of the analysis that reported the fusions.
#' @param fusions The fusions to upload.
#' @return The number of fusions uploaded to the db.
#' @examples
#' upload_fusions(con, 1438, up_fusions)
upload_fusions <- function(con, analysis_id, fusions) {
  fusions <- mutate(fusions, analysis_id = analysis_id, fusion_rank = Rank, gene1 = Gene1, 
    chrom1 = Chr1, break1 = Break1, strand1_plus = (Strand1 == "+"), gene2 = Gene2, 
    chrom2 = Chr2, break2 = Break2, strand2_plus = (Strand2 == "+"))
  fusions <- select(fusions, analysis_id, fusion_rank, gene1, chrom1, break1, strand1_plus, 
    gene2, chrom2, break2, strand2_plus)
  no_fusions_before <- dbGetQuery(con, "SELECT COUNT(*) as count from fusion;")
  output <- dbWriteTable(con, "fusion", fusions, append = TRUE, row.names = FALSE)
  no_fusions_after <- dbGetQuery(con, "SELECT COUNT(*) as count from fusion;")
  diff_fusions <- as.integer(no_fusions_after$count - no_fusions_before$count)
  return(diff_fusions)
}
