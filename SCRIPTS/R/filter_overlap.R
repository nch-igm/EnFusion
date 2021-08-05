#!/usr/bin/env Rscript

library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)
library(readr, quietly = TRUE, warn.conflicts = FALSE)

# Process command line input

usage <- "filter_overlap.R --overlap file.tsv"
description <- paste0("filter_overlap.R takes an overlap report and filters the results based on a set of criteria.")
option_list <- list(make_option("--overlap", type = "character", help = "Overlap file. Must be in the current directory ./ (required)"))
opt <- parse_args(OptionParser(usage = usage, description = description, option_list = option_list))
stopifnot(is.character(opt$overlap))
overlap_file <- opt$overlap

# Import overlap report

report <- read_tsv(overlap_file, comment = "#", na = "NA", col_types = cols(.default = col_character(), 
  NumTools = col_integer(), GenePairFrequency = col_double(), GenePairCount = col_integer(), 
  SampleCount = col_integer(), Rank = col_integer(), Total = col_integer(), SupportingReads = col_integer(), 
  Break1 = col_integer(), Break2 = col_integer()))

# Create filtered report

freport <- filter(report, GenePairFrequency < 0.1)  # Database frequency for fusion < 10%.
freport$Distance = ifelse(freport$Chr1 == freport$Chr2 & freport$Strand1 == freport$Strand2, 
  abs((as.numeric(freport$Break1)) - (as.numeric(freport$Break2))), NA)
freport <- filter(freport, Distance > 2e+05 | is.na(Distance))  # Filter out read through events
freport_read <- filter(freport, SupportingReads >= 4)  # At least 1 caller must have 4 reads of evidence.
fusions_with_enough_reads <- unique(freport_read$UnorderedFusion)
freport <- filter(freport, UnorderedFusion %in% fusions_with_enough_reads)
freport <- freport[c("UnorderedFusion", "OrderedFusion", "NumTools", "GenePairFrequency", 
  "GenePairCount", "SampleCount", "Tool", "Rank", "Total", "SupportingReads", "Gene1", 
  "Chr1", "Break1", "Strand1", "Gene2", "Chr2", "Break2", "Strand2", "Distance", 
  "FusionJunctionSequence")]

# Write filtered report

foutReport <- paste0("filtered_", overlap_file)
write_tsv(freport, foutReport)
