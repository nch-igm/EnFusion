#!/bin/env Rscript

library(DBI, quietly = TRUE, warn.conflicts = FALSE)
library(dbplyr, quietly = TRUE, warn.conflicts = FALSE)
library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
library(purrr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(readr, quietly = TRUE, warn.conflicts = FALSE)

srcdir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = T)))
source(file.path(srcdir, "utils_equalize_helpers.R"))
source(file.path(srcdir, "connect_to_db.R"))

# Make database connection

con <- connect_to_db()

# Get fusion and analysis tables

fusion_db <- tbl(con, "fusion")
analysis_db <- tbl(con, "analysis")

# Unordered gene pairs, no chromosomes

ordered_db <- left_join(fusion_db, analysis_db, "analysis_id") %>% distinct(sample_id, 
  gene1, gene2) %>% collect()

unordered_db <- ordered_db %>% mutate(Genes = get_unordered(gene1, gene2)) %>% distinct(sample_id, 
  Genes)

gene_counts <- unordered_db %>% group_by(Genes) %>% summarise(Count = n()) %>% arrange(desc(Count)) %>% 
  mutate(Samples = length(unique(unordered_db$sample_id)), Freq = Count/Samples)

write_tsv(gene_counts, paste0("GenePairCounts_", Sys.Date(), ".tsv"))

# Unordered gene pairs, include chromosomes

get_first <- function(Gene1, Gene2) {
  unordered_fusion <- purrr::map2_int(Gene1, Gene2, function(g1, g2) {
    if (g1 < g2) {
      first <- 1L
    } else {
      first <- 2L
    }
    return(first)
  })
  return(unordered_fusion)
}

match_chrom <- function(Chrom1, Chrom2, First) {
  unordered <- purrr::pmap_chr(list(Chrom1, Chrom2, First), function(c1, c2, first) {
    if (first == 1) {
      result <- paste(c1, c2, sep = "+")
    } else {
      result <- paste(c2, c1, sep = "+")
    }
    return(result)
  })
  return(unordered)
}

unordered_db <- fusion_db %>% left_join(analysis_db, "analysis_id") %>% distinct(sample_id, 
  gene1, chrom1, gene2, chrom2) %>% collect()

gene_counts <- unordered_db %>% mutate(Genes = get_unordered(gene1, gene2), First = get_first(gene1, 
  gene2), Chroms = match_chrom(chrom1, chrom2, First)) %>% distinct(sample_id, 
  Genes, Chroms) %>% group_by(Genes, Chroms) %>% summarise(Count = n()) %>% arrange(desc(Count))

samples_no <- length(unique(unordered_db$sample_id))

gene_counts <- gene_counts %>% mutate(Samples = samples_no, Freq = Count/Samples)

write_tsv(gene_counts, paste0("GenePairCounts_withChroms", Sys.Date(), ".tsv"))

# Ordered gene pairs, no chromosomes

ordered_db <- ordered_db %>% mutate(Fusions = get_ordered(gene1, gene2)) %>% distinct(sample_id, 
  Fusions)

gene_counts_ordered <- ordered_db %>% group_by(Fusions) %>% summarise(Count = n()) %>% 
  arrange(desc(Count)) %>% mutate(Samples = length(unique(ordered_db$sample_id)), 
  Freq = Count/Samples)

write_tsv(gene_counts_ordered, paste0("FusionPairCounts_", Sys.Date(), ".tsv"))

# Ordered gene pairs, include chromosomes

chrom_next <- function(Chrom1, Chrom2) {
  chroms <- paste(Chrom1, Chrom2, sep = ">>")
  return(chroms)
}

ordered_db <- fusion_db %>% left_join(analysis_db, "analysis_id") %>% distinct(sample_id, 
  gene1, chrom1, gene2, chrom2) %>% collect()

gene_counts <- ordered_db %>% mutate(Genes = get_ordered(gene1, gene2), Chroms = chrom_next(chrom1, 
  chrom2)) %>% distinct(sample_id, Genes, Chroms) %>% group_by(Genes, Chroms) %>% 
  summarise(Count = n()) %>% arrange(desc(Count))

samples_no <- length(unique(ordered_db$sample_id))

gene_counts_ordered <- gene_counts %>% mutate(Samples = samples_no, Freq = Count/Samples)

write_tsv(gene_counts_ordered, paste0("FusionPairCounts_withChroms_", Sys.Date(), 
  ".tsv"))
