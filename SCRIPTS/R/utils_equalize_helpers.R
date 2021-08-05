#' Get the ordered fusions given two equalized genes.
#'
#' @param Gene1 A character vector of the equalized first fusion gene partners.
#' @param Gene2 A character vector of the equalized second fusion gene partners.
#' @return A character vector with the ordered fusions 'Gene1->Gene2'.
#' @examples
#' get_ordered(star_fusion$Gene1, star_fusion$Gene2)
get_ordered <- function(Gene1, Gene2) {
  ordered_fusion <- paste(Gene1, Gene2, sep = ">>")
  return(ordered_fusion)
}

#' Get the unordered fusions given two equalized genes.
#'
#' @param Gene1 A character vector of the equalized first fusion gene partners.
#' @param Gene2 A character vector of the equalized second fusion gene partners.
#' @return A character vector with the unordered fusions 'GeneA+GeneB'.
#'   'GeneA' < 'GeneB' evaluates to TRUE so that 'GeneA->GeneB' and
#'   'GeneB->GeneA' are converted to the same unordered fusion.
#' @examples
#' get_unordered(Gene1, Gene2)
get_unordered <- function(Gene1, Gene2) {
  unordered_fusion <- purrr::map2_chr(Gene1, Gene2, function(g1, g2) {
    if (g1 < g2) {
      result <- paste(g1, g2, sep = "+")
    } else {
      result <- paste(g2, g1, sep = "+")
    }
    return(result)
  })
  return(unordered_fusion)
}

#' Convert chromosome records to a standard format.
#'
#' Convert chromosome records to a standard format. Common chromosome formats
#' include 'chrX' and 'X'. This function will use the format 'X', so that a
#' record of 'chrX' -> 'X'.
#'
#' @param chrom A character vector of the chromosomes.
#' @return A character vector with the chromosomes in a standard format.
#' @examples
#' equalize_chrom('chr3')
#' equalize_chrom('3')
#' equalize_chrom(c('chr3', '3'))
equalize_chrom <- function(chrom) {
  chr_opt <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
    "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M")
  pattern <- paste0("(", stringr::str_c(chr_opt, collapse = "|"), ")(?:_.+)?")
  chrom <- purrr::map_chr(chrom, function(chr) {
    if ((nchar(chr) >= 3) & (tolower(substr(chr, 1, 3)) == "chr")) {
      chr <- substr(chr, 4, nchar(chr))
    }
    if (stringr::str_detect(chr, pattern) == FALSE) {
      print(paste0("Chr '", chr, "' is not a valid chr (equalize_chrom)."))
      stop()
    }
    return(chr)
  })
  return(chrom)
}

#' Determine the HGNC IDs and gene symbols that an alias matches to
#' @param gene A string corresponding to an alias
#' @param breakpoint A string corresponding to gene breakpoint
#' @param alias_list A list containing HGNC IDs as keys and gene symbol, alias, 
#' locus_type, chr, and location as values
#' @return A dataframe listing all matching HGNC IDs, their associated symbols
#' and a boolean indicating if an alias maps to two or more HGNC IDs
#' @example alias_matching('RAF1', 'chr3:12600415', alias_list)

alias_matching <- function(gene, breakpoint, alias_list) {
  ids <- c()
  symbols <- c()
  
  for (i in 1:length(alias_list)) {
    if (gene %in% alias_list[[i]]$aliases | gene %in% alias_list[[i]]$symbol) {
      ids <- c(ids, names(alias_list)[i]) # Return IDS
      symbols <- c(symbols, alias_list[[i]]$symbol) # Return symbols
    }
  }
  
  if (gene %in% symbols) { # If alias is also a symbol, end search
    loc <- match(gene, symbols)
    return(symbols[loc])
  }
  
  multimatch <- ifelse(length(ids) > 1, "Yes", "No") # Determine if alias has multiple matches
  
  ids_collapsed <- paste0(ids, collapse = '|')
  symbols_collapsed <- paste0(symbols, collapse = '|')
  
  # Alias has no match
  if (is.null(ids)) {
    return(gene)
  }
  
  
  # If alias maps to one symbol
  if (multimatch == "No") {
    return(symbols_collapsed)
  }
  
  ids_to_keep <- c()
  symbols_to_keep <- c()
  ids_to_remove <- c()
  
  # Extract breakpoint
  breakpoint <- as.integer(unlist(strsplit(breakpoint, '[:]')[[1]][2]))
  
  # Assess breakpoint 
  for (i in 1:length(ids)) {
    if (breakpoint >= alias_list[[ids[i]]][['location']][1] & breakpoint <= 
        alias_list[[ids[i]]][['location']][2]) {
      ids_to_keep <- c(ids_to_keep, ids[i]) # Keep id
      symbols_to_keep <- c(symbols_to_keep, alias_list[[ids[i]]][['symbol']]) # Keep symbol
    }
    else {
      ids_to_remove <- c(ids_to_remove, ids[i])
    }
  }
  
  
  
  if (length(ids_to_keep) == 0) { # If multimatch is True but no match is found, return all IDs
    ids_to_keep <- ids
    symbols_to_keep <- symbols
    ids_to_remove <- c()
  }
  
  ids_collapsed <- paste0(ids_to_keep, collapse = '|')
  idsr_collapsed <- paste0(ids_to_remove, collapse = '|')
  symbols_collapsed <- paste0(symbols_to_keep, collapse = '|')
  
  # Return updated dataframe within dropped ids
  return(symbols_collapsed)
  
}