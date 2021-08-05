#!/bin/env Rscript

#' Connect to the fusion detection database.
#'
#' .dbconfig.R must set variables dbname, host, port, user, and passwd.
#'
#' @return A database connection.
#' @examples
#' con <- connect_to_db()
connect_to_db <- function() {
  srcdir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = T)))
  suppressWarnings({
    source(file.path(srcdir, ".dbconfig.R"))
  })
  library(RPostgres, quietly = TRUE, warn.conflicts = FALSE)
  
  drv <- Postgres()
  cat(paste0("Trying to connect to DB...\n"))
  cat(paste0("DB Name: ", dbname, "\n"))
  cat(paste0("Host: ", host, "\n"))
  cat(paste0("Port: ", port, "\n"))
  cat(paste0("User: ", user, "\n"))
  con <- dbConnect(drv, dbname = dbname, host = host, port = port, user = user, 
    password = passwd, connect_timeout = 60)
  cat(paste0("Connection to the DB has succeeded!", "\n"))
  cat(paste0("\n"))
  
  return(con)
}
