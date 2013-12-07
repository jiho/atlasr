#
#      Set global variables/functions for the Shiny app
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------


debug <- T
verbose <- F

# debug message
dmess <- function(...) {
  if (debug) {
    message("DEBUG: ", ...)
  }
  return(invisible(1))
}

vmess <- function(...) {
  if (verbose) {
    message("-> ", ...)
  }
  return(invisible(1))
}
