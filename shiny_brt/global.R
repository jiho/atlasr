#
#      Set global variables/functions for the Shiny app
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------


DEBUG <- TRUE
# DEBUG <- FALSE

# debug message
dmess <- function(...) {
  if (DEBUG) {
    message("DEBUG: ", ...)
  }
  return(invisible(1))
}
