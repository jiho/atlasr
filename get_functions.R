#
#     Source all available functions and install requirements
#
# (c) Copyright 2012 Jean-Olivier Irisson
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------

# Functions are in the folder "library", which is next to the current file; so we need to detect where we are
# This does it, for reasons I don't fully understand
# source: http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
path <- dirname(frame_files[[length(frame_files)]])

# Detect names of source files
sourceFiles = list.files(paste(path, "/library", sep=""), pattern="\\.(R|r)$", full=TRUE)

# Detect which packages are used in those source files
requiredPackages = c()
for (file in sourceFiles) {
    # read content
    content = scan(file, what="character", quiet=T, sep="\n")
    # extract require calls
    m = regexpr("require\\(.*?\\)", content)
    matched = regmatches(content, m)
    # extract the name of the package inside the require call
    m = regexpr("\\\".*?\\\"", matched)
    matched = regmatches(matched, m)
    # remove quotes
    pack = gsub("\\\"", "", matched)

    # store it as required
    requiredPackages = c(requiredPackages, pack)
}

# Install missing packages
requiredPackages = c(unique(requiredPackages), "mapproj")
installedPackages = row.names(installed.packages())
missingPackages = setdiff(requiredPackages, installedPackages)
if (length(missingPackages) > 0) {
    install.packages(missingPackages, repos="http://cran.at.r-project.org")

    # Install RCurl from the dedicated repository on windows
    if (Sys.info()["sysname"] == "Windows") {
      install.packages("RCurl", repos="http://www.stats.ox.ac.uk/pub/RWin/")
    }
}

# for (pack in requiredPackages) {
#   suppressPackageStartupMessages(library(pack, character.only=TRUE, quietly=TRUE))
# }

# Source required functions
for (file in sourceFiles) {
    source(file)
}

rm(file, content, installedPackages, m, matched, missingPackages, pack, requiredPackages, sourceFiles)

