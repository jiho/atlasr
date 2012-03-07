#
#     Source all available functions and install requirements
#
# (c) Copyright 2012 Jean-Olivier Irisson
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------

# Detect names of source files
sourceFiles = list.files("library", pattern=glob2rx("*.R|r"), full=TRUE)

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
}

# for (pack in requiredPackages) {
#   suppressPackageStartupMessages(library(pack, character.only=TRUE, quietly=TRUE))
# }

# Source required functions
for (file in sourceFiles) {
    source(file)
}

rm(file, content, installedPackages, m, matched, missingPackages, pack, requiredPackages, sourceFiles)

