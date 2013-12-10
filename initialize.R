#
#     Source all available functions and install requirements
#
# (c) Copyright 2012 J-O Irisson
#     http://creativecommons.org/licenses/by/3.0/
#
#-----------------------------------------------------------------------------


## List source files
#-----------------------------------------------------------------------------

# Functions are in the folder "library", which is next to the current file; so we need to detect where we are
# This does it, for reasons I don't fully understand
# source: http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
if (length(frame_files) == 0) {
	path <- getwd()
} else {
	path <- dirname(frame_files[[length(frame_files)]])
}

# Detect names of source files
sourceFiles = list.files(paste(path, "/library", sep=""), pattern="\\.(R|r)$", full=TRUE)


## Detect and install packages
#-----------------------------------------------------------------------------

# Detect which packages are used in those source files
requiredPackages = c()
for (file in sourceFiles) {
    # read content
    content = scan(file, what="character", quiet=T, sep="\n")
    # extract require or library calls
    m = regexpr("(require|library)\\(.*?\\)", content)
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
library("utils")
# NB: this is for installed.packages(). utils is usually loaded automatically with R but when we exectute the initialize.R script from .Rprofile, it seems this packages is not loaded yet.
installedPackages = row.names(installed.packages())
missingPackages = setdiff(requiredPackages, installedPackages)
if (length(missingPackages) > 0) {
    install.packages(missingPackages, repos="http://cran.at.r-project.org")
}


## Set up computation environment
#-----------------------------------------------------------------------------

# Source required functions
for (file in sourceFiles) {
    source(file)
}
message("-> All functions loaded")


# Check the status of the environment database

# location of environmental data
if (is.null(getOption("atlasr.env.data"))) {
  options(atlasr.env.data=paste(path, "/../env_data", sep=""))
}
# NB: it can be set by the user in a script or in the .Rprofile config file

# check and update if necessary
update.env.data()


# Cleanup the environment
rm(file, content, installedPackages, m, matched, missingPackages, pack, requiredPackages, sourceFiles, frame_files, path)
