#
#   Deal with data
#
#   - download the database if it does not exists
#   - inspect the netCDF database
#   - read netCDF files
#   - read observation data
#   - create a grid for the prediction space
#   - collate environmental data to a given dataset
#
# (c) Copyright 2011-2012 S Mormede, J-O Irisson
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------


## Environment database
#-----------------------------------------------------------------------------

dl.env.data <- function(url="http://dl.dropbox.com/u/1047321/antarctic_data_archive.zip", path="env_data", ...) {
  # Download the environmental database
  # = get an archive of all netCDF files and decompress it
  #
  # url   URL where the database can be downloaded
  # path  where to store the netCDF files of the database
  #
  
  message("-> Download environment database")
  
  temp <- tempfile() # download the zip as a temporary file
  download.file(url, destfile=temp)
  unzip(temp, exdir=path)
  unlink(temp)
  
  return(invisible(path))
}


check.get.env.data <- function(path="env_data") {
  # Check if the database of environmental data is present; if not, download it
  #
  # path  where to store the netCDF files of the database
  #
  
  if (! file.exists(path)) {
    warning("No environment database", call.=FALSE, immediate.=TRUE)
    dl.env.data(path=path)
  }
  return(invisible(path))
}


list.env.data <- function(variables="", path="env_data", full=FALSE, ...) {
  #
  # List available environment variables
  #
  # variables   only output variables matching this (matches all when "")
  # path        where to look for netCDF files
  # full        when TRUE, return the full path to the files, otherwise only return the variable name
  # ...         passed to match.vars() (in particular to use quiet)
  #

  suppressPackageStartupMessages(require("stringr", quietly=TRUE))

  # check that environment data exists
  check.get.env.data(path=path)

  # list all files and transform file names into variable names
  ncFiles <- list.files(path, pattern=glob2rx("*.nc"), full.names=TRUE)

  ncVariables <- str_replace(ncFiles, paste(path, "/", sep=""), "")
  ncVariables <- str_replace(ncVariables, "_0.1_0.1.nc", "")

  # possibly match and expand variable names
  ncVariablesMatched <- match.vars(variables, ncVariables, ...)
  ncFilesMatched <- ncFiles[ncVariables %in% ncVariablesMatched]

  if (full) {
    return(ncFilesMatched)
  } else {
    return(ncVariablesMatched)
  }
}


read.env.data <- function(variables="", path="env_data", ...) {
  # Read data from the netCDF files
  #
  # variables   only output variables matching this (matches all when "")
  # path        path to the location of netCDF files
  # ...         passed to list.env.data()
  #
  # NB: previously we saved the database in a RData file, but it is actually quicker to read it directly from the netCDF files

  # functions to access netCDF files
  suppressPackageStartupMessages(require("ncdf"))
  # functions to automate loops and split-apply functions
  suppressPackageStartupMessages(require("plyr"))

  # check that environment data exists
  check.get.env.data(path=path)

  message("-> Read environment data")

  # select which netCDF files to read
  ncFiles = list.env.data(variables, path, full=T, quiet=FALSE)
  ncVariables = list.env.data(variables, path, quiet=TRUE)
  # NB: inform about variable names expansion only for the first pass

  # read data inside each file
  database <- alply(ncFiles, 1, function(ncFile) {

    # open netCDF file
    ncid <- open.ncdf(ncFile)

    # get the values of the lon and lat
    lon <- get.var.ncdf(ncid, varid="lon")
    lat <- get.var.ncdf(ncid, varid="lat")

    # extract the variable name
    varName <- names(ncid$var)
    # usually there is only one variable per netCDF file; check this
    # if (length(varName) > 1) {
    #   stop("Several variables in dataset", ncFile)
    # }

    # now retrieve the actual data
    dat = get.var.ncdf(ncid, varName)
    # NB: get the full data for now, we will subsample it afterwards rather than at reading time

    # get the missing value code
    missingValue <- att.get.ncdf(ncid, varName, "_FillValue")
    # replace any such data with NA
    dat[which(dat==missingValue$value)] <- NA

    # close the file now that we're finished with it
    close(ncid)

    # store as an xy-coords list (suitable for persp, image, contour, etc.)
    dat <- list(x=lon, y=lat, z=dat)

    # display data as a reality check
    # image(dat, main=varName)
  
    return(dat)
  }, .progress=ifelse(length(ncFiles) > 5,"text","none"))  # get a nice progress bar when there are several files to read

  # remove plyr attributes (useless here)
  attributes(database) <- list()

  # name elements of the list
  # NB: sometimes the data names in the netCDF files are the same across several files, we will instead use a name derived from the file
  names(database) <- ncVariables
  
  return(database)
}


associate.env.data <- function(dataset, database, lon.name="lon", lat.name="lat") {
  # Associate the environmental data from `database` to the points specified in `dataset`
  #
  # dataset     data.frame with at least lon and lat columns
  # database    list resulting from read.env.data
  # lon/lat.name  names (or numbers) of the lon and lat columns in the dataset
  #

  message("-> Fetch environment data for points in ", deparse(substitute(dataset)))

  # bilinear interpolation
  # NB: if we extract points that are exactly on the grid, we do not interpolate anything and this method is quite fast
  # TODO test if we can be faster by skipping the interpolation if we are always exactly on the grid
  suppressPackageStartupMessages(require("fields"))
  # contains function rename
  suppressPackageStartupMessages(require("plyr"))   

  # rename lon and lat to ensure consistency in the following
  names(dataset)[which(names(dataset) == lon.name)] <- "lon"
  names(dataset)[which(names(dataset) == lat.name)] <- "lat"
  # extract lon and lat from the dataset
  coord <- dataset[,c("lon","lat")]
  
  # extract/interpolate the environmental data at the locations specified in the dataset
  envData <- as.data.frame(lapply(database, function(x) {
    interp.surface(x, coord)
  }))
  
  # store everthing in the original dataset
  dataset <- cbind(dataset, envData)
  
  # remove points that are on land
  dataset <- dataset[dataset$bathymetry < 0,]

  return(dataset)
}


## Observation data
#-----------------------------------------------------------------------------

read.data <- function(file, ...) {
  #
  # Detect the format of a data file and read it
  #
  # file  path to the data file
  # ...   passed to the method used for reading the file: gdata::read.xls, read.csv, read.txt
  #
  
  if (length(file) > 1) {
    stop("Read only one file at a time")
  }
  
  file <- win2unix(file, ...)
  
  fileName <- basename(file)
  message("-> Read data in ", file)
  
  # get file extension
  extension <- strsplit(fileName, split=".", fixed=T)[[1]]
  extension <- extension[length(extension)]
  
  if (extension == "xls") {
    suppressPackageStartupMessages(require("gdata"))
    d <- read.xls(xls=file, ...)
  } else if (extension == "csv") {
    d <- read.csv(file=file, ...)
  } else if (extension == "txt") {
    d <- read.table(file=file, header=TRUE, ...)
  } else {
    stop("Unknown file type")
  }
  
  # give some flexibility in naming lat and lon
  names(d)[tolower(names(d)) %in% c("latitude","lat")] = "lat"
  names(d)[tolower(names(d)) %in% c("longitude","lon","long")] = "lon"

  return(d)
}


## Prediction data
#-----------------------------------------------------------------------------

build.grid <- function(lat.min=-80, lat.max=-30, lat.step=0.1, lon.min=-180, lon.max=180, lon.step=0.5, ...) {
  # Define the grid on which we extract the environmental data
  #
  # lat/lon.lim   vector with the minimum and maximum coordinates
  # lat/lon.step  step in degrees
  #
  
  # compute coordinates
  lat = seq(from=lat.min, to=lat.max, by=lat.step)
  lon = seq(from=lon.min, to=lon.max, by=lon.step)

  # compute the full grid
  return(expand.grid(lon=lon, lat=lat))
}


## Utility functions
#-----------------------------------------------------------------------------

match.vars <- function(vars, choices, quiet=TRUE) {
  # Match abbreviated or partial variable names
  #
  # vars    variable names to match
  # choices list of matching possibilities
  # quiet   wether to provide a message about matches
  #
  
  res = c()
  for (i in seq(along=vars)) {
    # try exact matches first
    matches = grep(paste("^",vars[i],"$",sep=""), choices, value=TRUE)
    if (length(matches) == 0) {
      # then try partial matches if needed
      matches = grep(vars[i], choices, value=TRUE)
    }
    # store all matches for all variables
    res = c(res, matches)

    # issue a message when an expansion match occurred
    if (! quiet & any(matches != vars[i])) {
      messageText = paste("   ", vars[i], " expanded to ", sep="")
      # compute the amount of padding to get a nicely aligned list
      padding = paste(rep(" ", times=nchar(messageText)), collapse="")
      # inform about the expansion
      message(messageText, paste(matches, collapse=paste("\n", padding, sep="")))
    }
  }

  return(res)
}


rasterize <- function(x, vars, n=10, precisions=sapply(lapply(x[vars],range),diff)/n, fun=sum, ...) {
  #
  # Reduce the precision of certain columns of a data.frame to bins and summarize the rest of the information per bin
  # This is a bit like reducing the precise information on locations in a 2D plane to pixels of a given grey level, hence the name
  #
  # x           original data.frame
  # vars        columns to bin
  # n           scalar, number of bins
  # precisions  the precisions used to cut each column in vars
  #             (overrides n)
  # fun         function used to perform the summary
  # ...         further arguments to `fun`
  #


  # checks
  if (length(precisions) != length(vars)) {
    stop("The vector of precisions does not have as many elements as variables in vars")
  }

  message("-> Bin data with precisions : ", paste(vars, precisions, sep="=", collapse=" x "))

  # round columns to the given precision
  suppressPackageStartupMessages(require("plyr"))
  for (j in seq(along=vars)) {
    x[,vars[j]] <- plyr::round_any(x[,vars[j]], accuracy=precisions[j])
  }

  # summarize information per bin
  x <- ddply(x, .variables=vars, function(X, .vars, ...) {
    actualData <- X[! names(X) %in% .vars]
    sapply(actualData, fun, ...)
  }, .vars=vars, ...)

  return(x)
}


win2unix <- function(file, ...) {
  #
  # Convert between windows and unix file paths
  #
  # file    file path to convert
  if (grepl("\\", file, fixed=T)) {
    file <- gsub("\\", "/", file, fixed=T)
  }
  return(file)
}