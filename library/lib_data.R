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
# (c) Copyright 2011-2012 S Mormede
#                         J-O Irisson
#     http://creativecommons.org/licenses/by/3.0/
#
#-----------------------------------------------------------------------------


## Environment database
#-----------------------------------------------------------------------------

# Ben
# "http://webdav.data.aad.gov.au/data/environmental/derived/antarctic"
# Bruno
# "http://share.biodiversity.aq/GIS/antarctic"

update.env.data <- function(url="http://share.biodiversity.aq/GIS/antarctic", path=getOption("atlasr.env.data")) {
  #
  # Check if the database of environmental data is present and up-to-date
  # If not, download the missing bits
  #
  # url   URL where the database can be downloaded
  # path  where to store the netCDF files of the database
  #

  suppressPackageStartupMessages(require("stringr", quietly=TRUE))  # string manipulation
  suppressPackageStartupMessages(require("tools", quietly=TRUE))    # md5 sum
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))     # automation and loops

  # get remote checksum
  checksumsFile <- tempfile()
  download.file(str_c(url, "/checksums.txt"), destfile=checksumsFile, quiet=T)
  remoteMD5 <- read.table(checksumsFile, header=F, sep=",", col.names=c("path", "md5"), stringsAsFactors=FALSE)
  unlink(checksumsFile)

  # remove asc files, for now
  # TODO remove this restriction at some point
  remoteMD5 <- remoteMD5[!str_detect(remoteMD5$path, "asc\\.zip$"), ]

  # extract directory names
  paths <- str_split(remoteMD5$path, pattern=fixed("/"))
  remoteMD5$dir <- laply(paths, `[`, 1)
  dirs <- unique(remoteMD5$dir)

  # create directories if they don't exist yet
  dirs <- str_c(path, "/", dirs)
  a_ply(dirs, 1, dir.create, showWarnings=FALSE, recursive=TRUE)

  # prepare url and destination
  remoteMD5$url <- str_c(url, "/", remoteMD5$path)
  remoteMD5$destination <- str_c(path, "/", remoteMD5$path)


  # get missing files
  if ( ! file.exists(path) ) {
    # when the database does not exist yet, get everything
    warning("No environment database", call.=FALSE, immediate.=TRUE)
    message("-> Download entire environment database")
    message("   This will be long...")

    # split by directory
    d_ply(remoteMD5, ~dir, function(x) {
      # print directory identifier
      message(str_c("   ", x$dir[1]))

      # download all files, with a global progress bar
      a_ply(x, 1, function(x) {
        download.file(x$url, destfile=x$destination, quiet=TRUE)
        # unzip if it is a zipped netcdf file
        if (str_detect(x$destination, "nc\\.zip$")) {
          unzip(x$destination, exdir=path)
        }
      }, .progress="text")
    })


  } else {
    # when the database exists, compare its content to the online repository

    # get local checksums
    files <- list.files(path, pattern="(zip|png|csv|shp|dbf|shx)$", full=T, recursive=TRUE)
    localMD5 <- data.frame(
      path=str_replace(files, str_c(path, "/"), ""),
      md5=as.character(md5sum(files)),
      stringsAsFactors=FALSE
    )

    # get missing files
    missing <- remoteMD5[ ! remoteMD5$md5 %in% localMD5$md5,]
    additional <- localMD5[ ! localMD5$md5 %in% remoteMD5$md5,]

    if (nrow(missing) > 0 | nrow(additional) > 0) {
      message("-> Updating environment database")

      # delete additional files
      if (nrow(additional) > 0) {
        # issue a message
        many <- ( nrow(additional) > 10 )
        if (many) {
          additionalFiles <- str_c(additional$path[1:10], collapse="\n     ")
          additionalFiles <- str_c(additionalFiles, "\n      and ", nrow(additional) - 10, " more ...")
        } else {
          additionalFiles <- str_c(additional$path, collapse="\n     ")
        }
        message("   The files : \n     ", additionalFiles, "\n   are not on the server.")

        # delete files
        # message("   They will be deleted.")
        # unlink(str_c(path, "/", additional$path))
      }

      # download missing files
      if (nrow(missing) > 0) {
        # issue a message
        many <- ( nrow(missing) > 10 )
        if (many) {
          missingFiles <- str_c(missing$path[1:10], collapse="\n     ")
          missingFiles <- str_c(missingFiles, "\n      and ", nrow(missing) - 10, " more ...")
        } else {
          missingFiles <- str_c(missing$path, collapse="\n     ")
        }
        message("   The files : \n     ", missingFiles, "\n   are missing or have been updated. They will be downloaded.")

        # download the missing files
        a_ply(missing, 1, function(x) {
          download.file(x$url, destfile=x$destination, quiet=many)
          # unzip if it is a zipped netcdf file
          if (str_detect(x$destination, "nc\\.zip$")) {
            unzip(x$destination, exdir=path)
          }
        }, .progress=ifelse(many, "text", "none"))
      }

    } else {
      message("-> Environment database up to date")
    }

  }

  return(invisible(path))
}


list.env.data <- function(variables="", path=getOption("atlasr.env.data"), full=FALSE, ...) {
  #
  # List available environment variables
  #
  # variables   only output variables matching this (matches all when "")
  # path        where to look for netCDF files
  # full        when TRUE, return the full path to the files, otherwise only return the variable name
  # ...         passed to partial.match() (in particular to use quiet)
  #

  suppressPackageStartupMessages(require("stringr", quietly=TRUE))

  if (!file.exists(path)) {
    stop("Environment database not found in ", path)
  }

  # list all files and transform file names into variable names
  ncFiles <- list.files(path, pattern=glob2rx("*.nc"), full.names=TRUE)

  ncVariables <- str_replace(ncFiles, str_c(path, "/"), "")
  ncVariables <- str_replace(ncVariables, fixed("_0.1_0.1"), "")
  # NB: this string does not appear in most recent version of env filenames, but keep for backwards compatibility
  ncVariables <- str_replace(ncVariables, fixed(".nc"), "")

  # possibly match and expand variable names
  ncVariablesMatched <- partial.match(variables, ncVariables, ...)

  # get corresponding file names (in the same order!)
  ncFilesMatched <- ncFiles[match(ncVariablesMatched, ncVariables)]

  if (full) {
    return(ncFilesMatched)
  } else {
    return(ncVariablesMatched)
  }
}

read.env.data <- function(variables="", path=getOption("atlasr.env.data"), verbose=FALSE) {
  # Read data from the netCDF files
  #
  # variables   only output variables matching this (matches all when "")
  # path        path to the location of netCDF files
  #
  # NB: previously we saved the database in a RData file, but it is actually quicker to read it directly from the netCDF files

  suppressPackageStartupMessages(require("ncdf"))
  suppressPackageStartupMessages(require("plyr"))

  if ( ! file.exists(path) ) {
    stop("Environment database not found in ", path)
  }

  if ( verbose ) message("Read environment data")

  # select which netCDF files to read
  ncFiles = list.env.data(variables, path, full=T, verbose=verbose)
  ncVariables = list.env.data(variables, path, verbose=FALSE)
  # NB: possibly inform about variable names expansion only for the first pass

  # read data inside each file
  database <- llply(ncFiles, function(ncFile) {

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

  # # remove plyr attributes (useless here)
  # attributes(database) <- list()

  # name elements of the list
  # NB: sometimes the data names in the netCDF files are the same across several files, we will instead use a name derived from the file
  names(database) <- ncVariables

  # treat geomorphology as a factor
  if ( ! is.null(database$geomorphology) ) {
    z <- factor(database$geomorphology$z)
    dim(z) <- dim(database$geomorphology$z)
    database$geomorphology$z <- z
  }

  return(database)
}


mask.env.data <- function(database, ...) {
  #
  # Mask points which are on land in the environmental data
  #
  # database  a list resulting from read.env.data
  #

  suppressPackageStartupMessages(require("plyr", quietly=TRUE))

  # read bathymetry
  suppressMessages(bathy <- read.env.data(variables="bathymetry", ...)[[1]])

  # mask points on land (replace by NA)
  onLand <- bathy$z >= 0
  database <- llply(database, function(x, mask) {
    x$z[mask] <- NA
    return(x)
  }, mask=onLand)

  return(database)
}


get.env.data <- function(lon, lat, database, verbose=FALSE) {
  # Associate the environmental data from `database` to the points specified in `dataset`
  #
  # dataset     data.frame with at least lon and lat columns
  # database    list resulting from read.env.data
  # lon/lat.name  names (or numbers) of the lon and lat columns in the dataset
  #

  if ( verbose ) message("Fetch environment data for observed points")

  # bilinear interpolation
  # NB: if we extract points that are exactly on the grid, we do not interpolate anything and this method is quite fast
  # TODO test if we can be faster by skipping the interpolation if we are always exactly on the grid
  suppressPackageStartupMessages(require("fields"))

  # rename lon and lat to ensure consistency in the following
  coord <- data.frame(lon, lat)

  # extract/interpolate the environmental data at the locations specified in the dataset
  envData <- as.data.frame(lapply(database, function(x) {
    if (class(x$z[1,1]) == "numeric") {
      interp.surface(x, coord)
    } else {
      interp.nn(x, coord)
    }
  }))

  return(envData)
}

closest.index <- function(x, y) {
  #
  #	Find the indexes of x such as the corresponding elements are closest to the values in y
  #
	round(approx(x,1:length(x),y)$y)
}

interp.nn <- function(x, coord) {
  #
  # Nearest neigbour interpolation
  #
  # x     list with components, x, y and z
  # coord x,y coordinates to interpolate to
  #

  # find the closest grid points in x and y
  xIdx <- closest.index(x$x, coord[,1])
  yIdx <- closest.index(x$y, coord[,2])

  # fetch the value at those grid points
  idx <- cbind(xIdx, yIdx)
  apply(idx, 1, function(id, z) { z[id[1],id[2]] }, z=x$z)
}


## Observation data
#-----------------------------------------------------------------------------

read.data <- function(file, verbose=FALSE, ...) {
  #
  # Detect the format of a data file and read it
  #
  # file  path to the data file
  # ...   passed to the method used for reading the file: gdata::read.xls, read.csv, read.txt
  #

  suppressPackageStartupMessages(library("tools", quietly=TRUE))

  # checks
  if ( length(file) > 1 ) {
    warning("Can only read only one file at a time. Using the first element")
    file <- file[1]
  }

  # ensure windows/unix interoperability and more
  file <- clean.path(file)

  # check file existence
  if ( ! file.exists(file) ) {
    stop("Cannot find file: ", file)
  }


  if ( verbose ) {
    fileName <- basename(file)
    message("Read data in ", fileName)
  }

  # get file extension
  ext <- file_ext(file)

  # choose reading method according to extension
  if (ext == "xls") {
    suppressPackageStartupMessages(require("gdata"))
    d <- read.xls(xls=file, ...)
  } else if (ext == "csv") {
    d <- read.csv(file=file, ...)
  } else if (ext == "txt") {
    d <- read.table(file=file, header=TRUE, ...)
  } else {
    stop("Unknown file type")
  }

  # give some flexibility in naming lat and lon
  names(d)[tolower(names(d)) %in% c("latitude","lat")] = "lat"
  names(d)[tolower(names(d)) %in% c("longitude","lon","long")] = "lon"

  # check that lat and lon are present
  if (! all(c("lat","lon") %in% names(d)) ) {
    stop("Need latitude and longitude in the input data")
  }

  # check consistency of lat and lon
  bizarreCoord <- which(d$lat < -90 | d$lat > 90 | d$lon < -180 | d$lon > 180)
  if ( length(bizarreCoord) > 0 ) {
    warning("Some coordinates look funny in your input file. Please check you data")
    print(x[bizarreCoord, c("lon", "lat")])
  }

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

partial.match <- function(pattern, choices, verbose=FALSE) {
  # Match abbreviated or partial variable names
  #
  # pattern pattern to match in choices
  # choices list of matching possibilities
  # quiet   wether to provide a message about matches
  #

  # prepare storage
  res = c()
  var.index = c()
  choice.index = c()

  for (i in seq(along=pattern)) {
    # try exact matches first
    idx = grep(paste("^",pattern[i],"$",sep=""), choices)
    if (length(idx) == 0) {
      # then try partial matches if needed
      idx = grep(pattern[i], choices)

      # if no variable can be matched, issue a warning
      if (length(idx) == 0) {
        nChoices <- length(choices)
        if ( nChoices > 6 ) {
          possibilities <- paste(paste(choices[1:6], collapse="\n    "), "\n    and", nChoices - 6, "more...")
        } else {
          possibilities <- paste(choices, collapse="\n    ")
        }
        stop("No name matching \"", pattern[i], "\" could be found.\n  The possibilities were:\n    ", possibilities)
      }
    }
    # store all matches for all variables
    matches = choices[idx]
    res = c(res, matches)
    choice.index = c(choice.index, idx)
    var.index = c(var.index, rep(i, length(matches)))

    # issue a message when an expansion match occurred
    if ( verbose & any(matches != pattern[i])) {
      messageText = paste(pattern[i], " expanded to ", sep="")
      # compute the amount of padding to get a nicely aligned list
      padding = paste(rep(" ", times=nchar(messageText)), collapse="")
      # inform about the expansion
      message(messageText, paste(matches, collapse=paste("\n", padding, sep="")))
    }
  }

  # store the indexes of the matches as attributes
  attr(res, "var.index") <- var.index
  attr(res, "choice.index") <- choice.index

  return(res)
}


weight.data <- function(x, weights, warn=TRUE) {
  #
  # Apply weights to the columns of a data.frame
  #
  # x         data.frame to be weighted
  # weights   named vector of weights
  #           names will be expanded using partial.match and the final set of names must match the names of x
  #           numeric values will be scaled to a maximum of 1
  # warn      wether to warn when weights are missing for some of the columns
  #

  # expand each element of weights by name
  # this allows to specify something like weights=c(nox=2) and have *all* nox variables double weighted
  expandedNames <- partial.match(names(weights), names(x), quiet=TRUE)
  # NB: we expand based in what is in x. If the name of an element in weights does not match, partial.match will throw an error

  # if name expansion resulted in more weights, replicate the weight values appropriately
  if (length(expandedNames) > length(weights)) {
    weights <- weights[attr(expandedNames, "var.index")]
    # NB: we use the attributes regarding indexes assigned by match.var to do this
  }

  # assign the new, expanded names
  names(weights) <- expandedNames

  # warn when all elements of x do not have a weight associated with them and use the default weight (1)
  weights <- weights[names(x)]
  NAweights <- is.na(weights)
  if ( any(NAweights) ) {
    if (warn) warning("Weights are missing for ", paste(names(x)[NAweights], collapse=", ") ,"\n  Assuming weight(s) of 1")
    weights[NAweights] <- 1
    names(weights) <- names(x)
  }

  # normalize so that max weight is 1
  weights <- weights / max(weights)

  # apply weights to the columns of x
  for (i in 1:ncol(x)) {
    x[,i] <- x[,i] * weights[i]
  }

  # store the weights as attributes
  attr(x, "weights") <- weights

  return(x)
}

function.maker <- function(str) {
    #
    # Transform a character string into a function
    #
    # str   character string defining the function
    #

    # first consider the situation in which the string is just the name of a function ("log", "sqrt")
    f <- tryCatch(
      eval(parse(text=str)),
      # if the string does not represent a function, evaluating it can throw errors or warnings. We catch them and forget about them
      warning=function(x) { NULL },
      error=function(x) { NULL }
    )

    if ( ! is.function(f) ) {
      # when it is not, assume it is the body of a function definition in which the argument is x

      # start with an empty function
      f <- function(x) {}
      environment(f) <- baseenv()

      # make sure str is enclosed in curly brackets in case the function has several commands (i.e. spans several lines)
      str <- str_c("{",str,"}")

      # fill in the body of the function
      body(f) <- tryCatch(
        parse(text=str),
        error = function(e) { stop("Error parsing transformation function : ", str, call.=FALSE) }
      )
    }

    return(f)
}

transform.data <- function(x, transformations, warn=TRUE) {
  #
  # Apply transformations to the columns of a data.frame
  #
  # x                 data.frame to be transformed
  # transformations   named vector of transformations
  #                   names will be expanded using partial.match and the final set of names must match the names of x
  # warn              wether to warn when transformations are missing for some of the columns
  #

  # NB: see weight.data, which is very similar, for more detailed comments

  # expand each element of transformations by name
  expandedNames <- partial.match(names(transformations), names(x), quiet=TRUE)

  # replicate transformations appropriately
  if (length(expandedNames) > length(transformations)) {
    transformations <- transformations[attr(expandedNames, "var.index")]
  }

  # assign the new, expanded names
  names(transformations) <- expandedNames

  # warn when all elements of x do not have a weight associated with them and use the default weight (1)
  transformations <- transformations[names(x)]
  NAtransformations <- is.na(transformations)
  if ( any(NAtransformations) ) {
    if (warn) warning("Transformations are missing for ", paste(names(x)[NAtransformations], collapse=", ") ,"\n  Assuming no transformation")
    transformations[NAtransformations] <- "x"
    names(transformations) <- names(x)
    # TODO just skip these columns for improved speed
  }

  # apply transformations to the columns of x
  for (i in 1:ncol(x)) {

    # make a function from the textual description
    trans <- function.maker(transformations[i])

    # apply the function, catching errors and warnings
    x[,i] <- tryCatch(
      # try applying the function
      trans(x[,i]),
      # if everything goes smoothly, the result is stored in x[,i]

      # if there is a warning, print it and still return the value
      warning=function(w) {
        # custom message
        warning("Warning when applying transformation \"", transformations[i], "\" to \"", names(x)[i], "\" : \n  ", w$message, call.=FALSE)
        # # original message
        # warning(w)
        # still compute the result, suppressing the warning
        y <- suppressWarnings(trans(x[,i]))
        return(y)
      },

      # if there is an error, convert it into a warning and return NA
      error=function(e) {
        warning("Error when applying transformation \"", transformations[i], "\" to \"", names(x)[i], "\" : \n  ", e$message, "\nReturning NA for \"", names(x)[i], "\"", call.=FALSE)
        return(NA)
      }
    )

  }

  # store the transformations as attributes
  attr(x, "transformations") <- transformations

  return(x)
}


rasterize <- function(x, vars, n=10, precisions=NULL, fun=sum, ...) {
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

  suppressPackageStartupMessages(require("plyr", quietly=TRUE))

  # checks
  OKvars <- vars %in% names(x)
  if (! all(OKvars) ) {
    stop("Variable(s) ", vars[!OKvars], " not in ", deparse(substitute(x)))
  }

  # if the precisions are not specified, use n
  if (is.null(precisions)) {
    precisions <- sapply(lapply(x[vars], range, na.rm=T), diff) / n
  }
  # otherwise, check that it is correctly specified
  else {
    if (length(precisions) != length(vars)) {
      stop("The vector of precisions does not have as many elements as variables in vars")
    }
  }

  message("-> Bin data with precisions : ", paste(vars, round(precisions, 3), sep="=", collapse=" x "))

  # round columns to the given precision
  suppressPackageStartupMessages(require("plyr"))
  for (j in seq(along=vars)) {
    x[,vars[j]] <- plyr::round_any(x[,vars[j]], accuracy=precisions[j])
  }

  # when there are variables in addition to the binned ones, summarize their information per bin
  if (ncol(x) > length(vars)) {
    x <- ddply(x, .variables=vars, function(X, .vars, ...) {
      actualData <- X[! names(X) %in% .vars]
      summarizedData <- data.frame(t(sapply(actualData, fun, ...)))
      summarizedData$freq <- nrow(X)
      return(summarizedData)
    }, .vars=vars, ...)
  }
  # otherwise, make sure that the combinations of bins are specified only once and count the number of observations per bin
  else {
    x <- plyr::count(x)
  }

  return(x)
}

rasterise <- rasterize

weight.per.bin <- function(lat, lon, bin.size) {

   suppressPackageStartupMessages(library("plyr", quietly=TRUE))

   # round lat and lon on the binnin grid
   lon <- round_any(lon, bin.size)
   lat <- round_any(lat, bin.size)

   # compute number of data points per bin
   x <- data.frame(lon, lat)
   nb <- count(x)

   # assign that number to each original point
   nb <- join(x, nb, by=c("lon", "lat"))

   # compute weights
   w <- 1 / nb$freq

   return(w)
}


clean.path <- function(file, ...) {
  #
  # Clean file paths
  #
  # file    file path to convert

  file <- normalizePath(file, winslash="/", mustWork=FALSE)

  return(file)
}

too.many.na <- function(x, p, weights=NULL) {
   # Identify lines with a missing value proportion above p
   #
   # x         a matrix or data.frame
   # p         the proportion of missing values
   # weights   weights for the columns of x, which weight the missing values in them

   # compute relative weights
   if (is.null(weights)) {
      weights <- rep(1, ncol(x))
   }
   weights <- weights / sum(weights)

   # find NAs
   na <- is.na(x)

   # count NAs per line, weighted
   props <- apply(na, 1, function(y, w) {
      sum(y * w)
   }, w=weights)

   return(props >= p)
}


write.shapefile <- function(x, name, variables=NULL) {
  #
  # Write shapefiles from a data.frame
  #
  # x           data.frame or object coercible as such
  #             should contain columns lat and lon + variables of interest
  # name        base name of the shapefile
  # variables   names or column numbers of the variables of interest
  #             when NULL, use all columns but lon and lat
  #

  suppressPackageStartupMessages(require("maptools", quietly=TRUE))

  # make sure x is a data.frame
  x <- as.data.frame(x)
  # NB: SpatialPointsDataFrame only accepts pure data.frame objects

  # determine variables
  if (is.null(variables)) {
    variables <- setdiff(names(x), c("lon", "lat"))
  }

  # convert to spatial object and write shapefile
  xSp <- SpatialPointsDataFrame(x[c("lon", "lat")], x[variables], proj4string=CRS("+proj=longlat +datum=WGS84"))
  writeSpatialShape(xSp, name)

  # add the .prj file
  cat("GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]]\n", file=paste(name, ".prj", sep=""))

}

write.raster <- function(x, name, variables=NULL) {
   library("sp")

   # make sure x is a data.frame
   x <- as.data.frame(x)
   # NB: SpatialPointsDataFrame only accepts pure data.frame objects

   # determine which variables to keep
   if (is.null(variables)) {
     variables <- setdiff(names(x), c("lon", "lat"))
   }

   # create a spatial object
   x <- x[x$lon > -180,]
   # xSp <- SpatialPixelsDataFrame(points=x[,c("lon", "lat")], data=x[,variables])
   xSp <- SpatialPixelsDataFrame(points=x[,c("lon", "lat")], data=x[variables], proj4string=CRS("+proj=longlat +datum=WGS84"))

   # convert into raster
   library("raster")
   xR <- stack(xSp)

   # write it in a raster file
   writeRaster(xR, filename=name, format="raster")
   library("ncdf")
   writeRaster(xR, filename=name, format="CDF", overwrite=T)
   writeAsciiGrid(xSp, fname=paste(name, ".ag", sep=""))
   # writeRaster(xR1, filename=name, format="ascii")
   # writeRaster(xR, filename=name, format="EHdr")
   # writeRaster(xR, filename=name, format="GTiff")
   # writeRaster(xR, filename=name, format="netCDF")

}


}
