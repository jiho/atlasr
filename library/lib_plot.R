#
#     Plot data viewed from the South Pole
#
# The plot uses ggplot2 rather than base graphics, because it eases
# the projection and mapping of colours.
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------


## Toolbox functions
#-----------------------------------------------------------------------------

polar_proj <- function(projection="stereographic", orientation=c(-90,0,0)) {
  #
  # Easy access to a suitable polar projection
  # NB: view from south pole (-90)
  #
  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
  suppressPackageStartupMessages(require("mapproj", quietly=TRUE))
  c <- coord_map(projection=projection, orientation=orientation)
  return(c)
}

discrete.colourmap <- function(n=10) {
  #
  # Define some qualitatively different colours
  #
  # n     number of colours on the scale
  #

  # base colour map from http://colorbrewer2.org/
  suppressPackageStartupMessages(require("RColorBrewer", quietly=TRUE))
  cmap <- brewer.pal(12, "Set3")    # NB: 12 is the maximum

  # remove the grey, which is confusing the the following (results in white after saturation change, confusing compared to the default background in ggplot, etc.)
  cmap <- cmap[cmap != "#D9D9D9"]

  # prepare saturated and under-saturated versions of the colors for use when n > 122
  cmapHSV <- rgb2hsv(col2rgb(cmap))

  cmapSat <- cmapHSV
  # saturate
  cmapSat["s",] <- cmapSat["s",] + 0.3
  # ensure saturation is less than 1
  cmapSat["s",][cmapSat["s",] > 1] <- 1
  # convert to colors
  cmapSat <- hsv(cmapSat[1,], cmapSat[2,], cmapSat[3,])

  cmapUndersat <- cmapHSV
  # de-saturate
  cmapUndersat["s",] <- cmapUndersat["s",] - 0.15
  # ensure saturation is less than 1
  cmapUndersat["s",][cmapUndersat["s",] < 0] <- 0
  # convert to colors
  cmapUndersat <- hsv(cmapUndersat[1,], cmapUndersat[2,], cmapUndersat[3,])

  cmap <- c(cmap, cmapSat, cmapUndersat)
  # barplot(rep(1, times=36), col=cmap)

  # extract colours
  if (n > length(cmap)) {
    warning("Not enough colours to plot everything. It is unlikely that you will be able to discriminate between more than 33 colors on the plot anyway.")
  }
  cmap <- cmap[1:n]

  return(cmap)
}

continuous.colourmap <- function(n=10) {
  #
  # Define a sequence of colours on a gradient
  #

  # get nice colours
  suppressPackageStartupMessages(require("RColorBrewer", quietly=TRUE))
  colours <- brewer.pal(n=6, name="Spectral")

  # reverse the colours to get the colder first (and map those to the smaller values)
  colours <- rev(colours)

  # interpolate them
  colorRampPalette(colors=colours, space="Lab")(n)

}

plot.pred <- function(x, ...) {
  #
  # Generic for the plot of predictions from a model
  #
  UseMethod("plot.pred")
}

clip.polygon <- function(x, lon.min=-180, lon.max=180, lat.min=-90, lat.max=90) {
  #
  # Clip a polygon within given limits
  #
  # x       a data.frame with coordinates (as the first two columns, in the order lon then lat) defining a polygon, or several polygons separated by NAs
  # lon.*
  # lat.*   limits within which the polygon is cut
  #

  suppressPackageStartupMessages(require("gpclib", quietly=TRUE))
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))

  # reformat the data.frame x to be converted into polygons
  x <- x[,1:2]
  names(x) <- c("x", "y")

  # cut all polygons
  x$id <- cumsum(is.na(x$x))
  x <- na.omit(x)
  xL <- split(x[,1:2], f=x$id)

  # remove completely out of range polygons
  # this accelerates the next steps
  inRange <- laply(xL, function(x) {
    any(x$x >= lon.min & x$x <= lon.max & x$y >= lat.min & x$y <= lat.max)
  })
  xL <- xL[inRange]

  # convert each piece into a polygon object for gpclib
  xP <- as(xL[[1]], "gpc.poly")
  if (length(xL) > 1) {
    for (i in 2:length(xL)) {
      tempP <- as(na.omit(xL[[i]]), "gpc.poly")
      xP <- append.poly(xP, tempP)
    }
  }

  # prepare the mask
  mask <- data.frame(x=c(lon.min, lon.min, lon.max, lon.max), y=c(lat.min, lat.max, lat.max, lat.min))
  maskP <- as(mask, "gpc.poly")

  # compute clipping
  clippedP <- intersect(xP, maskP)

  # convert into a data.frame
  clipped <- clippedP@pts
  clipped <- ldply(clipped, function(x) {
    # extract elements
    X <- data.frame(lon=x$x, lat=x$y, hole=x$hole)
    # close polygon
    X <- rbind(X, X[1,])
    # separate from next polygon
    X <- rbind(X, NA)
    return(X)
  })

  return(clipped)
}

clip.to.data <- function(x, data, expand=1) {
  #
  # Clip a polygon within limits of the data argument
  #
  # x       a data.frame with coordinates (lon and lat) defining a polygon, or several polygons separated by NAs
  #         typically, a coastline
  # data    a data.frame with coordinates (lon and lat) where the data points are and to which the
  # expand  how much to expand x layer around the data (in lon/lat units)
  #

  # get data range
  lons <- range(data$lon, na.rm=TRUE) + c(-expand, +expand)
  lats <- range(data$lat, na.rm=TRUE) + c(-expand, +expand)

  # select portions of land which fit this data
  clipped <- clip.polygon(x[,c("lon", "lat")], lons[1], lons[2], lats[1], lats[2])

  return(clipped)
}


## Plotting functions
#-----------------------------------------------------------------------------

layer_land <- function(x, expand=1, path=getOption("atlasr.env.data"), ...) {
  #
  # Produce a ggplot layer with coastlines corresponding to a dataset
  #
  # x       dataset of interest, with columns lon and lat at least
  #         the boundaries of the land layer are taken from this
  # expand  how much to expand the land layer around the data (in lon/lat units)
  # path    path to environmental database
  #

  # read coordinates of land masses
  land <- read.csv(str_c(path, "/worldmap-below_30-rough-no_countries.csv"))

  # get data range
  lats <- range(x$lat, na.rm=TRUE) + c(-expand, +expand)
  lons <- range(x$lon, na.rm=TRUE) + c(-expand, +expand)

  # select portions of land which fit this data
  land <- land[land$lat >= lats[1] & land$lat <= lats[2] & land$lon >= lons[1] & land$lon <= lons[2],]
  # TODO use cut polygon in one of the mapping libraries

  landLayer <- geom_path(aes(x=lon, y=lat), size=0.3, data=land, na.rm=T)

  return(landLayer)
}


polar.ggplot <- function(data, mapping=aes(), geom=c("auto", "point", "tile"), lat.precision=NULL, lon.precision=NULL, draw.coast=TRUE, scale=1, ...) {
  #
  # data          data frame with columns lat, lon, and variables to plot
  # mapping       a call to `aes()` which maps a variable to a plotting
  #               aesthetic characteristic (fill, size, alpha, etc.)
  # geom          the type of plot ("geometry" in ggplot parlance) to produce
  #               = points or tiles (possibly better looking, longer to plot)
  #               auto chooses points or tiles depending on the size of the data
  # lat.precision
  # lon.precision the precision at which lat and lon are considered
  #               (in degrees). If they are larger than the original
  #               precision, the data is *subsampled* to those locations
  #               i.e. some data is actually dropped. If you want to average or sum
  #               the data per cell, use rasterize() in lib_data.R
  # draw.coast    wether to draw a basic coastline
  # scale         scale of the points plotted
  # ...           passed to the appropriate geom
  #

  suppressPackageStartupMessages(require("plyr", quietly=TRUE))
  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
  suppressPackageStartupMessages(require("stringr", quietly=TRUE))

  # Check arguments
  # geoms
  geom <- match.arg(geom)

  # allow lon/lat to be called more liberally
  names(data)[tolower(names(data)) %in% c("latitude","lat")] <- "lat"
  names(data)[tolower(names(data)) %in% c("longitude","lon","long")] <- "lon"
  # check that we have something that looks like lat and lon
  if (! all(c("lat","lon") %in% names(data)) ) {
    stop("Need two columns named lat and lon to be able to plot\nYou have ", paste(names(data), collapse=", "))
  }

  # if new precisions are specified for lat or lon, subsample the data
  if (!is.null(lat.precision)) {
    # compute the vector of latitudes
    lats <- sort(unique(data$lat))
    # recompute the precision to be the closest multiple of the current step
    step <- unique(round(diff(lats),5))
    # NB: use round to solve floating point precision issues
    if (lat.precision > step) {
      # proceed only if the given precision is actually larger than the current step (otherwise there is nothing to resample)
      lat.precision <- round(lat.precision/step)*step
      # reduce to the given precision
      lats <- seq(from=min(lats), to=max(lats), by=lat.precision)
      # select points at those latitudes only
      data <- data[round(data$lat,5) %in% round(lats,5),]
      # NB: use round to solve floating point precision issues
    }
  }
  if (!is.null(lon.precision)) {
    lons <- sort(unique(data$lon))
    step <- unique(round(diff(lons),5))
    if (lon.precision > step) {
      lon.precision <- round(lon.precision/step)*step
      lons <- seq(from=min(lons), to=max(lons), by=lon.precision)
      data <- data[round(data$lon,5) %in% round(lons,5),]
    }
  }

  # define the geom depending on the size of the data when geom="auto"
  if (geom == "auto") {
    if (nrow(data) <= 2000) {
      geom <- "tile"
    } else {
      geom <- "point"
    }
  }

  # remove colour mapping (one must use fill because we only use filled shapes)
  mapping <- mapping[names(mapping) != "colour",]

  # Plot
  # prepare plot
  p <- ggplot(data, aes(x=lon, y=lat)) +
        # stereographic projection
        polar_proj()

  # plot points or tiles depending on the geom argument
  if (geom == "point") {
    # add mapping of size to better cover the space (smaller points near the center)
    mapping = c(aes(size=lat), mapping)
    class(mapping) = "uneval"

    # try to guess the size of points based on how many data there is
    # (the more data, the smaller the points)
    nLats <- length(unique(data$lat))
    baseSize <- scale * 90 / nLats

    # plot
    # NB: shape: 21 = filled point, 22 = filled square, 23 = filled losange
    p <- p + geom_point(mapping=mapping, shape=21, colour=NA, ...) + scale_size(range=c(baseSize, baseSize*2.2), guide=FALSE)
  } else if (geom == "tile"){
    p <- p + geom_tile(mapping=mapping, ...)
  }

  # Get and re-cut coastline if need be
  if (draw.coast) {
    # extract the whole world
    suppressPackageStartupMessages(require("maps", quietly=TRUE))
    coast <- map("world", interior=FALSE, plot=FALSE)
    coast <- data.frame(lon=coast$x, lat=coast$y)
    # restrict the coastline info to what we need given the data
    expand <- 2       # add a little wiggle room
    # compute extent of data
    lats <- range(data$lat) + c(-expand, +expand)
    lons <- range(data$lon) + c(-expand, +expand)
    # re-cut the coastline
    # coast <- coast[coast$lat >= lats[1] & coast$lat <= lats[2] & coast$lon >= lons[1] & coast$lon <= lons[2],]
    coast <- coast[coast$lat <= lats[2] & coast$lon >= lons[1] & coast$lon <= lons[2],]

    # add the geom
    p <- p + geom_path(data=coast, na.rm=TRUE, colour="grey50")
    # NB: silently remove missing values which are inherent to coastline data
  }

  # use nice colours
  if ("fill" %in% names(mapping)) {
    fill.data <- data[,as.character(mapping$fill)]
    if (is.numeric(fill.data)) {
      # if the data is numeric, use a continuous, jet-like colour scale
      p <- p + scale_fill_gradientn(colours=continuous.colourmap())
    } else if (is.factor(fill.data)) {
      # if the data is a factor, use a discrete colour scape
      p <- p + scale_fill_manual(values=discrete.colourmap(n=nlevels(fill.data)))
    }
  }

  # no background
  p <- p + theme_bw()

  # nicer, simpler scales
  p <- p +
    # scale_x_continuous(name="", breaks=c(0)) +
    # NB: fails due to a bug in ggplot now, instead use
    opts(axis.text.x=theme_blank(), axis.title.x=theme_blank()) +
    scale_y_continuous(name="Latitude")

  return(p)
}


plot.env.data <- function(variables="", path=getOption("atlasr.env.data"), ...) {
  #
  # Plot all environmental data to PNG files
  #
  # variables   only output variables matching this (matches all when "")
  # path        path to environmental database
  #
  suppressPackageStartupMessages(require("stringr", quietly=TRUE))
  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
  suppressPackageStartupMessages(require("ncdf", quietly=TRUE))
  suppressPackageStartupMessages(require("reshape2", quietly=TRUE))

  # read all data files
  database <- read.env.data(variables=variables, path=path)

  # identify each element by its name and file of origin
  ncVariables <- names(database)
  ncFiles <- list.env.data(variables=variables, path=path, full=T)
  for (i in seq(along=database)) {
    database[[i]]$variable <- ncVariables[i]
    database[[i]]$file <- ncFiles[i]
  }

  # read coordinates of land masses
  land <- read.csv(str_c(path, "/worldmap-below_30-rough-no_countries.csv"))
  landLayer <- geom_polygon(aes(x=lon, y=lat), alpha=0.5, data=land)

  message("-> Plot variables")

  # loop on all files
  l_ply(database, function(x) {

    # convert into a data.frame, for ggplot
    d <- melt(x$z)
    names(d) <- c("x", "y", "z")
    d$x <- x$x[d$x]
    d$y <- x$y[d$y]

    # better variable name
    variable <- str_replace_all(x$variable, "_", "\n")

    # png file to plot into
    file <- str_replace(x$file, "\\.nc$", ".png")

    # plot in the png file
    png(file=file, width=1200, height=900, units="px", res=90)

    p <- ggplot(d) +
      # plot points
      geom_point(aes(x=x, y=y, colour=z), size=0.5) +
      # plot land
      landLayer +
      # nice colour gradient
      scale_colour_distiller(name=variable, palette="Spectral", guide="colourbar") +
      # blank theme
      opts(panel.background=theme_blank(),
           panel.grid.major=theme_blank(),
           axis.ticks=theme_blank(),
           axis.title.x=theme_blank(),
           axis.title.y=theme_blank(),
           axis.text.x=theme_blank(),
           axis.text.y=theme_blank(),
           plot.margins=c(0,0,0,0)
      ) +
      # polar view
      polar_proj()
    print(p)

    dev.off()

    # cleanup
    rm(d, p)

    return(NULL)

  }, .progress="text")

  return(invisible(NULL))
}

