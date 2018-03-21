#
#      Useful functions for plots of models from an Antarctic perspective
#
#  (c) Copyright 2011-2013 Jean-Olivier Irisson
#      http://creativecommons.org/licenses/by/3.0/
#
#-----------------------------------------------------------------------------


# Colour scales for ggplot
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


# Generics
plot.pred <- function(x, ...) {
  #
  # Generic for the plot of predictions from a model
  #
  UseMethod("plot.pred")
}

plot.effects <- function(x, ...) {
  #
  # Generic for the plot of predictions from a model
  #
  UseMethod("plot.effects")
}


# Polygon clipping
clip.polygon <- function(x, lon.min=-180, lon.max=180, lat.min=-90, lat.max=90) {
  #
  # Clip a polygon within given limits
  #
  # x       a data.frame with coordinates (lon and lat) defining a polygon, or several polygons separated by NAs
  # lon.*
  # lat.*   limits within which the polygon is cut
  #

  suppressPackageStartupMessages(require("rgeos", quietly=TRUE))
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))
  suppressPackageStartupMessages(require("sp", quietly=TRUE))
  suppressPackageStartupMessages(require("raster", quietly=TRUE))

  # convert data.frame x into SpatialPolygons object(s)

  # identify each polygon in the original data
  x$id <- cumsum(is.na(x$lon))
  x <- na.omit(x)
  xPolygon <- dlply(x, ~id, function(x) { Polygon(x[,c("lon","lat")]) } )

  # convert it into one big SpatialPolygons object
  # TODO we could actually keep several separate polygons as in the original data
  xPolygons <- Polygons(xPolygon, ID="foo")
  xSpatialPolygons <- SpatialPolygons(list(xPolygons), proj4string=CRS("+proj=longlat +datum=WGS84"))

  # do not cut the pole when the longitude spans its full range
  # if (lon.min <= -180 & lon.max >= 180) {
  #   if (lat.min < -50) {
  #     lat.min <- -91
  #   }
  #   if (lat.max > 50) {
  #     lat.max <- 91
  #   }
  # }
  # NB: this allows to get the full antarctic continent in polar view
  if ( lon.min < -180 ) { lon.min <- -180}
  if ( lon.max >  180 ) { lon.max <-  180}
  if ( lat.min < -90 ) { lat.min <- -90}
  if ( lat.max >  90 ) { lat.max <-  90}

  # create the clipping polygon
  clip.extent <- as(extent(lon.min, lon.max, lat.min, lat.max), "SpatialPolygons")
  proj4string(clip.extent) <- CRS(proj4string(xSpatialPolygons))

  # clip the map
  clipped <- gIntersection(xSpatialPolygons, clip.extent)

  # convert it back into a data.frame
  clipped <- ldply(clipped@polygons[[1]]@Polygons, function(x) {
    # extract coordinates
    X <- as.data.frame(x@coords)
    # separate from next polygon
    X <- rbind(X, NA)
    return(X)
  })
  names(clipped) <- c("lon", "lat")  
  # ggplot(clipped) + geom_polygon(aes(x=lon, y=lat))

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


# South pole plots in ggplot
south_pole_proj <- function(projection="stereographic", orientation=c(-90,0,0)) {
   #
   # Easy access to a suitable polar projection
   # NB: view from south pole (-90)
   #
   suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
   suppressPackageStartupMessages(require("mapproj", quietly=TRUE))

   out <- list(
      # polar projection
      coord_map(projection=projection, orientation=orientation),
      # simpler scales (because polar projection screws up longitude scales)
      theme(axis.text.x=element_blank(),axis.title.x=element_blank())#,
      # scale_y_continuous(name="Latitude")
   )

   return(out)
}

layer_land <- function(x, expand=1, path=getOption("atlasr.env.data"), ...) {
  #
  # Produce a ggplot layer with coastlines corresponding to a dataset
  #
  # x       dataset of interest, with columns lon and lat at least
  #         the boundaries of the land layer are taken from this
  # expand  how much to expand the land layer around the data (in lon/lat units)
  # path    path to environmental database
  #

  # compute data range
  lonRange <- diff(range(x$lon, na.rm=T))
  latRange <- diff(range(x$lat, na.rm=T))
  minRange <- min(lonRange, latRange) + expand * 2

  # choose the resolution based on the range of the data
  if (minRange < 35) {
    resolution <- 0.01
  } else if (minRange < 45){
    resolution <- 0.05
  } else {
    resolution <- 0.1
  }

  # read coordinates of land masses at the appropriate resolution
  land <- read.csv(str_c(path, "/coastline/worldmap-", resolution, ".csv"))

  # clip land masses to data range
  land <- clip.to.data(land, x, expand=expand)

  # draw the layer
  landLayer <- geom_polygon(aes(x=lon, y=lat), data=land, na.rm=T, ...)

  return(landLayer)
}

polar_ggplot <- function(data, mapping=aes(), geom=c("auto", "raster", "point", "tile"), path=getOption("atlasr.env.data"), land=NULL, scale=1, ...) {
  #
  # data          data frame with columns lat, lon, and variables to plot
  # mapping       a call to `aes()` which maps a variable to a plotting
  #               aesthetic characteristic (fill, size, alpha, etc.)
  # geom          the type of plot ("geometry" in ggplot parlance) to produce:
  #               raster  fast but does not allow projection
  #               points  reasonably fast, allows projection
  #               tiles   slow but better looking, allows projection
  #               auto chooses points or tiles depending on the size of the data
  # lat.precision
  # lon.precision the precision at which lat and lon are considered
  #               (in degrees). If they are larger than the original
  #               precision, the data is *subsampled* to those locations
  #               i.e. some data is actually dropped. If you want to average or
  #               sum the data per cell, use rasterize() in lib_data.R
  # land          layer (result of a geom call) to show the land on top of the data; NULL uses layer_land to produce it
  # scale         scale of the points plotted
  # ...           passed to the selected geom
  #

  suppressPackageStartupMessages(require("plyr", quietly=TRUE))
  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
  suppressPackageStartupMessages(require("stringr", quietly=TRUE))

  # Check arguments
  # geoms
  geom <- match.arg(geom)

  # check that we have something that looks like lat and lon
  if (! all(c("lat","lon") %in% names(data)) ) {
    stop("Need two columns named lat and lon to be able to plot\nYou have ", paste(names(data), collapse=", "))
  }

  # define the geom depending on the size of the data when geom="auto"
  n <- nrow(data)
  if (geom == "auto") {
    if (n <= 2000) {
      geom <- "tile"
    } else if (n <= 50000) {
      geom <- "point"
    } else {
      geom <- "raster"
    }
  }

  # remove colour mapping (one must use fill because we only use filled shapes)
  mapping <- mapping[names(mapping) != "colour",]

  # Plot
  # prepare plot
  p <- ggplot(data, aes(x=lon, y=lat))

  # plot points or tiles depending on the geom argument
  if (geom == "point") {
    # add mapping of size to better cover the space (smaller points near the center)
    mapping = c(aes(size=lat), mapping)
    class(mapping) = "uneval"

    # try to guess the size of points based on how many different latitudes there are
    # (the more data, the smaller the points)
    nLats <- length(unique(data$lat))
    baseSize <- scale * 90 / nLats

    # plot
    # NB: shape: 21 = filled point, 22 = filled square, 23 = filled losange
    p <- p + geom_point(mapping=mapping, shape=21, stroke=0, ...) + scale_size(range=c(baseSize, baseSize*2.2), guide=FALSE)
    } else if (geom == "tile") {
      p <- p + geom_tile(mapping=mapping, ...)

    } else if (geom == "raster") {
      p <- p + geom_raster(mapping=mapping, ...)
    }

  # use nice fill colours
  if ( "fill" %in% names(mapping) ) {
     # get the datat that is mapped to fill
    fill.data <- data[,as.character(mapping$fill)]

    if (is.numeric(fill.data)) {
      # if the data is numeric, use a continuous, jet-like colour scale
      p <- p + scale_fill_gradientn(colours=continuous.colourmap(), guide="colourbar")

    } else if (is.factor(fill.data)) {
      # if the data is a factor, use a discrete colour scape
      p <- p + scale_fill_manual(values=discrete.colourmap(n=nlevels(fill.data)))
    }
  }

  # no background
  p <- p + theme_bw()

  # remove boundaries
  p <- p + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))

  return(p)
}
