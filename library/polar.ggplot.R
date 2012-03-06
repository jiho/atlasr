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


polar.ggplot <- function(data, mapping=aes(), geom=c("points", "tiles"), lat.precision=NULL, lon.precision=NULL, coast=NULL, ...) {
  #
  # data          data frame with columns lat, lon, and variables to plot
  # mapping       a call to `aes()` which maps a variable to a plotting
  #               aesthetic characteristic (fill, colour, size, alpha, etc.)
  # geom          the type of plot ("geometry" in ggplot parlance) to produce
  #               = points (the default) or tiles
  # lat.precision
  # lon.precision the precision at which lat and lon are considered
  #               (in degrees). If they are larger than the original
  #               precision, the data is averaged within the new cells
  # coast         coastline geom, if none is provided, a basic coastline is drawn
  #

  suppressPackageStartupMessages(require("plyr", quietly=TRUE))
  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))

  # Check arguments
  geom <- match.arg(geom)

  # Allow lon/lat to be called more liberally
  names(data)[tolower(names(data)) %in% c("latitude","lat")] <- "lat"
  names(data)[tolower(names(data)) %in% c("longitude","lon","long")] <- "lon"
  # check that we have something that looks like lat and lon
  if (! all(c("lat","lon") %in% names(data)) ) {
    stop("Need two columns named lat and lon to be able to plot\nYou have ", paste(names(data), collapse=", "))
  }

  # If new precisions are specified for lon or lat, regrid the data (for speed purposes)
  if (! all(is.null(c(lat.precision, lon.precision))) ) {

    # round latitude and longitude to the given precision
    if (!is.null(lat.precision)) {
      data$lat <- round_any(data$lat, lat.precision)
    }
    if (!is.null(lon.precision)) {
      data$lon <- round_any(data$lon, lon.precision)
    }

    # compute the average within the new grid cells
    # NB: use only numeric columns to avoid issues with mean
    data <- ddply(data[, sapply(data, is.numeric)], ~lon+lat, colMeans, na.rm=T)
  }

  # Get and re-cut coastline if none is provided
  if (is.null(coast)) {
    # extract the whole world
    suppressPackageStartupMessages(require("maps", quietly=TRUE))
    coast <- map("world", interior=FALSE, plot=FALSE)
    coast <- data.frame(lon=coast$x, lat=coast$y)    
    # restrict the coastline info to what we need given the data
    expand <- 1       # add a little wiggle room
    # compute extent of data
    lats <- range(data$lat) + c(-expand, +expand)
    lons <- range(data$lon) + c(-expand, +expand)
    # re-cut the coastline
    # coast <- coast[coast$lat >= lats[1] & coast$lat <= lats[2] & coast$lon >= lons[1] & coast$lon <= lons[2],]
    coast <- coast[coast$lat <= lats[2] & coast$lon >= lons[1] & coast$lon <= lons[2],]
    
    # prepare the geom
    coast <- geom_path(data=coast, na.rm=TRUE, fill="grey50")
    # NB: silently remove missing values which are inherent to coastline data
  }

  # Plot
  # prepare plot
  p <- ggplot(data, aes(x=lon, y=lat)) +
        # stereographic projection
        # NB: view from south pole (-90)
        coord_map(projection="stereographic", orientation=c(-90,0,0))

  # plot points or tiles depending on the geom argument
  p <- p + switch(geom,
    points = geom_point(mapping=mapping),
    tiles  = geom_tile(mapping=mapping)
  )

  # plot the coastline
  p <- p + coast

  # use nicer colours
  if ("fill" %in% names(mapping)) {
    fill.data <- data[,as.character(mapping$fill)]
    if (is.numeric(fill.data)) {
      # if the data is numeric, use a yellow to red gradient
      p <- p + scale_fill_gradient(low="#FAF3A9", high="#F62B32")
    } else if (is.factor(fill.data)) {
      # if the data is discrete, use a colorbrewer scale if possible (less than 12 colours)
      if (nlevels(fill.data)<=12) {
        p <- p + scale_fill_brewer(palette="Set3")
      }
      # otherwise just use the default colours of ggplot
    }
  }
  # same for coulour
  if ("colour" %in% names(mapping)) {
    colour.data <- data[,as.character(mapping$colour)]
    if (is.numeric(colour.data)) {
      p <- p + scale_colour_continuous(low="#FAF3A9", high="#F62B32")
    } else if (is.factor(colour.data)) {
      if (nlevels(colour.data)<=12) {
        p <- p + scale_colour_brewer(palette="Set3")
      }
    }
  }

  # nicer, simpler scales
  p <- p + labs(x="Longitude", y="Latitude")

  # no background
  p <- p + theme_bw()

  return(p)
}

