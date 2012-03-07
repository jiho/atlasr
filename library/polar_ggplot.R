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
  c <- coord_map(projection=projection, orientation=orientation)
  return(c)
}


spectral_colours <- function(aesthetic=c("fill", "colour"), ...) {
  #
  # Spectral colour scale from http://colorbrewer2.org/ used as a continuous scale in ggplot2
  #
  # aesthetic   type of scale to generate (fill or colour)
  # ...         passed to scale_***_gradientn
  #
  aesthetic <- match.arg(aesthetic)

  # get the palette
  suppressPackageStartupMessages(require("RColorBrewer", quietly=TRUE))
  colours <- rev(brewer.pal(6, "Spectral"))

  if (aesthetic == "fill") {
    s <- scale_fill_gradientn(colours=colours, guide=guide_colorbar(), ...)
  } else {
    s <- scale_colour_gradientn(colours=colours, guide=guide_colorbar(), ...)
  }

  return(s)
}
spectral_colors <- spectral_colours


## Plotting functions
#-----------------------------------------------------------------------------

polar.ggplot <- function(data, mapping=aes(), geom=c("point", "tile"), lat.precision=NULL, lon.precision=NULL, coast=NULL, ...) {
  #
  # data          data frame with columns lat, lon, and variables to plot
  # mapping       a call to `aes()` which maps a variable to a plotting
  #               aesthetic characteristic (fill, colour, size, alpha, etc.)
  # geom          the type of plot ("geometry" in ggplot parlance) to produce
  #               = points (the default) or tiles (possibly better looking, longer)
  # lat.precision
  # lon.precision the precision at which lat and lon are considered
  #               (in degrees). If they are larger than the original
  #               precision, the data is averaged within the new cells
  # coast         coastline geom, if none is provided, a basic coastline is drawn
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

  # if we need to re-grid, we need both a lat and a lon precision
  if (! all(is.null(c(lat.precision, lon.precision))) ) {
    stop("Need to provide both lat and lon precision to be able to regrid the data")
  }


  # If new precisions are specified for lon or lat, regrid the data (for speed purposes)
  if (! all(is.null(c(lat.precision, lon.precision))) ) {

    # NB: use only numeric columns to avoid issues when rasterizing
    numColumns <- sapply(data, is.numeric)
    if (sum(numColumns) != ncols(data)) {
      warning("Columns", str_c(names(data)[!numColumns], collapse=", "), "were deleted to allow regriding")
    }

    # compute mean per bin
    data <- rasterize(data[numColumns], vars=c("lat", "lon"), precisions=c(lat.precision, lon.precision), fun=mean, na.rm=T)
  }


  # Get and re-cut coastline if none is provided
  if (is.null(coast)) {
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

    # prepare the geom
    coast <- geom_path(data=coast, na.rm=TRUE, colour="grey50")
    # NB: silently remove missing values which are inherent to coastline data
  }


  # Plot
  # prepare plot
  p <- ggplot(data, aes(x=lon, y=lat)) +
        # stereographic projection
        polar_proj()

  # plot points or tiles depending on the geom argument
  p <- p + switch(geom,
    point = geom_point(mapping=mapping),
    tile  = geom_tile(mapping=mapping)
  )

  # plot the coastline
  p <- p + coast

  # use nicer colours
  if ("fill" %in% names(mapping)) {
    fill.data <- data[,as.character(mapping$fill)]
    if (is.numeric(fill.data)) {
      # if the data is numeric, use a yellow to red gradient
      p <- p + scale_fill_continuous(low="#FAF3A9", high="#F62B32")
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


plot.env.data <- function(variables="", path="env_data", ...) {
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
      # nice colour gradient
      spectral_colours(aes="colour", name=variable) +
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