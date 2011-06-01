#
#     Plot data viewed from the South Pole
#
# The plot uses ggplot2 rather than base graphics, because it eases
# the projection and mapping of colours.
#
# data          data frame with columns lat, long, and variables to plot
# mapping       a call to `aes()` which maps a variable to a plotting
#               aesthetic characteristic (fill, colour, size, alpha, etc.)
# geom          the type of plot ("geometry" in ggplot parlance) to produce
#               = points (the default) or tiles
# lat.precision
# lon.precision the precision at which lat and lon are considered
#               (in degrees). If they are larger than the original
#               precision, the data is averaged within the new cells
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

polar.ggplot <- function(data, mapping=aes(), geom=c("points", "tiles"), lat.precision=NULL, lon.precision=NULL, ...) {

  suppressMessages(require(plyr, quietly=TRUE))
  suppressMessages(require(ggplot2, quietly=TRUE))
  suppressMessages(require(maps, quietly=TRUE))

  # Check arguments
  geom = match.arg(geom)

  # Allow long/lat to be called more liberally
  names(data)[tolower(names(data)) %in% c("latitude","lat")] = "lat"
  names(data)[tolower(names(data)) %in% c("longitude","lon","long")] = "long"
  # check that we have someting that looks like lat and long
  if (! all(c("lat","long") %in% names(data)) ) {
    stop("Need two columns named lat and long to be able to plot\nYou have ", paste(names(data), collapse=", "))
  }

  # If new precisions are specified for lon or lat, regrid the data (for speed purposes)
  if (! all(is.null(c(lat.precision, lon.precision))) ) {

    # round latitude and longitude to the given precision
    if (!is.null(lat.precision)) {
      data$lat = round_any(data$lat, lat.precision)
    }
    if (!is.null(lon.precision)) {
      data$long = round_any(data$long, lon.precision)
    }

    # compute the average within the new grid cells
    data = ddply(data, ~long+lat, mean)
  }

  # Plot
  # prepare plot
  p = ggplot(data, aes(x=long, y=lat)) +
        # stereographic projection
      	# NB: view from south pole (-90)
      	coord_map(projection="stereographic", orientation=c(-90,0,0))

	# plot points or tiles depending on the geom argument
  p = p + switch(geom,
    points = geom_point(mapping=mapping),
    tiles = geom_tile(mapping=mapping),
  )

	# plot the coastline
  coast = map("world", interior=FALSE, plot=FALSE)
  coast = data.frame(long=coast$x, lat=coast$y)
  coast = coast[coast$lat <= max(data$lat)+2,]
	p = p + geom_path(data=coast)

  # use nicer colours
  if ("fill" %in% names(mapping)) {
    fill.data = data[,as.character(mapping$fill)]
    if (is.numeric(fill.data)) {
      # if the data is numeric, use a yellow to red gradient
      p = p + scale_fill_gradient(low="#FAF3A9", high="#F62B32")
    } else if (is.factor(fill.data)) {
      # if the data is discrete, use a colorbrewer scale if possible (less than 12 colours)
      if (nlevels(fill.data)<=12) {
        p = p + scale_fill_brewer(palette="Set3")
      }
      # otherwise just use the default colours of ggplot
    }
  }
  # same for coulour
  if ("colour" %in% names(mapping)) {
    colour.data = data[,as.character(mapping$colour)]
    if (is.numeric(colour.data)) {
      p = p + scale_colour_gradient(low="#FAF3A9", high="#F62B32")
    } else if (is.factor(colour.data)) {
      if (nlevels(colour.data)<=12) {
        p = p + scale_colour_brewer(palette="Set3")
      }
    }
  }

	# nicer, simpler scales
	p = p + scale_x_continuous(breaks=180) + labs(x="Longitude", y="Latitude")

	# no background
	p = p +	theme_bw()

  return(p)
}

