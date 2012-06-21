#
#     Perform abiotic regionalisation
#
# (c) Copyright 2012 Ben Raymond, ben dot raymond at aad dot gov dot au
#                    Jean-Olivier Irisson
#     GNU General Public License v3
#     Last-Modified: <2012-Jun-21 17:04>
#
#-----------------------------------------------------------------------------

## Toolbox functions
#-----------------------------------------------------------------------------

flipud <- function(A) A[nrow(A):1,]

mcolor=function(x,y=NULL,z=NULL,interp=F,col=topo.colors(100),clim=NULL) {
    ## uses rasterImage() as a faster alternative to image()
    if (is.null(y) & is.null(z)) {
        z=x
        x=1:dim(z)[1]
        y=1:dim(z)[2]
    }
    ncolours=length(col)
    tempz=z
    if (!is.null(clim)) {
       tempz[tempz<clim[1]]=clim[1]
       tempz[tempz>clim[2]]=clim[2]
    }
    temp=round((flipud(t(tempz))-min(tempz,na.rm=T))*(ncolours-1)/(max(tempz,na.rm=T)-min(tempz,na.rm=T)))+1
    tempa=as.raster(col[temp],nrow=dim(temp)[1])

    ## show on figure
    xbin=mean(abs(diff(x)))
    ybin=mean(abs(diff(y)))
    plot(c(min(x,na.rm=T)-xbin/2,max(x,na.rm=T)+xbin/2),c(min(y,na.rm=T)-ybin/2,max(y,na.rm=T)+ybin/2),type="n",xlab="",ylab="",xlim=c(min(x,na.rm=T)-xbin/2-0.05,max(x,na.rm=T)+xbin/2+0.05),ylim=c(min(y,na.rm=T)-ybin/2-0.05,max(y,na.rm=T)+ybin/2+0.05),yaxs='i', xaxs='i')
    rasterImage(tempa,min(x,na.rm=T)-xbin/2,min(y,na.rm=T),max(x,na.rm=T)+xbin/2,max(y,na.rm=T),interp=F)
}


## Run bioregionalisation
#-----------------------------------------------------------------------------

compute.bioreg <- function(
  #
  # Perform bioregionalisation based on clustering
  #
  data,                       # data.frame containing only variables to cluster (after transforming, scaling, weighting, etc.)
  data.orig=data,             # original data.frame, before data processing (and usually with lon and lat); this is later used for plots etc.
  n.groups=12,                # either an integer number of groups (e.g. 10), or the height at which to cut the dendrogram (e.g. 0.13) in the hierarchical clustering step (hclust)
  n.groups.intermediate=200,  # number of groups in the non-hierarchical clustering step (clara); increasing this increases computation time significantly
  quick=TRUE,                 # TRUE produces less precise clustering, suitable for exploratory runs; FALSE is for high quality final analyses
  quiet=FALSE                 # output messages along the computation
)
{
  # load packages
  suppressPackageStartupMessages(require("cluster", quietly=TRUE))
  suppressPackageStartupMessages(require("vegan", quietly=TRUE))
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))

  # check arguments
  if (nrow(data) != nrow(data.orig)) {
    stop("Rows (i.e. locations) in data and data.orig do not match")
  }


  # remove lines (i.e. locations) NAs but keep track of them
  missing.mask <- rowSums(is.na(data)) > 0
  data.used <- data[!missing.mask,]


  if ( ! quiet ) cat("   Perform non-hierarchical clustering first\n")
  # For the later hierarchical clustering we will need to compte a distance matrix between all data points. This is obviously impossible on the full data set, so we reduce the information to a smaller number of similar clusters through non-hierarchical clustering

  # number of samples according to the quality argument (smaller numbers speed-up computation)
  if (quick) {
    samples <- 5
  } else {
    samples <- 50
  }

  # perform clustering
  cl <- clara(data.used, k=n.groups.intermediate, metric="manhattan", stand=FALSE, samples=samples)

  # associate non-hierachical cluster number with data
  data.used$clara.num <- cl$clustering


  if ( ! quiet ) cat("   Perform hierarchical clustering on the result\n")
  # Do a hierarchical clustering using the output of the nonhierarchical step. This defines the bioregions

  # first calculate mean properties of the non-hierarchical clusters
  data.mean <- ddply(data.used, ~clara.num, colMeans, na.rm=TRUE)

  # dissimilarities of these clusters
  D <- vegdist(data.mean[!names(data.mean) %in% "clara.num"], method="gower")

  # hierarchical clustering
  hcl <- hclust(D, method="ave")

  # now extract the desired number of groups from the dendrogram
  if (floor(n.groups)==n.groups) {
      # n.groups is integer, i.e. we specified a number of groups directly
      hcl$clustering <- cutree(hcl, k=n.groups)
      # work out the dissimilarity level (height) that corresponds to this number of groups
      hcl$cut.height <- mean(c(hcl$height[length(hcl$height)+2-n.groups], hcl$height[length(hcl$height)+2-n.groups-1]))
  } else {
      # we specified a height at which to cut the dendrogram
      # show on the dendrogram the height at which we are cutting
      hcl$cut.height <- n.groups
      hcl$clustering <- cutree(hcl, h=n.groups)
      n.groups <- length(unique(hcl$clustering))
  }
  # associate hierachical cluster number to each non-hierarchical cluster
  data.mean$hclust.num <- hcl$clustering


  # Store result
  # associate non-hierarchical cluster number to each data point
  data.orig$clara.num[!missing.mask] <- cl$clustering

  # associate hierarchical cluster number to each data point
  data.orig <- join(data.orig, data.mean[,c("clara.num", "hclust.num")], by="clara.num", type="full")
  data.orig <- rename(data.orig, c(hclust.num="cluster"))

  # transform into factors
  data.orig$clara.num <- factor(data.orig$clara.num)
  data.orig$cluster <- factor(data.orig$cluster)

  # store everything in a list
  res <- list(
    cl=cl,
    hcl=hcl,
    data=data.orig
  )
  class(res) <- c("bioreg", "list")

  return(res)
}

bioreg <- function(
  #
  # User-friendly interface to bioregionalisation
  #
  variables,                  # names (or abbreviations) of the variables to use for the bioregionalisation
  lat.min=-80, lat.max=-30,   # definition of the prediction grid
  lat.step=0.1,
  lon.min=-180, lon.max=180,
  lon.step=0.5,
  transformations=NULL,       # named vector of transformations applied to each variable (has to match variables)
  weights=NULL,               # named vector of weigths associated with each variable (has to match variables; will be scaled to a max of 1)
  quick=TRUE,                 # TRUE produces less precise clustering and faster plots;  FALSE is for high quality final analyses
  output.dir=NULL,            # destination for output files; if NULL, no output files will be saved and plots will be printed to the screen
  path=getOption("atlasr.env.data"),  # path to the environment database
  ...                         # passed to compute.bioreg()
)
{
    # load packages
    suppressPackageStartupMessages(require("plyr", quietly=TRUE))
    suppressPackageStartupMessages(require("reshape2", quietly=TRUE))
    suppressPackageStartupMessages(require("stringr", quietly=TRUE))
    suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))


    ## Get and transform data
    #-------------------------------------------------------------------------

    message("-> Get environmental variables")
    allVariables <- list.env.data(path=path)

    # expand variable names
    variables <- match.vars(variables, allVariables, quiet=FALSE)

    # check that there are enought variables
    if (length(variables) < 2) {
        stop("You must specify at least two input variables")
    }

    # get database
    database <- read.env.data(variables=variables, path=path)
    # remove information on land
    database <- mask.env.data(database, path=path)
    # build region of interest
    prediction_grid <- build.grid(
                          lat.min=lat.min, lat.max=lat.max, lat.step=lat.step,
                          lon.min=lon.min, lon.max=lon.max, lon.step=lon.step
    )
    # get environment data at the points of interest
    data.raw <- associate.env.data(prediction_grid, database)

    # transform data
    if (!is.null(transformations)) {
      data.transformed <- transform.data(data.raw[,! names(data.raw) %in% c("lon", "lat")], transformations)
    } else {
      data.transformed <- data.raw[,! names(data.raw) %in% c("lon", "lat")]
    }

    for (i in 1:ncol(data.transformed)) {
      # clean up any Inf values
      data.transformed[is.infinite(data.transformed[,i]),i] <- NA

      # normalise each column of x to 0-1 range
      data.transformed[,i] <- data.transformed[,i] - min(data.transformed[,i], na.rm=T)
      data.transformed[,i] <- data.transformed[,i] / max(data.transformed[,i], na.rm=T)
    }

    # weight data
    if (!is.null(weights)) {
      data.transformed <- weight.data(data.transformed, weights)
    }


    ## Run bioregionalisation
    #-------------------------------------------------------------------------
    message("-> Compute bioregions")
    bioregObj <- compute.bioreg(data=data.transformed, data.orig=data.raw, quick=quick, ...)


    ## Output data
    #-------------------------------------------------------------------------

    # should we save some data on disk ?
    output <- !is.null(output.dir)

    if (output) {
      output.dir <- clean.path(output.dir)
      message("-> Write output to ", output.dir)

      dir.create(output.dir, recursive=TRUE, showWarnings=FALSE)

      # use a date/time based suffix to differenciate between runs
      suffix <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      baseName <- str_c(output.dir,"/bioreg-", suffix)

      # Rdata with full object
      rdataFile <- str_c(baseName, ".Rdata")
      save(bioregObj, file=rdataFile)

      # Shapefiles of clusters
      write.shapefile(bioregObj$data, name=baseName, variables="cluster")

      # CSV file of clusters
      csvFile <- str_c(baseName, ".csv")
      write.table(bioregObj$data, file=csvFile, sep=",", row.names=FALSE)

      # TXT file with variables info
      infoFile <- str_c(baseName, "-info.txt")
      infoD <- data.frame(variables)
      if (!is.null(transformations)) {
        infoD$transformations <- attr(data.transformed, "transformations")
      } else {
        infoD$transformations <- ""
      }
      if (!is.null(weights)) {
        infoD$weights <- attr(data.transformed, "weights")
      } else {
        infoD$weights <- 1
      }
      write.table(infoD, file=infoFile, sep="\t", row.names=FALSE)
    }

    message("-> Produce plots")

    if (output) {
      # plot in a PDF file
      pdfFile <- str_c(baseName, ".pdf")
      pdf(pdfFile, width=11.7, height=8.3)
    }

    # dendrogram
    plot.dendro(bioregObj)

    # ask for further plots if we are using interactive output
    if (!output) devAskNewPage(TRUE)

    # plot of variables distributions within each cluster
    print(boxplotPlot <- plot.bioreg(bioregObj, geom="boxplot"))
    print(violinPlot <- plot.bioreg(bioregObj, geom="violin"))

    # map of clusters
    print(clusterPlot <- plot.pred.bioreg(bioregObj, quick=quick))

    if (!output) {
      # switch back to not asking for plots if we were working interactively and switched it on
      devAskNewPage(FALSE)
    } else {
      # close the PDF file otherwise
      dev.off()
    }

    message("Done")

    return(invisible(bioregObj))
    # NB: invisible() so that it doesn't get printed to console if not assigned to a variable
}


## Plots
#-----------------------------------------------------------------------------

plot.dendro <- function(x) {
  #
  # Plot the dendrogram of the final hierarchical clustering
  #
  # x     object of class "bioreg" resulting from a bioreg() or compute.bioreg() call
  #

  # re-extract variables from the object
  n.groups <- nlevels(x$data$cluster)
  n.groups.intermediate <- nlevels(x$data$clara.num)

  # plot the dendrogram
  plot(x$hcl, labels=F, hang=-1)

  # add the cutting level
  lines(c(1,n.groups.intermediate), c(x$hcl$cut.height, x$hcl$cut.height), lty=2, col="red")

  # add markers for group labels
  # get colour map
  cmap <- discrete.colourmap(n.groups)
  # reorder colours
  colours <- cmap[x$hcl$clustering][x$hcl$order]
  suppressPackageStartupMessages(require("scales", quietly=TRUE))   # for function alpha()
  points(1:n.groups.intermediate, rep(-0.004, n.groups.intermediate), col=alpha("black", 0.5), bg=colours, pch=22, cex=1)

  return(invisible(x))
}

plot.bioreg <- function(x, geom=c("violin", "boxplot"), ...) {
  #
  # Plot the effects in a bioregionalisation study:
  # the values of the variables in each cluster
  #
  # x     object of class "bioreg" resulting from a bioreg() or compute.bioreg() call
  # geom  type of plot: violin plot or boxplot
  #

  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))

  # extract the data.frame with predictions from the object
  x <- x$data

  geom <- match.arg(geom)

  # reformat data
  cols <- setdiff(names(x), c("lon", "lat", "clara.num"))
  xm <- melt(x[,cols], id.vars="cluster")
  xm <- na.omit(xm)

  # remove plotting box/violin plots for clusters with a very small number of points
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))
  counts <- plyr::count(x, "cluster")
  smallClusters <- na.omit(counts$cluster[which(counts$freq < 3)])

  xmB <- xm[!xm$cluster %in% smallClusters,]
  xmS <- xm[xm$cluster %in% smallClusters,]

  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
  if (geom == "violin") {
    p <- ggplot() + geom_violin(aes(x=cluster, y=value, fill=cluster), data=xmB)
  } else  {
    p <- ggplot() + geom_boxplot(aes(x=cluster, y=value, fill=cluster), data=xmB, outlier.size=1)
  }

  # plot points for small clusters, if any
  if (nrow(xmS) > 0) {
    p <- p + geom_point(aes(cluster, y=value, fill=cluster), data=xmS, shape=21, size=1.5)
  }

  # get colour map
  cmap <- discrete.colourmap(n=nlevels(x$cluster))
  p <- p +
        # modify scales to use nice colours and *not* drop empty clusters
        scale_fill_manual(values=cmap, drop=FALSE) +
        scale_colour_manual(values=cmap, drop=FALSE) +
        scale_x_discrete(drop=FALSE) +
        # split in facets
        facet_wrap(~variable, scales="free") +
        # orient horizontally
        coord_flip()

  return(p)
}

plot.pred.bioreg <- function(x, quick=FALSE, path=getOption("atlasr.env.data"), ...) {
  #
  # Plot a map of bioregionalisation clusters
  #
  # x     object of class "bioreg" resulting from a bioreg() or compute.bioreg() call
  # quick when TRUE, use a raster plot with no projection; when FALSE, use a vector plot with polar stereographic projection
  #

  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))

  # extract the data.frame with predictions from the object
  x <- x$data

  # get colours
  cmap <- discrete.colourmap(n=nlevels(x$cluster))

  if (quick) {
    # raster based plot

    clusterMap <- ggplot(x, aes(x=lon, y=lat)) +
      geom_raster(aes(fill=cluster)) +
      # no extra space
      scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
      # nice colours
      scale_fill_manual(values=c(cmap)) +
      # super-impose land
      layer_land(x) +
      # use theme_bw() otherwise holes (missing data) look very similar to the grey cluster
      theme_bw()

  } else {
    # polar projected plot

    clusterMap <- polar.ggplot(x, aes(fill=cluster)) +
      # nice colours
      scale_colour_manual(values=cmap)

  }

  return(clusterMap)
}


## GUI
#-----------------------------------------------------------------------------

do.bioreg <- function(lon.min=30, lon.max=60, lat.min=-62, lat.max=-45) {
  #
  # Open a GUI to select the arguments of the bioreg() function
  # First window: select variables and set all options
  #
  # Ben Raymond
  # Last-Modified: <2012-May-21 18:11>

  suppressPackageStartupMessages(require("rpanel"))

  # default window dimensions
  w <- 600          # width of the window
  windowH <- 700    # height of the window
  h <- 50           # height of elements
  spacer <- 10      # height of spacer

  # main window
  win <- rp.control(title="Run regionalisation analysis", size=c(w, windowH), aschar=F)

  # NB: positions are x, y, width, eight
  #     x, y are the coordinates of the top-left hand corner

  # destination directory for output files and then proceed to second GUI panel
  win$output.dir <- NULL
  rp.button(win, title="Choose output\ndirectory", pos=c(0, 0, w/4, h), action=function(win) {

    # choose the data file
    win$output.dir <- tk_choose.dir(getwd(), "Choose a suitable folder for the output files")
    # tkdestroy(win$window)

    # write the filename, for information
    rp.text(win, txt=win$output.dir, initval=win$output.dir, pos=c(w/4, 0, 3*w/4, h))

    return(win)
  })

  # variable list
  # available variables
  allVariables <- list.env.data()
  # dimensions
  listWcols <- 72   # compatible with 600px window
  listH <- 380
  listHrows <- 24   # fits in 380px height

  rp.listbox.mult(win, var=variables, vals=allVariables, title="Variables to use in the regionalisation", rows=listHrows, cols=listWcols, initval="", pos=c(0, h, w, listH), aschar=F, action=function(win) {
    return(win)
  })

  # current vertical dimension
  mid <- listH + h + spacer

  # options
  optionsFile <- tempfile()

  # variables transformation
  win <- rp.button(win, "Variables transformation", pos=c(0, mid, w/2, h) , action=function(win) {
    if (length(win$variables) < 2) {
      rp.messagebox("At least two variables must be selected", title="Warning")
    } else {
      do.bioreg.variables(variables=win$variables, file=optionsFile)
    }
    return(win)
  })


  # quality of run
  # better quality is slower, which can be painful at the exploratory stage
  rp.radiogroup(win, quick, values=c('Exploratory run (faster)','Final run (better quality)'), title="Analysis type", pos=c(0, mid+h, w/2, h*2))

  # number of clusters
  rp.slider(win, n.groups,  from=2, to=40,  resolution=1,   title="Number of clusters", initval=12 , showvalue=TRUE, pos=c(0, mid+3*h, w/2, h))

  # location
  rp.slider(win, lat.max,  from=-90, to=-30,  resolution=1,   title="North"   , initval=lat.max , showvalue=TRUE, pos=c(w/2+w/8, mid    , w/4, h))
  rp.slider(win, lon.min,  from=-180, to=180, resolution=1,   title="West"    , initval=lon.min,  showvalue=TRUE, pos=c(w/2, mid+h*1, w/4, h))
  rp.slider(win, lon.max,  from=-180, to=180, resolution=1,   title="East"    , initval=lon.max , showvalue=TRUE, pos=c(3*w/4, mid+h*1, w/4, h))
  rp.slider(win, lat.min,  from=-90, to=-30,  resolution=1,   title="South"   , initval=lat.min , showvalue=TRUE, pos=c(w/2+w/8, mid+h*2, w/4, h))
  rp.slider(win, lon.step, from=0.1, to=4,    resolution=0.1, title="Step lon", initval=0.1 ,     showvalue=TRUE, pos=c(w/2, mid+h*3, w/4, h))
  rp.slider(win, lat.step, from=0.1, to=4,    resolution=0.1, title="Step lat", initval=0.1 ,     showvalue=TRUE, pos=c(3*w/4, mid+h*3, w/4, h))

  # action buttons
  rowY <- mid+h*4+spacer
  rp.button(win, "Help", pos=c(0, rowY, w/4, h), action=function(win) {
    rp.messagebox("Someday... Maybe", title="Help")
    return(win)
  })
  rp.button(win, "Cancel", quitbutton=TRUE, pos=c(w/2, rowY, w/4, h) , action=function(win) {
    message("Aborting");
    return(win)
  })
  rp.button(win, "Run", pos=c(3*w/4, rowY, w/4, h), action=function(win) {

    # variables
    # print(win$variables)
    if (length(win$variables) < 2) {
      rp.messagebox("At least two variables must be selected", title="Warning")
    } else {

      # variables transformations
      if (file.exists(optionsFile)) {
        # read the options in the file
        load(optionsFile)
        transformations <- opts$transformations
        weights <- opts$weights
      } else {
        transformations <- NULL
        weights <- NULL
      }
      # print(transformations)
      # print(weights)


      # options
      quick <- TRUE
      if (grepl('^Final',win$quick)) {
        quick <- FALSE
      }
      # print(quick)

      # print other option values
      # print(win$output.dir)
      # print(win$lat.min)
      # print(win$lat.max)
      # print(win$lat.step)
      # print(win$lon.min)
      # print(win$lon.max)
      # print(win$lon.step)

      # build the function call
      message("Command:")
      call <- str_c("bioreg(",
        "variables=", str_c(deparse(win$variables, width=500), collapse=""),
        ", n.groups=", win$n.groups,
        ", weights=", str_c("c(",str_c(names(weights), unlist(weights), sep="=", collapse=", "), ")"),
        ", transformations=", str_c("c(",str_c(names(transformations), str_c("\"", unlist(transformations),"\"") , sep="=", collapse=", "), ")"),
        ", lat.min=", win$lat.min, ", lat.max=", win$lat.max, ", lat.step=", win$lat.step,
        ", lon.min=", win$lon.min, ", lon.max=", win$lon.max, ", lon.step=", win$lon.step,
        ", quick=", quick,
        ", output.dir=", deparse(win$output.dir),
        ")"
      )

      cat(call, "\n")

      # execute call
      b <- eval(parse(text=call))

    }

    return(win)
  })

}

do.bioreg.variables <- function(variables, file) {
  #
  # Second part of the GUI. Called from do.bioreg()
  # Select variable weights and transformations
  #
  # Ben Raymond
  # Last-Modified: <2012-05-02 12:34:54>

  suppressPackageStartupMessages(require("rpanel"))

  # defaults
  weights <- rep(1,length(variables))
  names(weights) <- variables
  transformations <- rep("",length(variables))
  names(transformations) <- variables

  # if the options file exists, read it to use the previously saved weights and transformations
  # TODO avoid using a file storage, we should be able to pass everything through functions and operate in memory
  if (file.exists(file)) {
    load(file)

    # transformations and weights are stored as lists in the options file
    # transform them back to vectors
    opts <- lapply(opts, function(x) {
      y <- as.character(x)
      names(y) <- names(x)
      return(y)
    })

    # the selected variables may have changed, use what we can from the options file (opts$***) and use defaults for the rest (***)
    commonNames <- intersect(names(transformations), names(opts$transformations))
    transformations[match(commonNames, names(transformations))] <- opts$transformations[match(commonNames, names(opts$transformations))]
    weights[match(commonNames, names(weights))] <- opts$weights[match(commonNames, names(opts$weights))]
  }

  # default dimensions (in px)
  w <- 1000         # width of the window
  windowH <- 300    # height of the window
  # TODO adapt it to the number of variables
  h <- 50           # height of elements
  spacer <- 10      # height of spacer

  # main window
  win <- rp.control(title="Variables", size=c(w, windowH), aschar=F)

  # transformations
  # TODO look into why with rp.textentry.immediate, the results are not all carried to the stage of the close button
  rp.textentry(win, var=transformsBox, labels=variables, title="Transformations", initval=transformations, pos=c(0,0,w/2,windowH/2), action=function(win) {
    return(win)
  })

  # provide example transformations
  example.transforms.labels=c("log10(x+1)","Square root","log10(-1*negative values only)")
  example.transforms.functions=list('"log10(x+1)"', '"sqrt(x)"', '"x[x>=0]=NA; log10(-x)"')
  rp.textentry(win, var=exampletransformBox, labels=example.transforms.labels, title="Example transformations", initval=example.transforms.functions, pos=c(0,windowH/2,w/2,windowH/2), action=function(win) {
    return(win)
  })

  # selection for weights associated with each variable
  rp.textentry(win, var=weightsBox, labels=variables, title="Weighting", initval=weights,pos=c(w/2, 0, w/2, windowH/2), action=function(win) {
    return(win)
  })

  # close button
  rp.button(win, title="Close", quitbutton=TRUE, pos=c(3/4*w, windowH-h, w/4, h), action=function(win) {
    # print("Transformations")
    # print(win$transformsBox)
    # print("Weights")
    # print(win$weightsBox)

    # extract transformations and weights and store them in named lists
    transformations <- as.list(win$transformsBox)
    weights <- as.list(win$weightsBox)
    weights <- lapply(weights, as.numeric)

    # store them in a single object
    opts <- list(transformations=transformations, weights=weights)

    # save this object to the options file
    save(opts, file=file)

    return(win)
  })

  return(invisible(win))
}
