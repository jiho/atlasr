#
#     Perform abiotic regionalisation
#
# (c) Copyright 2012 Ben Raymond, ben dot raymond at aad dot gov dot au
#     Last-Modified: <2012-05-14 14:16:34>
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

function.maker <- function(str) {
    #
    # Transform a character string into a function
    #
    # str   character string defining the function
    #

    suppressPackageStartupMessages(require("stringr", quietly=TRUE))

    # empty function
    f <- function(x) {}
    environment(f) <- baseenv()

    str <- str_c("{",str,"}") # make sure str is enclosed in curly brackets (will it matter if user also supplies these? - to check)

    # fill in the body of the function
    body(f) <- substitute(tryCatch(expr,
                                   error=function(e) "Error applying transformation function"),
                          list(expr=parse(text=str)[[1]])
    )

    return(f)
}


## Run bioregionalisation
#-----------------------------------------------------------------------------

bioreg <- function(variables, n.groups=12, lat.min=-80, lat.max=-30, lat.step=0.1, lon.min=-180, lon.max=180, lon.step=0.5, transformations=NULL, weights=NULL, quality=c("low","high"), path="env_data", output.dir=NULL)
{
    #
    # Perform bioregionalisation based on clustering
    #
    # variables     vector of names of environmental variables used in the bioregionalisation
    # n.groups      either an integer number of groups, or the height at which to cut the dendrogram (e.g. 0.13)
    # l**.min
    # l**.max
    # l**.step      definition of the grid on which the clustering will be done
    # transformations   list giving transformation function for each variable, or NULL for no transformations
    # weights       vector giving the weight for each variable
    # quality       "low" or "high"; low quality is faster, suitable for exploratory runs; high quality for final analyses
    # path          path where the environmental data is to be found
    # output.dir    destination for output files; if NULL, no output files will be saved

    # load packages
    suppressPackageStartupMessages(require("cluster", quietly=TRUE))
    suppressPackageStartupMessages(require("vegan", quietly=TRUE))
    suppressPackageStartupMessages(require("plyr", quietly=TRUE))
    suppressPackageStartupMessages(require("reshape2", quietly=TRUE))
    suppressPackageStartupMessages(require("stringr", quietly=TRUE))
    suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))


    # Check input arguments
    # variables
    if (length(variables) < 2) {
        stop("You must specify at least two input variables")
    }

    # quality
    quality <- match.arg(quality)

    # weights
    if (is.null(weights)) {
        weights <- rep(1,length(variables))
    }
    weights <- as.numeric(weights)
    weights <- weights / max(weights) # normalize so that max weight is 1

    # transformation functions
    if (!is.null(transformations)) {
        tfuncs <- list()
        for (i in seq(along=transformations)) {
            if (is.character(transformations[[i]]) & nchar(transformations[[i]])>0) {
                # convert from string expression into actual function
                tryCatch(
                    tfuncs[[i]] <- function.maker(transformations[[i]]),
                    error=function(e) stop("Error parsing transformation function : ", transformations[[i]], "\n  ", e)
                )

            } else if (is.function(transformations[[i]])) {
                # store actual functions
                tfuncs[[i]] <- transformations[[i]]

            } else if (is.null(transformations[[i]]) | is.na(transformations[[i]]) | (is.character(transformations[[i]]) & nchar(transformations[[i]])==0)) {
                # pass-through for the rest
                tfuncs[[i]] <- function(x){x}

            } else {
                stop("Supplied transformation is neither a string nor a function")
            }
        }
        transformations <- tfuncs
    }

    # Get and transform data
    # get database
    database <- read.env.data(variables=variables, path=path, match.names=F)
    # remove information on land
    database <- mask.env.data(database, path=path)
    # build region of interest
    prediction_grid <- build.grid(
                          lat.min=lat.min, lat.max=lat.max, lat.step=lat.step,
                          lon.min=lon.min, lon.max=lon.max, lon.step=lon.step
    )
    # get environment data at the points of interest
    data.raw <- associate.env.data(prediction_grid, database)

    data.transformed <- data.raw[,! names(data.raw) %in% c("lon", "lat")]
    for (i in 1:ncol(data.transformed)) {
        # apply user supplied transformations
        # TODO this is fragile : dependent on ordering of weights and transformations and data columns being the same: need to code it better
        # TODO: suggestion JO: use named lists/vectors for transformation and weights with names matching data columns? i.e. list(bathymetry="log(x)", floor_temperature=ceiling). But it's much more cumbersome to write then
        if (!is.null(transformations[[i]])) {
            data.transformed[,i] <- transformations[[i]](data.raw[,c(names(data.transformed)[i])]) # can't use i to index data.raw, because it has lon and lat cols and data.transformed does not!
        }

        # clean up any Inf values
        data.transformed[is.infinite(data.transformed[,i]),i] <- NA

        # normalise each column of x to 0-1 range
        data.transformed[,i] <- data.transformed[,i] - min(data.transformed[,i], na.rm=T)
        data.transformed[,i] <- data.transformed[,i] / max(data.transformed[,i], na.rm=T)

        # apply weighting
        # cat(sprintf("i=%d, weights[i-2]: %s\n",i,str(weights[[i-2]])))
        data.transformed[,i] <- data.transformed[,i] * weights[i]
    }

    # record which lines (i.e. locations) are masked out because of missing data (including land)
    missing.mask <- rowSums(is.na(data.transformed)) > 0
    data.trans.noNA <- na.omit(data.transformed)

    message("-> Perform non-hierarchical clustering first")
    # For the later hierarchical clustering we will need to compte a distance matrix between all data points. This is obviously impossible on the full data set, so we reduce the information to a smaller number of similar clusters through non-hierarchical clustering
    # number of clusters
    num.groups.intermediate <- 200

    # number of samples according to the quality argument (smaller numbers speed-up computation)
    samples <- switch(quality, low=5, high=50)

    # perform clustering
    cl <- clara(data.trans.noNA, k=num.groups.intermediate, metric="manhattan", stand=FALSE, samples=samples)

    # extract cluster numbers
    data.trans.noNA$clara.num <- cl$clustering


    message("-> Perform hierarchical clustering on the result")
    # Do a hierarchical clustering using the output of the nonhierarchical step. This defines the bioregions

    # first calculate mean properties of the non-hierarchical clusters
    xc <- ddply(data.trans.noNA, ~clara.num, colMeans, na.rm=TRUE)

    # dissimilarities of these clusters
    D <- vegdist(xc[!names(xc) %in% "cluster.num"], method="gower")

    # hierarchical clustering
    hcl <- hclust(D, method="ave")

    # now extract the desired number of groups from the dendrogram
    if (floor(n.groups)==n.groups) {
        # n.groups is integer, i.e. we specified a number of groups directly
        hclust.num <- cutree(hcl, k=n.groups)
        # work out the dissimilarity level (height) that corresponds to this number of groups
        temph <- mean(c(hcl$height[length(hcl$height)+2-n.groups], hcl$height[length(hcl$height)+2-n.groups-1]))
    } else {
        # we specified a height at which to cut the dendrogram
        # show on the dendrogram the height at which we are cutting
        temph <- n.groups
        hclust.num <- cutree(hcl, h=n.groups)
        n.groups <- length(unique(cn.new))
    }
    # associate hierachical cluster number to each non-hierarchical cluster
    xc$hclust.num <- hclust.num

    # associate hierarchical cluster number to each data point on the total grid
    # non-hierarchical cluster number
    data.raw$clara.num[!missing.mask] <- cl$clustering
    # hierarchical cluster number
    data.raw <- join(data.raw, xc[,c("clara.num", "hclust.num")], by="clara.num", type="full")
    data.raw <- rename(data.raw, c(hclust.num="cluster"))
    data.raw$cluster <- factor(data.raw$cluster)

    # define the output to be of class "bioreg"
    class(data.raw) <- c("bioreg", class(data.raw))

    message("-> Produce plots")

    # get a nice colour map
    cmap <- discrete.colourmap(n.groups)

    # Dendrogram
    dev.new()
    # dendrogram
    plot(hcl, labels=F, hang=-1)
    # cutting level
    lines(c(1,num.groups.intermediate), c(temph,temph), lty=2, col="red")
    # markers for group labels
    colours <- cmap[hclust.num][hcl$order]
    points(1:200, rep(-0.02, num.groups.intermediate), col=NA, bg=colours, pch=21, cex=1)

    # plot of variables distributions within each cluster
    dev.new()
    variablesPlot <- plot.bioreg(data.raw)
    print(variablesPlot)

    # Image map
    dev.new()
    if (quality=="low") {
        # read coordinates of land masses
        land <- read.csv(str_c(path, "/worldmap-below_30-rough-no_countries.csv"))
        if (F) {
            ## ggplot-based version: this has problems - the land layer does not line up with the raster, and the colours in the raster don't match the clustering results quite right
            landLayer <- geom_polygon(aes(x=lon, y=lat), alpha=0.5, data=land)
            clusterMap <- ggplot(data.raw, aes(x=lon, y=lat)) + geom_raster(aes(fill=cluster)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + scale_fill_manual(values=c(cmap)) + landLayer + theme_bw() # use theme_bw() otherwise missing data (grey) looks very similar to the grey cluster
            # TODO use polar.ggplot with subsampling and geom_tile instead? # BR- would prefer not, it will lose too much detail
            print(clusterMap)
        } else {
            temp <- join(prediction_grid, data.raw[,c("lon","lat","cluster")], by=c("lon","lat"))
            nrow=attr(prediction_grid,'out.attrs')$dim[1]
            temp$cluster=as.numeric(temp$cluster)
            mcolor(matrix(temp$lon,nrow=nrow),matrix(temp$lat,nrow=nrow),matrix(temp$cluster,nrow=nrow),col=cmap)
            lines(land$lon,land$lat)
        }
    } else {
        clusterMap <- polar.ggplot(data.raw, geom="point", aes(colour=cluster)) + scale_colour_manual(values=cmap)
        ## TODO fix error when longitudes of [-180,180] are used: gap in plot at lon==180
        print(clusterMap)
    }

    # Output data
    if (!is.null(output.dir)) {
        output.file <- normalizePath(str_c(output.dir,"/bioreg.Rdata"), winslash="/", mustWork=F)
        save(data.raw, data.transformed, file=output.file)
    }
    invisible(data.raw) # invisible() so that it doesn't get printed to console if not assigned to a variable
}

## Plots
#-----------------------------------------------------------------------------

plot.bioreg <- function(x, geom=c("violin", "boxplot"), ...) {
  #
  # Plot the effects in a bioregionalisation study:
  # the values of the variables in each cluster
  #
  # x     data.frame of class "bioreg" resulting from a bioreg() call
  # geom  type of plot: violin plot or boxplot
  #

  geom <- match.arg(geom)

  # reformat data
  cols <- setdiff(names(x), c("lon", "lat", "clara.num"))
  xm <- melt(x[,cols], id.vars="cluster")
  xm <- na.omit(xm)

  # remove plotting box/violin plots for clusters with a very small number of points
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))
  counts <- count(x, "cluster")
  smallClusters <- na.omit(counts$cluster[which(counts$freq < 3)])

  xmB <- xm[!xm$cluster %in% smallClusters,]
  xmS <- xm[xm$cluster %in% smallClusters,]

  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
  if (geom == "violin") {
    p <- ggplot() + geom_violin(aes(x=cluster, y=value, fill=cluster), data=xmB)
  } else  {
    p <- ggplot() + geom_boxplot(aes(x=cluster, y=value, fill=cluster), data=xmB)
  }

  # plot points for small clusters
  p <- p + geom_point(aes(cluster, y=value, fill=cluster), data=xmS, shape=21, size=1.5)

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

## GUI
#-----------------------------------------------------------------------------

do.bioreg <- function(...) {
    ##
    ## Open a GUI to select the arguments of the bioreg() function
    ##
    ## ...   passed to bioreg()
    ##
    ## Ben Raymond
    ## Last-Modified: <2012-05-10 10:22:33>

  suppressPackageStartupMessages(require("rpanel"))


  # default dimensions (in px)
  w <- 600      # width of the window
  h <- 50       # height of elements
  spacer <- 10  # height of spacer
  main.height=700

  # main window
  win <- rp.control(title="Regionalisation", size=c(w,main.height),aschar=F)


  nRows <- 30
  nBoxes <- 20
  hSel <- 380
  allVariables <- list.env.data()

  rp.text(win,txt="Choose the variables to use in the regionalisation analysis",pos=c(spacer,spacer,w-2*spacer,40))

  blah=rp.listbox.mult(win, var=availableVariables, vals=allVariables, title="Available variables",  rows=nRows, cols=min(50,max(sapply(allVariables,nchar))), initval="", pos=c(spacer, 40+spacer, w-2*spacer, main.height-40-h-2*spacer), aschar=F, action=function(win) {
      if (all(win$availableVariables == "")) {
          rp.messagebox("At least two variables must be selected", title="Warning")
      }
      return(win)
  })

  ## destination directory for output files and then proceed to second GUI panel
  rp.button(win, title="Choose directory for output files", pos=c(spacer,main.height-h-spacer,w-2*spacer,h), action=function(win) {
      if (length(win$availableVariables)<2) {
          rp.messagebox("At least two variables must be selected first", title="Warning")
      } else {
          outputDir=tk_choose.dir(getwd(), "Choose a suitable folder for the output files")
          tkdestroy(win$window)
          bioreg.secondpanel(win$availableVariables,output.dir=outputDir,...)
      }
      return(win) })
}

bioreg.secondpanel <- function(selectedVariables,varweights=rep(1,length(selectedVariables)),vartransforms=rep("",length(selectedVariables)),output.dir=getwd(),lon.min=30,lon.max=60,lat.min=-62,lat.max=-45) {
    ##
    ## Second part of the GUI. Called from do.bioreg()
    ##
    ## Ben Raymond
    ## Last-Modified: <2012-05-02 12:34:54>

  suppressPackageStartupMessages(require("rpanel"))


  # default dimensions (in px)
  w <- 1200      # width of the window
  h <- 50       # height of elements
  spacer <- 10  # height of spacer
  main.height=700

  # main window
  win <- rp.control(title="Regionalisation", size=c(w,main.height),aschar=F)


  nRows <- 20
  nBoxes <- 20
  hSel <- 380


  ## selection for weights associated with each variable
  rp.textentry.immediate(win, var=weightbox, labels=selectedVariables, title="Weighting", initval=varweights,pos=c(spacer,spacer,w/3-spacer,main.height/2-spacer),
               action=function(win) {
                   return(win) })

  rp.textentry.immediate(win, var=transformbox, labels=selectedVariables, title="Transformations", initval=vartransforms,pos=c(w/3+spacer,spacer,w/3-spacer,main.height/2-spacer),
               action=function(win) {
                   return(win) })

  ## provide example transformations
  example.transforms.labels=c('log10(-1*negative values only)','log10(x+1)','Square root')
  example.transforms.functions=list('"x[x>=0]=NA; log10(-x)"','"log10(x+1)"','"sqrt(x)"')
  rp.textentry(win, var=exampletransformbox, labels=example.transforms.labels, title="Example transformations",initval=example.transforms.functions,pos=c(2*w/3+spacer,spacer,w/3-spacer,main.height/2-spacer),
               action=function(win) {
                   return(win) })

  ## number of groups in final result
  ## to add

  ## quality of run? better quality is slower, which can be painful at the exploratory stage
  rp.radiogroup(win, quality, values=c('Exploratory run (faster)','Final run (better quality)'), title="Analysis type", pos=c(spacer,main.height/2+spacer, w/4-spacer, main.height/4-spacer))

  ## location
  rp.slider(win, lat.max,  from=-90, to=-30,  resolution=1,   title="North"   , initval=lat.max , showvalue=TRUE, pos=c(w/4+w/12+spacer,main.height/2+spacer,w/8-spacer,h))
  rp.slider(win, lon.min,  from=-180, to=180, resolution=1,   title="West"    , initval=lon.min, showvalue=TRUE, pos=c(w/4+spacer,main.height/2+spacer+h,w/8-spacer,h))
  rp.slider(win, lon.max,  from=-180, to=180, resolution=1,   title="East"    , initval=lon.max , showvalue=TRUE, pos=c(w/4+5*w/32+spacer,main.height/2+spacer+h,w/8-spacer,h))
  rp.slider(win, lat.min,  from=-90, to=-30,  resolution=1,   title="South"   , initval=lat.min , showvalue=TRUE, pos=c(w/4+w/12+spacer,main.height/2+spacer+2*h,w/8-spacer,h))

  rp.slider(win, n.groups,  from=2, to=40,  resolution=1,   title="Number of clusters"   , initval=12 , showvalue=TRUE, pos=c(3*w/4+spacer,main.height/2+spacer,w/4-spacer,h))

  ## hide output directory for now, since we aren't using it
  rp.text(win,txt=paste("Output directory:",output.dir),pos=c(spacer,main.height-2*h-2*spacer,w-spacer*2,h))

  rp.button(win,title="Run", pos=c(w-100,main.height-h-spacer,80,h), action=function(win, ...) {

      lat.step=0.1
      lon.step=0.1
      quality="low"
      if (grepl('^Final',win$quality)) {
          quality="high"
      }
      weights=lapply(win$weightbox,as.numeric)
#      cat(sprintf('gui weights:\n'))
#      cat(str(weights))

      b=bioreg(variables=selectedVariables,n.groups=win$n.groups,weights=weights, transformations=win$transformbox,lat.min=win$lat.min,lat.max=win$lat.max,lat.step=lat.step,lon.min=win$lon.min,lon.max=win$lon.max,lon.step=lon.step,quality=quality,output.dir=output.dir,...)
      return(win) })


  rp.button(win,title="Start again", pos=c(spacer,main.height-h-spacer,80,h), action=function(win) {
      tkdestroy(win$window)
      do.bioreg()
      return(win) })

  return(win)

}
