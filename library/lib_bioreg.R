#
#     Perform abiotic regionalisation
#
# (c) Copyright 2012 Ben Raymond, ben dot raymond at aad dot gov dot au
#     Last-Modified: <2012-05-03 13:23:16>
#
#-----------------------------------------------------------------------------

## Toolbox functions
#-----------------------------------------------------------------------------

get.bioreg.colourmap <- function(n=10) {
    #
    # Define some not-too-appalling colours to use
    #
    # n     number of colours on the scale
    #

    # base colour map (16 colours)
    cmap <- c("#C7D79EFF", "#FA9864FF", "#BEFFE8FF", "#D69DBCFF", "#F7ED59FF", "#A5F57AFF", "#FF3D4AFF", "#7AB6F5FF", "#369C5DFF", "#A80084FF", "#AA66CDFF", "#FFAA00FF", "#7AF5CAFF", "#FFBEBEFF", "#0070FFFF", "#E9FFBEFF")
    # repeat it if necessary (n > 16)
    cmap <- rep(cmap, ceiling(n/length(cmap)))
    # extract colours
    cmap <- cmap[1:n]

    return(cmap)
}


mcolor <- function(x, y=NULL, z=NULL, interp=FALSE, col=topo.colors(100), clim=NULL) {
    # plot a colour map
    # TODO deplace that with ggplot calls
    # uses rasterImage() as a faster alternative to image()

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
    temp=round((apply(t(tempz),2,rev)-min(tempz,na.rm=T))*(ncolours-1)/(max(tempz,na.rm=T)-min(tempz,na.rm=T)))+1
    tempa=as.raster(col[temp],nrow=dim(temp)[1])

    # show on figure
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
    # TODO why is that necessary? the function should get all its arguments passed to it?

    str <- str_c("{",str,"}") # make sure str is enclosed in curly brackets (will it matter if user also supplies these? - to check)

    # fill in the body of the function
    body(f) <- substitute(tryCatch(expr,
                                   error=function(e) "Error applying transformation function"),
                          list(expr=parse(text=str)[[1]])
    )

    return(f)
}


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


    # Check input arguments
    # variables
    # expand variable names
    variables <- list.env.data(variables, quiet=FALSE)
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
    # cat(deparse(transformations[[1]]))


    # Get and transform data
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

    data.transformed <- data.raw[,! names(data.raw) %in% c("lon", "lat")]
    for (i in 1:ncol(data.transformed)) {
        # apply user supplied transformations
        # TODO this is fragile : dependent on ordering of weights and transformations and data columns being the same: need to code it better
        # TODO: suggestion JO: use named lists/vectors for transformation and weights with names matching data columns? i.e. list(bathymetry="log(x)", floor_temperature=ceiling). But it's much more cumbersome to write then
        if (!is.null(transformations[[i]])) {
            data.transformed[,i] <- transformations[[i]](data.raw[,i])
        }

        # clean up any Inf values
        data.transformed[is.infinite(data.transformed[,i]),i] <- NA

        # normalise each column of x to 0-1 range
        # TODO: wouldn't scaling (0 mean, unit variance) be more appropriate?
        # data.transformed[,i] <- scale(data.transformed[,i])
        data.transformed[,i] <- data.transformed[,i] - min(data.transformed[,i], na.rm=T)
        data.transformed[,i] <- data.transformed[,i] / max(data.transformed[,i], na.rm=T)

        # apply weighting
        # cat(sprintf("i=%d, weights[i-2]: %s\n",i,str(weights[[i-2]])))
        data.transformed[,i] <- data.transformed[,i] * weights[i]
    }

    # record which lines (i.e. locations) are masked out because of missing data (including land)
    missing.mask <- rowSums(is.na(data.transformed)) > 0
    data.trans.noNA <- na.omit(data.transformed)

    message("-> Non-hierarchical clustering")
    # For the later hierarchical clustering we will need to compte a distance matrix between all data points. This is obviously impossible on the full data set, so we reduce the information to a smaller number of similar clusters through non-hierarchical clustering
    # number of clusters
    num.groups.intermediate <- 200

    # number of samples according to the quality argument (smaller numbers speed-up computation)
    samples <- switch(quality, low=5, high=50)

    # perform clustering
    cl <- clara(data.trans.noNA, k=num.groups.intermediate, metric="manhattan", stand=FALSE, samples=samples)

    # extract cluster numbers
    data.trans.noNA$clara.num <- cl$clustering


    message("-> Hierarchical clustering")
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


    message("-> Producing plots")

    # get a nice colour map
    cmap <- get.bioreg.colourmap(n.groups)

    # Dendrogram
    dev.new()
    # dendrogram
    plot(hcl, labels=F, hang=-1)
    # cutting level
    lines(c(1,num.groups.intermediate), c(temph,temph), lty=2, col="red")
    # markers for group labels
    colours <- cmap[hclust.num][hcl$order]
    points(1:200, rep(-0.02, num.groups.intermediate), col=NA, bg=colours, pch=21, cex=1)

    # Violin plot
    dev.new()
    datam <- melt(data.raw[,c(variables, "cluster")], id.vars="cluster")
    # ggplot(na.omit(datam)) + geom_violin(aes(x=cluster, y=value, fill=cluster), colour="black") + scale_fill_manual(values=cmap, guide="none") + coord_flip() + facet_wrap(~variable, scales="free")
    # boxplot for now
    ggplot(na.omit(datam)) + geom_boxplot(aes(x=cluster, y=value, fill=cluster), colour="black") + scale_fill_manual(values=cmap, guide="none") + coord_flip() + facet_wrap(~variable, scales="free")

    # Image map
    dev.new()
    polar.ggplot(data.raw, aes(colour=cluster))  + scale_colour_manual(values=cmap)

    # Output data
    if (!is.null(output.dir)) {
        output.file <- normalizePath(str_c(output.dir,"/bioreg.Rdata"), winslash="/", mustWork=F)
        save(data.raw, data.transformed, file=output.file)
    }
    invisible(data.raw) # invisible() so that it doesn't get printed to console if not assigned to a variable
}
