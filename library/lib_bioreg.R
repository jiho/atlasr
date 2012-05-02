#
#   Perform abiotic regionalisation
#
# (c) Copyright 2012 Ben Raymond, ben dot raymond at aad dot gov dot au
#     Last-Modified: <2012-05-02 12:33:41>
#
#-----------------------------------------------------------------------------


# some functions to make things easier
get.bioreg.colourmap=function(n=10) {
    # define some not-too-appalling colours to use
    cmap=c("#C7D79EFF","#FA9864FF","#BEFFE8FF","#D69DBCFF","#F7ED59FF","#A5F57AFF","#FF3D4AFF","#7AB6F5FF","#369C5DFF","#A80084FF","#AA66CDFF","#FFAA00FF","#7AF5CAFF","#FFBEBEFF","#0070FFFF","#E9FFBEFF")
    cmap=rep(cmap,ceiling(n/length(cmap)))
    cmap=cmap[1:n]
    return(cmap)
}

mcolor=function(x,y=NULL,z=NULL,interp=F,col=topo.colors(100),clim=NULL) {
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
  f <- function(x) {}
  environment(f) <- baseenv()
  str=paste('{',str,'}',sep='') # make sure str is enclosed in curly brackets (will it matter if user also supplies these? - to check)
  body(f) <- substitute(tryCatch(expr,
                                 error=function(e) "Error applying transformation function"),
                        list(expr=parse(text=str)[[1]]))
  f
}


bioreg <- function(variables, n.groups=12, lat.min=-80, lat.max=-30, lat.step=0.1, lon.min=-180, lon.max=180, lon.step=0.5, transformations=NULL, weights=NULL, path="env_data", quality=c("low","high"), output.dir=NULL)
{
    # n.groups: either an integer number of groups, or the height at which to cut the dendrogram (e.g. 0.13)
    # quality: "low" or "high" - low quality is faster, suitable for exploratory runs. Use high quality for final analyses.
    # weights: list giving the weight for each variable
    # transformations: list giving transformation function for each variable, or NULL for no transformations for any
    # output.dir: destination for output files - if NULL, no output files will be saved

    suppressPackageStartupMessages(require("cluster", quietly=TRUE))
    suppressPackageStartupMessages(require("vegan", quietly=TRUE))

    if (length(variables)<2) {
        stop('You must specify at least two input variables')
    }
    quality=match.arg(quality)


    # weights
    if (is.null(weights)) {
        # default to equal (unit) weighting
        weights=rep(1,length(variables))
    }
    weights=as.numeric(weights)
    #cat(str(weights))
    weights=weights/max(weights) # normalize so that max weight is 1

    # transformations
    if (!is.null(transformations)) {
        # convert each from string expression into actual function, if necessary
        tfuncs=list()
        for (i in (1:length(transformations))) {
            if (is.character(transformations[[i]]) & nchar(transformations[[i]])>0) {
                tryCatch(tfuncs[[i]]<-function.maker(transformations[[i]]),
                         error=function(e) {stop(sprintf('Error parsing transformation function \"%s\"',transformations[[i]])) }
                         )
            } else if (is.function(transformations[[i]])) {
                tfuncs[[i]]=transformations[[i]]
            } else if (is.null(transformations[[i]]) | (is.character(transformations[[i]]) & nchar(transformations[[i]])==0)) {
                tfuncs[[i]]=function(x){x}
            } else {
                stop(sprintf('Supplied transformation is neither a string nor a function'))
            }
        }
        transformations=tfuncs
    }
#    cat(deparse(transformations[[1]]))

    num.groups.intermediate=200 #number of clusters to produce in the non-hierarchical clustering step

    database=read.env.data(variables=variables, path=path)
    prediction_grid <- build.grid(
                          lat.min=lat.min, lat.max=lat.max, lat.step=lat.step,
                          lon.min=lon.min, lon.max=lon.max, lon.step=lon.step
    )
    data.raw <- associate.env.data(prediction_grid, database)


    data.transformed=data.raw
    for (i in 3:dim(data.transformed)[2]) { # skip first two cols: those are lon and lat
        # ** this is fragile (dependent on ordering of data with first two cols lon and lat; dependent on ordering of weights and transformations and data columns being the same): need to code it better
        if (!is.null(transformations[[i-2]])) {
            data.transformed[,i]=transformations[[i-2]](data.raw[,i])
        }
        # clean up any Inf values
        data.transformed[which(is.infinite(data.transformed[,i])),i]=NA

        # normalise each column of x to 0-1 range (but not lon or lat)
        data.transformed[,i]=data.transformed[,i]-min(data.transformed[,i],na.rm=T)
        data.transformed[,i]=data.transformed[,i]/(max(data.transformed[,i],na.rm=T))

        # apply weighting
#        cat(sprintf("i=%d, weights[i-2]: %s\n",i,str(weights[[i-2]])))
        data.transformed[,i]=data.transformed[,i]*weights[[i-2]]
    }

    # which records to mask out because of missing data (including land)
    missing.mask=which(rowSums(is.na(data.transformed))>0)
    not.missing.mask=which(rowSums(is.na(data.transformed))==0)

    # which columns are the data columns (i.e. not lon or lat)
    datcols=which(!(names(data.transformed) %in% c('lon','lat')))

    message(sprintf('-> Non-hierarchical clustering ... ')); flush.console()

    # non-hierarchical clustering step
    if (quality=="low") {
        cl=clara(data.transformed[not.missing.mask,datcols],num.groups.intermediate,metric="manhattan",stand=FALSE,samples=5)
    } else {
        cl=clara(data.transformed[not.missing.mask,datcols],num.groups.intermediate,metric="manhattan",stand=FALSE,samples=50)
    }
    cluster.num=cl$clustering

    message(sprintf('-> Hierarchical clustering ... ')); flush.console()

    # now do a hierarchical clustering using the output of the nonhierarchical step
    # first calculate mean properties of the nonhierarchical clusters
    xc=matrix(NA,nrow=num.groups.intermediate,ncol=length(datcols))
    u.cluster.num=unique(cluster.num)

    for (k in 1:length(u.cluster.num)) {
        tempidx=which(cluster.num==u.cluster.num[k])
        xc[k,]=colMeans(data.transformed[not.missing.mask[tempidx],datcols],na.rm=T)
    }

    # dissimilarities of these clusters
    D=vegdist(xc,method="gower")

    # hierarchical clustering
    hcl=hclust(D,method="ave")

    # now extract the desired number of groups from the dendrogram
    if (floor(n.groups)==n.groups) {
        # we specified a number of groups directly
        cn.new=cutree(hcl,k=n.groups)
        # work out the dissimilarity level (height) that corresponds to this number of groups
        temph=mean(c(hcl$height[length(hcl$height)+2-n.groups],hcl$height[length(hcl$height)+2-n.groups-1]))
    } else {
        # we specified a height at which to cut the dendrogram
        # show on the dendrogram the height at which we are cutting
        temph=n.groups
        cn.new=cutree(hcl,h=n.groups)
        n.groups=length(unique(cn.new))
    }

    message(sprintf('-> Producing plots ... ')); flush.console()
    cmap=get.bioreg.colourmap(n.groups)
    dev.new()
    plot(hcl,labels=F,hang=-1)
    lines(c(1,num.groups.intermediate),c(temph,temph),lty=2,col=2)
    # add markers for group labels
    dorder=order.dendrogram(as.dendrogram(hcl))
    for (k in 1:n.groups) {
        temp=which(cn.new[dorder]==k)
        points(temp,rep(-0.02,length(temp)),col=cmap[k],bg=cmap[k],pch=21,cex=5)
    }

    cluster.num.new=integer(length=dim(data.transformed)[1])
    cluster.num.new[missing.mask]=NA
    for (k in 1:length(u.cluster.num)) {
        tempidx=which(cluster.num==u.cluster.num[k])
        cluster.num.new[not.missing.mask[tempidx]]=cn.new[k]
    }
    temp=prediction_grid
    temp$cluster.num=cluster.num.new

    dev.new()
    mcolor(matrix(temp$lon,nrow=attr(temp,'out.attrs')$dim[1]),matrix(temp$lat,nrow=attr(temp,'out.attrs')$dim[1]),matrix(temp$cluster.num,nrow=attr(temp,'out.attrs')$dim[1]),col=cmap)
#    image(matrix(temp$lon,nrow=attr(temp,'out.attrs')$dim[1])[,1],matrix(temp$lat,nrow=attr(temp,'out.attrs')$dim[1])[1,],matrix(temp$cluster.num,nrow=attr(temp,'out.attrs')$dim[1]),col=cmap)

    # collate environmental (raw) data per cluster, and produce boxplot
    data.raw$cluster=cluster.num.new
    dev.new()
    par(mfcol=c(1,length(datcols)))
    for (k in datcols) {
        boxplot(as.formula(sprintf('%s ~ cluster',names(data.raw)[k])),data=data.raw,col=cmap,xlab=names(data.raw[k]),horizontal=T)
    }
    invisible(data.raw) # invisible() so that it doesn't get printed to console if not assigned to a variable
    if (!is.null(output.dir)) {
        output.file=normalizePath(paste(output.dir,'/','bioreg.Rdata',sep=''),winslash='/',mustWork=F)
        save(data.raw,data.transformed,file=output.file)
    }
}
