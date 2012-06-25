#
#     General Dissimilarity Modelling
#
# (c) Copyright 2011-2012 S Mormede, Ben Raymond, J-O Irisson
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------

# Ben Raymond, 2011:
# note that the training data cluster labels and the indicator species info uses simple spatial interpolation, and will not work with time-varying models


## Run GDM analysis
#-----------------------------------------------------------------------------

compute.gdm <- function(
  resp.vars,          # names of the response variables (taxa)
  pred.vars,          # names of the predictor variables (environment variables)
  data,               # data frame containing locations (lat, lon), presence/absence data and environment variables
  newdata,            # prediction grid coordinates and associated environment variables
  n.groups=NULL,      # number of clusters in the result; when NULL (the default) the optimal number of clusters is determined by a stepwise procedure
  pre.sample=2500,    # sub-sample the input data before feeding it to the GDM function. When NULL, no subsampling occurs. Subsampling reduces the number of pairwise operations needed, speeds up the function and determines the amount of memory required. Since the library is 32 bits, the amount of allocatable memory is limited and this should not be much larger than 2500
  intern.sample=NULL, # internally, the GDM function can also subsample the pairwise dissimilarities. When NULL, no subsampling occurs. The total number of pairwise dissimilarities is n * ( n - 1 ) / 2, where n is the number of rows in `data` (or `pre.sample` when subsampling a priori), so for this to be effective, intern.sample must be lower than that. As a rule of thumb, n=2500 gives over 3 million pairwise dissimilarities.
  ...                 # passed to the internal GDM function
)
{
  ## Run GDM
  #--------------------------------------------------------------------------
  message("-> Compute dissimilarities")

  # set the random number generator seed for replicable results
  set.seed(123)

  # remove lines with missing data
  data <- data[,c("lon", "lat", resp.vars, pred.vars)]
  data <- na.omit(data)

  # subsample input data ( and remove lon and lat)
  if (!is.null(pre.sample)) {
    sdata <- data[sample.int(nrow(data), pre.sample), -c(1,2)]
  } else {
    sdata <- data[,-c(1,2)]
  }

  # run GDM
  m.gdm <- gdm.fit(sdata[,pred.vars], sdata[,resp.vars], sample=intern.sample, ...)


  ## Cluster the data
  #--------------------------------------------------------------------------
  message("-> Cluster based on dissimilarities")
  suppressPackageStartupMessages(require("cluster", quietly=TRUE))

  # remove NAs from prediction data
  newdata <- na.omit(newdata[, c("lat", "lon", pred.vars)])

  # predict the environment according to the GDM model
  # NB: this is done by pieces of 1000 to limit the amount of memory required
  pred <- newdata[1, pred.vars]
  for (i in seq(from=1,to=nrow(newdata)%/%1000*1000+1, by=1000)) {
    pred <- rbind(pred, gdm.transform(m.gdm, newdata[i:min(i+999,nrow(newdata)), pred.vars]))
  }
  pred<-pred[-1,]
  # TODO re-implement the 1000 split more efficiently

  if ( is.null(n.groups) ) {
    # optimise the number of clusters
    # prepare storage
    res <- array(NA,c(10,8,4))

    # try several possibilities
    for (cl in 4:10) {
      for (samples in 3:8) {
        for (sampsz in 2:4) {
          temp <- clara(pred, k=cl, metric="manhattan", keep.data=FALSE, samples=samples, sampsize=min(nrow(pred), 40 + sampsz * cl))
          res[cl,samples,sampsz] <- temp$silinfo$avg.width
        }
      }
    }

    # select the best one(s)
    temp <- which(res==max(res,na.rm=T), arr.ind=T)
    # when there are several, choose the first (least number of clusters)
    if (nrow(temp)>1) temp<-temp[1,]

    # recompute clustering with this optimal parameters
    cl <- clara(pred, k=temp[1], metric="manhattan", keep.data=FALSE, samples=temp[2], sampsize=min(nrow(pred), 40 +temp[3]*temp[1]))

  } else {
    # compute the clustering with the given number of groups
    cl <- clara(pred, k=n.groups, metric="manhattan", keep.data=FALSE)
  }


  ## Return the result
  #-----------------------------------------------------------------------------

  # store the cluster number in the data
  newdata$cluster <- factor(cl$clustering)

  # store elements in the resulting object
  res <- list(
    data=data,
    model=m.gdm,
    cl=cl,
    prediction=newdata
  )

  # class it as gdm
  class(res) <- c("gdm", "list")

  return(invisible(res))
}

gdm <- function(
  file,               # name of the file where the presence/abundance data is
  taxa="",            # names (or abbreviations) of the taxa of interest (by default, all columns in file except for lat and lon)
  variables,          # names (or abbreviations) of the variables to use for the prediction
  lat.min=-80, lat.max=-30, lat.step=0.1,   # definition of the prediction grid
  lon.min=-180, lon.max=180, lon.step=0.5,
  path=getOption("atlasr.env.data"),        # path to the environmental database
  ...                 # passed to compute.brt()
)
{

  # read selected variables from the database
  database <- read.env.data(variables, path=path, quiet=FALSE)
  # get the full names of those variables
  pred.vars <- names(database)

  # remove information on land
  database <- mask.env.data(database, path=path)


  # read input dataset
  file <- clean.path(file)
  if (file.exists(file)) {
    input.data <- read.data(file)
  } else {
    stop("Cannot find file : ", file)
  }

  # get the names of the taxa of interest
  allTaxa <- names(input.data[,!names(input.data) %in% c("lat", "lon")])
  resp.vars <- match.vars(taxa, allTaxa)

  # get environment data for the observations
  input.data <- associate.env.data(input.data, database)

  # build prediction grid
  prediction.data <- build.grid(
    lat.min=lat.min, lat.max=lat.max, lat.step=lat.step,
    lon.min=lon.min, lon.max=lon.max, lon.step=lon.step
  )
  # get environment data on this grid
  prediction.data <- associate.env.data(prediction.data, database)

  # compute the GDM model and clustering
  gdmObject <- compute.gdm(resp.vars=resp.vars, pred.vars=pred.vars, data=input.data, newdata=prediction.data, ...)

  return(invisible(gdmObject))
}


## Plots
#-----------------------------------------------------------------------------




plot.pred.gdm <- function(x, ...) {
  polar.ggplot(x$prediction, aes(fill=cluster), ...)
}

## and use gdm.transform to create the curves

## make environmental ranges to plot response curves

# env.ranges<-dat[1:200,]
# for (i in pred.var.col) {
#     env.ranges[,i]<-seq(from=(Min(dat[,i])),to=Max((dat[,i])),length=200)
# }
# env.ranges<-env.ranges[,pred.var.col]
# no.curves <- gdm.transform(no.gdm,env.ranges)
#
# ## sort them by max to min
# tp<-names(rev(sort(apply(no.curves,2,max))))
# temp<-match(tp,names(no.curves))
#
#
# ## then plot them out
# par(mfrow=c(3,4))
# par(cex=0.9)
# j<-1
# for (i in temp) {
#     if(names(no.curves)[i] %in% names(lookup.names.variables)) {
#         tp<-as.character(lookup.names.variables[names(no.curves)[i]])
#     } else {
#         tp<-names(no.curves)[i]
#     }
#     plot(env.ranges[,i],no.curves[,i],type='l',cex=1.5,
#          xlab=tp, ylab = "transform")
#     rug(quantile(dat[,names(no.curves)[i]], probs = seq(0, 1, 0.1), na.rm = TRUE))
# }
#


# ## plot the results
# p = polar.ggplot(newdata, aes(colour=cluster), geom="point")
# # add a title
# p = p + opts(title="GDM")
# # display the plot
# print(p)
# # NB: suppress warnings about missing values: they are necessary to split the coastline in several bits

# ## calculate indicator species using dufrene-legendre method
# if (do.indicator.species) {
#     require(labdsv)
#     library(reshape)
#     temp <- as.matrix(cast(newdata, lat ~  long, value = "cluster",add.missing=T))
#     datlonidx <- round(approx(as.numeric(colnames(temp)),1:ncol(temp),dat$long)$y)
#     datlatidx <- round(approx(as.numeric(rownames(temp)),1:nrow(temp),dat$lat)$y)
#
#     dat$cluster <- NA
#     for (i in 1:length(datlonidx)){
#       if (!is.na(datlonidx) && !is.na(datlatidx)) {
#         dat[i,"cluster"] <- temp[datlatidx[i],datlonidx[i]]
#       }
#     }
#
#     nnanidx=which(!is.na(dat$cluster)) ## only include non-NA clusters
#     ## calculate indicator species stuff
#     frodo.baggins=indval(dat[nnanidx,resp.var.col],clustering=dat$cluster[nnanidx])
# }

# newdata$cluster<-factor(newdata$cluster)
# dat$cluster <- factor(dat$cluster)


# close the PDF file
# if (!is.null(plot.name)) {
#   dev.off()
# }
#
#   if (do.indicator.species) {
#     res=list('predicted'=newdata,'model'=no.gdm,'plot.obj'=p,'dat.cluster'=dat$cluster,'indval'=frodo.baggins)
#   } else {
