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


## Analyse GDM output
#-----------------------------------------------------------------------------

print.gdm <- function(x, ...) {
  suppressPackageStartupMessages(require("stringr", quietly=TRUE))

  species <- setdiff(names(x$data), c("lon", "lat", x$model$predictors))

  cat("\n     GDM model\n\n")

  cat("Species :\n")
  cat(" ", str_c(species, collapse="\n  "), "\n")
  cat("Predictors :\n")
  cat(" ", str_c(x$model$predictors, collapse="\n  "), "\n")

}

summary.gdm <- function(x, ...) {
  gdm.summary(x$model)
}


plot.gdm <- function(x, ...) {
  #
  # Plot "effects" in a GDM model
  #
  # x   object of class gdm
  #
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))
  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
  suppressPackageStartupMessages(require("reshape2", quietly=TRUE))

  # compute environmental data range (omitting lat and lon)
  env.ranges <- llply(x$data[, x$model$predictors], function(x) {
    seq(from=min(x, na.rm=T), to=max(x, na.rm=T), length.out=200)
  })
  env.ranges <- data.frame(do.call(cbind, env.ranges))

  # compute the transformation induced by the model on these ranges
  env.curves <- gdm.transform(x$model, env.ranges)

  # match original ranges to the transformations
  env.ranges.m <- melt(env.ranges, measure.vars=names(env.curves), value.name="value")
  env.curves.m <- melt(env.curves, measure.vars=names(env.curves), value.name="transform")
  # and join them
  env <- data.frame(env.curves.m, env.ranges.m)

  # compute quantiles of the original observations, to plot them as a rug plot
  quant <- data.frame(llply(x$data[,x$model$predictors], function(x) {
    # to easily mix them with the others, compute 200 quantiles
    as.numeric(quantile(x, probs=seq(0, 1, length.out=200)))
  }))
  quant.m <- melt(quant, measure.vars=names(quant), value.name="quantiles")
  # put this with the rest of the data
  env$quantiles <- quant.m$quantiles

  # sort the variables in decreasing order of transform value
  maxVal <- ldply(env.curves, max, na.rm=T)
  maxVal <- maxVal[order(maxVal$V1, decreasing=TRUE),]
  # reorder factor for display
  env$variable <- factor(env$variable, levels=maxVal$.id)

  # plot
  p <- ggplot(env) +
    # the transform result
    geom_path(aes(x=value, y=transform)) +
    # the distribution of data
    geom_rug(aes(x=quantiles), alpha=0.2) +
    # for each variable
    facet_wrap(~variable, scales="free_x")

  return(p)
}

plot.pred.gdm <- function(x, ...) {
  polar.ggplot(x$prediction, aes(fill=cluster), ...)
}


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
