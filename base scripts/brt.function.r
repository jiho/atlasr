###########################################################################
#
# Wrap around BRT function
# does model, effects, bootstraps, predictions in one function
#
# Sophie Mormede May 2011 s.mormede@niwa.co.nz
#
###########################################################################

do.brt <- function (dat, resp.vars, predvar, int = 2, distrib = "bernoulli", wghts = NULL, monotone = NULL, n.boot = NA, plot.name = NULL, n.pred = NA, pred.data, ...) {

# Arguments
#
# dat         dataset containing coords, species presence, environmental data
# resp.vars   names or indexes of response variables (species to model)
# predvar     names of the predictive variable (environmental data)
# int         number of interactions, or tree complexity
# distrib     distribution family: bernoulli (= binomial), poisson,
#             laplace or gaussian
# wgths       name of the variable in `dat` which corresponds to weights
# monotone    a vector of 0 and 1 of the lengths of predvar
#             1 is monotone increasing
#             -1 monotone decreasing
#             0 is arbitrary
# n.boot      number of bootstraps for the predictions, NA means single
#             point estimate. Beware, bootstraps are slow to run and
#             less than 50-100 bootstraps will fail.
# plot.name   prefix for the name of plots
#             can be a relative path or full path
#             plots will be saved as PDF files, one for each response variable
#             if NULL plots are only displayed to screen and not saved
# n.pred      number of bootstraps for predictions
#             1 if no bootstraps wanted and NA if no predictions wanted
# pred.data   predictive dataframe, needs lat and long and environmental
#             variables defined in predvar

# Value
#
# It will output a list of objects containing, for each response variable: the brt object, the responses, the interactions, the predictions, etc.

# get useful functions
suppressMessages(require("gbm", quietly=TRUE))

result <- list()

# loop on response variables
# BRT only deals with one species at a time
for (resp in resp.vars) {

  # create a new PDF file for the plots
  if (!is.null(plot.name)) {
    pdf(paste(plot.name, "-", resp, ".pdf" ,sep=""), width=10, height=8, useDingbats=FALSE, onefile=TRUE)
  }

  # start the output object
  result[[resp]] <- list()

  ###########################################################################
  # setup the dataset
  ###########################################################################
  # remove observations where the current species information is Non-Available
  dat <- dat[!is.na(dat[,resp]),]

  # make the data binomial if it needs to be
  if (distrib=="bernoulli") { dat[,resp] <- dat[,resp]>0 }

  # make uniform weights when they are not provided
  if (is.null(wghts)) {
    wghts <- rep(1,nrow(dat))
  } else {
    wgths <- dat$wghts
  }

  # make indifferent monotone argument when it is not provided
  if (is.null(monotone)) { monotone <- rep(0,length(predvar)) }

  ###########################################################################
  # run the BRT
  ###########################################################################
  cat(paste("  ", resp,":\n  -> optimising BRT model ",sep=""))

  # initial values before entering the while loop
  lr <- 0.05
  no.trees <- 0

  while ( no.trees < 1000 & lr > 0.0005 ) {
    cat(".")
    try( obj <-  gbm.step(
                           data = dat,
                           gbm.x = match(predvar, names(dat)),
                           gbm.y = match(resp, names(dat)),
                           learning.rate = lr,
                           tree.complexity = int,
                           site.weights = wghts,
                           family = distrib,
                           max.trees = 10000,
                           var.monotone = monotone,
                           n.trees = 50,
                           silent = T,
                           plot.main = F
                         )
    )

    # if the gbm does not converge, the return object is null or of size 0
    if ( ! is.null(obj) ) {
      if ( object.size(obj) > 0 ) {
        no.trees <- obj$gbm.call$best.trees
      }
    } else {
      no.trees <- 0
    }

    # decrease the learning rate
    lr <- lr / 2
  }
  cat("\n")

  # store the gbm object in the result object
  result[[resp]]$obj <- obj


  ###########################################################################
  # write additional results in the result object
  ###########################################################################
  cat("  -> writing results\n")

  temp <- list()
  # family
  temp$family <- obj$gbm.call$family
  # number of interaction = depth of tree
  temp$tree.complexity <- obj$gbm.call$tree.complexity
  # deviance explained
  temp$perc.deviance.explained <- base::round( obj$cv.statistics$deviance.mean / obj$self.statistics$mean.null, 2)
  # Area Under the receiver operative Curve
  # = quality of the prediction of presence/absence
  #   0.5 is indifferent
  temp$AUC <- round(min(obj$cv.roc.matrix), 2)
  result[[resp]]$deviance <- temp

  # write contributions
  temp <- obj$contributions$rel.inf
  names(temp) <- obj$contributions$var
  result[[resp]]$contributions <- temp

  # remove the temp object, to be clean
  rm(temp)

  ###########################################################################
  # plot effects, bootstrapped or not
  ###########################################################################
  cat("  -> plot effects\n")

  if (is.na(n.boot)) {
    # no bootstrap, just plot
    gbm.plot(obj)

  } else {
    # perform bootstrap
    cat("  -> running bootstrap on BRT model\n")
    boot <- NULL
    try(
      boot <- gbm.bootstrap(obj, n.reps=n.boot, verbose=F)
    )
    if (!is.null(boot)) {
      # when it runs correctly
      # store the output in the result object
      result[[resp]]$boot <- boot
      # and use bootstrapped plots
      gbm.plot.boot(obj, boot)

    } else {
      # when it does not (usually because of too few bootstraps), issue a warning and plot the simple object
      warning("Need more bootstraps to provide a confidence interval\n")
      gbm.plot(obj)
    }
  }

  ###########################################################################
  # run the predictions
  ###########################################################################

  if (!is.na(n.pred)) {
    # when n.pred is not NA, do the predictions
    cat("  -> make predictions\n")

    if (n.pred > 1) {
      # bootstrap predictions
      temp <- gbm.bootstrap(obj, pred.data = pred.data, return.pred.matrix = T, n.reps=n.pred, verbose=F)
      temp <- temp$pred.matrix
      # prediction and confidence zone
      pred = apply(temp,1,mean)
      CVpred = apply(temp,1,cv)

    } else {
      # do a simple prediction
      pred = predict(obj, newdata = pred.data, n.trees = obj$n.trees, type="response")
      CVpred = NA # CV can't be computed without bootstrap
    }
    # store it in the result object
    result[[resp]]$pred <- data.frame(pred.data[,c("lat", "long")], pred, CVpred)

    cat("  -> plot predictions\n")
    # plot the predictions
    if (all(is.na(CVpred))) {
      # if all uncertainty on prediction are NAs we probably did not do bootstrap and we are interested in plotting things quickly
      # => we use points
      p = polar.ggplot(result[[resp]]$pred, aes(colour=pred), geom="point")
    } else {
      # if we have uncertainty data we want to represent it (as transparency)
      # => we use tiles
      p = polar.ggplot(result[[resp]]$pred, aes(fill=pred, alpha=-CVpred), geom="tile")
    }
    # add a title
    p = p + opts(title=paste(resp, "- BRT"))
    # store the plot object in the result
    result$pred.plot = p
    # display the plot
    suppressWarnings(print(p))
    # NB: suppress warnings about missing values: they are necessary to split the coastline in several bits
  }

  # close the PDF file
  if (!is.null(plot.name)) {
    dev.off()
  }

}  # end of resp.vars loop

return(result)

} # end of function



############################################################################
# Some general functions
############################################################################

Sum <- function(x) {sum(x,na.rm=T)}
Max <-function(x) {max(x,na.rm=T)}
Min <- function(x) {min(x,na.rm=T)}
Median <- function(x) {median(x,na.rm=T)}
Mean <- function(x) {mean(x,na.rm=T)}
Length <- function(x) length(x[!is.na(x)])
Table <- function(...) {
  a <-table(...,useNA="ifany")
  b <- length(dim(a))
  if (b==2) {
    a<-cbind(a,Sum = margin.table(a,1))
  } else {
    a<-c(a,Sum = margin.table(a))
  }
  return(a)
}
