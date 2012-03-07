#
#     Boosted Regression Trees
#
#     - bootstrapping functions (complement those in package dismo)
#     - wrap around BRT function which does model, effects, bootstraps, predictions
#     - fully automated function
#
# (c) Copyright 2011-2012 S Mormede, J-O Irisson
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------



## Functions modified / absent from the dismo package
#-----------------------------------------------------------------------------

gbm.plot.boot <- function(
     gbm.object,                    # a gbm object - could be one from gbm.step
     boot.object,                   # a gbm bootstrapped object - could be one from bootstrap
     variable.no = 0,               # the var to plot - if zero then plots all
     smooth = FALSE,                # should we add a smoothed version of the fitted function
     rug = T,                       # plot a rug of deciles
     n.plots = length(pred.names),  # plot the first n most important preds
     common.scale = T,              # use a common scale on the y axis
     write.title = T,               # plot a title above the plot
     y.label = "fitted function",   # the default y-axis label
     x.label = NULL,                # the default x-axis label
     show.contrib = T,              # show the contribution on the x axis
     plot.layout = c(3,4),          # define the default layout for graphs on the page
     ...                            # other arguments to pass to the plotting
                                    # useful options include cex.axis, cex.lab, etc.
)
{
# function to plot gbm bootstrapped response variables, with the option
# of adding a smooth representation of the response if requested
# additional options in this version allow for plotting on a common scale
# note too that fitted functions are now centered by subtracting their mean
#
# version 1.0
#
# S Mormede Feb 2010
#

  require(gbm)
  require(splines)

  #cleaning us start stuff

  gbm.call <- boot.object$gbm.call
  gbm.x <- gbm.call$gbm.x
  pred.names <- gbm.call$predictor.names
  response.name <- gbm.call$response.name
  data <- boot.object$x
  boot<-boot.object$function.dataframe

  max.plots <- plot.layout[1] * plot.layout[2]
  plot.count <- 0
  n.pages <- 1

  if (length(variable.no) > 1) {stop("only one response variable can be plotted at a time")}

  if (variable.no > 0) {
    # we are plotting all vars in rank order of contribution
    n.plots <- 1
  }

  max.vars <- length(gbm.object$contributions$var)
  if (n.plots > max.vars) {
    n.plots <- max.vars
    cat("warning - reducing no of plotted predictors to maximum available (",max.vars,")\n",sep="")
  }


  # get the variable names here and min and max
  if (n.plots == 1) {
    varname <- as.character(pred.names[variable.no])
  } else varname <- as.character(gbm.object$contributions$var)

  name.vals<-paste(varname,".vals",sep="")
  name.mean<-paste(varname,".mean",sep="")
  name.lower<-paste(varname,".lower",sep="")
  name.upper<-paste(varname,".upper",sep="")

  ymin <- min (boot[,names(boot) %in% name.lower])
  ymax <- max (boot[names(boot) %in% name.upper])

  # now do the actual plots

  for (j in c(1:length(varname))) {

    if (plot.count == max.plots) {
      plot.count = 0
      n.pages <- n.pages + 1
    }

    if (plot.count == 0) {
      par(mfrow = plot.layout)
      # NB: this also creates a new page when there is none
    }

    plot.count <- plot.count + 1

    if (show.contrib) {
      x.label <- paste (varname[j],"  (",round(gbm.object$contributions[gbm.object$contributions$var == varname[j],2],1),"%)",sep="")
    } else {
      x.label <- varname[j]
    }

    if (!common.scale) {
      ymin <-min(boot[,name.lower[j]])
      ymax <-max(boot[,name.upper[j]])
    }

    if (is.factor(data[,varname[j]])) {
      num<-length(unique(data[,varname[j]]))
      plot(boot[1:num,name.vals[j]],boot[1:num,name.mean[j]],ylim=c(ymin,ymax), type="n", xlab = x.label, ylab = y.label)
      lines(boot[1:num,name.vals[j]],boot[c(1:num),name.lower[j]],col=8)
      lines(boot[1:num,name.vals[j]],boot[c(1:num),name.upper[j]],col=8)

    } else {
       plot(boot[,name.vals[j]],boot[,name.mean[j]],ylim=c(ymin,ymax), type='l', xlab = x.label, ylab = y.label)
       lines(boot[,name.vals[j]],boot[,name.lower[j]],col=8)
       lines(boot[,name.vals[j]],boot[,name.upper[j]],col=8)
    }


   if (smooth & is.vector(boot[,name.vals[j]]) ) {
      temp.lo <- loess(boot[,name.mean[j]] ~ boot[,name.vals[j]], span = 0.3)
      lines(boot[,name.vals[j]],fitted(temp.lo), lty = 2, col = 2)
   }
   if (plot.count == 1) {
     if (write.title) {
        title(paste(response.name," - page ",n.pages,sep=""))
     }
   } else {
      if (write.title & j == 1) {
        title(response.name)
      }
   }

    k <- match(varname[j],names(data))
    if (rug & is.vector(data[,k])) {
      rug(quantile(data[,k], probs = seq(0, 1, 0.1), na.rm = TRUE),lwd=2)
    }
  }
}


gbm.bootstrap <- function(
  gbm.object,                   # a gbm object describing sample intensity
  bootstrap.model = TRUE,       # bootstrap the model fitted and predicted values using 632+
  bootstrap.functions = TRUE,   # estimate confidence intervals for the fitted functions
  pred.data = NULL,             # an independent evaluation data set - leave here for now
  return.train.matrix = FALSE,  # return the full matrix of predictions
  return.pred.matrix = FALSE,   # return the full matrix of predictions
  return.tails = FALSE,         # return the tails of the predictions
  CI = 95,                      # the required confidence interval - should be 90, 95, 99
  n.reps = 200,       # number of bootstrap samples
  verbose = T                   # control reporting
)
{
# function to calculate bootstrap estimates of 95% confidence intervals for
# a brt model created using gbm.step
#
# version 2.9 - J. Leathwick/J. Elith - June 2007
#
# function takes an initial BRT model fitted with gbm.step and uses bootstrapping
# to estimate model uncertainty either in the model fitted values, the model functions
# or in predictions to new data
# Because model fitted values are only estimated for a subset of sites
# for each bootstrap sample, n reps will typically produce n * 0.632 estimates
# because this number fluctuates, we calculate the 95% confidence intervals
# for these using interpolation and dropping the "missing" values
# we also use interpolation to establish the 95% confidence limits for the fitted
# functions, allowing these to be estimate from smaller numbers of repetitions
# To allow bootstrap estimation at new predicted sites for large datasets,
# only the upper and lower tails of the predicted values are saved

require(gbm)

# first get the original analysis details..

  gbm.call <- gbm.object$gbm.call
  train.data <- gbm.object$x
  n.obs <- nrow(train.data)
  gbm.x <- gbm.call$gbm.x
  gbm.y <- gbm.call$gbm.y
  family <- gbm.call$family
  lr <- gbm.call$learning.rate
  tc <- gbm.call$tree.complexity
  response.name <- gbm.call$response.name
  predictor.names <- gbm.call$predictor.names
  n.preds <- length(gbm.x)
  n.trees <- gbm.call$best.trees
  weights <- gbm.object$weights

# now setup storage space for the model bootstrapping if required

  if (bootstrap.model) {

    train.fit.matrix <- matrix(NA, ncol = n.reps, nrow = n.obs)
    train.pred.matrix <- matrix(NA, ncol = n.reps, nrow = n.obs)

  }

  if (bootstrap.functions) {
    # then set up a dataframe of means for predicting them

    function.matrix <- array(0,dim=c(n.reps,n.preds,200))

    function.pred.frame <- as.data.frame(as.vector(rep(mean(train.data[,gbm.x[1]],na.rm=T),200)))

    for (i in 2:n.preds) {
      # step through the predictor set

      j <- gbm.x[i]

      if (is.vector(train.data[,j])) {
        temp.vector <- as.vector(rep(mean(train.data[,j],na.rm=T),200))
        function.pred.frame <- cbind(function.pred.frame,temp.vector)
      }
      if (is.factor(train.data[,j])) {
        temp.table <- sort(table(train.data[,j]),decreasing = TRUE)
        temp.factor <- factor(rep(names(temp.table)[1],200),levels = levels(train.data[,j]))
        function.pred.frame <- cbind(function.pred.frame,temp.factor)
      }
    }
    names(function.pred.frame) <- predictor.names
  }


# then check and see if a prediction dataframe is provided...

  if (is.null(pred.data)) {

    bootstrap.predictions <- FALSE

  } else {
    # we have a prediction dataset

    bootstrap.predictions <- TRUE

    function.pred.frame.name <- deparse(substitute(pred.data))   # get the dataframe name
    pred.data <- eval(pred.data)
    n.pred.obs <- nrow(pred.data)

# and setup storage space for the prediction data
# first calculating the number of columns needed for the tail matrices
# by first getting half of 100 - CI, multipling by n.reps, and taking the ceiling
# then add one to take the values added for each bootstrap model

    tail.cols <- ceiling(((100-CI)/200) * n.reps) + 1

    preds.mean <- rep(0, n.pred.obs)
    preds.lower.tail <- matrix(NA, ncol = tail.cols, nrow = n.pred.obs)
    preds.upper.tail <- matrix(NA, ncol = tail.cols, nrow = n.pred.obs)

    if (return.pred.matrix) boot.preds <- matrix(NA, ncol = n.reps, nrow = n.pred.obs)

  }

  cat("gbm.pred.bootstrap - version 2.9","\n\n")
  cat("bootstrap resampling gbm.step model for ",response.name,"\n",sep="")
  cat("with ",n.trees," trees and ",n.obs," observations\n\n",sep="")
  if (bootstrap.predictions) {
    cat("prediction dataset has ",n.pred.obs," rows\n\n",sep="")
  } else {
    cat("no prediction dataset provided...\n\n")
  }

# initiate timing call

  z1 <- unclass(Sys.time())

# create gbm.fixed function call

  gbm.call.string <- paste("gbm.fixed(data=boot.data,gbm.x=gbm.x,gbm.y=gbm.y,",sep="")
  gbm.call.string <- paste(gbm.call.string,"family=family,learning.rate=lr,tree.complexity=tc,",sep="")
  gbm.call.string <- paste(gbm.call.string,"n.trees = ",n.trees,", site.weights = weights,verbose=FALSE)",sep="")

# now start the main bootstrap loop

  for (i in 1:n.reps) {

    if (i == 6 & verbose) {

      z2 <- unclass(Sys.time())
      est.time <- (z2 - z1)/60  # time for five reps
      est.time <- est.time * (n.reps/5) * 2  # multiply by two as sorting takes time
      cat("five bootstrap samples processed \n"," estimated time for completion is ", round(est.time,1)," minutes \n",sep="")
    } else {
      if (verbose) cat(i,"\n")
    }

# create a vector with which to select the bootstrap sample

    i.select <- sample(1:n.obs, n.obs, replace = T)
    boot.data <- train.data[i.select,]

    boot.model.gbm <- eval(parse(text=gbm.call.string))

    if (bootstrap.model) {

      i.match <- 1:n.obs %in% i.select # reference back to the original training data

      train.preds <- predict.gbm(boot.model.gbm, train.data,"link", n.trees = n.trees)
      train.fit.matrix[i.match,i] <- train.preds[i.match]
      train.pred.matrix[!i.match,i] <- train.preds[!i.match]
    }

    if (bootstrap.functions) {

      for (j in 1:n.preds) {
        #cycle through the first time and get the range of the functions

        k <- gbm.x[j]
        temp.frame <- function.pred.frame

        if (is.vector(train.data[,k])) {
          temp.vector <- seq(from=min(train.data[,k],na.rm=T),to=max(train.data[,k],na.rm=T),length=200)
          temp.frame[,j] <- temp.vector
        }
        if (is.factor(train.data[,k])) {
          temp.table <- table(train.data[,k])
          temp.factor <- rep(as.character(names(temp.table)),length.out = 200)
          temp.factor <- factor(temp.factor)
          temp.frame[,j] <- temp.factor
        }
        function.pred <- predict.gbm(boot.model.gbm, temp.frame,"link", n.trees = n.trees)
        function.matrix[i,j,] <- function.pred
      }
    }

    if (bootstrap.predictions) {

      preds <- predict.gbm(boot.model.gbm, pred.data,"link", n.trees = n.trees)
      preds.mean <- preds.mean + preds
      if (return.pred.matrix) boot.preds[,i] <- preds

      if (i <= tail.cols) {
        preds.lower.tail[,i] <- preds
        preds.upper.tail[,i] <- preds
      } else {
        preds.lower.tail <- t(apply(preds.lower.tail,1,sort))
        preds.upper.tail <- t(apply(preds.upper.tail,1,sort))
        preds.lower.tail[,tail.cols] <- preds
        preds.upper.tail[,1] <- preds
      }
    }
  }  # end of main bootstrap loop

# now caclulate the mean fitted and predicted values within the training data

  cat("calculating final values \n")

  if (bootstrap.model) {

    train.fit.mean <- apply(train.fit.matrix, 1, mean, na.rm=T)

# interpolate the percentiles for fits

    train.fit.lower <- rep(0,n.obs)
    train.fit.upper <- rep(0,n.obs)
    train.fit.count <- rep(0,n.obs)

    for (i in 1:n.obs) {
      temp <- train.fit.matrix[i,]
      temp <- temp[!is.na(temp)]
      train.fit.count[i] <- length(temp)
      if (length(temp) >= 2) {
        train.fit.upper[i] <- approx(ppoints(temp),sort(temp),0.975)$y
        train.fit.lower[i] <- approx(ppoints(temp),sort(temp),0.025)$y
      } else {
        train.fit.upper[i] <- NA
        train.fit.lower[i] <- NA
      }
    }

    train.pred.mean <- apply(train.pred.matrix,1,mean,na.rm=T)

# interpolate the percentiles and predictions on withheld training obs

    train.pred.lower <- rep(0,n.obs)
    train.pred.upper <- rep(0,n.obs)
    train.pred.count <- rep(0,n.obs)

    for (i in 1:n.obs) {
      temp <- train.pred.matrix[i,]
      temp <- temp[!is.na(temp)]
      train.pred.count[i] <- length(temp)
      if (length(temp) >= 2) {
        train.pred.upper[i] <- approx(ppoints(temp),sort(temp),0.975)$y
        train.pred.lower[i] <- approx(ppoints(temp),sort(temp),0.025)$y
      } else {
        train.pred.upper[i] <- NA
        train.pred.lower[i] <- NA
      }
    }
    if (family == "bernoulli") {
      train.fit.mean <- exp(train.fit.mean)/(1 + exp(train.fit.mean))
      train.fit.upper <- exp(train.fit.upper)/(1 + exp(train.fit.upper))
      train.fit.lower <- exp(train.fit.lower)/(1 + exp(train.fit.lower))

      train.pred.mean <- exp(train.pred.mean)/(1 + exp(train.pred.mean))
      train.pred.upper <- exp(train.pred.upper)/(1 + exp(train.pred.upper))
      train.pred.lower <- exp(train.pred.lower)/(1 + exp(train.pred.lower))
    }
    if (family == "poisson") {
      train.fit.mean <- exp(train.fit.mean)
      train.fit.upper <- exp(train.fit.upper)
      train.fit.lower <- exp(train.fit.lower)

      train.pred.mean <- exp(train.pred.mean)
      train.pred.upper <- exp(train.pred.upper)
      train.pred.lower <- exp(train.pred.lower)
    }
  }

# now calculate values for the fitted functions

  if (bootstrap.functions) {

    function.dataframe <- as.data.frame(matrix(0,nrow=200,ncol=n.preds*4))
    for (i in 1:n.preds) {
      j <- (i * 4) - 3
      names(function.dataframe)[j:(j+3)] <- paste(predictor.names[i],c(".vals",".lower",".mean",".upper"), sep="")
    }

    for (i in 1:n.preds) {
      k <- (i * 4) - 3
      if (is.vector(train.data[,gbm.x[i]])) {
        function.dataframe[,k] <- seq(min(train.data[,gbm.x[i]], na.rm = T), max(train.data[,gbm.x[i]], na.rm = T),length = 200)
        function.dataframe[,k+2] <- apply(function.matrix[,i,],2,mean)
        for (j in 1:200) {
          temp <- function.matrix[,i,j]
          function.dataframe[j,k+1] <- approx(ppoints(temp),sort(temp),0.025)$y
          function.dataframe[j,k+3] <- approx(ppoints(temp),sort(temp),0.975)$y
        }
      }
      if (is.factor(train.data[,gbm.x[i]])) {
        temp.table <- table(train.data[,gbm.x[i]])
        n <- length(temp.table)
        function.dataframe[,k] <- as.factor(rep(unlist(labels(temp.table)),length=200))
        for (j in 1:n) {
          temp <- function.matrix[,i,j]
          function.dataframe[j,k+1] <- approx(ppoints(temp),sort(temp),0.025)$y
          function.dataframe[j,k+3] <- approx(ppoints(temp),sort(temp),0.975)$y
        }
      }
    }
  }

# now a final sort of the predictions data

  if (bootstrap.predictions) {

    if (return.pred.matrix) boot.preds <- t(apply(boot.preds,1,sort))
    preds.mean <- preds.mean/n.reps
    preds.lower.tail <- t(apply(preds.lower.tail,1,sort))
    preds.upper.tail <- t(apply(preds.upper.tail,1,sort))

# and then clip off the irrelevant working columns

    preds.lower.tail <- preds.lower.tail[,1:tail.cols - 1]
    preds.upper.tail <- preds.upper.tail[,2:tail.cols]

    if (family == "bernoulli") {

      if (return.pred.matrix) boot.preds <- exp(boot.preds)/(1 + exp(boot.preds))
      preds.mean <- round(exp(preds.mean)/(1 + exp(preds.mean)),4)
      preds.lower.tail <- round(exp(preds.lower.tail)/(1 + exp(preds.lower.tail)),4)
      preds.upper.tail <- round(exp(preds.upper.tail)/(1 + exp(preds.upper.tail)),4)
    }

    if (family == "poisson") {
      if (return.pred.matrix) boot.preds <- exp(boot.preds)
      preds.mean <- round(exp(preds.mean),4)
      preds.lower.tail <- round(exp(preds.lower.tail),4)
      preds.upper.tail <- round(exp(preds.upper.tail),4)
    }

# and assign the upper and lower CI estimates

    preds.upper.limit <- preds.upper.tail[,1]
    preds.lower.limit <- preds.lower.tail[,tail.cols - 1]
  }


# final timing call

  z2 <- unclass(Sys.time())

  elapsed.time <- round((z2 - z1)/60,2)

  cat("analysis took ",round(elapsed.time,1)," minutes \n\n")

  gbm.call$n.bootstrap.reps <- n.reps
  gbm.call$bootstrap.CI <- CI
  if (bootstrap.predictions) gbm.call$prediction.dataframe <- function.pred.frame.name
  gbm.call$bootstrap.time <- elapsed.time

  final.object <- list(gbm.call = gbm.call)
  if (bootstrap.model) {
    final.object$train.fit.stats <- data.frame(fit.mean = train.fit.mean, fit.lower = train.fit.lower, fit.upper = train.fit.upper, fit.count = train.fit.count)
    final.object$train.pred.stats <- data.frame(pred.mean = train.pred.mean, pred.lower = train.pred.lower, pred.upper = train.pred.upper, pred.count = train.pred.count)
  }

  if (bootstrap.functions) final.object$function.dataframe <- function.dataframe

  if (bootstrap.predictions) {
    final.object$prediction.stats = data.frame(pred.mean = preds.mean, upper.limit = preds.upper.limit, lower.limit = preds.lower.limit)
    if (return.tails) {
      final.object$preds.lower.tail <- preds.lower.tail
      final.object$preds.upper.tail <- preds.upper.tail
    }
    if(return.pred.matrix) final.object$pred.matrix <- boot.preds
  }

  return(final.object)
}


## Run BRT analysis 
#-----------------------------------------------------------------------------

brt <- function(resp.var, pred.vars, data, family = c("bernoulli", "gaussian", "poisson"), tree.complexity = 2, n.boot.effects=0, plot.layout = c(2,2), predict=FALSE, newdata=data, n.boot.pred=0, quiet=FALSE, ...) {
    #
    # Fit, evaluate and predict BRT (for one species)
    #
    # resp.var          response variable = column of data holding abundance or presence
    # pred.vars         predictor variables = columns of data holding environmental data
    # data              data.frame of presence/abundance and environmental data
    # family            distribution family
    #                   (bernoulli for presence-absence, gaussian or poisson for abundances)
    # tree.complexity   number of interactions
    # n.boot.effects    number of bootstraps for the estimation of the effects of the environment
    # plot.layout       dimension of the matrix of effects plots
    # predict           whether to perform the prediction or only fit the model
    # newdata           environmental conditions at locations where presence/abundance should be predicted
    # n.boot.effects    number of bootstraps for the prediction
    #                   (allows to estimate error through cross validation)
    # quiet             when TRUE, do not print messages when true
    # ...               passed to dismo::gbm.step, plot.pred.brt
    #

    suppressPackageStartupMessages(require("gbm", quietly=T))
    suppressPackageStartupMessages(require("dismo", quietly=T))

    # arguments checks
    if (!(resp.var %in% names(data))) {
        stop("Response variables not found in the original dataset")
    }
    if (!all(pred.vars %in% names(data))) {
        stop("Not all response variables are in the original dataset")
    }
    if (predict & !all(pred.vars %in% names(newdata))) {
        stop("Not all response variables are in the new dataset")
    }

    if ((n.boot.effects > 0 & n.boot.effects < 100)) {
        warning("Less than 100 bootstraps will fail.\nn.boot.effects was increased to 100", call.=FALSE)
        n.boot.effects = 100
    }
    if (predict & (n.boot.pred > 0 & n.boot.pred < 100)) {
        warning("Less than 100 bootstraps will fail.\nn.boot.pred was increased to 100", call.=FALSE)
        n.boot.pred = 100
    }

    family <- tolower(family)
    family <- match.arg(family)

    # start the output object
    result <- list()
    class(result) <- c("brt", "list")

    ## Setup data
    #-----------------------------------------------------------------------------
    if ( ! quiet ) cat("   setup data\n")

    # remove observations where the current species information is Non-Available
    data <- data[!is.na(data[,resp.var]),]

    # make the data binomial if it needs to be
    if (family == "bernoulli") { data[,resp.var] = data[,resp.var] > 0 }

    # store reformated data in the result object
    result$data <- data
    result$resp.var <- resp.var
    result$pred.vars <- pred.vars


    ## Run BRT model
    #-----------------------------------------------------------------------------
    if ( ! quiet ) cat("   optimising BRT model ")

    # initial values before entering the while loop
    lr = 0.05
    no.trees = 0
    
    # NB: gbm.plot uses eval() in the global environment to find the dataset which was used to fit the model
    #     this is highly flawed and only works when used interactively directly from the command prompt
    #     to bypass this, we write the dataset to the global environment with a funky name
    myDataSetForGbmPlot <- data
    assign("myDataSetForGbmPlot", myDataSetForGbmPlot, .GlobalEnv)

    while ( no.trees < 1000 & lr > 0.0005 ) {
        if ( ! quiet ) { cat(".") }
        try( obj <- gbm.step(
                    data = myDataSetForGbmPlot,
                    gbm.x = match(pred.vars, names(data)),
                    gbm.y = match(resp.var, names(data)),
                    learning.rate = lr,
                    silent = TRUE,
                    plot.main = FALSE,
                    #
                    family = family,
                    tree.complexity = tree.complexity,
                    ...
                    )
        )

        # if the GBM does not converge, the return object is NULL or of size 0
        if ( ! is.null(obj) ) {
            if ( object.size(obj) > 0 ) {
                no.trees = obj$gbm.call$best.trees
            }
        } else {
            no.trees = 0
        }

        # decrease the learning rate
        lr = lr / 2
    }
    if ( ! quiet ) { cat("\n") }

    # store the gbm object in the result object
    result$obj = obj


    ## Store additional results
    #-----------------------------------------------------------------------------
    if ( ! quiet ) cat("   write results\n")

    # Deviance related information
    temp = list()

    # deviance explained
    temp$perc.deviance.explained = base::round( obj$cv.statistics$deviance.mean / obj$self.statistics$mean.null, 2)

    # Area Under the receiver operative Curve (AUC)
    # = quality of the prediction of presence/absence
    #   0.5 is indifferent
    temp$AUC = base::round(min(obj$cv.roc.matrix), 2)
    result$deviance = temp

    # Contributions
    temp = obj$contributions$rel.inf
    names(temp) = obj$contributions$var
    result$contributions = temp

    # remove the temp object, to be clean
    rm(temp)


    # Plot effects
    #-----------------------------------------------------------------------------
    if ( ! quiet ) cat("   plot effects\n")
    
    
    if (n.boot.effects == 0) {
    
        # no bootstrap, just plot
        gbm.plot(obj, plot.layout = plot.layout)
    
    } else {
    
        # perform bootstrap
        if ( ! quiet ) cat("    running bootstrap on BRT model\n")
        boot = NULL
        try( boot = gbm.bootstrap(obj, n.reps=n.boot.effects, verbose=FALSE) )
        if ( ! is.null(boot) ) {
            # when it runs correctly
            # store the output in the result object
            result$boot = boot
            # and use bootstrapped plots
            gbm.plot.boot(obj, boot, plot.layout = plot.layout)
    
        } else {
            # when it does not (usually because of too few bootstraps), issue a warning and plot the simple object
            warning("Bootstrapping of effects failed, not enough information to provide a confidence interval\nTry increasing n.boot.effects")
            gbm.plot(obj, plot.layout = plot.layout)
        }
    }
    
    # NB: remove the funky variable for gbm.plot
    rm(myDataSetForGbmPlot, pos=.GlobalEnv)


    ## Predict presence / abundance
    #-----------------------------------------------------------------------------

    if (predict) {
        if ( ! quiet ) cat("   make predictions\n")

        if (n.boot.pred > 1) {

            # bootstrap predictions
            temp = gbm.bootstrap(obj, pred.data = newdata, return.pred.matrix = T, n.reps=n.boot.pred, verbose=F)
            temp = temp$pred.matrix

            # prediction and confidence zone
            pred = apply(temp,1,mean)
            CVpred = apply(temp,1,cv)

        } else {

            # do a simple prediction
            pred = predict(obj, newdata = newdata, n.trees = obj$n.trees, type="response")
            CVpred = NA     # CV can't be computed without bootstrap
        }

        # store it in the result object
        result$prediction = data.frame(taxon=resp.var, newdata[,c("lat", "lon")], pred, CVpred)


        if ( ! quiet ) cat("   plot predictions\n")
        p <- plot.pred.brt(result, ...)
        # store the plot in the result object
        result$plot.pred <- p

        print(p)

    }

    # print info about the fit
    summary(result)

    return(result)
}

brts <- function(file, taxa, variables, lat.min=-80, lat.max=-30, lat.step=0.1, lon.min=-180, lon.max=180, lon.step=0.5, predict=FALSE, bin=FALSE, path="env_data", plot.layout=c(2,2), ...) {

  suppressPackageStartupMessages(require("stringr", quietly=TRUE))

  # read dataset
  if (file.exists(file)) {
    observed_data <- read.data(file)    
  } else {
    stop("Cannot find file : ", file)
  }
  # get selected taxa
  allTaxa <- setdiff(names(observed_data), c("lat", "lon"))
  taxa <- match.vars(taxa, allTaxa, quiet=FALSE)
  observed_data <- observed_data[, c("lat", "lon", taxa)]

  # filter-out unwanted locations
  observed_data <- observed_data[observed_data$lat >= lat.min & observed_data$lat <= lat.max &
                                 observed_data$lon >= lon.min & observed_data$lon <= lon.max,]
                                 
  # bin observation data
  if (bin) {
    observed_data <- rasterize(observed_data, c("lat", "lon"), precisions=c(lat.step, lon.step), fun=sum)
    # NB: when the initial data is presence, we compute total number of observations per bin
    #     but this is OK, because the family should then be `bernoulli` and the data is reconverted to presence/absence
    
  }
  
  # read selected variables from the database
  database <- read.env.data(variables, path=path, quiet=FALSE)
  
  # get full, expanded names of selected environment variables
  variables <- names(database)

  # get environment data for the observations
  obsdata <- associate.env.data(observed_data, database)

  # build prediction grid if needed
  if (predict) {
    prediction_grid <- build.grid(lat.min=lat.min, lat.max=lat.max, lat.step=lat.step,
                                  lon.min=lon.min, lon.max=lon.max, lon.step=lon.step)
    preddata <- associate.env.data(prediction_grid, database)
  } else {
    preddata <- obsdata
  }

  # prepare storage for the results
  result <- list()
  class(result) <- c("brt.list", "list")

  for (i in seq(along=taxa)) {
    message("-> Running BRT for ", taxa[i])

    # prepare names of the files where the results will be written
    # taxon name
    cTaxon <- taxa[i]
    # data file without extension
    fileName <- str_replace(file, "\\.(csv|xls|txt)$", "")
    # prepare a basic name from this
    baseName <- str_c(fileName, cTaxon, "BRT", sep="-")

    # prepare specific filenames
    pdfFile <- str_c(baseName, ".pdf", sep="")
    csvFile <- str_c(baseName, ".csv", sep="")
    rdataFile <- str_c(baseName, ".Rdata", sep="")

    # open the PDF
    pdf(pdfFile, paper="a4r")

    # b <- brt(resp.var=taxa[i], pred.vars=variables, data=obsdata, predict=predict, newdata=preddata, plot.layout=plot.layout, ...)

    b <- tryCatch(
      brt(resp.var=taxa[i], pred.vars=variables, data=obsdata, predict=predict, newdata=preddata, plot.layout=plot.layout, ...),
      # do not stop on error
      error=function(e) {
        warning(e)
        return(NULL)
      }
    )

    # close PDF
    dev.off()

    if (is.null(b)) {
      # there was an error, just skip to the next species
      next

    } else {
      # write the results in files
      save(b, file=rdataFile)
      if (predict) {
        write.table(b$prediction, file=csvFile, sep=",", append=(i!=1), col.names=(i==1), row.names=FALSE)
      }
    }

    # store the result in the total object
    result[[taxa[i]]] <- b
  }

  return(invisible(result))
}



## Analyze BRT output
#-----------------------------------------------------------------------------

print.brt <- function(x, ...) {
  #
  # Print information about a brt model
  #
  # x   object of class brt
  #

  resp.var <- x$obj$gbm.call$response.name
  pred.vars <- paste(x$obj$gbm.call$predictor.names, collapse=" + ")
  cat(resp.var, "~", pred.vars, "\n")
  cat("Family =",x$obj$gbm.call$family)
  if (!is.null(x$boot)) {
    cat(", bootstrapped effects")
  }
  if (!is.null(x$prediction)) {
    cat(", ")
    if (!all(is.na(x$prediction$CVpred))) {
      cat("bootstraped ")
    }
    cat("prediction")
  }
  cat("\n")

}


summary.brt <- function(x, ...) {
  #
  # Summarize information about the fit of a brt model
  #
  # x   object of class brt
  #

  cat("\n")

  resp.var <- x$obj$gbm.call$response.name
  cat(resp.var, ", Family = ",x$obj$gbm.call$family, "\n", sep="")

  dev <- x$deviance
  cat(dev$perc.deviance.explained * 100, "% of deviance explained, AUC =", dev$AUC, "\n")

  print(data.frame("Contribution"=x$contributions), check.names=F)

  if (!is.null(x$prediction)) {
    latRange <- range(x$prediction$lat)
    lonRange <- range(x$prediction$lon)
    latStep <- diff(sort(unique(x$prediction$lat)))[1]
    lonStep <- diff(sort(unique(x$prediction$lon)))[1]
    cat("\nPrediction : lat = ", latRange[1], " ..", latStep,".. ", latRange[2], "\n")
    cat("             lon = ", lonRange[1], " ..", lonStep,".. ", lonRange[2], "\n")
    print(summary(x$prediction$pred))
  }

  cat("\n")

  return(invisible(x))
}

summary.brt.list <- function(x) {
  #
  # Same as summary brt for the output of several models
  #
  # x   an object of class brt.list
  #
  
  lapply(x, summary.brt)
  return(invisible(x))
}


plot.brt <- function(x, plot.layout=c(2,2), ...) {
  #
  # Plot effects in a brt model
  #
  # x             object of class brt
  # plot.layout   dimensions of the matrix of plots
  #
  
  # NB: gbm.plot uses eval() in the global environment to find the dataset which was used to fit the model
  #     this is highly flawed and only works when used interactively directly from the command prompt
  #     to bypass this, we write the dataset to the global environment with a funky name
  myDataSetForGbmPlot <- x$data
  assign(x="myDataSetForGbmPlot", value=myDataSetForGbmPlot, pos=.GlobalEnv)

  if (is.null(x$boot)) {
    gbm.plot(x$obj, plot.layout = plot.layout)
  } else {
    gbm.plot.boot(x$obj, boot, plot.layout = plot.layout)
  }

  rm(myDataSetForGbmPlot, pos=.GlobalEnv)

  return(invisible(x))
}


plot.pred <- function(x, ...) {
  UseMethod("plot.pred")
}

plot.pred.brt <- function(x, type=c("quick", "full"), overlay.stations=FALSE, ...) {
  #
  # Plot BRT predictions
  #
  # x                     object of class brt
  # type                  quick = point based plot (with subsampling)
  #                       full = tile based plot with the whole dataset
  # overlay.observations  add points at the location of observation in to original data
  # ...                   passed on to polar.ggplot
  #

  suppressPackageStartupMessages(require("ggplot2"))

  if (is.null(x$prediction)) {
    # check that predictions are there
    stop("No predictions in this brt object; re-run the analysis with predict=TRUE")

  } else {
    # check arguments
    type <- match.arg(type)

    # check wether CV error is computed
    if (all(is.na(x$prediction$CVpred))) {
      # map only prediction
      mapping <- aes_string(fill="pred", colour="pred")
    } else {
      # map prediction as colour and error as transparency
      mapping <- aes(fill=pred, colour=pred, alpha=CVpred)
    }

    # main plot
    if (type == "quick") {
      # use points (and possibly regrid)
      p = polar.ggplot(x$prediction, mapping=mapping, geom="point", ...)
    } else {
      # use tiles
      p = polar.ggplot(x$prediction, mapping=mapping, geom="tile", ...)
    }

    # overlay stations
    if (overlay.stations) {
      if (!all(c("lat", "lon") %in% names(x$data))) {
        warning("Cannot overlay data points because the coordinates were not in the original dataset", immediate.=TRUE)
      } else {
        p <- p + geom_point(aes_string(x="lon", y="lat", size=x$resp.var), data=x$data) + scale_size(range=c(1,2), breaks=seq(from=0, to=max(x$data[,x$resp.var], by=1)))
      }
    }
    
    # add a title
    p = p + opts(title=paste(x$obj$gbm.call$response.name, "- BRT"))
        
    return(p)
  }
}

