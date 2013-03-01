#
#     Boosted Regression Trees
#
#     - bootstrapping functions (complement those in package dismo)
#     - wrap around BRT function which does model, effects, bootstraps, predictions
#     - fully automated function
#     - GUI interface
#
# (c) Copyright 2011-2012 S Mormede
#                         J-O Irisson
#     http://creativecommons.org/licenses/by/3.0/
#
#-----------------------------------------------------------------------------



## Functions modified / absent from the dismo package
#-----------------------------------------------------------------------------

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
  response.name <- gbm.call$response.name
  predictor.names <- gbm.call$predictor.names
  dd <- data.frame(matrix(gbm.object$data$x, nrow=nrow(gbm.object$data$x.order)))
  classes <- attr(gbm.object$Terms, "dataClasses")[predictor.names]
  factorVars <- which(classes == "factor")
  if (length(factorVars) > 0) {
    for (i in as.numeric(factorVars)) {
      levels <- as.numeric(gbm.object$var.levels[[i]]) - 1
      dd[,i] <- factor(dd[,i], levels=levels, labels=gbm.object$var.levels[[i]])
    }
  }
  train.data <- data.frame(gbm.object$data$y, dd)
  names(train.data) <- c(response.name, predictor.names)
  n.obs <- nrow(train.data)
  gbm.y <- 1
  gbm.x <- 2:ncol(train.data)
  family <- gbm.call$family
  lr <- gbm.call$learning.rate
  tc <- gbm.call$tree.complexity
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

  # cat("gbm.pred.bootstrap - version 2.9","\n\n")
  # cat("bootstrap resampling gbm.step model for ",response.name,"\n",sep="")
  # cat("with ",n.trees," trees and ",n.obs," observations\n\n",sep="")
  # if (bootstrap.predictions) {
  #   cat("prediction dataset has ",n.pred.obs," rows\n\n",sep="")
  # } else {
  #   cat("no prediction dataset provided...\n\n")
  # }

# initiate timing call

  z1 <- unclass(Sys.time())

# create gbm.fixed function call

  gbm.call.string <- paste("gbm.fixed(data=boot.data,gbm.x=gbm.x,gbm.y=gbm.y,",sep="")
  gbm.call.string <- paste(gbm.call.string,"family=family,learning.rate=lr,tree.complexity=tc,",sep="")
  gbm.call.string <- paste(gbm.call.string,"n.trees = ",n.trees,", site.weights = weights,verbose=FALSE)",sep="")

# now start the main bootstrap loop

  for (i in 1:n.reps) {

    # if (i == 6 & verbose) {
    #
    #   z2 <- unclass(Sys.time())
    #   est.time <- (z2 - z1)/60  # time for five reps
    #   est.time <- est.time * (n.reps/5) * 2  # multiply by two as sorting takes time
    #   cat("five bootstrap samples processed \n"," estimated time for completion is ", round(est.time,1)," minutes \n",sep="")
    # } else {
    #   if (verbose) cat(i,"\n")
    # }

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

  # cat("calculating final values \n")

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

  # cat("analysis took ",round(elapsed.time,1)," minutes \n\n")

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


gbm.fixed.etc <- function(data, gbm.x, gbm.y, tree.complexity = 1, site.weights = rep(1, nrow(data)), verbose = TRUE, learning.rate = 0.001, n.trees = 2000, bag.fraction = 0.5, family = "bernoulli", keep.data = FALSE, var.monotone = rep(0, length(gbm.x)), ...) {
  #
  # Fix gbm.fixed to
  #
  # allow the ellipsis argument (map all other arguments but ...)
  res <- gbm.fixed(
    data=data,
    gbm.x=gbm.x,
    gbm.y=gbm.y,
    tree.complexity=tree.complexity,
    site.weights=site.weights,
    verbose=verbose,
    learning.rate=learning.rate,
    n.trees=n.trees,
    bag.fraction=bag.fraction,
    family=family,
    keep.data=keep.data,
    var.monotone=var.monotone
  )

  # return the species name
  res$gbm.call$response.name <- names(data)[gbm.y]

  return(res)
}


## Run BRT analysis
#-----------------------------------------------------------------------------

compute.brt <- function(
  #
  # Fit, evaluate and predict BRT (for one species)
  #
  resp.var,               # response variable = column of data holding abundance or presence
  pred.vars,              # predictor variables = columns of data holding environmental data
  data,                   # data.frame of presence/abundance and environmental data
  family="bernoulli",     # distribution family ("bernoulli" for presence-absence, "gaussian" or "poisson" for abundances)
  tree.complexity=2,      # number of interactions
  n.boot.effects=0,       # number of bootstraps for the estimation of the effects of the environment
  predict=FALSE,          # whether to perform the prediction or only fit the model
  newdata=data,           # environmental conditions at locations where presence/abundance should be predicted
  extrapolate.env=FALSE,  # whether to extrapolate outside of the environmental range for the prediction. FALSE removes the points, TRUE keeps them, NA replaces the values with NA
  n.boot.pred=0,          # number of bootstraps for the prediction (allows to estimate error through cross validation)
  n.trees.fixed=0,        # if > 0, specifies the fixed number of trees to use in the BRT model. Otherwise, the number of trees is estimated by a stepwise procedure
  quiet=FALSE,            # do not print messages when TRUE
  site.weights=rep(1, nrow(data)), # location weights
  ...                     # passed to dismo::gbm.step or dismo::gbm.fixed depending on n.trees.fixed
)
{
    suppressPackageStartupMessages(require("gbm", quietly=T))
    suppressPackageStartupMessages(require("dismo", quietly=T))

    ## Check arguments
    #-------------------------------------------------------------------------
    if (!(resp.var %in% names(data))) {
        stop("Response variable not found in the original dataset")
    }
    if (!all(pred.vars %in% names(data))) {
        stop("Not all explanatory variables are in the original dataset")
    }
    if (predict & !all(pred.vars %in% names(newdata))) {
        stop("Not all explanatory variables are in the prediction dataset")
    }

    if (n.trees.fixed > 0 & (n.boot.effects > 0 | (predict & n.boot.pred > 0))) {
      warning("Bootstrapping cannot be done with a fixed number of trees (quick computation). Using the full computation, with optimization of the number of trees instead.")
      n.trees.fixed <- 0
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
    family <- match.arg(family, c("bernoulli", "gaussian", "poisson"))

    # start the output object
    result <- list()
    class(result) <- c("brt", "list")

    ## Setup data
    #-------------------------------------------------------------------------
    if ( ! quiet ) cat("   setup data\n")

    # remove observations where either
    # - the current species information
    # - one of the explanatory variables
    # is Non-Available
    hasNA <- aaply(is.na(data[,c("lat", "lon", resp.var, pred.vars)]), 1, any, .expand=FALSE)
    data <- data[!hasNA,]
    site.weights <- site.weights[!hasNA]

    # remove NAs in the prediction dataset
    newdata <- na.omit(newdata)

    # make the data binomial if it needs to be
    if (family == "bernoulli") { data[,resp.var] = data[,resp.var] > 0 }

    # store reformated data in the result object
    result$data <- data
    result$resp.var <- resp.var
    result$pred.vars <- pred.vars


    ## Run BRT model
    #-------------------------------------------------------------------------

    # NB: gbm.plot uses eval() in the global environment to find the dataset which was used to fit the model
    #     to avoid conflicts, we use a funky name
    myFunkyDatasetNameForGbmPlot <- data

    if (n.trees.fixed <= 0) {
      # initial learning rate value
      lr = 0.05

      # estimate number of trees required, using gbm.step
      if ( ! quiet ) cat("   optimise BRT model ")

      no.trees = 0
      while ( no.trees < 1000 & lr > 0.0005 ) {
        if ( ! quiet ) { cat(".") }
        obj <- tryCatch(
          gbm.step(
            data = myFunkyDatasetNameForGbmPlot,
            gbm.x = match(pred.vars, names(data)),
            gbm.y = match(resp.var, names(data)),
            tree.complexity = tree.complexity,
            learning.rate = lr,
            family = family,
            plot.main = FALSE,
            silent = TRUE,
            site.weights = site.weights,
            ...
          ),
          # transform errors into warnings and make obj=NULL when that happens
          error=function(e) {
            warning(e)
            NULL
          }
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

    } else {
      # fit model with fixed number of trees using gbm.fixed
      if ( ! quiet ) cat("   fit BRT model\n")

      obj <- tryCatch(
        gbm.fixed.etc(
          data = myFunkyDatasetNameForGbmPlot,
          gbm.x = match(pred.vars, names(data)),
          gbm.y = match(resp.var, names(data)),
          tree.complexity = tree.complexity,
          verbose = FALSE,
          learning.rate = 2/n.trees.fixed,
          n.trees = n.trees.fixed,
          family = family,
          site.weights = site.weights,
          ...
        ),
        # transform errors into warnings and make obj=NULL when that happens
        error=function(e) {
          warning(e)
          NULL
        }
      )
    }

    # store the gbm object in the result object
    result$obj = obj


    ## Store additional diagnostics
    #-------------------------------------------------------------------------
    if ( ! quiet ) cat("   compute diagnostics\n")

    # Deviance related information
    temp = list()

    if (n.trees.fixed <= 0) {
        # deviance explained
        temp$perc.deviance.explained = 1 - base::round( obj$cv.statistics$deviance.mean / obj$self.statistics$mean.null, 2)
        # Area Under the receiver operative Curve (AUC)
        # = quality of the prediction of presence/absence
        #   0.5 is indifferent
        temp$AUC = base::round(min(obj$cv.roc.matrix), 2)
    } else {
        # deviance explained
        temp$perc.deviance.explained = 1 - base::round(obj$self.statistics$resid.deviance / obj$self.statistics$null.deviance, 2)
        # AUC is not defined if cross-validation was not used to select the number of trees (i.e. if gbm.fixed was used)
        temp$AUC = NA
    }
    result$deviance = temp

    # Contributions
    temp = obj$contributions$rel.inf
    names(temp) = obj$contributions$var
    result$contributions = temp

    # remove the temp object, to be clean
    rm(temp)


    # Bootstrap effects
    #-------------------------------------------------------------------------
    if (n.boot.effects != 0) {

        # perform bootstrap
        if ( ! quiet ) cat("   bootstrap model and effects\n")
        boot <- NULL
        boot <- tryCatch(
          gbm.bootstrap(obj, n.reps=n.boot.effects, verbose=FALSE),
          # convert errors into warnings
          error=function(e) {
            warning(e)
            NULL
          }
        )
        if ( ! is.null(boot) ) {
            # when it runs correctly
            # store the output in the result object
            result$boot = boot
        } else {
            # when it does not (usually because of too few bootstraps), issue a warning
            warning("Bootstrapping of effects failed, not enough information to provide a confidence interval\nTry increasing n.boot.effects")
        }
    }


    ## Predict presence / abundance
    #-------------------------------------------------------------------------
    if (predict) {

        if (n.boot.pred > 1) {
            if ( ! quiet ) cat("   bootstrap predictions\n")

            # bootstrap predictions
            temp = gbm.bootstrap(obj, pred.data = newdata, return.pred.matrix = T, n.reps=n.boot.pred, verbose=F)
            temp = temp$pred.matrix

            # prediction and confidence zone
            pred = apply(temp,1,mean)
            CVpred = apply(temp,1,cv)

        } else {
            if ( ! quiet ) cat("   make predictions\n")

            # do a simple prediction
            pred = predict(obj, newdata = newdata, n.trees = obj$n.trees, type="response")
            CVpred = NA     # CV can't be computed without bootstrap
        }

        # store it as a data.frame
        prediction = data.frame(taxon=resp.var, newdata[,c("lat", "lon")], pred, CVpred)

        # if we choose *not* to extrapolate, remove the data outside of the original environmental range
        if (is.na(extrapolate.env)) {
          # extrapolate.env = NA means we want to replace extrapolated points with NA; i.e. we don't want to extrapolate
          extrapolate <- FALSE
        } else {
          extrapolate <- extrapolate.env
        }
        if ( ! extrapolate ) {
          # detect points in the prediction range which are outside the original range
          outsideRange <- c()
          for (pred.var in pred.vars) {
            cData <- data[,pred.var]        # input data
            cNewData <- newdata[,pred.var]  # new, prediction data
            if (is.factor(cData)) {
              # for a factor, compute the available levels and only keep those where there is more that one observation (with only 1, the computation of effects in the BRT fails anyway)
              cLevels <- table(cData)
              cLevels <- names(cLevels[cLevels > 1])
              outsideRange <- c(outsideRange, which(! cNewData %in% cLevels))
            } else {
              # for a numeric value, compute the numerical range
              cRange <- range(cData, na.rm=T)
              outsideRange <- c(outsideRange, which(cNewData < cRange[1] | cNewData > cRange[2]))
            }
          }
          outsideRange <- unique(outsideRange)
          # inform the user about it
          message("   Removed ", length(outsideRange), " points outside of the original data environmental range")

          if (is.na(extrapolate.env)) {
            # replaced extrapolated points with NA
            prediction[outsideRange, c("pred", "CVpred")] <- NA
          } else {
            # remove extrapolated points
            prediction <- prediction[-outsideRange,]
          }
        }

        # store it in the resulting object
        result$prediction <- prediction
    }

    return(result)
}

brt <- function(
  #
  # User-friendly interface to BRT modelling
  #
  file,                   # name of the file where the presence/abundance data is
  taxa,                   # names (or abbreviations) of the taxa of interest
  variables,              # names (or abbreviations) of the variables to use for the prediction
  lat.min=-80, lat.max=-30, lat.step=0.1,   # definition of the prediction grid
  lon.min=-180, lon.max=180, lon.step=0.5,
  bin=FALSE,              # whether to bin the observation data on the prediction grid
  quick=TRUE,             # when TRUE, a fixed number of trees is used in the fit and the prediction plot is subsampled, to increase speed
  overlay.stations=FALSE, # when TRUE, stations in the observed data are overlaid on top of the prediction plot
  path=getOption("atlasr.env.data"),  # path to the environmental database
  # compute.brt() arguments
  predict=FALSE,          # whether to perform the prediction or only fit the model
  quiet=FALSE,            # do not print messages when TRUE
  n.trees.fixed=ifelse(quick, 1000, 0), # if > 0, specifies the fixed number of trees to use in the BRT model. Otherwise, the number of trees is estimated by a stepwise procedure
  save=TRUE,              # save output to files
  transformations=NULL,   # named vector of transformations applied to each variable (has to match variables)
  ...                     # passed to compute.brt()
)
{
  suppressPackageStartupMessages(require("stringr", quietly=TRUE))
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))

  # read dataset
  file <- clean.path(file)
  if (file.exists(file)) {
    input_data <- read.data(file)
  } else {
    stop("Cannot find file : ", file)
  }
  # get selected taxa
  allTaxa <- setdiff(names(input_data), c("lat", "lon"))
  taxa <- match.vars(taxa, allTaxa, quiet=FALSE)
  input_data <- input_data[, c("lat", "lon", taxa)]

  # bin observation data
  if (bin) {
    # round lat and lon on the grid
    input_data$lon <- round_any(input_data$lon, 0.1)
    input_data$lat <- round_any(input_data$lat, 0.1)

    # compute number of data points per bin
    nb <- ddply(input_data, ~lat+lon, function(x) { c(weights=nrow(x)) })

    # compute weight
    nb$weights <- 1 / nb$weights

    # associate weights to input_data
    input_data <- join(input_data, nb, by=c("lon", "lat"))
  } else {
    input_data$weights <- 1
  }

  # read selected variables from the database
  database <- read.env.data(variables, path=path, quiet=FALSE)
  # remove information on land
  database <- mask.env.data(database, path=path)

  # get full, expanded names of selected environment variables
  variables <- names(database)

  # get environment data for the observations
  obsdata <- associate.env.data(input_data, database)

  # build prediction grid if needed
  if (predict) {
    prediction_grid <- build.grid(
                          lat.min=lat.min, lat.max=lat.max, lat.step=lat.step,
                          lon.min=lon.min, lon.max=lon.max, lon.step=lon.step
    )
    preddata <- associate.env.data(prediction_grid, database)
  } else {
    preddata <- obsdata
  }

  # transform data if needed
  if (!is.null(transformations)) {
    obsdata[,variables] <- transform.data(obsdata[,variables], transformations)
    preddata[,variables] <- transform.data(preddata[,variables], transformations)
  }

  # prepare storage for the results
  result <- list()
  class(result) <- c("brt.list", "list")

  for (i in seq(along=taxa)) {
    message("-> Run BRT for ", taxa[i])

    brtObj <- tryCatch(
      compute.brt(resp.var=taxa[i], pred.vars=variables, data=obsdata, predict=predict, newdata=preddata, n.trees.fixed=n.trees.fixed, site.weights=obsdata$weights, ...),
      # do not stop on error
      error=function(e) {
        warning(e)
        return(NULL)
      }
    )

    if (is.null(brtObj)) {
      # there was an error, just skip to the next species
      next

    } else {

      if (save) {
        # prepare names of the files where the results will be written
        # taxon name
        cTaxon <- taxa[i]
        # data file without extension
        fileName <- str_replace(file, "\\.(csv|xls|txt)$", "")
        # replace "." by "_" because shapefile writing does not support "."
        # fileName <- str_replace_all(fileName, fixed("."), "_")
        cTaxon <- str_replace_all(cTaxon, fixed("."), "_")
        # prepare directory
        dirName <- str_c(fileName, "/", cTaxon, "-BRT")
        dir.create(dirName, showWarnings=FALSE, recursive=TRUE)
        if (!file.exists(dirName)) {
          stop("Could not produce output directory : ", dirName)
        }
        # prepare a basic name from this
        baseName <- str_c(dirName, "/", cTaxon, "-BRT")

        # prepare specific filenames
        pdfFile <- str_c(baseName, ".pdf", sep="")
        csvFile <- str_c(baseName, ".csv", sep="")
        rdataFile <- str_c(baseName, ".Rdata", sep="")
        infoFile <- str_c(baseName, "-info.txt", sep="")

        # open the PDF
        pdf(pdfFile, width=11.7, height=8.3)
      }

      # plot effects
      if ( ! quiet ) cat("   plot results\n")
      plot.brt(brtObj, quick=quick, overlay.stations=overlay.stations, path=path)

      # print info about the fit
      summary(brtObj)

      if (save) {
        # close PDF
        dev.off()

        # write the results in files
        message("   Write output to ", dirName)
        save(brtObj, file=rdataFile)
        capture.output(summary(brtObj), file=infoFile)
        if (predict) {
          # CSV file
          write.table(brtObj$prediction, file=csvFile, sep=",", row.names=FALSE)

          # Shapefiles
          write.shapefile(brtObj$prediction, baseName, c("pred", "CVpred"))
        }
      }
    }

    # store the result in the total object
    result[[taxa[i]]] <- brtObj
  }

  message("\nDone\n")

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

  return(invisible(x))
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


plot.brt <- function(x, ...) {
  #
  # Produce all plots for a BRT object
  #
  # x   object of class brt
  #

  if (dev.interactive() | names(dev.cur()) == "null device") devAskNewPage(TRUE)

  plot.effects.brt(x, ...)

  if (!is.null(x$prediction)) {
    print(plot.pred.brt(x, ...))
  }

  devAskNewPage(FALSE)
}

plot.effects.brt <- function(x, n=200, ...) {

  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))
  suppressPackageStartupMessages(require("reshape2", quietly=TRUE))
  suppressPackageStartupMessages(require("gridExtra", quietly=TRUE))
  suppressPackageStartupMessages(require("gbm", quietly=TRUE))

  # extract into from object
  vars <- x$obj$var.names
  data <- x$data[,vars]
  contrib <- sort(x$contributions, decreasing=T)

  # detect if the object has bootstraps of effects
  bootstrap <- ! is.null(x$boot)

  # compute effects for all variables
  if ( ! bootstrap ) {
    # without bootstrap
    effects <- alply(vars, 1, function(variable, x, n, data) {

      # compute effects
      out <- plot.gbm(x$obj, i.var=variable, continuous.resolution = n, return.grid=TRUE)
      out$y <- out$y - mean(out$y, na.rm=T)
      out <- rename(out, c(y="fitted"))

      if (is.factor(out[,variable])) {
        # remove unobserved factor levels when data is discrete
        counts <- table(data[,variable])
        zeroCounts <- counts[which(counts == 0)]
        out$fitted[out[,variable] %in% names(zeroCounts)] <- NA
      } else {
        # compute quantiles of the data when it is continuous
        out$quantile <- quantile(data[,variable], probs=seq(0, 1, length.out=n))
      }
      return(out)
    }, x=x, n=n, data=data)

  } else {
    # with bootstraps
    d <- x$boot$function.dataframe

    effects <- alply(vars, 1, function(variable, d) {
      # extract columns corresponding to the current variable
      out <- d[str_detect(names(d),variable)]
      # rename them
      names(out) <- c(variable, "lower", "fitted", "upper")

      if (is.factor(out[,variable])) {
        # the data is padded with zeros, only keep the interesting bit
        out <- out[1:nlevels(out[,variable]),]
        # reorder levels as in the original data
        out[,variable] <- factor(out[,variable], levels=levels(data[,variable]))
        counts <- table(data[,variable])
        zeroCounts <- counts[which(counts == 0)]
        out[out[,variable] %in% names(zeroCounts), c("lower", "fitted", "upper")] <- NA
      } else {
        # add quantiles for continuous variables
        out$quantile <- quantile(data[,variable], probs=seq(0, 1, length.out=nrow(out)))
      }
      return(out)
    }, d=d)
  }

  # rank variables by contribution
  names(effects) <- vars
  effects <- effects[names(contrib)]

  # compute y-limits (same for all plot panels)
  if ( bootstrap ) {
    ranges <- ldply(effects, function(X) {
      c(min(X$lower, na.rm=T), max(X$upper, na.rm=T))
    })
  } else {
    ranges <- ldply(effects, function(X) {
      range(X$fitted, na.rm=T)
    })
  }
  ylim <- scale_y_continuous(limits=c(min(ranges$V1), max(ranges$V2)))

  # prepare labels
  varLabels <- str_c(names(contrib), " (", round(contrib,1), "%)")
  names(varLabels) <- names(contrib)

  # prepare plots
  plots <- llply(effects, function(X) {
    # base plot
    p <- ggplot(X, aes_string(x=names(X)[1]))

    # add y limits
    p <- p + ylim

    # add x label
    p <- p + xlab(varLabels[names(X)[1]])

    if (is.factor(X[,1])) {
      # discrete variable
      if (is.null(x$boot)) {
        # without bootstrap
        p <- p + geom_point(aes(y=fitted), na.rm=T, shape=15)
      } else {
        # with bootstrap
        p <- p + geom_linerange(aes(ymin=lower, ymax=upper), na.rm=T)
      }
    } else {
      # continuous variable
      if (! is.null(x$boot)) {
        # with bootstrap
        p <- p + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)
      }
      # with or without bootstrap
      p <- p +
        geom_path(aes(y=fitted), na.rm=T) +
        geom_rug(aes(x=quantile), alpha=0.6)
    }
    return(p)
  })

  # print plots on a grid
  do.call(grid.arrange, plots)

  return(invisible(plots))
}

plot.pred.brt <- function(x, quick=FALSE, overlay.stations=FALSE, geom="auto", ...) {
  #
  # Plot BRT predictions
  #
  # x                     object of class brt
  # quick                 subsample to 1º x 2º in lat x lon to plot more quickly
  # overlay.observations  add points at the location of observation in to original data
  # ...                   passed on to polar.ggplot
  #

  suppressPackageStartupMessages(require("ggplot2"))

  if (is.null(x$prediction)) {
    # check that predictions are there
    stop("No predictions in this brt object; re-run the analysis with predict=TRUE")

  } else {

    # check wether CV error is computed
    if (all(is.na(x$prediction$CVpred))) {
      # map only prediction
      mapping <- aes(fill=pred)
    } else {
      # map prediction as colour and error as transparency
      mapping <- aes(fill=pred, alpha=CVpred)
    }

    # main plot
    if (quick) {
      # subsample the plot
      geom = "raster"
      lat.precision <- 1
      lon.precision <- 2
    } else {
      # do not subsample and use the default geom
      lat.precision <- NULL
      lon.precision <- NULL
    }
    p <- polar.ggplot(x$prediction, mapping=mapping, geom=geom, lat.precision=lat.precision, lon.precision=lon.precision, ...)

    # overlay stations
    if (overlay.stations) {
      if (!all(c("lat", "lon") %in% names(x$data))) {
        warning("Cannot overlay data points because the coordinates were not in the original dataset", immediate.=TRUE)
      } else {
        # make the plot
        if (is.numeric(x$data[,x$resp.var])) {
          # numerical values (abundances)
          # = use coloured points with white outline
          p <- p + geom_point(aes_string(x="lon", y="lat", size=x$resp.var), data=x$data, alpha=0.5) + scale_size_continuous(range=0.5, 4)
        } else {
          # presence-absence values
          # = use points (presence) and crosses (absence)
          p <- p + geom_point(aes_string(x="lon", y="lat", shape=x$resp.var), data=x$data, alpha=0.7) + scale_shape_manual(values=c(4, 16))
        }
      }
    }

    # add a title
    p = p + ggtitle(paste(x$obj$gbm.call$response.name, "- BRT"))

    return(p)
  }
}

plot.brt.list <- function(x, ...) {
  lapply(x, plot.brt)
  return(invisible(x))
}

plot.effects.brt.list <- function(x, ...) {
  if (dev.interactive() | names(dev.cur()) == "null device") devAskNewPage(TRUE)
  lapply(x, plot.effects.brt)
  devAskNewPage(FALSE)
  return(invisible(x))
}

plot.pred.brt.list <- function(x, ...) {
  if (dev.interactive() | names(dev.cur()) == "null device") devAskNewPage(TRUE)
  print(lapply(x, plot.pred.brt))
  devAskNewPage(FALSE)
  return(invisible(x))
}


## GUI
#-----------------------------------------------------------------------------

do.brt <- function() {
  #
  # Open a GUI to select the arguments of the brt() function
  #

  suppressPackageStartupMessages(require("rpanel"))


  # default dimensions (in px)
  w <- 600      # width of the window
  w.h <- 700    # height of the window
  h <- 50       # height of elements
  spacer <- 10  # height of spacer

  # main window
  win <- rp.control(title="Run BRT model", size=c(w,w.h))

  # NB: positions are x, y, width, eight
  #     x, y are the coordinates of the top-left hand corner

  rp.button(win, title="Choose file", pos=c(0, 0, w/4, h), action=function(win) {

    # choose the data file
    file <- file.choose()
    # file <- "../data/Austropallene.csv"
    win$file <- clean.path(file)


    # write the filename, for information
    rp.text(win, txt=file, initval=file, pos=c(w/4, 0, 3*w/4, h))

    # build species list
    suppressMessages(data <- read.data(file))

    # the selection panel should be 380px high for the whole window to fit in 700px height
    hSel <- 380
    # this is roughly equivalent to 20 rows in the listbox
    nRows <- 20

    # we can fit ~20 checkboxes into that space, so we need to split the checkbox list into columns when there are more than that
    nBoxes <- 20
    # however, it is not possible to do that programmatically in a for loop because variables with numbered names have to be defined in each case and this requires eval(parse(...)) constructs which mostly don't work
    # so we only accommodate for a situation with two columns
    allTaxa <- setdiff(names(data), c("lat", "lon"))
    nTaxa <- length(allTaxa)

    if (nTaxa <= nBoxes) {
      # when all taxa can fit into one column, just use this
      rp.checkbox(win, var=taxa, labels=allTaxa, title="Taxa", initval=c(TRUE, rep(FALSE, times=nTaxa-1)), pos=c(0, h, 2*w/3, hSel), action=function(win) {
        if (all(! win$taxa)) {
          rp.messagebox("At least one taxon must be selected", title="Warning")
        }
        return(win)
      })

    } else {
      # otherwise, split in half
      half <- ceiling(nTaxa/2)
      rest <- nTaxa - half
      allTaxa1 <- allTaxa[1:half]
      allTaxa2 <- allTaxa[(half+1):nTaxa]

      rp.checkbox(win, var=taxa1, labels=allTaxa1, title="Taxa", initval=c(TRUE, rep(FALSE, times=half-1)), pos=c(0, h, w/3, hSel), action=function(win) {
        if (all(! c(win$taxa1, win$taxa2))) {
          rp.messagebox("At least one taxon must be selected", title="Warning")
        }
        return(win)
      })
      rp.checkbox(win, var=taxa2, labels=allTaxa2, title="", initval=rep(FALSE, times=rest), pos=c(w/3, h+7, w/3, hSel-7), action=function(win) {
        if (all(! c(win$taxa1, win$taxa2))) {
          rp.messagebox("At least one taxon must be selected", title="Warning")
        }
        return(win)
      })

    }

    # build environment variables list
    allVariables <- list.env.data()
    rp.listbox.mult(win, var=variables, vals=allVariables, title="Variables",  rows=nRows, cols=22, pos=c(2*w/3, h, w/3, hSel-h), action=function(win) {
      if (all(win$variables == "")) {
        rp.messagebox("At least one variable must be selected", title="Warning")
      }
      return(win)
    })

    # variables transformations
    optionsFile <- tempfile()
    win <- rp.button(win, "Variables transformation", pos=c(2*w/3, hSel, w/3, h) , action=function(win) {
      if (all(win$variables == "")) {
        rp.messagebox("At least one variable must be selected", title="Warning")
      } else {
        do.brt.variables(variables=win$variables, file=optionsFile)
      }
      return(win)
    })


    # compute vertical coordinate of the middle of the window
    mid <- hSel + h + spacer


    # effects options
    checkH <- 4*h/5
    rp.radiogroup(win, family, values=c("bernoulli", "gaussian", "poisson"), title="Distribution", pos=c(0, mid, w/4, checkH*2))

    # options checkboxes
    rp.checkbox(win, bin, title="Bin input data\non 0.1 x 0.1 grid", initval=FALSE, pos=c(0, mid+checkH*2, w/4, checkH))
    rp.checkbox(win, extrapolate.env, title="Extrapolate envi-\nronmental range", initval=FALSE, pos=c(0, mid+checkH*3, w/4, checkH))
    rp.checkbox(win, quick, title="Quick computation\n(faster fit and plot)", initval=TRUE, pos=c(0, mid+checkH*4, w/4, checkH))

    rp.checkbox(win, bootstrap.effects, title="Bootstrap effects", initval=FALSE, pos=c(w/4, mid, w/4, checkH), action=function(win) {
      if (win$quick & win$bootstrap.effects) {
        rp.messagebox("Bootstrapping won't work with the quick computation.\nPlease disable one")
      }
      return(win)
    })
    rp.checkbox(win, predict, title="Predict", initval=FALSE, pos=c(w/4, mid+checkH, w/4, checkH))
    rp.checkbox(win, bootstrap.prediction, title="Bootstrap\nprediction", initval=FALSE, pos=c(w/4, mid+checkH*2, w/4, checkH), action=function(win) {
      if (win$quick & win$bootstrap.prediction) {
        rp.messagebox("Bootstrapping won't work with the quick computation.\nPlease disable one.")
      }
      return(win)
    })
    rp.checkbox(win, overlay.stations, title="Overlay stations\non prediction plot", initval=FALSE, pos=c(w/4, mid+checkH*3, w/4, checkH))
    rp.checkbox(win, save, title="Save output", initval=TRUE, pos=c(w/4, mid+checkH*4, w/4, checkH))

    # location
    rp.slider(win, lat.max,  from=-90, to=-30,  resolution=2,   title="North"   , initval=-30 , showvalue=TRUE, pos=c(w/2+w/8, mid    , w/4, h))
    rp.slider(win, lon.min,  from=-180, to=180, resolution=5,   title="West"    , initval=-180, showvalue=TRUE, pos=c(w/2, mid+h*1, w/4, h))
    rp.slider(win, lon.max,  from=-180, to=180, resolution=5,   title="East"    , initval=180 , showvalue=TRUE, pos=c(3*w/4, mid+h*1, w/4, h))
    rp.slider(win, lat.min,  from=-90, to=-30,  resolution=2,   title="South"   , initval=-80 , showvalue=TRUE, pos=c(w/2+w/8, mid+h*2, w/4, h))
    rp.slider(win, lon.step, from=0.1, to=4  ,  resolution=0.1, title="Step lon", initval=1 , showvalue=TRUE, pos=c(w/2, mid+h*3, w/4, h))
    rp.slider(win, lat.step, from=0.1, to=4   , resolution=0.1, title="Step lat", initval=1 , showvalue=TRUE, pos=c(3*w/4, mid+h*3, w/4, h))


    # action buttons
    rowY <- mid+h*4+spacer
    rp.button(win, "Help", pos=c(0, rowY, w/4, h), action=function(win) {
      browseURL("https://github.com/jiho/atlasr/blob/master/documentation/HOWTO%20create%20a%20BRT%20model.md")
      return(win)
    })
    rp.button(win, "Cancel", quitbutton=TRUE, pos=c(w/2, rowY, w/4, h) , action=function(win) {
      message("Aborting");
      return(win)
    })
    rp.button(win, "Run", pos=c(3*w/4, rowY, w/4, h), action=function(win) {
      # print(win$file)

      if (nTaxa <= nBoxes) {
        taxa <- names(win$taxa)[win$taxa]
      } else {
        taxa <- names(win$taxa1)[win$taxa1]
        taxa <- c(taxa, names(win$taxa2)[win$taxa2])
      }
      # print(taxa)

      variables <- win$variables
      # print(variables)

      # variables transformations
      if (file.exists(optionsFile)) {
        # read the options in the file
        load(optionsFile)
        transformations <- opts$transformations
      } else {
        transformations <- NULL
      }
      # print(transformations)

      # print(win$family)
      # print(win$bootstrap.effects)
      if(win$bootstrap.effects) {
        n.boot.effects <- 200
      } else {
        n.boot.effects <- 0
      }

      # print(win$predict)
      # print(win$bootstrap.prediction)
      if(win$bootstrap.prediction) {
        n.boot.pred <- 200
      } else {
        n.boot.pred <- 0
      }

      # print(win$lat.max)
      # print(win$lat.min)
      # print(win$lat.step)
      # print(win$lon.max)
      # print(win$lon.min)
      # print(win$lon.step)

      # print(win$quick)
      # print(win$extrapolate.env)
      # print(win$overlay.stations)
      # print(win$bin)
      # print(win$save)


      # build the function call
      message("Command:")
      call <- str_c("brt(",
        "file=", str_c(deparse(win$file, width=500), collapse=""),
        ", taxa=", str_c(deparse(taxa, width=500), collapse=""),
        ", variables=", str_c(deparse(variables, width=500), collapse=""),
        ifelse( is.null(transformations),
          ", transformations=NULL",
          str_c(", transformations=", str_c("c(",str_c(names(transformations), str_c("\"", unlist(transformations),"\"") , sep="=", collapse=", "), ")"))
        ),
        ", lat.min=", win$lat.min, ", lat.max=", win$lat.max, ", lat.step=", win$lat.step,
        ", lon.min=", win$lon.min, ", lon.max=", win$lon.max, ", lon.step=", win$lon.step,
        ", predict=", win$predict,
        ", bin=", win$bin,
        ", family=", deparse(win$family),
        ", n.boot.effects=", n.boot.effects,
        ", n.boot.pred=", n.boot.pred,
        ", quick=", win$quick,
        ", extrapolate.env=", win$extrapolate.env,
        ", overlay.station=", win$overlay.stations,
        ", save=", win$save,
        ")"
      )
      cat(call, "\n")

      # execute call
      b <- eval(parse(text=call))

      return(win)

    })


    return(win)
  })

}

do.brt.variables <- function(variables, file) {
  #
  # Second part of the GUI. Called from do.brt()
  # Select variable transformations
  #

  suppressPackageStartupMessages(require("rpanel"))

  # defaults
  nVars <- length(variables)
  transformations <- rep("x", nVars)
  names(transformations) <- variables

  # if the options file exists, read it to use the previously saved transformations
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
  }

  # default dimensions (in px)
  w <- 500          # width of the window
  h <- 50           # height of elements
  spacer <- 10      # height of spacer
  varH <- nVars * 25 + spacer
  exampleH <- 3 * 25 + spacer
  windowH <- varH + exampleH + h

  # main window
  win <- rp.control(title="Variables", size=c(w, windowH), aschar=F)

  # transformations
  # TODO look into why with rp.textentry.immediate, the results are not all carried to the stage of the close button
  rp.textentry(win, var=transformsBox, labels=variables, title="Transformations", initval=transformations, pos=c(0,0,w,varH), action=function(win) {
    return(win)
  })

  # provide example transformations
  example.transforms.labels=c("Remove small values", "Remove large values", "Remove  values based on quantiles")
  example.transforms.functions=list('"x[x<0.01]=NA; x"', '"x[x>300]=NA; x"', '"q=quantile(x,c(0.1,0.9)); x[x<q[1]|x>q[2]]=NA; x"')
  rp.textentry(win, var=exampletransformBox, labels=example.transforms.labels, title="Example transformations", initval=example.transforms.functions, pos=c(0,varH,w,exampleH), action=function(win) {
    return(win)
  })

  # close button
  rp.button(win, title="Close", quitbutton=TRUE, pos=c(3/4*w, windowH-h, w/4, h), action=function(win) {
    # print("Transformations")
    # print(win$transformsBox)

    # extract transformations and weights and store them in named lists
    transformations <- as.list(win$transformsBox)

    # store them in a single object
    opts <- list(transformations=transformations)

    # save this object to the options file
    save(opts, file=file)

    return(win)
  })

  return(invisible(win))
}


