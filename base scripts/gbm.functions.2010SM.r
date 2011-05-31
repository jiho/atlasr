# Modified by Mormede from Leathwick and Elith as stated below
# to allow for transformation of outputs
# currently only GBM plot has been modified. 10 Feb 2010


gbm.step <- function (
  data,                                     # the input dataframe
  gbm.x,                                    # the predictors
  gbm.y,                                    # and response
  offset = NULL,                            # allows an offset to be specified
  fold.vector = NULL,                       # allows a fold vector to be read in for CV with offsets,
  tree.complexity = 1,                      # sets the complexity of individual trees
  learning.rate = 0.01,                     # sets the weight applied to inidivudal trees
  bag.fraction = 0.75,                      # sets the proportion of observations used in selecting variables
  site.weights = rep(1, nrow(data)),        # allows varying weighting for sites
  var.monotone = rep(0, length(gbm.x)),     # restricts responses to individual predictors to monotone
  n.folds = 10,                             # number of folds
  prev.stratify = TRUE,                     # prevalence stratify the folds - only for p/a data
  family = "bernoulli",                     # family - bernoulli (=binomial), poisson, laplace or gaussian
  n.trees = 50,                             # number of initial trees to fit
  step.size = n.trees,                      # numbers of trees to add at each cycle
  max.trees = 10000,                        # max number of trees to fit before stopping
  tolerance.method = "auto",                # method to use in deciding to stop - "fixed" or "auto"
  tolerance = 0.001,                        # tolerance value to use - if method == fixed is absolute,
                                            # if auto is multiplier * total mean deviance
  keep.data = FALSE,                        # keep raw data in final model
  plot.main = TRUE,                         # plot hold-out deviance curve
  plot.folds = FALSE,                       # plot the individual folds as well
  verbose = TRUE,                           # control amount of screen reporting
  silent = FALSE,                           # to allow running with no output for simplifying model)
  keep.fold.models = FALSE,                 # keep the fold models from cross valiation
  keep.fold.vector = FALSE,                 # allows the vector defining fold membership to be kept
  keep.fold.fit = FALSE,                    # allows the predicted values for observations from CV to be kept
  ...                                       # allows for any additional plotting parameters
)
{
#
# j. leathwick/j. elith - 19th September 2005
#
# version 2.9
#
# function to assess optimal no of boosting trees using k-fold cross validation
#
# implements the cross-validation procedure described on page 215 of
# Hastie T, Tibshirani R, Friedman JH (2001) The Elements of Statistical Learning:
# Data Mining, Inference, and Prediction Springer-Verlag, New York.
#
# divides the data into 10 subsets, with stratification by prevalence if required for pa data
# then fits a gbm model of increasing complexity along the sequence from n.trees to n.trees + (n.steps * step.size)
# calculating the residual deviance at each step along the way
# after each fold processed, calculates the average holdout residual deviance and its standard error
# then identifies the optimal number of trees as that at which the holdout deviance is minimised
# and fits a model with this number of trees, returning it as a gbm model along with additional information
# from the cv selection process
#
# updated 13/6/05 to accommodate weighting of sites
#
# updated 19/8/05 to increment all folds simultaneously, allowing the stopping rule
# for the maxinum number of trees to be fitted to be imposed by the data,
# rather than being fixed in advance
#
# updated 29/8/05 to return cv test statistics, and deviance as mean
# time for analysis also returned via unclass(Sys.time())
#
# updated 5/9/05 to use external function calc.deviance
# and to return cv test stats via predictions formed from fold models
# with n.trees = target.trees
#
# updated 15/5/06 to calculate variance of fitted and predicted values across folds
# these can be expected to approximate the variance of fitted values
# as would be estimated for example by bootstrapping
# as these will underestimate the true variance
# they are corrected by multiplying by (n-1)2/n
# where n is the number of folds
#
# updated 25/3/07 tp allow varying of bag fraction
#
# requires gbm library from Cran
# requires roc and calibration scripts of J Elith
# requires calc.deviance script of J Elith/J Leathwick
#
#

  require(gbm)

  if (silent) verbose <- FALSE

  # initiate timing call

  z1 <- unclass(Sys.time())

  # setup input data and assign to position one

  dataframe.name <- deparse(substitute(data))   # get the dataframe name

  data <- eval(data)
  x.data <- eval(data[, gbm.x])                 #form the temporary datasets
  names(x.data) <- names(data)[gbm.x]
  y.data <- eval(data[, gbm.y])
  sp.name <- names(data)[gbm.y]
  if (family == "bernoulli") prevalence <- mean(y.data)

  assign("x.data", x.data, env = globalenv())               #and assign them for later use
  assign("y.data", y.data, env = globalenv())

  offset.name <- deparse(substitute(offset))   # get the dataframe name
  offset = eval(offset)

  n.cases <- nrow(data)
  n.preds <- length(gbm.x)

  if (!silent) {
    cat("\n","\n","GBM STEP - version 2.9","\n","\n")
    cat("Performing cross-validation optimisation of a boosted regression tree model \n")
    cat("for",sp.name,"with dataframe",dataframe.name,"and using a family of",family,"\n\n")
    cat("Using",n.cases,"observations and",n.preds,"predictors \n\n")
  }

  # set up the selector variable either with or without prevalence stratification

  if (is.null(fold.vector)) {

    if (prev.stratify & family == "bernoulli") {
      presence.mask <- data[,gbm.y] == 1
      absence.mask <- data[,gbm.y] == 0
      n.pres <- sum(presence.mask)
      n.abs <- sum(absence.mask)

      # create a vector of randomised numbers and feed into presences
      selector <- rep(0,n.cases)
      temp <- rep(seq(1, n.folds, by = 1), length = n.pres)
      temp <- temp[order(runif(n.pres, 1, 100))]
      selector[presence.mask] <- temp

      # and then do the same for absences
      temp <- rep(seq(1, n.folds, by = 1), length = n.abs)
      temp <- temp[order(runif(n.abs, 1, 100))]
      selector[absence.mask] <- temp
    } else {
      #otherwise make them random with respect to presence/absence
      selector <- rep(seq(1, n.folds, by = 1), length = n.cases)
      selector <- selector[order(runif(n.cases, 1, 100))]
    }
  } else {
    if (length(fold.vector) != n.cases) stop("supplied fold vector is of wrong length")
    cat("loading user-supplied fold vector \n\n")
    selector <- eval(fold.vector)
  }

  # set up the storage space for results

  pred.values <- rep(0, n.cases)

  cv.loss.matrix <- matrix(0, nrow = n.folds, ncol = 1)
  training.loss.matrix <- matrix(0, nrow = n.folds, ncol = 1)
  trees.fitted <- n.trees

  model.list <- list(paste("model",c(1:n.folds),sep=""))     # dummy list for the tree models

  # set up the initial call to gbm

  if (is.null(offset)) {
    gbm.call <- paste("gbm(y.subset ~ .,data=x.subset, n.trees = n.trees,
    interaction.depth = tree.complexity, shrinkage = learning.rate,
    bag.fraction = bag.fraction, weights = weight.subset,
    distribution = as.character(family), var.monotone = var.monotone,
    verbose = FALSE)", sep="")
    }
  else {
    gbm.call <- paste("gbm(y.subset ~ . + offset(offset.subset),
    data=x.subset, n.trees = n.trees,
    interaction.depth = tree.complexity, shrinkage = learning.rate,
    bag.fraction = bag.fraction, weights = weight.subset,
    distribution = as.character(family), var.monotone = var.monotone,
    verbose = FALSE)", sep="")
    }

  n.fitted <- n.trees

# calculate the total deviance

  y_i <- y.data

  u_i <- sum(y.data * site.weights) / sum(site.weights)
  u_i <- rep(u_i,length(y_i))

  total.deviance <- calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = FALSE)

  mean.total.deviance <- total.deviance/n.cases

  tolerance.test <- tolerance

  if (tolerance.method == "auto") {
     tolerance.test <- mean.total.deviance * tolerance
  }

# now step through the folds setting up the initial call

  if (!silent){

    cat("creating",n.folds,"initial models of",n.trees,"trees","\n")

    if (prev.stratify & family == "bernoulli") cat("\n","folds are stratified by prevalence","\n","\n")
      else cat("\n","folds are unstratified","\n","\n")

    cat ("total mean deviance = ",round(mean.total.deviance,4),"\n","\n")

    cat("tolerance is fixed at ",round(tolerance.test,4),"\n","\n")

    if (tolerance.method != "fixed" & tolerance.method != "auto") {
      cat("invalid argument for tolerance method - should be auto or fixed","\n")
      return()}
  }

  if (verbose) cat("ntrees resid. dev.","\n")

  for (i in 1:n.folds) {

    model.mask <- selector != i  #used to fit model on majority of data
    pred.mask <- selector == i   #used to identify the with-held subset

    y.subset <- y.data[model.mask]
    x.subset <- x.data[model.mask,]
    weight.subset <- site.weights[model.mask]

    if (!is.null(offset)) {
      offset.subset <- offset[model.mask]
    }
    else {
      offset.subset <- NULL
    }

    model.list[[i]] <- eval(parse(text = gbm.call))

    fitted.values <- model.list[[i]]$fit  #predict.gbm(model.list[[i]], x.subset, type = "response", n.trees = n.trees)
    if (!is.null(offset)) fitted.values <- fitted.values + offset[model.mask]
    if (family == "bernoulli") fitted.values <- exp(fitted.values)/(1 + exp(fitted.values))
    if (family == "poisson") fitted.values <- exp(fitted.values)

    pred.values[pred.mask] <- predict.gbm(model.list[[i]], x.data[pred.mask, ], n.trees = n.trees)
    if (!is.null(offset)) pred.values[pred.mask] <- pred.values[pred.mask] + offset[pred.mask]
    if (family == "bernoulli") pred.values[pred.mask] <- exp(pred.values[pred.mask])/(1 + exp(pred.values[pred.mask]))
    if (family == "poisson") pred.values[pred.mask] <- exp(pred.values[pred.mask])

# calc training deviance

    y_i <- y.subset
    u_i <- fitted.values
    weight.fitted <- site.weights[model.mask]
    training.loss.matrix[i,1] <- calc.deviance(y_i, u_i, weight.fitted, family = family)

# calc holdout deviance

    y_i <- y.data[pred.mask]
    u_i <- pred.values[pred.mask]
    weight.preds <- site.weights[pred.mask]
    cv.loss.matrix[i,1] <- calc.deviance(y_i, u_i, weight.preds, family = family)

  } # end of first loop

# now process until the change in mean deviance is =< tolerance or max.trees is exceeded

  delta.deviance <- 1

  cv.loss.values <- apply(cv.loss.matrix,2,mean)
  if (verbose) cat(n.fitted,"  ",round(cv.loss.values,4),"\n","\n")

  if (!silent) cat("")
  if (!silent) cat("now adding trees...","\n")

  j <- 1

  while (delta.deviance > tolerance.test & n.fitted < max.trees) {  # beginning of inner loop

# add a new column to the results matrice..

    training.loss.matrix <- cbind(training.loss.matrix,rep(0,n.folds))
    cv.loss.matrix <- cbind(cv.loss.matrix,rep(0,n.folds))

    n.fitted <- n.fitted + step.size
    trees.fitted <- c(trees.fitted,n.fitted)

    j <- j + 1

    for (i in 1:n.folds) {

      model.mask <- selector != i  #used to fit model on majority of data
      pred.mask <- selector == i   #used to identify the with-held subset

      y.subset <- y.data[model.mask]
      x.subset <- x.data[model.mask,]
      weight.subset <- site.weights[model.mask]
      if (!is.null(offset)) {
        offset.subset <- offset[model.mask]
      }

      model.list[[i]] <- gbm.more(model.list[[i]], weights = weight.subset, step.size)

      fitted.values <- model.list[[i]]$fit # predict.gbm(model.list[[i]],x.subset, type = "response", n.trees = n.fitted)
      if (!is.null(offset)) fitted.values <- fitted.values + offset[model.mask]
      if (family == "bernoulli") fitted.values <- exp(fitted.values)/(1 + exp(fitted.values))
      if (family == "poisson") fitted.values <- exp(fitted.values)

      pred.values[pred.mask] <- predict.gbm(model.list[[i]], x.data[pred.mask, ], n.trees = n.fitted)
      if (!is.null(offset)) pred.values[pred.mask] <- pred.values[pred.mask] + offset[pred.mask]
      if (family == "bernoulli") pred.values[pred.mask] <- exp(pred.values[pred.mask])/(1 + exp(pred.values[pred.mask]))
      if (family == "poisson") pred.values[pred.mask] <- exp(pred.values[pred.mask])

# calculate training deviance

      y_i <- y.subset
      u_i <- fitted.values
      weight.fitted <- site.weights[model.mask]
      training.loss.matrix[i,j] <- calc.deviance(y_i, u_i, weight.fitted, family = family)

# calc holdout deviance

      u_i <- pred.values[pred.mask]
      y_i <- y.data[pred.mask]
      weight.preds <- site.weights[pred.mask]
      cv.loss.matrix[i,j] <- calc.deviance(y_i, u_i, weight.preds, family = family)

    }  # end of inner loop

    cv.loss.values <- apply(cv.loss.matrix,2,mean)

    if (j < 5) {
      if (cv.loss.values[j] > cv.loss.values[j-1]) {
        if (!silent) cat("restart model with a smaller learning rate or smaller step size...")
        return()
      }
    }

      if (j >= 20) {   #calculate stopping rule value
        test1 <- mean(cv.loss.values[(j-9):j])
        test2 <- mean(cv.loss.values[(j-19):(j-9)])
        delta.deviance <- test2 - test1
      }

      if (verbose) cat(n.fitted," ",round(cv.loss.values[j],4),"\n")

  } # end of while loop

# now begin process of calculating optimal number of trees

  training.loss.values <- apply(training.loss.matrix,2,mean)

  cv.loss.ses <- rep(0,length(cv.loss.values))
  cv.loss.ses <- sqrt(apply(cv.loss.matrix,2,var)) / sqrt(n.folds)

# find the target holdout deviance

  y.bar <- min(cv.loss.values)

# plot out the resulting curve of holdout deviance

  if (plot.main) {

    y.min <- min(cv.loss.values - cv.loss.ses)  #je added multiplier 10/8/05
    y.max <- max(cv.loss.values + cv.loss.ses)  #je added multiplier 10/8/05 }

    if (plot.folds) {
      y.min <- min(cv.loss.matrix)
      y.max <- max(cv.loss.matrix) }

      plot(trees.fitted, cv.loss.values, type = 'l', ylab = "holdout deviance",
          xlab = "no. of trees", ylim = c(y.min,y.max), ...)
      abline(h = y.bar, col = 2)

      lines(trees.fitted, cv.loss.values + cv.loss.ses, lty=2)
      lines(trees.fitted, cv.loss.values - cv.loss.ses, lty=2)

      if (plot.folds) {
        for (i in 1:n.folds) {
          lines(trees.fitted, cv.loss.matrix[i,],lty = 3)
      }
    }
  }

# identify the optimal number of trees

  target.trees <- trees.fitted[match(TRUE,cv.loss.values == y.bar)]

  if(plot.main) {
    abline(v = target.trees, col=3)
    title(paste(sp.name,", d - ",tree.complexity,", lr - ",learning.rate, sep=""))
  }

# estimate the cv deviance and test statistics
# includes estimates of the standard error of the fitted values added 2nd may 2005

  cv.deviance.stats <- rep(0, n.folds)
  cv.roc.stats <- rep(0, n.folds)
  cv.cor.stats <- rep(0, n.folds)
  cv.calibration.stats <- matrix(0, ncol=5, nrow = n.folds)
  if (family == "bernoulli") threshold.stats <- rep(0, n.folds)

  fitted.matrix <- matrix(NA, nrow = n.cases, ncol = n.folds)  # used to calculate se's
  fold.fit <- rep(0,n.cases)

  for (i in 1:n.folds) {

    pred.mask <- selector == i   #used to identify the with-held subset
    model.mask <- selector != i  #used to fit model on majority of data

    fits <- predict.gbm(model.list[[i]], x.data[model.mask, ], n.trees = target.trees)
    if (!is.null(offset)) fits <- fits + offset[model.mask]
    if (family == "bernoulli") fits <- exp(fits)/(1 + exp(fits))
    if (family == "poisson") fits <- exp(fits)
    fitted.matrix[model.mask,i] <- fits

    fits <- predict.gbm(model.list[[i]], x.data[pred.mask, ], n.trees = target.trees)
    if (!is.null(offset)) fits <- fits + offset[pred.mask]
    fold.fit[pred.mask] <- fits  # store the linear predictor values
    if (family == "bernoulli") fits <- exp(fits)/(1 + exp(fits))
    if (family == "poisson") fits <- exp(fits)
    fitted.matrix[pred.mask,i] <- fits

    y_i <- y.data[pred.mask]
    u_i <- fitted.matrix[pred.mask,i]  #pred.values[pred.mask]
    weight.preds <- site.weights[pred.mask]

    cv.deviance.stats[i] <- calc.deviance(y_i, u_i, weight.preds, family = family)

    cv.cor.stats[i] <- cor(y_i,u_i)

    if (family == "bernoulli") {
      cv.roc.stats[i] <- roc(y_i,u_i)
      cv.calibration.stats[i,] <- calibration(y_i,u_i,"binomial")
      threshold.stats[i] <- approx(ppoints(u_i), sort(u_i,decreasing = T), prevalence)$y
    }

    if (family == "poisson") {
      cv.calibration.stats[i,] <- calibration(y_i,u_i,"poisson")
    }
  }

  fitted.vars <- apply(fitted.matrix,1, var, na.rm = TRUE)

# now calculate the mean and se's for the folds

  cv.dev <- mean(cv.deviance.stats, na.rm = TRUE)
  cv.dev.se <- sqrt(var(cv.deviance.stats)) / sqrt(n.folds)

  cv.cor <- mean(cv.cor.stats, na.rm = TRUE)
  cv.cor.se <- sqrt(var(cv.cor.stats, use = "complete.obs")) / sqrt(n.folds)

  cv.roc <- 0.0
  cv.roc.se <- 0.0

  if (family == "bernoulli") {
    cv.roc <- mean(cv.roc.stats,na.rm=TRUE)
    cv.roc.se <- sqrt(var(cv.roc.stats, use = "complete.obs")) / sqrt(n.folds)
    cv.threshold <- mean(threshold.stats, na.rm = T)
    cv.threshold.se <- sqrt(var(threshold.stats, use = "complete.obs")) / sqrt(n.folds)
  }

  cv.calibration <- 0.0
  cv.calibration.se <- 0.0

  if (family == "poisson" | family == "bernoulli") {
    cv.calibration <- apply(cv.calibration.stats,2,mean)
    cv.calibration.se <- apply(cv.calibration.stats,2,var)
    cv.calibration.se <- sqrt(cv.calibration.se) / sqrt(n.folds) }

# fit the final model

  if (is.null(offset)) {
    gbm.call <- paste("gbm(y.data ~ .,data=x.data, n.trees = target.trees,
    interaction.depth = tree.complexity, shrinkage = learning.rate,
    bag.fraction = bag.fraction, weights = site.weights,
    distribution = as.character(family), var.monotone = var.monotone,
    verbose = FALSE)", sep="")
    }
  else {
    gbm.call <- paste("gbm(y.data ~ . + offset(offset),data=x.data, n.trees = target.trees,
    interaction.depth = tree.complexity, shrinkage = learning.rate,
    bag.fraction = bag.fraction, weights = site.weights,
    distribution = as.character(family), var.monotone = var.monotone,
    verbose = FALSE)", sep="")
    }

  if (!silent) cat("fitting final gbm model with a fixed number of ",target.trees," trees for ",sp.name,"\n")

  gbm.object <- eval(parse(text = gbm.call))

  best.trees <- target.trees

#extract fitted values and summary table

  gbm.summary <- summary(gbm.object,n.trees = target.trees, plotit = FALSE)

  fits <- predict.gbm(gbm.object,x.data,n.trees = target.trees)
  if (!is.null(offset)) fits <- fits + offset
  if (family == "bernoulli") fits <- exp(fits)/(1 + exp(fits))
  if (family == "poisson") fits <- exp(fits)
  fitted.values <- fits

  y_i <- y.data
  u_i <- fitted.values
  resid.deviance <- calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = FALSE)

  self.cor <- cor(y_i,u_i)
  self.calibration <- 0.0
  self.roc <- 0.0

  if (family == "bernoulli") {  # do this manually as we need the residuals
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    self.roc <- roc(y_i,u_i)
    self.calibration <- calibration(y_i,u_i,"binomial")
  }

  if (family == "poisson") {   # do this manually as we need the residuals
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    self.calibration <- calibration(y_i,u_i,"poisson")
  }

  if (family == "gaussian" | family == "laplace") {
    residuals <- y_i - u_i
  }

  mean.resid.deviance <- resid.deviance/n.cases

  z2 <- unclass(Sys.time())
  elapsed.time.minutes <- round((z2 - z1)/ 60,2)  #calculate the total elapsed time

  if (verbose) {
    cat("\n")
    cat("mean total deviance =", round(mean.total.deviance,3),"\n")
    cat("mean residual deviance =", round(mean.resid.deviance,3),"\n","\n")
    cat("estimated cv deviance =", round(cv.dev,3),"; se =",
      round(cv.dev.se,3),"\n","\n")
    cat("training data correlation =",round(self.cor,3),"\n")
    cat("cv correlation = ",round(cv.cor,3),"; se =",round(cv.cor.se,3),"\n","\n")
    if (family == "bernoulli") {
      cat("training data ROC score =",round(self.roc,3),"\n")
      cat("cv ROC score =",round(cv.roc,3),"; se =",round(cv.roc.se,3),"\n","\n")
    }
    cat("elapsed time - ",round(elapsed.time.minutes,2),"minutes","\n")
  }

  if (n.fitted == max.trees & !silent) {
    cat("\n","########### warning ##########","\n","\n")
    cat("maximum tree limit reached - results may not be optimal","\n")
    cat("  - refit with faster learning rate or increase maximum number of trees","\n")
  }

# now assemble data to be returned

  gbm.detail <- list(dataframe = dataframe.name, gbm.x = gbm.x, predictor.names = names(x.data),
    gbm.y = gbm.y, response.name = sp.name, offset = offset.name, family = family, tree.complexity = tree.complexity,
    learning.rate = learning.rate, bag.fraction = bag.fraction, cv.folds = n.folds,
    prev.stratification = prev.stratify, max.fitted = n.fitted, n.trees = target.trees,
    best.trees = target.trees, train.fraction = 1.0, tolerance.method = tolerance.method,
    tolerance = tolerance, var.monotone = var.monotone, date = date(),
    elapsed.time.minutes = elapsed.time.minutes)

  training.stats <- list(null = total.deviance, mean.null = mean.total.deviance,
    resid = resid.deviance, mean.resid = mean.resid.deviance, correlation = self.cor,
    discrimination = self.roc, calibration = self.calibration)

  cv.stats <- list(deviance.mean = cv.dev, deviance.se = cv.dev.se,
    correlation.mean = cv.cor, correlation.se = cv.cor.se,
    discrimination.mean = cv.roc, discrimination.se = cv.roc.se,
    calibration.mean = cv.calibration, calibration.se = cv.calibration.se)

  if (family == "bernoulli") {
    cv.stats$cv.threshold <- cv.threshold
    cv.stats$cv.threshold.se <- cv.threshold.se
  }

  rm(x.data,y.data, envir = globalenv())           #finally, clean up the temporary dataframes

# and assemble results for return

  gbm.object$gbm.call <- gbm.detail
  gbm.object$x <- data
  gbm.object$fitted <- fitted.values
  gbm.object$fitted.vars <- fitted.vars
  gbm.object$residuals <- residuals
  gbm.object$contributions <- gbm.summary
  gbm.object$self.statistics <- training.stats
  gbm.object$cv.statistics <- cv.stats
  gbm.object$weights <- site.weights
  gbm.object$trees.fitted <- trees.fitted
  gbm.object$training.loss.values <- training.loss.values
  gbm.object$cv.values <- cv.loss.values
  gbm.object$cv.loss.ses <- cv.loss.ses
  gbm.object$cv.loss.matrix <- cv.loss.matrix
  gbm.object$cv.roc.matrix <- cv.roc.stats

  if (keep.fold.models) gbm.object$fold.models <- model.list
  else gbm.object$fold.models <- NULL

  if (keep.fold.vector) gbm.object$fold.vector <- selector
  else gbm.object$fold.vector <- NULL

  if (keep.fold.fit) gbm.object$fold.fit <- fold.fit
  else gbm.object$fold.fit <- NULL

  return(gbm.object)
}

gbm.fixed <- function(
  data,                                # the input dataframe
  gbm.x,                               # indices of the predictors in the input dataframe
  gbm.y,                               # index of the response in the input dataframe
  tree.complexity = 1,                 # the tree depth - sometimes referred to as interaction depth
  site.weights = rep(1, nrow(data)),   # by default set equal
  verbose = TRUE,                      # to control reporting
  learning.rate = 0.001,               # controls speed of the gradient descent
  n.trees = 2000,                      # default number of trees
  train.fraction = 1,
  bag.fraction = 0.5,                  # varies random sample size for each new tree
  family = "bernoulli",                # can be any of bernoulli, poisson, gaussian, laplace - note quotes
  keep.data = FALSE,                   # keep original data
  var.monotone = rep(0, length(gbm.x)) # constrain to positive (1) or negative monontone (-1)
)
{
#
# j leathwick, j elith - 6th May 2007
#
# version 2.9 - developed in R 2.0
#
# calculates a gradient boosting (gbm)object with a fixed number of trees
# with the number of trees identified using gbm.step or some other procedure
#
# mostly used as a utility function, e.g., when being called by gbm.simplify
#
# takes as input a dataset and args selecting x and y variables, learning rate and tree complexity
#
# updated 13/6/05 to accommodate weighting of sites when calculating total and residual deviance
#
# updated 10/8/05 to correct how site.weights are returned
#
# requires gbm
#
#
  require(gbm)

  # setup input data and assign to position one

  dataframe.name <- deparse(substitute(data))   # get the dataframe name

  x.data <- eval(data[, gbm.x])                 # form the temporary datasets
  names(x.data) <- names(data)[gbm.x]
  y.data <- eval(data[, gbm.y])
  sp.name <- names(data)[gbm.y]


  assign("x.data", x.data, pos = 1)             #and assign them for later use
  assign("y.data", y.data, pos = 1)

  # fit the gbm model

  z1 <- unclass(Sys.time())

  gbm.call <- paste("gbm(y.data ~ .,n.trees = n.trees, data=x.data, verbose = F, interaction.depth = tree.complexity, weights = site.weights, shrinkage = learning.rate, distribution = as.character(family), var.monotone = var.monotone, bag.fraction = bag.fraction, keep.data = keep.data)", sep="")

  if (verbose) {
    print(paste("fitting gbm model with a fixed number of ",n.trees," trees for ",sp.name,sep=""),quote=FALSE)
  }

  gbm.object <- eval(parse(text = gbm.call))

  best.trees <- n.trees

  #extract fitted values and summary table

  fitted.values <- predict.gbm(gbm.object,x.data,n.trees = n.trees,type="response")
  gbm.summary <- summary(gbm.object,n.trees = n.trees, plotit = FALSE)

  y_i <- y.data
  u_i <- fitted.values

  if (family == "poisson") {
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    resid.deviance <- 2 * sum(deviance.contribs * site.weights)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)

    u_i <- sum(y.data * site.weights) / sum(site.weights)
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    total.deviance <- 2 * sum(deviance.contribs * site.weights)
  }

  if (family == "bernoulli") {
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    resid.deviance <- -2 * sum(deviance.contribs * site.weights)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)

    u_i <- sum(y.data * site.weights) / sum(site.weights)
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    total.deviance <- -2 * sum(deviance.contribs * site.weights)
  }

  if (family == "laplace") {
    resid.deviance <- sum(abs(y_i - u_i))
    residuals <- y_i - u_i
    u_i <- mean(y.data)
    total.deviance <- sum(abs(y_i - u_i))
  }

  if (family == "gaussian") {
    resid.deviance <- sum((y_i - u_i) * (y_i - u_i))
    residuals <- y_i - u_i
    u_i <- mean(y.data)
    total.deviance <- sum((y_i - u_i) * (y_i - u_i))
  }

  if (verbose) {
    print(paste("total deviance = ",round(total.deviance,2),sep=""),quote=F)
    print(paste("residual deviance = ",round(resid.deviance,2),sep=""),quote=F)
  }

  # now assemble data to be returned

  z2 <- unclass(Sys.time())
  elapsed.time.minutes <- round((z2 - z1)/ 60,2)  # calculate the total elapsed time

  gbm.detail <- list(dataframe = dataframe.name, gbm.x = gbm.x, predictor.names = names(x.data), gbm.y = gbm.y, reponse.name = names(y.data), family = family, tree.complexity = tree.complexity, learning.rate = learning.rate, bag.fraction = bag.fraction, cv.folds = 0, n.trees = n.trees, best.trees = best.trees, train.fraction = train.fraction, var.monotone = var.monotone, date = date(), elapsed.time.minutes = elapsed.time.minutes)

  gbm.object$gbm.call <- gbm.detail
  gbm.object$fitted <- fitted.values
  gbm.object$residuals <- residuals
  gbm.object$contributions <- gbm.summary
  gbm.object$self.statistics <- list(null.deviance = total.deviance, resid.deviance = resid.deviance)
  gbm.object$weights <- site.weights

  rm(x.data,y.data, pos=1)           # finally, clean up the temporary dataframes

  return(gbm.object)
}

gbm.cv <- function (
  data,                                # the input data frame
  gbm.x,                               # the predictors as numeric indices
  gbm.y,                               # the response
  n.trees = 500,                       # number of trees to fit
  tree.complexity = 1,                 # or interaction depth
  learning.rate = 0.001,               # controls rate of convergence - see gbm
  cv.folds = 5,                        # number of folds for cross validation
  family = "bernoulli",                # options as specified for gdm
  var.monotone = rep(0, length(gbm.x)),# see GDM documentation
  site.weights = rep(1, nrow(data)),
  verbose = TRUE,
  keep.data = FALSE
)
{
#
# j leathwick, j elith - October 2006
#
# version 2.9 - developed in R 2.3.1
#
# calculates a gradient boosting (gbm)object using Ridgeway's cross validation method
# this should be used with care, as performance is sometimes idisyncratic
# a safer option is to use gbm.step and manually select the required number of trees
#
# takes as input a dataset and args selecting x and y variables, and degree of interaction depth
#
# requires gbm
#
#
  require(gbm)

  # setup input data and assign to position one

  dataframe.name <- deparse(substitute(data))  # get the dataframe name

  x.data <- eval(data[, gbm.x])                 #form the temporary datasets
  names(x.data) <- names(data)[gbm.x]
  y.data <- eval(data[, gbm.y])
  sp.name <- names(data)[gbm.y]
  train.fraction <- 1

  assign("x.data", x.data, pos = 1)               #and assign them for later use
  assign("y.data", y.data, pos = 1)

  # fit the gbm model

  gbm.call <- paste("gbm(y.data ~ .,n.trees = n.trees, data=x.data, verbose = F, interaction.depth = tree.complexity, weights = site.weights, shrinkage = learning.rate, cv.folds = cv.folds, distribution = as.character(family), var.monotone = var.monotone, keep.data = keep.data)", sep="")

  if (verbose) {
    print(paste("fitting cross validation gbm model of ",n.trees," trees for ",sp.name, " with ",cv.folds," folds", sep=""),quote = FALSE)
  }

  gbm.object <- eval(parse(text = gbm.call))

  # identify the best number of trees using method appropriate to model

  best.trees <- gbm.perf(gbm.object, method = 'cv', plot.it = FALSE)

  if (best.trees == n.trees) {
    # we've failed to reach a satisfactory number...
    print(" WARNING - cross validation model failed to identify optimal number of trees",quote = FALSE)
    print(" re-fit with increased value for n.trees ", quote = FALSE)
  } else {
    print (paste("model successfully fitted, with ",best.trees," used out of a total of ", n.trees,sep=""), quote = FALSE)
  }

#extract fitted values and summary table

  fitted.values <- predict.gbm(gbm.object,x.data,n.trees = best.trees,type="response")
  gbm.summary <- summary(gbm.object,n.trees = best.trees, plotit = FALSE)

  y_i <- y.data
  u_i <- fitted.values

  if (family == "poisson") {
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    resid.deviance <- 2 * sum(deviance.contribs)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    u_i <- mean(y.data)
    total.deviance <- 2 * sum(ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i))
  }

  if (family == "bernoulli") {
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    resid.deviance <- -2 * sum(deviance.contribs)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    u_i <- mean(y.data)
    total.deviance <- -2 * sum((y_i * log(u_i)) + ((1-y_i) * log(1 - u_i)))
  }

  if (verbose) {
    print(paste("total deviance = ",round(total.deviance,2),sep=""),quote=F)
    print(paste("residual deviance = ",round(resid.deviance,2),sep=""),quote=F)
  }

  # now assemble data to be returned

  gbm.detail <- list(dataframe = dataframe.name, gbm.x = gbm.x, predictor.names = names(x.data),
    gbm.y = gbm.y, response.names = sp.name, degree= tree.complexity, n.trees = n.trees,
    learning.rate = learning.rate, best.trees = best.trees, cv.folds = cv.folds,
    family = family, train.fraction = train.fraction, var.monotone = var.monotone)

  gbm.object$gbm.call <- gbm.detail
  gbm.object$fitted <- fitted.values
  gbm.object$residuals <- residuals
  gbm.object$contributions <- gbm.summary
  gbm.object$deviances <- list(null.deviance = total.deviance, resid.deviance = resid.deviance)
  gbm.object$weights <- weights

  rm(x.data,y.data, pos=1)           #finally, clean up the temporary dataframes

  return(gbm.object)
}

gbm.oob <- function (data, gbm.x, gbm.y,interaction.depth = 1, site.weights = rep(1, nrow(data)), max.trees = 20000, verbose = TRUE, learning.rate = 0.001, n.trees = 200, add.trees = n.trees, cv.folds = 0, distribution = "bernoulli", train.fraction = 1, var.monotone = rep(0, length(gbm.x)), keep.data = TRUE)
{
#
# j leathwick, j elith - 23rd March 2005
#
# version 1.4 - developed in R 2.0
#
# calculates a gradient boosting (gbm)object in which an initial model is fitted
# and then trees are progressively added in sets of 10 testing performance
# along the way, using gbm.perf to test performance using out-of-bag estimation until the optimal
# number of trees is identified
#
# this method is not particularly recommended, as it is known to be very conservatice in its choice
# of number of trees - gbm.step should be used for small to moderate sized datasets or
# gbm.holdout can be used for larger datasets
#
# takes as input a dataset and args selecting x and y variables
#
# requires gbm
#
#
  require(gbm)

# setup input data and assign to position one

  dataframe.name <- deparse(substitute(data))  # get the dataframe name

  x.data <- eval(data[, gbm.x])                 #form the temporary datasets
  names(x.data) <- names(data)[gbm.x]
  y.data <- eval(data[, gbm.y])
  sp.name <- names(data)[gbm.y]


  assign("x.data", x.data, pos = 1)               #and assign them for later use
  assign("y.data", y.data, pos = 1)

# fit the gbm model

  if (verbose) {
    print(paste("fitting initial gbm model of ",n.trees," trees for ",sp.name, " and expanding using OOB method",sep=""),quote=FALSE)
  }

  gbm.call <- paste("gbm(y.data ~ .,n.trees = n.trees, data=x.data, verbose = F, interaction.depth = interaction.depth, weights = site.weights, shrinkage = learning.rate, distribution = as.character(distribution), var.monotone = var.monotone, keep.data = keep.data)", sep="")

  gbm.object <- eval(parse(text = gbm.call))

# identify the best number of trees using method appropriate to model

  best.trees <- gbm.perf(gbm.object, method == 'OOB', plot.it = FALSE)

  n.fitted <- n.trees

  if (verbose) {
    print("expanding model to find optimal no of trees...",quote=FALSE)
  }

  while(gbm.object$n.trees - best.trees < n.trees & n.fitted < max.trees) {
    gbm.object <- gbm.more(gbm.object, add.trees)
    n.fitted <- n.fitted + add.trees
    best.trees <- gbm.perf(gbm.object, plot.it = FALSE)
    if (n.fitted %% 100 == 0) {
      #report times along the way
      if (verbose) {
        print(paste("fitted trees = ", n.fitted, sep = ""), quote = FALSE)
      }
    }
  }

  if (verbose) {
    print(paste("fitting stopped at ",best.trees," trees",sep=""),quote=FALSE)
  }

  # extract fitted values and summary table

  fitted.values <- predict.gbm(gbm.object,x.data,n.trees = n.trees,type="response")
  y_i <- y.data
  u_i <- fitted.values

  if (distribution == "poisson") {
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    resid.deviance <- 2 * sum(deviance.contribs)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    u_i <- mean(y.data)
    total.deviance <- 2 * sum(ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i))
  }

  if (distribution == "bernoulli") {
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    resid.deviance <- -2 * sum(deviance.contribs)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    u_i <- mean(y.data)
    total.deviance <- -2 * sum((y_i * log(u_i)) + ((1-y_i) * log(1 - u_i)))
  }

  if (verbose) {
    print(paste("total deviance = ",round(total.deviance,2),sep=""),quote=F)
    print(paste("residual deviance = ",round(resid.deviance,2),sep=""),quote=F)
  }

  gbm.summary <- summary(gbm.object,n.trees = best.trees, plotit = FALSE)

  # now assemble data to be returned

  gbm.detail <- list(dataframe = dataframe.name, gbm.x = gbm.x, x.names = names(x.data),
      gbm.y = gbm.y, y.names = names(y.data), degree= interaction.depth, n.trees = n.fitted,
      learning.rate = learning.rate, best.trees = best.trees, cv.folds = 0,
      distribution = distribution, train.fraction = train.fraction, var.monotone = var.monotone)

  rm(x.data,y.data, pos=1)           #finally, clean up the temporary dataframes

  return(list(gbm.object = gbm.object, fitted.values = fitted.values, residuals = residuals, summary = gbm.summary, deviance = deviance, deviances = list(null.deviance = total.deviance, resid.deviance = resid.deviance), weights = weights, gbm.call = gbm.detail))
}

gbm.holdout <- function (
   data,                               # the input data frame
   gbm.x,                              # indices of predictor variables
   gbm.y,                              # index of response variable
   learning.rate = 0.001,              # typically varied between 0.1 and 0.001
   tree.complexity = 1,                # sometimes called interaction depth
   family = "bernoulli",               # "bernoulli","poisson", etc. as for gbm
   n.trees = 200,                      # initial number of trees
   add.trees = n.trees,                # number of trees to add at each increment
   max.trees = 20000,                  # maximum number of trees to fit
   verbose = TRUE,                     # controls degree of screen reporting
   train.fraction = 0.8,               # proportion of data to use for training
   permute = TRUE,                     # reorder data to start with
   prev.stratify = TRUE,               # stratify selection for p/a data
   var.monotone = rep(0, length(gbm.x)),# allows constraining of response to monotone
   site.weights = rep(1, nrow(data)),  # set equal to 1 by default
   refit = TRUE,                       # refit the model with the full data but id'd no of trees
   keep.data = TRUE                    # keep copy of the data
)
{
#
# j leathwick, j elith - October 2006
#
# version 2.9 - developed in R 2.3.1
#
# calculates a gradient boosting (gbm)object in which model complexity is
# determined using a training set with predictions made to a withheld set
# an initial set of trees is fitted, and then trees are progressively added
# testing performance # along the way, using gbm.perf until the optimal
# number of trees is identified
#
# as any structured ordering of the data should be avoided, a copy of the data set
# BY DEFAULT is randomly reordered each time the function is run
#
# takes as input a dataset and args selecting x and y variables, and degree of interaction depth
#
# requires gbm
#
#
  require(gbm)

# setup input data and assign to position one

  dataframe.name <- deparse(substitute(data))  # get the dataframe name
  cv.folds <- 0

  if (permute) {
    print("",quote=FALSE)
    print("WARNING - data is being randomly reordered to avoid confounding effects",quote=FALSE)
    print("of inherent structure as submitted - use permute = FALSE to turn off this option",quote=FALSE)
    n.rows <- nrow(data)

    if (prev.stratify == TRUE & family == "bernoulli") {

      presence.mask <- data[,gbm.y] == 1
      absence.mask <- data[,gbm.y] == 0
      n.pres <- sum(presence.mask)
      n.abs <- sum(absence.mask)

      selector <- seq(1,n.rows)

      temp <- sample(selector[presence.mask],size = n.pres * train.fraction)
      selector[temp] <- 0

      temp <- sample(selector[absence.mask],size = n.abs * train.fraction)
      selector[temp] <- 0

      sort.vector <- sort(selector,index.return = TRUE)[[2]]
    } else {
      sort.vector <- sample(seq(1,n.rows),n.rows,replace=FALSE)
    }

    sort.data <- data[sort.vector,]

    x.data <- eval(sort.data[, gbm.x])                 #form the temporary datasets
    y.data <- eval(sort.data[, gbm.y])
  } else {
    x.data <- eval(data[, gbm.x])                 #form the temporary datasets
    y.data <- eval(data[, gbm.y])
  }

  names(x.data) <- names(data)[gbm.x]
  sp.name <- names(data)[gbm.y]

  assign("x.data", x.data, pos = 1)               #and assign them for later use
  assign("y.data", y.data, pos = 1)

# fit the gbm model

  print(paste("fitting initial gbm model of ",n.trees," trees for ",sp.name,sep=""),quote=FALSE)
  print(" and expanding using withheld data for evaluation",quote=FALSE)

  gbm.call <- paste("gbm(y.data ~ .,n.trees = n.trees, data=x.data, verbose = F, interaction.depth = tree.complexity, weights = site.weights, shrinkage = learning.rate, cv.folds = 0, distribution = as.character(family), train.fraction = train.fraction, var.monotone = var.monotone, keep.data = keep.data)", sep="")

  gbm.object <- eval(parse(text = gbm.call))

# identify the best number of trees using method appropriate to model

  best.trees <- gbm.perf(gbm.object, method = 'test', plot.it = FALSE)

  n.fitted <- n.trees

  if (verbose) print("expanding model to find optimal no of trees...",quote=FALSE)

  while(gbm.object$n.trees - best.trees < n.trees & n.fitted < max.trees){

    gbm.object <- gbm.more(gbm.object, add.trees)
    best.trees <- gbm.perf(gbm.object, method = 'test', plot.it = FALSE)
    n.fitted <- n.fitted + add.trees

    if (n.fitted %% 100 == 0){ 
      #report times along the way
      if (verbose) print(paste("fitted trees = ", n.fitted, sep = ""), quote = FALSE)
    }
  }

  if (verbose) print(paste("fitting stopped at ",best.trees," trees",sep=""),quote=FALSE)

  if (refit) {
    # we are refitting the model with fixed tree size
    print(paste("refitting the model to the full dataset using ",best.trees," trees",sep=""),quote=FALSE)

    x.data <- eval(data[, gbm.x])                 #form the temporary datasets
    y.data <- eval(data[, gbm.y])

    gbm.call <- eval(paste("gbm(y.data ~ .,n.trees = best.trees, data=x.data, verbose = F, interaction.depth = tree.complexity, weights = site.weights, shrinkage = learning.rate, cv.folds = 0, distribution = as.character(family), var.monotone = var.monotone, keep.data = keep.data)", sep=""))

    gbm.object <- eval(parse(text = gbm.call))

  }

#extract fitted values and summary table

  fitted.values <- predict.gbm(gbm.object,x.data,n.trees = best.trees,type="response")
  gbm.summary <- summary(gbm.object,n.trees = best.trees, plotit = FALSE)

  y_i <- y.data
  u_i <- fitted.values

  if (family == "poisson") {
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    resid.deviance <- 2 * sum(deviance.contribs)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    u_i <- mean(y.data)
    total.deviance <- 2 * sum(ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i))
  }

  if (family == "bernoulli") {
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    resid.deviance <- -2 * sum(deviance.contribs)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    u_i <- mean(y.data)
    total.deviance <- -2 * sum((y_i * log(u_i)) + ((1-y_i) * log(1 - u_i)))
  }

  if (verbose) {
    print(paste("total deviance = ",round(total.deviance,2),sep=""),quote=F)
    print(paste("residual deviance = ",round(resid.deviance,2),sep=""),quote=F)
  }

# now assemble data to be returned

  gbm.detail <- list(dataframe = dataframe.name, gbm.x = gbm.x, predictor.names = names(x.data),
    gbm.y = gbm.y, response.name = sp.name, tree.complexity = tree.complexity, n.trees = best.trees,
    learning.rate = learning.rate, best.trees = best.trees, cv.folds = cv.folds,
    family = family, train.fraction = train.fraction, var.monotone = var.monotone )

  gbm.object$fitted <- fitted.values
  gbm.object$residuals <- residuals
  gbm.object$contributions <- gbm.summary
  gbm.object$deviances <- list(null.deviance = total.deviance, resid.deviance = resid.deviance)
  gbm.object$weights <- weights
  gbm.object$gbm.call <- gbm.detail

  rm(x.data,y.data, pos=1)           #finally, clean up the temporary dataframes

  return(gbm.object)
}

gbm.kfold <- function (
   gbm.object,
   n.folds = 10,
   prev.stratify = F,
   fold.vector = NULL,
   strict = T,
   verbose = TRUE,
   seed = NULL
)
{
#
# j. leathwick/j. elith - 24th March 2004
#
# version 1.4
#
# function to perform k-fold cross validation
# with optional full model perturbation for each fold
#
# requires gbm library from Cran
# requires functions roc and calibration of j. elith
#
# takes a gbm object
# and first assesses the full model, and then
# randomly subsets the dataset into n.folds folds and drops
# each fold in turn, fitting on remaining data
# and predicting for withheld data
#
# calculate rocs and calibration on folds as well as full data
# returning the mean and se of the ROC scores
# and the mean calibration statistics
#
# modified 8/10/04 to
#   1. add prevalence stratification
#
# modified 28 Jan 05 to:
#   1. accommodate cross validation fitting
#   2. fit folds model using fit.gbm
#   3. use only the original data set, obviating need for storage of temporary data sets
#
# modified 2nd March to accommodate different data distributions but still requires residual deviance for normal data
#
# modified 16th June to accommodate site weights

    require(gbm)

#    gbm.terms <- gbm.object$Terms       #and the gbm call details

    gbm.call <- gbm.object$gbm.call
    n.preds <- length(gbm.call$pred.names)
    data <- eval(parse(text = gbm.call$dataframe))
    gbm.x <- gbm.call$gbm.x
    gbm.y <- gbm.call$gbm.y
    sp.name <- names(data)[gbm.y]

    tree.complexity <- gbm.call$tree.complexity
    learning.rate <- gbm.call$learning.rate
    n.trees <- gbm.call$n.trees
    best.trees <- gbm.call$best.trees
    family <- gbm.call$family
    train.fraction <- gbm.call$train.fraction
    site.weights <- eval(gbm.object$weights)
    var.monotone <- gbm.call$var.monotone

    n.cases <- nrow(data)

    if (verbose) cat("Calculating statistics from full model for",sp.name,"\n")
    if (verbose) cat("using family of",family,"\n")

    if (strict) {
        if (verbose) cat(" using re-estimation of optimal tree size","\n")
    }
    else {
        if (verbose) cat("  but with fixed numbers of trees","\n")
    }

    y_i <- data[, gbm.y]
    u_i <- gbm.object$fitted

    self.resid.deviance <- calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = TRUE)
    self.correlation <- cor(y_i, u_i)

    if (family=="bernoulli") {
        self.discrimination <- roc(y_i, u_i)
        self.calibration <- calibration(y_i, u_i)
	}

    if (family=="poisson") {
        self.calibration <- calibration(y_i, u_i, family = "poisson")
	}

# set up for results storage

    fold.resid.deviance <- rep(0,n.folds)
    fold.discrimination <- rep(0,n.folds)
    fold.correlation <- rep(0,n.folds)
    fold.calibration <- as.data.frame(matrix(0,ncol=5,nrow=n.folds))
    names(fold.calibration) <- c("intercept","slope","test1","test2","test3")

    predictions <- rep(0,n.cases)

# now setup for withholding random folds

    pred.values <- rep(0, n.cases)

    if (is.null(fold.vector)) {
      if (prev.stratify) {

        presence.mask <- data[,gbm.y] == 1
        absence.mask <- data[,gbm.y] == 0
        n.pres <- sum(presence.mask)
        n.abs <- sum(absence.mask)

        selector <- rep(0,n.cases)
        #create a vector of randomised numbers and feed into presences
        temp <- rep(seq(1, n.folds, by = 1), length = n.pres)
        temp <- temp[order(runif(n.pres, 1, 100))]
        selector[presence.mask] <- temp

        #and then do the same for absences
        temp <- rep(seq(1, n.folds, by = 1), length = n.abs)
        temp <- temp[order(runif(n.abs, 1, 100))]
        selector[absence.mask] <- temp

      } else {
		#otherwise make them random with respect to presence/absence
        selector <- rep(seq(1, n.folds, by = 1), length = n.cases)
        selector <- selector[order(runif(n.cases, 1, 100))]
      }
    }
    if (!is.null(fold.vector)) {
      selector <- fold.vector
      if (verbose) cat("found saved fold vector and using this for all new models","\n","\n")
    }
    if (verbose) cat("Now processing folds...","\n")

    for (i in 1:n.folds) {
        if (verbose) cat(i," ")
        model.mask <- selector != i  #used to fit model on majority of data
        pred.mask <- selector == i   #used to identify the with-held fold

    # fit new gbm model

        if (!is.null(seed)) set.seed(seed)

        if (!strict) {
			# we have a fixed model
            gbmmod.new <- gbm.fixed(data = data[model.mask, ], gbm.x = gbm.x, gbm.y = gbm.y, tree.complexity = tree.complexity, n.trees = best.trees, learning.rate = learning.rate, var.monotone = var.monotone, family = family, verbose = FALSE)
        }

        best.trees <- gbmmod.new$gbm.call$best.trees  #only update this for strict cross validation with full model refit

        # print(paste("step ",i," model fitted with ",best.trees," trees",sep=""),quote = FALSE)

        predictions[pred.mask] <- predict.gbm(gbmmod.new, data[pred.mask, gbm.x], type = "response", n.trees = best.trees)

        y_i <- data[pred.mask, gbm.y]
        u_i <- predictions[pred.mask]

        weights.fold <- site.weights[pred.mask]

        fold.resid.deviance[i] <- calc.deviance(y_i, u_i, weights = weights.fold, family = family, calc.mean = TRUE)
        fold.correlation[i] <- cor(y_i,u_i)
        fold.calibration[i,] <- calibration(y_i, u_i, family = family)

        if (family=="bernoulli") {
           fold.discrimination[i] <- roc(y_i, u_i)
        }
       # fold.calibration[i,] <- calibration(y_i, u_i)}
       # 
       # if (family=="poisson"){
       #     fold.calibration[i,] <- calibration(y_i, u_i, family = "poisson")
       # }
    }


# and assemble results for return

    if (verbose) cat("","\n")   #send a return to the screen first

    self.stats <- list(self.resid.deviance = self.resid.deviance, self.correlation = self.correlation, self.discrimination = self.discrimination, self.calibration = self.calibration)

    y_i <- data[, gbm.y]
    u_i <- predictions

    pooled.resid.deviance <- calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = TRUE)
    pooled.correlation <- cor(u_i, u_i)
    pooled.calibration <- calibration(y_i, u_i, family = family)

    if (family=="bernoulli") {
      pooled.discrimination <- roc(y_i, u_i)
    } else {
      pooled.discrimination <- NA
    }
	
    pooled.stats <- list(pooled.deviance = pooled.resid.deviance, pooled.correlation = pooled.correlation, pooled.discrimination = pooled.discrimination, pooled.calibration = pooled.calibration)

    fold.deviance.mean <- mean(fold.resid.deviance)
    fold.deviance.se <- sqrt(var(fold.resid.deviance))/sqrt(n.folds)

    fold.correlation.mean <- mean(fold.correlation)
    fold.correlation.se <- sqrt(var(fold.correlation))/sqrt(n.folds)

    if (family == "bernoulli") {
      fold.discrimination.mean <- mean(fold.discrimination)
      fold.discrimination.se <- sqrt(var(fold.discrimination))/sqrt(n.folds)
    } else {
      fold.discrimination.mean <- NA
      fold.discrimination.se <- NA
    }

    fold.calibration.mean <- apply(fold.calibration,2,mean)
    names(fold.calibration.mean) <- names(fold.calibration)

    fold.calibration.se <- apply(fold.calibration,2,var)
    fold.calibration.se <- sqrt(fold.calibration.se) / sqrt(n.folds)

    cv.stats <- list(deviance.mean = fold.deviance.mean, deviance.se = fold.deviance.se, correlation.mean = fold.correlation.mean, correlation.se = fold.correlation.se, discrimination.mean = fold.discrimination.mean, discrimination.se = fold.discrimination.se, calibration.mean = fold.calibration.mean, calibration.se = fold.calibration.se )

    fold.values = list(fold.deviance = fold.resid.deviance, fold.correlation = fold.correlation, fold.discrimination = fold.discrimination, fold.calibration = fold.calibration)

    return(list(gbm.call = gbm.call, self.statistics = self.stats, pooled.statistics = pooled.stats, cv.statistics = cv.stats, fold.values = fold.values))
}

fit.gbm <- function(data, gbm.x, gbm.y,interaction.depth = 1, site.weights = rep(1, nrow(data)), max.trees = 20000, verbose = TRUE, learning.rate = 0.001, init.trees = 200, cv.folds = 1, prev.stratify = FALSE, distribution = "bernoulli",keep.data = TRUE)
{
#
# j leathwick, j elith - 7th January 2005
#
# version 1.2 - developed in R 2.0
#
# calculates a gradient boosting (gbm)object in which an initial model is fitted
# and then trees are progressively added in sets of 10 testing performance
# along the way, using gbm.perf to test performance until the optimal
# number of trees is identified
#
# takes as input a dataset and args selecting x and y variables, and degree of interaction depth
#
# requires gbm
#
# modified 28 Jan 05 to allow:
#   1. fitting using cross validation
#   2. more agressive forward searching when using OOB to determine optimum number of trees
#
# modified 2nd March 2005 to accommodate different error distributions
#
# replaced by gbm.step, gbm.fixed, gbm.cv, and gbm.oob which should be used in preference - 23rd March 2005
#
    require(gbm)

# setup input data and assign to position one

    dataframe.name <- deparse(substitute(data))  # get the dataframe name

    x.data <- eval(data[, gbm.x])                 #form the temporary datasets
    names(x.data) <- names(data)[gbm.x]
    y.data <- eval(data[, gbm.y])
    sp.name <- names(data)[gbm.y]


    assign("x.data", x.data, pos = 1)               #and assign them for later use
    assign("y.data", y.data, pos = 1)

# fit the gbm model

    gbm.call <- paste("gbm(y.data ~ .,n.trees = init.trees, data=x.data, verbose = F, interaction.depth = interaction.depth, weights = site.weights, shrinkage = learning.rate, cv.folds = cv.folds, distribution = as.character(distribution), keep.data = keep.data)", sep="")

    if (verbose) {
      if (cv.folds == 1) {
        print(paste("fitting initial gbm model of ",init.trees," trees for ",sp.name, " and expanding using OOB method",sep=""),quote=FALSE)
      } else {
        print(paste("fitting cross validation gbm model of ",init.trees," trees for ",sp.name, " with ",cv.folds," folds",sep=""),quote = FALSE)
      }
    }

    gbm.object <- eval(parse(text = gbm.call))

# identify the best number of trees using method appropriate to model

    if (cv.folds == 1) {
       best.trees <- gbm.perf(gbm.object, method == 'OOB', plot.it = FALSE)
    } else {
       best.trees <- gbm.perf(gbm.object, method = 'cv', plot.it = FALSE)
    }

    if (cv.folds == 1) {
      # a non cross-validation model

        n.fitted <- init.trees

        if (verbose) print("expanding model to find optimal no of trees...",quote=FALSE)

        while(gbm.object$n.trees - best.trees < init.trees & n.fitted < max.trees) {
            gbm.object <- gbm.more(gbm.object, init.trees)
            best.trees <- gbm.perf(gbm.object, plot.it = FALSE)
            n.fitted <- n.fitted + init.trees
            if (n.fitted %% 100 == 0) {
              # report times along the way
              if (verbose) print(paste("fitted trees = ", n.fitted, sep = ""), quote = FALSE)
            }
        }
        if (verbose) print(paste("fitting stopped at ",n.fitted," trees",sep=""),quote=FALSE)
    } else {
      # we're doing a cross validation model

        if (best.trees == init.trees) {
            # we've failed to reach a satisfactory number...
            print(" WARNING - cross validation model failed to identify optimal number of trees",quote = FALSE)
            print(" re-fit with increased value for init.trees ", quote = FALSE)
        } else {
            if (verbose) print (paste("model successfully fitted, with ",best.trees," used out of a total of ", init.trees,sep=""), quote = FALSE)
        }
    }

#extract fitted values and summary table

    fitted.values <- predict.gbm(gbm.object,x.data,n.trees = init.trees,type="response")
    gbm.summary <- summary(gbm.object,n.trees = init.trees, plotit = FALSE)

# now assemble data to be returned

    gbm.detail <- list(dataframe = dataframe.name, gbm.x = gbm.x, x.names = names(x.data),
        gbm.y = gbm.y, y.names = names(y.data), degree= interaction.depth, init.trees = init.trees,
        learning.rate = learning.rate, best.trees = best.trees, cv.folds = cv.folds, distribution = distribution)

    rm(x.data,y.data, pos=1)           #finally, clean up the temporary dataframes

    return(list(gbm.object = gbm.object, fitted.values = fitted.values, summary = gbm.summary, weights = weights, gbm.call = gbm.detail))
}

gbm.anz <- function(formula = formula(data), distribution = "bernoulli", data = list(), weights, offset = NULL, var.monotone = NULL, n.trees = 100, interaction.depth = 1, n.minobsinnode = 10, shrinkage = 0.001, bag.fraction = 0.5, train.fraction = 1, cv.folds = 0, keep.data = TRUE, verbose = TRUE)
{
    call <- match.call()
    m <- match.call(expand = FALSE)
    m$distribution <- m$offset <- m$var.monotone <- m$n.trees <- NULL
    m$interaction.depth <- m$n.minobsinnode <- m$shrinkage <- NULL
    m$bag.fraction <- m$train.fraction <- m$keep.data <- m$verbose <- NULL
    m$cv.folds <- NULL
    m[[1]] <- as.name("model.frame")
    m$na.action <- na.pass
    m.keep <- m
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    a <- attributes(Terms)
    y <- model.extract(m, response)
    x <- model.frame(delete.response(Terms), data, na.action = na.pass)
    w <- model.extract(m, weights)
    if (!is.null(model.extract(m, offset))) {
        stop("In gbm the offset term needs to be specified using the offset parameter.")
    }
    var.names <- a$term.labels
    response.name <- dimnames(attr(terms(formula), "factors"))[[1]][1]
    cv.error <- NULL
    if (cv.folds > 1) {
        i.train <- 1:floor(train.fraction * length(y))

        presence.mask <- y[1:length(i.train)] == 1  #next 10 rows stratify by prevalence - J Leathwick - 28 Jan 05
        absence.mask  <- y[1:length(i.train)] == 0
        n.pres <- sum(presence.mask)
        n.abs <- sum(absence.mask)
        cv.group <- rep(0,length(i.train))

        #insert randomised numbers into presence rows
        cv.group[presence.mask] <- sample(rep(1:cv.folds, length = n.pres))

        #and then do the same for absences
        cv.group[absence.mask] <-  sample(rep(1:cv.folds, length = n.abs))

        if (verbose) print(table(cv.group,y))
        cv.error <- rep(0, n.trees)
        for (i.cv in 1:cv.folds) {
            if (verbose)
                cat("CV:", i.cv, "\n")
            i <- order(cv.group == i.cv)
            gbm.obj <- gbm.fit(x[i.train, ][i, ], y[i.train][i],
                offset = offset[i.train][i], distribution = distribution,
                w = ifelse(w == NULL, NULL, w[i.train][i]), var.monotone = var.monotone,
                n.trees = n.trees, interaction.depth = interaction.depth,
                n.minobsinnode = n.minobsinnode, shrinkage = shrinkage,
                bag.fraction = bag.fraction, train.fraction = mean(cv.group != i.cv),
                keep.data = FALSE, verbose = verbose,
                var.names = var.names, response.name = response.name
            )
            cv.error <- cv.error + gbm.obj$valid.error * sum(cv.group == i.cv)
        }
        cv.error <- cv.error/length(i.train)
    }
    gbm.obj <- gbm.fit(x, y, offset = offset, distribution = distribution,
        w = w, var.monotone = var.monotone, n.trees = n.trees,
        interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode,
        shrinkage = shrinkage, bag.fraction = bag.fraction, train.fraction = train.fraction,
        keep.data = keep.data, verbose = verbose, var.names = var.names,
        response.name = response.name
    )
    gbm.obj$Terms <- Terms
    gbm.obj$cv.error <- cv.error
    if (!keep.data) {
        gbm.obj$m <- m.keep
    }
    return(gbm.obj)
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
  n.reps = 200,                 # number of bootstrap samples
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
  train.data <- eval(parse(text=gbm.call$dataframe))
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

  if (bootstrap.functions) {  # then set up a dataframe of means for predicting them

    function.matrix <- array(0,dim=c(n.reps,n.preds,200))

    function.pred.frame <- as.data.frame(as.vector(rep(mean(train.data[,gbm.x[1]],na.rm=T),200)))

    for (i in 2:n.preds) {  # step through the predictor set

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
  if (bootstrap.predictions) cat("prediction dataset has ",n.pred.obs," rows\n\n",sep="")
    else cat("no prediction dataset provided...\n\n")

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

      for (j in 1:n.preds) {  #cycle through the first time and get the range of the functions

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
      }
      else {
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
      train.fit.upper[i] <- approx(ppoints(temp),sort(temp),0.975)$y
      train.fit.lower[i] <- approx(ppoints(temp),sort(temp),0.025)$y
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
      train.pred.upper[i] <- approx(ppoints(temp),sort(temp),0.975)$y
      train.pred.lower[i] <- approx(ppoints(temp),sort(temp),0.025)$y
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
      names(function.dataframe)[j:(j+3)] <-
        paste(predictor.names[i],c(".vals",".lower",".mean",".upper"),sep="")
    }

    for (i in 1:n.preds) {
      k <- (i * 4) - 3
      if (is.vector(train.data[,gbm.x[i]])) {
        function.dataframe[,k] <- seq(min(train.data[,gbm.x[i]], na.rm = T),
          max(train.data[,gbm.x[i]], na.rm = T),length = 200)
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
    final.object$train.fit.stats <- data.frame(fit.mean = train.fit.mean, fit.lower = train.fit.lower,
     fit.upper = train.fit.upper, fit.count = train.fit.count)
    final.object$train.pred.stats <- data.frame(pred.mean = train.pred.mean, pred.lower = train.pred.lower,
     pred.upper = train.pred.upper, pred.count = train.pred.count)
    }

  if (bootstrap.functions) final.object$function.dataframe <- function.dataframe

  if (bootstrap.predictions) {
    final.object$prediction.stats = data.frame(pred.mean = preds.mean, upper.limit = preds.upper.limit,
     lower.limit = preds.lower.limit)
    if (return.tails) {
      final.object$preds.lower.tail <- preds.lower.tail
      final.object$preds.upper.tail <- preds.upper.tail
    }
    if(return.pred.matrix) final.object$pred.matrix <- boot.preds
  }

  return(final.object)
}

gbm.plot <- function(
     gbm.object,                    # a gbm object - could be one from gbm.step
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
     trans = NA,                    # if the results need to be transformed prior to plotting
     ...                            # other arguments to pass to the plotting
                                    # useful options include cex.axis, cex.lab, etc.
)
{
# function to plot gbm response variables, with the option
# of adding a smooth representation of the response if requested
# additional options in this version allow for plotting on a common scale
# note too that fitted functions are now centered by subtracting their mean
#
# version 2.9
#
# j. leathwick/j. elith - March 2007
#

require(gbm)
require(splines)

gbm.call <- gbm.object$gbm.call
gbm.x <- gbm.call$gbm.x
pred.names <- gbm.call$predictor.names
response.name <- gbm.call$response.name
data <- gbm.object$x

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

predictors <- list(rep(NA,n.plots)) # matrix(0,ncol=n.plots,nrow=100)
responses <- list(rep(NA,n.plots)) # matrix(0,ncol=n.plots,nrow=100)

#cycle through the first time and get the range of the functions
for (j in c(1:n.plots)) {  
  if (n.plots == 1) {
    k <- variable.no
  } else {
    k <- match(gbm.object$contributions$var[j],pred.names)
  }

  if (is.null(x.label)) {
    var.name <- gbm.call$predictor.names[k]
  } else {
    var.name <- x.label
  }
  
  pred.data <- data[,gbm.call$gbm.x[k]]

  response.matrix <- plot.gbm(gbm.object, k, return.grid = TRUE)

  predictors[[j]] <- response.matrix[,1]
  if (is.factor(data[,gbm.call$gbm.x[k]])) {
    predictors[[j]] <- factor(predictors[[j]],levels = levels(data[,gbm.call$gbm.x[k]]))
  }

  ###############
  # here is the SM change
  #responses[[j]] <- response.matrix[,2] - mean(response.matrix[,2])
  if (is.na(trans)) {
    responses[[j]] <- response.matrix[,2] - mean(response.matrix[,2])
  } else {
    responses[[j]] <- trans(response.matrix[,2])
  }

  #end of SM change
  ###############


  if(j == 1) {
    ymin = min(responses[[j]])
    ymax = max(responses[[j]])
  } else {
    ymin = min(ymin,min(responses[[j]]))
    ymax = max(ymax,max(responses[[j]]))
  }
}

# now do the actual plots

  for (j in c(1:n.plots)) {

   if (plot.count == max.plots) {
     plot.count = 0
     n.pages <- n.pages + 1
   }

   if (plot.count == 0) {
     dev.new(width = 11, height = 8)
     par(mfrow = plot.layout)
   }

    plot.count <- plot.count + 1

    if (n.plots == 1) {
      k <- match(pred.names[variable.no],gbm.object$contributions$var)
      if (show.contrib) {
         x.label <- paste(var.name,"  (",round(gbm.object$contributions[k,2],1),"%)",sep="")
      }
    } else {
      k <- match(gbm.object$contributions$var[j],pred.names)
      var.name <- gbm.call$predictor.names[k]
      if (show.contrib) {
         x.label <- paste(var.name,"  (",round(gbm.object$contributions[j,2],1),"%)",sep="")
      } else {
        x.label <- var.name
      }
    }

    if (common.scale) {
      plot(predictors[[j]],responses[[j]],ylim=c(ymin,ymax), type='l', xlab = x.label, ylab = y.label, ...)
    } else {
      plot(predictors[[j]],responses[[j]], type='l', xlab = x.label, ylab = y.label, ...)
    }
    if (smooth & is.vector(predictors[[j]])) {
      temp.lo <- loess(responses[[j]] ~ predictors[[j]], span = 0.3)
      lines(predictors[[j]],fitted(temp.lo), lty = 2, col = 2)
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
    if (rug & is.vector(data[,gbm.call$gbm.x[k]])) {
      rug(quantile(data[,gbm.call$gbm.x[k]], probs = seq(0, 1, 0.1), na.rm = TRUE),lwd=2)
    }
  }
}

gbm.unipred <- function(
  gbm.object,
  x = NULL,           # the predictor variable to plot - NULL to force a value
  x.label = NULL,     # allows manual specification of the x label
  x.range = NULL,     # manual range specification for the x variable
  y.label = NULL,     # allows manual specification of the y label
  y.range = NULL,     # manual range specification for the predicted response
  pred.means = NULL,  # allows specification of values for other variables in a list
  smooth = FALSE,     # controls smoothing of the response function
  mask = FALSE,       # allows a mask object to be used to restrict the plotting range
  mask.object = NULL,
  add = FALSE,        # controls whether a new plot will be drawn
  return.data = FALSE,# should the predicted values be returned
  ...                 # additional plot parameters
)
{
#
# gbm.unipred version 2.9 August 2006
# J Leathwick/J Elith
#
# takes a gbm boosted regression tree object produced by gbm.step and
# plots a perspective plot showing predicted values for two predictors
# as specified by number using x and y
# values for all other variables are set at their mean by default
# but values can be specified by giving a list consisting of the variable name
# and its desired value, e.g., c(name1 = 12.2, name2 = 57.6)
#
  require(gbm)
  require(splines)

# then get the model parameters

  gbm.call <- gbm.object$gbm.call
  family <- gbm.call$family
  gbm.x <- gbm.call$gbm.x
  n.preds <- length(gbm.x)
  n.trees <- gbm.object$gbm.call$best.trees
  data <- eval(parse(text=gbm.call$dataframe))[,gbm.x]
  pred.names <- gbm.call$predictor.names

  if (is.null(x.label)) {
    x.label <- gbm.call$predictor.names[x]
  }

  if (is.null(y.label)) {
    y.label <- "fitted response"
  }

#first check variables specified in pred.means

  if (!is.null(pred.means)) {
    for (k in 1:length(pred.means)) {
      m <- match(names(pred.means)[k],pred.names)
      if (is.na(m)) {
        stop(paste("invalid variable name",names(pred.means)[k],"specified in argument pred.means"))
      }
    }
  }

#now check that we have a vector variable

  if (is.null(x) | is.factor(data[,x]))  {
    cat("STOP - the number of a vector input variable must be specified...\n")
    return()
  }

# set up storage for the predictors to enable a prediction dataframe to be made
  pred.list <- list(rep(NA,n.preds))

# now set the range of x

  if (is.null(x.range)) {
    x.var <- seq(min(data[,x],na.rm=T),max(data[,x],na.rm=T),length = 100)
  } else {
    x.var <- seq(x.range[1],x.range[2],length = 100)
  }

  pred.list[[1]] <- x.var
  names(pred.list)[1] <- names(data)[x]

#now cycle through and add other variables

  j <- 2
  for (i in 1:n.preds) {
    if (i != x) {
      if (is.vector(data[,i])) {
        m <- match(pred.names[i],names(pred.means))
        if (is.na(m)) {
          pred.list[[j]] <- mean(data[,i],na.rm=T)
        } else pred.list[[j]] <- pred.means[[m]]
      }
      if (is.factor(data[,i])) {
        temp.table <- sort(table(data[,i]),decreasing = TRUE)
        m <- match(pred.names[i],names(pred.means))
        if (is.na(m)) {
          pred.list[[j]] <- names(temp.table)[1]
        } else {
          n <- match(pred.means[[m]],names(temp.table))
          if (is.na(n)){
            stop(paste("invalid level",pred.means[[m]],"provided for ",pred.names[i]))
          }
          pred.list[[j]] <- pred.means[m]
        }
        paste(i,pred.names[i],"\n")
        pred.list[[j]] <- factor(pred.list[[j]],levels=levels(data[,i]))
      }
      names(pred.list)[[j]] <- pred.names[i]
      j <- j + 1
    }
  }

pred.frame <- expand.grid(pred.list)
#
# form the prediction
#
  prediction <- predict.gbm(gbm.object,pred.frame,n.trees = n.trees, type="response")
#
# model smooth if required
#
  if (smooth) {
    pred.glm <- loess(prediction ~ x.var, span = 0.2, data=pred.frame)
#    pred.glm <- glm(prediction ~ ns(AvgDepth, df = 8), data=pred.frame,family=poisson)
    prediction <- fitted(pred.glm)
  }
#
# report the maximum value
#
  max.pred <- max(prediction)
  print(paste("maximum value = ",round(max.pred,2)),quote=FALSE)
#
# mask out values inside hyper-rectangle but outside of sample space
#
  if (mask) {
    mask.trees <- mask.object$gbm.call$best.trees
    point.prob <- predict.gbm(mask.object[[1]],pred.frame, n.trees = mask.trees, type="response")
    prediction[point.prob < 0.5] <- NA
  }
#
# and finally plot the result
#
  if (is.null(y.range)) {
    if (family == "bernoulli") ylim = c(0,1)
    else ylim <- c(min(prediction),max(prediction))
  } else ylim <- y.range

  if (add) {
    lines(x.var, prediction, ...)
  } else {
    plot(x.var, prediction, xlab = x.label, ylab = y.label, type = "l", ylim = ylim,...)
  }

  if (return.data){
    return.data <- as.data.frame(cbind(x.var,prediction))
    names(return.data) <- c(names(data)[x],"richness")
    return(return.data)
  }
}

gbm.perspec <- function(
     gbm.object,
     x = 1,                # the first variable to be plotted
     y = 2,                # the second variable to be plotted
     pred.means = NULL,    # allows specification of values for other variables
     x.label = NULL,       # allows manual specification of the x label
     x.range = NULL,       # manual range specification for the x variable
     y.label = NULL,       # and y label
     y.range = NULL,       # and the y
     z.label = NULL,       # and z label
     z.range = NULL,       # allows control of the vertical axis
     ticktype = "detailed",# specifiy detailed types - otherwise "simple"
     theta = 55,           # rotation
     phi=40,               # and elevation
     smooth = "none",      # controls smoothing of the predicted surface
     mask = FALSE,         # controls masking using a sample intensity model
     perspective = TRUE,   # controls whether a contour or perspective plot is drawn
     ...                   # allows the passing of additional arguments to plotting routine
                           # useful options include shade, ltheta, lphi for controlling illumination
                           # and cex for controlling text size - cex.axis and cex.lab have no effect
)
{
#
# gbm.perspec version 2.9 April 2007
# J Leathwick/J Elith
#
# takes a gbm boosted regression tree object produced by gbm.step and
# plots a perspective plot showing predicted values for two predictors
# as specified by number using x and y
# values for all other variables are set at their mean by default
# but values can be specified by giving a list consisting of the variable name
# and its desired value, e.g., c(name1 = 12.2, name2 = 57.6)

  require(gbm)
  require(splines)

#get the boosting model details

  gbm.call <- gbm.object$gbm.call
  gbm.x <- gbm.call$gbm.x
  n.preds <- length(gbm.x)
  gbm.y <- gbm.call$gbm.y
  pred.names <- gbm.call$predictor.names
  family = gbm.call$family

  x.name <- gbm.call$predictor.names[x]

  if (is.null(x.label)) {
    x.label <- gbm.call$predictor.names[x]
  }

  y.name <- gbm.call$predictor.names[y]

  if (is.null(y.label)) {
    y.label <- gbm.call$predictor.names[y]
  }

  if (is.null(z.label)) {
    z.label <- "fitted value"
  }

  data <- eval(parse(text=gbm.call$dataframe))[,gbm.x]
  n.trees <- gbm.call$best.trees

  if (is.null(x.range)) {
    x.var <- seq(min(data[,x],na.rm=T),max(data[,x],na.rm=T),length = 50)
  } else {x.var <- seq(x.range[1],x.range[2],length = 50)}

  if (is.null(y.range)) {
    y.var <- seq(min(data[,y],na.rm=T),max(data[,y],na.rm=T),length = 50)
  } else {y.var <- seq(y.range[1],y.range[2],length = 50)}

  pred.frame <- expand.grid(list(x.var,y.var))
  names(pred.frame) <- c(x.name,y.name)

  j <- 3
  for (i in 1:n.preds) {
    if (i != x & i != y) {
      if (is.vector(data[,i])) {
        m <- match(pred.names[i],names(pred.means))
        if (is.na(m)) {
          pred.frame[,j] <- mean(data[,i],na.rm=T)
        } else pred.frame[,j] <- pred.means[m]
      }
      if (is.factor(data[,i])) {
        m <- match(pred.names[i],names(pred.means))
        temp.table <- table(data[,i])
        if (is.na(m)) {
          pred.frame[,j] <- rep(names(temp.table)[2],2500)
        } else pred.frame[,j] <- pred.means[m]
        pred.frame[,j] <- factor(pred.frame[,j],levels=names(temp.table))
      }
      names(pred.frame)[j] <- pred.names[i]
      j <- j + 1
    }
  }
#
# form the prediction
#
  prediction <- predict.gbm(gbm.object,pred.frame,n.trees = n.trees, type="response")

# model smooth if required

  if (smooth == "model") {
    pred.glm <- glm(prediction ~ ns(pred.frame[,1], df = 8) * ns(pred.frame[,2], df = 8), data=pred.frame,family=poisson)
    prediction <- fitted(pred.glm)
  }

# report the maximum value and set up realistic ranges for z

  max.pred <- max(prediction)
  cat("maximum value = ",round(max.pred,2),"\n")

  if (is.null(z.range)) {
    if (family == "bernoulli") {
      z.range <- c(0,1)
    } else if (family == "poisson") {
      z.range <- c(0,max.pred * 1.1)
    } else {
      z.min <- min(data[,y],na.rm=T)
      z.max <- max(data[,y],na.rm=T)
      z.delta <- z.max - z.min
      z.range <- c(z.min - (1.1 * z.delta), z.max + (1.1 * z.delta))
    }
  }
# form the matrix

  pred.matrix <- matrix(prediction,ncol=50,nrow=50)

# kernel smooth if required

  if (smooth == "average") {  
    # apply a 3 x 3 smoothing average
    pred.matrix.smooth <- pred.matrix
    for (i in 2:49) {
      for (j in 2:49) {
        pred.matrix.smooth[i,j] <- mean(pred.matrix[c((i-1):(i+1)),c((j-1):(j+1))])
      }
    }
  pred.matrix <- pred.matrix.smooth
  }

# mask out values inside hyper-rectangle but outside of sample space

  if (mask) {
    mask.trees <- mask.object$gbm.call$best.trees
    point.prob <- predict.gbm(mask.object[[1]],pred.frame, n.trees = mask.trees, type="response")
    point.prob <- matrix(point.prob,ncol=50,nrow=50)
    pred.matrix[point.prob < 0.5] <- 0.0
  }
#
# and finally plot the result
#
  if (!perspective) {
    image(x = x.var, y = y.var, z = pred.matrix, zlim = z.range)
  } else {
    persp(x=x.var, y=y.var, z=pred.matrix, zlim= z.range,      # input vars
      xlab = x.label, ylab = y.label, zlab = z.label,          # labels
      theta=theta, phi=phi, r = sqrt(10), d = 3,               # viewing pars
      ticktype = ticktype, mgp = c(4,1,0), ...                 #
    )
  }
}

calibration <- function(obs, preds, family = "binomial")
{
#
# j elith/j leathwick 17th March 2005
# calculates calibration statistics for either binomial or count data
# but the family argument must be specified for the latter
# a conditional test for the latter will catch most failures to specify
# the family
#

if (family == "bernoulli") family <- "binomial"
pred.range <- max(preds) - min(preds)
if(pred.range > 1.2 & family == "binomial") {
  print(paste("range of response variable is ", round(pred.range, 2)), sep = "", quote = F)
  print("check family specification", quote = F)
  return()
}
if(family == "binomial") {
  pred <- preds + 1e-005
  pred[pred >= 1] <- 0.99999
  mod <- glm(obs ~ log((pred)/(1 - (pred))), family = binomial)
  lp <- log((pred)/(1 - (pred)))
  a0b1 <- glm(obs ~ offset(lp) - 1, family = binomial)
  miller1 <- 1 - pchisq(a0b1$deviance - mod$deviance, 2)
  ab1 <- glm(obs ~ offset(lp), family = binomial)
  miller2 <- 1 - pchisq(a0b1$deviance - ab1$deviance, 1)
  miller3 <- 1 - pchisq(ab1$deviance - mod$deviance, 1)
}
if(family == "poisson") {
  mod <- glm(obs ~ log(preds), family = poisson)
  lp <- log(preds)
  a0b1 <- glm(obs ~ offset(lp) - 1, family = poisson)
  miller1 <- 1 - pchisq(a0b1$deviance - mod$deviance, 2)
  ab1 <- glm(obs ~ offset(lp), family = poisson)
  miller2 <- 1 - pchisq(a0b1$deviance - ab1$deviance, 1)
  miller3 <- 1 - pchisq(ab1$deviance - mod$deviance, 1)
}
calibration.result <- c(mod$coef, miller1, miller2, miller3)
names(calibration.result) <- c("intercept", "slope", "testa0b1", "testa0|b1", "testb1|a")
return(calibration.result)
}

roc <- function(obsdat, preddat)
{
# code adapted from Ferrier, Pearce and Watson's code, by J.Elith
#
# see:
# Hanley, J.A. & McNeil, B.J. (1982) The meaning and use of the area
# under a Receiver Operating Characteristic (ROC) curve.
# Radiology, 143, 29-36
#
# Pearce, J. & Ferrier, S. (2000) Evaluating the predictive performance
# of habitat models developed using logistic regression.
# Ecological Modelling, 133, 225-245.
# this is the non-parametric calculation for area under the ROC curve,
# using the fact that a MannWhitney U statistic is closely related to
# the area
#
  if (length(obsdat) != length(preddat)) {
    stop("obs and preds must be equal lengths")  
  }
  n.x <- length(obsdat[obsdat == 0])
  n.y <- length(obsdat[obsdat == 1])
  xy <- c(preddat[obsdat == 0], preddat[obsdat == 1])
  rnk <- rank(xy)
  wilc <- ((n.x * n.y) + ((n.x * (n.x + 1))/2) - sum(rnk[1:n.x]))/(n.x * n.y)
  return(round(wilc, 4))
}

calc.deviance <- function(obs.values, fitted.values, weights = rep(1,length(obs.values)), family="binomial", calc.mean = TRUE)
{
# j. leathwick/j. elith
#
# version 2.1 - 5th Sept 2005
#
# function to calculate deviance given two vectors of raw and fitted values
# requires a family argument which is set to binomial by default
#
#

if (length(obs.values) != length(fitted.values))
   stop("observations and predictions must be of equal length")

y_i <- obs.values

u_i <- fitted.values

if (family == "binomial" | family == "bernoulli") {

   deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
   deviance <- -2 * sum(deviance.contribs * weights)

}

if (family == "poisson" | family == "Poisson") {

    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    deviance <- 2 * sum(deviance.contribs * weights)

}

if (family == "laplace") {
    deviance <- sum(abs(y_i - u_i))
}

if (family == "gaussian") {
    deviance <- sum((y_i - u_i) * (y_i - u_i))
}



if (calc.mean) deviance <- deviance/length(obs.values)

return(deviance)

}

gbm.interactions <- function(
   gbm.object,
   use.weights = FALSE,     # use weights for samples
   mask.object              # a gbm object describing sample intensity
)
{
#
# gbm.interactions version 2.9
#
# j. leathwick, j. elith - May 2007
#
# functions assesses the magnitude of 2nd order interaction effects
# in gbm models fitted with interaction depths greater than 1
# this is achieved by:
#   1. forming predictions on the linear scale for each predictor pair;
#   2. fitting a linear model that relates these predictions to the predictor
#        pair, with the the predictors fitted as factors;
#   3. calculating the mean value of the residuals, the magnitude of which
#        increases with the strength of any interaction effect;
#   4. results are stored in an array;
#   5. finally, the n most important interactions are identified,
#        where n is 25% of the number of interaction pairs;

  require(gbm)

  gbm.call <- gbm.object$gbm.call
  n.trees <- gbm.call$best.trees
  depth <- gbm.call$interaction.depth
  gbm.x <- gbm.call$gbm.x
  n.preds <- length(gbm.x)
  pred.names <- gbm.object$gbm.call$predictor.names
  cross.tab <- matrix(0,ncol=n.preds,nrow=n.preds)
  dimnames(cross.tab) <- list(pred.names,pred.names)

  if (use.weights) mask.trees <- mask.object$gbm.call$best.trees

  cat("gbm.interactions - version 2.9 \n")
  cat("Cross tabulating interactions for gbm model with ",n.preds," predictors","\n",sep="")

  data <- eval(parse(text=gbm.call$dataframe))[,gbm.x]

  for (i in 1:(n.preds - 1)) {
    # step through the predictor set

    if (is.vector(data[,i])) {  
      # create a sequence through the range
      x.var <- seq(min(data[,i],na.rm=T),max(data[,i],na.rm=T),length = 20)
    } else {
      # otherwise set up simple factor variable
      x.var <- factor(names(table(data[,i])),levels = levels(data[,i]))
    }
    x.length <- length(x.var)

    cat(i,"\n")

    for (j in (i+1):n.preds) {
      #create vector or factor data for second variable

      if (is.vector(data[,j])) {
        y.var <- seq(min(data[,j],na.rm=T),max(data[,j],na.rm=T),length = 20)
      } else {
        y.var <- factor(names(table(data[,j])),levels = levels(data[,j]))
      }
      y.length <- length(y.var)

# and now make a temporary data frame

      pred.frame <- expand.grid(list(x.var,y.var))
      names(pred.frame) <- c(pred.names[i],pred.names[j])

      n <- 3 # and add the balance of the variables to it

      for (k in 1:n.preds) {
        if (k != i & k != j) {
          if (is.vector(data[,k])) {
            # either with the mean
            pred.frame[,n] <- mean(data[,k],na.rm=T)
          } else {
            # or the most common factor level
            temp.table <- sort(table(data[,k]),decreasing = TRUE)
            pred.frame[,n] <- rep(names(temp.table)[1],x.length * y.length)
            pred.frame[,n] <- as.factor(pred.frame[,n])
          }
          names(pred.frame)[n] <- pred.names[k]
          n <- n + 1
        }
      }
#
# form the prediction
#
      prediction <- predict.gbm(gbm.object,pred.frame,n.trees = n.trees, type="link")

      if (use.weights) {
        point.prob <- predict.gbm(mask.object[[1]],pred.frame, n.trees = mask.trees, type="response")
        interaction.test.model <- lm(prediction ~ as.factor(pred.frame[,1]) + as.factor(pred.frame[,2]), weights = point.prob)
        
      } else {
        interaction.test.model <- lm(prediction ~ as.factor(pred.frame[,1]) + as.factor(pred.frame[,2]))
      }

      interaction.flag <- round(mean(resid(interaction.test.model)^2) * 1000,2)

      cross.tab[i,j] <- interaction.flag

    }   # end of j loop
  }  # end of i loop

# create an index of the values in descending order

  search.index <- ((n.preds^2) + 1) - rank(cross.tab, ties.method = "first")

  n.important <- max(2,round(0.1 * ((n.preds^2)/2),0))
  var1.names <- rep(" ",n.important)
  var1.index <- rep(0,n.important)
  var2.names <- rep(" ",n.important)
  var2.index <- rep(0,n.important)
  int.size <- rep(0,n.important)

  for (i in 1:n.important) {

    index.match <- match(i,search.index)

    j <- trunc(index.match/n.preds) + 1
    var1.index[i] <- j
    var1.names[i] <- pred.names[j]

    k <- index.match%%n.preds
    if (k > 0) {
      #only do this if k > 0 - otherwise we have all zeros from here on
      var2.index[i] <- k
      var2.names[i] <- pred.names[k]

      int.size[i] <- cross.tab[k,j]
    }

  }

  rank.list <- data.frame(var1.index,var1.names,var2.index,var2.names,int.size)

  return(list(rank.list = rank.list, interactions = cross.tab, gbm.call = gbm.object$gbm.call))
}

bu.gbm <- function(in.list=c("gbm.step","gbm.fixed","gbm.cv","gbm.oob","gbm.holdout","gbm.kfold","fit.gbm","gbm.anz","gbm.bootstrap","gbm.plot","gbm.unipred", "gbm.perspec","calibration","roc","calc.deviance","gbm.interactions","bu.gbm","map.predictions","gbm.plot.fits", "gbm.simplify","gbm.bootstrap","gbm.predict.grids"), out.file="gbm.functions.2.9.R", drive = "h:/gbm/") {
  print(paste("backing up functions to", paste(drive,out.file,sep="")),quote=FALSE)
  dump(in.list,paste(drive,out.file,sep=""))
}

map.predictions <- function(gbm.preds, x = pred.data$col, y = pred.data$row, plot = TRUE, export = FALSE, filename = NULL)
{
  nrows <- 695     #is x axis
  ncols <- 688     #is y axis
  vec.length <- nrows * ncols
  big.vector <- rep(-0.1, vec.length)
  big.vector[((y - 1) * nrows) + x] <- round(gbm.preds,2)

  map.data <- matrix(big.vector, ncol = ncols, nrow = nrows) #fills by columns

  if (export) {
    if (length(filename) == 0) {return("filename required...")}

    write("ncols         695",filename)
    write("nrows         688",filename,append=T)
    write("xllcorner     4420000",filename,append=T)
    write("yllcorner     -5000000",filename,append=T)
    write("cellsize      4000",filename,append=T)
    write("NODATA_value  -0.1",filename,append=T)
    for (i in 1:ncols) {          #first column is first row
      write(map.data[,i],filename,ncol=nrows,append=T)
      if ((i %% 10) == 0) cat(".")
    }

  }

  if (plot) {

   #reverse column order to get map displaying correctly
    map.data <- map.data[,seq(from=ncols, to = 1, by = -1)]

    image(z = map.data, zlim = c(0,35), col = topo.colors(12))
  }
}

gbm.plot.fits <- function(
  gbm.object,
  mask.presence = FALSE,
  use.factor = FALSE,
  plot.layout = c(3,4)          # define the default layout for graphs on the page
)
{
#
# j leathwick, j elith - 7th January 2005
#
# version 2.0 - developed in R 2.0
#
# to plot distribution of fitted values in relation to ydat from mars or other p/a models
# allows masking out of absences to enable focus on sites with high predicted values
# fitted values = those from model; raw.values = original y values
# label = text species name; ydat = predictor dataset
# mask.presence forces function to only plot fitted values for presences
# use.factor forces to use quicker printing box and whisker plot
# file.name routes to a pdf file of this name
#
# SM change 2010: changed dat to dats to avoid conflict with external name.

  max.plots <- plot.layout[1] * plot.layout[2]
  plot.count <- 0

  dats <- gbm.object$gbm.call$dataframe    #get the dataframe name
  dats <- as.data.frame(eval(parse(text=dats)))   #and now the data

  n.cases <- nrow(dats)

  gbm.call <- gbm.object$gbm.call #and the mars call details
  gbm.x <- gbm.call$gbm.x
  gbm.y <- gbm.call$gbm.y
  family <- gbm.call$family

  xdat <- as.data.frame(dats[,gbm.x])
  ydat <- as.data.frame(dats[,gbm.y])

  n.preds <- ncol(xdat)

  fitted.values <- gbm.object$fitted

  pred.names <- names(dats)[gbm.x]
  sp.name <- names(dats)[gbm.y]

  if (mask.presence) {
  mask <- ydat == 1 
  } else {
    mask <- rep(TRUE, length = n.cases)
  }

  robust.max.fit <- approx(ppoints(fitted.values[mask]), sort(fitted.values[mask],na.last=T), 0.99) #find 99%ile value

  for (j in 1:n.preds) {

    if (plot.count == max.plots) {
     plot.count = 0
    }

    if (plot.count == 0) {
     dev.new(width = 11, height = 8)
     par(mfrow = plot.layout)
    }

    plot.count <- plot.count + 1

    if (is.numeric(xdat[mask,j])) {
      wt.mean <- zapsmall(mean((xdat[mask, j] * fitted.values[mask]^5)/mean(fitted.values[mask]^5),na.rm=TRUE),2)
    } else {
      wt.mean <- "na"
    }
    if (use.factor) {
      temp <- factor(cut(xdat[mask, j], breaks = 12))
      if (family == "binomial") {
        plot(temp, fitted.values[mask], xlab = pred.names[j], ylab = "fitted values", ylim = c(0, 1))
      } else {
        plot(temp, fitted.values[mask], xlab = pred.names[j], ylab = "fitted values")
      }
    } else {
      if (family == "binomial") {
        plot(xdat[mask, j], fitted.values[mask], xlab = pred.names[j], ylab = "fitted values", ylim = c(0, 1))
      } else {
        plot(xdat[mask, j], fitted.values[mask], xlab = pred.names[j], ylab = "fitted values")
      }
    }
    abline(h = (0.333 * robust.max.fit$y), lty = 2.)
    if (j == 1) {
      title(paste(sp.name, ", wtm = ", wt.mean))
    } else {
      title(paste("wtm = ", wt.mean))
    }
  }
}

gbm.simplify <- function(
  gbm.object,                 # a gbm object describing sample intensity
  n.folds = 10,               # number of times to repeat the analysis
  n.drops = "auto",           # can be automatic or an integer specifying the number of drops to check
  alpha = 1,                  # controls stopping when n.drops = "auto"
  prev.stratify = TRUE,       # use prevalence stratification in selecting evaluation data
  eval.data = NULL,           # an independent evaluation data set - leave here for now
  plot = TRUE                 # plot results
)
{
# function to simplify a brt model fitted using gbm.step
#
# version 2.9 - J. Leathwick/J. Elith - June 2007
#
# starts with an inital cross-validated model as produced by gbm.step
# and then assesses the potential to remove predictors using k-fold cv
# does this for each fold, removing the lowest contributing predictor,
# and repeating this process for a set number of steps
# after the removal of each predictor, the change in predictive deviance
# is computed relative to that obtained when using all predictors
# it returns a list containing the mean change in deviance and its se
# as a function of the number of variables removed
# having completed the cross validation, it then identifies the sequence
# of variable to remove when using the full data set, testing this
# up to the number of steps used in the cross-validation phase of the analysis
# with results reported to the screen - it then returns
# a table containing the order in which variables are to be removed
# and a list of vectors, each of which specifies the predictor col numbers
# in the original dataframe  - the latter can be used as an argument to gbm.step
# e.g., gbm.step(data = data, gbm.x = simplify.object$pred.list[[4]]...)
# would implement a new analysis with the original predictor set, minus its
# four lowest contributing predictors
#

require(gbm)

# first get the original analysis details..

  gbm.call <- gbm.object$gbm.call
  data <- eval(parse(text=gbm.call$dataframe))
  n.cases <- nrow(data)
  gbm.x <- gbm.call$gbm.x
  gbm.y <- gbm.call$gbm.y
  family <- gbm.call$family
  lr <- gbm.call$learning.rate
  tc <- gbm.call$tree.complexity
  start.preds <- length(gbm.x)
  max.drops <- start.preds - 2
  response.name <- gbm.call$response.name
  predictor.names <- gbm.call$predictor.names
  n.trees <- gbm.call$best.trees
  pred.list <- list(initial = gbm.x)
  weights <- gbm.object$weights

  if (n.drops == "auto") {
    auto.stop <- TRUE
  } else {
    auto.stop <- FALSE
  }

# take a copy of the original data and starting predictors

  orig.data <- data
  orig.gbm.x <- gbm.x

#  if (!is.null(eval.data)) independent.test <- TRUE
#    else independent.test <- FALSE

# extract original performance statistics...

  original.deviance <- round(gbm.object$cv.statistics$deviance.mean,4)
  original.deviance.se <- round(gbm.object$cv.statistics$deviance.se,4)

  cat("gbm.simplify - version 2.9","\n\n")
  cat("simplifying gbm.step model for ",response.name," with ",start.preds," predictors",sep="")
  cat(" and ",n.cases," observations \n",sep="")
  cat("original deviance = ",original.deviance,"(",original.deviance.se,")\n\n",sep="")

# check that n.drops is less than n.preds - 2 and update if required

  if (auto.stop) {
    cat("variable removal will proceed until average change exceeds the original se\n\n")
    n.drops <- 1 
  } else{
    if (n.drops > start.preds - 2) {
      cat("value of n.drops (",n.drops,") is greater than permitted","\n",
        "resetting value to ",start.preds - 2,"\n\n",sep="")
      n.drops <- start.preds - 2
    } else {
      cat("a fixed number of",n.drops,"drops will be tested\n\n")
    }
  }

# set up storage for results

  dev.results <- matrix(0, nrow = n.drops, ncol = n.folds)
  dimnames(dev.results) <- list(paste("drop.",1:n.drops,sep=""),
   paste("rep.",1:n.folds,sep=""))

  drop.count <- matrix(NA, nrow = start.preds, ncol = n.folds)
  dimnames(drop.count) <- list(predictor.names,paste("rep.",1:n.folds,sep=""))

  original.deviances <- rep(0,n.folds)

  model.list <- list(paste("model",c(1:n.folds),sep=""))     # dummy list for the tree models

# create gbm.fixed function call

  gbm.call.string <- paste("try(gbm.fixed(data=train.data,gbm.x=gbm.new.x,gbm.y=gbm.y,",sep="")
  gbm.call.string <- paste(gbm.call.string,"family=family,learning.rate=lr,tree.complexity=tc,",sep="")
  gbm.call.string <- paste(gbm.call.string,"n.trees = ",n.trees,", site.weights = weights.subset,verbose=FALSE))",sep="")

# now set up the fold structure

  if (prev.stratify & family == "bernoulli") {
    presence.mask <- data[,gbm.y] == 1
    absence.mask <- data[,gbm.y] == 0
    n.pres <- sum(presence.mask)
    n.abs <- sum(absence.mask)

    # create a vector of randomised numbers and feed into presences
    selector <- rep(0,n.cases)
    temp <- rep(seq(1, n.folds, by = 1), length = n.pres)
    temp <- temp[order(runif(n.pres, 1, 100))]
    selector[presence.mask] <- temp

    # and then do the same for absences
    temp <- rep(seq(1, n.folds, by = 1), length = n.abs)
    temp <- temp[order(runif(n.abs, 1, 100))]
    selector[absence.mask] <- temp
  } else {  
    
    #otherwise make them random with respect to presence/absence
    selector <- rep(seq(1, n.folds, by = 1), length = n.cases)
    selector <- selector[order(runif(n.cases, 1, 100))]
  }

# now start by creating the intial models for each fold

  cat("creating initial models...\n\n")

  gbm.new.x <- orig.gbm.x

  for (i in 1:n.folds) {

# create the training and prediction folds

    train.data <- orig.data[selector!=i,]
    weights.subset <- weights[selector != i]
    eval.data <- orig.data[selector==i,]

    model.list[[i]] <- eval(parse(text=gbm.call.string))  # create a fixed size object

# now make predictions to the withheld fold

    u_i <- eval.data[,gbm.y]
    y_i <- predict.gbm(model.list[[i]], eval.data, n.trees, "response")

    original.deviances[i] <- round(calc.deviance(u_i,y_i, family = family, calc.mean = TRUE),4)

  } # end of creating initial models

  n.steps <- 1

  while (n.steps <= n.drops & n.steps <= max.drops) {

    cat("dropping predictor",n.steps,"\n")

    for (i in 1:n.folds) {

# get the right data

    train.data <- orig.data[selector!=i,]
    eval.data <- orig.data[selector==i,]
    weights.subset <- weights[selector != i]

# get the current model details

    gbm.x <- model.list[[i]]$gbm.call$gbm.x
    n.preds <- length(gbm.x)
    these.pred.names <- model.list[[i]]$gbm.call$predictor.names
    contributions <- model.list[[i]]$contributions

# get the index number in pred.names of the last variable in the contribution table

    last.variable <- match(as.character(contributions[n.preds,1]),these.pred.names)
    gbm.new.x <- gbm.x[-last.variable]

# and keep a record of what has been dropped

    last.variable <- match(as.character(contributions[n.preds,1]),predictor.names)
    drop.count[last.variable,i] <- n.steps

    model.list[[i]] <- eval(parse(text=gbm.call.string))  # create a fixed size object

    u_i <- eval.data[,gbm.y]
    y_i <- predict.gbm(model.list[[i]],eval.data,n.trees,"response")

    deviance <- round(calc.deviance(u_i,y_i, family = family, calc.mean = TRUE),4)

# calculate difference between intial and new model by subtracting new from old because we want to minimise deviance

    dev.results[n.steps,i] <- round(deviance - original.deviances[i] ,4)

    }

  if (auto.stop) {
    # check to see if delta mean is less than original deviance error estimate

    delta.mean <- mean(dev.results[n.steps,])

    if (delta.mean < (alpha * original.deviance.se)) {
      n.drops <- n.drops + 1
      dev.results <- rbind(dev.results, rep(0,n.folds))
    }
  }
  n.steps <- n.steps + 1
  }

# now label the deviance matrix

  dimnames(dev.results) <- list(paste("drop.",1:n.drops,sep=""),
   paste("rep.",1:n.folds,sep=""))

# calculate mean changes in deviance and their se

  mean.delta <- apply(dev.results,1,mean)
  se.delta <- sqrt(apply(dev.results,1,var))/sqrt(n.folds)

###########################

  if (plot) {
    y.max <- 1.5 * max(mean.delta + se.delta)
    y.min <- 1.5 * min(mean.delta - se.delta)
    plot(seq(0,n.drops),c(0,mean.delta),xlab="variables removed",
      ylab = "change in predictive deviance",type='l',ylim=c(y.min,y.max))
    lines(seq(0,n.drops),c(0,mean.delta) + c(0,se.delta),lty = 2)
    lines(seq(0,n.drops),c(0,mean.delta) - c(0,se.delta),lty = 2)
    abline(h = 0 , lty = 2, col = 3)
    min.y <- min(c(0,mean.delta))
    min.pos <- match(min.y,c(0,mean.delta)) - 1 # subtract one because now zero base
    abline(v = min.pos, lty = 3, col = 2)
    abline(h = original.deviance.se, lty = 2, col = 2)
    title(paste("RFE deviance - ",response.name," - folds = ",n.folds,sep=""))
  }

# and do a final backwards drop sequence from the original model

  cat("\nnow processing final dropping of variables with full data \n\n")

  gbm.call.string <- paste("try(gbm.fixed(data=orig.data,gbm.x=gbm.new.x,gbm.y=gbm.y,",sep="")
  gbm.call.string <- paste(gbm.call.string,"family=family,learning.rate=lr,tree.complexity=tc,",sep="")
  gbm.call.string <- paste(gbm.call.string,"n.trees = ",n.trees,", site.weights = weights,verbose=FALSE))",sep="")

  n.steps <- n.steps - 1 #decrement by one to reverse last increment in prev loop

  final.model <- gbm.object  # restore the original model and data
  train.data <- orig.data

# and set up storage

  final.drops <- matrix(NA, nrow = start.preds, ncol = 1)
  dimnames(final.drops) <- list(predictor.names,"step")

  for (i in 1:n.steps) {

# get the current model details

    gbm.x <- final.model$gbm.call$gbm.x
    n.preds <- length(gbm.x)
    these.pred.names <- final.model$gbm.call$predictor.names
    contributions <- final.model$contributions

    cat(i,"-",as.character(contributions[n.preds,1]),"\n")

# get the index number in pred.names of the last variable in the contribution table

    last.variable <- match(as.character(contributions[n.preds,1]),these.pred.names)
    gbm.new.x <- gbm.x[-last.variable]

# and keep a record of what has been dropped

    last.variable <- match(as.character(contributions[n.preds,1]),predictor.names)
    final.drops[last.variable] <- i

    final.model <- eval(parse(text=gbm.call.string))  # create a fixed size object

  }

#and then the corresponding numbers

  removal.list <- dimnames(final.drops)[[1]]
  removal.list <- removal.list[order(final.drops)]
  removal.list <- removal.list[1:n.drops]

  removal.numbers <- rep(0,n.steps)

# construct predictor lists to faciliate final model fitting

  for (i in 1:n.steps) {
    removal.numbers[i] <- match(removal.list[i],predictor.names)
    pred.list[[i]] <- orig.gbm.x[0-removal.numbers[1:i]]
    names(pred.list)[i] <- paste("preds.",i,sep="")
  }

  deviance.summary <- data.frame(mean = round(mean.delta,4), se = round(se.delta,4))

  final.drops <- data.frame("preds" = dimnames(final.drops)[[1]][order(final.drops)],
     "order" = final.drops[order(final.drops)])

  return(list(deviance.summary = deviance.summary, deviance.matrix = dev.results, drop.count = drop.count, final.drops = final.drops, pred.list = pred.list, gbm.call = gbm.call))

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
  train.data <- eval(parse(text=gbm.call$dataframe))
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
  if (bootstrap.predictions) cat("prediction dataset has ",n.pred.obs," rows\n\n",sep="")
  else cat("no prediction dataset provided...\n\n")

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
      train.fit.upper[i] <- approx(ppoints(temp),sort(temp),0.975)$y
      train.fit.lower[i] <- approx(ppoints(temp),sort(temp),0.025)$y
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
      train.pred.upper[i] <- approx(ppoints(temp),sort(temp),0.975)$y
      train.pred.lower[i] <- approx(ppoints(temp),sort(temp),0.025)$y
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


gbm.predict.grids <- function(
       model,
       new.dat,
       export.grids = F,
       preds2R = T,
       sp.name = "preds",
       filepath = NULL,
       num.col = NULL,
       num.row = NULL,
       row.index = NULL,
       col.index = NULL,
       xll = NULL,
       yll = NULL,
       cell.size = NULL,
       no.data = NULL,
       plot=F,
       full.grid=T,
       header = T,
       north.limit = NULL,
       south.limit = NULL
)
{
# J.Elith / J.Leathwick, March 07
# to make predictions to sites or grids. If to sites, the
# predictions are written to the R workspace. If to grid,
# the grids are written to a nominated directory and optionally also
# plotted in R
#
# new data (new.dat) must be a data frame with column names identical
# to names for all variables in the model used for prediction
#
# pred.vec is a vector of -9999's, the length of the scanned full grid
# (i.e. without nodata values excluded)
#
# filepath must specify the whole path as a character vector,but without the final file
# name - eg "c:/gbm/"

  require(gbm)

  pred.vec <- rep(no.data,num.col * num.row)

  cat("\nForming predictions for gridded data...\n")

  temp <- predict.gbm(model, new.dat, n.trees=model$gbm.call$best.trees, type="response")

  if (!is.null(north.limit)) {
     temp[new.dat$y < north.limit] <- no.data
  }

  if (!is.null(south.limit)) {
     temp[new.dat$y > south.limit] <- no.data
  }

  if (is.null(row.index)) {
    pred.vec <- temp
  } else {
    cell.pointer <- ((row.index - 1) * num.col) + (col.index)
    pred.vec[cell.pointer] <- temp
  }

  cat("\n and exporting the results...\n")
  if(export.grids) {
    newname <- paste(filepath, sp.name,".asc", sep="")

    #full.pred <- pred.vec
    #full.pred[as.numeric(row.names(new.dat))] <- temp

    if(header){
      write(paste("ncols          ",num.col,sep=""),newname)
      write(paste("nrows          ",num.row,sep=""),newname,append=T)
      write(paste("xllcorner      ",xll,sep=""),newname,append=T)
      write(paste("yllcorner      ",yll,sep=""),newname,append=T)
      write(paste("cellsize       ",cell.size,sep=""),newname,append=T)
      write(paste("NODATA_value ",no.data,sep=""),newname,append=T)
    }

    full.pred.mat <- matrix(pred.vec, nrow=num.row, ncol=num.col, byrow=T)

    write.table(full.pred.mat, newname, sep=" ", append=T, row.names=F, col.names=F)

    if (plot) {
      image(z = t(full.pred.mat)[, nrow(full.pred.mat):1], zlim =  c(0,1), col = rev(topo.colors(12)))
    }
  }

#also write to R directory, if required:

  if(preds2R){
    assign(sp.name,temp, pos=1)
  }
}



#######################################################
# SM function to plot bootstrapped info


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
      dev.new(width = 11, height = 8)
      par(mfrow = plot.layout)
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
