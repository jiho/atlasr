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


brt.fit <- function(x, y, n.trees=NULL, shrinkage=0.01, min.n.trees=3000, n.boot=0, verbose=FALSE, ...) {
   #
   # Fit a Boosted Regression Tree model with optimization of number of trees and bootstraps
   #
   # x            matrix or data.frame of predictors
   # y            vector of responses
   # n.trees      number of trees
   #              if NULL, the number of trees and shrinkage are estimated
   # shrinkage    shrinkage value for fixed number of trees or starting value for
   #              the optimization of shrinkage and number of trees
   # min.n.trees  minimum number of trees to be reached in the optimization of
   #              shrinkage and number of trees
   # n.boot       number of bootstraps (should be at least 100)
   # continuous.resolution number of equally space points at which to
   #              evaluate continuous predictors

   suppressPackageStartupMessages(library("gbm", quietly=TRUE))
   suppressPackageStartupMessages(library("plyr", quietly=TRUE))

   # sanity checks
   #--------------------------------------------------------------------------
   if ( ! is.vector(y) ) {
      stop("Need only one response variable, y should be a vector")
   }

   if (nrow(x) != length(y)) {
      stop("Predictors and response variable must be the same length")
   }

   if (is.matrix(x)) {
      x <- data.frame(x)
   }

   if ( n.boot != 0 & n.boot < 100 ) {
      warning("Less than 100 bootstraps is useless. Forcing n.boot=100")
      n.boot <- 100
   }

   # compute good CV fold (at least n.per.fold observations per fold)
   n.per.fold <- 30
   cv.fold <- length(y) %/% n.per.fold
   if (cv.fold < 2) {
      stop("Not enough data to correctly fit a model. Need at least ", n.per.fold * 2, " data points")
   }
   cv.fold <- min(cv.fold, 10)

   # add response variable to the data.frame
   d <- x
   d$response <- y


   # fit model
   #--------------------------------------------------------------------------

   # optimise shrinkage and number of trees if the number of trees is not forced
   if ( is.null(n.trees) ) {
      if ( verbose ) message("Optimising shrinkage and number of trees")

      # start with some defaults
      shrinkage <- shrinkage * 1/0.75   # NB: will be reduced before the first model
      best.iter <- 0
      step.n.trees <- 100

      while ( best.iter < min.n.trees ) {

         # if the number of trees is too low, try reducing the shrinkage
         shrinkage <- shrinkage * 0.75

         # prime the model
         m <- gbm(response ~ . , data=d, shrinkage=shrinkage, n.trees=step.n.trees, train.fraction=0.9, ...)
         best.iter <- suppressWarnings(gbm.perf(m, method="OOB", plot.it=FALSE))
         # NB: at this point best.iter is probably the maximum number of trees in the model

         # increase number of trees until the best number of trees (best.iter) is well below the total number of trees computed. This means enough trees will have been computed.
         # The cross-validation estimation of the optimal number of trees is the best but is quite long. The estimations using OOB or test methods are not as good as the one using cross validation, especially for small shrinkage.
         # => we take an additional buffer in the number of trees, inversely proportional to shrinkage
         while ( m$n.trees - best.iter < 2 / shrinkage ) {
            # add some trees
            m <- gbm.more(m, n.new.trees=step.n.trees)

            # OOB underestimates
            # test over estimates
            # => we take the average
            best.iter.OOB <- suppressWarnings(gbm.perf(m, method="OOB", plot.it=FALSE))
            best.iter.test <- suppressWarnings(gbm.perf(m, method="test", plot.it=FALSE))
            best.iter <- round((best.iter.OOB + best.iter.test) / 2)

            if (verbose) message("  shrinkage : ", round(shrinkage, 4), " | trees : ", best.iter, " / ", m$n.trees)
         }
      }
      n.trees <- m$n.trees
   }

   # fit the final model with cross validation and get a better, final, estimate of the actual optimal number of trees
   if ( verbose ) message("Cross-validating model (", cv.fold, " folds)")
   m <- gbm(response ~ . , data=d, shrinkage=shrinkage, n.trees=n.trees, cv.fold=cv.fold, ...)
   m$best.iter <- gbm.perf(m, method="cv", plot.it=verbose)
   if (verbose) message("  shrinkage : ", round(shrinkage, 4), " | trees : ", m$best.iter, " / ", m$n.trees)


   # bootstraps
   #--------------------------------------------------------------------------

   if (n.boot > 0) {

      if ( verbose ) message("Bootstrapping model (", n.boot, " bootstraps)")

      boot.brt <- function(i, d, n.trees, shrinkage, ...) {
         # i            bootstrap number
         # d            original data
         # n.trees      number of trees to use: optimal number of trees in the original model
         # shrinkage    shrinkage of the original model

         # resampe original data with replacement
         dd <- d[sample.int(nrow(d), replace=TRUE),]

         # fit model on this resampled dataset
         library("gbm")
         mb <- gbm(response ~ . , data=d, shrinkage=shrinkage, n.trees=n.trees, keep.data=FALSE, ...)

         return(mb)
      }

      # parallel computation takes too much memory, probably because of a problem with gbm
      # disable it for now because it takes almost as much time and way more memory
      parallel <- FALSE
      if ( parallel ) {
         # try to parallelise bootstraps
         suppressPackageStartupMessages(library("parallel", quietly=TRUE))
         n <- detectCores()

         # parallel
         cl <- makeCluster(n)
         if (verbose) message("  parallel computation on ", n, " cores")
         m$boot <- clusterApply(cl=cl, x=1:n.boot, fun=boot.brt, d=d, n.trees=m$best.iter, shrinkage=shrinkage, ...)
         stopCluster(cl)

         # # foreach (a bit worse regarding memory management)
         # library("doParallel")
         # registerDoParallel(cores=n)
         # m$boot <- alply(1:n.boot, 1, .fun=boot.brt, d=d, n.trees=m$best.iter, shrinkage=shrinkage, .progress="text", .parallel=T, ...)

      } else {
         m$boot <- alply(1:n.boot, 1, .fun=boot.brt, d=d, n.trees=m$best.iter, shrinkage=shrinkage, .progress="text", ...)
      }
   }

   # create a specific class
   class(m) <- c("brt", class(m))

   return(m)
}


AUC <- function(y, p) {
   #
   # Code derived from the .roc function by Elith
   # This is the non-parametric calculation for area under the ROC curve, using the fact that a MannWhitney U statistic is closely related to the area
   #
   # y   observed data (pres=1, abs=0)
   # p   predicted probabilities
   #

   if ( length(y) != length(p) ) {
      stop("y and p must be the same length")
   }
   n.0 <- sum(y==0)
   n.1 <- sum(y==1)
   pOrder <- c(p[y==0], p[y==1])
   pRank <- rank(pOrder)

   U <- ( n.0 * n.1  +  n.0 * (n.0 + 1) / 2  -  sum(pRank[1:n.0]) ) / ( n.0 * n.1 )

   return(U)
}

ilogit <- function(x) {
   # Inverse logit function
   1 / (1 + exp(-x))
}

summary.brt <- function(m, n.trees=m$best.iter, ...) {

   if (m$distribution$name != "bernoulli") {
      stop("Summary statistics only implemented for Bernoulli distribution")
   }

   # percentage of deviance explained
   y <- m$data$y
   w <- m$data$w
   # null hypothesis proba and associated null deviance
   p <- sum(y * w) / sum(w)
   logLikelihood <- sum( (y * log(p) + (1-y) * log(1 - p)) * w )
   nullDeviance <- -2 * logLikelihood
   nullDeviance <- nullDeviance / length(y)
   # explained deviance (remaining deviance is m$train.error)
   explainedDevianceTrain <- 1 - m$train.error[n.trees] / nullDeviance
   explainedDevianceCV    <- 1 - m$cv.error[n.trees]    / nullDeviance

   # AUC
   AUC <- AUC(y, ilogit(m$cv.fitted))

   # performance
   if (AUC < 0.7) {
      performance <- "poor :-("
   } else if ( AUC < 0.8) {
      performance <- "acceptable :-/"
   } else if ( AUC < 0.9) {
      performance <- "good :-)"
   } else if ( AUC > 0.9) {
      performance <- "great :-D !"
   }

   # relative influence of variables
   rel.inf <- relative.influence(m, m$best.iter)
   rel.inf <- 100 * rel.inf/sum(rel.inf)
   rel.inf <- sort(rel.inf, decreasing=TRUE)

   cat("A gradient boosted model with bernoulli loss function.\n")
   cat(m$n.trees,"iterations were performed\n")
   cat("The best cross-validated iteration was",m$best.iter,"\n")
   cat("Predictors contributions in % of total influence :\n")
   print(rel.inf)
   cat("Explained deviance (on training set) =", round(explainedDevianceTrain * 100, 1), "%\n")
   cat("Cross-validated explained deviance =", round(explainedDevianceCV * 100, 1), "%, AUC =", round(AUC, 2), "\n")
   cat("Predictive performance of the model is", performance, "\n")

   return(invisible(rel.inf))
}

summary(m)



