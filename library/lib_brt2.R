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

library("dismo")
data(Anguilla_train)
d <- Anguilla_train
d <- d[,-1]
d <- d[1:200,]
# d <- d[,-which(names(d) %in% "Method")]
data <- d
y <- d$Angaus
x <- d[,-which(names(d)=="Angaus")]


brt.fit <- function(x, y, n.trees=NULL, min.n.trees=3000, n.boot=0, verbose=FALSE, ...) {
   #
   # Fit a Boosted Regression Tree model with optimization of number of trees and bootstraps
   #
   # x            matrix or data.frame of predictors
   # y            vector of responses
   # n.trees      number of trees
   #              if NULL, the number of trees and shrinkage are estimated
   # min.n.trees  minimum number of trees to be reached in the optimization of shrinkage and number of trees
   # n.boot       number of bootstraps (should be at least 100)
   #
   
   library("gbm")
   
   # sanity checks
   if ( ! is.vector(y) ) {
      stop("Need only one response variable, y should be a vector")
   }
   
   if (nrow(x) != length(y)) {
      stop("Predictors and response variable must be the same length")
   }

   if (is.matrix(x)) {
      x <- data.frame(x)
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
   
   # optimise shrinkage and number of trees if the number of trees is not forced
   if ( is.null(n.trees) ) {
      if ( verbose ) message("Optimising shrinkage and number of trees")

      # start with some defaults
      shrinkage <- 0.02    # NB: will be reduced before the first model
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
   if ( verbose ) message("Cross-validating final model (", cv.fold, " folds)")
   m <- gbm(response ~ . , data=d, shrinkage=shrinkage, n.trees=n.trees, cv.fold=cv.fold, ...)
   best.iter <- gbm.perf(m, method="cv", plot.it=verbose)
   if (verbose) message("  shrinkage : ", round(shrinkage, 4), " | trees : ", best.iter, " / ", m$n.trees)

   return(m)
}

m <- brt.fit(x, y, distribution="bernoulli", interaction.depth=4, verbose=TRUE)






