#
#      Test access to data
#
#  (c) Copyright 2012 Jean-Olivier Irisson
#      http://creativecommons.org/licenses/by/3.0/
#
#-----------------------------------------------------------------------------

suppressPackageStartupMessages(require("testthat", quietly=TRUE))


context("Access to environment data")

test_that("error when matching non existent variables", {
  expect_that(
    partial.match(vars="foo", choices=c("bar", "bob"), quiet=TRUE),
    throws_error("No variable matching")
  )
})

test_that("variable name expansion works", {
  expect_that( partial.match("foo", "foobar", quiet=TRUE), matches("foobar") )
})

test_that("variable name expansion issues a message when asked for", {
  expect_that( partial.match("foo", "foobar", quiet=FALSE), shows_message("foobar") )
})

# make sure that the path to the env_data folder works when it is relative
path <- getOption("atlasr.env.data")
if (substr(path, 1, 1) != "/") {
 path <- paste("../", path, sep="")
}

test_that("listing environment returns a non-empty character vector", {
  x <- list.env.data(path=path)

  expect_that( x, is_a("character") )
  expect_that( length(x) > 1, is_true() )
})

test_that("masking land works", {
  x <- read.env.data("bathymetry", path=path)
  x <- mask.env.data(x, path=path)

  expect_that( all(na.omit(x[[1]]$z) < 0), is_true() )
})


context("Access to species data files")

test_that("latitude and longitude names are flexible", {
  temp <- tempfile(fileext=".csv")
  write.csv(data.frame(latitude=1, longitude=1), temp)
  d <- read.data(temp)
  expect_that( all(c("lat", "lon") %in% names(d)), is_true() )

  write.csv(data.frame(Latitude=1, LonGItude=1), temp)
  d <- read.data(temp)
  expect_that( all(c("lat", "lon") %in% names(d)), is_true() )

  write.csv(data.frame(lat=1, long=1), temp)
  d <- read.data(temp)
  expect_that( all(c("lat", "lon") %in% names(d)), is_true() )

  write.csv(data.frame(LAT=1, LONG=1), temp)
  d <- read.data(temp)
  expect_that( all(c("lat", "lon") %in% names(d)), is_true() )

  unlink(temp)
})


context("Writing shapefiles")

test_that("writing shapefile to a directory with accents in the name works", {
  temp <- tempfile("test_accent_éü")

  d <- data.frame(lon=180, lat=-30, x=1)

  write.shapefile(d, temp, variables="x")

  expect_that( file.exists(paste(temp, ".shp", sep="")), is_true() )
  expect_that( file.exists(paste(temp, ".shx", sep="")), is_true() )
  expect_that( file.exists(paste(temp, ".dbf", sep="")), is_true() )

  unlink(temp)
})

test_that("writing shapefile to a directory with spaces in the name works", {
  temp <- tempfile("test_space_ _")

  d <- data.frame(lon=180, lat=-30, x=1)

  write.shapefile(d, temp, variables="x")

  expect_that( file.exists(paste(temp, ".shp", sep="")), is_true() )
  expect_that( file.exists(paste(temp, ".shx", sep="")), is_true() )
  expect_that( file.exists(paste(temp, ".dbf", sep="")), is_true() )

  unlink(temp)
})


context("Weighting data")

test_that("weighting data warns or errors appropriately when names do not match", {
  x <- data.frame(foo=1:3, bar=4:6)
  # additional weight
  expect_that(weight.data(x, c(foo=1, bar=2, bob=3)), throws_error())
  # missing weight and additional weight
  expect_that(weight.data(x, c(foo=1, bob=2)), throws_error())
  # missing weight
  expect_that(weight.data(x, c(foo=1)), gives_warning("Weights are missing"))
})

test_that("weighting data is not influenced by column order", {
  data <- data.frame(foo=1:3, bar=4:6)
  # precompute what weighting should give
  dataW <- data
  dataW$bar <- dataW$bar * 0.5
  expect_equivalent(weight.data(data, c(foo=1, bar=0.5)), dataW)
  expect_equivalent(weight.data(data, c(bar=0.5, foo=1)), dataW)
  expect_equivalent(weight.data(data[,c(2,1)], c(foo=1, bar=0.5)), dataW[,c(2,1)])
  # NB: weight.data sets some attributes so we use "equivalent" instead of "equal"
})

test_that("omitting weights defaults to weight=1", {
  data <- data.frame(foo=1:3, bar=4:6)
  # all weights are 1
  expect_that(weight.data(data, weights=c(foo=1)  , warn=FALSE)$bar, equals(data$bar))
  # after scaling, weights=c(2,1) become c(1, 0.5)
  expect_that(weight.data(data, weights=c(foo=2)  , warn=FALSE)$bar, equals(data$bar*0.5))
  # since weight of foo is 0.5, weight of bar=1 becomes the maximum and is not affected by scaling
  expect_that(weight.data(data, weights=c(foo=0.5), warn=FALSE)$bar, equals(data$bar))
})

test_that("name expansion cause weights to be replicated", {
  data <- data.frame(foo=1:3, fii=4:6, bar=4:6)
  dataW <- weight.data(data, weights=c(f=1, bar=2))
  # weight of 1 should be replicated and becomes 0.5 after scaling by maximum weight
  expect_that(dataW$foo, equals(data$foo * 0.5))
  expect_that(dataW$fii, equals(data$fii * 0.5))
})

test_that("weighting saves weights as attributes", {
  x <- data.frame(foo=1:3, bar=4:6)

  xW <- weight.data(x, c(foo=2, bar=1))
  w <- attr(xW, "weights")
  expect_that(names(w), equals(names(x)))

  # no matter the order
  xW <- weight.data(x, c(bar=1, foo=2))
  w2 <- attr(xW, "weights")
  expect_that(w2, equals(w))

  # or whether some are ommitted
  xW <- weight.data(x, c(foo=2), warn=F)
  w2 <- attr(xW, "weights")
  expect_that(w2, equals(w))
})


context("Transforming data")

test_that("transforming data warns or errors appropriately when names do not match", {
  x <- data.frame(foo=1:3, bar=4:6)
  # additional transformation
  expect_that(transform.data(x, c(foo="x", bar="x", bob="x")), throws_error())
  # missing transformation and additional transformation
  expect_that(transform.data(x, c(foo="x", bob="x")), throws_error())
  # missing transformation
  expect_that(transform.data(x, c(foo="x")), gives_warning("Transformations are missing"))
})

test_that("transforming data is not influenced by column order", {
  data <- data.frame(foo=1:3, bar=4:6)
  # precompute what transforming should give
  dataW <- data
  dataW$foo <- sqrt(dataW$foo)
  dataW$bar <- log(dataW$bar)
  expect_equivalent(transform.data(data, c(foo="sqrt(x)", bar="log(x)")), dataW)
  expect_equivalent(transform.data(data, c(bar="log(x)", foo="sqrt(x)")), dataW)
  expect_equivalent(transform.data(data[,c(2,1)], c(foo="sqrt(x)", bar="log(x)")), dataW[,c(2,1)])
})

test_that("omitting transformation defaults to no transformation", {
  data <- data.frame(foo=1:3, bar=4:6)
  expect_that(transform.data(data, c(foo="log(x)"), warn=FALSE)$bar, equals(data$bar))
})

test_that("name expansion cause transformations to be replicated", {
  data <- data.frame(foo=1:3, fii=4:6, bar=4:6)
  dataW <- transform.data(data, c(f="log(x)", bar="x"))
  # log transformation is replicated
  expect_that(dataW$foo, equals(log(data$foo)))
  expect_that(dataW$fii, equals(log(data$fii)))
})

test_that("transforming saves transforms as attributes", {
  x <- data.frame(foo=1:3, bar=4:6)

  xW <- transform.data(x, c(foo="log(x)", bar="x"))
  w <- attr(xW, "transformations")
  expect_that(names(w), equals(names(x)))

  # no matter the order
  xW <- transform.data(x, c(bar="x", foo="log(x)"))
  w2 <- attr(xW, "transformations")
  expect_that(w2, equals(w))

  # or whether some are ommitted
  xW <- transform.data(x, c(foo="log(x)"), warn=F)
  w2 <- attr(xW, "transformations")
  expect_that(w2, equals(w))
})

test_that("defining a function from a string accepts various inputs", {
  # explicit function definitions
  expect_that(is.function(function.maker("log(x)")), is_true())
  # pre-existing functions
  expect_that(is.function(function.maker("log")), is_true())
  # custom defined functions
  myFunc <- function(x) {x}
  expect_that(is.function(function.maker("myFunc")), is_true())
})

test_that("unparsable functions produce an error", {
  expect_that( function.maker("$"),  throws_error() )
  expect_that( transform.data(data.frame(foo=1:2), c(foo="$")),  throws_error() )
})

test_that("errors and warnings in transformation functions are handled", {
  # TODO split that into several tests
  data <- data.frame(foo=-1:1)

  # errors are transformed into warnings
  expect_that(transform.data(data, c(foo="stop('my error')")), gives_warning("Error when applying transformation"))
  # and return NA
  res <- suppressWarnings(transform.data(data, c(foo="stop('my error')"))$foo)
  expect_that(all(is.na(res)), is_true())

  # warnings (here : produce NaNs) are issued in context
  expect_that(transform.data(data, c(foo="log(x)")), gives_warning("Warning when applying transformation"))
  # but still computes the result
  res <- suppressWarnings(transform.data(data, c(foo="log(x)"))$foo)
  expect_that(res, equals(suppressWarnings(log(data$foo))))

  # applying a meaningless "function" throws an error which is converted into a warning
  expect_that(transform.data(data, c(foo="aMeaningLessVariable39825")), gives_warning("Error when applying transformation"))
})
