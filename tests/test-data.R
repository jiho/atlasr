suppressPackageStartupMessages(require("testthat", quietly=TRUE))


context("Access to environment data")

test_that("error when matching non existent variables", {
  expect_that(
    match.vars(vars="foo", choices=c("bar", "bob"), quiet=TRUE),
    throws_error("No variable matching")
  )
})

test_that("variable name expansion works", {
  expect_that( match.vars("foo", "foobar", quiet=TRUE), matches("foobar") )
})

test_that("variable name expansion issues a message when asked for", {
  expect_that( match.vars("foo", "foobar", quiet=FALSE), shows_message("foobar") )
})

test_that("listing environment returns a non-empty character vector", {
  x <- list.env.data()

  expect_that( x, is_a("character") )
  expect_that( length(x) > 1, is_true() )
})

test_that("masking land works", {
  x <- read.env.data("bathymetry")
  x <- mask.env.data(x)

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
})


context("Data processing")

test_that("writing shapefile to a directory with accents in the name works", {
  temp <- tempfile("test_accent_éü")

  d <- data.frame(lon=180, lat=-30, x=1)

  write.shapefile(d, temp, variables="x")

  expect_that( file.exists(paste(temp, ".shp", sep="")), is_true() )
  expect_that( file.exists(paste(temp, ".shx", sep="")), is_true() )
  expect_that( file.exists(paste(temp, ".dbf", sep="")), is_true() )
})

test_that("writing shapefile to a directory with spaces in the name works", {
  temp <- tempfile("test_space_ _")

  d <- data.frame(lon=180, lat=-30, x=1)

  write.shapefile(d, temp, variables="x")

  expect_that( file.exists(paste(temp, ".shp", sep="")), is_true() )
  expect_that( file.exists(paste(temp, ".shx", sep="")), is_true() )
  expect_that( file.exists(paste(temp, ".dbf", sep="")), is_true() )
})


