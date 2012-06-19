#
#      Test polygon clipping functions
#
#  (c) Copyright 2012 Jean-Olivier Irisson
#      GNU General Public License v3
#
#-----------------------------------------------------------------------------

suppressPackageStartupMessages(require("testthat", quietly=TRUE))

context("Polygon clipping")

# generate data
set.seed(123)
x <- data.frame(lon=runif(100, min=-180, max=180), lat=runif(100, min=-90, max=90))

# set some limits
lon.min <- -90
lon.max <- 20
lat.min <- -70
lat.max <- -30

# cut polygon
xC <- clip.polygon(x, lon.min, lon.max, lat.min, lat.max)

test_that("only data in the appropriate range is kept", {
  # compute range of result
  lonR <- range(xC$lon, na.rm=T) 
  latR <- range(xC$lat, na.rm=T)
  expect_that(
    c(lonR, latR), equals(c(lon.min, lon.max, lat.min, lat.max))
  )
})

test_that("order of lon and lat does not matter", {
  # clip polygon inverting lat and lon
  xC2 <- clip.polygon(x[,c(2,1)], lon.min, lon.max, lat.min, lat.max)
  expect_that(
    xC2, equals(xC)
  )
})

# prepare a fake data data.frame
data <- data.frame(lon=c(lon.min, lon.min, lon.max, lon.max), lat=c(lat.min, lat.max, lat.max, lat.min))

test_that("clipping to data is equivalent to clipping to limits", {
  xCd <- clip.to.data(x, data, expand=0)
  expect_that(
    xCd, equals(xC)
  )
})

test_that("expand works", {
  expand <- 3.75
  xCd <- clip.to.data(x, data, expand=expand)
  lonR <- range(xCd$lon, na.rm=T) 
  latR <- range(xCd$lat, na.rm=T)
  expect_that(
    c(lonR, latR), equals(c(lon.min-expand, lon.max+expand, lat.min-expand, lat.max+expand))
  )
})
