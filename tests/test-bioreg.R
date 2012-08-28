#
#      Test bioregionalisation functions
#
#  (c) Copyright 2012 Jean-Olivier Irisson
#      http://creativecommons.org/licenses/by/3.0/
#
#-----------------------------------------------------------------------------

suppressPackageStartupMessages(require("testthat", quietly=TRUE))
suppressPackageStartupMessages(require("stringr", quietly=TRUE))

# make sure that the path to the env_data folder works when it is relative
path <- getOption("atlasr.env.data")
if (substr(path, 1, 1) != "/") {
 path <- paste("../", path, sep="")
}

context("Bioregionalisation")

# # prepare fake data
# set.seed(123)
# n <- 20
# lons <- seq(30, 50, n)
# lats <- seq(-50, -30, n)
# data <- 
# n <- 100
# data <- data.frame(lon=runif(n, -180, 180), lat=runif(n, -80, -30), sp1=round(runif(100)), sp2=round(runif(100)))
# temp <- tempfile(fileext=".csv")
# write.table(data, file=temp, row.names=FALSE, sep=",")
# 

temp <- tempdir()

# run bioregionalisation
capture.output(suppressMessages(
  b <- bioreg(c("bathymetry", "sst_summer_climatology"), n.groups=2, n.groups.intermediate=5, lat.min=-40, lat.max=-30, lat.step=2, lon.min=0, lon.max=20, lon.step=2, output=temp, path=path)
), file=tempfile())


test_that("clara clustering is stable", {
  expect_that(
    as.numeric(b$data$clara.num),
    equals(c(1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 4, 5, 3, 3, 3, 3, 4, 4, 5, 5, 4, 4, 3, 3, 3, 3, 4, 5, 5))
  )
})

test_that("hclust clustering is stable", {
  expect_that(
    as.numeric(b$data$cluster),
    equals(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2))
  )
})


# test_that("missing weights induce warning", {
#   expect_that(
#     bioreg(c("bathymetry", "sst_summer_climatology"), weights=c(bathy=1), n.groups=2, n.groups.intermediate=5, lat.min=-40, lat.max=-30, lat.step=2, lon.min=0, lon.max=20, lon.step=2, output=temp),
#     gives_warning()
#   )
# })
# 
# test_that("supplementary weights induce error", {
#   expect_that(
#     bioreg(c("bathymetry", "sst_summer_climatology"), weights=c(nox=1), n.groups=2, n.groups.intermediate=5, lat.min=-40, lat.max=-30, lat.step=2, lon.min=0, lon.max=20, lon.step=2, output=temp),
#     throws_error("No variable matching")
#   )
# })
# 
# test_that("missing transformations induce warning", {
#   expect_that(
#     bioreg(c("bathymetry", "sst_summer_climatology"), transformations=c(bathy="x"), n.groups=2, n.groups.intermediate=5, lat.min=-40, lat.max=-30, lat.step=2, lon.min=0, lon.max=20, lon.step=2, output=temp),
#     gives_warning("Transformations are missing")
#   )
# })
# 
# test_that("supplementary transformations induce error", {
#   expect_that(
#     bioreg(c("bathymetry", "sst_summer_climatology"), transformations=c(nox="x"), n.groups=2, n.groups.intermediate=5, lat.min=-40, lat.max=-30, lat.step=2, lon.min=0, lon.max=20, lon.step=2, output=temp),
#     throws_error("No variable matching")
#   )
# })


unlink(temp)