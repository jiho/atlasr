#
#      Test bioregionalisation functions
#
#  (c) Copyright 2012 Jean-Olivier Irisson
#      GNU General Public License v3
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
  b <- bioreg(c("bathymetry", "bathymetry_slope"), n.groups=5, n.groups.intermediate=100, lat.min=-50, lat.max=-30, lat.step=2, lon.min=0, lon.max=90, lon.step=4, output=temp, path=path)
), file=tempfile())


test_that("clara clustering is stable", {
  expect_that(
    as.numeric(b$data$clara), equals(c(1, 2, 3, 4, 5, 6, 7, 8, 7, 9, 10, 5, 11, 12, 10, 13, 14, 15, 16, 17, 18, 19, 19, 20, 21, 22, 23, 24, 25, 25, 26, 27, 28, 4, 29, 1, 30, 31, 10, 32, 33, 34, 35, 2, 36, 37, 38, 39, 28, 40, 23, 41, 42, 43, 23, 17, 44, 45, 46, 47, 4, 48, 49, 50, 51, 52, 20, 37, 20, 5, 53, 28, 7, 54, 7, 55, 42, 9, 56, 57, 58, 44, 59, 12, 60, 39, 61, 37, 59, 62, 63, 59, 64, 60, 65, 66, 9, 67, 68, 69, 70, 6, 71, 32, 72, 29, 12, 65, 60, 21, 4, 72, 32, 73, 44, 60, 54, 74, 75, 13, 70, 11, 2, 76, 77, 1, 35, 51, 78, 25, 42, 67, 19, 4, 76, 72, 62, 72, 79, 67, 67, 80, 60, 75, 74, 76, 31, 81, 25, 44, 82, 83, 75, 41, 39, 31, 10, 2, 37, 84, 52, 60, 41, 65, 55, 12, 33, 85, 25, 5, 80, 70, 2, 11, 56, 6, 55, 10, 21, 19, 51, 44, 86, 4, 48, 67, 41, 60, 84, 1, 69, 58, 80, 57, 87, 3, 88, 89, 48, 23, 39, 62, 2, 76, 1, 11, 48, 54, 5, 90, 91, 92, 67, 93, 91, 13, 69, 85, 23, 48, 61, 94, 18, 18, 95, 96, 9, 67, 86, 33, 97, 98, 67, 47, 99, 80, 48, 53, 13, 100, 61, 73, 1, 61, 57))
  )
})

test_that("hclust clustering is stable", {
  expect_that(
    as.numeric(b$data$cluster), equals(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2))  
  )
})

unlink(temp)