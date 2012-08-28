#
#      Test BRT functions
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

context("BRT computation")

# prepare fake data
set.seed(123)
n <- 100
data <- data.frame(lon=runif(n, -180, 180), lat=runif(n, -80, -30), sp1=round(runif(100)), sp2=round(runif(100)))
temp <- tempfile(fileext=".csv")
write.table(data, file=temp, row.names=FALSE, sep=",")

# run BRT models
# optimised number of trees
capture.output(suppressMessages(
  b <- brt(file=temp, taxa="sp1", variables=c("bathymetry", "ssh"), predict=TRUE, lat.step=2, lon.step=4, quick=FALSE, path=path, quiet=TRUE, save=T)
), file=tempfile())
# fixed number of trees
capture.output(suppressMessages(
  bf <- brt(file=temp, taxa="sp1", variables=c("bathymetry", "ssh"), predict=TRUE, lat.step=2, lon.step=4, quick=TRUE, path=path, quiet=TRUE, save=T)
), file=tempfile())


test_that("contributions are replicable", {
  expect_that(
    as.numeric(b$sp1$contributions), equals(c(61.9345447875574, 38.0654552124426), tolerance=10^-13)
  )
  expect_that(
    as.numeric(bf$sp1$contributions), equals(c(52.7459199221002, 47.2540800778998), tolerance=10^-13)
  )
})

test_that("deviance diagnostics are replicable", {
  expect_that(
    as.numeric(b$sp1$deviance), equals(c(0.02, 0.25), tolerance=10^-2)
  )
  expect_that(
    as.numeric(bf$sp1$deviance$perc.deviance.explained), equals(0.62, tolerance=10^-2)
  )
})

test_that("predictions are replicable", {
  pred <- b$sp1$prediction$pred
  expect_that(
    as.numeric(summary(pred)), equals(c(0.3744, 0.4597, 0.5039, 0.502, 0.5754, 0.5942), tolerance=10^-5)
  )
  pred <- bf$sp1$prediction$pred
  expect_that(
    as.numeric(summary(pred)), equals(c(0.004359, 0.2462, 0.4592, 0.4981, 0.7714, 0.9959), tolerance=10^-5)
  )
})

outDir <- str_replace(temp, ".csv", "")
outDir <- str_c(outDir, "/sp1-BRT/")
outPrefix <- str_c(outDir, "sp1-BRT")

test_that("produces all output files", {
  # plots
  expect_that( file.exists(str_c(outPrefix, ".pdf")), is_true() )
  # shapefiles
  expect_that( file.exists(str_c(outPrefix, ".shp")), is_true() )
  expect_that( file.exists(str_c(outPrefix, ".dbf")), is_true() )
  expect_that( file.exists(str_c(outPrefix, ".shx")), is_true() )
  expect_that( file.exists(str_c(outPrefix, ".prj")), is_true() )
  # Rdata archive
  expect_that( file.exists(str_c(outPrefix, ".Rdata")), is_true() )
  # csv file
  expect_that( file.exists(str_c(outPrefix, ".csv")), is_true() )
  # info file
  expect_that( file.exists(str_c(outPrefix, "-info.txt")), is_true() )
  
  # count number of files
  expect_that( length(list.files(outDir)), equals(8) )
})

unlink(temp)
unlink(outDir)
