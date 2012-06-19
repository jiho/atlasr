#
#      Test BRT functions
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

context("BRT computation")

# prepare fake data
set.seed(123)
n <- 100
data <- data.frame(lon=runif(n, -180, 180), lat=runif(n, -80, -30), sp1=round(runif(100)), sp2=round(runif(100)))
temp <- tempfile(fileext=".csv")
write.table(data, file=temp, row.names=FALSE, sep=",")

# run a BRT model
capture.output(suppressMessages(
  b <- brts(file=temp, taxa="sp1", variables=c("bathymetry", "ssh"), path=path, quiet=T, predict=TRUE, lat.step=2, lon.step=4)
), file=tempfile())


test_that("contributions are replicable", {
  expect_that(
    as.numeric(b$sp1$contributions), equals(c(61.9345447875574, 38.0654552124426), tolerance=10^-13)
  )
})

test_that("deviance diagnostics are replicable", {
  expect_that(
    as.numeric(b$sp1$deviance), equals(c(0.02, 0.25), tolerance=10^-4)
  )
})

test_that("predictions are replicable", {
  pred <- b$sp1$prediction$pred
  expect_that(
    as.numeric(summary(pred)), equals(c(0.3744, 0.4597, 0.5039, 0.502, 0.5754, 0.5942), tolerance=10^-4)
  )
})

outDir <- str_replace(temp, ".csv", "")
outDir <- str_c(outDir, "/sp1-BRT/")
outPrefix <- str_c(outDir, "sp1-BRT.")

test_that("produces all output files", {
  # plots
  expect_that( file.exists(str_c(outPrefix, "pdf")), is_true() )
  # shapefiles
  expect_that( file.exists(str_c(outPrefix, "shp")), is_true() )
  expect_that( file.exists(str_c(outPrefix, "dbf")), is_true() )
  expect_that( file.exists(str_c(outPrefix, "shx")), is_true() )
  expect_that( file.exists(str_c(outPrefix, "prj")), is_true() )
  # Rdata archive
  expect_that( file.exists(str_c(outPrefix, "Rdata")), is_true() )
  # csv file
  expect_that( file.exists(str_c(outPrefix, "csv")), is_true() )
  
  # count number of files
  expect_that( length(list.files(outDir)), equals(7) )
})

unlink(temp)
unlink(outDir)
