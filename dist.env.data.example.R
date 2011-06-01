############################################################################
# This script will calculate and plot the environmental distance between the
# data used for the model and the prediction dataset
# need to use the appropriate dataset and env.data below
# Sophie Mormede, May 2011, s.mormede@niwa.co.nz
############################################################################

# Adapt the path to your computer setup
# setwd("C:\\Projects\\Specific projects\\2011\\Bioregionalisation\\R")


############################################################################
# Load environmental data for the full grid
#  = grid on which the predictions will be made
############################################################################

env.data.file = "env.data.RData"

# either read it from the saved copy (this is faster than re-reading it every time)
if (file.exists(env.data.file)) {
  load(env.data.file)

# or extract it from the files
} else {
  source("base scripts/load env data.r")
  # total.env <- load.env(dataset=NA, env.dat=T, path="../Antarctic/")
  # NB !!!!! For now we subsample with a low precision
  total.env <- load.env(dataset=NA, env.dat=T, path="../Antarctic/", lat.step=1, long.step=2)
  env.data <- total.env[["env.dat"]]

  # remove points on land
  env.data <- env.data[env.data$bathymetry < 0, ]

  # save it for later
  save(env.data, file=env.data.file)
}

############################################################################
# Load environmental data for your data
#  = data points on which the model is built
############################################################################

dataset <- read.csv("../Data/data_butch_update.csv")

  # load the corresponding environmental data
  source("base scripts/load env data.r")
  #dataset.env <- load.env(dataset=dataset, env.dat=NA, path="../Antarctic/")
  dataset.env <- load.env(dataset=dataset, env.dat=NA, path="C:\\Projects\\Specific projects\\ANT and Arctic env layers\\Antarctic\\")
  dataset <- dataset.env[["dataset"]]
  # optionnally inspect the content of the new dataset, with environmental data added
  # names(dataset)
  # summary(dataset)

  # remove (completely) points which are on land
  # those points are probably coastal points for which the grid cell
  # is mostly land, so the average altitude is positive and the rest
  # of the environemental data is probably representing mostly land also
  dataset <- dataset[dataset$bathymetry<0,]


############################################################################
# do the mahalanobis distance and plot
############################################################################

# chose the same predictive variables as you will use in the model
pred.vars <- c("bathymetry", "bathymetry_slope", "floor_temperature")
Mean <- function(x) {mean(x,na.rm=T)}

dat <- dataset[,names(dataset) %in% pred.vars]
dat <- na.omit(dat)
env <- env.data[,names(env.data) %in% pred.vars]
env <- na.omit(env)
env.dist <- env.data[,names(env.data) %in% c(pred.vars,"lat","long")]
env.dist <- na.omit(env.dist)
env.dist <- env.dist[,c("lat","long")]

temp <- mahalanobis(env,center = apply(dat,2,Mean),cov=cov(dat))
env.dist$dist <- temp
rm(temp,dat,env)


# Get suppport functions
source("base scripts/polar.ggplot.R")

# Plot and save
polar.ggplot(env.dist, mapping=aes(colour=dist), geom="points")
savePlot("env.dist.butch.png",type="png")

