# example for a GDM analysis
# Sophie Mormede, May 2011, s.mormede@niwa.co.nz
############################################################################


############################################################################
# Set working directory
############################################################################

# Adapt the path to your computer setup
# setwd("C:\\Projects\\Specific projects\\2011\\Bioregionalisation\\R")
source("base scripts/load env data.r")


############################################################################
# load the predictive environmental dataset
############################################################################

# if you are using the general predictive dataset
load("env.data.RData")

# OR if you want to load your own predictive grid
# env.data <- load.env(dataset=NA, env.dat=T, path="../Antarctic/",lat.lim=c(-80,-30),lat.step=1,long.lim=c(-180,180),long.step=2)
# env.data <- env.data[["env.dat"]]


############################################################################
# load your own dataset and add environmental data
############################################################################

dataset <- read.csv("../Data/data_butch_update.csv")

# load the corresponding environmental data
dataset.env <- load.env(dataset=dataset, env.dat=NA, path="../Antarctic/")
#dataset.env <- load.env(dataset=dataset, env.dat=NA, path="C:\\Projects\\Specific projects\\ANT and Arctic env layers\\Antarctic\\")
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
# Chose environmental variables and variables to be predicted
############################################################################

# Environmental variables used for the prediction
# all of them
# pred.vars <- names(env.data)[! names(env.data) %in% c("lat","long")]

# OR a few
# pick in the following (except lat and long)
names(env.data)
pred.vars <- c("bathymetry", "bathymetry_slope", "floor_temperature")

# Response variables = taxa of interest
# all of them
# resp.vars <-names(dataset)[! names(dataset) %in% c(names(env.data),"lat","long")]

# OR a few
# pick in the following
names(dataset)[! names(dataset) %in% names(env.data)]
resp.vars <- c("Electrona_carlsbergi","Electrona_antarctica", "Electrona_paucirastra", "Electrona_risso", "Electrona_subaspera")


############################################################################
# Do GDM on important species
# at the moment ggplot doesn't plot well the cluster
############################################################################

# Get suppport functions
library(MASS)
require(cluster)
source("base scripts/gdmfuncs.1.1.R")
source("base scripts/lookup.names.variables.R")
source("base scripts/gdm.function.r")

result <- do.gdm(dat=dataset,resp.vars=resp.vars,predvar=pred.vars,samp=10000,pred.data=env.data, plot.name = "mycto.gdm",n.clust=NA,do.indicator.species=T)

  ## Base case GDM running
  ## dat has your dataset of responses
  ## resp.vars are the names of different variables you want to predict
  ## predvar are your predictive variable names
  ## samp is your sample size, if too big the function will fall over. The function will randomly pick that number of samples from the dataset. Should be 10000 or more unless the computer is overloaded
  ## plot.name = what name will the plot be saved under, can include a relative path or full path, ignore extension, it will be a pdf to allow multiple pages and mac users
  ## pred.data: predictive dataframe, needs lat and long, and environmental variables as in predvar
  ## n.clust = NA if automatically chose cluster number, or give a number
  ## do.indicator.species: calculate indicator species for clusters?
  ## image.name superseded, not saved anymore, left in for code running elsewhere


# Inspect the resulting object
str(result, 2)
# deviance explained and a few other informations
gdm.summary(result$model)
# predictions
head(result$predicted)
table(result$predicted$cluster) #check if not most in one cluster
# cluster composition
result$indval

# Redo the plot manually
polar.ggplot(result$predicted, mapping=aes(colour=cluster), geom="points")
# plot tiles and subsample the predicted data
polar.ggplot(result$predicted, mapping=aes(fill=cluster), geom="tile", lat.precision=1, lon.precision=2)


# save the R Data file for future work
save.image("mycto.RData")
