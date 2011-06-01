############################################################################
# Set working directory
############################################################################

# Adapt the path to your computer setup
# setwd("C:\\Projects\\Specific projects\\2011\\Bioregionalisation\\R")

# Name of the file containing your dataset
# It must be a CSV (comma separated values) file
# with lat and long as the first and second columns
dataset.filename = "../Data/Austropallene.csv"

############################################################################
# Load environmental data for the full grid
#  = grid on which the predictions will be made
############################################################################

env.data.saved.copy = "env.data.RData"

# when the saved copy exists, read it
# (this is faster than re-reading it every time)
if (file.exists(env.data.saved.copy)) {
	load(env.data.saved.copy)

# if there is no saved copy, extract the data from the original files
} else {
	source("base scripts/load env data.r")
  # total.env <- load.env(dataset=NA, env.dat=T, path="../Antarctic/")
  # NB !!!!! For now we subsample with a low precision
	total.env <- load.env(dataset=NA, env.dat=T, path="../Antarctic/", lat.step=1, long.step=2)
	env.data <- total.env[["env.dat"]]

	# remove points on land
	env.data <- env.data[env.data$bathymetry < 0, ]

	# save it for later
	save(env.data, file=env.data.saved.copy)
}

############################################################################
# Load environmental data for your data
#  = data points on which the model is built
############################################################################

dataset.saved.copy = "spider.comm.RData"

# when the saved copy exists, read it
# (this is faster than re-reading it every time)
if (file.exists(dataset.saved.copy)) {
	load(dataset.saved.copy)

# if there is no saved copy, extract the data from the original files
} else {
	# read your dataset
	dataset <- read.csv(dataset.filename)
	# optionnally, get a few informations about the dataset
	# names(dataset)
	# summary(dataset)
	# make sure the first two columns are named "lat" and "long"
	names(dataset)[1:2] <- c("lat","long")

	# load the corresponding environmental data
	source("base scripts/load env data.r")
	dataset.env <- load.env(dataset=dataset, env.dat=NA, path="../Antarctic/")
	dataset <- dataset.env[["dataset"]]
	# optionnally inspect the content of the new dataset, with environmental data added
	# names(dataset)
	# summary(dataset)

	# remove (completely) points which are on land
	# those points are probably coastal points for which the grid cell
	# is mostly land, so the average altitude is positive and the rest
	# of the environemental data is probably representing mostly land also
	dataset <- dataset[dataset$bathymetry<0,]

	# save the dataset for later
	save(dataset, file=dataset.saved.copy)
}


############################################################################
# Do BRT on important species
############################################################################

# Environmental variables used for the prediction
# all of them
# pred.vars <- names(env.data)[! names(env.data) %in% c("lat","long")]
# a few
# pick in the following (except lat and long)
names(env.data)
pred.vars <- c("bathymetry", "bathymetry_slope", "floor_temperature")

# Response variables = taxa of interest
# all of them
# resp.vars <-names(dataset)[! names(dataset) %in% c(names(env.data),"lat","long")]
# a few
# pick in the following
names(dataset)[! names(dataset) %in% names(env.data)]
resp.vars <- c("Austropallene.cornigera")

# Get suppport functions
source("base scripts/gbm.functions.2010SM.r")
source("base scripts/brt.function.r")
source("base scripts/polar.ggplot.R")

result <- do.brt(dat=dataset, resp.vars=resp.vars, predvar=pred.vars, int=2, distrib="bernoulli", wghts=NULL, monotone=NULL, n.boot=NA, plot.name="BRT", image.name="spider.brt.file", n.pred=1, pred.data=env.data)
# Where:
# dat         dataset containing coords, species presence, environmental data
# resp.vars   names or indexes of response variables (species to model)
# predvar     names of the predictive variable (environmental data)
# int         number of interactions, or tree complexity
# distrib     distribution family: bernoulli (= binomial), poisson,
#             laplace or gaussian
# wgths       name of the variable in `dat` which corresponds to weights
# monotone    a vector of 0 and 1 of the lengths of predvar
#             1 is monotone increasing
#             -1 monotone decreasing
#             0 is arbitrary
# n.boot      number of bootstraps for the predictions, NA means single
#             point estimate. Beware, bootstraps are slow to run and
#             less than 50-100 bootstraps will fail.
# plot.name   prefix for the name of plots
#             can be a relative path or full path
#             plots will be saved as PDF files, one for each response variable
#             if NULL plots are only displayed to screen and not saved
# image.name  name the R workspace is saved under
#             ignore the extension, it will be RData
#             can be a relative path or full path
# n.pred      number of bootstraps for predictions
#             1 if no bootstraps wanted and NA if no predictions wanted
# pred.data   predictive dataframe, needs lat and long and environmental
#             variables defined in predvar

# Inspect the resulting object
str(result, 2)
# deviance explained and a few other informations
result[[1]]$deviance
# contribution of predictive variables
result[[1]]$obj$contributions
# predictions
head(result[[1]]$pred)

# Redo the plot manually
polar.ggplot(result[[1]]$pred, mapping=aes(colour=pred), geom="points")
# plot tiles and subsample the predicted data
polar.ggplot(result[[1]]$pred, mapping=aes(fill=pred), geom="tile", lat.precision=1, lon.precision=2)
