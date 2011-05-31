############################################################################
# Set working directory
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
	total.env <- load.env(dataset=NA, env.dat=T, path="../Antarctic/")
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

dataset.file = "spider.comm.RData"

# either read it from the saved copy (this is faster than re-reading it every time)
if (file.exists(dataset.file)) {
	load(dataset.file)

# or extract it from the files
} else {
	# read your dataset
	dataset <- read.csv("../Data/Austropallene.csv")
	# optionnally, get a few informations about the dataset
	# names(dataset)
	# summary(dataset)
	# make sure the first two collumns are named "lat" and "long"
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
	save(dataset, file=dataset.file)
}



############################################################################
# Do BRT on important species
############################################################################

# Environmental variables used for the prediction
# all of them
pred.vars <- names(env.data)[! names(env.data) %in% c("lat","long")]
# a few
# pick in the following (except lat and long)
# names(env.data)
# pred.vars <- c("bathymetry", "bathymetry_slope", "floor_temperature")

# Response variables = taxa of interest
# all of them
# resp.vars <-names(dataset)[! names(dataset) %in% c(names(env.data),"lat","long")]
# a few
# pick in the following
names(dataset)[! names(dataset) %in% names(env.data)]
resp.vars <- c("Austropallene.cornigera")

# Subsample one in 5 rows for now, for speed purposes
env.data <- env.data[seq(1, nrow(env.data), by=5),]

# Get suppport functions
source("base scripts/gbm.functions.2010SM.r")
source("base scripts/brt.function.r")

result <- do.brt(dat=dataset, resp.vars=resp.vars, predvar=pred.vars, int=2, distrib="bernoulli", wghts=NULL, monotone=NULL, n.boot=NA, plotname="spider.brt.effects", image.name="spider.brt.file", n.pred=1, pred.data=env.data)
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
# plotname    name the plots are saved under
#             ignore the extension, it will be PDF
#             can be a relative path or full path
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
