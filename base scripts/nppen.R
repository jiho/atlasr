#
#     Non-Parametric Probabilitic Ecological Niche
#
#   Compute the probability of presence of a given taxon based on
#   the environmental parameters
#
#   Beaugrand, G, Lenoir, S, Ibañez, F, and Manté, C
#   A new model to assess the probability of occurrence of a species,
#   based on presence-only data
#   Marine Ecology Progress Series, 424:175--190, 2011
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

# nppen.stats <- function(X, Y) {
# 
#   # vector sizes
#   n = nrow(X)
#   m = nrow(Y)
# 
#   # compute distance between all points of Y and the centroid of X
#   d2 = mahalanobis(Y, mean(X), cov(X))
# 
#   # compute probabilities of occurrence of those distances
#   # to do this compute the distance between each element of X and the rest of X+each element of Y. For each element of Y count how many of these permutations lead to a distance larger than the one computed.
#   p = c()
#   for (iy in 1:m) {
#     y = Y[iy,]
#     d2m = c()
#     for (ix in 1:n) {
#       Z = rbind(y, X[-ix,])
#       d2m[ix] = mahalanobis(X[ix,], mean(Z), cov(Z))
#     }
#     p[iy] = sum(d2[iy] > d2m)/n
#   }
# 
#   return(p)
# }

nppen <- function(X, Y)
# X   matrix or data.frame containing the environmental data at
#     locations of presence
# Y   matrix or data.frame containing the environmental data at
#     the points where the probability of presence needs to be
#     predicted (i.e., usually on a grid)
{

	# Note about Mahalanobis computation:
	#
	#   mahalanobis(Z[1,], mean(Z), cov(Z))
	# is equivalent to
	#   Zs = scale(Z)
	#   k = Zs[1,]
	#   k %*% solve(cov(Zs)) %*% k
	# which is a form close to what Ibanez uses

	# The computation here is the Ibanez form, which gives slightly different distances but exactly the same probabilities

	# vector sizes
	n = nrow(X)
	m = nrow(Y)

  d2 = c()
	p = c()
	for (iy in 1:m) {
		# distance from the current element of Y to centroid of X
		Z = rbind(Y[iy,], X)
		Zs = scale(Z)
		k = Zs[1,] - apply(Zs[-1,], 2, mean)
		invCov = solve(cov(Zs))
		e0 = as.numeric(k %*% invCov %*% k)

		d2m = c()
		for (ix in 1:n) {
			# distance from the current element of X and the centroid of (the rest of X + the current y)
			k = Zs[1+ix,] - apply(Zs[-(1+ix),], 2, mean)
			d2m[ix] = k %*% invCov %*% k
		}

    # store the mahalanobis distance
    d2[iy] = e0
    
		# compute the probability
		p[iy] = sum(d2m >= e0)/n
	}
	return(data.frame(d2, p))
}

do.nppen <- function(dat, resp.vars, pred.data, pred.vars, plot.name=NULL)
# dat         dataset containing coords, species presence, environmental data
# resp.vars   names or indexes of response variables (species to model)
# pred.data   predictive dataframe, needs lat and long and environmental
#             variables defined in predvar
# pred.vars   names of the predictive variable (environmental data)
# plot.name   prefix for the name of plots
#             can be a relative path or full path
#             plots will be saved as PDF files, one for each response variable
#             if NULL plots are only displayed to screen and not saved
{
  
  result = list()

  for (resp in resp.vars) {
    
    # Keep environemental data only where there are presences
    X = dataset[dataset[resp] > 0 , pred.vars]

    # Extract the predictive variables we are interested in
    Y = env.data[c("lat", "long", pred.vars)]

    # Remove NAs
    X = na.omit(X)
    Y = na.omit(Y)

    # Compute the environmental distance and probability of presence
    pred = nppen(X, Y[pred.vars])

    # Add coordinates
    pred = data.frame(Y[,c("lat","long")], pred)

    # Plot the results
    # compute a measure of 'confidence' based on Mahalanobis distance
    pred$confidence = sqrt(pred$d2)
    pred$confidence = 1 - pred$confidence/max(pred$confidence)
    # set up plots
    suppressMessages(require("ggplot2"))
    prob.plot = polar.ggplot(pred, aes(colour=p)) + opts(title=paste(resp, " - presence"))
    conf.plot = polar.ggplot(pred, aes(colour=confidence)) + opts(title=paste(resp, " - confidence"))
    
    if (!is.null(plot.name)) {
      # if a plot.name is provided, plot to a file
      pdf(paste(plot.name, "-", resp, ".pdf", sep=""), width=12, height=6)
    }
    # set up the plots side by side
    grid.newpage()
    suppressMessages(require("grid"))
    pushViewport(viewport(layout = grid.layout(1, 2)))
    print(prob.plot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
    print(conf.plot, vp=viewport(layout.pos.row=1, layout.pos.col=2))

    if (!is.null(plot.name)) {
      dev.off()
    }
    
    # Store probabilities and plots in the result object
    result[[resp]] = list(pred, prob.plot, conf.plot)
  }
  
  return(result)
}