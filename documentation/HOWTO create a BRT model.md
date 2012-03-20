
# HOWTO create a BRT model

Boosted Regression Trees use the [bgm](http://cran.r-project.org/web/packages/gbm/) and [dismo](http://cran.r-project.org/web/packages/dismo/index.html) R packages.


## Formatting your data

The data should be formatted in a table with stations as lines and one column per taxon, plus two additional columns with latitude and longitude of stations. Further columns will not fail but will show up as additional taxa.

The first line should contain column headers. The order of columns does not mater, only their header does. The header of the latitude column should be "lat" or "latitude". The header of the longitude column should be "lon", "long", or "longitude"

The file can be formatted as :

*	`.csv`: Comma Separated Values (French users should make sure an actual comma "," is used, and not a semicolon ";" as is customary in France) ;

*	`.txt`: Text file, with columns separated by tabulations (copy-paste from a spreadsheet into a text file should give this format) ;

*	`.xls`: Excel file (on OS X / Linux only) with data in the first sheet.


## Using the GUI

A custom graphical interface should ease the use of the underlying functions, which are documented further down. To display this GUI, follow the instructions in the README document to get started and then type

	do.brt()

in the console.

This presents a window with just one button labelled `Choose file`. Clicking on the button opens the standard file chooser on your system and should allow you to navigate to your data file.

Once the data file is chosen, the interface will present two lists (taxa on the left, environmental variables on the right), several options, a geographic grid picker, and, at the bottom, a row of buttons.

Check the taxa you want to study. BRTs will be run independently and successively for all of them. Results will also be save independently for each taxon.

Select the environmental variables you want to include in the model. This list works as a file browser: press CTRL (COMMAND on Mac OS X) to select several variables, press SHIFT to select a range. Ticking the checkbox below the list does not visually change the selection (this is an unsolvable bug) but it actually does select all variables. This can be useful as a first step but beware that several variables are repeated (in raw and interpolated forms) in the default environmental list.

The options (lists and checkboxes) map to the following arguments of the command line functions, which are explained in the next paragraph :

*	Distribution : `family`

*	Prediction : `predict` (and `n.boot.pred=200` when bootstrapping)

*	Bootstrap effects : `n.boot.effects=200` when checked

*	Bin original data : `bin=TRUE` when checked

*	Extrapolate environmental range : `extrapolate.env=TRUE` when checked

*	Subsample prediction plot : `quick=TRUE` when checked

*	Overlay stations : `overlay=TRUE` when checked

The geographic range determines both which original data points are used and where the prediction is to be made. They map to the arguments `lat/lon.min` and `lat/lon.max`. The steps in lon and lat define the precision of the prediction grid and the bin size when binning the original data (when `bin=TRUE`). They map to the arguments `lat/lon.step`.

The `Help` button currently does not do much. The `Cancel` button closes the interface (but does not stop the analysis if one was started). The `Run` button runs the analysis with the currently selected argument values.

The R console should then print information along the computation of the model(s). It will write "Done" when the run is over. The GUI should still be open, so parameter values can be tweaked and another run started, if need be.


## Using the command line

The function `brts()` is the direct equivalent of the GUI. It accepts the following arguments and default values:

*	`file`	the path and name of the data file.

*	`taxa`	a vector of taxa names as character strings. Those will be matched against the actual names in the data file and can therefore be abbreviated. For example, if the taxa in the data file are "calanus", "calanus_propinquus" and "calanus_simillimus", `taxa="pro"` will match "calanus_propinquus", `taxa="cala"` will match all three, but  `taxa="calanus"` will match only "calanus" (because the name is specified completely).

*	`variables`	a vector of environment variables names as character strings. Like `taxa` they can be abbreviated and will be matched against names of variables in the the environmental dataset.

*	`family=c("bernoulli", "gaussian", "poisson")`	the distribution of data for each taxon. "bernoulli" is for presence/absence, "guassian" and "poisson" are for abundances.

*	`n.boot.effects=0`	number of bootstraps in the estimation of the effects of environmental data on the targeted taxa. This allows to judge the robustness and "significativity" of the predicted effect. Less than 50-100 will probably fail. Beware, these are long to run.

*	`bin=FALSE`	wether to bin the original observations spatially, on the prediction grid. Some locations might be sampled more often than others for practical reasons. They will be associated with the same values of the environmental data when those will be extracted from the environmental database, and those values will have increased weight in the definition of the environmental niche for the taxa. This is likely to be purely a bias of the sampling. Binning the data summarizes information per grid cell: with several presence/absence observations, a presence is recorded for the cell when at least one of the observation records a presence; with several abundance observations, the mean abundance is used. 

*	`predict=FALSE`	whether to predict the probability of presence / abundance of the species on the prediction grid. When this is `FALSE`, only the effects of environmental parameters are estimated.

*	`n.boot.pred=0`	number of bootstraps for the prediction. This allows to estimate a cross-validation error on the prediction. Less than 50-100 will fail. Beware, these are long to run.

*	`extrapolate.env=FALSE`	whether to extrapolate the prediction outside of the *environmental* range of the original data. The environmental range of the original data is the range of values associated with the observation points for each environmental variable. This includes points where both non-zero and zero abundances are recorded, assuming that points that were sampled could have contained the taxon, even though it was not actually detected. For example, samples might have been taken in temperatures ranging from 2 to 4ºC only. From this data alone, it might be difficult to extrapolate to what should happen at 5 or 6ºC.  
    When `TRUE`, predictions are made for all points of the predictions grid. When `FALSE`, the predictions are deleted for points whose associated environmental data is outside the range observed in the original data (and not plotted). When `NA`, the predictions at such points are replaced by `NA` (and plotted as grey).

*	`quick=TRUE`	wether to subsample the predicted data on a 1 x 2º lat/lon grid before plotting, to speed things up.

*	`overlay.stations=FALSE`	wether to overlay stations in the original dataset on top of the prediction plot.

*	`lat.min/max/step`	definition of the latitude coordinate of the prediction grid (minimum, maximum and step). By default, from -30º to -80º with a 0.1º step.

*	`lon.min/max/step`	definition of the longitude coordinate of the prediction grid (minimum, maximum and step). By default, from -180º to -180º with a 0.5º step.

*	`path="env_data"`	path to the environmental data location (see `HOWTO deal with environmental data.txt` for more information on the environmental data).

*	`plot.layout=c(2,2)`	layout of the matrix of effects plots (number of lines, number of columns).

*	`quiet=FALSE`	wether to remove the messages printed to the console along the analysis.


It returns a named list with one element per taxa studied. Each element is an object of class `brt`, i.e. a list with components:

*	`data`  the original `data.frame`, with associated environmental data.

*   `resp.var`  the name of the response variable, i.e., the taxon.

*   `pred.vars` the names of the prediction variables, i.e., the environmental variables.

*   `obj`  the object resulting from the model fit. See the help of `gbm.step`, `gbm`, `gbm.object` in the `dismo` package for more information.

*   `deviance`  a list with information about the deviance explained, namely

    *   `perc.deviance.explained`  the percentage of deviance explained by the model;

    *   `AUC`  the Area Under the receiver operating characteristic Curve.

*   `contributions`  named vector of the percentage of contributions of all prediction variables to the model (sum to 100)

*   `prediction`  `data.frame` with the prediction grid (`lat`, `long`), the predicted probability of presence / abundance (`pred`) at the grid points and possibly the cross validated error (`CVpred`).

*   `plot.pred`  a ggplot2 object with the prediction plot.


## Analysis routine

do not predict with many parameters
trim down parameters
do not predict and boot strap effects
trim down parameters
predict with quick plot
check
predict with bootstrap and full plot


## Viewing and interpreting results


## Inspecting the results again
