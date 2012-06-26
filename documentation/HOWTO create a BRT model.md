
# HOWTO create a BRT model

Boosted Regression Trees use the [gbm](http://cran.r-project.org/web/packages/gbm/) and [dismo](http://cran.r-project.org/web/packages/dismo/index.html) R packages.


## Formatting your data

The input data should be formatted in a table with stations as lines and one column per taxon, plus two additional columns with latitude and longitude of stations. Further columns will not fail but will show up as additional taxa.

The first line should contain column headers. The order of columns does not mater, only their header does. The header of the latitude column should be "lat" or "latitude". The header of the longitude column should be "lon", "long", or "longitude"

The file can be formatted as :

*	`.csv`: Comma Separated Values (French users should make sure an actual comma "," is used, and not a semicolon ";" as is customary in France) ;

*	`.txt`: Text file, with columns separated by tabulations (copy-paste from a spreadsheet into a text file should give this format) ;

*	`.xls`: Excel file (readable on OS X / Linux only) with data in the first sheet.


## Using the GUI

A custom graphical interface should ease the use of the underlying functions, which are documented in the next paragraph. To display this GUI, follow the instructions in the README document to get started and then type

	do.brt()

in the console.

This presents a window with just one button labelled `Choose file`. Clicking on the button opens the standard file chooser on your system and should allow you to navigate to your data file.

Once the data file is chosen, the interface will present two lists (taxa on the left, environmental variables on the right), several options, a geographic grid picker, and, at the bottom, a row of buttons.

Check the taxa you want to study. BRTs will be run independently and successively for all of them. Results will also be saved independently for each taxon.

Select the environmental variables you want to include in the model. This list works as a file browser: press CTRL (COMMAND on Mac OS X) to select several variables, press SHIFT to select a range. Ticking the checkbox below the list does not visually change the selection (this is an unsolvable bug) but it actually does select all variables. This can be useful as a first step but beware that several variables are repeated (in raw and interpolated forms) in the default environmental list and using variables as correlated as those makes the interpretation of effects more difficult.

The options (lists and checkboxes) map to the following arguments of the command line functions, which are explained in the next paragraph :

*	Distribution : `family`

*	Bin original data : `bin=TRUE` when checked

*	Extrapolate environmental range : `extrapolate.env=TRUE` when checked

*	Quick computation : `quick=TRUE` and `n.trees.fixed=1000` when checked

*	Bootstrap effects : `n.boot.effects=200` when checked

*	Predict : `predict=TRUE` when checked

*	Bootstrap prediction : `n.boot.pred=200` when checked

*	Overlay stations : `overlay.stations=TRUE` when checked

*   Save output : `save=TRUE` when checked

The geographic range determines where the prediction is to be made. Options map to the arguments `lat/lon.min` and `lat/lon.max`. The steps in lon and lat define the precision of the prediction grid and the bin size when binning the original data (when `bin=TRUE`). They map to the arguments `lat/lon.step`.

The `Help` button currently does not do much. The `Cancel` button closes the interface (but does not stop the analysis if one was started). The `Run` button runs the analysis with the currently selected argument values.

The R console should then print information along the computation of the model(s). In particular, it shows the underlying command that is run, allowing you to copy-paste it into a file and re-run the exact same model at any time. It will write "Done" when the run is over. The GUI should still be open, so parameter values can be tweaked and another run started, if need be.


## Using the command line

The function `brt()` is the direct equivalent of the GUI. It formats the data, sends it to the function `compute.brt()`, collect the result and produces all plots and output. It accepts the following arguments and default values:

*	`file`	the path and name of the data file.

*	`taxa`	a vector of taxa names as character strings. Those will be matched against the actual names in the input data file and can therefore be abbreviated. For example, if the taxa in the input data are "calanus", "calanus_propinquus" and "calanus_simillimus", `taxa="pro"` will match "calanus_propinquus", `taxa="cala"` will match all three, but  `taxa="calanus"` will match only "calanus" (because the name is specified completely).

*	`variables`	a vector of environment variables names as character strings. Like `taxa` they can be abbreviated and will be matched against names of variables in the the environmental dataset.

*	`family=c("bernoulli", "gaussian", "poisson")`	the distribution of data for each taxon. "bernoulli" is for presence/absence, "gaussian" and "poisson" are for abundances.

*	`predict=FALSE`	whether to predict the probability of presence / abundance of the species on the prediction grid. When this is `FALSE`, only the effects of environmental parameters are estimated.

*	`lat.min/max/step`	definition of the latitude coordinate of the prediction grid (minimum, maximum and step). By default, from -30º to -80º with a 0.1º step.

*	`lon.min/max/step`	definition of the longitude coordinate of the prediction grid (minimum, maximum and step). By default, from -180º to -180º with a 0.5º step.

*	`n.boot.pred=0`	number of bootstraps for the prediction. This allows to estimate a cross-validation error on the prediction. Less than 50-100 will fail. Beware, these are long to run.

*	`n.boot.effects=0`	number of bootstraps in the estimation of the effects of environmental data on the targeted taxa. This allows to judge the robustness and "significativity" of the predicted effect. Less than 50-100 will probably fail. Beware, these are long to run.

*	`bin=FALSE`	whether to associate weights to spatial bins of the input data, with the same precision as the prediction grid. Some locations might be sampled more often than others for practical reasons. They will be associated with similar values of the environmental data and those values will therefore have increased weight in the definition of the environmental niche for the taxa. This is likely to be purely a bias of the sampling. By binning, we compute the number of data points per grid cell and associate weights inversely proportional to this number, compensating for the increased weight do to repetition.

*	`extrapolate.env=FALSE`	whether to extrapolate the prediction outside of the *environmental* range of the input data. The environmental range of the input data is the range of values associated with the observed points for each environmental variable. This includes points where both presence and absence, assuming that points that were sampled could have contained the taxon, even though it was not actually detected. For example, samples might have been taken in temperatures ranging from 2 to 4ºC only. From this data alone, it might be risky to extrapolate to what should happen at 5 or 6ºC. This is whay `extrapolate.env` is `FALSE` by default.  
    When `TRUE`, predictions are made for all points of the predictions grid. When `FALSE`, the predictions are deleted for points whose associated environmental data is outside the range observed in the input data (and not plotted). When `NA`, the predictions at such points are replaced by `NA` (and plotted as grey).

*	`quick=TRUE`	when `TRUE`, computation is made faster:
    
    -   the model is fitted with a fixed number of trees (`n.trees.fixed=1000`) instead of estimating the number of trees by a stepwise procedure; this is faster but prone to over/under-fitting.
    -   the predicted data is subsampled on a 1 x 2º lat/lon grid before plotting, to speed up plotting and reduce the size of the resulting plot file.

*   `n.trees.fixed=0`   when 0, the optimal number of trees is estimated through a stepwise procedure. When > 0 (either specified explicity or with `quick=TRUE`), a fixed numnber of trees is used.

*	`overlay.stations=FALSE`	wether to overlay stations in the input dataset on top of the prediction plot.

*   `save=TRUE` wether to save the output to files or just print info and the console and produce plots

*	`path="env_data"`	path to the environmental data location (see `HOWTO deal with environmental data.txt` for more information on the environmental data).

*	`plot.layout=c(2,2)`	layout of the matrix of effects plots (number of lines, number of columns).

*	`quiet=FALSE`	wether to remove the messages printed in the console along the analysis.


It returns a named list with one element per taxa studied. Each element is an object of class `brt`, i.e. a list with components:

*	`data`  the original `data.frame`, with associated environmental data.

*   `resp.var`  the name of the response variable, i.e., the taxon.

*   `pred.vars` the names of the prediction variables, i.e., the environmental variables.

*   `obj`  the object resulting from the model fit. See the help of `gbm.step`, `gbm`, `gbm.object` in the `dismo` package for more information.

*   `deviance`  a list with information about the deviance explained, namely

    *   `perc.deviance.explained`  the percentage of deviance explained by the model;

    *   `AUC`  the Area Under the receiver operating characteristic Curve.

*   `contributions`  named vector of the percentage of contributions of all prediction variables to the model (sum to 100)

*   `prediction`  `data.frame` with the prediction grid (`lat`, `lon`), the predicted probability of presence / abundance (`pred`) at the grid points and possibly the cross validated error (`CVpred`).


## Viewing and interpreting results

Results for each taxon are saved in a folder hierachy, next to your data. Let us say you have a file named `data_fish.csv` with presence data for species `Electrona_antarctica`. The file is on your `Desktop`. Running the a BRT model on this species will result in this hierachy

    Desktop/
        data_fish.csv
        data_fish/
            Electrona_antarctica-BRT/
                Electrona_antarctica-BRT-info.txt
                Electrona_antarctica-BRT.csv
                Electrona_antarctica-BRT.dbf
                Electrona_antarctica-BRT.pdf
                Electrona_antarctica-BRT.prj
                Electrona_antarctica-BRT.Rdata
                Electrona_antarctica-BRT.shp
                Electrona_antarctica-BRT.shx

The `info.txt` file contains basic information about the model, which was previously also printed to the screen.

The PDF file contains the plots of the effect of each variable on the probability of presence/abundance of *Electrona antarctica* and possibly the map of the prediction of its probablity of presence/abundance (if `predict=TRUE`).

The effects plots are ordered according to their contribution to the model. The curve is interpretable as a measure of the probablity of presence of the species according to the value of the variable (its ecological niche). So, when the curve is flat, the variable has little to no influence. Bootstrapping effects adds confidence intervals around theses curves, allowing to better assess which variables really have a significant effect.

The prediction probability is mapped to a colour scale. When bootstrapping predictions, the confidence in the prediction is mapped to transparency (the less confident, the more transparent). Every point of the grid would be represented unless:

*   `quick=TRUE`, in which case the data is subsampled to a 1ºx2º latxlon grid
*   `extrapolate.env=FALSE` or `NA`, in which case, points where the values of the environment are outside the environmental range of the original data are dropped or represented in grey (respectively).

The other files contain the data represented in this PDF file.

The `Rdata` file contains the R object returned by the analysis. It is called `brtObj` and loading the `Rdata` file in R (Right-click, Open With > RStudio for example) should put you in a R workspace where this object exists (check in the "Workspace" tab, on the right hand side).

The `csv` file and the shape files (`dbf`, `prj`, `shp`, `shx`) contain the predicted probability of presence at every point of the prediction grid.


## Analysis steps

A suggested analysis routine is:

1. run an initial model with no prediction, quick computation and many variables
2. remove variables which have little to no influence in the model
3. run a new model with the selected variables, full computation (quick unchecked) and bootstrapping of effects, but still with no prediction
4. with the new, bootstrapped, effects plots, trim further fown the explanatory variables
5. run a new model with the final set of variables, with prediction, on a **0.1 x 0.1 degree grid** and, possibly, bootstrapping of the prediction

It is paramount to use the same grid to fit the model and to predict. The resolution of the grid used to fit the model is the resolution of the environmental data layers (0.1 degree), therefore the final prediction has to be on a 0.1 degree grid.


## Inspecting the results again

As explained above, every element of the analysis is saved in the `Rdata` file. This file can be loaded in an R session with the `load()` command, the `Load` button in the workspace tab of RStudio or by double clicking the `Rdata` file itself, thus opening it in your R interface of choice (RStudio for example).

This gives you a `brtObj` object, which can be inspected with the usual R functions:

*   `print(brtObj)` (or simply `brtObj`) gives a few facts about the model

*   `summary(brtObj)` gives information about deviance, contributions and predictions, when available

*   `plot(brtObj)` produces the effects plot

*   `plot.pred(brtObj)` produces the prediction map

The object itself is simply a list, can be inspected using `str(brtObj)` and elements can be extracted from it with the usual `$` operator (Cf. above for the constituting elements of this object).
