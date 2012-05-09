
# HOWTO run a regionalisation analysis

## Using the GUI

A custom graphical interface should ease the use of the underlying functions, which are documented in the next paragraph. To display this GUI, follow the instructions in the README document to get started and then type

	do.bioreg()

in the console.

This presents a window with a list of variables. Select the ones you want to use as input variables to the regionalisation analysis (use control-leftclick to select multiple variables). Once you have selected the variables of interest, click the button at the bottom of the window and select a suitable directory for storing the output files.

A new window will then open, with a number of options.

### Options

#### Weighting

This allows you to control the relative influence of each variable. By default, all variables have equal weight (weight==1). If you wish to increase (decrease) the influence of one or more variables, you can increase (decrease) their associated weight values. During processing, weights are always scaled relative to the maximum weight specified (i.e. the maximum weight is scaled to 1, and the others are scaled accordingly).

### Transformations

This allows you to specify transformations that will be applied to each variable. Transformation functions should be written as valid R expressions, using `x` as the variable name. For example:

*     `log(x)` : will take a simple logarithm (base e)
*     `log10(x)` : logarithm (base 10)
*     `sqrt(x)` : square root
*     `x^0.25`	: fourth-root (root-root) transform

Compound expressions can also be used. For example, the bathymetry variable gives the elevation of each pixel, with negative values indicating lower than sea level. Land will have elevation values above zero. Thus, we probably want to discard positive values, and log-transform the negative values. This can be done like so:
*	 `x[x>=0]=NA; log10(-x)` : first set non-negative values to NA, then change sign and log10-transform the negative values

Some example transformations are included in the right-hand list. You can copy and paste these into the "Transformations" list.


### Analysis type

*   `Exploratory` : runs more quickly and produces a map that is not quite so pretty
*   `Final` : uses more rigorous settings in the clustering step, and produces nicer looking maps (but is much slower)

Use the `exploratory` setting to decide on the variables to use and other settings, then run a `final` quality analysis. Note that the clustering results will change slightly with `final` quality selected.

### Geographic range

The geographic range determines the spatial extent of the analyses. They map to the arguments `lat/lon.min` and `lat/lon.max`.

### Number of clusters

Use this slider to select the number of clusters (the level at which to cut the dendrogram from the hierarchical clustering step).

### Running the analysis

The `Start again` button returns you to the previous (variable selection) screen.

The `Run` button will run the analysis using your chosen settings. The R console should then print information as the analysis progresses. The GUI should still be open, so parameter values can be tweaked and another run started, if need be.


## Using the command line

The function `bioreg()` is the direct equivalent of the GUI. It accepts the following arguments and default values:

*	`variables`	a vector of variable names to use (e.g. "c('bathymetry','sst_summer_climatology')")

*	`n.groups=12`	the number of clusters to choose (the level at which to cut the dendrogram from the hierarchical clustering step). By default 12

*	`lat.min/max/step`	definition of the latitude coordinate of the prediction grid (minimum, maximum and step). By default, from -30º to -80º with a 0.1º step.

*	`lon.min/max/step`	definition of the longitude coordinate of the prediction grid (minimum, maximum and step). By default, from -180º to -180º with a 0.5º step.

*	`path="env_data"`	path to the environmental data location (see `HOWTO deal with environmental data.txt` for more information on the environmental data).

*	`transformations`	a list of transformations (strings) to apply to each variable. Use an empty string if no transformation is required for a given variable. e.g. "transformations=list('log10(x)','')" will apply a log10 transformation to the first variable and no transformation to the second.

*	`weights`		a vector of weights to control the relative influence of each variable. By default, all weights are equal.

*	`quality="low"`		`quality="low"` is faster and suitable for exploratory analyses; use `quality="high"` for final analyses.

*	`output.dir`		destination for output files; if NULL, no output files will be saved. Default is NULL.


It returns a data frame with the cluster number appended as an additional variable. The variable `clara.num` is the cluster number from the intermediate (non-hierarchical) clustering step. The final cluster number is indicated by the variable `cluster`.

