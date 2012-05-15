
# HOWTO run a regionalisation analysis

## Using the GUI

A custom graphical interface should ease the use of the underlying function, which are documented in the next paragraph. To display this GUI, follow the instructions in the README document to get started and then type

	do.bioreg()

in the console. This presents a window with a several elements.

### Main window

The first button allows to choose an output directory. When you do, all output will be saved there. When you don't, nothing is saved but plots are displayed after the analysis run.

Below, a list of available environmental variables is presented. Select the ones you want to use. This list works as a file browser: press CTRL (COMMAND on Mac OS X) to select several variables, press SHIFT to select a range.

Variables can be used as is, but data can also be transformed and some variables can be given more weight than others. This is done through another panel which you can open using the "Variables transformation" button. See below for its description.

The following options map directly to command-line options. See the next paragraph to get more information on their meaning.

* Analysis type: "Exploratory" sets `quick=TRUE`, "Final" sets `quick=FALSE`

* Number of clusters: `n.groups`

The geographic range determines the spatial extent of the analyses. They map to the arguments `lat/lon.min` and `lat/lon.max`. The steps in lon and lat define the precision of the prediction grid. They map to the arguments `lat/lon.step`.

The `Help` button currently does not do much. The `Cancel` button closes the interface (but does not stop the analysis if one was started). The `Run` button runs the analysis using  your chosen settings.

The R console should then print information as the analysis progresses. It will write "Done" when the analysis is over. The GUI should still be open, so parameter values can be tweaked and another run started, if need be.

### Transformation panel

#### Transformations

This allows you to specify transformations that will be applied to each variable. Transformation functions should be written as valid R expressions, using `x` as the variable name. For example:

*     `log(x)` : will take a simple logarithm (base e)
*     `log10(x)` : logarithm (base 10)
*     `sqrt(x)` : square root
*     `x^0.25`	: fourth-root (root-root) transform

Compound expressions can also be used. For example, the bathymetry variable gives the elevation of each pixel, with negative values indicating lower than sea level. Land will have elevation values above zero. Thus, we probably want to discard positive values, and log-transform the negative values. This would be done by using `x[x>=0]=NA; log10(-x)` : first set non-negative values to NA, then change sign and log10-transform the now-positive values. Actually, all data on land is automatically masked so in that particular example `log10(-x)` would be enough.

Some example transformations are included in the bottom list. You can copy and paste these into the "Transformations" list.

**WARNING** Make sure to press `Return` after typing anything in a field, otherwise the input is not recorded.

#### Weighting

This allows you to control the relative influence of each variable. By default, all variables have equal weight (weight=1). If you wish to increase (decrease) the influence of one or more variables, you can increase (decrease) their associated weight values. During processing, weights are always scaled relative to the maximum weight specified (i.e. the maximum weight is scaled to 1, and the others are scaled accordingly).

**WARNING** Similarly, make sure to press `Return` after typing anything in a field, otherwise the input is not recorded.



## Using the command line

The function `bioreg()` is the direct equivalent of the GUI. It accepts the following arguments and default values:

*	`variables`	a vector of variable names to use (e.g. `c("bathymetry","sst_summer_climatology")`)

*	`n.groups=12`	the number of clusters to choose (the level at which to cut the dendrogram from the hierarchical clustering step). By default 12

*	`lat.min/max/step`	definition of the latitude coordinate of the prediction grid (minimum, maximum and step). By default, from -30º to -80º with a 0.1º step.

*	`lon.min/max/step`	definition of the longitude coordinate of the prediction grid (minimum, maximum and step). By default, from -180º to -180º with a 0.5º step.

*	`path="env_data"`	path to the environmental data location (see `HOWTO deal with environmental data.txt` for more information on the environmental data).

*	`transformations`	a list of transformations (strings) to apply to each variable. Use an empty string if no transformation is required for a given variable. e.g. `transformations=list("log10(x)","")` will apply a log10 transformation to the first variable and no transformation to the second.

*	`weights`		a list or vector of weights to control the relative influence of each variable. By default, all weights are equal.

*   `quick=TRUE` runs more quickly and produces a map that is not quite so pretty; `quick=FALSE` uses more rigorous settings in the clustering step, and produces nicer looking maps (but is much slower).
    Use the "Exploratory" setting to decide on the variables to use and other settings, then run a "Final" quality analysis. Note that the clustering results will change slightly with "Final" quality selected.

*	`output.dir`	destination for output files; if NULL, no output files will be saved and graphics will be displayed interactively. Default is NULL.


The function produces several plots

*   the dendrogram of the clustering, showing the colour-coded clusters and the cutoff value

*   the distribution of environmental values in each cluster shown as a boxplot and a violin plot. Boxplot display, as usual, the median (thick line), quartiles (box), 1.58 * interquartile range / sqrt(n) (whiskers), and further points are considered outliers. Violin plots represent an estimate of the probability density of the data distribution; i.e., the wider the violin, the mode data.

*   a map of the clusters either as a flat image with no projection (when `quick=TRUE`) or as a true map with seterographic projection (when `quick=FALSE`)


The function returns a data frame with the cluster number appended as an additional variable. The variable `clara.num` is the cluster number from the intermediate (non-hierarchical) clustering step. The final cluster number is indicated by the variable `cluster`.


When `output.dir` is set, several files are returned:

*   `Rdata` data archive in R format which contains the data.frame returned by the function and the transformed data.
*   `pdf`   all plots saved to a PDF file
*   `shp`, `shx`, `dbf`, `prj`  shapefiles and projection (WGS 84) with the grid locations and the corresponding cluster value.
