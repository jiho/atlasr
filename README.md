
# How to get started

The scripts are all written in R and should work on all platforms, provided that you have the right setup.


## Pre-requisites

### Install R and associated tools

Go to the [R download site](http://cran.at.r-project.org/) and download R for your system.

*	Mac OS X users, you should get the `R-2.**.pkg` as well as `tcltk-**.pkg` in the `tools` subdirectory further down the page. Installing is just a matter of double-clicking and saying `OK` to everything.

*	Windows users, you should select the `base` subdirectory and then the main `Download R 2.** for Windows` link. Install with the default choices everywhere.

*   Linux users, you should get it from your package management system rather than from the above site. The pacckage is usually called R-cran.

We recommend using a specialized editor on top of R, which will become your only interface to R. Download and install [RStudio](http://rstudio.org/); on the download page, select `RStudio Desktop` and it should detect your operating system and present the appropriate file for you on the next page.


### Download the latest version of the scripts

These scripts are rapidly evolving (at times...) and the latest version is always kept on this page (https://github.com/jiho/atlasr). You can download a zip file using the `ZIP` button or directly at the address:

https://github.com/jiho/atlasr/zipball/master

Unzipping the file should create an `atlasr` directory with all the scripts.


## Run an analysis

### With RStudio

Double-click the `atlasr.Rproj` project file inside the `atlasr` directory. This will start RStudio, R, load all functions and update the environemental data repository. The first time you do this, it will be quite long (over an hour) because it will download all necessary R packages and all environment data layers (well over a Gb of data).

In R's console (left side of the RStudio interface) type the command you want to run. These are for example

    do.bioreg()
    do.brt()
    do.gdm()

### Without RStudio

Open R as you normally would.

Tell it where to fetch and store the environemnt data layers:

    options(atlasr.env.data="/path/to/somewhere")

You can write this line in your personal `.Rprofile` file to avoid having to repeat it everytime you start R.

Load all functions and update the environment data by running

    source("/path/to/atlasr/intitialize.R")

The first time you do this, it will be quite long (over an hour) because it will download all necessary R packages and all environment data layers (well over a Gb of data).

Next type the command you want to run, for example

    do.bioreg()
    do.brt()
    do.gdm()

### Read documentation

For each analysis, documentation is provided in the  [`documentation`](https://github.com/jiho/atlasr/tree/master/documentation) directory.
