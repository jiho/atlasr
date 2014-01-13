# License and citation

This code is licensed under a Creative Commons Attribution 3.0 Unported license (CC BY 3.0, http://creativecommons.org/licenses/by/3.0/). It allows you to use, share and remix the code as long as you cite us in the following way:

> J-O Irisson, S Mormede, B Raymond (alphabetical order). atlasr: R code for ecoregionalisation. https://github.com/jiho/atlasr.


The data layers downloaded automatically by this code are under a Polar Information commons license http://polarcommons.org/ethics-and-norms-of-data-sharing.php and copyright Ben Raymond. The PIC license have similar terms to the CC BY license.

The GDM computation functions (files `gdmfuncs.1.1.R` and `gdmlib.dll`) are downloaded from http://www.biomaps.net.au/gdm/, retain their original license and should be cited as advised on the original page.


# How to get started

The scripts are all written in R and should work on all platforms, provided that you have the right setup.


## Pre-requisites

### Install R and associated tools

Go to the [R download site](http://cran.at.r-project.org/) and download R for your system.

*	Mac OS X users, you should get the `R-**.pkg` as well as `tcltk-**.pkg` in the `tools` subdirectory further down the page. Installing is just a matter of double-clicking and saying `OK` to everything.

*	Windows users, you should select the `base` subdirectory and then the main `Download R ** for Windows` link. Install with the default choices everywhere.

*   Linux users, you should get it from your package management system rather than from the above site. The package is usually called R-cran.

We recommend using a specialized editor on top of R, which will become your only interface to R. Download and install [RStudio](http://rstudio.org/); on the download page, select `RStudio Desktop` and it should detect your operating system and present the appropriate file for you on the next page.


### Download the latest version of the scripts

These scripts are rapidly evolving (at times...) and the latest version is always kept on this page (https://github.com/jiho/atlasr). You can download a zip file using the `ZIP` button or directly at the address:

https://github.com/jiho/atlasr/zipball/master

Unzipping the file should create an `atlasr` directory with all the scripts.


## Run an analysis

Double-click the `atlasr.Rproj` project file inside the `atlasr` directory. This will start RStudio and load R in the appropriate location. Load all functions and update the environemental data repository by typing the following command in the "Console" section:

    source("initialize.R")
    
The first time you do this, it will be quite long (over an hour) because it will download all necessary R packages and all environment data layers (well over a Gb of data). In R's console type the command you want to run. These are for example:

    do.bioreg()
    do.brt()
    do.gdm()

### Read documentation

For some analyses, documentation is provided in the  [`documentation`](https://github.com/jiho/atlasr/tree/master/documentation) directory.
