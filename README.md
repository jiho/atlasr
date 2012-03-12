
# How to get started

The scripts are all written in R and should work on all platforms, provided that you have the right setup.


## Pre-requisites

### Install R and associated tools

Go to the [R download site](http://cran.at.r-project.org/) and download R for your system.

*	Mac OS X users you should get the `R-2.**.pkg` as well as `tcltk-**.pkg` in the `tools` subdirectory further down the page. Installing is just a matter of double-clicking and saying `OK` to everything.

*	Windows users you should select the `base` subdirectory and then the main `Download R 2.** for Windows` link. Install with the default choices  everywhere.

*   Linux users, you should get it from your package management system.

Optionally (but this is recommended as it makes things easier in the following), you can install a Graphical User Interface on top of R. We recommend [RStudio](http://rstudio.org/). On the download page, select `RStudio Desktop` and it should detect your operating system and present the recommended file for you on the next page.


### Download the latest version of the scripts

These scripts are rapidly evolving (at times...) and the latest version is always kept on this page (https://github.com/jiho/atlasr). You can download a zip file using the `ZIP` button or directly at the address:

	https://github.com/jiho/atlasr/zipball/master

Unzipping the file should create an `atlasr` directory with all the scripts.


## Running an analysis

### Start R and tell it were to work

In the just-created `atlasr` directory, open `get_functions.R` with RStudio (right-click and `Open With...` if necessary). You should get a four panes window with the script on the top-left, the R console on the bottom-left and some stuff we don't care about yet on the right. Use the menu `Tools > Set Working Directory > To Source File Location`. This should print something starting with `setwd(` in the console. It tells R where to find the functions. Your data can be anywhere though.

Other interfaces of R, such as the default R GUI on Windows or Mac, also provide similar functionality to set the working directory. Alternatively you can use the `setwd()` command directly. You should set it to where `get_functions.R` is.


### Get all supporting functions

In RStudio, click on the `Source` button at the top of the file pane. This should write something starting with `source(` in the console. This command loads all functions stored in the `library` folder as well as install all the "packages" that they use to get additional functionality. This is long the first time you run it because it needs to download all the packages, but should be instantaneous afterwards.

With any R interface (not just RStudio), you could just write

	source("get_functions.R")

in the console.

You can now close the `get_functions.R` file.


### Run an analysis

In the console, type the command to run your analysis. For example, to run a BRT model using the automatic graphical interface, type

	do.brt()

Read the additional documentation for each analysis in the `documentation` directory. When you are done, quit RStudio. It asks wether you want to save an image of your session. You do not need to save anything (all results are saved appropriately when you run each analysis).
