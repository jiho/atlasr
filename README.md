
# How to get started

The scripts are all written in R and should work on all platforms, provided that you have the right setup.


## Pre-requisites

### Install R and associated tools

Go to the [R download site](http://cran.at.r-project.org/) and download R for your system.

*	Mac OS X users, you should get the `R-2.**.pkg` as well as `tcltk-**.pkg` in the `tools` subdirectory further down the page. Installing is just a matter of double-clicking and saying `OK` to everything.

*	Windows users, you should select the `base` subdirectory and then the main `Download R 2.** for Windows` link. Install with the default choices everywhere.

*   Linux users, you should get it from your package management system rather than from the above site.

Optionally (but this is recommended as it makes things easier in the following), you can install a specialized editor on top of R, which can become your only interface to R. We recommend [RStudio](http://rstudio.org/). On the download page, select `RStudio Desktop` and it should detect your operating system and present the recommended file for you on the next page. If you already use RStudio, **ensure that you have the latest version installed**.


### Download the latest version of the scripts

These scripts are rapidly evolving (at times...) and the latest version is always kept on this page (https://github.com/jiho/atlasr). You can download a zip file using the `ZIP` button or directly at the address:

https://github.com/jiho/atlasr/zipball/master

Unzipping the file should create an `atlasr` directory with all the scripts.


## Running an analysis

### Start R and get all supporting functions

In the just-created `atlasr` directory, open `initialize.R` with RStudio (right-click and `Open With...` if necessary). You should get a four panes window with the script on the top-left, the R console on the bottom-left and some stuff we don't care about yet on the right. 

Then, click on the `Source` button at the top of the file pane. This should write something starting with `source(` in the console. This command executes `initialize.R` which loads all functions stored in the `library` folder, installs all the R "packages" that they use, and fetches the environment data layers. This is long the first time you run it because it needs to download everything, but it should be instantaneous afterwards.

With any kind of R interface (not just RStudio), you could just write

	source("initialize.R")

in the console to achieve this.


By default, R works in your home directory (`Documents and Settings/bla bla/something` on Windows, `/Users/yourself` on Mac OS X, `/home/yourself` on Linux). If you would rather work somewhere else, you can use the menu `Tools > Set Working Directory > Choose Directory` in RStudio to tell it where to work. It should print something starting with `setwd(` in the console. You have to do this every time you open R.

If you want to work in the `atlasr` directory directly, keeping the code and the data together, you can even use the shortcut `Tools > Set Working Directory >  To Source File Location` in RStudio, which saves you from having to browse through your computer.

You can now close the `initialize.R` file.

### Run an analysis

In the console, type the command to run your analysis. For example, to run a BRT model using the automatic graphical interface, type

	do.brt()

Read the additional documentation for each analysis in the [`documentation`](https://github.com/jiho/atlasr/tree/master/documentation) directory. When you are done, quit RStudio. It asks wether you want to save an image of your session. You do not need to save anything (all results are saved appropriately when you run each analysis).
