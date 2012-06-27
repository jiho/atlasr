#
#     General Dissimilarity Modelling
#
# (c) Copyright 2011-2012 S Mormede, Ben Raymond, J-O Irisson
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------

# Ben Raymond, 2011:
# note that the training data cluster labels and the indicator species info uses simple spatial interpolation, and will not work with time-varying models


## Run GDM analysis
#-----------------------------------------------------------------------------

compute.gdm <- function(
  resp.vars,          # names of the response variables (taxa)
  pred.vars,          # names of the predictor variables (environment variables)
  data,               # data frame containing locations (lat, lon), presence/absence data and environment variables
  newdata,            # prediction grid coordinates and associated environment variables
  n.groups=NULL,      # number of clusters in the result; when NULL (the default) the optimal number of clusters is determined by a stepwise procedure
  min.n.groups=3,
  max.n.groups=10,    # minimum and maximum number tested when searching for the optimal number of clusters; testing a larger range takes more time
  intern.sample=NULL, # internally, the GDM function can also subsample the pairwise dissimilarities. When NULL, no subsampling occurs. The total number of pairwise dissimilarities is n * ( n - 1 ) / 2, where n is the number of rows in `data` (or `pre.sample` when subsampling a priori), so for this to be effective, intern.sample must be lower than that. As a rule of thumb, n=2500 gives over 3 million pairwise dissimilarities
  ...                 # passed to the internal GDM function
)
{
  ## Run GDM
  #--------------------------------------------------------------------------
  message("-> Compute dissimilarities")

  # remove lines with missing data
  data <- data[,c("lon", "lat", resp.vars, pred.vars)]
  data <- na.omit(data)

  # run GDM
  m.gdm <- gdm.fit(data[,pred.vars], data[,resp.vars], sample=intern.sample, ...)


  ## Cluster the data
  #--------------------------------------------------------------------------
  message("-> Cluster based on dissimilarities")
  suppressPackageStartupMessages(require("cluster", quietly=TRUE))

  # remove NAs from prediction data
  newdata <- na.omit(newdata[, c("lat", "lon", pred.vars)])

  # predict the environment according to the GDM model
  # NB: this is done by pieces of 1000 to limit the amount of memory required
  pred <- newdata[1, pred.vars]
  for (i in seq(from=1,to=nrow(newdata)%/%1000*1000+1, by=1000)) {
    pred <- rbind(pred, gdm.transform(m.gdm, newdata[i:min(i+999,nrow(newdata)), pred.vars]))
  }
  pred<-pred[-1,]
  # TODO re-implement the 1000 split more efficiently

  if ( is.null(n.groups) ) {
    # optimise the number of clusters
    # prepare storage
    res <- array(NA,c(max.n.groups,10,4))

    # try several possibilities
    for (cl in min.n.groups:max.n.groups) {
      for (samples in 5:10) {
        for (sampsz in 2:4) {
          temp <- clara(pred, k=cl, metric="manhattan", keep.data=FALSE, samples=samples, sampsize=min(nrow(pred), 40 + sampsz * cl))
          res[cl,samples,sampsz] <- temp$silinfo$avg.width
        }
      }
    }

    # select the best one(s)
    temp <- which(res==max(res,na.rm=T), arr.ind=T)
    # when there are several, choose the first (least number of clusters)
    if (nrow(temp)>1) temp<-temp[1,]

    # recompute clustering with this optimal parameters
    cl <- clara(pred, k=temp[1], metric="manhattan", keep.data=FALSE, samples=temp[2], sampsize=min(nrow(pred), 40 +temp[3]*temp[1]))

  } else {
    # compute the clustering with the given number of groups
    cl <- clara(pred, k=n.groups, metric="manhattan", keep.data=FALSE)
  }


  ## Return the result
  #-----------------------------------------------------------------------------

  # store the cluster number in the data
  newdata$cluster <- factor(cl$clustering)

  # store elements in the resulting object
  res <- list(
    data=data,
    model=m.gdm,
    cl=cl,
    prediction=newdata
  )

  # class it as gdm
  class(res) <- c("gdm", "list")

  return(invisible(res))
}

gdm <- function(
  file,               # name of the file where the presence/abundance data is
  taxa="",            # names (or abbreviations) of the taxa of interest (by default, all columns in file except for lat and lon)
  variables,          # names (or abbreviations) of the variables to use for the prediction
  transformations=NULL,                     # named vector of transformations applied to each variable (has to match variables)
  lat.min=-80, lat.max=-30, lat.step=0.1,   # definition of the prediction grid
  lon.min=-180, lon.max=180, lon.step=0.5,
  pre.sample=2500,    # sub-sample the input data before feeding it to the GDM function. When NULL, no subsampling occurs. Subsampling reduces the number of pairwise operations needed, speeds up the function and determines the amount of memory required. Since the library is 32 bits, the amount of allocatable memory is limited and this should not be much larger than 2500
  taxa.min=1,         # minimum number of taxa with a non-null abundance or presence needed to consider a given location in the analysis
  save=TRUE,          # whether to save output to files or just print info on the console and the screen
  quick=TRUE,         # quick plot
  intern.sample=ifelse(quick, 10^5, NULL),  # set subsampling of dissimilarity matrix
  path=getOption("atlasr.env.data"),        # path to the environmental database
  ...                 # passed to compute.gdm()
)
{

  suppressPackageStartupMessages(require("plyr", quietly=TRUE))

  ## Prepare data
  #--------------------------------------------------------------------------

  # read selected variables from the database
  database <- read.env.data(variables, path=path, quiet=FALSE)
  # get the full names of those variables
  pred.vars <- names(database)

  # remove information on land
  database <- mask.env.data(database, path=path)


  # read input dataset
  file <- clean.path(file)
  if (file.exists(file)) {
    input_data <- read.data(file)
  } else {
    stop("Cannot find file : ", file)
  }

  # get the names of the taxa of interest
  allTaxa <- names(input_data[,!names(input_data) %in% c("lat", "lon")])
  resp.vars <- match.vars(taxa, allTaxa)

  # keep only those taxa
  input_data <- input_data[,c("lat", "lon", resp.vars)]

  # bin lat and lon at the same resolution as the environmental data grid
  input_data <- rasterize(input_data, vars=c("lat", "lon"), precisions=c(lat=0.1, lon=0.1), fun=function(x) {as.numeric(any(x>0))})
  # NB: for each species, this records a presence as soon as at least one presence is recorded within the grid cell. The number of absences does not count and so 1 presence for 100 absences will still be recorded as a presence. It's not perfect but it's probably the best we can do.

  # keep only locations with a given number of species present
  message("-> Compute number of species per location (min=", taxa.min,")")
  nbSp <- aaply(input_data, 1, function(x) {sum(x>0)}, .expand=FALSE)
  notEnough <- nbSp < taxa.min
  if (any(notEnough)) {
    message("   ", sum(notEnough), " / ", nrow(input_data), " locations removed")
    input_data <- input_data[!notEnough, ]
  }

  # subsample input data when required
  n <- nrow(input_data)
  if (!is.null(pre.sample) & pre.sample < n) {
    message("-> Subsample input data (choose ", pre.sample, " lines randomly)")
    input_data <- input_data[sample.int(n, pre.sample),]
  } else {
    if (n > 2600) {
      warning("Your input dataset still had ", n, " data points, which is a lot and may lead to memory errors. Consider using the argument `pre.sample` to subsample it randomly.")
    }
  }

  # get environment data for the observations
  input_data <- associate.env.data(input_data, database)

  # build prediction grid
  prediction_grid <- build.grid(
    lat.min=lat.min, lat.max=lat.max, lat.step=lat.step,
    lon.min=lon.min, lon.max=lon.max, lon.step=lon.step
  )
  # get environment data on this grid
  prediction_grid <- associate.env.data(prediction_grid, database)


  # transform environmental data if necessary
  if (!is.null(transformations)) {
    message("-> Apply transformation to environmental data")
    input_data[,pred.vars] <- transform.data(input_data[,pred.vars], transformations)
    suppressWarnings(prediction_grid[,pred.vars] <- transform.data(prediction_grid[,pred.vars], transformations))
  }

  ## Compute the GDM model and clustering
  #--------------------------------------------------------------------------
  gdmObj <- compute.gdm(resp.vars=resp.vars, pred.vars=pred.vars, data=input_data, newdata=prediction_grid, intern.sample=intern.sample, ...)


  ## Store output
  #--------------------------------------------------------------------------

  if (save) {
    # prepare names of the files where the results will be written
    # data file without extension
    fileName <- str_replace(file, "\\.(csv|xls|txt)$", "")
    baseName <- basename(fileName)

    # prepare directory
    dirName <- str_c(fileName, "-GDM")
    dir.create(dirName, showWarnings=FALSE, recursive=TRUE)
    if (!file.exists(dirName)) {
      stop("Could not produce output directory : ", dirName)
    }
    # prepare a basic name from this
    baseName <- str_c(dirName, "/", baseName, "-GDM")

    # prepare specific filenames
    pdfFile <- str_c(baseName, ".pdf")
    csvFile <- str_c(baseName, ".csv")
    rdataFile <- str_c(baseName, ".Rdata")
    infoFile <- str_c(baseName, "-info.txt")

    # print info about the fit
    # summary(brtObj)
    # TODO shorten the summary

    # write the results in files
    message("-> Write output to ", dirName)

    # RData file with the object
    save(gdmObj, file=rdataFile)

    # CSV file
    write.table(gdmObj$prediction, file=csvFile, sep=",", row.names=FALSE)

    # Shapefiles
    write.shapefile(gdmObj$prediction, baseName, c("cluster"))

    # Info file
    capture.output(summary(gdmObj), file=infoFile)

    # open the PDF
    pdf(pdfFile, width=11.7, height=8.3)
  }

  message("-> Plot results")
  plot.gdm(gdmObj, quick=quick, ...)

  if (save) {
    # close PDF
    dev.off()
  }

  return(invisible(gdmObj))
}


## Analyse GDM output
#-----------------------------------------------------------------------------

print.gdm <- function(x, ...) {
  suppressPackageStartupMessages(require("stringr", quietly=TRUE))

  species <- setdiff(names(x$data), c("lon", "lat", x$model$predictors))

  cat("\n     GDM model\n\n")

  cat("Species :\n")
  cat(" ", str_c(species, collapse="\n  "), "\n")
  cat("Predictors :\n")
  cat(" ", str_c(x$model$predictors, collapse="\n  "), "\n")

}

summary.gdm <- function(x, ...) {
  gdm.summary(x$model)
}

plot.gdm <- function(x, ...) {
  #
  # Produce all plots for a GDM object
  #
  # x   object of class gdm
  #

  if (dev.interactive() | names(dev.cur()) == "null device") devAskNewPage(TRUE)

  print(plot.effects.gdm(x, ...))

  print(plot.pred.gdm(x, ...))

  devAskNewPage(FALSE)
}

plot.effects.gdm <- function(x, ...) {
  #
  # Plot "effects" in a GDM model
  #
  # x   object of class gdm
  #
  suppressPackageStartupMessages(require("plyr", quietly=TRUE))
  suppressPackageStartupMessages(require("ggplot2", quietly=TRUE))
  suppressPackageStartupMessages(require("reshape2", quietly=TRUE))

  # compute environmental data range (omitting lat and lon)
  env.ranges <- llply(x$data[, x$model$predictors], function(x) {
    seq(from=min(x, na.rm=T), to=max(x, na.rm=T), length.out=200)
  })
  env.ranges <- data.frame(do.call(cbind, env.ranges))

  # compute the transformation induced by the model on these ranges
  env.curves <- gdm.transform(x$model, env.ranges)

  # match original ranges to the transformations
  env.ranges.m <- melt(env.ranges, measure.vars=names(env.curves), value.name="value")
  env.curves.m <- melt(env.curves, measure.vars=names(env.curves), value.name="transform")
  # and join them
  env <- data.frame(env.curves.m, env.ranges.m)

  # compute quantiles of the original observations, to plot them as a rug plot
  quant <- data.frame(llply(x$data[,x$model$predictors], function(x) {
    # to easily mix them with the others, compute 200 quantiles
    as.numeric(quantile(x, probs=seq(0, 1, length.out=200)))
  }))
  quant.m <- melt(quant, measure.vars=names(quant), value.name="quantiles")
  # put this with the rest of the data
  env$quantiles <- quant.m$quantiles

  # sort the variables in decreasing order of transform value
  maxVal <- ldply(env.curves, max, na.rm=T)
  maxVal <- maxVal[order(maxVal$V1, decreasing=TRUE),]
  # reorder factor for display
  env$variable <- factor(env$variable, levels=maxVal$.id)

  # plot
  p <- ggplot(env) +
    # the transform result
    geom_path(aes(x=value, y=transform)) +
    # the distribution of data
    geom_rug(aes(x=quantiles), alpha=0.2) +
    # for each variable
    facet_wrap(~variable, scales="free_x")

  return(p)
}

plot.pred.gdm <- function(x, quick=FALSE, overlay.stations=FALSE, geom="auto", ...) {
  #
  # Plot GDM predictions
  #
  # x                     object of class gdm
  # quick                 subsample to 1 x 2 degree in lat x lon to plot more quickly
  # overlay.observations  add points at the location of observation in to original data
  # ...                   passed on to polar.ggplot
  #

  suppressPackageStartupMessages(require("ggplot2"))

  if (quick) {
    # subsample the plot
    geom = "raster"
    lat.precision <- 1
    lon.precision <- 2
  } else {
    # do not subsample and use the default geom
    lat.precision <- NULL
    lon.precision <- NULL
  }

  # main plot
  p <- polar.ggplot(x$prediction, mapping=aes(fill=cluster), geom=geom, lat.precision=lat.precision, lon.precision=lon.precision, ...)

  # overlay stations
  if (overlay.stations) {
    p <- p + geom_point(aes_string(x=lon, y=lat), data=x$data, size=1, alpha=0.7)
  }

  return(p)
}



# ## calculate indicator species using dufrene-legendre method
# if (do.indicator.species) {
#     require(labdsv)
#     library(reshape)
#     temp <- as.matrix(cast(newdata, lat ~  long, value = "cluster",add.missing=T))
#     datlonidx <- round(approx(as.numeric(colnames(temp)),1:ncol(temp),dat$long)$y)
#     datlatidx <- round(approx(as.numeric(rownames(temp)),1:nrow(temp),dat$lat)$y)
#
#     dat$cluster <- NA
#     for (i in 1:length(datlonidx)){
#       if (!is.na(datlonidx) && !is.na(datlatidx)) {
#         dat[i,"cluster"] <- temp[datlatidx[i],datlonidx[i]]
#       }
#     }
#
#     nnanidx=which(!is.na(dat$cluster)) ## only include non-NA clusters
#     ## calculate indicator species stuff
#     frodo.baggins=indval(dat[nnanidx,resp.var.col],clustering=dat$cluster[nnanidx])
# }

# newdata$cluster<-factor(newdata$cluster)
# dat$cluster <- factor(dat$cluster)


# close the PDF file
# if (!is.null(plot.name)) {
#   dev.off()
# }
#
#   if (do.indicator.species) {
#     res=list('predicted'=newdata,'model'=no.gdm,'plot.obj'=p,'dat.cluster'=dat$cluster,'indval'=frodo.baggins)
#   } else {


## GUI
#-----------------------------------------------------------------------------


do.gdm <- function() {
  #
  # Open a GUI to select the arguments of the gdm() function
  #

  suppressPackageStartupMessages(require("rpanel"))


  # default dimensions (in px)
  w <- 600      # width of the window
  w.h <- 700    # height of the window
  h <- 50       # height of elements
  spacer <- 10  # height of spacer

  # main window
  win <- rp.control(title="Run GDM model", size=c(w,w.h))

  # NB: positions are x, y, width, eight
  #     x, y are the coordinates of the top-left hand corner

  rp.button(win, title="Choose file", pos=c(0, 0, w/4, h), action=function(win) {

    # choose the data file
    file <- file.choose()
    # file <- "../data/euphau-clean.txt"
    win$file <- clean.path(file)


    # write the filename, for information
    rp.text(win, txt=file, initval=file, pos=c(w/4, 0, 3*w/4, h))

    # build species list
    suppressMessages(data <- read.data(file))

    # the selection panel should be 380px high for the whole window to fit in 700px height
    hSel <- 380
    # this is roughly equivalent to 20 rows in the listbox
    nRows <- 20

    # we can fit ~20 checkboxes into that space, so we need to split the checkbox list into columns when there are more than that
    nBoxes <- 20
    # however, it is not possible to do that programmatically in a for loop because variables with numbered names have to be defined in each case and this requires eval(parse(...)) constructs which mostly don't work
    # so we only accommodate for a situation with two columns
    allTaxa <- setdiff(names(data), c("lat", "lon"))
    nTaxa <- length(allTaxa)

    if (nTaxa <= nBoxes) {
      # when all taxa can fit into one column, just use this
      rp.checkbox(win, var=taxa, labels=allTaxa, title="Taxa", initval=rep(TRUE, times=nTaxa), pos=c(0, h, 2*w/3, hSel), action=function(win) {
        if (sum(win$taxa) < 2) {
          rp.messagebox("At least two taxa must be selected", title="Warning")
        }
        return(win)
      })

    } else {
      # otherwise, split in half
      half <- ceiling(nTaxa/2)
      rest <- nTaxa - half
      allTaxa1 <- allTaxa[1:half]
      allTaxa2 <- allTaxa[(half+1):nTaxa]

      rp.checkbox(win, var=taxa1, labels=allTaxa1, title="Taxa", initval=rep(TRUE, times=half), pos=c(0, h, w/3, hSel), action=function(win) {
        if (all(! c(win$taxa1, win$taxa2))) {
          rp.messagebox("At least one taxon must be selected", title="Warning")
        }
        return(win)
      })
      rp.checkbox(win, var=taxa2, labels=allTaxa2, title="", initval=rep(TRUE, times=rest), pos=c(w/3, h+7, w/3, hSel-7), action=function(win) {
        if (all(! c(win$taxa1, win$taxa2))) {
          rp.messagebox("At least one taxon must be selected", title="Warning")
        }
        return(win)
      })

    }

    # build environment variables list
    allVariables <- list.env.data()
    rp.listbox.mult(win, var=variables, vals=allVariables, title="Variables",  rows=nRows, cols=22, pos=c(2*w/3, h, w/3, hSel-h), action=function(win) {
      return(win)
    })

    # transform environmental variables
    # options
    optionsFile <- tempfile()
    win <- rp.button(win, "Variables transformation", pos=c(2*w/3, hSel, w/3, h) , action=function(win) {
      if (length(win$variables) < 2) {
        rp.messagebox("At least two variables must be selected", title="Warning")
      } else {
        do.gdm.variables(variables=win$variables, file=optionsFile)
      }
      return(win)
    })


    # compute vertical coordinate of the middle of the window
    mid <- hSel + h + spacer

    # clustering options
    rp.slider(win, n.groups, from=0, to=15, resolution=1, title="Number of regions", initval=0 , showvalue=TRUE, pos=c(0, mid , w/4, h))
    rp.slider(win, min.n.groups, from=2, to=15, resolution=1, title="Min nb of regions", initval=3, showvalue=TRUE, pos=c(0, mid+h, w/4, h))
    rp.slider(win, max.n.groups, from=2, to=15, resolution=1, title="Max nb of regions", initval=10, showvalue=TRUE, pos=c(0, mid+h*2, w/4, h))
    rp.slider(win, taxa.min, from=1, to=5, resolution=1, title="Min taxa per location", initval=1, showvalue=TRUE, pos=c(0, mid+h*3, w/4, h))

    # checkboxes range
    rp.checkbox(win, pre.sample, title="Subsample data", initval=TRUE, pos=c(w/4, mid, w/4, h))
    rp.checkbox(win, intern.sample, title="Subsample\ndissimilarity matrix", initval=FALSE, pos=c(w/4, mid+h, w/4, h))
    rp.checkbox(win, save, title="Save output", initval=FALSE, pos=c(w/4, mid+h*2, w/4, h))

    # location
    rp.slider(win, lat.max,  from=-90, to=-30,  resolution=2,   title="North"   , initval=-30 , showvalue=TRUE, pos=c(w/2+w/8, mid    , w/4, h))
    rp.slider(win, lon.min,  from=-180, to=180, resolution=5,   title="West"    , initval=-180, showvalue=TRUE, pos=c(w/2, mid+h*1, w/4, h))
    rp.slider(win, lon.max,  from=-180, to=180, resolution=5,   title="East"    , initval=180 , showvalue=TRUE, pos=c(3*w/4, mid+h*1, w/4, h))
    rp.slider(win, lat.min,  from=-90, to=-30,  resolution=2,   title="South"   , initval=-80 , showvalue=TRUE, pos=c(w/2+w/8, mid+h*2, w/4, h))
    rp.slider(win, lon.step, from=0.1, to=4  ,  resolution=0.1, title="Step lon", initval=1 , showvalue=TRUE, pos=c(w/2, mid+h*3, w/4, h))
    rp.slider(win, lat.step, from=0.1, to=4   , resolution=0.1, title="Step lat", initval=1 , showvalue=TRUE, pos=c(3*w/4, mid+h*3, w/4, h))


    # action buttons
    rowY <- mid+h*4+spacer
    rp.button(win, "Help", pos=c(0, rowY, w/4, h), action=function(win) {
      rp.messagebox("Someday... Maybe", title="Help")
      return(win)
    })
    rp.button(win, "Cancel", quitbutton=TRUE, pos=c(w/2, rowY, w/4, h) , action=function(win) {
      message("Aborting");
      return(win)
    })

    rp.button(win, "Run", pos=c(3*w/4, rowY, w/4, h), action=function(win) {

      # print(win$file)

      if (nTaxa <= nBoxes) {
        taxa <- names(win$taxa)[win$taxa]
      } else {
        taxa <- names(win$taxa1)[win$taxa1]
        taxa <- c(taxa, names(win$taxa2)[win$taxa2])
      }
      # print(taxa)

      if (length(taxa) < 2) {
        rp.messagebox("At least two taxa must be selected", title="Warning")
        return(win)   # return early
      }

      # simplify the taxa argument if all taxa are selected
      if (all(allTaxa %in% taxa)) {
        taxa <- ""
      }

      variables <- win$variables
      if (length(win$variables) < 2) {
        rp.messagebox("At least two variables must be selected", title="Warning")
        return(win)   # return early
      }
      # print(variables)

      # variables transformations
      if (file.exists(optionsFile)) {
        # read the options in the file
        load(optionsFile)
        transformations <- opts$transformations
      } else {
        transformations <- NULL
      }
      # print(transformations)

      print(win$n.groups)
      if (win$n.groups == 0) {
        n.groups <- NULL
      } else {
        n.groups <- win$n.groups
      }
      # print(n.groups)
      # print(win$min.n.groups)
      # print(win$max.n.groups)

      if (win$pre.sample) {
        pre.sample=2500
      } else {
        pre.sample=NULL
      }
      # print(pre.sample)

      if (win$intern.sample) {
        intern.sample=100000
      } else {
        intern.sample=NULL
      }
      # print(intern.sample)

      # build the function call
      message("Command:")
      call <- str_c("gdm(",
        "file=", str_c(deparse(win$file, width=500), collapse=""),
        ", taxa=", str_c(deparse(taxa, width=500), collapse=""),
        ", variables=", str_c(deparse(variables, width=500), collapse=""),
        ifelse( is.null(transformations),
          ", transformations=NULL",
          str_c(", transformations=", str_c("c(",str_c(names(transformations), str_c("\"", unlist(transformations),"\"") , sep="=", collapse=", "), ")"))
        ),
        ", lat.min=", win$lat.min, ", lat.max=", win$lat.max, ", lat.step=", win$lat.step,
        ", lon.min=", win$lon.min, ", lon.max=", win$lon.max, ", lon.step=", win$lon.step,
        # ", quick=", win$quick,
        # ", extrapolate.env=", win$extrapolate.env,
        # ", overlay.station=", win$overlay.stations,
        ", n.groups=", deparse(dput(n.groups)),
        ", min.n.groups=", win$min.n.groups, ", max.n.groups=", win$max.n.groups,
        ", taxa.min=", win$taxa.min,
        ", pre.sample=", pre.sample, ", intern.sample=", intern.sample,
        ", save=", win$save,
        ")"
      )
      cat(call, "\n")

      # execute call
      b <- eval(parse(text=call))

      return(win)

    })


    return(win)
  })

}

do.gdm.variables <- function(variables, file) {
  #
  # Second part of the GUI. Called from do.gdm()
  # Select variables transformations
  #

  suppressPackageStartupMessages(require("rpanel"))

  # defaults
  transformations <- rep("x",length(variables))
  names(transformations) <- variables

  # if the options file exists, read it to use the previously saved weights and transformations
  # TODO avoid using a file storage, we should be able to pass everything through functions and operate in memory
  if (file.exists(file)) {
    load(file)

    # transformations are stored as a list in the options file
    # transform them back to vectors
    opts <- lapply(opts, function(x) {
      y <- as.character(x)
      names(y) <- names(x)
      return(y)
    })

    # the selected variables may have changed, use what we can from the options file (opts$***) and use defaults for the rest (***)
    commonNames <- intersect(names(transformations), names(opts$transformations))
    transformations[match(commonNames, names(transformations))] <- opts$transformations[match(commonNames, names(opts$transformations))]
  }

  # default dimensions (in px)
  w <- 600          # width of the window
  h <- 50           # height of elements
  h.var <- 25       # height of a variable in the list
  spacer <- 10      # height of spacer


  n <- length(variables)
  varH <- (n+1) * h.var
  examplesH <- (3+1) * h.var
  windowH <-  varH + examplesH + spacer + h

  # main window
  win <- rp.control(title="Variables", size=c(w, windowH), aschar=F)

  # transformations
  # TODO look into why with rp.textentry.immediate, the results are not all carried to the stage of the close button
  rp.textentry(win, var=transformsBox, labels=variables, title="Transformations", initval=transformations, pos=c(0,0,w,varH), action=function(win) {
    return(win)
  })
  # TODO adapt height to the number of variables

  # provide example transformations
  example.transforms.labels=c("log10(x+1)","Square root","log10(-1*negative values only)")
  example.transforms.functions=list('"log10(x+1)"', '"sqrt(x)"', '"x[x>=0]=NA; log10(-x)"')
  rp.textentry(win, var=exampletransformBox, labels=example.transforms.labels, title="Example transformations", initval=example.transforms.functions, pos=c(0,varH,w,examplesH), action=function(win) {
    return(win)
  })

  # close button
  rp.button(win, title="Close", quitbutton=TRUE, pos=c(1/2*w, varH + examplesH + spacer, w/2, h), action=function(win) {
    # print("Transformations")
    # print(win$transformsBox)
    # print("Weights")
    # print(win$weightsBox)

    # extract transformations and weights and store them in named lists
    transformations <- as.list(win$transformsBox)

    # store them in a single object
    opts <- list(transformations=transformations)

    # save this object to the options file
    save(opts, file=file)

    return(win)
  })

  return(invisible(win))
}

