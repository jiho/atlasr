#
#      Shiny app server, which does most of the work
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

dmess("run server.R")

library("gbm")
library("plyr")
library("stringr")
library("ggplot2")

source("../library/lib_data.R")
source("../library/lib_plot.R")

shinyServer(function(input, output) {

  dmess("execute shinyServer")

  # Reactive functions
  # Their output is used by several other functions and depends on the input values

  # read user supplied data
  get.data <- reactive({
    dmess("run get.data")

    if (is.null(input$dataFile)) {

      # if no data file is uploaded, return NULL
      dmess("no data file selected")
      d <- NULL

    } else {

      # otherwise read the data and return the data.frame
      dmess("read data file ", input$dataFile$name)
      dataFile <- input$dataFile$datapath
      d <- read.data(dataFile, filetype="csv")
    }
    d
  })

  # fit the model and compute predictions
  fit.model <- reactive({
    dmess("fit model")

    # read data
    d <- get.data()

    # if there is no data, stop and give an error
    if (is.null(d)) {
      stop("Upload a CSV file to start")
    }

    # if no explanatory variables are selected, stop and give an error
    dmess("test if options are correctly set")
    if (is.null(input$vars)) {
      stop("Select at least one explanatory variable")
    }

    # interpolate environmental data at data points
    dmess("read and mask env data")
    env <- read.env.data(input$vars, verbose=FALSE)
    env <- mask.env.data(env)

    dmess("associate env data with input data")
    d <- associate.env.data(d, env)

    # set up the model call
    dmess("set up model")
    formula <- str_c(input$species, " ~ ", str_c(input$vars, collapse=" + "))
    dmess(formula)

    call <- str_c("gbm(", formula, ", data=d, distribution='bernoulli', n.trees=", input$n.trees, ", interaction.depth=", input$interaction.depth, ", shrinkage=", input$shrinkage, ", bag.fraction=", input$bag.fraction, ", cv.fold=", input$cv.fold, ", verbose=FALSE)")
    dmess(call)
    # TODO optimize number of trees

    m <- eval(parse(text=call))

    # find best number of trees
    m$best.iter <- gbm.perf(m, method="cv", plot.it=FALSE)

    # predict
    if (input$predict) {
      # prepare prediction grid
      lon <- seq(input$lon[1], input$lon[2], by=input$lonStep)
      lat <- seq(input$lat[1], input$lat[2], by=input$latStep)
      pred.d <- expand.grid(lon=lon, lat=lat)
      pred.d <- associate.env.data(pred.d, env)

      # predict probability of presence using the optimal number of trees
      pred.d$proba <- predict(m, newdata=pred.d, ntrees=m$best.iter, type="response")

      # store this in the object
      m$prediction <- pred.d
    }

    # return the object
    m
  })


  # UI functions
  # Create UI which depends on input values

  # list of species in input file
  output$speciesList <- renderUI({
    dmess("render species list")
    d <- get.data()
    sp <- names(d)[-c(1,2)]
    if (!is.null(sp)) selectInput("species", "Species to model", sp)
  })

  # list of variables depending on options chosen above in the UI
  output$variablesList <- renderUI({
    dmess("render variables list")
    vars <- list.env.data()

    # selected interpolated or non-interpolated variables
    interpolatedVersion <- vars[str_detect(vars, "interpolated")]
    nonInterpolatedVersion <- str_replace(interpolatedVersion, "_interpolated", "")
    neverInterpolated <- vars[! vars %in% c(interpolatedVersion, nonInterpolatedVersion)]
    if ( input$interpolated ) {
      vars <- vars[vars %in% c(neverInterpolated, interpolatedVersion)]
    } else {
      vars <- vars[vars %in% c(neverInterpolated, nonInterpolatedVersion)]
    }

    # use distance based variables
    if ( ! input$distance ) {
      vars <- vars[!str_detect(vars, "distance_")]
    }

    # select depths of interest
    depths <- str_replace(str_extract(vars, "_[0-9]+"), "_", "")
    hasDepth <- !is.na(depths)
    if (length(input$depths) > 0) {
      # detect the presence of any of the selected depths
      hasSelectedDepth <- depths %in% input$depths

      # keep depth-independent variables and depth-dependant variables at selected depths
      vars <- vars[ !hasDepth | hasSelectedDepth ]
    }

    # select season of interest
    if (length(input$season) == 1) {
      hasSeason <- str_detect(vars, "summer|winter")

      # detect the presence of the selected season
      hasSelectedSeason <- str_detect(vars, input$season)

      # keep season-independent variables and season-dependant variables for the selected season
      vars <- vars[ !hasSeason | hasSelectedSeason ]
    }

    checkboxGroupInput("vars", "Select variables:", vars)
  })


  # Output functions
  # Return text or plots to the user, destined to be displayed in the main panel

  # give a summary of the model
  output$modelSummary <- renderPrint({
    if (input$run > 0) {
      isolate({
        m <- fit.model()
        print(m)
        sum <- summary(m, plotit=FALSE)
        row.names(sum) <- NULL
        print(sum)
      })
    }
  })

  # plot the effects
  output$modelPlot <- renderPlot({
    if (input$run > 0) {
      isolate({
        m <- fit.model()

        n <- length(m$var.names)
        nrows <- ceiling(sqrt(n))
        ncols <- ceiling(n/nrows)
        par(mfrow=c(nrows, ncols))
        for (i in 1:n) {
          plot(m, i.var=i, n.trees=m$best.iter)
        }

        # # extract marginal effects for partial dependance plots
        # pdp <- adply(m$var.names, 1, function(name) {
        #   d <- plot(m, i.var=name, return.grid=TRUE)
        #   names(d) <- c("value", "marginal effect")
        #   return(d)
        # })
        #
        # # get quantiles of the data
        #
        #
        # # cleanup variable names
        # names(pdp)[1] <- "variable"
        # pdp$variable <- m$var.names[pdp$variable]
        # pdp$variable <- str_replace_all(pdp$variable, "_", " ")
        #
        # # plot effects
        # library("ggplot2")
        # ggplot(pdp) + geom_path(aes(x=value, y=`marginal effect`)) + facet_wrap(~variable, scales="free_x")
      })
    }
  })

  # plot the predictions
  output$predPlot <- renderPlot({
    if (input$run > 0) {
      isolate({

        m <- fit.model()

        source("library/lib_plot.R")
        if (input$quick) {
          geom="raster"
        } else {
          geom="tile"
        }

        print(polar.ggplot(m$prediction, aes(fill=proba), geom=geom) + scale_fill_gradientn(colours=continuous.colourmap(), limits=c(0,1)))
      })
    }
  })

})
