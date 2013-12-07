#
#         Shiny app server, which does most of the work
#
#   (c) Copyright 2013 Jean-Olivier Irisson
#         GNU General Public License v3
#
#--------------------------------------------------------------------------

dmess("run server.R")

library("plyr")

shinyServer(function(input, output) {

   dmess("run shinyServer")

   # UI functions
   #--------------------------------------------------------------------------
   # Create UI which depends on input values
   
   # list of species in input file
   output$speciesList <- renderUI({
      dmess("render species list")
      d <- read_input_data()
      sp <- setdiff(names(d), c("lat", "lon"))
      if ( ! is.null(sp) ) selectInput("species", "Species to model", sp)
   })

   # list of variables depending on options chosen above in the UI
   output$variablesList <- renderUI({
      dmess("render variables list")
      vars <- list.env.data(path="../../env_data")

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


   # Reactive functions
   #--------------------------------------------------------------------------
   # Their output is used by several other functions and depends on the input values

   # read user supplied data
   read_input_data <- reactive({
      dmess("run read_input_data")

      if (is.null(input$dataFile)) {

         # if no data file is uploaded, return NULL
         dmess("no data file selected")
         d <- NULL

      } else {

         # otherwise read the data and return the data.frame
         vmess("Read input data file ", input$dataFile$name)
         dataFile <- input$dataFile$datapath
         d <- read.data(dataFile, filetype="csv")
      }
      d
   })

   # read env data
   read_env_data <- reactive({
      dmess("run read_env_data")
      vmess("Read environmental database")

      # read selected variables from the database
      database <- read.env.data(input$vars, path="../../env_data", verbose=TRUE)
      # remove information on land
      database <- mask.env.data(database, path="../../env_data")
      
      database
   })

   # fit the model and compute predictions
   fit_model <- reactive({
      dmess("run fit_model")
      
      d <- read_input_data()
      # if there is no data, stop and give an error
      if ( is.null(d) ) {
         stop("Upload a CSV file to start")
      }
      database <- read_env_data()

      # if no explanatory variables are selected, stop and give an error
      if ( is.null(input$vars) ) {
         stop("Select at least one explanatory variable")
      }

      # weight observations within each cell of the environmental grid
      if ( input$bin ) {
         vmess("Weight data per bin")
         weights <- weight.per.bin(lon=d$lon, lat=d$lat, bin=0.1)         
      } else {
         weights <- rep(1, nrow(d))
      }

      # get environment data for the observations
      vmess("Associate env data with input data")
      x <- get.env.data(lon=d$lon, lat=d$lat, database)

      # remove points with only missing data
      onlyNA <- which(too.many.na(x, p=1))
      vmess("Removing ", length(onlyNA), " observations (over ", nrow(x), ") because of missing environmental data")
      x <- x[-onlyNA,]
      d <- d[-onlyNA,]
      weights <- weights[-onlyNA]

      # set up the model call
      vmess("Fitting BRT model for ", input$species)
      call <- str_c("brt.fit(x=x, y=d[,\"", input$species, "\"]",
                    # ", distribution=\"", input$distribution, "\"",
                    ", distribution=\"bernoulli\"",
                    ", n.trees=", input$n.trees,
                    ", interaction.depth=", input$interaction.depth,
                    ", shrinkage=", input$shrinkage,
                    ", bag.fraction=", input$bag.fraction,
                    ", max.cv.fold=", input$max.cv.fold,
                    # ", min.n.trees=", input$min.n.trees,
                    ", n.boot=", input$n.boot,
                    ", weights=weights",
                    ", verbose=", verbose,
                    ")"
      )
      dmess(call)
      m <- eval(parse(text=call))

      # return the object
      m
   })

   # generate prediction grid
   generate_pred_grid <- reactive({
      dmess("run generate_pred_grid")

      database <- read_env_data()
      m <- fit_model()
      variables <- m$var.names

      vmess("Generate prediction grid")
      predGrid <- build.grid(
         lat.min=input$lat[1], lat.max=input$lat[2], lat.step=input$latStep,
         lon.min=input$lon[1], lon.max=input$lon[2], lon.step=input$lonStep
      )
      xPred <- get.env.data(lon=predGrid$lon, lat=predGrid$lat, database)

      # remove data outside original data range
      if ( ! input$extrapolate ) {
         vmess("Remove data outside input environmental range")
         x <- data.frame(matrix(m$data$x, nrow=nrow(m$data$x.order)))
         names(x) <- variables
         ranges <- llply(x, range, na.rm=T)
         for (var in variables) {
            xPred[,var][which(xPred[,var] < ranges[[var]][1] | xPred[,var] > ranges[[var]][2]) ] <- NA
         }
      }

      # remove locations with only NA
      predData <- cbind(predGrid, xPred)
      onlyNA <- which(too.many.na(predData[,variables], p=1))
      if ( length(onlyNA) > 0 ) {
         vmess("Removing ", length(onlyNA), " points of the ", nrow(predData), " points prediction grid\n(on land or outside the range of the training data)")
         predData <- predData[-onlyNA,]
      }

      # remove locations with less that min.var.prop% of variance explained
      infl <- relative.influence(m, n.trees=m$best.iter)
      w <- infl / sum(infl)
      notEnoughData <- which(too.many.na(predData[,variables], p=input$min.var.prop/100, weights=w))
      if ( length(notEnoughData) > 0 ) {
         vmess("Further removing ", length(notEnoughData), " points of the ", nrow(predData), " points prediction grid\n(available environmental data does not allow to explain ", input$min.var.prop,"% of variance)")
         predData <- predData[-notEnoughData,]
      }
      
      predData
   })

   # predict distribution
   predict_distrib <- reactive({
      dmess("run predict_distrib") 

      m <- fit_model()
      predData <- generate_pred_grid()
      
      vmess("Predict distribution")
      pred <- predict(m, predData)
      # pred <- pred[,-which(names(pred) %in% m$var.names)]

      pred
   })


   # Output functions
   #--------------------------------------------------------------------------
   # Return text or plots to the user, destined to be displayed in the main panel

   # give a summary of the model
   output$modelSummary <- renderPrint({
      dmess("run modelSummary")

      if (input$run > 0) {
         isolate({
            m <- fit_model()
            
            vmess("Compute summary of model effects")
            summary.brt(m, plotit=FALSE)
         })
      }
   })

   # plot the effects
   output$modelPlot <- renderPlot({
      dmess("run modelPlot")
      if (input$run > 0) {
         isolate({
            m <- fit_model()
            
            vmess("Compute and plot model effects")
            plot.effects(m)
         })
      }
   })

   # plot the predictions
   output$predPlot <- renderPlot({
      dmess("run predPlot")
      if (input$run > 0) {
         isolate({
            if ( ! input$predict ) {
               stop("No predictions")
            } else {
               pred <- predict_distrib()
               
               vmess("Plot habitat suitability")
               print(plot.pred.brt(x=pred, quick=input$quick, scale=0.7))
            }
         })
      }
   })
 
   dlFile <- function(extension, suffix="") {
      dmess("generate filename for ", extension, " file")
      paste("brt-", format(Sys.time(), format="%Y%m%d_%H%M%S"), ".", extension, sep="")
   }
 
   output$downloadCSV <- downloadHandler(
      filename = dlFile(ext="csv"),
      content = function(file) {
         dmess("generate CSV file")         
         pred <- predict_distrib()
         write.csv(pred, file=file)
      },
      contentType = "text/csv"
   )
    
   output$downloadNetCDF <- downloadHandler(
      filename = dlFile(ext="nc"),
      content = function(file) {
         dmess("generate netCDF file")         
         pred <- predict_distrib()
         write.netcdf(pred, file=file, dimensions=c("lon","lat"))
      },
      contentType = "application/octet-stream"
   )

   output$downloadEffectsPlot <- downloadHandler(
      filename = dlFile(ext="pdf", suffix="effects"),
      content = function(file) {
         dmess("plot effects to PDF")
         m <- fit_model()
         pdf(file, width=8, height=6)
         plot.effects(m)
         dev.off()
      }
   )
   

})
