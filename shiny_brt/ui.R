#
#         Shiny server GUI, which presents controls and results
#
#   (c) Copyright 2013 Jean-Olivier Irisson
#         GNU General Public License v3
#
#--------------------------------------------------------------------------


dmess("run ui.R")

library("shiny")

shinyUI(pageWithSidebar(

   # Title
   headerPanel("Model data with Boosted Regression Trees"),


   # User input
   sidebarPanel(

      div( id="input",
         h5("Input data"),
                  
         # upload data file
         helpText("Upload a csv file with columns for latitude and longitude and each species of interest"),
         fileInput("dataFile", label="", multiple=FALSE, accept="text/csv"),
   
         # list species in the loaded data file
         # conditionalPanel(
         #    condition="typeof(input.dataFile) != 'undefined'",
            uiOutput("speciesList"),
         # ),
   
         checkboxInput("bin", "Bin input data on 0.1º grid", TRUE),
         helpText("All data points will be kept but points in the same grid cell will be weighted down")
      ),

      div( id="options",
         h5("Options - see ?gbm in R"),
      
         # selectInput("distribution", "Distribution", list(
         #    "bernoulli",
         #    "gaussian",
         #    "poisson"
         # )),
         # helpText("When bernoulli is selected, the data will be converted to presence and absence only, even if it contains abundances"),

         # checkboxInput("optim", "Choose shrinkage and number of trees automatically", FALSE),
         # conditionalPanel(
         #    condition="input.optim == false",
         #    wellPanel(
         #       helpText("Ideally, you want a large number of trees with a very small learning rate"),
               sliderInput("n.trees", "Maximum number of trees (the actual number will be estimated through cross validation)", min=500, max=10000, step=500, value=2000),
               sliderInput("shrinkage", "Shrinkage per tree", min=0.001, max=0.05, step=0.001, value=0.01),
         #    )
         # ),
         # conditionalPanel(
         #    condition="input.optim == true",
         #    wellPanel(
         #       sliderInput("min.n.trees", "Minimum number of trees", min=500, max=10000, step=500, value=2000),
         #       sliderInput("shrinkage", "Starting value for the shrinkage per tree", min=0.001, max=0.1, step=0.001, value=0.05)
         #    )
         # ),
      
         checkboxInput("advanced", "Advanced gbm settings", value=FALSE),
         conditionalPanel(
            condition="input.advanced == true",
            wellPanel(
               sliderInput("interaction.depth", "interaction.depth", min=1, max=10, step=1, value=3),
               sliderInput("bag.fraction", "bag.fraction", min=0, max=1, step=0.1, value=0.5),
               sliderInput("max.cv.fold", "cv.fold", min=0, max=20, step=2, value=6)
            )
         ),

         sliderInput("n.boot", "Number of bootstraps", min=0, max=500, step=50, value=0)
      ),

      div( id="prediction",
         h5("Prediction"),
         checkboxInput("predict", "Predict distribution", value=FALSE),

         # Domain for prediction
         conditionalPanel(
            condition="input.predict == true",
            wellPanel(
               checkboxInput("quick", "Quick, unprojected, prediction plot", FALSE),

               checkboxInput("extrapolate", "Predict beyond observed environmental range", FALSE),
               sliderInput("min.var.prop", "Minimum variance for prediction", min=0, max=100, step=10, value=50),
               helpText("When some environmental data is missing, prediction will only be made if the available environmental data allows to capture the specified percentage of the variance"),

               # checkboxInput("overlay", "Overlay stations on prediction map"),

               sliderInput("lat", "Latitudinal limits and step", min=-80, max=-30, step=1, value=c(-80,-30)),
               sliderInput("latStep", "", min=0.1, max=4, step=0.1, value=2),
               sliderInput("lon", "Longitudinal limits and step", min=-180, max=180, step=5, value=c(-180,180)),
               sliderInput("lonStep", "", min=0.1, max=4, step=0.1, value=2)
            )
         )
      ),

      div( id="predictors",
         h5("Explanatory variables"),

         checkboxInput("interpolated", "Prefer interpolated variables", TRUE),
         checkboxInput("distance", "Allow 'distance from ...' variables", FALSE),
         checkboxGroupInput("depths", "Restrict to depths:", c(0,50,200,250,500)),
         checkboxGroupInput("season", "Restrict to season:", c("summer", "winter"), "summer"),

         uiOutput("variablesList")
      ),

      # Run button
      div(
         style="text-align: right",
         actionButton("run", "Run !")
         # NB: actionButton tracks the number of time it is pressed
         #     this allows to not try plotting the model on the first run
      )
   ),


   # Render the result
   mainPanel(

      h5("Model summary"),
      verbatimTextOutput("modelSummary"),

      h5("Model effects"),
      plotOutput("modelPlot"),
      # TODO dynamic height based on number of variables
      downloadButton("downloadEffectsPlot", "Download plot as PDF"),


      h5("Model prediction"),
      plotOutput("predPlot"),
      # TODO dynamic height depending on wether we do a quick or not quick plot
      downloadButton('downloadNetCDF', 'Download predictions as netCDF'),
      downloadButton('downloadCSV', 'Download predictions as CSV')
   )

))
