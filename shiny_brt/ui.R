#
#      Shiny server GUI, which presents controls and results
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------


dmess("run ui.R")

library("shiny")
library("shinyIncubator")   # for action button


shinyUI(pageWithSidebar(

  # Title
  headerPanel("Model data with Boosted Regression Trees"),


  # User input
  sidebarPanel(
    h5("Input data"),

    # input data file
    helpText("Upload a csv file with columns for latitude and longitude and each species of interest"),
    fileInput("dataFile", label="", multiple=FALSE, accept="text/csv"),

    # List of species in the loaded data file
    conditionalPanel(
      condition="typeof(input.dataFile) != 'undefined'",
      uiOutput("speciesList")
    ),
    # TODO try to hide this with conditional panel based on input.dataFile
    #      this would allow to avoid the conditions which check wether the data in NULL in the corresponding server code

    h5("Options"),

    # Settings for the run
    # selectInput("distrib", "Distribution", list(
    #   "bernouilli",
    #   "gaussian",
    #   "poisson"
    # )),
    # helpText("When bernouilli is selected, the data will be converted to presence and absence only, even if it contains abundances"),

    # checkboxInput("bin", "Bin input data on 0.1º grid", FALSE),
    # helpText("All data points will be kept but points in the same grid cell will be weighted down"),

    # checkboxInput("boot", "Bootstrap model", value=FALSE),

    checkboxInput("optim", "Choose shrinkage (i.e. learning rate) and number of trees automatically", FALSE),
    conditionalPanel(
      condition="input.optim == false",
      wellPanel(
        helpText("Ideally, you want a large number of trees with a very small learning rate"),
        sliderInput("n.trees", "Maximum number of trees", min=500, max=10000, step=500, value=5000),
        sliderInput("shrinkage", "Shrinkage (i.e. learning rate) per tree", min=0.001, max=0.01, step=0.001, value=0.005)
      )
    ),

    checkboxInput("advanced", "Advanced gbm settings", value=FALSE),
    conditionalPanel(
      condition="input.advanced == true",
      wellPanel(
        helpText("See ?gbm in R"),
        sliderInput("interaction.depth", "interaction.depth", min=1, max=10, step=1, value=2),
        sliderInput("bag.fraction", "bag.fraction", min=0, max=1, step=0.1, value=0.5),
        sliderInput("cv.fold", "cv.fold", min=0, max=20, step=1, value=5)
      )
    ),


    checkboxInput("predict", "Predict distribution", value=FALSE),

    # Domain for prediction
    conditionalPanel(
      condition="input.predict == true",
      wellPanel(
        checkboxInput("quick", "Quick, unprojected, prediction plot", TRUE),

        checkboxInput("extrapolate", "Predict beyond observed environmental range"),
        # helpText("not recommended"),

        # checkboxInput("overlay", "Overlay stations on prediction map"),

        # checkboxInput("bootpred", "Bootstrap predictions"),

        sliderInput("lat", "Latitudinal limits and step", min=-80, max=-30, step=1, value=c(-80,-30)),
        sliderInput("latStep", "", min=0.1, max=4, step=0.1, value=1),
        sliderInput("lon", "Longitudinal limits and step", min=-180, max=180, step=5, value=c(-180,180)),
        sliderInput("lonStep", "", min=0.1, max=4, step=0.1, value=1)
      )
    ),

    # Environmental variables
    h5("Explanatory variables"),
    checkboxInput("interpolated", "Prefer interpolated variables", TRUE),
    checkboxInput("distance", "Allow 'distance from ...' variables", FALSE),
    checkboxGroupInput("depths", "Restrict to depths:", c(0,50,200,250,500)),
    checkboxGroupInput("season", "Restrict to season:", c("summer", "winter"), "summer"),
    uiOutput("variablesList"),
    # checkboxGroupInput("vars", "", allVariables),

    # Run button
    # submitButton("Run model")
    div(
      style="text-align: right",
      actionButton("run", "Run !")
    )

  ),


  # Render the result
  mainPanel(
    tabsetPanel(
      tabPanel("Summary", verbatimTextOutput("modelSummary")),
      tabPanel("Effects", plotOutput("modelPlot")),
      tabPanel("Prediction", plotOutput("predPlot"))
    )
    # h5("Model summary"),
    # verbatimTextOutput("modelSummary"),
    # h5("Model effects"),
    # plotOutput("modelPlot"),
    # h5("Model prediction"),
    # plotOutput("predPlot")
  )

))
