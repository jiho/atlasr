#
#     Custom GUIs to run high level functions
#
#     - modification of rpanel widgets
#     - BRT GUI
#
# (c) Copyright 2012 Jean-Olivier Irisson
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------


rp.listbox.mult <- function (panel, var, vals, labels = vals, rows = length(vals),
    initval = vals[1], parent = window, pos = NULL, title = deparse(substitute(var)),
    action = I, ...)
{
    #
    # List box allowing multiple selection
    #

    suppressPackageStartupMessages(require("rpanel"))

    varname <- deparse(substitute(var))
    ischar <- is.character(panel)
    if (ischar) {
        panelname <- panel
        panel <- rpanel:::.geval(panel)
    }
    else {
        panelname <- panel$intname
        panelreturn <- deparse(substitute(panel))
        rpanel:::.gassign(panel, panelname)
    }
    pos = rpanel:::.newpos(pos, ...)
    inittclvalue <- rpanel:::.rp.initialise(panelname, varname, initval = initval)
    if (rpanel:::.checklayout(pos)) {
        # if (is.null(pos$grid)) {
            gd = panel$window
        # }
        # else {
        #     gd = rpanel:::.geval(panelname, "$", pos$grid)
        # }
        newlistbox <- tkwidget(gd, "labelframe", text = title)
        # if ((is.null(pos$row)) && (is.null(pos$column))) {
            rpanel:::.rp.layout(newlistbox, pos)
        # }
        # else {
        #     if (is.null(pos$sticky)) {
        #         pos$sticky <- "w"
        #     }
        #     if (is.null(pos$rowspan)) {
        #         pos$rowspan = 1
        #     }
        #     if (is.null(pos$columnspan)) {
        #         pos$columnspan = 1
        #     }
        #     tkgrid(newlistbox, row = pos$row, column = pos$column,
        #         sticky = pos$sticky, `in` = gd, rowspan = pos$rowspan,
        #         columnspan = pos$columnspan)
        # }
        if (rows != length(vals)) {
            scr <- tkscrollbar(newlistbox, repeatinterval = 5,
                command = function(...) tkyview(listBox, ...))
            # if ((is.null(pos$width)) && (is.null(pos$height))) {
                listBox <- tklistbox(newlistbox, height = rows,
                  selectmode = "extended", yscrollcommand = function(...) tkset(scr,
                    ...), background = "white")
            # }
            # else {
            #     listBox <- tklistbox(newlistbox, height = rows,
            #       selectmode = "extended", yscrollcommand = function(...) tkset(scr,
            #         ...), background = "white", width = pos$width,
            #       height = pos$height)
            # }
        }
        else {
            # if ((is.null(pos$width)) && (is.null(pos$height))) {
                listBox <- tklistbox(newlistbox, height = rows,
                  selectmode = "extended", background = "white")
            # }
            # else {
            #     listBox <- tklistbox(newlistbox, height = rows,
            #       selectmode = "extended", background = "white",
            #       width = pos$width, height = pos$height)
            # }
        }
        if (rows != length(vals)) {
            tkgrid(listBox, scr)
            tkgrid.configure(scr, rowspan = length(vals), sticky = "nsw")
        }
        else {
            tkgrid(listBox)
        }
        selected <- 0
        for (i in (1:length(vals))) {
            tkinsert(listBox, "end", labels[i])
            if (inittclvalue == vals[i])
                selected <- i
        }
        tkselection.set(listBox, selected - 1)
        tkbind(listBox, "<ButtonRelease-1>", function(...) {
            rpanel:::.geval(panelname, "$", varname, " <- c('", paste(vals[as.numeric(tkcurselection(listBox)) + 1], collapse="','"),"')")
            panel <- action(rpanel:::.geval(panelname))
            if (!is.null(panel$intname)) {
                rpanel:::.gassign(panel, panelname)
            }
            else {
                stop("The panel was not passed back from the action function.")
            }
        })
    }
    if (ischar)
        invisible(panelname)
    else assign(panelreturn, rpanel:::.geval(panelname), envir = parent.frame())
}

do.brt <- function(...) {
  #
  # Open a GUI to select the arguments of the brts() function
  #
  # ...   passed to brts()
  #

  suppressPackageStartupMessages(require("rpanel"))


  # default dimensions (in px)
  w <- 600      # width of the window
  h <- 50       # height of elements
  spacer <- 10  # height of spacer

  # main window
  win <- rp.control(title="Run BRT model", size=c(w,700))


  rp.button(win, title="Choose file", pos=c(0, 0, w/4, h), action=function(win) {

    # choose the data file
    file <- file.choose()
    # file <- "Austropallene.csv"
    file <- win2unix(file, ...)


    # write the filename, for information
    rp.textentry(win, file, title="Dataset", initval=file, pos=c(w/4, h/4, 3*w/4, h))

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
      rp.checkbox(win, var=taxa, labels=allTaxa, title="Taxa", initval=c(TRUE, rep(FALSE, times=nTaxa-1)), pos=c(0, h, 2*w/3, hSel), action=function(win) {
        if (all(! win$taxa)) {
          rp.messagebox("At least one taxon must be selected", title="Warning")
        }
        return(win)
      })

    } else {
      # otherwise, split in half
      half <- ceiling(nTaxa/2)
      rest <- nTaxa - half
      allTaxa1 <- allTaxa[1:half]
      allTaxa2 <- allTaxa[(half+1):nTaxa]

      rp.checkbox(win, var=taxa1, labels=allTaxa1, title="Taxa", initval=c(TRUE, rep(FALSE, times=half-1)), pos=c(0, h, w/3, hSel), action=function(win) {
        if (all(! c(win$taxa1, win$taxa2))) {
          rp.messagebox("At least one taxon must be selected", title="Warning")
        }
        return(win)
      })
      rp.checkbox(win, var=taxa2, labels=allTaxa2, title="", initval=rep(FALSE, times=rest), pos=c(w/3, h+7, w/3, hSel-7), action=function(win) {
        if (all(! c(win$taxa1, win$taxa2))) {
          rp.messagebox("At least one taxon must be selected", title="Warning")
        }
        return(win)
      })

    }

    # build environment variables list
    allVariables <- list.env.data()
    rp.listbox.mult(win, var=variables, vals=allVariables, title="Variables",  rows=nRows, pos=c(2*w/3, h, w/3, hSel-h), action=function(win) {
      if (all(win$variables == "")) {
        rp.messagebox("At least one variable must be selected", title="Warning")
      }
      return(win)
    })
    rp.checkbox(win, var=selAllVariables, title="Select all", initval=FALSE, pos=c(2*w/3, hSel, w/3, h))


    # compute vertical coordinate of the middle of the window
    mid <- hSel + h + spacer


    # effects options
    rp.radiogroup(win, family, values=c("bernoulli", "gaussian", "poisson"), title="Effects", pos=c(0, mid, w/4, h*2))

    # prediction options
    rp.radiogroup(win, prediction, values=c("no", "yes", "yes + bootstrap"), initval="yes", title="Prediction", pos=c(0, mid+h*2, w/4, h*2))

    # checkboxes range
    checkH <- 4*h/5
    rp.checkbox(win, bootstrap.effects, title="Bootstrap effects", initialval=FALSE, pos=c(w/4, mid, w/4, checkH))
    rp.checkbox(win, bin, title="Bin original data\non prediction grid", initval=FALSE, pos=c(w/4, mid+checkH, w/4, checkH))
    rp.checkbox(win, extrapolate.env, title="Extrapolate envi-\nronmental range", initval=FALSE, pos=c(w/4, mid+checkH*2, w/4, checkH))
    rp.checkbox(win, quick.plot, title="Subsample predict-\nion plot (faster)", initialval=FALSE, pos=c(w/4, mid+checkH*3, w/4, h))
    rp.checkbox(win, overlay.stations, title="Overlay stations\non prediction plot", initialval=FALSE, pos=c(w/4, mid+checkH*4, w/4, checkH))

    # location
    rp.slider(win, lat.max,  from=-90, to=-30,  resolution=1,   title="North"   , initval=-30 , showvalue=TRUE, pos=c(w/2+w/8, mid    , w/4, h))
    rp.slider(win, lon.min,  from=-180, to=180, resolution=1,   title="West"    , initval=-180, showvalue=TRUE, pos=c(w/2, mid+h*1, w/4, h))
    rp.slider(win, lon.max,  from=-180, to=180, resolution=1,   title="East"    , initval=180 , showvalue=TRUE, pos=c(3*w/4, mid+h*1, w/4, h))
    rp.slider(win, lat.min,  from=-90, to=-30,  resolution=1,   title="South"   , initval=-80 , showvalue=TRUE, pos=c(w/2+w/8, mid+h*2, w/4, h))
    rp.slider(win, lon.step, from=0.1, to=4  ,  resolution=0.1, title="Step lon", initval=0.5 , showvalue=TRUE, pos=c(w/2, mid+h*3, w/4, h))
    rp.slider(win, lat.step, from=0.1, to=4   , resolution=0.1, title="Step lat", initval=0.1 , showvalue=TRUE, pos=c(3*w/4, mid+h*3, w/4, h))


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
    rp.button(win, "Run", pos=c(3*w/4, rowY, w/4, h), action=function(win, ...) {
      # print(win$file)

      if (nTaxa <= nBoxes) {
        taxa <- names(win$taxa)[win$taxa]
      } else {
        taxa <- names(win$taxa1)[win$taxa1]
        taxa <- c(taxa, names(win$taxa2)[win$taxa2])
      }
      # print(taxa)

      # print(win$variables)
      # print(win$selAllVariables)
      if (win$selAllVariables) {
        variables <- allVariables
      } else {
        variables <- win$variables
      }

      # print(win$family)
      # print(win$bootstrap.effects)
      if(win$bootstrap.effects) {
        n.boot.effects <- 200
      } else {
        n.boot.effects <- 0
      }

      # print(win$prediction)
      if (win$prediction == "no") {
        predict <- FALSE
        n.boot.pred <- 0
      } else if (win$prediction == "yes") {
        predict <- TRUE
        n.boot.pred <- 0
      } else if (win$prediction == "bootstrap") {
        predict <- TRUE
        n.boot.pred <- 200
      }

      # print(win$lat.max)
      # print(win$lat.min)
      # print(win$lat.step)
      # print(win$lon.max)
      # print(win$lon.min)
      # print(win$lon.step)

      # print(win$quick.plot)
      # print(win$extrapolate.env)
      # print(win$overlay.stations)
      # print(win$bin)


      b <- brts(
                  file=win$file,
                  taxa=taxa, variables=variables,
                  lat.min=win$lat.min, lat.max=win$lat.max, lat.step=win$lat.step,
                  lon.min=win$lon.min, lon.max=win$lon.max, lon.step=win$lon.step,
                  predict=predict,
                  bin=win$bin,
                  #
                  family=win$family,
                  n.boot.effects=n.boot.effects,
                  n.boot.pred=n.boot.pred,
                  quick=win$quick.plot,
                  extrapolate.env=win$extrapolate.env,
                  overlay.station=win$overlay.stations,
                  ...
              )
      return(win)
    }, ...)


    return(win)
  })

}


# Experimentation with gwidgets
# library("gWidgets")
# options("guiToolkit"="tcltk")
# # options("guiToolkit"="RGtk2")
#
#
# defaultMessage="Choose a dataset file first"
#
# # Window
# win <- gwindow(title="Run a BRT model", visible=FALSE, width=500, height=720)
#
# # Layout
# # global container
# globalG <- ggroup(horizontal=FALSE, container=win)
#
# tab <- glayout(container=globalG, spacing=5)
#
# datasetNameW <- glabel(defaultMessage, container=tab, fill="x")
# datasetSelectW <-  gbutton("Choose file", container=tab, fill="both", expand=T, handler=function(h, ...) {
#   # datasetName <- file.choose()
#   datasetName <- "Austropallene.csv"
#
#   # read the data and detect species lists
#   observed_data <- read.data(datasetName)
#   species <- setdiff(names(observed_data), c("lat", "lon"))
#
#   # write the filename, for information
#   svalue(datasetNameW) <- datasetName
#
#   # update the choice of species
#   speciesW[] <- data.frame(Species=species)
# })
# tab[1,1] <- datasetSelectW
# tab[1,2:4] <- datasetNameW
#
# # species choice
# speciesW <- gtable(defaultMessage, multiple=TRUE, container=tab)
# tab[2,1:2] <- speciesW
#
# # variables choice
# envVariables <- list.env.data()
# envVariables <- data.frame(Variables=envVariables)
# variablesW <- gtable(envVariables, multiple=TRUE, container=tab)
# allVariablesW <- gcheckbox("Select All", checked=FALSE, container=tab, handler=function(h, ...) {
#   if (svalue(h$obj)) {
#     # cat("Select all variables\n")
#     variablesW[] <- data.frame(Variables = "All")
#   } else {
#     variablesW[] <- envVariables
#   }
# })
# tab[2,3:4] <- variablesW
# tab[3,3:4] <- allVariablesW
#
# effectsG <- gframe("Effects", container=tab, horiz=F, expand=T)
# familyW <- gcombobox(c("Bernoulli", "Gaussian", "Poisson"), container=effectsG)
# bootEffectsW <- gcheckbox("Bootstrap", checker=FALSE, container=effectsG)
# tab[4,1:2] <- effectsG
#
# locationG <- gframe("Location", container=tab, horiz=F, expand=T)
# locs <- glayout(container=locationG)
#
# locs[1,2] <- "North"
# locs[2,2] <- gspinbutton(from=-90, to=0, by=1, value=-30, container=locs)
#
# locs[2,1] <- "West"
# locs[3,1] <- gspinbutton(from=-180, to=180, by=1, value=-180, container=locs)
#
# locs[3,2] <- "South"
# locs[4,2] <- gspinbutton(from=-90, to=0, by=1, value=-80, container=locs)
#
# locs[2,3] <- "East"
# locs[3,3] <- gspinbutton(from=-180, to=180, by=1, value=180, container=locs)
#
# tab[4,3:4] <- locationG
#
# visible(win) <- TRUE
#
#
#
# buttonCancelW <- gbutton("Cancel", container=buttonsG, handler <- function(h, ...) {
#   message("OK, aborting ...")
#   dispose(h$obj)
# })
# buttonOKW <- gbutton("Run", container=buttonsG, handler <- function(h, ...) {
#   selectedSpecies <- svalue(speciesW)
#   message(selectedSpecies)
#   dispose(h$obj)
# })
#
#