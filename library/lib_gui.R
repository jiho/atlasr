#
#     Custom GUIs to run high level functions
#
#     - modification of rpanel widgets
#     - BRT GUI
#
# (c) Copyright 2012 Jean-Olivier Irisson, Ben Raymond
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------


rp.listbox.mult <- function (panel, var, vals, labels = vals, rows = length(vals), cols=NULL,
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
            if (is.null(cols)) {
                listBox <- tklistbox(newlistbox, height = rows,
                  selectmode = "extended", yscrollcommand = function(...) tkset(scr,
                    ...), background = "white")
            } else {
                listBox <- tklistbox(newlistbox, height = rows, width=cols,
                  selectmode = "extended", yscrollcommand = function(...) tkset(scr,
                    ...), background = "white")
            }
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
            if (is.null(cols)) {
                listBox <- tklistbox(newlistbox, height = rows,
                  selectmode = "extended", background = "white")
            } else {
                listBox <- tklistbox(newlistbox, height = rows, width=cols,
                  selectmode = "extended", background = "white")
            }
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


rp.text <- function (panel, txt="", parent = window, pos = NULL, ...) {
  #
  # Simple text-displaying panel
  #
  # BR April 2012

  suppressPackageStartupMessages(require("rpanel"))

  ischar <- is.character(panel)
  if (ischar) {
    panelname <- panel
    panel <- rpanel:::.geval(panel)
  } else {
    panelname <- panel$intname
    panelreturn <- deparse(substitute(panel))
    rpanel:::.gassign(panel, panelname)
  }
  pos = rpanel:::.newpos(pos, ...)
  f <- function(...) {
    panel <- action(rpanel:::.geval(panelname))
    if (!is.null(panel$intname)) {
      rpanel:::.gassign(panel, panelname)
    } else {
      stop("The panel was not passed back from the action function.")
    }
  }
  if (rpanel:::.checklayout(pos)) {
    if ((!is.list(pos)) || (is.null(pos$row))) {
      newbutton <- tktext(panel$window, borderwidth=0, fg="black", bg="grey85", wrap="char")
      tkinsert(newbutton,"end",txt)
      rpanel:::.rp.layout(newbutton, pos)
    } else {
      if (is.null(pos$grid)) {
        gd = panel$window
      } else {
        gd = rpanel:::.geval(panelname, "$", pos$grid)
      }
      if ((is.null(pos$width)) && (is.null(pos$height))) {
        newbutton <- tktext(panel$window)
      } else {
        newbutton <- tktext(panel$window, width = pos$width, height = pos$height)
      }
      tkinsert(newbutton,"end",txt)
      if (is.null(pos$sticky)) {
        pos$sticky <- "w"
      }
      if (is.null(pos$rowspan)) {
        pos$rowspan = 1
      }
      if (is.null(pos$columnspan)) {
        pos$columnspan = 1
      }
      tkgrid(newbutton, row = pos$row, column = pos$column,
             sticky = pos$sticky, `in` = gd, rowspan = pos$rowspan,
             columnspan = pos$columnspan
      )
    }
  }
  if (ischar)
      invisible(panelname)
  else assign(panelreturn, rpanel:::.geval(panelname), envir = parent.frame())
}

rp.textentry.immediate <- function (panel, var, action = I, labels = NA, names = labels,
    title = NA, initval = NA, parent = window, pos = NULL, ...)
{
    ## Same as rp.textentry, but the "var" entries are modified immediately (as the user enters them)
    ##  rather than waiting until the user presses enter
    ##
    ## BR April 2012

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
    if (rpanel:::.checklayout(pos)) {
        varname <- deparse(substitute(var))
        if ((varname %in% names(panel)) && (all(is.na(initval)))) {
            initval <- rpanel:::.geval(panelname, "$", varname)
        }
        if ((length(initval) == 1) && (is.na(initval))) {
            if ((length(labels) == 1) && (is.na(labels))) {
                nboxes <- 1
                if (is.na(title))
                  title <- varname
                labels <- varname
            }
            else {
                nboxes <- length(labels)
                if (is.na(title) & (nboxes == 1))
                  title <- labels
            }
            initval <- rep(NA, nboxes)
        }
        else {
            nboxes <- length(initval)
            if ((length(labels) == 1) && (is.na(labels)))
                if (nboxes != 1) {
                  labels <- paste(varname, 1:nboxes, sep = "")
                }
                else {
                  labels <- varname
                }
            else if (length(labels) != nboxes)
                stop("lengths of labels and initval do not match.")
        }
        if ((nboxes == 1) & (!is.na(title)))
            labels <- title
        rpanel:::.geval(panelname, "$", varname, " <- vector(length=",
            nboxes, ")")
        if (nboxes > 1)
            if (is.na(title))
                title <- varname
        if ((!is.list(pos)) || (is.null(pos$grid))) {
            gd = panel$window
        }
        else {
            gd = rpanel:::.geval(panelname, "$", pos$grid)
        }
        if (nboxes > 1) {
            frame <- tkwidget(gd, "labelframe", text = title,
                padx = 2, pady = 2)
        }
        else {
            frame <- tkframe(gd)
        }
        if ((!is.list(pos)) || ((is.null(pos$row)) && (is.null(pos$column)))) {
            rpanel:::.rp.layout(frame, pos)
        }
        else {
            if (is.null(pos$sticky)) {
                pos$sticky <- "w"
            }
            if (is.null(pos$rowspan)) {
                pos$rowspan = 1
            }
            if (is.null(pos$columnspan)) {
                pos$columnspan = 1
            }
            tkgrid(frame, row = pos$row, column = pos$column,
                sticky = pos$sticky, `in` = gd, rowspan = pos$rowspan,
                columnspan = pos$columnspan)
        }
        for (i in 1:nboxes) {
            if (is.na(initval[i]))
                initval[i] <- "NA"
            inittclvalue <- rpanel:::.rp.initialise(panelname, paste(varname,
                i, sep = ""), initval[i])
            tclvariable <- rpanel:::.geval(panelname, "$", varname, i,
                ".tcl <- tclVar(", deparse(inittclvalue), ")")
            if (is.numeric(inittclvalue)) {
                rpanel:::.geval(panelname, "$", varname, "[", i, "] <- deparse(",
                  inittclvalue, ")")
            }
            else {
                rpanel:::.geval(panelname, "$", varname, "[", i, "] <- '",
                  inittclvalue, "'")
            }
            if (!any(is.na(names))) {
                rpanel:::.geval("names(", panelname, "$", varname, ")[",
                  i, "] <- '", names[i], "'")
            }
            f <- function() {
                for (i in 1:nboxes) {
                  rpanel:::.geval(panelname, "$", varname, "[", i, "] <- tclvalue(",
                    panelname, "$", varname, i, ".tcl)")
                  if (!any(is.na(names))) {
                    rpanel:::.geval("names(", panelname, "$", varname,
                      ")[", i, "] <- '", names[i], "'")
                  }
                }
                panel <- action(rpanel:::.geval(panelname))
                if (!is.null(panel$intname)) {
                  rpanel:::.gassign(panel, panelname)
                }
                else {
                  stop("The panel was not passed back from the action function.")
                }
            }
            if (!any(is.na(names))) {
                rpanel:::.geval("names(", panelname, "$", varname, ")[",
                  i, "] <- '", names[i], "'")
            }
            label <- tklabel(frame, text = labels[i], height = "1")
            if ((!is.list(pos)) || ((is.null(pos$width)) && (is.null(pos$height)))) {
                entry <- tkentry(frame, textvariable = tclvariable)
            }
            else {
                entry <- tkentry(frame, textvariable = tclvariable,
                  width = pos$width)
            }
            tkgrid(label, entry)
            tkgrid.configure(label, sticky = "w")
            tkgrid.configure(entry, sticky = "e")
            tkbind(entry, "<Enter>", f)
        }
    }
    if (ischar)
        invisible(panelname)
    else assign(panelreturn, rpanel:::.geval(panelname), envir = parent.frame())
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
