#
#     Customization of rpanel widgets
#
# (c) Copyright 2012 Jean-Olivier Irisson, Ben Raymond
#     GNU General Public License v3
#
#-----------------------------------------------------------------------------


rp.listbox.mult <- function (panel, var, vals, labels = vals, rows = length(vals), cols=NULL, initval = vals[1], parent = window, pos = NULL, title = deparse(substitute(var)), action = I, ...) {
  #
  # List box allowing multiple selection
  #

  suppressPackageStartupMessages(require("rpanel"))

  varname <- deparse(substitute(var))
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
  inittclvalue <- rpanel:::.rp.initialise(panelname, varname, initval = initval)
  if (rpanel:::.checklayout(pos)) {
    # if (is.null(pos$grid)) {
      gd = panel$window
    # } else {
    #   gd = rpanel:::.geval(panelname, "$", pos$grid)
    # }
    newlistbox <- tkwidget(gd, "labelframe", text = title)
    # if ((is.null(pos$row)) && (is.null(pos$column))) {
      rpanel:::.rp.layout(newlistbox, pos)
    # } else {
    #   if (is.null(pos$sticky)) {
    #     pos$sticky <- "w"
    #   }
    #   if (is.null(pos$rowspan)) {
    #     pos$rowspan = 1
    #   }
    #   if (is.null(pos$columnspan)) {
    #     pos$columnspan = 1
    #   }
    #   tkgrid(newlistbox, row = pos$row, column = pos$column,
    #     sticky = pos$sticky, `in` = gd, rowspan = pos$rowspan,
    #     columnspan = pos$columnspan
    #   )
    # }
    if (rows != length(vals)) {
      scr <- tkscrollbar(newlistbox, repeatinterval = 5, command = function(...) tkyview(listBox, ...))
      # if ((is.null(pos$width)) && (is.null(pos$height))) {
      if (is.null(cols)) {
        listBox <- tklistbox(newlistbox, height = rows, selectmode = "extended", yscrollcommand = function(...) tkset(scr, ...), background = "white")
      } else {
        listBox <- tklistbox(newlistbox, height = rows, width=cols, selectmode = "extended", yscrollcommand = function(...) tkset(scr, ...), background = "white")
      }
      # } else {
      #   listBox <- tklistbox(newlistbox, height = rows, selectmode = "extended", yscrollcommand = function(...) tkset(scr, ...), background = "white", width = pos$width, height = pos$height)
      # }
    } else {
      # if ((is.null(pos$width)) && (is.null(pos$height))) {
      if (is.null(cols)) {
        listBox <- tklistbox(newlistbox, height = rows, selectmode = "extended", background = "white")
      } else {
        listBox <- tklistbox(newlistbox, height = rows, width=cols, selectmode = "extended", background = "white")
      }
      # } else {
      #   listBox <- tklistbox(newlistbox, height = rows, selectmode = "extended", background = "white", width = pos$width, height = pos$height)
      # }
    }
    if (rows != length(vals)) {
      tkgrid(listBox, scr)
      tkgrid.configure(scr, rowspan = length(vals), sticky = "nsw")
    } else {
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
      } else {
        stop("The panel was not passed back from the action function.")
      }
    })
  }
  if (ischar) {
    invisible(panelname)
  } else {
    assign(panelreturn, rpanel:::.geval(panelname), envir = parent.frame())
  }
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
  if (ischar) {
    invisible(panelname)
  } else {
    assign(panelreturn, rpanel:::.geval(panelname), envir = parent.frame())
  }
}

rp.textentry.immediate <- function (panel, var, action = I, labels = NA, names = labels, title = NA, initval = NA, parent = window, pos = NULL, ...) {
  #
  # Same as rp.textentry, but the "var" entries are modified immediately (as the user enters them)
  #  rather than waiting until the user presses enter
  #
  # BR April 2012

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
      } else {
        nboxes <- length(labels)
        if (is.na(title) & (nboxes == 1))
          title <- labels
      }
      initval <- rep(NA, nboxes)
    } else {
      nboxes <- length(initval)
      if ((length(labels) == 1) && (is.na(labels)))
        if (nboxes != 1) {
          labels <- paste(varname, 1:nboxes, sep = "")
        }
        else {
          labels <- varname
        } else if (length(labels) != nboxes)
        stop("lengths of labels and initval do not match.")
    }
    if ((nboxes == 1) & (!is.na(title)))
      labels <- title
    rpanel:::.geval(panelname, "$", varname, " <- vector(length=", nboxes, ")")
    if (nboxes > 1) {
      if (is.na(title)) {
        title <- varname
      }
    }
    if ((!is.list(pos)) || (is.null(pos$grid))) {
      gd = panel$window
    } else {
      gd = rpanel:::.geval(panelname, "$", pos$grid)
    }
    if (nboxes > 1) {
      frame <- tkwidget(gd, "labelframe", text = title, padx = 2, pady = 2)
    } else {
      frame <- tkframe(gd)
    }
    if ((!is.list(pos)) || ((is.null(pos$row)) && (is.null(pos$column)))) {
      rpanel:::.rp.layout(frame, pos)
    } else {
      if (is.null(pos$sticky)) {
        pos$sticky <- "w"
      }
      if (is.null(pos$rowspan)) {
        pos$rowspan = 1
      }
      if (is.null(pos$columnspan)) {
        pos$columnspan = 1
      }
      tkgrid(frame, row = pos$row, column = pos$column, sticky = pos$sticky, `in` = gd, rowspan = pos$rowspan, columnspan = pos$columnspan)
    }
    for (i in 1:nboxes) {
      if (is.na(initval[i])) {
        initval[i] <- "NA"
      }
      inittclvalue <- rpanel:::.rp.initialise(panelname, paste(varname, i, sep = ""), initval[i])
      tclvariable <- rpanel:::.geval(panelname, "$", varname, i, ".tcl <- tclVar(", deparse(inittclvalue), ")")
      if (is.numeric(inittclvalue)) {
        rpanel:::.geval(panelname, "$", varname, "[", i, "] <- deparse(", inittclvalue, ")")
      } else {
        rpanel:::.geval(panelname, "$", varname, "[", i, "] <- '", inittclvalue, "'")
      }
      if (!any(is.na(names))) {
        rpanel:::.geval("names(", panelname, "$", varname, ")[", i, "] <- '", names[i], "'")
      }
      f <- function() {
        for (i in 1:nboxes) {
          rpanel:::.geval(panelname, "$", varname, "[", i, "] <- tclvalue(", panelname, "$", varname, i, ".tcl)")
          if (!any(is.na(names))) {
            rpanel:::.geval("names(", panelname, "$", varname, ")[", i, "] <- '", names[i], "'")
          }
        }
        panel <- action(rpanel:::.geval(panelname))
        if (!is.null(panel$intname)) {
          rpanel:::.gassign(panel, panelname)
        } else {
          stop("The panel was not passed back from the action function.")
        }
      }
      if (!any(is.na(names))) {
        rpanel:::.geval("names(", panelname, "$", varname, ")[", i, "] <- '", names[i], "'")
      }
      label <- tklabel(frame, text = labels[i], height = "1")
      if ((!is.list(pos)) || ((is.null(pos$width)) && (is.null(pos$height)))) {
        entry <- tkentry(frame, textvariable = tclvariable)
      } else {
        entry <- tkentry(frame, textvariable = tclvariable, width = pos$width)
      }
      tkgrid(label, entry)
      tkgrid.configure(label, sticky = "w")
      tkgrid.configure(entry, sticky = "e")
      tkbind(entry, "<Enter>", f)
    }
  }
  if (ischar) {
    invisible(panelname)
  } else {
    assign(panelreturn, rpanel:::.geval(panelname), envir = parent.frame())
  }
}

