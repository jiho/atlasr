#
#     Customization of rpanel widgets
#
# (c) Copyright 2012 J-O Irisson
#                    Ben Raymond
#     http://creativecommons.org/licenses/by/3.0/
#
#-----------------------------------------------------------------------------


w.listbox.mult <- function(parent, title = NA, labels, rows = length(labels), initval = labels[1], action=I, pos=NULL, foreground=NULL, background="white", font=NULL, sleep = 0.01) {
  if (is.na(title))
     widget <- rpanel:::w.createwidget(parent, pos, background)
  else
     widget <- rpanel:::w.createwidget(parent, pos, background, title)
  widget$.type <- "listbox"

  if (rows < length(labels)) {
     widget$.widget <- rpanel:::handshake(tklistbox, parent$.handle, height=rows, selectmode = "extended",
                                 yscrollcommand = function(...) rpanel:::handshake(tkset, scr, ...))
     scr <- rpanel:::handshake(tkscrollbar, parent$.handle, repeatinterval = 5,
                      command=function(...) rpanel:::handshake(tkyview, widget$.widget,...))
     rpanel:::w.appearancewidget(widget, font, foreground, background, scr)
  }
  else {
     widget$.widget <- rpanel:::handshake(tklistbox, parent$.handle, height=rows, selectmode = "extended")
     rpanel:::w.appearancewidget(widget, font, foreground, background)
  }

  selection <- 1

#  tkinsert(widget$.widget, "end", "test")

  for (i in (1:length(labels))) {
    Sys.sleep(sleep)
    rpanel:::handshake(tkinsert, widget$.widget, "end", as.character(labels[[i]]))
    if (labels[[i]] == initval) selection <- i - 1
  }

  rpanel:::handshake(tkselection.set, widget$.widget, selection)

  f <- function(...) action(labels[as.numeric(rpanel:::handshake(tkcurselection, widget$.widget)) + 1])

  rpanel:::handshake(tkbind, widget$.widget, "<ButtonRelease-1>", f)

  invisible(widget)
}


rp.listbox.mult <- function(panel, variable, vals, labels = vals, rows = length(labels), initval = vals[1], pos = NULL, title = deparse(substitute(variable)), action = I,  foreground = NULL, background = NULL, font = NULL, parentname = deparse(substitute(panel)), sleep = 0.01, name = paste("listbox", rpanel:::.nc(), sep = ""), ...) {

  if (!exists(panel$panelname, rpanel:::.rpenv, inherits = FALSE)) { # if the panelname is not set then
     panelname <- deparse(substitute(panel)) # the panel name should be the panel deparse subst'ed
     # 13/03/2012 these lines are not commented out in previous version
     #    panel <- rp.control.get(panelname, panel) # now get the panel
     #    panel$panelname = panelname # now set the panelname properly
     #    assign(panelname, panel, envir=.rpenv) # now send back the panel
  }
  else
    panelname <- panel$panelname
    # 13/03/2012 these lines are not commented out in previous version
    #    panel <- rp.control.get(panelname, panel) # now get the panel

  varname <- deparse(substitute(variable))
  if (!rpanel:::rp.isnull(panelname, varname))
     variable = rpanel:::rp.var.get(panelname, varname)
  else
     variable = initval; rpanel:::rp.var.put(panelname, varname, variable)

  if (is.null(pos) && length(list(...)) > 0) pos <- list(...)

  f <- function(val) {
     rpanel:::rp.var.put(panelname, varname, val)
     panel <- rpanel:::rp.control.get(panelname)
     panel <- action(panel)
     rpanel:::rp.control.put(panelname, panel)
  }

  if (rpanel:::rp.widget.exists(panelname, parentname))
     parent <- rpanel:::rp.widget.get(panelname, parentname)
  else
     parent <- panel
  if (is.list(pos) && !is.null(pos$grid)) parent <- rp.widget.get(panelname, pos$grid)

  widget <- w.listbox.mult(parent, title, labels, rows, initval = variable, action=f, pos,
                      foreground, background, font, sleep)
  rpanel:::rp.widget.put(panelname, name, widget)
  if (rpanel:::.rpenv$savepanel) rp.control.put(panelname, panel) # put the panel back into the environment
  invisible(panelname)
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
      # initial value in rp.textentry
      # tkbind(entry, "<Key-Return>", f)
      # make the entry immediate
      tkbind(entry, "<Enter>", f)
    }
  }
  if (ischar) {
    invisible(panelname)
  } else {
    assign(panelreturn, rpanel:::.geval(panelname), envir = parent.frame())
  }
}

