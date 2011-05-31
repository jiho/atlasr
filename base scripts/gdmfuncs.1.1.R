"gdm.fit" <- 
function (x, y, geo=FALSE, wtype=c("equal", "standard", "custom"),w=NULL, sample=NULL) 
{
    options(warn.FPU = FALSE)
    if ( nrow(x) != nrow(y) )
	stop( "x and y have different row counts" )

    if ( missing(wtype) )
        wtype <- "equal"

    switch(wtype,
           equal=weights <- rep(1,times=nrow(x)),
           standard=weights <- apply(y, 1, sum),
           custom=weights <- w)

    if ( is.null(weights) )
	stop( "ERROR: NULL weights here, have you supplied a vector for argument w?" )

    if ( geo == TRUE ) {
	  DoGeo <- 1
	  nPredSplines <- ( ncol(x) - 1 ) * 3
	  predlist <- c("Geographic", names(x)[3:length(names(x))])
    }
    else {
        DoGeo <- 0
        nPredSplines <- ncol(x) * 3
        predlist <- names(x)
    }

    npairs <- nrow(x) * (nrow(x)-1) / 2
    if ( is.null(sample) ) {        
        samplevec <- 0
    }
    else {
        if ( sample >= npairs ) {  # this is really an argument error but we'll fix it
	    samplevec <- 0         # so treat as not sub-sampled...
        }
        else {
            totalpairs <- npairs
	    npairs <- sample
            samplevec <- sort(sample(seq(from=0,to=totalpairs,by=1),npairs))
        }
    }

    dyn.load("gdmlib.dll")
    p1 <- 0
    p2 <- 0
    p3 <- 0
    p4 <- 0
    p5 <- rep(0,times=nPredSplines)
    p6 <- rep(0,times=nPredSplines)
    z <- .C( "DoRStatModel", 
             gdmdev = as.double(p1),
	     nulldev = as.double(p2),
	     expdev = as.double(p3),
             intercept = as.double(p4),		
             coeffs = as.double(p5),
	     quants = as.double(p6),
             as.double(as.matrix(t(x))), 
             as.integer(as.matrix(t(y))), 
             as.integer(nrow(x)), 
             as.integer(ncol(x)), 
	     as.integer(ncol(y)), 
             as.double(weights),
             as.integer(DoGeo),
             as.integer(npairs),
             as.integer(samplevec))
    dyn.unload("gdmlib.dll")
    call <- match.call()
    m <- match.call(expand = F)
    if ( wtype == "equal" ) {
	weightReturn <- "equal"
    }
    else if ( wtype == "standard" ) {
	weightReturn <- "standard"
    }
    else {
	weightReturn <- m[[5]];
    }
	
    return(structure(list(xdataname = m[[2]],
                          ydataname = m[[3]],
                          geo = geo,
                          sample = npairs,
                          gdmdeviance = z$gdmd,
                          nulldeviance = z$nulldev,
                          explained = z$expdev,
                          intercept = z$intercept,
                          predictors = predlist,
                          coefficients = z$coeffs,
                          quantiles = z$quants,
                          weights = weightReturn,
                          creationdate=date())))
}



"gdm.plot" <-
function (model, xdata, ydata, plot.layout = c(1,1), plot.color = rgb(0.5,0.5,0.5), plot.linewidth=2.0) 
{
    options(warn.FPU = FALSE)
    PSAMPLE <- 200
    dyn.load("gdmlib.dll")
    pObserved <- rep(0,times=model$sample)
    pPredicted <- rep(0,times=model$sample)
    nPreds <- length(model$predictors)
    pTmpMatrix <- matrix(0,PSAMPLE,nPreds)
    for ( i in 1:nPreds ) {
        pTmpMatrix[,i] <- seq( from=model$quantiles[[(i*3)-2]],to=model$quantiles[[i*3]], length=PSAMPLE )
    }
    z <- .C( "DoRStatModelPlot", 
              observed = as.double(pObserved),
              predicted = as.double(pPredicted),
              as.double(as.matrix(t(xdata))), 
              as.integer(as.matrix(t(ydata))),
              tmpmat = as.double(as.matrix(pTmpMatrix)),  
              as.double(model$intercept),
	      as.double(model$coefficients),
              as.double(model$quantiles),	
              as.integer(nrow(xdata)), 
              as.integer(ncol(xdata)), 
	      as.integer(ncol(ydata)),
              as.integer(model$geo)
              )
    dyn.unload("gdmlib.dll")
    pn <- model$predictors
    s1 <- pn[[1]]
    nthis <- 0
    for ( i in 1:nPreds ) {            
        if ( sum(model$coefficients[((i*3)-2):(i*3)]) > 0 ) {
	    if ( nthis == 0 ) 
                s1 <- pn[[i]]
            else
	        s1 <- paste(s1,"+",pn[[i]])
	    nthis <- nthis + 1
	}
    }      	
              
    totalplots <- plot.layout[1] * plot.layout[2]
    if ( totalplots == 1 ) {
        singlepage <- T
        par(cex = 1.0)
    }
    else {
        singlepage <- F
        par(cex = 0.5)
    }

    ## apply the link function and plot.....
    plot(z$predicted, z$observed, 
         xlab="Predicted Ecological Distance", 
         ylab="Observed Compositional Dissimilarity", type="n" )
    points( z$predicted, z$observed, pch=20, cex=0.25, col=plot.color )
    overlayX <- seq( from=min(z$predicted), to=max(z$predicted), length=PSAMPLE )
    overlayY <- 1 - exp( - overlayX )
    lines( overlayX, overlayY, lwd=plot.linewidth ) #title(sub=s1)

    ## use the raw data and plot.....
    x11()
    dev.next()
    plot(1.0 - exp(-1.0*z$predicted), z$observed, 
         xlab="Predicted Compositional Dissimilarity", 
         ylab="Observed Compositional Dissimilarity", type="n" )
    points( 1.0 - exp(-1.0*z$predicted), z$observed, pch=20, cex=0.25, col=plot.color )
    overlayX <- overlayY <- seq( from=min(1.0 - exp(-1.0*z$predicted)), to=max(1.0 - exp(-1.0*z$predicted)), length=PSAMPLE )
    lines( overlayX, overlayY, lwd=plot.linewidth ) #title(sub=s1)

    ## only plot the predictors with non-zero coeffiecients.....
    par(cex = 0.5)
    par( mfcol=plot.layout )
    finalMatrix <- matrix(0,PSAMPLE,nPreds)
    pmin <- 1
    pmax <- PSAMPLE
    nstart <- 0
    nthis <- nstart
    totalmin <- min(z$tmpmat)
    totalmax <- max(z$tmpmat)
    for ( i in 1:nPreds ) {            
        finalMatrix[,i] <- z$tmpmat[pmin:pmax]
        pmin <- pmin + PSAMPLE
        pmax <- pmax + PSAMPLE
        ## only if the sum of the coefficients associated with this predictor is > 0.....
	if ( sum(model$coefficients[((i*3)-2):(i*3)]) > 0 ) {
            if ( singlepage ) {
	        x11()
      	        dev.next()
                plot( seq(from=model$quantiles[[(i*3)-2]],to=model$quantiles[[(i*3)]], length=PSAMPLE),
                      finalMatrix[,i], xlab=pn[i], ylab=paste("f(", pn[i], ")", sep="" ), ylim=c(totalmin,totalmax), type="l" )
	    }
	    else {
                plot( seq(from=model$quantiles[[(i*3)-2]],to=model$quantiles[[(i*3)]], length=PSAMPLE),
                      finalMatrix[,i], xlab=pn[i], ylab=paste("f(", pn[i], ")", sep="" ), ylim=c(totalmin,totalmax), type="l" )
		nthis <- nthis + 1
		if ( nthis >= totalplots ) {
                    if ( i == nPreds )
                        break				    
                    else {  
                        ## create a new page..... 
                        x11()
                        dev.next()
                        par( mfcol=plot.layout )
                        nthis <- nstart
                    }                               
                }
            }
        }
    }
}




"gdm.predict" <-
function (model, matrixA, matrixB) 
{
    options(warn.FPU = FALSE)
    PSAMPLE <- 200
   
    ## sanity check of matrix dimensions
    if ( nrow(matrixA) != nrow(matrixB) ) stop( "matrixA and matrixB have different row counts!" )
    if ( ncol(matrixA) != ncol(matrixB) ) stop( "matrixA and matrixB have different column counts!" )

    pPredicted <- rep(0,times=nrow(matrixA))
    dyn.load("gdmlib.dll")
    z <- .C( "DoRStatModelPredict", 
             predicted = as.double(pPredicted),
             as.double(as.matrix(t(matrixA))), 
             as.double(as.matrix(t(matrixB))),
             as.double(model$intercept),
             as.double(model$coefficients),
             as.double(model$quantiles),	
             as.integer(nrow(matrixA)), 
             as.integer(ncol(matrixB)), 
             as.integer(model$geo)
             )
    dyn.unload("gdmlib.dll")

    ## apply the link function before returning
    return( 1 - exp( - z$predicted ) )
}



"gdm.summary" <-
function (model) 
{
    print( "", quote=F )	
    print( "", quote=F )	
    print( "GDM Modelling Summary", quote=F );
    print( paste( "Creation Date: ", model$creationdate ), quote=F );
    print( "", quote=F )	
    call <- match.call()
    m <- match.call(expand = F)
    print( paste( "Name: ", m[[2]] ), quote=F )
    print( "", quote=F )	
    print( paste( "X Data: ", model$xdataname ), quote=F )
    print( "", quote=F )	
    print( paste( "Y Data: ", model$ydataname ), quote=F )
    print( "", quote=F )
    print( paste( "Weighting: ", model$weights ), quote=F )
    print( "", quote=F )	
    print( paste( "Samples: ", model$sample ), quote=F )
    print( "", quote=F )	
    print( paste( "Use Geographical Distance: ", model$geo ), quote=F )
    print( "", quote=F )	
    print( paste( "NULL Deviance: ", model$nulldeviance ), quote=F )
    print( paste( "GDM Deviance: ", model$gdmdeviance ), quote=F )	
    print( paste( "Deviance Explained: ", model$explained ), quote=F )
    print( "", quote=F )	
    print( paste( "Intercept: ", model$intercept ), quote=F )
    print( "", quote=F )	
    thiscoeff <- 1
    thisquant <- 1
    for ( i in 1:length(model$predictors) ) {
	print( paste( "Predictor ",i,": ",model$predictors[[i]], sep="" ), quote=F )            
        for ( j in 1:3 ) {
	    if ( j == 1 ) print( paste( "Min Quantile: ",model$quantiles[[thisquant]], sep="" ), quote=F )
	    else if ( j == 2 ) print( paste( "50% Quantile: ",model$quantiles[[thisquant]], sep="" ), quote=F )
	    else if ( j == 3 ) print( paste( "Max Quantile: ",model$quantiles[[thisquant]], sep="" ), quote=F )
	    thisquant <- thisquant + 1
	}
	for ( j in 1:3 ) {
	    print( paste( "Coefficient[",j,"]: ",model$coefficients[[thiscoeff]], sep="" ), quote=F )
	    thiscoeff <- thiscoeff + 1
	}
	print( "", quote=F )	            
    }	
}



"gdm.transform" <-
function (model, xdata) 
{
    options(warn.FPU = FALSE)
    PSAMPLE <- nrow(xdata)
    dyn.load("gdmlib.dll")
    nPreds <- length(model$predictors)
    pTmpMatrix <- matrix(0,PSAMPLE,nPreds)
##    for ( i in 1:nPreds ) {
##        pTmpMatrix[,i] <- seq( from=model$quantiles[[(i*3)-2]],to=model$quantiles[[i*3]], length=PSAMPLE )
##    }
    z <- .C( "DoRStatModelTransform", 
             as.double(as.matrix(xdata)), 
             tmpmat = as.double(as.matrix(pTmpMatrix)),  
             as.double(model$intercept),
	       as.double(model$coefficients),
             as.double(model$quantiles),	
             as.integer(nrow(xdata)), 
             as.integer(ncol(xdata)), 
	       as.integer(model$geo)
             )
    dyn.unload("gdmlib.dll")
    pn <- model$predictors
      
    finalMatrix <- as.data.frame(matrix(0,PSAMPLE,nPreds))
    pmin <- 1
    pmax <- PSAMPLE
    for ( i in 1:nPreds ) {            
        finalMatrix[,i] <- z$tmpmat[pmin:pmax]
        pmin <- pmin + PSAMPLE
        pmax <- pmax + PSAMPLE
        names(finalMatrix)[[i]] <- pn[i]
    }
    return(finalMatrix)
}


