#############################################################################
TDboost.profile.real <- function(formula, p.vec=NULL, xi.vec=NULL, link.power = 0, 
                            data, weights, offset, fit.glm=FALSE, 
                            do.smooth=TRUE, do.plot=FALSE, 
                            do.ci=do.smooth, eps=1/6,
                            control=list( epsilon=1e-09, maxit=glm.control()$maxit, trace=glm.control()$trace ),
                            do.points=do.plot, method="inversion", conf.level=0.95, 
                            phi.method=ifelse(method=="saddlepoint","saddlepoint","mle"), verbose=FALSE, add0=FALSE,verbose1=FALSE,...) {
  # verbose gives feedback on screen:
  #    0 : minimal (FALSE)
  #    1 : small amount (TRUE)
  #    2 : lots
  # Determine the value of  p  to use to fit a Tweedie distribution
  # A profile likelihood is used (can can be plotted with do.plot=TRUE)
  
  # The plot can be (somewhat) recovered using (where  out  is the
  # returned object, with smoothing):
  #     plot( out$x, out$y)
  # The actual points can then be added:
  #     points( out$p, out$L)
  # (Note that out$L and out$y are the same if smoothing is  not requested.  
  # Likewise for out$p and out$x.)
  
  # Peter Dunn
  # 07 Dec 2004
  
  
  if ( is.logical( verbose ) ) {
    verbose <- as.numeric(verbose)
  }
  
  if (verbose >= 1 ) {
    cat("---\n This function may take some time to complete;\n")
    cat(" Please be patient.  If it fails, try using  method=\"series\"\n")
    cat(" rather than the default  method=\"inversion\"\n")
    cat(" Another possible reason for failure is the range of p:\n")
    cat(" Try a different input for  p.vec\n---\n")
  }
  
  
  cl <- match.call()
  mf <- match.call()
  m <- match(c("formula", "data", "weights","offset"), names(mf), 0L)
  mf <- mf[c(1, m)]
  formula_input <- mf[[2]]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  if (!is.null(offset)) {
    if (length(offset) == 1)
      offset <- rep(offset, NROW(Y))
    else if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  
  
  
  ### NOW some notation stuff
  xi.notation <- TRUE
  if ( is.null(xi.vec) & !is.null(p.vec) ){ # If p.vec given, but not xi.vec
    xi.vec <- p.vec
    xi.notation <- FALSE
  }
  
  if ( is.null(p.vec) & !is.null(xi.vec) ){ # If xi.vec given, but not p.vec
    p.vec <- xi.vec
    xi.notation <- TRUE
  }
  
  if ( is.null( p.vec ) & is.null(xi.vec)) { # Neither xi.vec or p.vec are given
    if ( any(Y==0) ){
      p.vec <- seq(1.2, 1.8, by=0.1)
    } else {
      p.vec <- seq(1.5, 5, by=0.5)
    }
    xi.vec <- p.vec
    xi.notation <- TRUE
  }
  
  ### AT THIS POINT, we have both xi.vec and p.vec declared, and both are the same
  ### but we stick with using xi.vec hereafter
  
  
  # Determine notation to use in output (consistent with what was supplied by the user)
  index.par <- ifelse( xi.notation, "xi","p")
  
  # Find if values of p/xi have been specified that are not appropriate
  xi.fix <- any( xi.vec<=1 & xi.vec!=0, na.rm=TRUE)
  
  if ( xi.fix ) {
    xi.fix.these <- (xi.vec<=1 & xi.vec!=0)
    xi.vec[ xi.fix.these ] <- NA 
    if ( length(xi.vec) == 0 ) {
      stop(paste("No values of",index.par,"between 0 and 1, or below 0, are possible.\n"))
    } else {
      cat("Values of",index.par,"between 0 and 1 and less than zero have been removed: such values are not possible.\n")
    }
  }
  
  # Warnings
  if ( any( Y == 0 ) & any( (xi.vec >= 2) | (xi.vec <=0) ) ) {
    xi.fix.these <- (xi.vec>=2 | xi.vec<=0)
    xi.vec[ xi.fix.these ] <- NA 
    if ( length(xi.vec) == 0 ) {
      stop(paste("When the response variale contains exact zeros, all values of",index.par,"must be between 1 and 2.\n"))
    } else {
      cat("When the response variale contains exact zeros, all values of",index.par,"must be between 1 and 2; other values have been removed.\n")
    }
  }
  
  
  # Remove any NAs in the p/xi values
  xi.vec <- xi.vec[ !is.na(xi.vec) ] 
  
  if ( do.smooth & (length(xi.vec) < 5) ) {
    warning(paste("Smoothing needs at least five values of",index.par,".") )
    do.smooth <- FALSE
    do.ci <- FALSE
  }
  if ( (conf.level >= 1) | (conf.level <=0)  ){
    stop("Confidence level must be between 0 and 1.")
  }
  
  if ( !do.smooth & do.ci ) {
    do.ci <- FALSE
    warning("Confidence intervals only computed if  do.smooth=TRUE\n")
  }
  
  if ( add0 ) {
    p.vec  <- c( 0, p.vec )
    xi.vec <- c( 0, xi.vec )
  }
  
  
  cat(xi.vec,"\n")
  
  
  
  
  # Some renaming
  ydata <- Y
  model.x <- X
  
  
  
  # Now, fit models!  Need the Tweedie class
  
  # First define the function to *minimize*;
  # since we want a *maximum* likeihood estimator,
  # define the negative.
  dtweedie.nlogl <- function(phi, y, mu, power) {
    ans <- - sum( log( dtweedie( y=y, mu=mu, phi=phi, power=power ) ) )
    if ( is.infinite( ans ) ) {
      # If infinite, replace with saddlepoint estimate?
      ans <- sum( tweedie.dev(y=y, mu=mu, power=power) )/length( y )
    }
    attr(ans, "gradient") <- dtweedie.dldphi(y=y, mu=mu, phi=phi, power=power)
    
    ans
  }
  
  
  # Set up some parameters
  xi.len <- length(xi.vec)
  phi <- NaN
  
  L <- array( dim = xi.len )
  phi.vec <- L
  
  # Now most of these are for debugging:
  b.vec <- L
  c.vec <- L
  mu.vec <- L
  b.mat <- array( dim=c(xi.len, length(ydata) ) )
  
  for (i in (1:xi.len)) {
    
    if ( verbose>0) {
      cat( paste(index.par," = ",xi.vec[i],"\n", sep="") )
    } else {
      cat(".")
    }
    
    
    # We find the mle of phi.
    # We try to bound it as well as we can, 
    # then use  uniroot to solve for it.
    
    # Set things up
    p <- xi.vec[i]
    
    phi.pos <- 1.0e2
    bnd.neg <- -Inf
    bnd.pos <- Inf
    
    if (verbose==2) cat("* Fitting initial model:")
    
    # Fit the model with given p
    catch.possible.error <- try(
      fit.model <- TDboost(formula_input,
                           data = mf,
                           distribution = list(name="EDM",alpha=p),
                           shrinkage = 0.005,               
                           interaction.depth = 7,
                           bag.fraction = 0.5,
                           n.minobsinnode = 20,
                           keep.data = TRUE,
                           cv.folds = 5,
                           verbose = verbose1,...),
      
      silent = TRUE
    )
    
    
    skip.obs <- FALSE
    if ( class( catch.possible.error )=="try-error" ) {
      skip.obs <- TRUE 
    }
    
    if( skip.obs ) {
      warning(paste("  Problem near ",index.par," = ",p,"; this error reported:\n     ",
                    catch.possible.error,
                    " Examine the data and function inputs carefully.") )
      
      # NOTE:  We need epsilon to be very small to make a
      #        smooth likelihood plot
      mu <- rep(NA, length(ydata) )
      
    } else {
      best.iter <- TDboost.perf(fit.model,method="cv", plot.it = FALSE)
      print(best.iter)	
      tmp <- predict.TDboost(fit.model,mf, n.trees = best.iter)      
      mu <- tmp # be careful prediction is for mu directly
    }
    if (verbose==2) cat(" Done\n")
    
    ### -- Start:  Estimate phi
    if (verbose>=1) cat("* Phi estimation, method: ", phi.method)
    
    if( skip.obs ) {
      if (verbose>=1) cat("; but skipped for this obs\n")
      phi.est <- phi <- phi.vec[i] <- NA
    } else {
      if ( phi.method=="mle"){
        
        if (verbose>=1) cat(" (using optimize): ")
        
        # Saddlepoint approx of phi:
        phi.saddle <- sum( tweedie.dev(y=ydata, mu=mu, power=p) )/length( ydata )
        
        if ( is.nan(phi) ) {
          
          # NOTE:  This is only used for first p.
          #        The remainder use the previous estimate of phi as a starting point
          
          # Starting point for phi:  use saddlepoint approx/EQL.
          phi.est <- phi.saddle
          
        } else {
          phi.est <- phi
        }
        low.limit <- min( 0.001, phi.saddle/2)
        
        #    ans <- nlm(p=phi.est, f=dtweedie.nlogl, 
        #                hessian=FALSE,
        #                power=p, mu=mu, y=data)
        if ( p!= 0 ) {
          ans <- optimize(f=dtweedie.nlogl, maximum=FALSE, interval=c(low.limit, 10*phi.est),
                          power=p, mu=mu, y=ydata )
          phi <- phi.vec[i] <- ans$minimum
        } else {
          phi <- phi.vec[i] <- sum( (ydata-mu)^2 ) / length(ydata)
        }
        if (verbose>=1) cat(" Done (phi =",phi.vec[i],")\n")
        
      } else{ # phi.method=="saddlepoint")
        
        if (verbose>=1) cat(" (using mean deviance/saddlepoint): ")
        phi <- phi.est <- phi.vec[i] <- sum( tweedie.dev(y=ydata, mu=mu, power=p) )/length( ydata )
        if (verbose>=1) cat(" Done (phi =",phi,")\n")
      }
    }
    
    ### -- Done:  estimate phi
    
    # Now determine log-likelihood at this p
    # Best to use the same density evaluation for all:  series or inversion
    
    if (verbose>=1) cat("* Computing the log-likelihood ")
    
    # Now compute the log-likelihood
    if (verbose>=1) cat("(method =",method,"):")
    
    if ( skip.obs ) {
      if (verbose>=1) cat(" but skipped for this obs\n")
      L[i] <- NA
    } else {   
      if ( method=="saddlepoint") {
        L[i] <- dtweedie.logl.saddle(y=ydata, mu=mu, power=p, phi=phi, eps=eps)
      } else {
        if (p==2) {
          L[i] <- sum( log( dgamma( rate=1/(phi*mu), shape=1/phi, x=ydata ) ) )
        } else {
          if ( p == 1 ) {
            if ( phi==1 ){
              # The phi==1 part added 30 Sep 2010.
              L[i] <- sum( log( dpois(x=ydata/phi, lambda=mu/phi ) ) )
            } else { # p=1, but phi is not equal to zero
              # As best as we can, we determine if the values
              # of  ydata  are multiples of  phi
              y.on.phi <- ydata/phi
              close.enough <- array( dim=length(y.on.phi))
              for (i in (1:length(y.on.phi))){
                if (isTRUE(all.equal(y.on.phi, as.integer(y.on.phi)))){
                  L[i] <- sum( log( dpois(x=round(y/phi), lambda=mu/phi ) ) )
                } else {
                  L[i] <- 0
                }
              }
            }
          } else {
            if ( p == 0 ) {
              L[i] <- sum( dnorm(x=ydata, mean=mu, sd=sqrt(phi), log=TRUE) )
            } else {
              L[i] <- switch(
                pmatch(method, c("interpolation","series", "inversion"), 
                       nomatch=2),
                "1" = dtweedie.logl( mu=mu, power=p, phi=phi, y=ydata),
                "2" = sum( log( dtweedie.series( y=ydata, mu=mu, power=p, phi=phi) ) ),
                "3" = sum( log( dtweedie.inversion( y=ydata, mu=mu, power=p, phi=phi) ) ) )
            }
          }
        }
      }
    }
    
    if (verbose>=1) {
      cat(" Done: L =",L[i],"\n")
    } 
  }
  
  if ( verbose == 0 ) cat("Done.\n")
  ### Smoothing
  # If there are infs etc in the log-likelihood,
  # the smooth can't happen.  Perhaps this can be worked around.
  # But at present, we note when it happens to ensure no future errors
  # (eg in computing confidence interval for p).
  
  # y and x are the smoothed plot produced; here we set them
  # as NA so they are defined when the function is returned to
  # the workspace
  y <- NA
  x <- NA
  
  if ( do.smooth ) {
    L.fix <- L
    xi.vec.fix <- xi.vec
    phi.vec.fix <- phi.vec
    if ( any( is.nan(L) ) | any( is.infinite(L) ) | any( is.na(L) ) ) {
      retain.these <- !( ( is.nan(L) ) | ( is.infinite(L) ) | ( is.na(L) ) )
      xi.vec.fix <- xi.vec.fix[ retain.these ]
      phi.vec.fix <- phi.vec.fix[ retain.these ]
      L.fix <- L.fix[ retain.these ]
      
      # 	    p.vec.fix <- p.vec.fix[ !is.infinite(L.fix) ]
      # 	    phi.vec.fix <- phi.vec.fix[ !is.infinite(L.fix) ]
      # 	    L.fix <- L.fix[ !is.infinite(L.fix) ]
      
      if (verbose>=1) cat("Smooth perhaps inaccurate--log-likelihood contains  Inf  or  NA.\n")
    }
    #else {
    
    if ( length( L.fix ) > 0 ) {
      if (verbose>=1) cat(".")
      if (verbose>=1) cat(" --- \n")
      if (verbose>=1) cat("* Smoothing: ")
      # Smooth the points
      # - get smoothing spline
      ss <- splinefun( xi.vec.fix, L.fix )
      
      # Plot smoothed data
      xi.smooth <- seq(min(xi.vec.fix), max(xi.vec.fix), length=50 )
      L.smooth <- ss(xi.smooth )
      
      if ( do.plot) {
        keep.these <- is.finite(L.smooth) & !is.na(L.smooth)
        L.smooth <- L.smooth[ keep.these ] 
        xi.smooth <- xi.smooth[ keep.these ] 
        if ( verbose>=1 & any( !keep.these ) ) {
          cat(" (Some values of L are infinite or NA for the smooth; these are ignored)\n")
        }
        
        yrange <- range( L.smooth, na.rm=TRUE )
        
        plot( yrange ~ range(xi.vec),
              type="n",
              las=1,
              xlab=ifelse(xi.notation, expression(paste( xi," index")), expression(paste( italic(p)," index")) ),
              ylab=expression(italic(L)))
        lines( xi.smooth, L.smooth,
               lwd=2)
        rug( xi.vec )
        if (do.points) {
          points( L ~ xi.vec, pch=19)
        }
        if (add0) lines(xi.smooth[xi.smooth<1], L.smooth[xi.smooth<1], col="gray", lwd=2)
      }
      x <- xi.smooth
      y <- L.smooth
      
    } else {
      cat("  No valid values of the likelihood computed: smooth aborted\n",
          "  Consider trying another value for the input  method.\n")
    }
  }
  else {
    if ( do.plot) {
      
      keep.these <- is.finite(L) & !is.na(L)
      xi.vec <- xi.vec[ keep.these ] 
      L <- L[ keep.these ] 
      phi.vec <- phi.vec[ keep.these ]
      if ( verbose>=1 & any( keep.these ) ) {
        cat(" Some values of L are infinite or NA, and hence ignored\n")
      }
      
      yrange <- range( L, na.rm=TRUE )
      # Plot the data we have
      plot( yrange ~ range(xi.vec),
            type="n",     
            las=1,
            xlab=ifelse( xi.notation, expression(paste(xi," index")), expression(paste(italic(p)," index")) ),
            ylab=expression(italic(L)))
      lines( L ~ xi.vec, lwd=2)
      
      rug( xi.vec )
      if (do.points) {
        points( L ~ xi.vec, pch=19)
      }
      
    }
    x <- xi.vec
    y <- L
  }
  if (verbose>=2) cat(" Done\n")
  
  
  ### Maximum likelihood estimates
  if ( do.smooth ){
    # WARNING:  This spline does not necessarily pass
    #           through the given points.  This should
    #           be seen as a deficiency in the method.
    
    if (verbose>=2) cat(" Estimating phi:  ")
    
    # Find maximum from profile
    
    L.max <- max(y, na.rm=TRUE)
    xi.max <- x[ y==L.max ]
    
    # In some unusual cases, when add0=TRUE, the maximum of p/xi
    # may be between 0 and 1.
    # In this situation, we... set the mle to either 0 or 1,
    # whichever gives the larger values of the log-likelihood
    if ( (xi.max > 0)  & ( xi.max < 1 ) ) {
      L.max <- max( c( y[xi.vec==0], y[xi.vec==1]) )
      xi.max <- xi.vec[ L.max == y ]
      cat("MLE of",index.par,"is between 0 and 1, which is impossible.",
          "Instead, the MLE of",index.par,"has been set to",xi.max,
          ".  Please check your data and the call to tweedie.profile().")
    }
    
    # Now need to find mle of  phi  at this very value of  p
    # - Find bounds
    
    phi.1 <-   2 * max( phi.vec.fix, na.rm=TRUE )
    phi.2 <- 0.5 * min( phi.vec.fix, na.rm=TRUE )
    
    if ( phi.1 > phi.2 ) {
      phi.hi <- phi.1
      phi.lo <- phi.2 
    } else {
      phi.hi <- phi.2
      phi.lo <- phi.1
    }
    
    if (verbose>=2) cat(" Done\n")
    
    # Now, if the maximum happens to be at an endpoint,
    # we have to do things differently:
    if ( (xi.max==xi.vec[1]) | (xi.max==xi.vec[length(xi.vec)]) ) {
      
      if ( xi.max==xi.vec[1]) phi.max <- phi.vec[1]
      if ( xi.max==xi.vec[length(xi.vec)]) phi.max <- phi.vec[length(xi.vec)]
      #      phi.max <- phi.lo # and is same as phi.max
      warning("True maximum possibly not detected.")
      
    } else {
      
      # - solve
      
      if ( phi.method=="saddlepoint"){
        
        fit.model1 <- TDboost(formula_input,
                              data = mf,
                              distribution = list(name="EDM",alpha=xi.max),
                              shrinkage = 0.005,               
                              interaction.depth = 7,
                              bag.fraction = 0.5,
                              n.minobsinnode = 20,
                              keep.data = TRUE,
                              cv.folds = 5,
                              verbose = verbose1,...)
        best.iter1 <- TDboost.perf(fit.model1,method="cv", plot.it = FALSE)
        print(best.iter1)
        tmp1 <- predict.TDboost(fit.model1,mf, n.trees = best.iter1)
        mu <-  tmp1 # be careful prediction is for mu directly
        
        phi.max <- sum( tweedie.dev(y=ydata, mu=mu, power=xi.max) )/length( ydata )
        
      } else {
        #        phi.max <- nlm(p=phi.est, f=dtweedie.nlogl, 
        #                    hessian=FALSE,
        #                    power=p, mu=mu, y=data, 
        #                    gradient=dtweedie.dldphi)$estimate
        
        fit.model1 <- TDboost(formula_input,
                              data = mf,
                              distribution = list(name="EDM",alpha=xi.max),
                              shrinkage = 0.005,               
                              interaction.depth = 7,
                              bag.fraction = 0.5,
                              n.minobsinnode = 20,
                              keep.data = TRUE,
                              cv.folds = 5,
                              verbose = verbose1,...)
        best.iter1 <- TDboost.perf(fit.model1,method="cv", plot.it = FALSE)
        print(best.iter1)
        tmp1 <- predict.TDboost(fit.model1,mf, n.trees = best.iter1)
        mu <- tmp1 # be careful prediction is for mu directly
        phi.max <- optimize( f=dtweedie.nlogl, maximum=FALSE, interval=c(phi.lo, phi.hi ), 
                             # set lower limit phi.lo???
                             power=xi.max, mu=mu, y=ydata)$minimum
        
        #       Note that using the Hessian often produces computational problems
        #       for little (no?) advantage.
        
      }
    }
  } else {
    # Find maximum from computed data
    if (verbose>=2) cat(" Finding maximum likelihood estimates:  ")
    
    L.max <- max(L)
    
    xi.max   <- xi.vec  [ L == L.max ]
    phi.max <- phi.vec[ L == L.max ]
    
    if (verbose>=2) cat(" Done\n")
  }
  
  # Now report
  if ( verbose >= 2 ) {
    cat( "ML Estimates:  ",index.par,"=",xi.max," with phi=",phi.max," giving L=",L.max,"\n")
    cat(" ---\n")
  }
  
  # Now find approximate, nominal 95% CI info
  ht <- L.max - ( qchisq(conf.level, 1) / 2 )
  ci <- array( dim=2, NA )
  
  
  # if ( do.ci ) {
  #     if (verbose==2) cat("* Finding confidence interval:")
  #    if ( !do.smooth ) {
  #       warning("Confidence interval may be very inaccurate without smoothing.\n")
  # 		y <- L
  # 		x <- xi.vec
  #    }
  # 
  #    if ( do.plot ) {
  #       abline(h=ht, lty=2)
  #       title( sub=paste("(",100*conf.level,"% confidence interval)", sep="") )
  #    }
  # 
  #    # Now find limits on p:
  #    # - fit a smoothing spline the other way:  based on L predicting p
  #    #   so we can ensure that the limits on  p  are found OK
  # 
  #    # --- Left side ---
  #    cond.left <- (y < ht ) & (x < xi.max )
  #    if ( all(cond.left==FALSE) ) {
  #       warning("Confidence interval cannot be found: insufficient data to find left CI.\n")
  #    }else{
  # 		# Now find it!
  #       approx.left <- max( x[cond.left] )
  #       index.left <- seq(1, length(cond.left) )
  #       index.left <- index.left[x==approx.left]
  #       left.left <- max( 1, index.left - 5)
  #       left.right <- min( length(cond.left), index.left + 5)
  # 
  #       # Left-side spline
  #       ss.left <- splinefun( y[left.left:left.right], x[left.left: left.right] )
  #       ci.new.left <- ss.left( ht )
  # 
  #       ci[1] <- ci.new.left
  #    }
  # 
  #    # --- Right side ---
  #    cond.right <- (y < ht ) & (x > xi.max )
  #    if ( all( cond.right==FALSE ) ) {
  #       warning("Confidence interval cannot be found: insufficient data to find right CI.\n")
  #    }else{
  #      approx.right <- min( x[cond.right] )
  #       index.right <- seq(1, length(cond.right) )
  #       index.right <- index.right[x==approx.right]
  #       right.left <- max( 1, index.right - 5)
  #       right.right <- min( length(cond.left), index.right + 5)
  # 
  #       # Right-side spline
  #       ss.right <- splinefun( y[right.left:right.right], x[right.left: right.right] )
  #       ci.new.right <- ss.right(ht )
  # 
  #       ci[2] <- ci.new.right
  #    }
  #     if (verbose==2) cat(" Done\n")
  # }
  
  if ( fit.glm ) {
    out.glm <- glm.fit( x=model.x, y=ydata, weights=weights, offset=offset,
                        family=tweedie(var.power=xi.max, link.power=link.power) )
    
    if ( xi.notation){
      out <- list( y=y, x=x, ht=ht, L=L, xi=xi.vec, xi.max=xi.max, L.max=L.max, 
                   phi=phi.vec, phi.max=phi.max, ci=ci, method=method, phi.method=phi.method,
                   glm.obj = out.glm)
    } else {
      out <- list( y=y, x=x, ht=ht, L=L, p=p.vec, p.max=xi.max, L.max=L.max, 
                   phi=phi.vec, phi.max=phi.max, ci=ci, method=method, phi.method=phi.method,
                   glm.obj = out.glm)
    }
  } else {
    if (xi.notation ){
      out <- list( y=y, x=x, ht=ht, L=L, xi=xi.vec, xi.max=xi.max, L.max=L.max, 
                   phi=phi.vec, phi.max=phi.max, ci=ci, method=method, phi.method=phi.method)
      
    } else {
      out <- list( y=y, x=x, ht=ht, L=L, p=p.vec, p.max=xi.max, L.max=L.max, 
                   phi=phi.vec, phi.max=phi.max, ci=ci, method=method, phi.method=phi.method)
    }
    
  }
  
  if ( verbose ) cat("\n")
  
  invisible( out )
  
}



