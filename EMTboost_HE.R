# since we want a *maximum* likeihood estimator,
# define the negative.
dtweedie.nlogl.weight <- function(phi, y, mu, power, weight) {
  ans <- - sum( weight*log( dtweedie( y=y, mu=mu, phi=phi, power=power ) ) )
  #if ( is.infinite( ans ) ) {
  # If infinite, replace with saddlepoint estimate?
  #  ans <- sum( tweedie.dev(y=y, mu=mu, power=power) )/length( y )
  #}
  
  #attr(ans, "gradient") <- dtweedie.dldphi(y=y, mu=mu, phi=phi, power=power)
  #derivatives of the log-likelihood w.r.t. phi
  ans
}

EMTboost_HE <- function(train.dat, 
                        test.dat, 
                        rho=1.5, 
                        maxsteps=100, 
                        verbose=FALSE, plot.it=TRUE,
                        shrinkage = 0.001,
                        n.trees = 15000,             
                        interaction.depth = 2,         
                        bag.fraction = 1,           
                        n.minobsinnode = 20,	
                        cv.folds = 5)
{
  
  #initial should in (0,1)
  step <- 0
  Y <- train.dat[,1]
  n <- length(Y)
  zero.index <- which(Y==0) ; nonzero.index <- which(Y!=0)
  zero.size <- length(zero.index) ; nonzero.size <- length(nonzero.index)
  
  ####################
  ## initialization ##
  #  initialize mu
  mu <- sum(Y)/nonzero.size
  pred_mu <- rep(mu,times=n)
  print(paste("initial mu:",mu))
  
  #  initialize phi
  phi.saddle <- sum( tweedie.dev(y=Y, mu=pred_mu, power=rho) )/n
  phi.est <- phi.saddle
  low.limit <- min( 0.001, phi.saddle/2)
  high.limit <- 10*phi.est
  ans <- optimize(f=dtweedie.nlogl.weight, maximum=FALSE, interval=c(low.limit, high.limit),
                  power=rho, mu=pred_mu, y=Y, weight=as.numeric(Y!=0) )
  phi <- ans$minimum
  print(paste("initial phi:",phi))
  
  #  initialize q
  q <- zero.size/n
  print(paste("initial q:",q))
  
  #  initialize delta
  delta <- delta.new <- rep(1,times=n)
  numerator <- (1-q)*exp((1/phi)*(-(pred_mu[zero.index]^(2-rho))/(2-rho)))
  delta.zero <- numerator/(numerator+q)
  delta[zero.index] <- delta.zero
  
  ## log-likelihood
  L <- sum( log( (1-q)*dtweedie( y=Y, mu=pred_mu, phi=phi, power=rho )+q*as.numeric(Y==0) ) )
  if (verbose==TRUE){
    print(paste("initial log-likelihood:",L))
  }

  
  ## Outputs
  phi.list <- q.list <- mad.list <- L.list <- rep(NA,times=maxsteps)
  pred_test <- array(NA,dim=c(maxsteps,nrow(test.dat)))
  
  ###################
  ## EM Algorithm  ##
  while(step<maxsteps){

    step <- step+1
    if (verbose==TRUE){
      print(paste('step',step))
    }
    
    ############
    ## M-step ##
    #  update F
    TDboost.obj <- TDboost(formula = formula(train.dat),
                           weights = delta,
                           data = train.dat,
                           distribution = list(name="EDM",alpha=rho), 
                           shrinkage = shrinkage,              
                           n.trees = n.trees,             
                           interaction.depth = interaction.depth,         
                           bag.fraction = bag.fraction,           
                           n.minobsinnode = n.minobsinnode,		
                           keep.data = TRUE,             
                           cv.folds = cv.folds,
                           verbose = FALSE)  
    best.iter<-TDboost.perf(TDboost.obj,method="cv", plot.it = FALSE)
    pred_mu.new <- predict.TDboost(TDboost.obj, train.dat, best.iter)
    
    #  update phi
    phi.saddle <- sum( tweedie.dev(y=Y, mu=pred_mu.new, power=rho) )/length( Y )
    phi.est <- phi #use the former phi as the current estimation
    low.limit <- min( 0.001, phi.saddle/2)
    high.limit <- 10*phi.est
    ans <- optimize(f=dtweedie.nlogl.weight, maximum=FALSE, interval=c(low.limit, high.limit),
                    power=rho, mu=pred_mu.new, y=Y, weight=delta )
    phi.new<-ans$minimum
    if (verbose==TRUE){
      print(paste("current phi",phi.new))
    }
    
    #  update q (prob. of zeros)
    q.new<-sum(1-delta)/n
    print(paste("current q",q.new))
    
    ##############
    ## E-step   ##
    numerator<-(1-q.new)*exp((1/phi.new)*(-(pred_mu.new[zero.index]^(2-rho))/(2-rho)))
    delta.zero.new<-numerator/(numerator+q.new)
    delta.new[zero.index]<-delta.zero.new
    

    ## Calculate log-likelihood
    L.new<- sum( log( (1-q.new)*dtweedie( y=Y, mu=pred_mu.new, phi=phi.new, power=rho )+q.new*as.numeric(Y==0) ) )
    if (verbose==TRUE){
      print(paste("current log-likelihood:",L.new))
    }
    
    ## Update
    delta[zero.index]<-delta.zero.new
    pred_mu<-pred_mu.new
    phi<-phi.list[step]<-phi.new
    q<-q.list[step]<-q.new
    L<-L.list[step]<-L.new
    #predict test
    pred_test[step,] <- predict.TDboost(TDboost.obj, test.dat, best.iter, type='link')
    
    

  }
  if (plot.it==TRUE){
    par(mfrow=c(2,2))
    plot(L.list~seq(1,maxsteps,1),
         type='b',
         las=1,
         xlab=expression("iteration step"),
         ylab=expression(italic(L)))
    plot(q.list~seq(1,maxsteps,1),
         type='b',
         las=1,
         xlab=expression("iteration step"),
         ylab=expression(italic(q)))
    plot(phi.list~seq(1,maxsteps,1),
         type='b',
         las=1,
         xlab=expression("iteration step"),
         ylab=expression(italic("phi")))
  }
  
  optimal.item<-which(L.list[(.5*maxsteps):maxsteps]==max(L.list[(.5*maxsteps):maxsteps]))+(.5*maxsteps)-1
  
  return(list(L.list=L.list, q.list=q.list, phi.list=phi.list, optimal.step=optimal.item, 
              pred_test=pred_test[optimal.item,], L=L.list[optimal.item], 
              phi=phi.list[optimal.item], q=q.list[optimal.item]))

}


EMTboost.profile <- function(train.dat, 
                             test.dat,
                             rho.list,
                             maxsteps=100,
                             verbose=TRUE, plot.it=FALSE,
                             shrinkage = 0.001,
                             n.trees = 15000,             
                             interaction.depth = 2,         
                             bag.fraction = 1,           
                             n.minobsinnode = 20,	
                             cv.folds = 5)
{
  EMTboost.object <- list()
  num.rho <- length(rho.list)
  L.optimal <- rep(NA, times=num.rho)
  for (m in 1:num.rho){
    rho <- rho.list[m]
    print(paste("rho:",rho))
    EMTboost <- EMTboost_HE(train.dat=train.dat, 
                            test.dat=test.dat, 
                            rho=rho, 
                            maxsteps=maxsteps,
                            verbose=verbose, 
                            plot.it=plot.it, 
                            shrinkage = shrinkage,
                            n.trees = n.trees,             
                            interaction.depth = interaction.depth,         
                            bag.fraction = bag.fraction,           
                            n.minobsinnode = n.minobsinnode,	
                            cv.folds = cv.folds)
    EMTboost.object[[m]] <- EMTboost
    L.optimal[m] <- EMTboost$L
  }
  
  ## Choose optimal rho through log-likelihood
  opt.index <- which(L.optimal==max(L.optimal))
  #  Optimal rho
  rho.star <- rho.list[opt.index]
  
  return(list(rho = rho.star))
  
}
