library(tweedie)
library(mgcv)
library(TDboost)
library(statmod)
library(data.table)
library(cplm)
library(ggplot2)
library(latex2exp)
source("TDBoost_source.R")
source("EMTboost_HE.R")


rand_tweedie<- function(mu,...) {
  Y=rtweedie(1, mu = mu,...)
  Y
}

h1 <- function(x) {
  exp( -5 * (1-x[,1])^2 + (x[,2])^2 ) + exp( -5 * (x[,1])^2 + (1-x[,2])^2 )
}

rho <- 1.5 # True rho
phi <- 1 # True phi
P <- 2 # feature dimension
q <- 0.5 # c(1, 0.85, 0.75, 0.5, 0.25, 0.1)
print(paste("rho=",rho,"phi=",phi,"q=",q))

## Test data features
#  generate grid X data
T1 = seq(0,1,0.03)
T2 = seq(0,1,0.03)
TT = expand.grid(T1,T2)
colnames(TT) <- c("V1","V2")
truef<-h1(TT)

N<-500 # number of train data 
nobs <- N + nrow(TT) # number of train + test data


nrep <- 20 # 1

## Save the parameters
rho_TD_list <- phi_TD_list <- rho_EMT_list <- phi_EMT_list <- q_EMT_list <- array(NA,dim=c(nrep,2))


## For performance comparison 
MAD.list <- GINI.list <- array(NA,dim=c(nrep,2))
gini_matrix.list <- array(NA,dim=c(nrep,2,2))
gini_select.list <- rep(NA,times=nrep)
model.list <- c("TDboost", "EMTboost")

for(i in 1:nrep){
  
  print(paste("Replication",i))
  
  ## Train features
  #  train features are generated randomly
  XX <- matrix(runif(N * P, 0, 1), N, P)
  
  ## Zero-inflated Tweedie data
  #  Combine train and test data
  X = rbind(XX,TT)
  Fx <- h1(X)
  mu <- exp(Fx)
  #  Tweedie
  Y <- sapply(mu, rand_tweedie, xi = rho, phi=phi)
  #  Zero-inflated
  u <- runif(nobs,0,1)
  Y[which(u>q)] <- 0
  
  zero.index<-which(Y[1:N]==0); zero.size<-length(zero.index)
  nonzero.index<-which(Y[1:N]!=0); nonzero.size<-length(nonzero.index)
  dat<-data.frame(cbind(Y,X))
  colnames(dat)<-c("Y","V1","V2")
  
  ## Split train and test datasets
  train.dat<-dat[1:N,]
  test.dat<-dat[(N+1):nobs,]
  #  true loss
  loss<-Y[(N+1):nobs]
  
  ###############
  ## TDboost   ##
  ################################################
  ## estimating phi, rho via profile likelihood ##
  TD_out <- TDboost.profile(formula = Y~V1+V2,
                            data=train.dat,
                            p.vec=seq(1.1,1.9,0.03)[6:18], 
                            verbose=2,
                            phi.method="mle",
                            n.trees=15000, 
                            do.plot=FALSE,
                            method="series", 
                            verbose1=FALSE, 
                            do.smooth=FALSE)
  rho_TD.star<-TD_out$p.max
  print(paste('The estimated rho is:',rho_TD.star))
  rho_TD_list[i] <- rho_TD.star

  phi_TD.star<-TD_out$phi.max
  print(paste('The estimated phi is:',phi_TD.star))
  phi_TD_list[i] <- phi_TD.star

  #######################################
  ## TDboost: using the estimated rho  ##  
  TD_m <- TDboost(Y~V1+V2,                                         
                  data = train.dat,
                  distribution = list(name="EDM",alpha=rho_TD.star), 
                  shrinkage = 0.001,                               
                  n.trees = 15000,                                 
                  interaction.depth = 2,                           
                  bag.fraction = 1,                                
                  n.minobsinnode = 20,		                         
                  keep.data = TRUE,                                
                  cv.folds = 5,
                  verbose = FALSE)                                 
  best.iter <- TDboost.perf(TD_m, method="cv", plot.it = TRUE)
  pred_mu_TD <- predict.TDboost(TD_m, test.dat, best.iter)
  pred_f_TD <- predict.TDboost(TD_m, test.dat, best.iter, type='link')
 
  
  
  ################
  ## EMTboost   ##
  ###########################################################
  ## Senerio 3: estimating phi, rho via profile likelihood ##
  rho.list <- seq(1.3,1.7,0.1)
  EMT_out <- EMTboost.profile(train.dat=train.dat, 
                              test.dat=test.dat, 
                              rho.list=rho.list, 
                              maxsteps=20,
                              verbose=TRUE, 
                              plot.it=FALSE,
                              shrinkage = 0.001,
                              n.trees = 15000,             
                              interaction.depth = 2,         
                              bag.fraction = 1,           
                              n.minobsinnode = 20,	
                              cv.folds = 5)
  
  #  Optimal rho
  rho_EMT.star <- EMT_out$rho
  print(paste("EMTboost optimal rho:",rho_EMT.star))
  rho_EMT_list[i] <- rho_EMT.star
  
  ########################################
  ## EMTboost: using the estimated rho  ##  
  EMT_m <- EMTboost_HE(train.dat=train.dat, 
                       test.dat=test.dat, 
                       rho=rho_EMT.star, 
                       maxsteps=20,
                       verbose=TRUE, 
                       plot.it=FALSE,
                       shrinkage = 0.001,
                       n.trees = 15000,             
                       interaction.depth = 2,         
                       bag.fraction = 1,           
                       n.minobsinnode = 20,	
                       cv.folds = 5)
  
  #  Optimal q (prob. of Tweedie)
  q_EMT.star<- 1-EMT_m$q
  print(paste('q.hat:',q_EMT.star))
  q_EMT_list[i] <- q_EMT.star
  
  #  Optimal phi
  phi_EMT.star<-EMT_m$phi
  print(paste("phi.hat:",phi_EMT.star))
  phi_EMT_list[i] <- phi_EMT.star
  
  #  Predicted mu
  pred_f_EMT <- EMT_m$pred_test
  pred_mu_EMT<- q_EMT.star *exp(EMT_m$pred_test)

  ###############
  ## Performance
  ## TDboost
  #  MAD
  mad_TD <- mean(abs(pred_f_TD - truef))
  print(paste('Round',i,'TDboost MAD is',mad_TD))
  MAD.list[i,1] <- mad_TD
  #  GINI
  GINI.list[i,1] <- ((sum(loss*rank(pred_mu_TD,ties.method="last"))/sum(loss))-((N+1)/2))/((sum(loss*rank(loss,ties.method="last"))/sum(loss))-((N+1)/2))
  print(paste("TDboost GINI:",GINI.list[i,1]))
  
  ## EMTboost
  #  MAD
  mad_EMT<-mean(abs(pred_f_EMT - truef))
  print(paste('Round',i,'EMTboost MAD is',mad_EMT))
  MAD.list[i,2] <- mad_EMT
  #  GINI
  GINI.list[i,2] <- ((sum(loss*rank(pred_mu_EMT,ties.method="last"))/sum(loss))-((N+1)/2))/((sum(loss*rank(loss,ties.method="last"))/sum(loss))-((N+1)/2))
  print(paste("EMTboost GINI:",GINI.list[i,2]))
  
  ###############################
  #  Model comparison with GINI
  da <- data.frame(Loss = loss, TDboost = pred_mu_TD, EMTboost = pred_mu_EMT)
  gg <- gini(loss = "Loss", score  = c("TDboost","EMTboost"), base=NULL, data = da)
  gini_matrix.list[i,,]<-gg@gini
  print(gg@gini)
  gini_select.list[i]<-which.min(apply(gg@gini, 1, max))
  print(paste("The model selected through GINI matrix is:", model.list[gini_select.list[i]]))
}


gini.list
apply(MAD.list,2,mean)
apply(MAD.list,2,sd)
apply(GINI.list,2,mean)
apply(GINI.list,2,sd)
mean(rho_EMT_list);sd(rho_EMT_list)
mean(phi_EMT_list);sd(phi_EMT_list)
mean(q_EMT_list);sd(q_EMT_list)
mean(rho_TD_list);sd(rho_TD_list)
mean(phi_TD_list);sd(phi_TD_list)


