library(tweedie)
library(mgcv)
library(TDboost)
library(statmod)
library(data.table)
library(cplm)
library(ggplot2)
library(latex2exp)
source("generator.R")
source("TDBoost_source.R")
source("EMTboost_HE.R")

#################################
## Random Function Generator   ##
#################################


rand_tweedie<- function(mu,...) {
  Y=rtweedie(1, mu = mu,...)
  Y
}
phi <-  1
q <- 0.1 #c(0, 0.15, 0.25, 0.5, 0.75, 0.9)
P <- 10
rho <- 1.5
mid <- 1000
N <- 2000

nrep <- 20 # 1
## For saving the parameters
rho_TD_list <- phi_TD_list <- rho_EMT_list <- phi_EMT_list <- q_EMT_list <- array(NA,dim=c(nrep,2))


## For performance comparison 

MAD.list <- GINI.list <- array(NA,dim=c(nrep,2))
gini_matrix.list <- array(NA,dim=c(nrep,2,2))
gini_select.list <- rep(NA,times=nrep)
model.list <- c("TDboost", "EMTboost")

for (i in 1:nrep)
{
  
  print(paste("Replication",i))
  
  ###############################
  ## Random Function Generator ##
  X <- matrix(rnorm(N * P, 0, 1), N, P) 
  Fx <- Fx.generator(X,N,P)
  Y = sapply(exp(Fx), rand_tweedie, xi = rho, phi=phi)
  u <- runif(N,0,1)
  Y[which(u>q)] <- 0
  train_x <- X[1:mid, ]
  train_y <- Y[1:mid]
  test_x <- X[(mid+1):N, ]
  loss <- test_y <- Y[(mid+1):N]
  true_f <- Fx[(mid+1):N]
  train.dat = data.frame(cbind(train_y,train_x))
  test.dat = data.frame(cbind(test_y,test_x))
  colnames(train.dat) <- c("train_y","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")
  colnames(test.dat) <- c("test_y","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")

  
  ###############
  ## TDboost   ##
  ################################################
  ## estimating phi, rho via profile likelihood ##
  TD_out <- TDboost.profile(train_y~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10,
                         data=train.dat,
                         p.vec=seq(1.2,1.6,0.03), 
                         verbose=2,
                         phi.method="mle",
                         n.trees=3000, 
                         do.plot=TRUE,
                         method="series", 
                         verbose1=FALSE, 
                         do.smooth=FALSE)
  
  rho_TD.star <- TD_out$p.max
  print(paste('The estimated rho is:',rho_TD.star))
  rho_TD_list[i] <- rho_TD.star
  
  phi_TD.star <- TD_out$phi.max
  print(paste('The estimated phi is:',phi_TD.star))
  phi_TD_list[i] <- phi_TD.star
  
  #######################################
  ## TDboost: using the estimated rho  ##  
  TD_m <- TDboost(train_y~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10,    
                  data = train.dat,
                  distribution = list(name="EDM",alpha=rho_TD.star), 
                  shrinkage = 0.005,               
                  n.trees = 3000,              
                  interaction.depth = 6,         
                  bag.fraction = 0.5,            
                  n.minobsinnode = 20,		 
                  keep.data = FALSE,              
                  cv.folds = 5,
                  verbose = FALSE)
  
  best.iter <- TDboost.perf(TD_m, method="cv", plot.it = TRUE)
  pred_mu_TD <- predict.TDboost(TD_m, test.dat, best.iter)
  pred_f_TD <- predict.TDboost(TD_m, test.dat, best.iter, type='link')
  
  
  
  ################
  ## EMTboost   ##
  ###########################################################
  ## Senerio 3: estimating phi, rho via profile likelihood ##
  rho.list <- seq(1.44,1.56,0.03)
  EMT_out <- EMTboost.profile(train.dat=train.dat, 
                              test.dat=test.dat, 
                              rho.list=rho.list, 
                              maxsteps=20,
                              verbose=TRUE, 
                              plot.it=FALSE,
                              shrinkage = 0.005,
                              n.trees = 3000,             
                              interaction.depth = 6,         
                              bag.fraction = 0.5,           
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
                       shrinkage = 0.005,
                       n.trees = 3000,             
                       interaction.depth = 6,         
                       bag.fraction = 0.5,           
                       n.minobsinnode = 20,	
                       cv.folds = 5)
  
  
  #  Optimal q
  q_EMT.star<- 1 - EMT_m$q
  print(paste('q.hat:',q_EMT.star))
  q_EMT_list[i] <- q_EMT.star
  
  #  Optimal phi
  phi_EMT.star<-EMT_m$phi
  print(paste("phi.hat:",phi_EMT.star))
  phi_EMT_list[i] <- phi_EMT.star
  
  #  Predicted mu
  pred_f_EMT <- EMT_m$pred_test
  pred_mu_EMT <- q_EMT.star * exp(EMT_m$pred_test)
  
  ##################
  ## Performance  ##
  ## TDboost
  #  MAD
  mad_TD <- mean(abs(pred_f_TD - true_f))
  print(paste('Round',i,'TDboost MAD is',mad_TD))
  MAD.list[i,1] <- mad_TD
  #  GINI
  GINI.list[i,1] <- ((sum(loss*rank(pred_mu_TD,ties.method="last"))/sum(loss))-((N+1)/2))/((sum(loss*rank(loss,ties.method="last"))/sum(loss))-((N+1)/2))
  print(paste("TDboost GINI:",GINI.list[i,1]))
  
  ## EMTboost
  #  MAD
  mad_EMT <- mean(abs(pred_f_EMT - true_f)) 
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




