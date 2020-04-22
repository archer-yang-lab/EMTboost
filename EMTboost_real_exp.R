library(cplm)
library(tweedie)
library(TDboost)
library(reshape2)
library(ggplot2)
source("TDBoost_source_real.R")
source("EMTboost_HE_real.R")



#############################################
## Undersampling the nonzeros in real data ##
## to generate a new dataset               ##
load("dat_clean.rda")
zero.index <- which(dat$CLM_AMT5==0)
zero.num <- length(zero.index)
nonzero.index <- which(dat$CLM_AMT5!=0)
nonzero.num <- length(nonzero.index)
zero.data <- dat[zero.index,]
nonzero.data <- dat[nonzero.index,]
over.para <- 1 
zero.over.num <- ceiling(zero.num*over.para)
under.para <- 0.15 # Undersampling fraction
nonzero.under.num <- ceiling(nonzero.num*under.para)
zero.percentage <- zero.over.num/(zero.over.num+nonzero.under.num)



nrep <- 20 # 1
## For saving the parameters
rho_TD_list <- phi_TD_list <- rho_EMT_list <- phi_EMT_list <- q_EMT_list <- array(NA,dim=c(nrep,2))


## For performance comparison 

MAE.list <- GINI.list <- array(NA,dim=c(nrep,2))
gini_matrix.list <- array(NA,dim=c(nrep,2,2))
gini_select.list <- rep(NA,times=nrep)
model.list <- c("TDboost", "EMTboost")
for (i in 1:nrep)
{
  
  print(paste("Replication",i))
  
  #####################################
  ## Generate dataset from real data ##
  under.nonzero.index <- sample(nonzero.index,nonzero.under.num,replace=TRUE)
  over.zero.index <- sample(zero.index,zero.over.num,replace=TRUE)
  U <- sample(c(0,1),size=nonzero.under.num,replace=TRUE,prob=c(0.5,0.5))
  train.nonzero.index <- under.nonzero.index[which(U==0)]
  test.nonzero.index <- under.nonzero.index[which(U==1)]
  V <- sample(c(0,1),size=zero.over.num,replace=TRUE,prob=c(0.5,0.5))
  train.zero.index <- over.zero.index[which(V==0)]
  test.zero.index <- over.zero.index[which(V==1)]
  train.data <- dat[c(train.nonzero.index,train.zero.index),]
  test.data <- dat[c(test.nonzero.index,test.zero.index),]
  
  loss <- test.data$CLM_AMT5
  N <- length(loss)
  
  
  ###############
  ## TDboost   ##
  ################################################
  ## estimating phi, rho via profile likelihood ##
  TD_out <- TDboost.profile.real(CLM_AMT5~KIDSDRIV+ 
                                 TRAVTIME+CAR_USE+BLUEBOOK+RETAINED+NPOLICY+CAR_TYPE+REVOLKED+MVR_PTS+AGE+HOMEKIDS+
                                 INCOME+GENDER+MARRIED+JOBCLASS+MAX_EDUC+AREA,    # formula
                                 data=train.data,
                                 p.vec=seq(1.2,1.6,0.05), 
                                 verbose=2,
                                 phi.method="mle",
                                 n.trees=250, 
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
  TD_m <- TDboost(CLM_AMT5~KIDSDRIV+ 
                  TRAVTIME+CAR_USE+BLUEBOOK+RETAINED+NPOLICY+CAR_TYPE+REVOLKED+MVR_PTS+AGE+HOMEKIDS+
                  INCOME+GENDER+MARRIED+JOBCLASS+MAX_EDUC+AREA,    
                  data = train.data,
                  distribution = list(name="EDM",alpha=rho_TD.star),
                  shrinkage = 0.005,               
                  n.trees = 250,              
                  interaction.depth = 7,        
                  bag.fraction = 0.5,            
                  n.minobsinnode = 20,		 
                  keep.data = FALSE,             
                  cv.folds = 5,
                  verbose = FALSE)
  
  best.iter <- TDboost.perf(TD_m, method="cv", plot.it = TRUE)
  pred_mu_TD <- predict.TDboost(TD_m, test.data, best.iter)
  pred_f_TD <- predict.TDboost(TD_m, test.data, best.iter, type='link')
  
  
  ################
  ## EMTboost   ##
  ###########################################################
  ## Senerio 3: estimating phi, rho via profile likelihood ##
  rho.list <- seq(1.3,1.7,0.05)
  EMT_out <- EMTboost.profile_real(train.dat=train.data, 
                                   test.dat=test.data, 
                                   rho.list=rho.list, 
                                   maxsteps=20,
                                   verbose=TRUE, 
                                   plot.it=FALSE,
                                   shrinkage = 0.005,
                                   n.trees = 250,             
                                   interaction.depth = 7,         
                                   bag.fraction = 0.5,           
                                   n.minobsinnode = 20,	
                                   cv.folds = 5)
  
  #  Optimal rho
  rho_EMT.star <- EMT_out$rho
  print(paste("EMTboost optimal rho:",rho_EMT.star))
  rho_EMT_list[i] <- rho_EMT.star
  
  ########################################
  ## EMTboost: using the estimated rho  ##  
  EMT_m <- EMTboost_HE_real(train.dat=train.data, 
                            test.dat=test.data, 
                            rho=rho_EMT.star, 
                            maxsteps=20,
                            verbose=TRUE, 
                            plot.it=FALSE,
                            shrinkage = 0.005,
                            n.trees = 250,             
                            interaction.depth = 7,         
                            bag.fraction = 0.5,           
                            n.minobsinnode = 20,	
                            cv.folds = 5)
  
  
  #  Optimal q
  q_EMT.star <- 1 - EMT_m$q
  print(paste('q.hat:',q_EMT.star))
  q_EMT_list[i] <- q_EMT.star
  
  #  Optimal phi
  phi_EMT.star <- EMT_m$phi
  print(paste("phi.hat:",phi_EMT.star))
  phi_EMT_list[i] <- phi_EMT.star
  
  #  Predicted mu
  pred_f_EMT <- EMT_m$pred_test
  pred_mu_EMT <- q_EMT.star * exp(pred_f_EMT)
  
  ##################
  ## Performance  ##
  ## TDboost
  #  MAE
  mae_TD <- mean(abs(pred_mu_TD - loss))
  print(paste('Round',i,'TDboost MAE is',mae_TD))
  MAE.list[i,1] <- mae_TD
  #  GINI
  GINI.list[i,1] <- ((sum(loss*rank(pred_mu_TD,ties.method="last"))/sum(loss))-((N+1)/2))/((sum(loss*rank(loss,ties.method="last"))/sum(loss))-((N+1)/2))
  print(paste("TDboost GINI:",GINI.list[i,1]))
  
  ## EMTboost
  #  MAD
  mae_EMT <- mean(abs(pred_mu_EMT - loss)) 
  print(paste('Round',i,'EMTboost MAE is',mae_EMT))
  MAE.list[i,2] <- mae_EMT
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
