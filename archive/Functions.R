#################################################################################################################
### Title: Adjusted EWAS Using Compositional Data Transformations - Functions V1
### Current Analyst: Alexander Alsup
### Last Updated: (02/09/2026) 
### Notes: 
#################################################################################################################
# PACKAGES -----
packages <- c("ggplot2","here","MASS","MCMCpack","limma","betareg","tidyverse","writexl","readxl","easyCODA","limma","dplyr")
lapply(packages, require, character.only = TRUE)
#################################################################################################################
# PATHS AND OUTPUTS -----
Path_Parent = "" 

#################################################################################################################
# FUNCTIONS -----
##################################################
## Simulation: DATA GENERATION ----
# Data Generation Function 
generate.data <- function(data_sim=data_sim,N=N,depth=depth,B1=B1,cpg_num=cpg_num){
  # Cell Count Preparation
  n=N/2
  B0=data_sim$B0
  # Cell Counts, Proportions, & Main Effect (X)
  pois0 <- sapply(1:length(data_sim$Lambda0), function(x) rpois(n=n, lambda=data_sim$Lambda0[x]))
  pois1 <- sapply(1:length(data_sim$Lambda1), function(x) rpois(n=n, lambda=data_sim$Lambda1[x]))
  data.cells <- rbind(pois0,pois1)
  W <- data.frame(as.matrix(data.cells/rowSums(data.cells)));colnames(W)=names(data_sim$Lambda0)
  X=data.frame("X"=c(rep(0,times=N/2),rep(1,times=N/2)))
  # CpG Generation Temp Objects
  w_mat=as.matrix(W)
  x_vec = X$X
  # True methylation; Alphas and Betas
  Cell.Methylation = w_mat %*% B0
  Total.Methylation <- Cell.Methylation+(B1*x_vec)
  alpha = sapply(1:length(Total.Methylation), function(x) Total.Methylation[x]*depth)
  beta = sapply(1:length(Total.Methylation), function(x) (1-Total.Methylation[x])*depth)
  # CpG Data
  Y = matrix(NA,nrow=cpg_num,ncol=N)
  for(k in 1:cpg_num){
    Y_i=NULL
    for(i in 1:N){
      Y_i[i]<-rbeta(1,alpha[i],beta[i])
    }
    Y[k,] = Y_i
  }
  CpG_Variances = apply(Y,1,var,na.rm=TRUE)
  # Return
  result <- list("X"=X,"W"=W,"Y"=Y)
  return(result) 
}

## Simulation: SIMULATION FUNCTION -----
Sim.Function <- function(data_sim=data_sim,N=N,B1=B1,depth=depth,cpg_num=cpg_num){
  Data <- generate.data(data_sim=data_sim,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
  X_var = as.matrix(Data$X) ; W_var = as.matrix(Data$W) ; Y_var = as.matrix(Data$Y)
  # LR Selection 
  LR_Selection_Result = LR_Selection.v3(X=Data$X,W=W_var,var.total=(ncol(W_var)-1))
  lr_covariate = LR_Selection_Result[[1]]
  # CLARE
  limma_Covars = as.matrix(cbind(rep(1,times=length(lr_covariate)),X_var,lr_covariate))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,2]
  Detected_CpGs = Test_orig[which(Test_orig<=0.05)]
  CLARE = length(Detected_CpGs)/nrow(Y_var)
  # Proportion
  limma_Covars = as.matrix(cbind(X_var,W_var))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,1]
  Detected_CpGs = Test_orig[which(Test_orig<=0.05)]
  Proportion = length(Detected_CpGs)/nrow(Y_var)
  # Unadjusted
  limma_Covars = as.matrix(cbind(rep(1,times=length(lr_covariate)),X_var))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,2]
  Detected_CpGs = Test_orig[which(Test_orig<=0.05)]
  Unadjusted = length(Detected_CpGs)/nrow(Y_var)
  # Return
  results = data.frame("Effect"=B1,"Sample Size"=N,"CLARE"=CLARE,"Proportion"=Proportion,"Unadjusted"=Unadjusted)
  return(results)
}


##################################################
## Functional: Zero Proportion Scaling ----
# Compostional Scaling 
Scale_W <- function(W=W){
  W_lr=W
  if(any(W_lr==0)){
    W_lr[W_lr==0]=W_lr[W_lr==0]+1e-5
    W_lr=W_lr/rowSums(W_lr)
  } else {
    W_lr = as.matrix(W_lr)
  }
  return(W_lr)
}
## Functional: Matrix of All Unique ALRs ----
ALR_Matrix = function(W=W){
  W_lr=Scale_W(W=W)
  alr_covars <- ALR(W_lr,denom=1)$LR
  for(denom in 2:ncol(W_lr)){
    exclude=(denom-1)
    alr_covars_i <- ALR(W_lr,denom=denom)$LR[,-(1:exclude)]
    alr_covars <- cbind(alr_covars,alr_covars_i);
  }
  colnames(alr_covars)[ncol(alr_covars)] = paste0(names(W)[ncol(W)],"/",names(W)[ncol(W)-1])
  colnames(alr_covars)=gsub("\\/","\\.",colnames(alr_covars))
  return(alr_covars)
}
## Functional: LR Sparse Selection -----
LR_Selection.v3 <- function(Y=Y,X=X,W=W,var.total=8){
  W_lr = Scale_W(W=W)
  X_var = X$X
  #======= ALR Screen 
  LR_combos <- combn(ncol(W),2,simplify=FALSE)
  LR_specs <- lapply(LR_combos,function(pair) {
    list(numer=pair[1],denom=pair[2])
  })
  LR_matrix = sapply(LR_specs,function(spec) {
    log(W_lr[,spec$numer]/W_lr[,spec$denom])
  })
  LR_MSE=NULL
  for(i in 1:ncol(LR_matrix)){
    mod.fit1 <- glm(X_var~LR_matrix[,i])
    LR_MSE[i] <- mean(mod.fit1$residuals^2)
  }
  best_index=which(LR_MSE==min(LR_MSE,na.rm=TRUE))
  Best_LR_MSE = LR_MSE[best_index]
  Best_LR_Spec = list(LR_specs[[best_index]])
  Best_LR = as.matrix(LR_matrix[,best_index])
  base_spec = Best_LR_Spec[[1]]
  numer_base = base_spec$numer ; denom_base = base_spec$denom
  #======= SLR Screen Loop
  if(var.total>2){
    for(k in 1:(var.total-1)){
      remaining = setdiff(1:ncol(W_lr),c(numer_base,denom_base))
      LR_specs_numer = lapply(remaining, function(j){list(numer=c(numer_base,j),denom=c(denom_base)) })
      LR_specs_denom = lapply(remaining, function(j){list(numer=c(numer_base),denom=c(denom_base,j)) })
      LR_specs = c(LR_specs_numer,LR_specs_denom)
      LR_matrix <- sapply(LR_specs,function(spec){SLR(W_lr,numer=spec$numer,denom=spec$denom)$LR})
      LR_MSE=NULL
      for(i in 1:ncol(LR_matrix)){
        mod.fit1 <- glm(X_var~LR_matrix[,i])
        LR_MSE[i] <- mean(mod.fit1$residuals^2)
      }
      #
      challenger_index=min(which(LR_MSE==min(LR_MSE,na.rm=TRUE)))
      Best_LR_MSE = c(Best_LR_MSE,LR_MSE[challenger_index])
      Best_LR_Spec[[(k+1)]] = LR_specs[[challenger_index]]
      Best_LR = cbind(Best_LR,LR_matrix[,challenger_index]) 
      #
      base_spec = Best_LR_Spec[[(k+1)]]
      numer_base = base_spec$numer ; denom_base = base_spec$denom
    }
  }
  #======= Selecting Best LR
  # Specify Best LR
  best_index = min(which(Best_LR_MSE==min(Best_LR_MSE,na.rm=TRUE)))
  Final_LR_MSE = Best_LR_MSE[best_index]
  Final_LR_Spec = Best_LR_Spec[[best_index]]
  Final_LR = Best_LR[,best_index]
  # Last Check 
  base_spec = Final_LR_Spec
  numer_base = base_spec$numer ; denom_base = base_spec$denom
  SLR_Check = SLR(W_lr,numer=c(numer_base),denom=c(denom_base))$LR
  if(length(Final_LR == SLR_Check)!=length(Final_LR)){cat("Error: Output SLR Not Equal to Checked SLR")}
  # Packaging Results
  Numer_Cells = list(numer_base)
  Denom_Cells = list(denom_base)
  return(list("SLR"=Final_LR,"Numerator"=Numer_Cells,"Denominator"=Denom_Cells,"MSE"=Final_LR_MSE))
}
##################################################
## Method: CLARE  ----
CLARE_Test <- function(Y=Y,X=X,W=W,var.tot=var.tot){
  X_var = X$X
  #
  lr_covar <- LR_Selection.v3(Y=Y,X=X,W=W,var.total=var.tot)
  lr_covariate = lr_covar[[1]]
  #
  fit <- lm(Y~X_var+lr_covariate)
  result_detection <- ifelse(summary(fit)$coefficients["X_var",4] <= 0.05,1,0)
  return(result_detection)
}
## Method: Cell Proportion Adjusted EWAS ---- 
ProportionEWAS_Test<- function(Y=Y,X=X,W=W){
  df=data.frame("Y"=Y,"X"=X,W)
  model <- lm(Y~X+.,data=df)
  test_result <- ifelse(summary(model)$coefficients[2,4] <= 0.05,1,0) 
  return(test_result)
}

## Method: Unadjusted EWAS ----
# Model with No Covariates
Unadjusted_Test <- function(Y=Y,X=X,W=W){
  test_result<-NULL
  df=data.frame("Y"=Y,"X"=X,W)
  model <- lm(Y~X,data=df)
  test_result <- ifelse(summary(model)$coefficients[2,4] <= 0.05,1,0) 
  return(test_result)
}


#################################################################################################################
# PARAMETERS -----
##################################################
