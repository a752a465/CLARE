#################################################################################################################
### Title: Adjusted EWAS Using Compositional Data Transformations - Supplementary Simulation #1 ----
### Current Analyst: Alexander Alsup (a752a465@kumc.edu) ----
### Last Updated: (08/14/2026) ----
### Notes: ----
#################################################################################################################
# Source Functions & Parent Filepath ----
# Replace with your own directory where applicable
Parent_Directory = ""
source(paste0(Parent_Directory,"Supplementary Simulation Prep.R"))
#################################################################################################################
# 1. SUPP. SIMULATION 1: CONTINUOUS EXPOSURE -----
## Functions  ----
# Generalized CLARE LR Search
LR_Search_General <- function(W, fit_fun, var.total, minimize = TRUE, early_stop = FALSE, verbose = FALSE, ...) {
  
  W_lr <- Scale_W(W)
  p <- ncol(W_lr)
  
  # Conditional Check for the var.total input
  if (var.total > p) stop("var.total cannot exceed the number of cell types in W")
  
  # Functions to control the selection process
  pick_best <- function(scores) if (minimize) which.min(scores) else which.max(scores)
  is_better <- function(new, old) if (minimize) new < old else new > old
  
  # Evaluate fit_fun on every candidate spec at a given step; returns scores + a record data.frame
  eval_candidates <- function(specs, lr_matrix, step) {
    fits <- lapply(seq_len(ncol(lr_matrix)), function(i) fit_fun(lr_matrix[, i], ...))
    scores <- sapply(fits, function(f) f$score)
    
    has_extra <- !is.null(fits[[1]]$extra)
    extra_df <- if (has_extra) {
      do.call(rbind, lapply(fits, function(f) as.data.frame(as.list(f$extra))))
    } else NULL
    
    numer_lab <- sapply(specs, function(s) paste(colnames(W)[s$numer], collapse = "+"))
    denom_lab <- sapply(specs, function(s) paste(colnames(W)[s$denom], collapse = "+"))
    
    record <- data.frame(size = step, numerator = numer_lab, denominator = denom_lab,
                         score = scores, stringsAsFactors = FALSE)
    if (has_extra) record <- cbind(record, extra_df)
    
    list(scores = scores, record = record)
  }
  
  All_Evaluated <- list()
  Best_Path <- list()
  #======= Step 1: screen all pairwise (single-cell) log-ratios
  LR_combos <- combn(p, 2, simplify = FALSE)
  LR_specs  <- lapply(LR_combos, function(pair) list(numer = pair[1], denom = pair[2]))
  LR_matrix <- sapply(LR_specs, function(spec) log(W_lr[, spec$numer] / W_lr[, spec$denom]))
  
  ev <- eval_candidates(LR_specs, LR_matrix, step = 1)
  All_Evaluated[[1]] <- ev$record
  
  best_index <- pick_best(ev$scores)
  Best_Score <- ev$scores[best_index]
  Best_Spec  <- list(LR_specs[[best_index]])
  Best_LR    <- as.matrix(LR_matrix[, best_index])
  Best_Path[[1]] <- ev$record[best_index, ]
  
  base_spec  <- Best_Spec[[1]]
  numer_base <- base_spec$numer
  denom_base <- base_spec$denom
  
  if (verbose) {
    cat(sprintf("Step 1 | numer=%s denom=%s | score=%.5f\n",
                paste(colnames(W)[numer_base], collapse = "+"),
                paste(colnames(W)[denom_base], collapse = "+"),
                Best_Score[1]))
  }
  
  #======= Stepwise growth: add one more cell type (to either numerator or denominator) per step
  if (var.total > 2) {
    for (k in 1:(var.total - 1)) {
      remaining <- setdiff(1:p, c(numer_base, denom_base))
      if (length(remaining) == 0) break
      
      LR_specs_numer <- lapply(remaining, function(j) list(numer = c(numer_base, j), denom = denom_base))
      LR_specs_denom <- lapply(remaining, function(j) list(numer = numer_base, denom = c(denom_base, j)))
      LR_specs_k <- c(LR_specs_numer, LR_specs_denom)
      
      LR_matrix_k <- sapply(LR_specs_k, function(spec) SLR(W_lr, numer = spec$numer, denom = spec$denom)$LR)
      
      ev_k <- eval_candidates(LR_specs_k, LR_matrix_k, step = k + 1)
      All_Evaluated[[k + 1]] <- ev_k$record
      
      challenger_index <- pick_best(ev_k$scores)
      challenger_score  <- ev_k$scores[challenger_index]
      
      if (early_stop && !is_better(challenger_score, Best_Score[k])) {
        if (verbose) cat(sprintf("Step %d | no improvement (score=%.5f vs %.5f) -> stopping\n",
                                 k + 1, challenger_score, Best_Score[k]))
        break
      }
      
      Best_Score <- c(Best_Score, challenger_score)
      Best_Spec[[k + 1]] <- LR_specs_k[[challenger_index]]
      Best_LR <- cbind(Best_LR, LR_matrix_k[, challenger_index])
      Best_Path[[k + 1]] <- ev_k$record[challenger_index, ]
      
      base_spec  <- Best_Spec[[k + 1]]
      numer_base <- base_spec$numer
      denom_base <- base_spec$denom
      
      if (verbose) {
        cat(sprintf("Step %d | numer=%s denom=%s | score=%.5f\n",
                    k + 1,
                    paste(colnames(W)[numer_base], collapse = "+"),
                    paste(colnames(W)[denom_base], collapse = "+"),
                    Best_Score[k + 1]))
      }
    }
  }
  
  #======= Select best step along the path (matches original CLARE logic when early_stop = FALSE)
  best_index_final <- pick_best(Best_Score)
  Final_Score <- Best_Score[best_index_final]
  Final_Spec  <- Best_Spec[[best_index_final]]
  Final_LR    <- Best_LR[, best_index_final]
  
  Best_Path_df <- do.call(rbind, Best_Path)
  rownames(Best_Path_df) <- NULL
  Best_Path_df$selected <- seq_len(nrow(Best_Path_df)) == best_index_final
  
  All_Evaluated_df <- do.call(rbind, All_Evaluated)
  rownames(All_Evaluated_df) <- NULL
  
  list(
    LR          = Final_LR,
    Numerator   = colnames(W)[Final_Spec$numer],
    Denominator = colnames(W)[Final_Spec$denom],
    Numerator_idx   = Final_Spec$numer,
    Denominator_idx = Final_Spec$denom,
    Score       = Final_Score,
    Best_Path   = Best_Path_df # one row per size: the winner at that size
    # ,All_Evaluated = All_Evaluated_df # every candidate log-ratio evaluated, every size
  )
}

# Fit function for Continuous Exposure
fit_function_suppsim1 <- function(lr, X, covariates = NULL, ...) {
  # Collecting the log-ratio and other data into a single data frame
  df <- data.frame(lr = lr, X = X)
  rhs <- "lr"
  # Conditionally adding the covariates if present
  if (!is.null(covariates)) {
    df <- cbind(df, covariates)
    rhs <- paste(c("lr", colnames(covariates)), collapse = " + ")
  }
  
  # Fitting the model: does the log-ratio (plus covariates) explain the exposure?
  form <- as.formula(paste("X ~", rhs))
  mod <- glm(form, data = df)
  s <- summary(mod)
  
  ## Results
  list(
    # Our "score" output is what the greedy selection uses -- MSE, minimized
    score = mean(mod$residuals^2),
    
    # The "extra" output is additional information we can look at
    extra = c(
      lr_coef = unname(coef(mod)["lr"]),
      lr_p    = unname(s$coefficients["lr", "Pr(>|t|)"]),
      AIC     = AIC(mod),
      R2      = 1 - (sum(mod$residuals^2) / sum((X - mean(X))^2))
    )
  )
}

# Modified Sim Function
Sim.Function.SuppSim1 <- function(Data_sim=Data){
  X_var = as.matrix(Data_sim$X) ; W_var = as.matrix(Data_sim$W) ; Y_var = as.matrix(Data_sim$Y)
  # LR Selection -- X is now passed explicitly so it reaches fit_function_suppsim1 via ...
  # instead of silently falling back to a same-named global variable.
  LR_Selection_Result = LR_Search_General(W=W_var, X=X_var, fit_fun=fit_function_suppsim1, var.total=length(lambda0))
  lr_covariate = LR_Selection_Result$LR
  # CLARE
  limma_Covars = as.matrix(cbind(rep(1,times=length(lr_covariate)),X_var,lr_covariate))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,2]
  Est_orig  <- fit$coefficients[,2]
  Detected_CpGs = Test_orig[which(Test_orig<=0.05)]
  CLARE = length(Detected_CpGs)/nrow(Y_var)
  CLARE_Bias = mean(Est_orig - B1, na.rm=TRUE)
  CLARE_RMSE = sqrt(mean((Est_orig - B1)^2, na.rm=TRUE))
  CLARE_Variance = mean((Est_orig - mean(Est_orig, na.rm=TRUE))^2, na.rm=TRUE)
  # Proportion
  limma_Covars = as.matrix(cbind(X_var,W_var))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,1]
  Est_orig  <- fit$coefficients[,1]
  Detected_CpGs = Test_orig[which(Test_orig<=0.05)]
  Proportion = length(Detected_CpGs)/nrow(Y_var)
  Proportion_Bias = mean(Est_orig - B1, na.rm=TRUE)
  Proportion_RMSE = sqrt(mean((Est_orig - B1)^2, na.rm=TRUE))
  Proportion_Variance = mean((Est_orig - mean(Est_orig, na.rm=TRUE))^2, na.rm=TRUE)
  # Unadjusted
  limma_Covars = as.matrix(cbind(rep(1,times=length(lr_covariate)),X_var))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,2]
  Est_orig  <- fit$coefficients[,2]
  Detected_CpGs = Test_orig[which(Test_orig<=0.05)]
  Unadjusted = length(Detected_CpGs)/nrow(Y_var)
  Unadjusted_Bias = mean(Est_orig - B1, na.rm=TRUE)
  Unadjusted_RMSE = sqrt(mean((Est_orig - B1)^2, na.rm=TRUE))
  Unadjusted_Variance = mean((Est_orig - mean(Est_orig, na.rm=TRUE))^2, na.rm=TRUE)
  # Return
  results = data.frame("Effect"=B1,"Sample Size"=N,
                       "CLARE"=CLARE,"Proportion"=Proportion,"Unadjusted"=Unadjusted,
                       "CLARE_Bias"=CLARE_Bias,"CLARE_RMSE"=CLARE_RMSE,"CLARE_Variance"=CLARE_Variance,
                       "Proportion_Bias"=Proportion_Bias,"Proportion_RMSE"=Proportion_RMSE,"Proportion_Variance"=Proportion_Variance,
                       "Unadjusted_Bias"=Unadjusted_Bias,"Unadjusted_RMSE"=Unadjusted_RMSE,"Unadjusted_Variance"=Unadjusted_Variance)
  return(results)
}

# Data Generation for the simulation
generate.data.suppsim1 <- function(N=N,depth=depth,B1=B1,cpg_num=cpg_num){
  # Data Generation
  ## Lambda Matrix
  # X = floor(runif(N,20,71))
  X = floor(rnorm(N,60,20))
  X_tau <- outer(X, tau)
  lambda=matrix(rep(lambda0,times=N),nrow=N,ncol=length(lambda0),byrow=TRUE)
  lambda=lambda+X_tau
  # Cell Counts, Proportions, & Main Effect (X)
  data.cells <- matrix(
    rpois(n = N * ncol(lambda), lambda = as.vector(lambda)),
    nrow = N, ncol = ncol(lambda)
  ) ; colnames(data.cells) = colnames(lambda)
  W <- data.frame(as.matrix(data.cells/rowSums(data.cells)));colnames(W)=colnames(lambda)
  # True methylation; Alphas and Betas
  w_mat=as.matrix(W)
  Cell.Methylation = w_mat %*% B0
  Total.Methylation <- Cell.Methylation+(B1*X)
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

## Sim Parameters -----

## Global Parameters
Nsims= 50 # Outer Replicates - Number of Unique Cell Mixture Samples
cpg_num=500 # Inner Replicates - Number of CpGs generated per Cell Mixture Sample
total.cells=1000 # The total cells used for measuring DNAm
depth=200 # The probe depth

## Simulation Parameters
sim_params_suppsim1 <- list(
  # DMP Power / True X Effect
  sim_params_DMP_effect = list(
    lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6),
    tau = c("Neut"=-0.1,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+0.05,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0),
    B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647),
    B1 = -0.0004 # Main Effect
  ),
  # Type I Error / Null X Effect
  sim_params_DMP_noeffect = list(
    lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6),
    tau = c("Neut"=-0.1,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+0.05,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0),
    B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647),
    B1 = 0 # Main Effect
  ),
  # ncDMP Power / True X Effect
  sim_params_ncDMP_effect = list(
    lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6),
    tau = c("Neut"=-0.1,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+0.05,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0),
    B0=c("Neut"=0.5,"Mono"=0.5,"Eos"=0.5,"Baso"=0.5,"CD4_Naive"=0.5,"CD4_Mem"=0.5,"CD8_Naive"=0.5,"CD8_Mem"=0.5,"B_Naive"=0.5,"B_Mem"=0.5,"Treg"=0.5,"NK"=0.5),
    B1 = -0.0004 # Main Effect
  ),
  # ncDMP Type I Error / Null X Effect
  sim_params_ncDMP_noeffect = list(
    lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6),
    tau = c("Neut"=-0.1,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+0.05,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0),
    B0=c("Neut"=0.5,"Mono"=0.5,"Eos"=0.5,"Baso"=0.5,"CD4_Naive"=0.5,"CD4_Mem"=0.5,"CD8_Naive"=0.5,"CD8_Mem"=0.5,"B_Naive"=0.5,"B_Mem"=0.5,"Treg"=0.5,"NK"=0.5),
    B1 = 0 # Main Effect
  )
  #
)

#####################
## Simulation ----
### Set.Seed() ----
set.seed(999848)
### Simulation : cDMP : Type I Error / Null X Effect -----

# Scenario Parameters
sim_params = sim_params_suppsim1$sim_params_DMP_noeffect
lambda0= sim_params$lambda0 ;tau=sim_params$tau ; B0=sim_params$B0 ; B1 = sim_params$B1
rho=total.cells/sum(lambda0) ; lambda0=lambda0*rho ; tau=tau*rho
# Simulation Loop
Sim_Results = data.frame()
# Outer Loop for Sample Size
for(N_i in 1:length(sample_size)){
  results_i <- data.frame()
  N=sample_size[N_i]
  # Inner Loop for sims
  for(sim in 1:Nsims){
    Data <- generate.data.suppsim1(N=N,depth=depth,B1=B1,cpg_num=cpg_num)
    results_init <- Sim.Function.SuppSim1(Data_sim=Data)
    results_i = bind_rows(results_i,results_init)
  }
  # Outer Loop Message Update
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}
Sim_Results_cDMP = Sim_Results 

### Simulation : cDMP : Power / True X Effect -----

# Scenario Parameters
sim_params = sim_params_suppsim1$sim_params_DMP_effect
lambda0= sim_params$lambda0 ;tau=sim_params$tau ; B0=sim_params$B0 ; B1 = sim_params$B1
rho=total.cells/sum(lambda0) ; lambda0=lambda0*rho ; tau=tau*rho
# Simulation Loop
Sim_Results = data.frame()
# Outer Loop for Sample Size
for(N_i in 1:length(sample_size)){
  results_i <- data.frame()
  N=sample_size[N_i]
  # Inner Loop for sims
  for(sim in 1:Nsims){
    Data <- generate.data.suppsim1(N=N,depth=depth,B1=B1,cpg_num=cpg_num)
    results_init <- Sim.Function.SuppSim1(Data_sim=Data)
    results_i = bind_rows(results_i,results_init)
  }
  # Outer Loop Message Update
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}
Sim_Results_cDMP = bind_rows(Sim_Results_cDMP,Sim_Results) 




### Simulation : ncDMP : Type I Error / Null X Effect -----

# Scenario Parameters
sim_params = sim_params_suppsim1$sim_params_ncDMP_noeffect
lambda0= sim_params$lambda0 ;tau=sim_params$tau ; B0=sim_params$B0 ; B1 = sim_params$B1
rho=total.cells/sum(lambda0) ; lambda0=lambda0*rho ; tau=tau*rho
# Simulation Loop
Sim_Results = data.frame()
# Outer Loop for Sample Size
for(N_i in 1:length(sample_size)){
  results_i <- data.frame()
  N=sample_size[N_i]
  # Inner Loop for sims
  for(sim in 1:Nsims){
    Data <- generate.data.suppsim1(N=N,depth=depth,B1=B1,cpg_num=cpg_num)
    results_init <- Sim.Function.SuppSim1(Data_sim=Data)
    results_i = bind_rows(results_i,results_init)
  }
  # Outer Loop Message Update
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}
Sim_Results_ncDMP = Sim_Results 

### Simulation : ncDMP : Power / True X Effect -----

# Scenario Parameters
sim_params = sim_params_suppsim1$sim_params_ncDMP_effect
lambda0= sim_params$lambda0 ;tau=sim_params$tau ; B0=sim_params$B0 ; B1 = sim_params$B1
rho=total.cells/sum(lambda0) ; lambda0=lambda0*rho ; tau=tau*rho
# Simulation Loop
Sim_Results = data.frame()
# Outer Loop for Sample Size
for(N_i in 1:length(sample_size)){
  results_i <- data.frame()
  N=sample_size[N_i]
  # Inner Loop for sims
  for(sim in 1:Nsims){
    Data <- generate.data.suppsim1(N=N,depth=depth,B1=B1,cpg_num=cpg_num)
    results_init <- Sim.Function.SuppSim1(Data_sim=Data)
    results_i = bind_rows(results_i,results_init)
  }
  # Outer Loop Message Update
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}
Sim_Results_ncDMP = bind_rows(Sim_Results_ncDMP,Sim_Results) 

#####################
### Prepping Results ----

# Notes for Individual Scenarios
Sim_Results_cDMP$Scenario = "cDMP"
Sim_Results_ncDMP$Scenario = "ncDMP"

# Results Simple-ish
res_wide = bind_rows(Sim_Results_cDMP,Sim_Results_ncDMP)%>%
  dplyr::filter(!is.na(Effect))

# Complete Data Frame
Results_Final_Sim1 = bind_rows(Sim_Results_cDMP
                               ,Sim_Results_ncDMP
)%>%
  mutate(Effect=case_when(Effect==0~"Type I Error",TRUE~"Power"))%>%
  rename(CLARE_Detection=CLARE,Proportion_Detection=Proportion,Unadjusted_Detection=Unadjusted)%>%
  dplyr::filter(!is.na(Sample.Size))%>%
  pivot_longer(cols = matches("_Bias$|_RMSE$|_Detection$|_Variance$"), names_to = c("Method","Metric"), names_sep = "_")

### Double Checking RMSE Capture ----
max(abs(as.numeric(res_wide$CLARE_RMSE)^2 - (as.numeric(res_wide$CLARE_Bias)^2 + as.numeric(res_wide$CLARE_Variance))))
max(abs(as.numeric(res_wide$Unadjusted_RMSE)^2 - (as.numeric(res_wide$Unadjusted_Bias)^2 + as.numeric(res_wide$Unadjusted_Variance))))
max(abs(as.numeric(res_wide$Proportion_RMSE)^2 - (as.numeric(res_wide$Proportion_Bias)^2 + as.numeric(res_wide$Proportion_Variance))))


### Summaries -----
### Bias/RMSE/Variance

metric_means <- Results_Final_Sim1 %>%
  mutate(value=as.numeric(value))%>%
  group_by(Effect,Scenario,Method,Metric)%>%
  summarise(Metric_Mean=mean(value,na.rm=TRUE))%>%
  arrange(Metric,Scenario,Effect,Method)

# metric_means2 <- Results_Final_Sim1 %>%
#   mutate(value=as.numeric(value))%>%
#   group_by(Effect,Method,Metric)%>%
#   summarise(Metric_Mean=mean(value,na.rm=TRUE))%>%
#   arrange(Metric,Effect,Method)

#####################
## Quick Figures ----
## Overall prep ----
#
colors <- c("#15bf3a",
            "#4075d6",
            "#631307")
shapes <- c(15,1,4)
sample.size.breaks=unique(Results_Final_Sim1$Sample.Size)
Simulation.title = "Supplementary Simulation 1:"
#
df <-  Results_Final_Sim1 %>% 
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Unadjusted")))%>%
  mutate(Effect=factor(Effect,levels=c("Type I Error","Power")))

## Mean Power & Type I Error Rate with SD Ribbons ----
#
df_graph = df %>%
  dplyr::filter(Method != "Unadjusted" & !is.na(Scenario) & Metric == "Detection")%>%
  group_by(Effect,Sample.Size,Scenario,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE),
            sd=sd(value,na.rm=TRUE))%>%
  ungroup()
#
plot1 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method))+
  geom_line(size=1,alpha=0.5)+
  geom_point(size=2,alpha=0.9)+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,fill=Method),alpha=0.1,color=NA)+
  #stat_summary(fun=mean,geom="line",size=1,alpha=0.5)+
  #stat_summary(fun.data=mean_cl_normal,geom="ribbon",alpha=0.2)
  labs(title=Simulation.title,
       subtitle="Power & Type I Error Rate",
       x="Sample Size",
       y="Detection Rate")+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="fixed")+
  scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1),minor_breaks=NULL)+
  geom_hline(yintercept=0.05,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plot1

## Mean Bias with SD Ribbons ----
#
df_graph = df %>%
  mutate(Effect=case_when(Effect=="Type I Error"~"No Effect",
                          Effect=="Power"~"True Effect",TRUE~NA))%>%
  mutate(Effect=factor(Effect,levels=c("No Effect","True Effect")))%>%
  dplyr::filter(!is.na(Scenario) & Metric == "Bias"
                #& Method != "Unadjusted"
  )%>%
  group_by(Effect,Sample.Size,Scenario,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE),
            sd=sd(value,na.rm=TRUE))%>%
  ungroup()
#
plot2 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method))+
  geom_line(size=1,alpha=0.5)+
  geom_point(size=2,alpha=0.9)+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,fill=Method),alpha=0.1,color=NA)+
  labs(title=Simulation.title,
       subtitle="Mean Bias in Estimated Effect",
       x="Sample Size",
       y=expression(paste("Mean Bias (", hat(beta)[1], " - ", beta[1], ")")))+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="free_y")+
  geom_hline(yintercept=0,color="red",alpha=0.3,linetype="dashed")+   # unbiasedness reference
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plot2

## Mean RMSE with SD Ribbons ----
df_graph = df %>%
  mutate(Effect=case_when(Effect=="Type I Error"~"No Effect",
                          Effect=="Power"~"True Effect",TRUE~NA))%>%
  mutate(Effect=factor(Effect,levels=c("No Effect","True Effect")))%>%
  dplyr::filter(!is.na(Scenario) & Metric == "RMSE"
                #& Method != "Unadjusted"
  )%>%
  group_by(Effect,Sample.Size,Scenario,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE),
            sd=sd(value,na.rm=TRUE))%>%
  ungroup()
#
plot3 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method))+
  geom_line(size=1,alpha=0.5)+
  geom_point(size=2,alpha=0.9)+
  geom_ribbon(aes(ymin=pmax(mean-sd,0),ymax=mean+sd,fill=Method),alpha=0.1,color=NA)+
  labs(title=Simulation.title,
       subtitle="Mean RMSE of Estimated Effect",
       x="Sample Size",
       y="Mean RMSE")+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="free_y")+
  # no zero-reference line here -- RMSE has no "unbiased" analog, lower is simply better
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plot3

## Mean Variance with SD Ribbons ----
df_graph = df %>%
  mutate(Effect=case_when(Effect=="Type I Error"~"No Effect",
                          Effect=="Power"~"True Effect",TRUE~NA))%>%
  mutate(Effect=factor(Effect,levels=c("No Effect","True Effect")))%>%
  dplyr::filter(!is.na(Scenario) & Metric == "Variance"
                #& Method != "Unadjusted"
  )%>%
  group_by(Effect,Sample.Size,Scenario,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE),
            sd=sd(value,na.rm=TRUE))%>%
  ungroup()
#
plot4 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method))+
  geom_line(size=1,alpha=0.5)+
  geom_point(size=2,alpha=0.9)+
  geom_ribbon(aes(ymin=pmax(mean-sd,0),ymax=mean+sd,fill=Method),alpha=0.1,color=NA)+
  labs(title=Simulation.title,
       subtitle="Mean Variance of Estimated Effect",
       x="Sample Size",
       y="Mean Variance")+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="free_y")+
  # no zero-reference line here -- RMSE has no "unbiased" analog, lower is simply better
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plot4


