#################################################################################################################
### Title: Adjusted EWAS Using Compositional Data Transformations - Supplementary Simulation #2 ----
### Current Analyst: Alexander Alsup (a752a465@kumc.edu) ----
### Last Updated: (08/14/2026) ----
### Notes: ----
#################################################################################################################
# Source Functions & Parent Filepath ----
# Replace with your own directory where applicable
Parent_Directory = ""
source(paste0(Parent_Directory,"Supplementary Simulation Prep.R"))
#################################################################################################################

## SUPPLEMENTARY SIMULATION 2.1 - Deconvolution Error - Neutrophils/CD4 Naive -----

## Functions -----
# Function for adding Deconvolution Error
add_deconvolution_error <- function(W, floor = 1e-4, seed = NULL, scale = 100, type = "standard", err_std = 0.5) {
  if (type == "standard") {
    rmse_vector <- c(
      Neut = err_std, Mono = err_std, Eos = err_std, Baso = err_std,
      CD4_Naive = err_std, CD4_Mem = err_std, CD8_Naive = err_std, CD8_Mem = err_std,
      B_Naive = err_std, B_Mem = err_std, Treg = err_std, NK = err_std
    ) / scale
  } else if (type == "per_cell") {
    rmse_vector <- c(
      Neut = 2.85, Mono = 1.71, Eos = 0.95, Baso = 0.69,
      CD4_Naive = 2.47, CD4_Mem = 2.47, CD8_Naive = 1.98, CD8_Mem = 1.98,
      B_Naive = 0.54, B_Mem = 0.54, Treg = 0.59, NK = 2.27
    ) / scale
  }
  
  if (!is.null(seed)) set.seed(seed)
  W <- as.matrix(W)
  if (!all(colnames(W) %in% names(rmse_vector))) {
    stop("Every column of W must have a matching entry in rmse_vector")
  }
  rmse_vector <- rmse_vector[colnames(W)]  # align order to W's columns exactly
  
  N <- nrow(W); K <- ncol(W)
  
  noise <- matrix(
    rnorm(N * K, mean = 0, sd = rep(rmse_vector, each = N)),
    nrow = N, ncol = K
  )
  
  W_noisy <- W + noise
  W_noisy[W_noisy < floor] <- floor        # clip negative / near-zero perturbed proportions
  W_noisy <- W_noisy / rowSums(W_noisy)    # renormalize so every row is a valid composition (sums to 1)
  colnames(W_noisy) <- colnames(W)
  
  W_noisy
}

# Modified Simulation Function
Sim.Function.suppsim2 <- function(data_sim=data_sim,N=N,B1=B1,depth=depth,cpg_num=cpg_num,std_err_i=std_err_i){
  Data <- generate.data(data_sim=data_sim,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
  X_var = as.matrix(Data$X) ; W_var = as.matrix(Data$W) ; Y_var = as.matrix(Data$Y)
  W_noisy = add_deconvolution_error(W_var,err_std=std_err_i)
  # LR Selection
  LR_Selection_Result = LR_Selection.v3(X=Data$X,W=W_noisy,var.total=(ncol(W_noisy)-1))
  lr_covariate = LR_Selection_Result[[1]]
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
  limma_Covars = as.matrix(cbind(X_var,W_noisy))
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
  results = data.frame("Effect"=B1,"Sample Size"=N,"Error"=std_err_i,
                       "CLARE"=CLARE,"Proportion"=Proportion,"Unadjusted"=Unadjusted,
                       "CLARE_Bias"=CLARE_Bias,"CLARE_RMSE"=CLARE_RMSE,"CLARE_Variance"=CLARE_Variance,
                       "Proportion_Bias"=Proportion_Bias,"Proportion_RMSE"=Proportion_RMSE,"Proportion_Variance"=Proportion_Variance,
                       "Unadjusted_Bias"=Unadjusted_Bias,"Unadjusted_RMSE"=Unadjusted_RMSE,"Unadjusted_Variance"=Unadjusted_Variance)
  return(results)
}

## Sim Parameters ----

### Basic Parameters
sample_size=c(50,100,150,200,250,300,400,500) # Sample Size
Nsims= 50 # Outer Replicates - Number of Unique Cell Mixture Samples
cpg_num=500 # Inner Replicates - Number of CpGs generated per Cell Mixture Sample
total.cells=1000 # The total cells used for measuring DNAm
depth=200 # The probe depth
B1_null = 0 # Null Main Effect
B1_nonnull = -0.025 # Non-Null Main Effect
err_std_vector = c(0.5,1.0,1.5) # Vector of Deconvolution Errors (uniform, percentage-point scale)

#### cDMP Parameters, Neut/CD4_Naive confounding (W correlated with Main Effect)
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-10,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+5,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
cDMP_neut_cd4_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#### cDMP Parameters, Baso/B_Naive confounding (W correlated with Main Effect)
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=0,"Mono"=0,"Eos"=0,"Baso"=-0.2,"CD4_Naive"=0,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=+1,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
cDMP_baso_bnaive_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#### ncDMP Parameters, Neut/CD4_Naive confounding (W correlated with Main Effect, no methylation effect)
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-10,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+5,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.50,"Mono"=0.50,"Eos"=0.50,"Baso"=0.50,"CD4_Naive"=0.50,"CD4_Mem"=0.50,"CD8_Naive"=0.50,"CD8_Mem"=0.50,"B_Naive"=0.50,"B_Mem"=0.50,"Treg"=0.50,"NK"=0.50)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
ncDMP_neut_cd4_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#### ncDMP Parameters, Baso/B_Naive confounding (W correlated with Main Effect, no methylation effect)
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=0,"Mono"=0,"Eos"=0,"Baso"=-0.2,"CD4_Naive"=0,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=+1,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.50,"Mono"=0.50,"Eos"=0.50,"Baso"=0.50,"CD4_Naive"=0.50,"CD4_Mem"=0.50,"CD8_Naive"=0.50,"CD8_Mem"=0.50,"B_Naive"=0.50,"B_Mem"=0.50,"Treg"=0.50,"NK"=0.50)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
ncDMP_baso_bnaive_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#####################
## Simulation: Simulation 1 -----
### Set.Seed
set.seed(999848)

## 1. DMPs (W correlated with X)

### Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = cDMP_neut_cd4_parameters   # swap to cDMP_baso_bnaive_parameters to run the other scenario

for(err_i in 1:length(err_std_vector)){
  std_err_i = err_std_vector[err_i]
  cat(paste0("Error=", std_err_i, "\n"))
  ### Type I Error Rate
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_null,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
  ### Statistical Power
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_nonnull,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
}

### Packing
Results_DMP_Corr = Sim_Results
rm(Sim_Results,Scenario_params)

## 2. Non-DMPs (W correlated with X) 

### Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = ncDMP_neut_cd4_parameters   # swap to ncDMP_baso_bnaive_parameters to run the other scenario


for(err_i in 1:length(err_std_vector)){
  std_err_i = err_std_vector[err_i]
  cat(paste0("Error=", std_err_i, "\n"))
  ### Type I Error Rate
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_null,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
  
  ### Statistical Power
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_nonnull,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
}

### Packing
Results_NonDMP_Corr = Sim_Results
rm(Sim_Results)

## 3. Packaging & Saving

# Notes for Individual Scenarios
Results_DMP_Corr$Scenario = "DMP Correlated"
Results_NonDMP_Corr$Scenario = "Non-DMP Correlated"

# Results Simple-ish
res_wide = bind_rows(Results_DMP_Corr,Results_NonDMP_Corr)%>%
  dplyr::filter(!is.na(Effect))

# Complete Data Frame
Results_Final_SuppSim2 = bind_rows(Results_DMP_Corr
                                   ,Results_NonDMP_Corr
)%>%
  mutate(Effect=case_when(Effect==0~"Type I Error",TRUE~"Power"))%>%
  rename(CLARE_Detection=CLARE,Proportion_Detection=Proportion,Unadjusted_Detection=Unadjusted)%>%
  dplyr::filter(!is.na(Sample.Size))%>%
  pivot_longer(cols = matches("_Bias$|_RMSE$|_Detection$|_Variance$"), names_to = c("Method","Metric"), names_sep = "_")

##################### 
### Summaries -----
### Double Checking RMSE Capture
max(abs(as.numeric(res_wide$CLARE_RMSE)^2 - (as.numeric(res_wide$CLARE_Bias)^2 + as.numeric(res_wide$CLARE_Variance))))
max(abs(as.numeric(res_wide$Unadjusted_RMSE)^2 - (as.numeric(res_wide$Unadjusted_Bias)^2 + as.numeric(res_wide$Unadjusted_Variance))))
max(abs(as.numeric(res_wide$Proportion_RMSE)^2 - (as.numeric(res_wide$Proportion_Bias)^2 + as.numeric(res_wide$Proportion_Variance))))



### Bias/RMSE/Variance

metric_means <- Results_Final_SuppSim2 %>%
  mutate(value=as.numeric(value))%>%
  group_by(Effect,Scenario,Error,Method,Metric)%>%
  summarise(Metric_Mean=mean(value,na.rm=TRUE))%>%
  arrange(Metric,Scenario,Error,Effect,Method)

#####################
## Quick Figures ----
## Overall prep ----
#
colors <- c("#15bf3a",
            "#4075d6",
            "#631307")
shapes <- c(15,1,4)
sample.size.breaks=unique(Results_Final_SuppSim2$Sample.Size)
Simulation.title = "Supplementary Simulation 2: "
# subtitle Prefix
if(Scenario_params$Lambda0[1]==Scenario_params$Lambda1[1]){
  Simulation.title.suffix="Basophil/B Naive"
}else if(Scenario_params$Lambda0[4]==Scenario_params$Lambda1[4]){
  Simulation.title.suffix="Neutrophils/CD4 Naive"
}
Simulation.title = paste0(Simulation.title,Simulation.title.suffix)
#
unique(Results_Final_SuppSim2$Scenario)
df <-  Results_Final_SuppSim2 %>%
  mutate(Scenario=case_when(Scenario=="DMP Correlated"~"cDMP",Scenario=="Non-DMP Correlated"~"ncDMP",TRUE~NA))%>%
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Unadjusted")))%>%
  mutate(Effect=factor(Effect,levels=c("Type I Error","Power")))%>%
  mutate(Error_Label=paste0("Error=",Error))%>%
  mutate(Error_Label=factor(Error_Label,levels=paste0("Error=",sort(unique(Error)))))

## B1. Mean Power & Type I Error Rate ----
df_graph = df %>%
  dplyr::filter(Method != "Unadjusted" & !is.na(Scenario) & Metric == "Detection")%>%
  group_by(Effect,Sample.Size,Scenario,Error_Label,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE))%>%
  ungroup()
#
plotB1 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,linetype=Error_Label))+
  geom_line(size=1,alpha=0.6)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Power & Type I Error Rate, by Deconvolution Error Level",
       x="Sample Size",
       y="Detection Rate",
       linetype="Error")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="fixed")+
  scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1),minor_breaks=NULL)+
  geom_hline(yintercept=0.05,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB1

## B2. Mean Bias ----
df_graph = df %>%
  mutate(Effect=case_when(Effect=="Type I Error"~"No Effect",
                          Effect=="Power"~"True Effect",TRUE~NA))%>%
  mutate(Effect=factor(Effect,levels=c("No Effect","True Effect")))%>%
  dplyr::filter(!is.na(Scenario) & Metric == "Bias")%>%
  group_by(Effect,Sample.Size,Scenario,Error_Label,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE))%>%
  ungroup()
#
plotB2 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,linetype=Error_Label))+
  geom_line(size=1,alpha=0.6)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Mean Bias in Estimated Effect, by Deconvolution Error Level",
       x="Sample Size",
       y=expression(paste("Mean Bias (", hat(beta)[1], " - ", beta[1], ")")),
       linetype="Error")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="free_y")+
  geom_hline(yintercept=0,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB2
#################################################################################################################
## SUPPLEMENTARY SIMULATION 2.2 - Deconvolution Error - Basophils/B Naive -----

## Functions -----
# Function for adding Deconvolution Error
add_deconvolution_error <- function(W, floor = 1e-4, seed = NULL, scale = 100, type = "standard", err_std = 0.5) {
  if (type == "standard") {
    rmse_vector <- c(
      Neut = err_std, Mono = err_std, Eos = err_std, Baso = err_std,
      CD4_Naive = err_std, CD4_Mem = err_std, CD8_Naive = err_std, CD8_Mem = err_std,
      B_Naive = err_std, B_Mem = err_std, Treg = err_std, NK = err_std
    ) / scale
  } else if (type == "per_cell") {
    rmse_vector <- c(
      Neut = 2.85, Mono = 1.71, Eos = 0.95, Baso = 0.69,
      CD4_Naive = 2.47, CD4_Mem = 2.47, CD8_Naive = 1.98, CD8_Mem = 1.98,
      B_Naive = 0.54, B_Mem = 0.54, Treg = 0.59, NK = 2.27
    ) / scale
  }
  
  if (!is.null(seed)) set.seed(seed)
  W <- as.matrix(W)
  if (!all(colnames(W) %in% names(rmse_vector))) {
    stop("Every column of W must have a matching entry in rmse_vector")
  }
  rmse_vector <- rmse_vector[colnames(W)]  # align order to W's columns exactly
  
  N <- nrow(W); K <- ncol(W)
  
  noise <- matrix(
    rnorm(N * K, mean = 0, sd = rep(rmse_vector, each = N)),
    nrow = N, ncol = K
  )
  
  W_noisy <- W + noise
  W_noisy[W_noisy < floor] <- floor        # clip negative / near-zero perturbed proportions
  W_noisy <- W_noisy / rowSums(W_noisy)    # renormalize so every row is a valid composition (sums to 1)
  colnames(W_noisy) <- colnames(W)
  
  W_noisy
}

# Modified Simulation Function
Sim.Function.suppsim2 <- function(data_sim=data_sim,N=N,B1=B1,depth=depth,cpg_num=cpg_num,std_err_i=std_err_i){
  Data <- generate.data(data_sim=data_sim,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
  X_var = as.matrix(Data$X) ; W_var = as.matrix(Data$W) ; Y_var = as.matrix(Data$Y)
  W_noisy = add_deconvolution_error(W_var,err_std=std_err_i)
  # LR Selection
  LR_Selection_Result = LR_Selection.v3(X=Data$X,W=W_noisy,var.total=(ncol(W_noisy)-1))
  lr_covariate = LR_Selection_Result[[1]]
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
  limma_Covars = as.matrix(cbind(X_var,W_noisy))
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
  results = data.frame("Effect"=B1,"Sample Size"=N,"Error"=std_err_i,
                       "CLARE"=CLARE,"Proportion"=Proportion,"Unadjusted"=Unadjusted,
                       "CLARE_Bias"=CLARE_Bias,"CLARE_RMSE"=CLARE_RMSE,"CLARE_Variance"=CLARE_Variance,
                       "Proportion_Bias"=Proportion_Bias,"Proportion_RMSE"=Proportion_RMSE,"Proportion_Variance"=Proportion_Variance,
                       "Unadjusted_Bias"=Unadjusted_Bias,"Unadjusted_RMSE"=Unadjusted_RMSE,"Unadjusted_Variance"=Unadjusted_Variance)
  return(results)
}

## Sim Parameters ----

### Basic Parameters
sample_size=c(50,100,150,200,250,300,400,500) # Sample Size
Nsims= 50 # Outer Replicates - Number of Unique Cell Mixture Samples
cpg_num=500 # Inner Replicates - Number of CpGs generated per Cell Mixture Sample
total.cells=1000 # The total cells used for measuring DNAm
depth=200 # The probe depth
B1_null = 0 # Null Main Effect
B1_nonnull = -0.025 # Non-Null Main Effect
err_std_vector = c(0.5,1.0,1.5) # Vector of Deconvolution Errors (uniform, percentage-point scale)

#### cDMP Parameters, Neut/CD4_Naive confounding (W correlated with Main Effect)
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-10,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+5,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
cDMP_neut_cd4_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#### cDMP Parameters, Baso/B_Naive confounding (W correlated with Main Effect)
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=0,"Mono"=0,"Eos"=0,"Baso"=-0.2,"CD4_Naive"=0,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=+1,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
cDMP_baso_bnaive_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#### ncDMP Parameters, Neut/CD4_Naive confounding (W correlated with Main Effect, no methylation effect)
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-10,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+5,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.50,"Mono"=0.50,"Eos"=0.50,"Baso"=0.50,"CD4_Naive"=0.50,"CD4_Mem"=0.50,"CD8_Naive"=0.50,"CD8_Mem"=0.50,"B_Naive"=0.50,"B_Mem"=0.50,"Treg"=0.50,"NK"=0.50)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
ncDMP_neut_cd4_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#### ncDMP Parameters, Baso/B_Naive confounding (W correlated with Main Effect, no methylation effect)
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=0,"Mono"=0,"Eos"=0,"Baso"=-0.2,"CD4_Naive"=0,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=+1,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.50,"Mono"=0.50,"Eos"=0.50,"Baso"=0.50,"CD4_Naive"=0.50,"CD4_Mem"=0.50,"CD8_Naive"=0.50,"CD8_Mem"=0.50,"B_Naive"=0.50,"B_Mem"=0.50,"Treg"=0.50,"NK"=0.50)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
ncDMP_baso_bnaive_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#####################
# Simulation: Simulation 1 -----

### Set.Seed()
set.seed(999848)

## 1. DMPs (W correlated with X)

### Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = cDMP_baso_bnaive_parameters   # swap to cDMP_baso_bnaive_parameters to run the other scenario
#
for(err_i in 1:length(err_std_vector)){
  std_err_i = err_std_vector[err_i]
  cat(paste0("Error=", std_err_i, "\n"))
  ### Type I Error Rate
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_null,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
  ### Statistical Power
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_nonnull,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
}

### Packing
Results_DMP_Corr = Sim_Results
rm(Sim_Results,Scenario_params)

## 2. Non-DMPs (W correlated with X)

### Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = ncDMP_baso_bnaive_parameters   # swap to ncDMP_baso_bnaive_parameters to run the other scenario
#
for(err_i in 1:length(err_std_vector)){
  std_err_i = err_std_vector[err_i]
  cat(paste0("Error=", std_err_i, "\n"))
  ### Type I Error Rate
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_null,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
  
  ### Statistical Power
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_nonnull,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
}

### Packing
Results_NonDMP_Corr = Sim_Results
rm(Sim_Results)

## 3. Packaging & Saving

# Notes for Individual Scenarios
Results_DMP_Corr$Scenario = "DMP Correlated"
Results_NonDMP_Corr$Scenario = "Non-DMP Correlated"

# Results Simple-ish
res_wide = bind_rows(Results_DMP_Corr,Results_NonDMP_Corr)%>%
  dplyr::filter(!is.na(Effect))

# Complete Data Frame
Results_Final_SuppSim2 = bind_rows(Results_DMP_Corr
                                   ,Results_NonDMP_Corr
)%>%
  mutate(Effect=case_when(Effect==0~"Type I Error",TRUE~"Power"))%>%
  rename(CLARE_Detection=CLARE,Proportion_Detection=Proportion,Unadjusted_Detection=Unadjusted)%>%
  dplyr::filter(!is.na(Sample.Size))%>%
  pivot_longer(cols = matches("_Bias$|_RMSE$|_Detection$|_Variance$"), names_to = c("Method","Metric"), names_sep = "_")

#####################
### Double Checking RMSE Capture ----
max(abs(as.numeric(res_wide$CLARE_RMSE)^2 - (as.numeric(res_wide$CLARE_Bias)^2 + as.numeric(res_wide$CLARE_Variance))))
max(abs(as.numeric(res_wide$Unadjusted_RMSE)^2 - (as.numeric(res_wide$Unadjusted_Bias)^2 + as.numeric(res_wide$Unadjusted_Variance))))
max(abs(as.numeric(res_wide$Proportion_RMSE)^2 - (as.numeric(res_wide$Proportion_Bias)^2 + as.numeric(res_wide$Proportion_Variance))))


### Summaries -----
### Bias/RMSE/Variance

metric_means <- Results_Final_SuppSim2 %>%
  mutate(value=as.numeric(value))%>%
  group_by(Effect,Scenario,Error,Method,Metric)%>%
  summarise(Metric_Mean=mean(value,na.rm=TRUE))%>%
  arrange(Metric,Scenario,Error,Effect,Method)

#####################
## Quick Figures ----
## Overall prep ----
#
colors <- c("#15bf3a",
            "#4075d6",
            "#631307")
shapes <- c(15,1,4)
sample.size.breaks=unique(Results_Final_SuppSim2$Sample.Size)
Simulation.title = "Supplementary Simulation 2: "
# subtitle Prefix
if(Scenario_params$Lambda0[1]==Scenario_params$Lambda1[1]){
  Simulation.title.suffix="Basophil/B Naive"
}else if(Scenario_params$Lambda0[4]==Scenario_params$Lambda1[4]){
  Simulation.title.suffix="Neutrophils/CD4 Naive"
}
Simulation.title = paste0(Simulation.title,Simulation.title.suffix)
#
unique(Results_Final_SuppSim2$Scenario)
df <-  Results_Final_SuppSim2 %>%
  mutate(Scenario=case_when(Scenario=="DMP Correlated"~"cDMP",Scenario=="Non-DMP Correlated"~"ncDMP",TRUE~NA))%>%
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Unadjusted")))%>%
  mutate(Effect=factor(Effect,levels=c("Type I Error","Power")))%>%
  mutate(Error_Label=paste0("Error=",Error))%>%
  mutate(Error_Label=factor(Error_Label,levels=paste0("Error=",sort(unique(Error)))))

## B1. Mean Power & Type I Error Rate ----
df_graph = df %>%
  dplyr::filter(Method != "Unadjusted" & !is.na(Scenario) & Metric == "Detection")%>%
  group_by(Effect,Sample.Size,Scenario,Error_Label,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE))%>%
  ungroup()
#
plotB1 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,linetype=Error_Label))+
  geom_line(size=1,alpha=0.6)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Power & Type I Error Rate, by Deconvolution Error Level",
       x="Sample Size",
       y="Detection Rate",
       linetype="Error")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="fixed")+
  scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1),minor_breaks=NULL)+
  geom_hline(yintercept=0.05,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB1

## B2. Mean Bias ----
df_graph = df %>%
  mutate(Effect=case_when(Effect=="Type I Error"~"No Effect",
                          Effect=="Power"~"True Effect",TRUE~NA))%>%
  mutate(Effect=factor(Effect,levels=c("No Effect","True Effect")))%>%
  dplyr::filter(!is.na(Scenario) & Metric == "Bias")%>%
  group_by(Effect,Sample.Size,Scenario,Error_Label,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE))%>%
  ungroup()
#
plotB2 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,linetype=Error_Label))+
  geom_line(size=1,alpha=0.6)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Mean Bias in Estimated Effect, by Deconvolution Error Level",
       x="Sample Size",
       y=expression(paste("Mean Bias (", hat(beta)[1], " - ", beta[1], ")")),
       linetype="Error")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="free_y")+
  geom_hline(yintercept=0,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB2



#################################################################################################################
## SUPPLEMENTARY SIMULATION 2.3 - Deconvolution Error - Even Cell Types -----
## Functions -----
# Function for adding Deconvolution Error
add_deconvolution_error <- function(W, floor = 1e-4, seed = NULL, scale = 100, type = "standard", err_std = 0.5) {
  if (type == "standard") {
    rmse_vector <- c(
      Neut = err_std, Mono = err_std, Eos = err_std, Baso = err_std,
      CD4_Naive = err_std, CD4_Mem = err_std, CD8_Naive = err_std, CD8_Mem = err_std,
      B_Naive = err_std, B_Mem = err_std, Treg = err_std, NK = err_std
    ) / scale
  } else if (type == "per_cell") {
    rmse_vector <- c(
      Neut = 2.85, Mono = 1.71, Eos = 0.95, Baso = 0.69,
      CD4_Naive = 2.47, CD4_Mem = 2.47, CD8_Naive = 1.98, CD8_Mem = 1.98,
      B_Naive = 0.54, B_Mem = 0.54, Treg = 0.59, NK = 2.27
    ) / scale
  }
  
  if (!is.null(seed)) set.seed(seed)
  W <- as.matrix(W)
  if (!all(colnames(W) %in% names(rmse_vector))) {
    stop("Every column of W must have a matching entry in rmse_vector")
  }
  rmse_vector <- rmse_vector[colnames(W)]  # align order to W's columns exactly
  
  N <- nrow(W); K <- ncol(W)
  
  noise <- matrix(
    rnorm(N * K, mean = 0, sd = rep(rmse_vector, each = N)),
    nrow = N, ncol = K
  )
  
  W_noisy <- W + noise
  W_noisy[W_noisy < floor] <- floor        # clip negative / near-zero perturbed proportions
  W_noisy <- W_noisy / rowSums(W_noisy)    # renormalize so every row is a valid composition (sums to 1)
  colnames(W_noisy) <- colnames(W)
  
  W_noisy
}

# Modified Simulation Function
Sim.Function.suppsim2 <- function(data_sim=data_sim,N=N,B1=B1,depth=depth,cpg_num=cpg_num,std_err_i=std_err_i){
  Data <- generate.data(data_sim=data_sim,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
  X_var = as.matrix(Data$X) ; W_var = as.matrix(Data$W) ; Y_var = as.matrix(Data$Y)
  W_noisy = add_deconvolution_error(W_var,err_std=std_err_i)
  # LR Selection
  LR_Selection_Result = LR_Selection.v3(X=Data$X,W=W_noisy,var.total=(ncol(W_noisy)-1))
  lr_covariate = LR_Selection_Result[[1]]
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
  limma_Covars = as.matrix(cbind(X_var,W_noisy))
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
  results = data.frame("Effect"=B1,"Sample Size"=N,"Error"=std_err_i,
                       "CLARE"=CLARE,"Proportion"=Proportion,"Unadjusted"=Unadjusted,
                       "CLARE_Bias"=CLARE_Bias,"CLARE_RMSE"=CLARE_RMSE,"CLARE_Variance"=CLARE_Variance,
                       "Proportion_Bias"=Proportion_Bias,"Proportion_RMSE"=Proportion_RMSE,"Proportion_Variance"=Proportion_Variance,
                       "Unadjusted_Bias"=Unadjusted_Bias,"Unadjusted_RMSE"=Unadjusted_RMSE,"Unadjusted_Variance"=Unadjusted_Variance)
  return(results)
}

## Sim Parameters ----

### Basic Parameters
sample_size=c(50,100,150,200,250,300,400,500) # Sample Size
Nsims= 50 # Outer Replicates - Number of Unique Cell Mixture Samples
cpg_num=500 # Inner Replicates - Number of CpGs generated per Cell Mixture Sample
total.cells=1000 # The total cells used for measuring DNAm
depth=200 # The probe depth
B1_null = 0 # Null Main Effect
B1_nonnull = -0.025 # Non-Null Main Effect
err_std_vector = c(0.5,1.0,1.5) # Vector of Deconvolution Errors (uniform, percentage-point scale)


#### cDMP Parameters, 2 DA
lambda0=c("Neut"=10,"Mono"=10,"Eos"=10,"Baso"=10,"CD4_Naive"=10,"CD4_Mem"=10,"CD8_Naive"=10,"CD8_Mem"=10,"B_Naive"=10,"B_Mem"=10,"Treg"=10,"NK"=10)
tau = c("Neut"=-1,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+1,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.75,"Mono"=0.5,"Eos"=0.25,"Baso"=0.5,"CD4_Naive"=0.25,"CD4_Mem"=0.25,"CD8_Naive"=0.25,"CD8_Mem"=0.75,"B_Naive"=0.25,"B_Mem"=0.5,"Treg"=0.5,"NK"=0.5)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
cDMP_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

# Checking Expected Base Methylation
# Emu0 = sum((lambda0/sum(lambda0))*B0)
# Emu1 = sum((lambda1/sum(lambda1))*B0)

#### ncDMP Parameters, 2 DA
lambda0=c("Neut"=10,"Mono"=10,"Eos"=10,"Baso"=10,"CD4_Naive"=10,"CD4_Mem"=10,"CD8_Naive"=10,"CD8_Mem"=10,"B_Naive"=10,"B_Mem"=10,"Treg"=10,"NK"=10)
tau = c("Neut"=-1,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+1,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.50,"Mono"=0.50,"Eos"=0.50,"Baso"=0.50,"CD4_Naive"=0.50,"CD4_Mem"=0.50,"CD8_Naive"=0.50,"CD8_Mem"=0.50,"B_Naive"=0.50,"B_Mem"=0.50,"Treg"=0.50,"NK"=0.50)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
ncDMP_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#####################
# Simulation: Standard Cell Abundances -----
### Set.Seed()
set.seed(999848)

## 1. DMPs (W correlated with X)

### Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = cDMP_parameters   # swap to cDMP_baso_bnaive_parameters to run the other scenario

for(err_i in 1:length(err_std_vector)){
  std_err_i = err_std_vector[err_i]
  cat(paste0("Error=", std_err_i, "\n"))
  ### Type I Error Rate
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_null,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
  ### Statistical Power
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_nonnull,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
}

### Packing
Results_DMP_Corr = Sim_Results
rm(Sim_Results,Scenario_params)

## 2. Non-DMPs (W correlated with X)

### Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = ncDMP_parameters   # swap to ncDMP_baso_bnaive_parameters to run the other scenario


for(err_i in 1:length(err_std_vector)){
  std_err_i = err_std_vector[err_i]
  cat(paste0("Error=", std_err_i, "\n"))
  ### Type I Error Rate
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_null,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
  
  ### Statistical Power
  for(N_i in 1:length(sample_size)){
    results_i = data.frame("Effect"=NA,"Sample Size"=NA,"Error"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
    for(sim in 1:Nsims){
      results_i = bind_rows(results_i,Sim.Function.suppsim2(data_sim=Scenario_params,
                                                            N=sample_size[N_i],
                                                            B1=B1_nonnull,
                                                            depth=depth,
                                                            cpg_num=cpg_num,
                                                            std_err_i=std_err_i))
    }
    cat(paste0(sample_size[N_i], ","))
    results_i = results_i[!is.na(results_i$Effect),]
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat("\n")
}

### Packing
Results_NonDMP_Corr = Sim_Results
rm(Sim_Results)

## 3. Packaging & Saving

# Notes for Individual Scenarios
Results_DMP_Corr$Scenario = "DMP Correlated"
Results_NonDMP_Corr$Scenario = "Non-DMP Correlated"

# Results Simple-ish
res_wide = bind_rows(Results_DMP_Corr,Results_NonDMP_Corr)%>%
  dplyr::filter(!is.na(Effect))

# Complete Data Frame
Results_Final_SuppSim2 = bind_rows(Results_DMP_Corr
                                   ,Results_NonDMP_Corr
)%>%
  mutate(Effect=case_when(Effect==0~"Type I Error",TRUE~"Power"))%>%
  rename(CLARE_Detection=CLARE,Proportion_Detection=Proportion,Unadjusted_Detection=Unadjusted)%>%
  dplyr::filter(!is.na(Sample.Size))%>%
  pivot_longer(cols = matches("_Bias$|_RMSE$|_Detection$|_Variance$"), names_to = c("Method","Metric"), names_sep = "_")

#####################
### Double Checking RMSE Capture ----
max(abs(as.numeric(res_wide$CLARE_RMSE)^2 - (as.numeric(res_wide$CLARE_Bias)^2 + as.numeric(res_wide$CLARE_Variance))))
max(abs(as.numeric(res_wide$Unadjusted_RMSE)^2 - (as.numeric(res_wide$Unadjusted_Bias)^2 + as.numeric(res_wide$Unadjusted_Variance))))
max(abs(as.numeric(res_wide$Proportion_RMSE)^2 - (as.numeric(res_wide$Proportion_Bias)^2 + as.numeric(res_wide$Proportion_Variance))))


### Summaries -----
### Bias/RMSE/Variance

metric_means <- Results_Final_SuppSim2 %>%
  mutate(value=as.numeric(value))%>%
  group_by(Effect,Scenario,Error,Method,Metric)%>%
  summarise(Metric_Mean=mean(value,na.rm=TRUE))%>%
  arrange(Metric,Scenario,Error,Effect,Method)

#####################
## Quick Figures ----
## Overall prep ----
#
colors <- c("#15bf3a",
            "#4075d6",
            "#631307")
shapes <- c(15,1,4)
sample.size.breaks=unique(Results_Final_SuppSim2$Sample.Size)
Simulation.title = "Supplementary Simulation 2: "
# subtitle Prefix
Simulation.title.suffix="Even Cell Type Abundances"
Simulation.title = paste0(Simulation.title,Simulation.title.suffix)
#
unique(Results_Final_SuppSim2$Scenario)
df <-  Results_Final_SuppSim2 %>%
  mutate(Scenario=case_when(Scenario=="DMP Correlated"~"cDMP",Scenario=="Non-DMP Correlated"~"ncDMP",TRUE~NA))%>%
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Unadjusted")))%>%
  mutate(Effect=factor(Effect,levels=c("Type I Error","Power")))%>%
  mutate(Error_Label=paste0("Error=",Error))%>%
  mutate(Error_Label=factor(Error_Label,levels=paste0("Error=",sort(unique(Error)))))

## B1. Mean Power & Type I Error Rate ----
df_graph = df %>%
  dplyr::filter(!is.na(Scenario) & Metric == "Detection" 
                & Method != "Unadjusted" 
  )%>%
  group_by(Effect,Sample.Size,Scenario,Error_Label,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE))%>%
  ungroup()
#
plotB1 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,linetype=Error_Label))+
  geom_line(size=1,alpha=0.6)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Power & Type I Error Rate, by Deconvolution Error Level",
       x="Sample Size",
       y="Detection Rate",
       linetype="Error")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="fixed")+
  scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1),minor_breaks=NULL)+
  geom_hline(yintercept=0.05,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB1

## B2. Mean Bias ----
df_graph = df %>%
  mutate(Effect=case_when(Effect=="Type I Error"~"No Effect",
                          Effect=="Power"~"True Effect",TRUE~NA))%>%
  mutate(Effect=factor(Effect,levels=c("No Effect","True Effect")))%>%
  dplyr::filter(!is.na(Scenario) & Metric == "Bias")%>%
  group_by(Effect,Sample.Size,Scenario,Error_Label,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE))%>%
  ungroup()
#
plotB2 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,linetype=Error_Label))+
  geom_line(size=1,alpha=0.6)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Mean Bias in Estimated Effect, by Deconvolution Error Level",
       x="Sample Size",
       y=expression(paste("Mean Bias (", hat(beta)[1], " - ", beta[1], ")")),
       linetype="Error")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Effect,scales="free_y")+
  geom_hline(yintercept=0,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB2