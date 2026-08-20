#################################################################################################################
### Title: Adjusted EWAS Using Compositional Data Transformations - Supplementary Simulation #3 ----
### Current Analyst: Alexander Alsup (a752a465@kumc.edu) ----
### Last Updated: (08/14/2026) ----
### Notes: ----
#################################################################################################################
# Source Functions & Parent Filepath ----
# Replace with your own directory where applicable
Parent_Directory = ""
source(paste0(Parent_Directory,"Supplementary Sim Prep.R"))
#################################################################################################################
## SUPPLEMENTARY SIMULATION 3 - VARIABLE CONFOUNDING EFFECTS -----
## Functions ----
Sim.Function.suppsim3 <- function(data_sim=data_sim,N=N,B1=B1,depth=depth,cpg_num=cpg_num){
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
  results = data.frame("Effect"=B1,"Sample Size"=N,"Theta_M"=data_sim$Theta_M,
                       "CLARE"=CLARE,"Proportion"=Proportion,"Unadjusted"=Unadjusted,
                       "CLARE_Bias"=CLARE_Bias,"CLARE_RMSE"=CLARE_RMSE,"CLARE_Variance"=CLARE_Variance,
                       "Proportion_Bias"=Proportion_Bias,"Proportion_RMSE"=Proportion_RMSE,"Proportion_Variance"=Proportion_Variance,
                       "Unadjusted_Bias"=Unadjusted_Bias,"Unadjusted_RMSE"=Unadjusted_RMSE,"Unadjusted_Variance"=Unadjusted_Variance)
  return(results)
}

## Sim Parameters ----

# Beta_1s
B1_null = 0 # Null Main Effect
B1_nonnull = -0.025 # Non-Null Main Effect


### cDMP: ΘM ≈ -0.02 ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-10,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+5,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.6,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.3,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
Theta_M = sum((lambda1/sum(lambda1))*B0)-sum((lambda0/sum(lambda0))*B0);Theta_M
#
cDMP_Theta02 <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0,"Theta_M"=Theta_M)

### cDMP: ΘM ≈ -0.04 ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-10,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+5,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
Theta_M = sum((lambda1/sum(lambda1))*B0)-sum((lambda0/sum(lambda0))*B0)
#
cDMP_Theta04 <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0,"Theta_M"=Theta_M)

### cDMP: ΘM ≈ -0.06 ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-10,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+5,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.93,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.07,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
Theta_M = sum((lambda1/sum(lambda1))*B0)-sum((lambda0/sum(lambda0))*B0);Theta_M
#
cDMP_Theta06 <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0,"Theta_M"=Theta_M)

##############################################
# Simulation:   -----
### Set.Seed()
set.seed(999848)
Sim_Results = data.frame()

## Simulation: ΘM ≈ -0.02 ----
# Scenario Parameters
sim_params = cDMP_Theta02
lambda0= sim_params$lambda0 ;lambda1=sim_params$lambda1 ; B0=sim_params$B0 

# Simulation Loop
for(N_i in 1:length(sample_size)){
  results_i <- data.frame()
  N=sample_size[N_i]
  # Inner Loop for Type I Error
  B1=B1_null
  for(sim in 1:Nsims){
    results_init <- Sim.Function.suppsim3(data_sim=sim_params,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
    results_i = bind_rows(results_i,results_init)
  }
  # Inner Loop for Power
  B1=B1_nonnull
  for(sim in 1:Nsims){
    results_init <- Sim.Function.suppsim3(data_sim=sim_params,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
    results_i = bind_rows(results_i,results_init)
  }
  # Outer Loop Message Update
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}
cat("\n")
## Simulation: ΘM ≈ -0.04 ----
# Scenario Parameters
sim_params = cDMP_Theta04
lambda0= sim_params$lambda0 ;lambda1=sim_params$lambda1 ; B0=sim_params$B0 

# Simulation Loop
for(N_i in 1:length(sample_size)){
  results_i <- data.frame()
  N=sample_size[N_i]
  # Inner Loop for Type I Error
  B1=B1_null
  for(sim in 1:Nsims){
    results_init <- Sim.Function.suppsim3(data_sim=sim_params,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
    results_i = bind_rows(results_i,results_init)
  }
  # Inner Loop for Power
  B1=B1_nonnull
  for(sim in 1:Nsims){
    results_init <- Sim.Function.suppsim3(data_sim=sim_params,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
    results_i = bind_rows(results_i,results_init)
  }
  # Outer Loop Message Update
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}
cat("\n")
## Simulation: ΘM ≈ -0.06 ----
# Scenario Parameters
sim_params = cDMP_Theta06
lambda0= sim_params$lambda0 ;lambda1=sim_params$lambda1 ; B0=sim_params$B0 

# Simulation Loop
for(N_i in 1:length(sample_size)){
  results_i <- data.frame()
  N=sample_size[N_i]
  # Inner Loop for Type I Error
  B1=B1_null
  for(sim in 1:Nsims){
    results_init <- Sim.Function.suppsim3(data_sim=sim_params,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
    results_i = bind_rows(results_i,results_init)
  }
  # Inner Loop for Power
  B1=B1_nonnull
  for(sim in 1:Nsims){
    results_init <- Sim.Function.suppsim3(data_sim=sim_params,N=N,B1=B1,depth=depth,cpg_num=cpg_num)
    results_i = bind_rows(results_i,results_init)
  }
  # Outer Loop Message Update
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}
cat("\n")
##############################################
## 3. Packaging & Saving  -----

# Results Simple-ish
res_wide = Sim_Results%>%
  dplyr::filter(!is.na(Effect))

# Complete Data Frame
Results_Final_SuppSim = Sim_Results %>%
  mutate(Theta_M=round(Theta_M,2))%>%
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

metric_means <- Results_Final_SuppSim %>%
  mutate(value=as.numeric(value))%>%
  group_by(Effect,Theta_M,Method,Metric)%>%
  summarise(Metric_Mean=mean(value,na.rm=TRUE))%>%
  arrange(Metric,Theta_M,Effect,Method)

#####################
## Quick Figures ----
## Overall prep ----
#
colors <- c("#15bf3a",
            "#4075d6",
            "#631307")
shapes <- c(15,1,4)
sample.size.breaks=unique(Results_Final_SuppSim$Sample.Size)
Simulation.title = "Supplementary Simulation 3: "
#
df <-  Results_Final_SuppSim %>%
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Unadjusted")))%>%
  mutate(Effect=factor(Effect,levels=c("Type I Error","Power")))%>%
  mutate(Theta_M_Label=paste0("ΘM ≈",Theta_M))%>%
  mutate(Theta_M_Label=factor(Theta_M_Label,levels=paste0("ΘM ≈",sort(unique(Theta_M),decreasing=TRUE))))

## B1. Mean Power & Type I Error Rate ----
df_graph = df %>%
  dplyr::filter(!is.na(Theta_M_Label) & Metric == "Detection" 
                & Method != "Unadjusted" 
  )%>%
  group_by(Effect,Sample.Size,Theta_M_Label,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE))%>%
  ungroup()
#
plotB1 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,linetype=Theta_M_Label))+
  geom_line(size=1,alpha=0.6)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Power & Type I Error Rate, by Confounding Magnitude",
       x="Sample Size",
       y="Detection Rate",
       linetype="Confounding Magnitude")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Theta_M_Label~Effect,scales="fixed")+
  scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1),minor_breaks=NULL)+
  geom_hline(yintercept=0.05,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB1

## B2. Mean Bias ----
df_graph = df %>%
  dplyr::filter(!is.na(Theta_M_Label) & Metric == "Bias" 
                #& Method != "Unadjusted" 
  )%>%
  group_by(Effect,Sample.Size,Theta_M_Label,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE))%>%
  ungroup()
#
plotB2 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,linetype=Theta_M_Label))+
  geom_line(size=1,alpha=0.6)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Mean Bias in Estimated Effect, by Confounding Magnitude",
       x="Sample Size",
       y=expression(paste("Mean Bias (", hat(beta)[1], " - ", beta[1], ")")),
       linetype="Confounding Magnitude")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Theta_M_Label~Effect,scales="fixed")+
  geom_hline(yintercept=0,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB2