#################################################################################################################
### Title: Adjusted EWAS Using Compositional Data Transformations - SIMULATION 1 ----
### Current Analyst: Alexander Alsup (a752a465@kumc.edu) ----
### Last Updated: (01/12/2026) ----
### Notes: ----
#################################################################################################################
# PACKAGES -----
packages <- c("ggplot2","here","MASS","MCMCpack","limma","betareg","tidyverse","writexl","readxl","easyCODA","limma","dplyr")
lapply(packages, require, character.only = TRUE)
# Clear Environment ----
rm(list=ls())
# Source Functions & Parent Filepath ----
Parent_Directory = ""
source(paste0(Parent_Directory,"Functions.R"))
#################################################################################################################
# GLOBAL PARAMETERS -----
## Outputs ----
Outputs=TRUE
### Basic Parameters ----
sample_size=c(50,100,150,200,250,300,400,500) # Sample Size
Nsims= 50 # Outer Replicates - Number of Unique Cell Mixture Samples
cpg_num=500 # Inner Replicates - Number of CpGs generated per Cell Mixture Sample

#### Set.Seed() ----
set.seed(999848)


### Data Generation Parameters ----
total.cells=1000 # The total cells used for measuring DNAm
depth=200 # The probe depth
B1_null = 0 # Null Main Effect
B1_nonnull = -0.025 # Non-Null Main Effect

### DMP/non-DMP Simulation Parameters ----

#### DMP Parameters (W correlated with Main Effect) ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-10,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+5,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
DMP_corr_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#### DMP Parameters (W uncorrelated with Main Effect) ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=0,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=0,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
DMP_uncorr_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#### Non-DMP Parameters (W correlated with Main Effect) ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-10,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+5,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.50,"Mono"=0.50,"Eos"=0.50,"Baso"=0.50,"CD4_Naive"=0.50,"CD4_Mem"=0.50,"CD8_Naive"=0.50,"CD8_Mem"=0.50,"B_Naive"=0.50,"B_Mem"=0.50,"Treg"=0.50,"NK"=0.50)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
NonDMP_corr_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#### Non-DMP Parameters (W uncorrelated with Main Effect) ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=0,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=0,"CD4_Mem"=0,"CD8_Naive"=0,"CD8_Mem"=0,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.50,"Mono"=0.50,"Eos"=0.50,"Baso"=0.50,"CD4_Naive"=0.50,"CD4_Mem"=0.50,"CD8_Naive"=0.50,"CD8_Mem"=0.50,"B_Naive"=0.50,"B_Mem"=0.50,"Treg"=0.50,"NK"=0.50)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
NonDMP_uncorr_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)

#
rm(lambda0,lambda1,tau,B0,rho)
#################################################################################################################
# SIMULATIONS ----
## 1. DMPs (W correlated with X) -----

# Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = DMP_corr_parameters
# Type I Error Rate
for(N_i in 1:length(sample_size)){
  results_i = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
  for(sim in 1:Nsims){
    results_i = bind_rows(results_i,Sim.Function(data_sim=Scenario_params,
                                                 N=sample_size[N_i],
                                                 B1=B1_null,
                                                 depth=depth,
                                                 cpg_num=cpg_num))
  }
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}

# Statistical Power
for(N_i in 1:length(sample_size)){
  results_i = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
  for(sim in 1:Nsims){
    results_i = bind_rows(results_i,Sim.Function(data_sim=Scenario_params,
                                                 N=sample_size[N_i],
                                                 B1=B1_nonnull,
                                                 depth=depth,
                                                 cpg_num=cpg_num))
  }
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}

# Packing
Results_DMP_Corr = Sim_Results
rm(Sim_Results,Scenario_params)

## 2. Non-DMPs (W correlated with X) ----

# Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = NonDMP_corr_parameters

# Type I Error Rate
for(N_i in 1:length(sample_size)){
  results_i = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
  for(sim in 1:Nsims){
    results_i = bind_rows(results_i,Sim.Function(data_sim=Scenario_params,
                                                 N=sample_size[N_i],
                                                 B1=B1_null,
                                                 depth=depth,
                                                 cpg_num=cpg_num))
  }
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}

# Statistical Power
for(N_i in 1:length(sample_size)){
  results_i = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
  for(sim in 1:Nsims){
    results_i = bind_rows(results_i,Sim.Function(data_sim=Scenario_params,
                                                 N=sample_size[N_i],
                                                 B1=B1_nonnull,
                                                 depth=depth,
                                                 cpg_num=cpg_num))
  }
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}

# Packing
Results_NonDMP_Corr = Sim_Results
rm(Sim_Results,Scenario_params)

## **Not Currently in Use** 3. DMPs (W uncorrelated with X) -----

# Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = DMP_uncorr_parameters

# Type I Error Rate
for(N_i in 1:length(sample_size)){
  results_i = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
  for(sim in 1:Nsims){
    results_i = bind_rows(results_i,Sim.Function(data_sim=Scenario_params,
                                                 N=sample_size[N_i],
                                                 B1=B1_null,
                                                 depth=depth,
                                                 cpg_num=cpg_num))
  }
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}

# Statistical Power
for(N_i in 1:length(sample_size)){
  results_i = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
  for(sim in 1:Nsims){
    results_i = bind_rows(results_i,Sim.Function(data_sim=Scenario_params,
                                                 N=sample_size[N_i],
                                                 B1=B1_nonnull,
                                                 depth=depth,
                                                 cpg_num=cpg_num))
  }
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}

# Packing
Results_DMP_UnCorr = Sim_Results
rm(Sim_Results,Scenario_params)

## **Not Currently in Use** 4. Non-DMPs (W uncorrelated with X) -----

# Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
Scenario_params = NonDMP_uncorr_parameters

# Type I Error Rate
for(N_i in 1:length(sample_size)){
  results_i = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
  for(sim in 1:Nsims){
    results_i = bind_rows(results_i,Sim.Function(data_sim=Scenario_params,
                                                 N=sample_size[N_i],
                                                 B1=B1_null,
                                                 depth=depth,
                                                 cpg_num=cpg_num))
  }
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}

# Statistical Power
for(N_i in 1:length(sample_size)){
  results_i = data.frame("Effect"=NA,"Sample Size"=NA,"CLARE"=NA,"Proportion"=NA,"Unadjusted"=NA)
  for(sim in 1:Nsims){
    results_i = bind_rows(results_i,Sim.Function(data_sim=Scenario_params,
                                                 N=sample_size[N_i],
                                                 B1=B1_nonnull,
                                                 depth=depth,
                                                 cpg_num=cpg_num))
  }
  cat(paste0(sample_size[N_i],","))
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
}

# Packing
Results_NonDMP_UnCorr = Sim_Results
rm(Sim_Results,Scenario_params)

## 5. Packaging & Saving ----

# Notes for Individual Scenarios
Results_DMP_Corr$Scenario = "DMP Correlated"
# Results_DMP_UnCorr$Scenario = "DMP Uncorrelated"
Results_NonDMP_Corr$Scenario = "Non-DMP Correlated"
# Results_NonDMP_UnCorr$Scenario = "Non-DMP Uncorrelated"

# Complete Data Frame
Results_Final_Sim1 = bind_rows(Results_DMP_Corr
                               ,Results_NonDMP_Corr
                               #,Results_DMP_UnCorr
                               #,Results_NonDMP_UnCorr
                               
                               )%>%
  mutate(Effect=case_when(Effect==0~"Type I Error",TRUE~"Power"))%>%
  dplyr::filter(!is.na(Sample.Size))%>%
  pivot_longer(names_to="Method",cols=c("CLARE","Proportion","Unadjusted"))

#################################################################################################################
# QUICK FIGURES ----
## PLOT 1: ----
colors <- c("#15bf3a",
            "#4075d6",
            "#631307")
shapes <- c(15,1,4)
#
sample.size.breaks=unique(Results_Final_Sim1$Sample.Size)

df = Results_Final_Sim1 %>% 
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Unadjusted")))%>%
  mutate(Effect=factor(Effect,levels=c("Type I Error","Power")))%>%
  mutate(Scenario_2= case_when(Scenario=="DMP Correlated"~"DMP",
                               Scenario=="DMP Uncorrelated"~NA,
                               Scenario=="Non-DMP Correlated"~"Non-DMP",
                               Scenario=="Non-DMP Uncorrelated"~NA,TRUE~NA))%>%
  dplyr::filter(Method != "Unadjusted" & !is.na(Scenario_2))
#
title= "Simulation 1 Results"
plot1 <- ggplot(df,aes(x=Sample.Size,y=value,color=Method,shape=Method))+
  geom_line(size=1,alpha=0.5)+
  geom_point(size=2,alpha=0.9)+
  labs(title=title,
       x="Sample Size",
       y="Detection Rate")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario_2~Effect,scales="fixed")+
  scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1),minor_breaks=NULL)+
  geom_hline(yintercept=0.05,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plot1

## PLOT 2: Mean Power & FDR with SD Ribbons ----
colors <- c("#15bf3a",
            "#4075d6",
            "#631307")
shapes <- c(15,1,4)
#
sample.size.breaks=unique(Results_Final_Sim1$Sample.Size)

df = Results_Final_Sim1 %>% 
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Unadjusted")))%>%
  mutate(Effect=factor(Effect,levels=c("Type I Error","Power")))%>%
  group_by(Effect,Sample.Size,Scenario,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE),
            sd=sd(value,na.rm=TRUE))%>%
  mutate(Scenario_2= case_when(Scenario=="DMP Correlated"~"DMP",
                              Scenario=="Non-DMP Correlated"~"Non-DMP"
                               ,TRUE~NA))%>%
  dplyr::filter(Method != "Unadjusted" & !is.na(Scenario_2))%>%
  ungroup()
#
title= "Simulation 1 Results"
plot1 <- ggplot(df,aes(x=Sample.Size,y=mean,color=Method,shape=Method))+
  geom_line(size=1,alpha=0.5)+
  geom_point(size=2,alpha=0.9)+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,fill=Method),alpha=0.1,color=NA)+
  #stat_summary(fun=mean,geom="line",size=1,alpha=0.5)+
  #stat_summary(fun.data=mean_cl_normal,geom="ribbon",alpha=0.2)
  labs(title=title,
       x="Sample Size",
       y="Detection Rate")+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario_2~Effect,scales="fixed")+
  scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1),minor_breaks=NULL)+
  geom_hline(yintercept=0.05,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plot1

#################################################################################################################
# EXPORTING RESULTS ----
if(Outputs==TRUE){
  filename <- "Simulation 1.xlsx"
  filename <- paste0(Parent_Directory,"Results/",filename)
  data_list <- list(Data=Results_Final_Sim1)
  write_xlsx(data_list,filename)
}
