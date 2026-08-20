#################################################################################################################
### Title: Adjusted EWAS Using Compositional Data Transformations - Supplementary Simulation #5 ----
### Current Analyst: Alexander Alsup (a752a465@kumc.edu) ----
### Last Updated: (08/14/2026) ----
### Notes: ----
#################################################################################################################
# Source Functions & Parent Filepath ----
# Replace with your own directory where applicable
Parent_Directory = ""
source(paste0(Parent_Directory,"Supplementary Sim Prep.R"))
#################################################################################################################
## SUPPLEMENTARY SIMULATION 5.1 - Additional DA Cell Types - Effects Similar to Simulation 2 -----
## Sim Parameters ----
### Basic Parameters ----
sample_size=c(50,150,250,400,500) # Sample Size
Nsims= 100 # Outer Replicates - Number of Unique Cell Mixture Samples
cpg_num=1000 # Inner Replicates - Number of CpGs generated per Cell Mixture Sample

### Data Generation Parameters ----
total.cells=1000 # The total cells used for measuring DNAm
depth=200 # The probe depth

#### 4 Cell Type DA (Indirect Effect: -0.04248105) ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-8,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+3,"CD4_Mem"=0,"CD8_Naive"=+2,"CD8_Mem"=-2,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
Cell_4_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)
Theta_M=sum((lambda1/sum(lambda1))*B0)-sum((lambda0/sum(lambda0))*B0);Theta_M
#
rm(lambda0,lambda1,tau,B0,rho)

#### 5 Cell Type DA (Indirect Effect: -0.04152458) ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-8,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+3,"CD4_Mem"=0,"CD8_Naive"=+2,"CD8_Mem"=-2,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=+1)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
Cell_5_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)
Theta_M=sum((lambda1/sum(lambda1))*B0)-sum((lambda0/sum(lambda0))*B0);Theta_M
#
rm(lambda0,lambda1,tau,B0,rho)

#### 6 Cell Type DA (Indirect Effect: -0.04150184) ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-7,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+3.5,"CD4_Mem"=0,"CD8_Naive"=+2,"CD8_Mem"=-2,"B_Naive"=0,"B_Mem"=0,"Treg"=+0.5,"NK"=+1)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
Cell_6_parameters <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)
Theta_M=sum((lambda1/sum(lambda1))*B0)-sum((lambda0/sum(lambda0))*B0);Theta_M
#
rm(lambda0,lambda1,tau,B0,rho)

###############################
# SIMULATIONS ----

#### Set.Seed() ----
set.seed(999848)

## 4 Cell Types -----

### Prepare Results Dataframe & Scenario Parameters
Scenario_params = Cell_4_parameters
Sim_Results = data.frame()
Scenario_name = "4 DA Cells"

# Loop
cat(paste0("4 DA Cells: ","\n"))
for(N_i in 1:length(sample_size)){
  results_i = data.frame()

  # Statistical Power
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=sample_size[N_i],B1=B1_nonnull,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
  
  # Type I Error Rate
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=sample_size[N_i],B1=B1_null,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
  
  # Message
  cat(paste0(sample_size[N_i],"|"))
}
cat(paste0("\n"))
rm(Scenario_params)

## 5 Cell Types -----

### Prepare Results Dataframe & Scenario Parameters
Scenario_params = Cell_5_parameters
Scenario_name = "5 DA Cells"

# Loop
cat(paste0("5 DA Cells: ","\n"))
for(N_i in 1:length(sample_size)){
  results_i = data.frame()

  # Statistical Power
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=sample_size[N_i],B1=B1_nonnull,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
  
  # Type I Error Rate
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=sample_size[N_i],B1=B1_null,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
  
  # Message
  cat(paste0(sample_size[N_i],"|"))
}
cat(paste0("\n"))
rm(Scenario_params)

## 6 Cell Types -----

### Prepare Results Dataframe & Scenario Parameters
Scenario_params = Cell_6_parameters
Scenario_name = "6 DA Cells"

# Loop
cat(paste0("6 DA Cells: ","\n"))
for(N_i in 1:length(sample_size)){
  results_i = data.frame()
  
  # Statistical Power
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=sample_size[N_i],B1=B1_nonnull,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
  
  # Type I Error Rate
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=sample_size[N_i],B1=B1_null,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
  
  # Message
  cat(paste0(sample_size[N_i],"|"))
}
cat(paste0("\n"))
rm(Scenario_params)

## 3. Packaging & Saving ----

# Complete Data Frame
Results_Final_Sim1 = Sim_Results%>%
  dplyr::filter(!is.na(Sample.Size))%>%
  rename(CLARE_Detection=CLARE,Proportion_Detection=Proportion,Unadjusted_Detection=Unadjusted)%>%
  pivot_longer(cols = matches("_Power$|_TypeI$|_Bias$|_RMSE$|_Detection$|_Variance$"), names_to = c("Method","Metric"), names_sep = "_")


###############################
# SUMMARIES -----
### Double Checking RMSE Capture ----
max(abs(as.numeric(Sim_Results$CLARE_RMSE)^2 - (as.numeric(Sim_Results$CLARE_Bias)^2 + as.numeric(Sim_Results$CLARE_Variance))))
max(abs(as.numeric(Sim_Results$Unadjusted_RMSE)^2 - (as.numeric(Sim_Results$Unadjusted_Bias)^2 + as.numeric(Sim_Results$Unadjusted_Variance))))
max(abs(as.numeric(Sim_Results$Proportion_RMSE)^2 - (as.numeric(Sim_Results$Proportion_Bias)^2 + as.numeric(Sim_Results$Proportion_Variance))))

## Bias/RMSE/Variance ----

metric_means <- Results_Final_Sim1 %>%
  mutate(value=as.numeric(value))%>%
  group_by(B1_split,Scenario,Method,Metric)%>%
  summarise(Metric_Mean=mean(value,na.rm=TRUE))%>%
  arrange(Metric,Scenario,B1_split,Method)

metric_means2 <- Results_Final_Sim1 %>%
  mutate(value=as.numeric(value))%>%
  group_by(B1_split,Method,Metric)%>%
  summarise(Metric_Mean=mean(value,na.rm=TRUE))%>%
  arrange(Metric,B1_split,Method)

###############################
# QUICK FIGURES ----
## Overall prep ----
#
colors <- c("#15bf3a",
            "#4075d6",
            #"#674A05",
            "#631307")
shapes <- c(15,
            1,
            #12,
            4)
#
unique(Results_Final_Sim1$Effect)
sample.size.breaks=unique(Results_Final_Sim1$Sample.Size)
Simulation.title = "Supplementary Simulation 5: "
df <-  Results_Final_Sim1 %>%
  mutate(Metric=case_when(Effect==0 & Metric=="Detection"~"Type I Error",
                          Effect==-0.025 & Metric=="Detection"~"Power",
                          TRUE~Metric))%>%
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Unadjusted")))

## PLOT 1: Mean Power & Type I Error Rate with Linetype ----

df_graph = df %>%
  dplyr::filter(Metric %in% c("Type I Error","Power") 
                & Method != "Unadjusted" 
  )%>%
  group_by(Scenario,Metric,Sample.Size,,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE),
            sd=sd(value,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(Metric=factor(Metric,levels=c("Type I Error","Power")))
#
plotB1 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,fill=Method,linetype=Scenario))+
  geom_line(size=1,alpha=0.6)+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,fill=Method),alpha=0.1,color=NA)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Power & Type I Error Rate, by Proportion of CpGs where P(B1≠0)=p",
       x="Sample Size",
       y="Detection Rate",
       linetype="DA Cell Types")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Metric,scales="fixed")+
  scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1),minor_breaks=NULL)+
  geom_hline(yintercept=0.05,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB1



#################################################################################################################
## SUPPLEMENTARY SIMULATION 5.2 - Distributions of Type I Error Rates at N=500 -----
## Sim Parameters ----
### Basic Parameters ----
sample_size=c(50,150,250,400,500) # Sample Size
Nsims=200 # Outer Replicates - Number of Unique Cell Mixture Samples
cpg_num=1000 # Inner Replicates - Number of CpGs generated per Cell Mixture Sample
sims_it = seq((Nsims/20),Nsims,by=((Nsims/20)))

### Data Generation Parameters ----
total.cells=1000 # The total cells used for measuring DNAm
depth=200 # The probe depth

#### 4 Cell Type DA (Simulation #2 Configuration) (Indirect Effect: -0.04248105) ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-8,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+3,"CD4_Mem"=0,"CD8_Naive"=+2,"CD8_Mem"=-2,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
Cell_4_parameters_1 <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)
Relative_Shift = round(((lambda1-lambda0)/lambda0)*100,2);Relative_Shift
sum(lambda1);sum(lambda0)
Theta_M=sum((lambda1/sum(lambda1))*B0)-sum((lambda0/sum(lambda0))*B0);Theta_M
#
rm(lambda0,lambda1,tau,B0,rho)

#### 4 Cell Type DA (Effects Scaled Relative to Baseline Abundance) (Indirect Effect: -0.04247735) ----
lambda0=c("Neut"=50,"Mono"=5,"Eos"=3,"Baso"=1,"CD4_Naive"=12,"CD4_Mem"=7,"CD8_Naive"=3,"CD8_Mem"=4,"B_Naive"=4,"B_Mem"=3,"Treg"=2,"NK"=6)
tau = c("Neut"=-11.35,"Mono"=0,"Eos"=0,"Baso"=0,"CD4_Naive"=+2.724,"CD4_Mem"=0,"CD8_Naive"=+0.681,"CD8_Mem"=-0.908,"B_Naive"=0,"B_Mem"=0,"Treg"=0,"NK"=0)
B0=c("Neut"=0.8,"Mono"=0.518,"Eos"=0.236,"Baso"=0.647,"CD4_Naive"=0.2,"CD4_Mem"=0.383,"CD8_Naive"=0.199,"CD8_Mem"=0.811,"B_Naive"=0.236,"B_Mem"=0.405,"Treg"=0.429,"NK"=0.647)
rho=total.cells/sum(lambda0)
lambda1 = lambda0+tau
lambda0=lambda0*rho
lambda1=lambda1*rho
#
Cell_4_parameters_2 <- list("Lambda0"=lambda0,"Lambda1"=lambda1,"B0"=B0)
Relative_Shift = round(((lambda1-lambda0)/lambda0)*100,2);Relative_Shift
(lambda0/rho*0.227)[tau!=0]
sum(lambda1);sum(lambda0)
Theta_M=sum((lambda1/sum(lambda1))*B0)-sum((lambda0/sum(lambda0))*B0);Theta_M
#
rm(lambda0,lambda1,tau,B0,rho)

###############################
# SIMULATIONS ----

#### Set.Seed() ----
set.seed(999848)

## 4 Cell Type DA (Simulation #2 Configuration) -----

### Prepare Results Dataframe & Scenario Parameters
Scenario_params = Cell_4_parameters_1
Sim_Results = data.frame()
Scenario_name = "Simulation 2"

# Loop
cat(paste0("Simulation 2: ","\n"))
  results_i = data.frame()
  # Statistical Power
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=500,B1=B1_nonnull,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
    if(sim %in% sims_it){
      # Message
      cat(paste0(sim,"|"))
    }
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
  # Type I Error Rate
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=500,B1=B1_null,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
    # Message
    if(sim %in% sims_it){
      cat(paste0(sim,"|"))
    }
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
cat(paste0("\n"))
rm(Scenario_params)

## 4 Cell Type DA (Effects Scaled Relative to Baseline Abundance) -----

### Prepare Results Dataframe & Scenario Parameters
Scenario_params = Cell_4_parameters_2
Scenario_name = "Scaled"

# Loop
cat(paste0("Concentrated: ","\n"))
  results_i = data.frame()
  
  # Statistical Power
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=500,B1=B1_nonnull,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
    # Message
    if(sim %in% sims_it){
      cat(paste0(sim,"|"))
    }
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
  
  # Type I Error Rate
  for(sim in 1:Nsims){
    results_it = Sim.Function(data_sim=Scenario_params,N=500,B1=B1_null,depth=depth,cpg_num=cpg_num)
    results_it$Scenario = Scenario_name
    results_i = bind_rows(results_i,results_it)
    # Message
    if(sim %in% sims_it){
      cat(paste0(sim,"|"))
    }
  }
  results_i = results_i[!is.na(results_i$Effect),]
  Sim_Results = bind_rows(Sim_Results,results_i)
  
cat(paste0("\n"))
rm(Scenario_params)

## 3. Packaging & Saving ----

# Complete Data Frame
Results_Final_Sim1 = Sim_Results%>%
  dplyr::filter(!is.na(Sample.Size))%>%
  rename(CLARE_Detection=CLARE,Proportion_Detection=Proportion,Unadjusted_Detection=Unadjusted)%>%
  pivot_longer(cols = matches("_Power$|_TypeI$|_Bias$|_RMSE$|_Detection$|_Variance$"), names_to = c("Method","Metric"), names_sep = "_")

###############################
# SUMMARIES -----
### Double Checking RMSE Capture ----
max(abs(as.numeric(Sim_Results$CLARE_RMSE)^2 - (as.numeric(Sim_Results$CLARE_Bias)^2 + as.numeric(Sim_Results$CLARE_Variance))))
max(abs(as.numeric(Sim_Results$Unadjusted_RMSE)^2 - (as.numeric(Sim_Results$Unadjusted_Bias)^2 + as.numeric(Sim_Results$Unadjusted_Variance))))
max(abs(as.numeric(Sim_Results$Proportion_RMSE)^2 - (as.numeric(Sim_Results$Proportion_Bias)^2 + as.numeric(Sim_Results$Proportion_Variance))))

###############################
# QUICK FIGURES ----
## Overall prep ----
#
colors <- c("#15bf3a",
            "#4075d6",
            #"#674A05",
            "#631307")
shapes <- c(15,
            1,
            #12,
            4)
#
sample.size.breaks=unique(Results_Final_Sim1$Sample.Size)
Simulation.title = "Supplementary Simulation 5.2: "
#
df <-  Results_Final_Sim1 %>%
  mutate(Metric=case_when(Effect==0 & Metric=="Detection"~"Type I Error",
                          Effect==-0.025 & Metric=="Detection"~"Power",
                          TRUE~Metric))%>%
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Unadjusted")))

## PLOT 1: Histogram with Both ----

df_graph = df %>%
  dplyr::filter(Metric %in% c("Type I Error") 
                & Method != "Unadjusted" 
  )%>%
  mutate(Method = factor(Method, levels = c("CLARE","Proportion")))
# Plot
plotB1 <- ggplot(df_graph, aes(x = value, fill = Method))+
  geom_histogram(binwidth = 0.005, alpha = 0.7, color = NA, position = "identity") +
  geom_vline(xintercept = 0.05, color = "red", alpha = 0.5, linetype = "dashed") +
  facet_grid(Scenario ~ Method, scales = "free_y") +
  scale_fill_manual(values = c("CLARE" = "#15bf3a", "Proportion" = "#4075d6")) +
  labs(title = "Supplementary Simulation 5.2:",
       subtitle = "Distribution of Per-Replicate Type I Error Rates, by Configuration",
       x = "Type I Error Rate (per replicate)",
       y = "Count") +
  theme(legend.position = "none",   # redundant with facet columns already labeling Method
        strip.text = element_text(size = 10))
plotB1

## PLOT 2: Histogram with Simulation 4 Type I error rate -----
df_graph = df %>%
  dplyr::filter(Metric %in% c("Type I Error") 
                & Method != "Unadjusted"
                & Scenario == "Simulation 2"
  )%>%
  mutate(Method = factor(Method, levels = c("CLARE","Proportion")))
# Plot
plotB1 <- ggplot(df_graph, aes(x = value, fill = Method))+
  geom_histogram(binwidth = 0.01, alpha = 1, color="black", position = "identity") +
  geom_vline(xintercept = 0.05, color = "red", alpha = 0.5, linetype = "dashed") +
  facet_wrap(~Method, ncol=1, scales = "fixed") +
  scale_x_continuous(limits=c(0,1))+
  scale_fill_manual(values = c("CLARE" = "#15bf3a", "Proportion" = "#4075d6")) +
  labs(title = "Supplementary Simulation 5.2:",
       subtitle = "Distribution of Per-Replicate Type I Error Rates, by Configuration",
       x = "Type I Error Rate (per replicate)",
       y = "Count") +
  theme(legend.position = "none",   # redundant with facet columns already labeling Method
        strip.text = element_text(size = 10))
plotB1

###############################
# EXPLORATION -----

## Find a problematic seed
seed_it <- ceiling(runif(1,1,10000))
set.seed(seed_it)
Scenario_params = Cell_4_parameters_3
Scenario_name = "Dispersed"
Sim.Function(data_sim=Scenario_params,N=500,B1=B1_null,depth=depth,cpg_num=cpg_num)

## Generate Data from the seed and check manually
N=500
Scenario_params = Cell_4_parameters_3
set.seed(6225)
Data <- generate.data(data_sim=Scenario_params,N=N,B1=B1_null,depth=depth,cpg_num=cpg_num)
X_var = as.matrix(Data$X) ; W_var = as.matrix(Data$W) ; Y_var = as.matrix(Data$Y)
# LR Selection 
LR_Selection_Result = LR_Selection.v3(X=Data$X,W=W_var,var.total=(ncol(W_var)-1))
lr_covariate = LR_Selection_Result[[1]]
colnames(W_var)[c(7,8)]
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

