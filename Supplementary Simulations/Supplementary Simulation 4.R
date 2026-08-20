#################################################################################################################
### Title: Adjusted EWAS Using Compositional Data Transformations - Supplementary Simulation #4 ----
### Current Analyst: Alexander Alsup (a752a465@kumc.edu) ----
### Last Updated: (08/14/2026) ----
### Notes: ----
#################################################################################################################
# Source Functions & Parent Filepath ----
# Replace with your own directory where applicable
Parent_Directory = ""
source(paste0(Parent_Directory,"Supplementary Sim Prep.R"))
#################################################################################################################
## SUPPLEMENTARY SIMULATION 4 - Comparison to Reference-Free Methods -----
## Functions ----
# ReFACTor Function
refactor_matrix <- function(O, k, covars = NULL, t = 500, numcomp = NULL, stdth = 0.02,
                            write_output = FALSE, out = "refactor") {
  
  O <- as.matrix(O)
  storage.mode(O) <- "numeric"
  cpgnames <- rownames(O)
  if (is.null(cpgnames)) cpgnames <- paste0("cpg", seq_len(nrow(O)))
  
  # Exclude low-variance sites (unchanged from original)
  sds <- apply(t(O), 2, sd)
  m_before <- length(sds)
  include <- which(sds >= stdth)
  O <- O[include, , drop = FALSE]
  cpgnames <- cpgnames[include]
  #message(sprintf("%d sites excluded due to low variance (< %.3g)", m_before - length(include), stdth))
  
  if (is.null(numcomp) || is.na(numcomp)) numcomp <- k
  
  # Adjust for covariates (unchanged logic, just no file read)
  if (!is.null(covars)) {
    covars <- as.matrix(covars)
    storage.mode(covars) <- "numeric"
    O_adj <- O
    for (site in 1:nrow(O)) {
      model <- lm(O[site, ] ~ covars)
      O_adj[site, ] <- residuals(model)
    }
    O <- O_adj
  }
  
  # Standard PCA (unchanged)
  pcs <- prcomp(scale(t(O)))
  coeff <- pcs$rotation
  score <- pcs$x
  
  # Low-rank approximation + site ranking (unchanged)
  x <- score[, 1:k] %*% t(coeff[, 1:k])
  An <- scale(t(O), center = TRUE, scale = FALSE)
  Bn <- scale(x, center = TRUE, scale = FALSE)
  An <- t(t(An) * (1 / sqrt(apply(An^2, 2, sum))))
  Bn <- t(t(Bn) * (1 / sqrt(apply(Bn^2, 2, sum))))
  
  distances <- apply((An - Bn)^2, 2, sum)^0.5
  dsort <- sort(distances, index.return = TRUE)
  ranked_list <- dsort$ix
  
  # ReFACTor components from the top-t ranked sites (unchanged)
  sites <- ranked_list[1:t]
  pcs2 <- prcomp(scale(t(O[sites, , drop = FALSE])))
  first_score <- score[, 1:k]
  score2 <- pcs2$x
  
  if (write_output) {
    write(t(cpgnames[ranked_list]), file = paste0(out, ".out.rankedlist.txt"), ncol = 1)
    write(t(score2[, 1:numcomp]), file = paste0(out, ".out.components.txt"), ncol = numcomp)
  }
  
  list(
    refactor_components = score2[, 1:numcomp, drop = FALSE],
    ranked_list = ranked_list,
    standard_pca = first_score
  )
}
# Data Generation Function 
generate.data.suppsim4 <- function(data_sim=data_sim,N=N,depth=depth,B1_effects=list("B1_null"=0,"B1_nonnull"=-0.025),cpg_num=cpg_num,B1_prop_sim=NULL){
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
  ## B1=True -----
  B1trueindex=1:(cpg_num*B1_prop_sim)
  B1_data=B1_effects$B1_nonnull
  Cell.Methylation = w_mat %*% B0
  Total.Methylation <- Cell.Methylation+(B1_data*x_vec)
  alpha = sapply(1:length(Total.Methylation), function(x) Total.Methylation[x]*depth)
  beta = sapply(1:length(Total.Methylation), function(x) (1-Total.Methylation[x])*depth)
  # CpG Data
  num_B1T=(cpg_num*B1_prop_sim)
  Y = matrix(NA,nrow=num_B1T,ncol=N)
  for(k in 1:num_B1T){
    Y_i=NULL
    for(i in 1:N){
      Y_i[i]<-rbeta(1,alpha[i],beta[i])
    }
    Y[k,] = Y_i
  }
  ## B1=Null -----
  B1falseindex=(length(B1trueindex)+1):(cpg_num)
  B1_data=B1_effects$B1_null
  Cell.Methylation = w_mat %*% B0
  Total.Methylation <- Cell.Methylation+(B1_data*x_vec)
  alpha = sapply(1:length(Total.Methylation), function(x) Total.Methylation[x]*depth)
  beta = sapply(1:length(Total.Methylation), function(x) (1-Total.Methylation[x])*depth)
  # CpG Data
  num_B1F=cpg_num-num_B1T
  Y_2 = matrix(NA,nrow=num_B1F,ncol=N)
  for(k in 1:num_B1F){
    Y_i=NULL
    for(i in 1:N){
      Y_i[i]<-rbeta(1,alpha[i],beta[i])
    }
    Y_2[k,] = Y_i
  }
  Y=rbind(Y,Y_2)
  #
  CpG_Variances = apply(Y,1,var,na.rm=TRUE)
  # Return
  result <- list("X"=X,"W"=W,"Y"=Y,"Effect_indices"=list("B1true"=B1trueindex,"B1false"=B1falseindex))
  return(result) 
}
# Simulation Function
Sim.Function.suppsim4 <- function(data_sim=data_sim,N=N,depth=depth,cpg_num=cpg_num,B1_prop_sim=NULL){
  B1_effects = list("B1_null"=0,"B1_nonnull"=-0.025)
  # Data Generation
  num_B1T=(cpg_num*B1_prop_sim)
  Data <- generate.data.suppsim4(data_sim=data_sim,N=N,depth=depth,B1_effects=B1_effects,cpg_num=cpg_num,B1_prop_sim=B1_prop_sim)
  X_var = as.matrix(Data$X) ; W_var = as.matrix(Data$W) ; Y_var = as.matrix(Data$Y) ; Effect_indices=Data$Effect_indices
  B1_vector=c(rep(B1_effects$B1_nonnull,times=length(Effect_indices$B1true)),
              rep(B1_effects$B1_null,times=length(Effect_indices$B1false)))
  # ReFACTor
  Y_var_uncorr = Y_var[Effect_indices$B1false,]
  t_sim = cpg_num/10 ; t_sim = 30
  refactor_covars = refactor_matrix(O=Y_var, k=12, covars = NULL, t = t_sim, 
                                    numcomp = NULL, stdth = 0.02,write_output = FALSE, out = "refactor")$refactor_components
  limma_Covars = as.matrix(cbind(rep(1,times=nrow(refactor_covars)),X_var,refactor_covars))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,2]
  Est_orig  <- fit$coefficients[,2]
  True_CpGs = Test_orig[Effect_indices$B1true] ; False_CpGs = Test_orig[Effect_indices$B1false]
  Refactor = length(True_CpGs[which(True_CpGs<=0.05)])/length(True_CpGs) ; Refactor_TypeI = length(False_CpGs[which(False_CpGs<=0.05)])/length(False_CpGs)
  Refactor_Bias = mean(Est_orig - B1_vector, na.rm=TRUE) ; Refactor_RMSE = sqrt(mean((Est_orig - B1_vector)^2, na.rm=TRUE)) ; Refactor_Variance = mean((Est_orig - mean(Est_orig, na.rm=TRUE))^2, na.rm=TRUE) 
  # CLARE
  LR_Selection_Result = LR_Selection.v3(X=Data$X,W=W_var,var.total=(ncol(W_var)-1))
  lr_covariate = LR_Selection_Result[[1]]
  limma_Covars = as.matrix(cbind(rep(1,times=length(lr_covariate)),X_var,lr_covariate))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,2]
  Est_orig  <- fit$coefficients[,2]
  True_CpGs = Test_orig[Effect_indices$B1true] ; False_CpGs = Test_orig[Effect_indices$B1false]
  CLARE = length(True_CpGs[which(True_CpGs<=0.05)])/length(True_CpGs) ; CLARE_TypeI = length(False_CpGs[which(False_CpGs<=0.05)])/length(False_CpGs)
  CLARE_Bias = mean(Est_orig - B1_vector, na.rm=TRUE) ; CLARE_RMSE = sqrt(mean((Est_orig - B1_vector)^2, na.rm=TRUE)) ; CLARE_Variance = mean((Est_orig - mean(Est_orig, na.rm=TRUE))^2, na.rm=TRUE) 
  # Proportion
  limma_Covars = as.matrix(cbind(X_var,W_var))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,1]
  Est_orig  <- fit$coefficients[,1]
  True_CpGs = Test_orig[Effect_indices$B1true] ; False_CpGs = Test_orig[Effect_indices$B1false]
  Proportion = length(True_CpGs[which(True_CpGs<=0.05)])/length(True_CpGs) ; Proportion_TypeI = length(False_CpGs[which(False_CpGs<=0.05)])/length(False_CpGs)
  Proportion_Bias = mean(Est_orig - B1_vector, na.rm=TRUE) ; Proportion_RMSE = sqrt(mean((Est_orig - B1_vector)^2, na.rm=TRUE)) ; Proportion_Variance = mean((Est_orig - mean(Est_orig, na.rm=TRUE))^2, na.rm=TRUE)
  # Unadjusted
  limma_Covars = as.matrix(cbind(rep(1,times=length(lr_covariate)),X_var))
  fit <- lmFit(Y_var,limma_Covars)
  fit <- eBayes(fit)
  Test_orig <- fit$p.value[,2]
  Est_orig  <- fit$coefficients[,2]
  True_CpGs = Test_orig[Effect_indices$B1true] ; False_CpGs = Test_orig[Effect_indices$B1false]
  Unadjusted = length(True_CpGs[which(True_CpGs<=0.05)])/length(True_CpGs) ; Unadjusted_TypeI = length(False_CpGs[which(False_CpGs<=0.05)])/length(False_CpGs)
  Unadjusted_Bias = mean(Est_orig - B1_vector, na.rm=TRUE) ; Unadjusted_RMSE = sqrt(mean((Est_orig - B1_vector)^2, na.rm=TRUE)) ; Unadjusted_Variance = mean((Est_orig - mean(Est_orig, na.rm=TRUE))^2, na.rm=TRUE)
  # Return
  results = data.frame("Sample Size"=N,"B1_split"=B1_prop_sim,
                       "CLARE_Power"=CLARE,"CLARE_TypeI"=CLARE_TypeI,"CLARE_Bias"=CLARE_Bias,"CLARE_RMSE"=CLARE_RMSE,"CLARE_Variance"=CLARE_Variance,
                       "Refactor_Power"=Refactor,"Refactor_TypeI"=Refactor_TypeI,"Refactor_Bias"=Refactor_Bias,"Refactor_RMSE"=Refactor_RMSE,"Refactor_Variance"=Refactor_Variance,
                       "Proportion_Power"=Proportion,"Proportion_TypeI"=Proportion_TypeI,"Proportion_Bias"=Proportion_Bias,"Proportion_RMSE"=Proportion_RMSE,"Proportion_Variance"=Proportion_Variance,
                       "Unadjusted_Power"=Unadjusted,"Unadjusted_TypeI"=Unadjusted_TypeI,"Unadjusted_Bias"=Unadjusted_Bias,"Unadjusted_RMSE"=Unadjusted_RMSE,"Unadjusted_Variance"=Unadjusted_Variance
  )
  return(results)
}

## Sim Parameters ----
### Basic Parameters ----
sample_size=c(50,150,250,400,500) # Sample Size
Nsims= 50 # Outer Replicates - Number of Unique Cell Mixture Samples
cpg_num=1000 # Inner Replicates - Number of CpGs generated per Cell Mixture Sample

### Data Generation Parameters ----
total.cells=1000 # The total cells used for measuring DNAm
depth=200 # The probe depth
# Vector of the proportion of CpGs with a true effect
B1_prop_vec = c(0.1,0.2,0.3)

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


#################################################################################################################
# SIMULATIONS ----

#### Set.Seed() ----
set.seed(999848)
## 1. DMPs (W correlated with X) -----

### Prepare Results Dataframe & Scenario Parameters
Sim_Results = data.frame()
Scenario_params = DMP_corr_parameters
Scenario_name = "cDMP"

### All Effect Proportions = c(0.1,0.2,0.3,0.4)
for(split in 1:length(B1_prop_vec)){
  B1_prop_sim=B1_prop_vec[split]
  cat(paste0("Proportion X-Associated CpGs: ",B1_prop_sim))
  cat(paste0("\n"))
  for(N_i in 1:length(sample_size)){
    results_i = data.frame()
    for(sim in 1:Nsims){
      results_it = Sim.Function.suppsim4(data_sim=Scenario_params,N=sample_size[N_i],depth=depth,cpg_num=cpg_num,B1_prop_sim=B1_prop_sim)
      results_it$Scenario = Scenario_name
      results_i = bind_rows(results_i,results_it)
    }
    cat(paste0(sample_size[N_i],"|"))
    Sim_Results = bind_rows(Sim_Results,results_i)
    
  }
  cat(paste0("\n"))
}
rm(Scenario_params)

## 2. Non-DMPs (W correlated with X) ----

### All Effect Proportions = c(0.1,0.2,0.3,0.4)
Scenario_params = NonDMP_corr_parameters
Scenario_name = "ncDMP"

for(split in 1:length(B1_prop_vec)){
  B1_prop_sim=B1_prop_vec[split]
  cat(paste0("Proportion X-Associated CpGs: ",B1_prop_sim))
  cat(paste0("\n"))
  for(N_i in 1:length(sample_size)){
    results_i = data.frame()
    for(sim in 1:Nsims){
      results_it = Sim.Function.suppsim4(data_sim=Scenario_params,N=sample_size[N_i],depth=depth,cpg_num=cpg_num,B1_prop_sim=B1_prop_sim)
      results_it$Scenario = Scenario_name
      results_i = bind_rows(results_i,results_it)
    }
    cat(paste0(sample_size[N_i],"|"))
    Sim_Results = bind_rows(Sim_Results,results_i)
  }
  cat(paste0("\n"))
}


## 3. Packaging & Saving ----

# Wide Results
Sim_Results = Sim_Results %>% dplyr::filter(!is.na(CLARE_TypeI))

# Complete Data Frame
Results_Final_Sim1 = Sim_Results%>%
  dplyr::filter(!is.na(Sample.Size))%>%
  pivot_longer(cols = matches("_Power$|_TypeI$|_Bias$|_RMSE$|_Detection$|_Variance$"), names_to = c("Method","Metric"), names_sep = "_")


#################################################################################################################
# SUMMARIES -----
### Double Checking RMSE Capture ----
max(abs(as.numeric(Sim_Results$CLARE_RMSE)^2 - (as.numeric(Sim_Results$CLARE_Bias)^2 + as.numeric(Sim_Results$CLARE_Variance))))
max(abs(as.numeric(Sim_Results$Unadjusted_RMSE)^2 - (as.numeric(Sim_Results$Unadjusted_Bias)^2 + as.numeric(Sim_Results$Unadjusted_Variance))))
max(abs(as.numeric(Sim_Results$Proportion_RMSE)^2 - (as.numeric(Sim_Results$Proportion_Bias)^2 + as.numeric(Sim_Results$Proportion_Variance))))
max(abs(as.numeric(Sim_Results$Refactor_RMSE)^2 - (as.numeric(Sim_Results$Refactor_Bias)^2 + as.numeric(Sim_Results$Refactor_Variance))))

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

#################################################################################################################
# QUICK FIGURES ----
## Overall prep ----
#
colors <- c("#15bf3a",
            "#4075d6",
            "#674A05",
            "#631307")
shapes <- c(15,1,12,4)
sample.size.breaks=unique(Results_Final_Sim1$Sample.Size)
Simulation.title = "Supplementary Simulation 4: "
df <-  Results_Final_Sim1 %>%
  mutate(Method=factor(Method,levels=c("CLARE","Proportion","Refactor")))%>%
  mutate(Metric=case_when(Metric=="TypeI"~"Type I Error",TRUE~Metric))%>%
  mutate(`Proportion B1≠0 Label`=paste0("p=",B1_split))%>%
  mutate(`Proportion B1≠0 Label`=factor(`Proportion B1≠0 Label`,levels=paste0("p=",sort(unique(B1_split),decreasing=FALSE))))

## PLOT 1: Mean Power & Type I Error Rate with Linetype ----
df_graph = df %>%
  dplyr::filter(Metric %in% c("Type I Error","Power") 
                & Method != "Unadjusted" 
  )%>%
  group_by(Scenario,Metric,Sample.Size,`Proportion B1≠0 Label`,Method)%>%
  summarise(mean=mean(value,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(Metric=factor(Metric,levels=c("Type I Error","Power")))
#
plotB1 <- ggplot(df_graph,aes(x=Sample.Size,y=mean,color=Method,shape=Method,linetype=`Proportion B1≠0 Label`))+
  geom_line(size=1,alpha=0.6)+
  geom_point(size=2,alpha=0.9)+
  labs(title=Simulation.title,
       subtitle="Power & Type I Error Rate, by Proportion of CpGs where P(B1≠0)=p",
       x="Sample Size",
       y="Detection Rate",
       linetype="P(B1≠0)")+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  facet_grid(Scenario~Metric,scales="fixed")+
  scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1),minor_breaks=NULL)+
  geom_hline(yintercept=0.05,color="red",alpha=0.3,linetype="dashed")+
  scale_x_continuous(breaks=sample.size.breaks,minor_breaks=NULL)+
  theme(legend.position="top",strip.text=element_text(size=12),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text=element_text(size=9))
plotB1


