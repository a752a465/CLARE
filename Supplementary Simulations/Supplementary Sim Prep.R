#################################################################################################################
### Title: Adjusted EWAS Using Compositional Data Transformations - Supplementary Simulation Prep ----
### Current Analyst: Alexander Alsup (a752a465@kumc.edu) ----
### Last Updated: (08/12/2026) ----
### Notes: ----
#################################################################################################################
# PACKAGES -----
packages <- c("ggplot2","here","MASS","MCMCpack","limma","betareg","tidyverse","writexl","readxl","easyCODA","limma","dplyr")
# Filter out packages that are already installed
missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# Install missing packages
if(length(missing_packages)) {install.packages(missing_packages)}
# Load all packages
lapply(packages, require, character.only = TRUE)
# If limma needs to be installed 
# BiocManager::install("limma", force = TRUE)

# Clear Environment ----
rm(list=ls())
# Source Functions & Parent Filepath ----
Parent_Directory = "S:/Biostats/BIO-STAT/Koestler Devin/GRAs/Alex Alsup/Project 3/GitHub Repo/update 08_12/"
source(paste0(Parent_Directory,"Functions_v2.R"))

#################################################################################################################
# GLOBAL PARAMETERS -----
## Outputs ----
Outputs=TRUE
### Basic Parameters ----
sample_size=c(50,100,150,200,250,300,400,500) # Sample Size
Nsims= 50 # Outer Replicates - Number of Unique Cell Mixture Samples
cpg_num=500 # Inner Replicates - Number of CpGs generated per Cell Mixture Sample

### Data Generation Parameters ----
total.cells=1000 # The total cells used for measuring DNAm
depth=200 # The probe depth
B1_null = 0 # Null Main Effect
B1_nonnull = -0.025 # Non-Null Main Effect
