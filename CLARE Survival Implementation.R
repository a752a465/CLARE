#################################################################################################################
### Title: CLARE - Implementation for Survival Biomarker Discovery
### Current Analyst: Alexander Alsup
### Last Updated: (08/04/2026) 
### Notes: 
#################################################################################################################
# Packages and Sources Scripts -----
library(survival)
library(easyCODA)

# Function Source: Re-write the file path based on your download location
file_path = "S:/Biostats/BIO-STAT/Koestler Devin/GRAs/Alex Alsup/Project 3/GitHub Repo/update 08_12/CLARE demos/"
source(paste0(file_path,"CLARE_Functions.R"))

##################################################
# ReadMe -----

# CLARE is a simple greedy stepwise search that builds Log-Ratios (LRs) as potential predictors and selects the LR which MOST improves the model fit.
# CLARE accepts fit statistics for assessment, and model specification through a fit_fun() argument. 
# The contribution that CLARE uniquely makes is building the Log-Ratios to be assessed and selecting the best LR. 

# Below, we explore an example of how CLARE could be used to identify an SLR (Summated Log-Ratio) that is pre-specified to be associated with survival. 

##################################################
## Data Simulation Function ----

## N patients, K immune cell types 
## (Poisson counts -> compositional proportions W)
simulate_survival_data <- function(N, cell_lambda, ratio_numerator, ratio_denominator, beta_ratio,
                                   covariates = NULL, beta_cov = NULL,
                                   baseline_hazard = 0.02, admin_censor_time = 10,
                                   extra_censor_rate = NULL, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if (!all(ratio_numerator %in% names(cell_lambda)) || !all(ratio_denominator %in% names(cell_lambda))) {
    stop("ratio_numerator and ratio_denominator must reference names in cell_lambda")
  }
  if (length(intersect(ratio_numerator, ratio_denominator)) > 0) {
    stop("ratio_numerator and ratio_denominator cannot share cell types")
  }
  
  K <- length(cell_lambda)
  counts <- sapply(cell_lambda, function(lam) rpois(N, lam))
  colnames(counts) <- names(cell_lambda)
  W <- as.data.frame(counts / rowSums(counts))
  
  W_lr <- Scale_W(W)
  true_ratio <- rowMeans(log(W_lr[, ratio_numerator, drop = FALSE])) -
    rowMeans(log(W_lr[, ratio_denominator, drop = FALSE]))
  true_composition_effect <- beta_ratio * true_ratio
  
  lp <- true_composition_effect
  if (!is.null(covariates)) {
    if (is.null(beta_cov) || length(beta_cov) != ncol(covariates)) {
      stop("beta_cov must have one entry per column of covariates")
    }
    lp <- lp + as.matrix(covariates) %*% beta_cov
  }
  
  # Exponential event times under Cox PH: T = -log(U) / (h0 * exp(lp))
  U <- runif(N)
  time_event <- -log(U) / (baseline_hazard * exp(lp))
  
  censor_time <- rep(admin_censor_time, N)
  if (!is.null(extra_censor_rate)) {
    censor_time <- pmin(censor_time, rexp(N, rate = extra_censor_rate))
  }
  
  time <- pmin(time_event, censor_time)
  status <- as.numeric(time_event <= censor_time)
  
  list(W = W, time = time, status = status,
       true_ratio = true_ratio, true_composition_effect = true_composition_effect,
       active_cells = c(ratio_numerator, ratio_denominator),
       covariates = covariates, pct_censored = mean(status == 0))
}

##################################################
## Demo ----

# Sample Size of the Example
N=600

### Cell Composition Effects -----

## 12 immune cell types measured for each patient. 
## Three of them (CD4 Naive, NK, Monocytes) truly drive survival.
## The other 9 cell types are pure noise with respect to the outcome.
cell_lambda <- c(Neut = 50, Mono = 5, Eos = 3, Baso = 1,
                 CD4_Naive = 12, CD4_Mem = 7, CD8_Naive = 3, CD8_Mem = 4,
                 B_Naive = 4, B_Mem = 3, Treg = 2, NK = 6)


## True biomarker: (Mono) / (CD8_Naive + CD8_Mem) -- myeloid burden relative to CD8 Abundance. 
## so beta_ratio > 0 means a HIGHER ratio -> HIGHER hazard -> WORSE (decreased) survival.
ratio_numerator   <- c("Mono")
ratio_denominator <- c("CD8_Naive", "CD8_Mem")
beta_ratio <- 0.6

### Covariate Effects -----

# We also include Age and Sex, with their own effects independent of composition
covariates <- data.frame(
  Age = rnorm(N, mean = 60, sd = 10),
  Sex = rbinom(N, 1, 0.5)   # 1 = male
)
beta_cov = c(0.02, 0.3)

### Data Simulation ----
sim <- simulate_survival_data(
  N = N,
  cell_lambda = cell_lambda,
  ratio_numerator = ratio_numerator,
  ratio_denominator = ratio_denominator,
  beta_ratio = beta_ratio,
  covariates = covariates,
  beta_cov = c(0.02, 0.3),   # Age, Sex effects on log-hazard: older and male -> higher hazard
  baseline_hazard = 0.02,
  admin_censor_time = 10,
  seed = 8976
)

cat("True biomarker: (", paste(ratio_numerator, collapse = "+"), ") / (",
    paste(ratio_denominator, collapse = "+"), ")\n", sep = "")
cat("Percent censored:", round(100 * sim$pct_censored, 1), "%\n")

##################################################
## Step 1: Look at each immune cell type on its own ----
## The true signal here is a RATIO among CD8, NK, and Mono, not any single cell type's level.

univariate_results <- lapply(names(cell_lambda), function(cell) {
  df <- data.frame(time = sim$time, status = sim$status, prop = sim$W[[cell]])
  mod <- coxph(Surv(time, status) ~ prop, data = df)
  s <- summary(mod)
  data.frame(cell = cell,
             coef = s$coefficients[1, "coef"],
             p = s$coefficients[1, "Pr(>|z|)"])
})
cat("\n--- Step 1: univariate association of each cell type's proportion with survival ---\n")
print(do.call(rbind, univariate_results))

##################################################
## Step 2: Covariates on their own ----
cat("\n--- Step 2: covariate-only model (Age, Sex; no immune composition) ---\n")
cov_mod <- coxph(Surv(time, status) ~ Age + Sex, data = cbind(time = sim$time, status = sim$status, sim$covariates))
print(summary(cov_mod)$coefficients)

# Our estimates for Age (0.0178) and Sex (0.2017) 

##################################################
## Step 3: Run CLARE to discover the composition biomarker, adjusting for Age/Sex ----

### 3.1 Define the Testing Function ----

# We define a Testing function to accept the covariates and generated Log-Ratios for Testing
# The function can then output some additional results for us to look at, and is flexible. 

fit_function_cox_concordance <- function(lr, time, status, covariates = NULL, ...) {
  # Collecting the log-ratio and other data into a single data frame
  df <- data.frame(lr = lr, time = time, status = status)
  rhs <- "lr"
  # Conditionally adding the covariates if present
  if (!is.null(covariates)) {
    df <- cbind(df, covariates)
    rhs <- paste(c("lr", colnames(covariates)), collapse = " + ")
  }
  
  # Fitting the Model for all included covariates
  form <- as.formula(paste("survival::Surv(time, status) ~", rhs))
  mod <- survival::coxph(form, data = df)
  s <- summary(mod)
  
  ## Results
  list(
    # Our "score" output is what the greedy selection uses
    score = unname(s$concordance["C"]),
    
    # The "extra" output is additional information we can look at
    extra = c(
      lr_coef       = unname(coef(mod)["lr"]),
      lr_hr         = unname(exp(coef(mod)["lr"])),
      lr_p          = unname(s$coefficients["lr", "Pr(>|z|)"]),
      AIC = AIC(mod),
      concordance   = unname(s$concordance["C"])
    )
  )
}


### 3.2 Run the Test ----

result <- CLARE(
  W = sim$W,
  fit_fun_select = fit_function_cox_concordance,
  zero_handling_fun=Scale_W,
  var.total = 5,
  minimize = FALSE, # Whether to find the minimal value for the selection. Since we are using Concordance, we want the highest value 
  verbose = TRUE,
  time = sim$time,
  status = sim$status,
  covariates = sim$covariates
)

cat("\n--- Step 3: CLARE search ---\n")
cat("\nSelected numerator cells: ", paste(result$selection$Numerator, collapse = " + "), "\n")
cat("Selected denominator cells:", paste(result$selection$Denominator, collapse = " + "), "\n")
cat("True active Ratio was: ", paste0(ratio_numerator,"/(",paste(ratio_denominator,collapse="+"),")","\n"))

# For this seed, Mono/(CD8_Naive+CD8_Mem) is selected as the best LR predictor. 
# Feel free to try other seeds, other fit function specifications, or whatever you like. 

cat("\n--- Best log-ratio at each size (Best_Path) ---\n")
print(result$selection$Best_Path)

### 3.3 Fit the Full Survival Model ----
lr_covar <- result$selection$LR
final_mod <- coxph(Surv(time, status) ~ Age + Sex + lr_covar, data = cbind(time = sim$time, status = sim$status, sim$covariates, lr_covar))
summary(final_mod)

# We detect effects in Age, Sex, and our selected SLR: ln[Mono/(CD8_Naive+CD8_Mem)]
# Compared to the covariates-only model, our estimate for Sex is now less biased (0.2897 vs 0.2017 ; True=0.3)
# While our estimate for Age has become slightly more biased (0.0171 vs 0.0178)
# Our estimate for the SLR is slightly inflated compared to ground truth (0.735 vs 0.6)




