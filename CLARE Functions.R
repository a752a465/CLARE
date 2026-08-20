#################################################################################################################
### Title: CLARE - Functions
### Current Analyst: Alexander Alsup
### Last Updated: (08/20/2026) 
### Notes: Greedy compositional log-ratio selection, decoupled from any specific model
### CLARE - Generalized Version
### Core idea: the greedy search only needs a function that takes a candidate log-ratio vector and returns a scalar score to optimize.
#################################################################################################################
# Packages: Only One Required ----
packages <- c("easyCODA")
lapply(packages, require, character.only = TRUE)

##################################################
## Default Zero Handling Function ----
# The default zero handling function used by CLARE.
# This simply adds a small constant to zero proportions and rescales to compositions to 1
Scale_W <- function(W = W) {
  W_lr <- W
  if (any(W_lr == 0)) {
    W_lr[W_lr == 0] <- W_lr[W_lr == 0] + 1e-5
    W_lr <- W_lr / rowSums(W_lr)
  } else {
    W_lr <- as.matrix(W_lr)
  }
  return(as.matrix(W_lr))
}

##################################################
## Generic Greedy Log-Ratio Search ----

### Function Argument Notes -----

# fit_fun: function(lr, ...) -> a LIST with:

#     $score  (required) single numeric fit statistic for this candidate log-ratio

#     $extra  (optional) a named vector of anything else worth keeping per candidate (coefficients, p-values, concordance, whatever your model produces).
#             LR_Search does not interpret $extra at all -- it just stacks it into the output tables so you can inspect the selection process afterward.

# zero_handling_fun

# minimize: TRUE if lower score = better fit (MSE, AIC, deviance),
#           FALSE if higher score = better fit (R^2, concordance, logLik)

# var.total: max number of cell types allowed to enter the log-ratio (numerator+denominator combined)

# early_stop: if TRUE, stop growing the log-ratio once a step fails to improve on the previous best

# zero_handling_fun: Allows the user to supply their own zero-replacement strategy (pseudocount, detection limits, etc.)
# instead of CLARE's built-in default (Scale_W: floor + renormalize if any zeroes present).
# If zero_handling_fun is set to NULL and there are NO zeroes in the proportion data, CLARE will not rescale at all
# If zero_handling_fun is set to NULL and there ARE zeroes, CLARE will default to Scale_W()

# ...: passed through to fit_fun at every evaluation (e.g. outcome, covariates)

### Function Specification ----
LR_Search <- function(W, fit_fun,zero_handling_fun=Scale_W, var.total, minimize = TRUE, early_stop = TRUE,
                       verbose = FALSE, ...) {

  # Zero Handling
  if (!is.null(zero_handling_fun)) {
    W_lr <- zero_handling_fun(W)
  } else if (any(W == 0)) {
    if(verbose==TRUE){message("Zero proportions detected in W; no zero_handling_fun supplied -- defaulting to Scale_W() (floor + renormalize).")}
    W_lr <- Scale_W(W)
  } else {
    if(verbose==TRUE){message("No zero proportions detected in W; proceeding without zero handling.")}
    W_lr <- as.matrix(W)
  }
  
  # Conditional Check for the var.total input
  p <- ncol(W_lr)
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
    Best_Path   = Best_Path_df,      # one row per size: the winner at that size
    All_Evaluated = All_Evaluated_df # every candidate log-ratio evaluated, every size
  )
}

##################################################
