
rm(list = ls())
library(methods)
library(dplyr)
library(ggplot2)
source("R/generics.R")

## ---- 1) Source core + designs (order matters)
source("R/randPar.R")
source("R/randSeq.R")

source("R/mwudPar.R")
source("R/mwudSeq.R")

source("R/dludPar.R")
source("R/dludSeq.R")

source("R/bsdPar.R")
source("R/bsdSeq.R")
source("R/Evaluation.R")

source("R/Crd.R")
source("R/MaxEnt.R")
source("R/plots.R")



## ---- 2) Fixed parameters
N      <- 200
Nsim   <- c(10,100,500)
ratios <- list(c(1,1,1,1), c(2,1,1,2), c(4,3,2,1), c(37,21,21,21))
mwud_alphas <- c(2,6,8)
dlud_as     <- c(2,6,8)
bsd_mtis    <- c(2,4,6)
maxent_etas <- c(1)
store_assignments  <- TRUE
store_ratio_string <- "4:3:2:1"
store_Nsim         <- 100
store_tuning       <- 6
store_designs      <- c("BSD", "MWUD", "DLUD")

set.seed(123)

## ---- 3) Helpers
safe_metric <- function(out, key) {
  if (!is.null(out) && !is.null(out[[key]]) && length(out[[key]]) == 1) out[[key]] else NA_real_
}

#convert an S4 seq object (with @M and getProbMatrix()) into the list format
# that evaluate_assignments() can consume (assignments + probabilities).
seq_to_eval_list <- function(seq_obj, groups) {
  trt_num <- as.integer(seq_obj@M[1, ])  # one sequence (r=1)
  trt_fac <- factor(trt_num, levels = seq_along(groups), labels = groups)

  P_kn <- getProbMatrix(seq_obj)         # K x N
  list(
    assignments   = data.frame(treatment = trt_fac),
    probabilities = t(P_kn)              # N x K (Evaluation.R expects list-based designs this way)
  )
}
append_assignments <- function(assignments_all, aix, out, design, ratio, Ns, tuning) {
  ratio_string <- paste(ratio, collapse = ":")

  keep_this <- store_assignments &&
    ratio_string == store_ratio_string &&
    Ns == store_Nsim &&
    tuning == store_tuning &&
    design %in% store_designs

  if (!keep_this || is.null(out$Assignments)) {
    return(list(assignments_all = assignments_all, aix = aix))
  }

  tmp <- out$Assignments %>%
    mutate(
      Design      = design,
      Ratio       = ratio_string,
      Nsim        = Ns,
      TuningParam = tuning
    )

  if ("Simulation" %in% names(tmp)) {
    tmp <- tmp %>% rename(simulation = Simulation)
  }

  assignments_all[[aix]] <- tmp
  aix <- aix + 1L

  list(assignments_all = assignments_all, aix = aix)
}
## ---- 4) Wrapper factories (friend-style: pass function + params to evaluator)

# CRD wrapper (returns list(assignments, probabilities))
crd_wrapper <- function(p) {
  crd_allocation_function(N = p$N, ratio = p$ratio, K = p$K, groups = p$groups)
}

# MaxEnt wrapper factory because eta varies
maxent_wrapper_factory <- function(eta) {
  function(p) {
    maxent_allocation_function(N = p$N, ratio = p$ratio, K = p$K, groups = p$groups, eta = eta)
  }
}

# MWUD wrapper factory because alpha varies (uses genSeq)
mwud_wrapper_factory <- function(alpha) {
  function(p) {
    par_obj <- mwudPar(N = p$N, ratio = p$ratio, groups = p$groups, alpha = alpha)
    seq_obj <- genSeq(par_obj, r = 1, seed = as.numeric(sample(.Machine$integer.max, 1)))
    seq_to_eval_list(seq_obj, p$groups)
  }
}

# DLUD wrapper factory because a varies
dlud_wrapper_factory <- function(a) {
  function(p) {
    par_obj <- dludPar(N = p$N, ratio = p$ratio, groups = p$groups,
                       immigration_rate = a, init_seed = 1, stochastic_round = FALSE)
    seq_obj <- genSeq(par_obj, r = 1, seed = as.numeric(sample(.Machine$integer.max, 1)))
    seq_to_eval_list(seq_obj, p$groups)
  }
}

# BSD wrapper factory because mti varies
bsd_wrapper_factory <- function(mti) {
  function(p) {
    par_obj <- bsdPar(N = p$N, mti = mti, K = p$K, ratio = p$ratio, groups = p$groups)
    seq_obj <- genSeq(par_obj, r = 1, seed = as.numeric(sample(.Machine$integer.max, 1)))
    seq_to_eval_list(seq_obj, p$groups)
  }
}
afi_steps <- data.frame()
## ---- storage for selected assignment-level allocations
assignments_all <- list()
aix <- 1L
## ---- 5) Main run loop -> results dataframe
results <- data.frame()
for (ratio in ratios) {

  K      <- length(ratio)
  groups <- LETTERS[1:K]

  # params object passed into evaluate_assignments
  params <- list(N = N, K = K, ratio = ratio, groups = groups)
  for (Ns in Nsim) {

    ## ---- CRD
    out <- evaluate_assignments(
      allocation_function = crd_wrapper,
      params = params,
      Nsim   = Ns,
      ratio  = ratio,
      groups = groups,
      method = "crd",
      plot   = FALSE
    )

    results <- rbind(results, data.frame(
      Design      = "CRD",
      Ratio       = paste(ratio, collapse=":"),
      Nsim        = Ns,
      TuningParam = NA_real_,
      Imbalance   = safe_metric(out, "Imbalance"),
      MPM         = safe_metric(out, "MPM"),

      FI_n        = safe_metric(out, "FI_n"),
      AFI         = safe_metric(out, "AFI"),
      UI          = safe_metric(out, "UI"),
      UR          = safe_metric(out, "UR"),
      G           = safe_metric(out, "G")
    ))

    ## ---- MaxEnt (eta grid)
    for (eta in maxent_etas) {
      out <- evaluate_assignments(
        allocation_function = maxent_wrapper_factory(eta),
        params = params,
        Nsim   = Ns,
        ratio  = ratio,
        groups = groups,
        method = "maxent",
        plot   = FALSE
      )

      results <- rbind(results, data.frame(
        Design      = "MaxEnt",
        Ratio       = paste(ratio, collapse=":"),
        Nsim        = Ns,
        TuningParam = eta,
        Imbalance   = safe_metric(out, "Imbalance"),
        MPM         = safe_metric(out, "MPM"),

        FI_n        = safe_metric(out, "FI_n"),
        AFI         = safe_metric(out, "AFI"),
        UI          = safe_metric(out, "UI"),
        UR          = safe_metric(out, "UR"),
        G           = safe_metric(out, "G")
      ))
    }

    ## ---- MWUD (alpha grid)
    for (alpha in mwud_alphas) {
      out <- evaluate_assignments(
        allocation_function = mwud_wrapper_factory(alpha),
        params = params,
        Nsim   = Ns,
        ratio  = ratio,
        groups = groups,
        method = "mwud",
        plot   = FALSE
      )
      tmp_assign <- append_assignments(assignments_all, aix, out, "MWUD", ratio, Ns, alpha)
      assignments_all <- tmp_assign$assignments_all
      aix <- tmp_assign$aix
      results <- rbind(results, data.frame(
        Design      = "MWUD",
        Ratio       = paste(ratio, collapse=":"),
        Nsim        = Ns,
        TuningParam = alpha,
        Imbalance   = safe_metric(out, "Imbalance"),
        MPM         = safe_metric(out, "MPM"),
        FI_n        = safe_metric(out, "FI_n"),
        AFI         = safe_metric(out, "AFI"),
        UI          = safe_metric(out, "UI"),
        UR          = safe_metric(out, "UR"),
        G           = safe_metric(out, "G")
      ))

    }

    ## ---- DLUD (a grid)
    for (a in dlud_as) {
      out <- evaluate_assignments(
        allocation_function = dlud_wrapper_factory(a),
        params = params,
        Nsim   = Ns,
        ratio  = ratio,
        groups = groups,
        method = "dlud",
        plot   = FALSE
      )
      tmp_assign <- append_assignments(assignments_all, aix, out, "DLUD", ratio, Ns, a)
      assignments_all <- tmp_assign$assignments_all
      aix <- tmp_assign$aix

      results <- rbind(results, data.frame(
        Design      = "DLUD",
        Ratio       = paste(ratio, collapse=":"),
        Nsim        = Ns,
        TuningParam = a,
        Imbalance   = safe_metric(out, "Imbalance"),
        MPM         = safe_metric(out, "MPM"),

        FI_n        = safe_metric(out, "FI_n"),
        AFI         = safe_metric(out, "AFI"),
        UI          = safe_metric(out, "UI"),
        UR          = safe_metric(out, "UR"),
        G           = safe_metric(out, "G")
      ))



    }
    ## ---- BSD (mti grid)
    for (mti in bsd_mtis) {
      out <- evaluate_assignments(
        allocation_function = bsd_wrapper_factory(mti),
        params = params,
        Nsim   = Ns,
        ratio  = ratio,
        groups = groups,
        method = "bsd",
        plot   = FALSE
      )
      tmp_assign <- append_assignments(assignments_all, aix, out, "BSD", ratio, Ns, mti)
      assignments_all <- tmp_assign$assignments_all
      aix <- tmp_assign$aix
      results <- rbind(results, data.frame(
        Design      = "BSD",
        Ratio       = paste(ratio, collapse=":"),
        Nsim        = Ns,
        TuningParam = mti,
        Imbalance   = safe_metric(out, "Imbalance"),
        MPM         = safe_metric(out, "MPM"),

        FI_n        = safe_metric(out, "FI_n"),
        AFI         = safe_metric(out, "AFI"),
        UI          = safe_metric(out, "UI"),
        UR          = safe_metric(out, "UR"),
        G           = safe_metric(out, "G")
      ))


    }

  }



}


print(results)

write.csv(results,
          file = "R/CRD_MaxEnt_MWUD_DLUD_BSD_results.csv",
          row.names = FALSE)


if (length(assignments_all) > 0) {
  assignments_long_all <- dplyr::bind_rows(assignments_all)

  write.csv(
    assignments_long_all,
    file = "R/assignments_long_all.csv",
    row.names = FALSE
  )

  cat("Saved assignments_long_all.csv with", nrow(assignments_long_all), "rows.\n")
  print(head(assignments_long_all))
} else {
  cat("No assignment-level data stored. Check store_ratio_string, store_Nsim, store_tuning, and store_designs.\n")
}


#### plots


results <- read.csv(
  "R/CRD_MaxEnt_MWUD_DLUD_BSD_results.csv",
  stringsAsFactors = FALSE
)

assignments_df <- read.csv(
  "R/assignments_long_all.csv",
  stringsAsFactors = FALSE
)

p_imb <- plot_imbalance_across_designs(results)
print(p_imb)

p_G <- plot_G_across_designs(results, chosen_param = 2)
print(p_G)

p_tradeoff <- plot_ui_ur_tradeoff(
  results,
  ratio0 = "37:21:21:21",
  Nsim0 = 100
)
print(p_tradeoff)

p_bsd_arp <- plot_bsd_arp(
  K = 3,
  N = 200,
  Nsim = 10000,
  ratio = c(1, 1, 1),
  groups = as.character(1:3),
  mti = 2,
  seed = 123
)
print(p_bsd_arp)

p_imb_time <- plot_imbalance_over_time(
  assignments_df,
  ratio0 = "4:3:2:1",
  Nsim0 = 100,
  tuning0 = 6
)
print(p_imb_time)


p_afi <- plot_afi_curves(
  ratio_vec = c(4, 3, 2, 1),
  tuning_vals = c(2, 8),
  N_total = 200,
  Nsim_plot = 100
)
print(p_afi)
