source("R/Crd.R")
source("R/MaxEnt.R")

# Metric Calculation Functions



# Average Standard Deviation

calculate_asd <- function(assignments, Nsim, N, K, ratio, groups) {
  if (Nsim == 1) return(rep(NA, N))
  asd_values <- numeric(N)
  all_proportions <- array(0, dim = c(Nsim, N, K))
  
  for (sim in 1:Nsim) {
    sim_data <- assignments[assignments$Simulation == sim, ]
    counts <- integer(K)
    for (m in 1:N) {
      idx <- which(groups == as.character(sim_data$treatment[m]))
      counts[idx] <- counts[idx] + 1
      all_proportions[sim, m, ] <- counts / m
    }
  }
  
  for (m in 1:N) {
    group_props <- matrix(all_proportions[, m, ], nrow = Nsim, ncol = K)
    sd_values <- apply(group_props, 2, sd)
    asd_values[m] <- sqrt(m * sum(sd_values^2))
  }
  return(asd_values)
}


# Imbalance

calculate_imbalance <- function(assignments, ratio) {
  n <- nrow(assignments)
  Ni <- table(assignments$treatment)
  target <- n * ratio / sum(ratio)
  sqrt(sum((Ni - target)^2))
}

# Randomness

calculate_FI_n <- function(P, target_proportions) {
  if (nrow(P) != length(target_proportions)) stop("P row count mismatch")
  FIs <- apply(P, 2, function(P_j) sqrt(sum((P_j - target_proportions)^2)))
  mean(FIs)
}

calculate_pi_matrix <- function(assignments, N, K, groups) {
  pi_matrix <- matrix(0, nrow = K, ncol = N)
  for (j in 1:N) {
    for (i in 1:K) {
      # number of patients allocated to trt i after jth allocation
      count_ij <- sum(assignments$treatment[1:j] == groups[i])
      pi_matrix[i, j] <- count_ij / j
    }
  }
  return(pi_matrix)
}

# Balance/Imbalance Ratio

calculate_performance_measures <- function(
    MPM_n_xi, MPM_n_MaxEnt, MPM_n_CRD,
    AFI_n_xi, AFI_n_MaxEnt, AFI_n_CRD,
    wI = 1, wR = 1
) {
  UI <- (MPM_n_xi - MPM_n_MaxEnt) / (MPM_n_CRD - MPM_n_MaxEnt)
  UR <- (AFI_n_xi) / (AFI_n_MaxEnt)
  
  G <- sqrt((wI * UI)^2 + (wR * UR)^2) / sqrt(wI^2 + wR^2)
  
  return(list(UI = UI, UR = UR, G = G))
}

get_param <- function(params, name) {
  if (isS4(params)) slot(params, name) else params[[name]]
}



# Reference Simulation Wrapper


simulate_reference_metrics <- function(
    Nsim, ratio, groups, N, K,
    ref_method, allocation_func,
    eta_maxent = 1  
) {
  cat("simulate_reference_metrics called for", ref_method, "\n")
  params <- list(N = N, ratio = ratio, K = K, groups = groups)
  class(params) <- paste0(ref_method, "Par")
  
  if (ref_method == "maxent") {
    allocation_wrapper <- function(p) {
      allocation_func(p$N, p$ratio, p$K, p$groups, eta = eta_maxent)
    }
  } else {
    allocation_wrapper <- function(p) {
      allocation_func(p$N, p$ratio, p$K, p$groups)
    }
  }
  
  result <- evaluate_assignments(allocation_wrapper, params, Nsim, ratio, groups, ref_method, plot = FALSE)
  list(MPM = result$MPM, AFI = result$AFI)
  
}


# Main Evaluation Function

evaluate_assignments <- function(allocation_function,
                                 params,
                                 Nsim,
                                 ratio,
                                 groups,
                                 method,
                                 eta_maxent_ref = 1,
                                 plot = TRUE) {

  # Basic setup

  N   <- get_param(params, "N")        # trial size n
  K   <- get_param(params, "K")        # number of treatments
  rho <- ratio / sum(ratio)            # target allocation proportions
  

  # Reference designs (CRD and MaxEnt(η = 1))
  # Only simulate them if the current design is NOT itself CRD or MaxEnt

  if (!(method %in% c("crd", "maxent"))) {
    ref_crd <- simulate_reference_metrics(
      Nsim         = Nsim,
      ratio        = ratio,
      groups       = groups,
      N            = N,
      K            = K,
      ref_method   = "crd",
      allocation_func = crd_allocation_function
    )
    
    ref_maxent <- simulate_reference_metrics(
      Nsim         = Nsim,
      ratio        = ratio,
      groups       = groups,
      N            = N,
      K            = K,
      ref_method   = "maxent",
      allocation_func = maxent_allocation_function,
      eta_maxent   = eta_maxent_ref
    )
    
  } else {
    ref_crd   <- NULL
    ref_maxent <- NULL
    
  }
  

  # Containers for all simulations

  # Imb_s(m) and FI_s(m) for each simulation s and step m
  Imb_all <- matrix(NA_real_, nrow = Nsim, ncol = N)  # rows: sims, cols: m
  FI_all  <- matrix(NA_real_, nrow = Nsim, ncol = N)
  
  # final-step values Imb_s(n), FI_s(n)
  Imb_final <- numeric(Nsim)
  FI_final  <- numeric(Nsim)
  
  assignments_list <- vector("list", Nsim)
  
  last_prb_mat <- NULL
  

  # Main simulation

  for (sim in seq_len(Nsim)) {
    
    # 1) Generate one randomization sequence for this design
    p <- list(
      N      = get_param(params, "N"),
      K      = get_param(params, "K"),
      ratio  = ratio,
      groups = groups
    )
    seq_obj <- allocation_function(p)
    
    # Extract assignments (handles S4 or list)
    assignments_data <- if (isS4(seq_obj)) slot(seq_obj, "assignments") else seq_obj$assignments
    assignments <- data.frame(
      patient   = 1:N,
      treatment = factor(assignments_data$treatment, levels = groups)
    )
    assignments_list[[sim]] <- assignments
    
    # 2) Imbalance by step: Imb_sim(m) for m = 1..N
    Imb_m <- numeric(N)
    for (m in seq_len(N)) {
      partial <- assignments[1:m, ]
      Ni_m <- table(factor(partial$treatment, levels = groups))
      
      if (K == 2 && ratio[1] == ratio[2]) {
        Imb_m[m] <- as.numeric(Ni_m[1] - Ni_m[2])
      } else {
        Imb_m[m] <- sqrt(sum((Ni_m - m * rho)^2))
      }
    }
    Imb_all[sim, ]  <- Imb_m
    Imb_final[sim]  <- Imb_m[N]  
    
    # 3) Conditional probabilities P(j) and FI_sim(m)

    # For FI/AFI use the design's conditional randomization
    # probabilities P(j), not empirical proportions.
    if (isS4(seq_obj)) {
      # prb_mat: K x N, column j = P(j) = (P_1(j), ..., P_K(j))
      prb_mat <- getProbMatrix(seq_obj)
      
    } else if (!is.null(seq_obj$probabilities)) {
      # CRD / MaxEnt / other list-based designs:
      # probabilities is N x K (rows = patients, cols = groups)
      # we need K x N (rows = groups, cols = patients)
      prb_mat <- t(seq_obj$probabilities)
      
    } else {
      stop(
        "No probability matrix available: expected either an S4 sequence with getProbMatrix(), ",
        "or a list with a 'probabilities' matrix (e.g. crdPar, maxentPar)."
      )
    }
    
    last_prb_mat <- prb_mat
    
    # Build matrix with rho repeated for each column j
    
    rho_mat <- matrix(rho, nrow = nrow(prb_mat), ncol = ncol(prb_mat), byrow = FALSE)
    
    # FI_j = sqrt(sum_i (P_i(j) - rho_i)^2)  
    diff_mat <- prb_mat - rho_mat          # K x N
    FI_j     <- sqrt(colSums(diff_mat^2))        # length N
    
    # FI_sim(m) = (1/m) Σ_{j=1}^m FI_j  for m = 1..N
    FI_m <- cumsum(FI_j) / seq_len(N)
    
    FI_all[sim, ] <- FI_m
    FI_final[sim] <- FI_m[N]      # FI(n) for this trial
  }
  
  # Step-wise expectations: MPM(m) and AFI(m)
  MPM_by_step <- colMeans(Imb_all)  
  
  AFI_by_step <- colMeans(FI_all)   
  
  
  MPM_design <- mean(MPM_by_step)
  
  AFI_design <- AFI_by_step[N]
  
  Imb_mean  <- mean(Imb_final)
  FI_n_mean <- mean(FI_final)
  

  # Unified imbalance UI, unified randomness UR, and overall G

  if (!(method %in% c("crd", "maxent"))) {
    perf <- calculate_performance_measures(
      MPM_design,       # MPM(n)_ξ
      ref_maxent$MPM,   # MPM(n)_MaxEnt
      ref_crd$MPM,      # MPM(n)_CRD
      AFI_design,       # AFI(n)_ξ
      ref_maxent$AFI,   # AFI(n)_MaxEnt
      ref_crd$AFI       # AFI(n)_CRD
    )
    UI <- perf$UI
    UR <- perf$UR
    G  <- perf$G
  } else if (method=="crd"){
    UI =1 
    UR =0 
    G <- NA_real_
  } else{
    UI =0 
    UR =1 
    G <- NA_real_
  }
  
  # ASD and plotting 
  assignments_df <- do.call(rbind, assignments_list)
  assignments_df$Simulation <- rep(seq_len(Nsim), each = N)
  
  asd_values <- calculate_asd(assignments_df, Nsim, N, K, ratio, groups)
  
  if (plot) {
    par(mfrow = c(1, 3))
    
    # Use last pi_matrix as an example trajectory
    matplot(t(last_prb_mat), type = "l", lty = 1, col = rainbow(K),
            main = paste("Treatment Allocation for", method),
            xlab = "Sample Size (m)", ylab = "Proportion Assigned", lwd = 2)
    legend("topright", legend = groups, col = rainbow(K), lty = 1, title = "Groups")
    
    barplot(table(assignments_df$treatment) / (N * Nsim),
            main = paste("Final Treatment Distribution for", method),
            col = rainbow(K), names.arg = groups, ylab = "Proportion")
    
    plot(1:N, asd_values, type = "l", col = "blue", lwd = 2,
         main = paste("ASD for", method),
         xlab = "Sample Size (m)", ylab = "ASD")
  }
  

  # Return all relevant metrics
  list(
    Assignments   = assignments_df,
    ASD           = asd_values,
    Imbalance     = Imb_mean,       # mean Imb(n) over simulations
    FI_n          = FI_n_mean,      # mean FI(n) over simulations
    MPM           = MPM_design,     # MPM(n)_ξ
    AFI           = AFI_design,     # AFI(n)_ξ
    MPM_by_step   = MPM_by_step,    # MPM(m), m = 1..N
    AFI_by_step   = AFI_by_step,    # AFI(m), m = 1..N
    UI            = UI,
    UR            = UR,
    G             = G,
    ref_crd     = ref_crd,
    ref_maxent  = ref_maxent
  )
}