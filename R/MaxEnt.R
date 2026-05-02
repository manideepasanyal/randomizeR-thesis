# -----------------------------------------------------------
# Function: allocation_prb_maxent
# Purpose:  Calculate adaptive allocation probabilities for MaxEnt
# -----------------------------------------------------------
# Arguments:
#   target   : numeric vector of target allocation weights (e.g. c(3,1,1))
#   N        : integer vector of current treatment counts (length K)
#   eta      : numeric ∈ [0, 1], entropy–balance control parameter
# Returns:
#   Numeric vector of probabilities summing to 1
# -----------------------------------------------------------

# Allocation probabilities for MaxEnt
allocation_prb_maxent <- function(rnd, N_counts, tol_var = 1e-16) {
  # rnd$list: list(target = ratio, parameter = eta)
  w   <- as.numeric(rnd$target)      # target allocation counts
  eta <- as.numeric(rnd$parameter)   # η in [0,1]
  
  ntrt <- length(w)
  N    <- as.numeric(N_counts)
  
  # current subject index
  j <- sum(N) + 1L
  
  # target allocation proportions ρ
  rho <- w / sum(w)
  
  # ---- compute hypothetical lack of balance B(k) ----
  B  <- numeric(ntrt)
  N1 <- N
  
  for (k in seq_len(ntrt)) {
    N1[k] <- N1[k] + 1
    B[k]  <- max(abs(N1 / j - rho))
    N1[k] <- N1[k] - 1
  }
  
  # handle var(B) edge case (length 1 → NA)
  varB <- stats::var(B)
  if (is.na(varB) || varB <= tol_var) {
    # if (almost) all B are equal → just target proportions
    return(rho)
  }
  
  # f(μ) as in the Julia version
  f_mu <- function(mu) {
    exp_term <- exp(-mu * B)
    num <- sum(B * rho * exp_term)
    den <- sum(rho * exp_term)
    
    min(B) * eta + (1 - eta) * sum(rho * B) - num / den
  }
  
  # solve f(μ) = 0
  mu <- tryCatch(
    {
      stats::uniroot(
        f_mu,
        lower = -1,
        upper =  1,
        extendInt = "yes",
        tol = 1e-10
      )$root
    },
    error = function(e) {
      # fallback: μ = 0 if solver misbehaves
      0
    }
  )
  
  # final probabilities
  p_num <- rho * exp(-mu * B)
  p <- p_num / sum(p_num)
  p
}


maxent_allocation_function <- function(N, ratio, K, groups, eta) {
  # N: total number of subjects
  # ratio: target allocation counts, length K
  # K: number of arms
  # groups: group labels (length K)
  # eta: MaxEnt parameter in [0,1]
  
  # current counts
  counts <- rep(0L, K)
  names(counts) <- groups
  
  # store results
  assignments <- character(N)
  prob_matrix <- matrix(
    NA_real_,
    nrow = N,
    ncol = K,
    dimnames = list(NULL, groups)
  )
  
  for (m in seq_len(N)) {
    rnd <- list(target = ratio, parameter = eta)
    
    # probabilities given current counts
    probs <- allocation_prb_maxent(rnd, counts)
    
    # draw treatment for subject m
    trt <- sample(groups, size = 1, prob = probs)
    
    assignments[m] <- trt
    prob_matrix[m, ] <- probs
    
    # update counts
    counts[trt] <- counts[trt] + 1L
  }
  
  structure(
    list(
      assignments   = data.frame(treatment = assignments),
      probabilities = prob_matrix
    ),
    class = "maxentPar"
  )
}