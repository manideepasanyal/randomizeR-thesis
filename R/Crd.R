# ------------------------------------------------------------------
# Completely Randomized Design (CRD)
# ------------------------------------------------------------------

# Allocation probability function for CRD
allocation_prb_crd <- function(rnd) {
  # rnd: a CRD object, containing 'target' (allocation weights)
  w <- rnd$target
  prb <- w / sum(w)
  return(prb)
}

# Allocation function for CRD
crd_allocation_function <- function(N, ratio, K, groups) {
  probs <- allocation_prb_crd(list(target = ratio))
  
  # Assignments
  assignments <- sample(groups, N, replace = TRUE, prob = probs)
  
  # Constant probability matrix for FI_n check
  prob_matrix <- matrix(rep(probs, each = N),
                        nrow = N, ncol = K,
                        byrow = FALSE,
                        dimnames = list(NULL, groups))
  
  structure(
    list(assignments   = data.frame(treatment = assignments),
         probabilities = prob_matrix),
    class = "crdPar"
  )
}