source("R/randSeq.R")  
#' @include randSeq.R
NULL

# Class bsdSeq                                


# Representation of sequences for the Big Stick Design
# @description Stores BSD randomization sequences with design parameters.
# @slot N total number of patients included in the trial.
# @slot mti Maximum tolerated imbalance (MTI).
# @slot M matrix containing randomization sequences of length \code{N} in rows.

setClass("bsdSeq",
         slots = c(mti = "numeric"),
         contains = "randSeq")


setClass("rBsdSeq",
         contains = c("rRandSeq", "bsdSeq"))


# Methods for bsdSeq


#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "bsdSeq"),
          function(obj) {
            M    <- obj@M
            K    <- obj@K
            N    <- obj@N
            rho  <- obj@ratio / sum(obj@ratio)
            MTI  <- obj@mti
            
            # Map 0/1 to 1/2 for K=2 if needed (back-compat with binary sequences)
            map01 <- function(x) if (K == 2 && all(x %in% c(0, 1))) x + 1 else x
            
            rowProb <- function(seqRow) {
              s <- as.integer(map01(seqRow))
              if (length(s) != N) return(NA_real_)   # basic guard
              
              alloc <- rep(0, K)
              p <- 1
              
              for (m in seq_len(N)) {
                chosen <- s[m]
                if (is.na(chosen) || chosen < 1 || chosen > K) return(NA_real_)
                
                # Current deviation from target just before assigning patient m
                expected <- rho * (m - 1)
                diff <- alloc - expected
                
                if (max(abs(diff)) < MTI) {
                  # Random zone: assign with target probabilities
                  p <- p * rho[chosen]
                } else {
                  # Forcing zone: choose among most underrepresented arms
                  underfull <- which(diff == min(diff))
                  if (!(chosen %in% underfull)) return(0)  # violates BSD rule
                  p <- p * (1 / length(underfull))
                }
                
                # Update allocations after assigning m
                alloc[chosen] <- alloc[chosen] + 1
              }
              p
            }
            
            apply(M, 1, rowProb)
          }
)


setMethod("genSeq", signature(obj = "bsdPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rBsdSeq",
                M      = t(sapply(seq_len(r), function(i)
                  bsdRand(N = obj@N, mti = obj@mti, K = obj@K, ratio = obj@ratio))),
                N      = obj@N,
                mti    = obj@mti,
                K      = obj@K,
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          })

setMethod("genSeq", signature(obj = "bsdPar", r = "numeric", seed = "missing"),
          function(obj, r) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rBsdSeq",
                M      = t(sapply(seq_len(r), function(i)
                  bsdRand(N = obj@N, mti = obj@mti, K = obj@K, ratio = obj@ratio))),
                N      = obj@N,
                mti    = obj@mti,
                K      = obj@K,
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          })

setMethod("genSeq", signature(obj = "bsdPar", r = "missing", seed = "numeric"),
          function(obj, seed) {
            set.seed(seed)
            new("rBsdSeq",
                M      = t(bsdRand(N = obj@N, mti = obj@mti, K = obj@K, ratio = obj@ratio)),
                N      = obj@N,
                mti    = obj@mti,
                K      = obj@K,
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          })

setMethod("genSeq", signature(obj = "bsdPar", r = "missing", seed = "missing"),
          function(obj) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rBsdSeq",
                M      = t(bsdRand(N = obj@N, mti = obj@mti, K = obj@K, ratio = obj@ratio)),
                N      = obj@N,
                mti    = obj@mti,
                K      = obj@K,
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          })
methods::setMethod("getProbMatrix", signature = "rBsdSeq", function(obj) {
  M    <- obj@M
  K    <- obj@K
  N    <- obj@N
  rho  <- obj@ratio / sum(obj@ratio)
  MTI  <- obj@mti
  r    <- nrow(M)
  
  P_avg <- matrix(0, nrow = K, ncol = N)
  
  for (sim in seq_len(r)) {
    alloc <- rep(0L, K)
    
    for (m in seq_len(N)) {
      expected <- rho * (m - 1)
      diff <- alloc - expected
      
      if (max(abs(diff)) < MTI) {
        p <- rho
      } else {
        under <- which(diff == min(diff))
        p <- rep(0, K)
        p[under] <- 1 / length(under)
      }
      
      P_avg[, m] <- P_avg[, m] + p
      
      chosen <- as.integer(M[sim, m])
      alloc[chosen] <- alloc[chosen] + 1L
    }
  }
  
  P_avg / r
})

