#' @include randPar.R
NULL


# Validity check for dludPar

validateDludPar <- function(object) {
  errors <- character()
  
  if (length(object@immigration_rate) != 1 || object@immigration_rate <= 0) {
    errors <- c(errors, "immigration_rate must be a single positive number.")
  }
  
  if (length(object@init_seed) != 1 || object@init_seed < 0) {
    errors <- c(errors, "init_seed must be a non-negative integer.")
  }
  
  if (length(object@ratio) != object@K) {
    errors <- c(errors, "Length of 'ratio' must match K.")
  }
  
  if (length(object@groups) != object@K) {
    errors <- c(errors, "Length of 'groups' must match K.")
  }
  
  if (length(errors) == 0) TRUE else errors
}


# dludPar class

#' @rdname dludPar-class
#' @export
setClass("dludPar",
         slots = c(N = "numeric", K = "numeric", ratio = "numeric", groups = "character",
                   immigration_rate = "numeric", init_seed = "numeric",
                   stochastic_round = "logical"),
                  validity = validateDludPar)

# Constructor
#' Create a Drop-the-Loser urn design (dludPar)
#'
#' @param N Number of patients
#' @param ratio Allocation ratio (vector of K values)
#' @param groups Group labels (optional, defaults to LETTERS)
#' @param immigration_rate Immigration parameter (a)
#' @param init_seed Initial count for urn balls
#' @param stochastic_round Whether to use stochastic rounding (TRUE/FALSE)
#'
#' @return An object of class \code{dludPar}
#' @export
dludPar <- function(N, ratio, groups = LETTERS[1:length(ratio)],
                    immigration_rate = 1, init_seed = 1, stochastic_round = FALSE) {
  K <- length(ratio)
  new("dludPar", N = N, K = K, ratio = ratio, groups = groups,
      immigration_rate = immigration_rate, init_seed = init_seed,
      stochastic_round = stochastic_round)
}

# Core DLUD function
dludRand <- function(N, K, w, immigration_rate, init_seed = 0L, stochastic_round = FALSE) {
  add_int <- function(x) {
    if (stochastic_round) {
      floor(x) + rbinom(1, 1, x - floor(x))
    } else {
      as.integer(round(x))
    }
  }
  c0 <- 1L
  c  <- rep(add_int(init_seed), K)
  trt <- integer(N)
  P <- matrix(NA_real_, nrow = K, ncol = N)
  for (i in seq_len(N)) {
    repeat {
      total <- c0 + sum(c)
      u <- runif(1)
      p0 <- c0 / total
      if (u < p0) {
        for (j in 1:K) c[j] <- c[j] + add_int(immigration_rate * w[j])
      } else {
        P[, i] <- c / sum(c)
        v <- (u - p0) / (1 - p0)
        cum <- cumsum(c) / sum(c)
        arm <- which(v <= cum)[1]
        trt[i] <- arm
        c[arm] <- c[arm] - 1L
        break
      }
    }
  }
  list(trt = trt, P = P)
}


# genSeq methods


setMethod("genSeq", signature(obj = "dludPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            
            sims <- lapply(seq_len(r), function(x) {
              dludRand(obj@N, obj@K, obj@ratio,
                       obj@immigration_rate, obj@init_seed, obj@stochastic_round)
            })
            
            simMatrix <- t(sapply(sims, `[[`, "trt"))
                        if (r == 1) {
              Pstore <- sims[[1]]$P
            } else {
              Pstore <- Reduce(`+`, lapply(sims, `[[`, "P")) / r
            }
            
            new("rDludSeq", M = simMatrix, N = obj@N, K = obj@K,
                ratio = obj@ratio, groups = obj@groups,
                immigration_rate = obj@immigration_rate,
                init_seed = obj@init_seed,
                stochastic_round = obj@stochastic_round,
                probMatrix = Pstore,     
                seed = seed)
          })


setMethod("genSeq", signature(obj = "dludPar", r = "numeric", seed = "missing"),
          function(obj, r) {
            
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            
            sims <- lapply(seq_len(r), function(x) {
              dludRand(obj@N, obj@K, obj@ratio,
                       obj@immigration_rate, obj@init_seed, obj@stochastic_round)
            })
            
            simMatrix <- t(sapply(sims, `[[`, "trt"))  # r x N
            
            if (r == 1) {
              Pstore <- sims[[1]]$P
            } else {
              Pstore <- Reduce(`+`, lapply(sims, `[[`, "P")) / r
            }
            
            new("rDludSeq",
                M = simMatrix, N = obj@N, K = obj@K,
                ratio = obj@ratio, groups = obj@groups,
                immigration_rate = obj@immigration_rate,
                init_seed = obj@init_seed,
                stochastic_round = obj@stochastic_round,
                probMatrix = Pstore,       # NEW
                seed = seed)
          })


#' @rdname getDesign


