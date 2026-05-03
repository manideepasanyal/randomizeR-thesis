source("R/randPar.R")
#' @include randPar.R
NULL

# Validation function for mwudPar
validateMwudPar <- function(object) {
  errors <- character()

  if (length(object@alpha) != 1 || object@alpha <= 0 || object@alpha != as.integer(object@alpha)) {
    errors <- c(errors, "alpha must be a single positive integer.")
  }

  if (length(object@ratio) != object@K) {
    errors <- c(errors, "Length of 'ratio' must match 'K'.")
  }

  if (length(object@groups) != object@K) {
    errors <- c(errors, "Length of 'groups' must match 'K'.")
  }

  if (length(errors) == 0) TRUE else errors
}

# mwudPar class
#' @rdname mwudPar-class
#' @export
#'
setClass("mwudPar",
         slots = c(N = "numeric", K = "numeric", ratio = "numeric", groups = "character", alpha = "numeric"),
         validity = validateMwudPar)

# Constructor for mwudPar

#' Create a mwudPar object
#'
#' @param N Total sample size
#' @param ratio Allocation ratio
#' @param groups Group labels
#' @param alpha Weighting parameter (default = 1)
#'
#' @export
mwudPar <- function(N, ratio, groups = LETTERS[1:length(ratio)], alpha = 1) {
  K <- length(ratio)
  new("mwudPar", N = N, K = K, ratio = ratio, groups = groups, alpha = alpha)
}

# MWUD Randomization Function (core)
mwudRand <- function(N, K, ratio, alpha) {
  rho <- ratio / sum(ratio)
  trt <- integer(N)
  counts <- integer(K)

  for (m in 1:N) {
    p <- pmax(alpha * rho - counts + sum(counts) * rho, 0)
    if (sum(p) <= 0) p <- rho else p <- p / sum(p)   # safety fallback
    trt[m] <- sample.int(K, size = 1, prob = p)
    counts[trt[m]] <- counts[trt[m]] + 1L
  }
  trt
}



#' @rdname getDesign

setMethod("genSeq", signature(obj = "mwudPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            simMatrix <- t(sapply(seq_len(r), function(x) {
              mwudRand(obj@N, obj@K, obj@ratio, obj@alpha)
            }))
            new("rMwudSeq", M = simMatrix, N = obj@N, K = obj@K,
                ratio = obj@ratio, groups = obj@groups, alpha = obj@alpha, seed = seed)
          })

setMethod("genSeq", signature(obj = "mwudPar", r = "numeric", seed = "missing"),
          function(obj, r) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            simMatrix <- t(sapply(seq_len(r), function(x) {
              mwudRand(obj@N, obj@K, obj@ratio, obj@alpha)
            }))
            new("rMwudSeq", M = simMatrix, N = obj@N, K = obj@K,
                ratio = obj@ratio, groups = obj@groups, alpha = obj@alpha, seed = seed)
          })
