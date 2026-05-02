#' @include randPar.R
NULL

validatebsdPar <- function(object) {
  errors <- character()
  mti <- object@mti
  if (length(mti) != 1L) errors <- c(errors, sprintf("mti has length %d. Should 
                                                be length one.", length(mti)))
  if (round(mti[1]) != mti) errors <- c(errors, sprintf("First element of mti 
                                          is %s. Should be an integer.", mti))
  if (mti[1] < 0) errors <- c(errors, "mti must be a non-negative integer")
  if (length(errors) == 0) TRUE else errors
}

setClass("bsdPar",
         slots    = c(mti = "numeric"),
         contains = "randPar",
         validity = validatebsdPar
)

#' @export
bsdPar <- function(N, mti, K = 2, ratio = rep(1, K), groups = LETTERS[1:K]) {
  new("bsdPar", N = N, mti = mti, K = K, ratio = ratio, groups = groups)
}

#Generalized BSD sampler 
bsdRand <- function(N, mti, K, ratio) {
  R <- integer(N)
  allocations <- integer(K)
  rho <- ratio / sum(ratio)
  
  for (m in 1:N) {
    expected <- rho * (m - 1)
    diff <- allocations - expected
    
    if (max(abs(diff)) < mti) {
      R[m] <- sample.int(K, 1, prob = rho)
    } else {
      under <- which(diff == min(diff))
      R[m] <- if (length(under) == 1L) under else sample(under, 1)
    }
    allocations[R[m]] <- allocations[R[m]] + 1L
  }
  R
}





