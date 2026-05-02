#' @include randSeq.R
NULL


# Class mwudSeq for deterministic sequences

#' @rdname mwudSeq-class
setClass("mwudSeq", slots = c(alpha = "numeric"), contains = "randSeq")

# Class rMwudSeq for random sequences
#' @rdname rMwudSeq-class
setClass("rMwudSeq", contains = c("rRandSeq", "mwudSeq"))


methods::setMethod(
  f = "getProbMatrix",
  signature = "rMwudSeq",
  definition = function(obj) {
    M <- obj@M
    N <- obj@N
    K <- obj@K
    rho <- obj@ratio / sum(obj@ratio)
    alpha <- obj@alpha
    
    P <- matrix(NA_real_, nrow = K, ncol = N)
    
    for (m in seq_len(N)) {
      n_k <- if (m == 1) rep(0L, K) else tabulate(M[1:(m-1)], nbins = K)
      p <- pmax(alpha * rho - n_k + sum(n_k) * rho, 0)
      P[, m] <- p / sum(p)
    }
    P
  }
)




