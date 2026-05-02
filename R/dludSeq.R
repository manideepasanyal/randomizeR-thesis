
#' @include randSeq.R
NULL


#' @rdname dludSeq-class
setClass("dludSeq",
         slots = c(
           immigration_rate = "numeric",
           init_seed = "numeric",
           stochastic_round = "logical",
           probMatrix = "matrix"   
         ),
         contains = "randSeq")

#' @rdname rDludSeq-class
setClass("rDludSeq",
         contains = c("rRandSeq", "dludSeq"))

# Method to extract design label
methods::setMethod(
  f = "getProbMatrix",
  signature = "rDludSeq",
  definition = function(obj) {
    obj@probMatrix   
  }
)





