set.seed(123)


ratio_list <- list(
  c(1,1,1,1),
  c(2,1,1,1),
  c(4,3,2,1),
  c(37,21,21,21)
  
)
ratio_names <- c("1:1:1:1","2:1:1:1","4:3:2:1","37:21:21:21")

# 2 and 6
tuning_vals <- c(2, 6)

#n and nsim for this block 
N_eval    <- 200
Nsim_eval <- 300
alpha_sig <- 0.05

# Effect under H1 
mu_H1_base <- c(0, 0.1, 0.3, 0.7)

# Extract allocation as integer 1..K 
# (works if treatment is "A/B/C" or 1/2/3)

get_alloc_int_from_wrapper <- function(allocation_function, p, groups) {
  seq_obj <- allocation_function(p)
  
  trt <- seq_obj$assignments$treatment
  if (is.factor(trt)) trt <- as.character(trt)
  
  if (is.numeric(trt) || is.integer(trt)) {
    alloc_int <- as.integer(trt)
  } else {
    alloc_int <- as.integer(factor(trt, levels = groups))
  }
  
  if (anyNA(alloc_int)) stop("Treatment coding produced NA (labels don't match groups).")
  alloc_int
}


# One design evaluation (type I + power) for a given ratio/tuning

evaluate_design_simple <- function(allocation_function, ratio, N, Nsim,
                                   model_type = c("without_drift","with_drift"),
                                   alpha = 0.05,
                                   mu_H0 = NULL, mu_H1 = NULL) {
  model_type <- match.arg(model_type)
  
  K <- length(ratio)
  groups <- LETTERS[1:K]
  p <- list(N = N, K = K, ratio = ratio, groups = groups)
  
  if (is.null(mu_H0)) mu_H0 <- rep(0, K)
  if (is.null(mu_H1)) mu_H1 <- mu_H1_base[seq_len(K)]
  
  run_test <- function(mu_vec) {
    reject <- logical(Nsim)
    
    for (s in seq_len(Nsim)) {
      alloc_int <- get_alloc_int_from_wrapper(allocation_function, p, groups)
      
      j <- seq_len(N)
      drift <- if (model_type == "with_drift") j / N else 0
      
      y <- mu_vec[alloc_int] + drift + rnorm(N, 0, 1)
      
      # IMPORTANT: treat alloc as a factor (groups irrelevant now, since alloc_int is 1..K)
      pval <- summary(aov(y ~ factor(alloc_int)))[[1]][["Pr(>F)"]][1]
      reject[s] <- isTRUE(pval < alpha)
    }
    mean(reject)
  }
  
  c(type1 = run_test(mu_H0), power = run_test(mu_H1))
}


# Driver: BSD/MWUD/DLUD only, tuning {2,6}, local ratios
type1_power_results <- data.frame()

for (r_idx in seq_along(ratio_list)) {
  ratio <- ratio_list[[r_idx]]
  ratio_label <- if (!is.null(ratio_names) && length(ratio_names) >= r_idx) ratio_names[r_idx] else paste(ratio, collapse=":")
  
  for (tval in tuning_vals) {
    
    # MWUD
    out_wo <- evaluate_design_simple(mwud_wrapper_factory(tval), ratio, N_eval, Nsim_eval, "without_drift", alpha_sig)
    out_w  <- evaluate_design_simple(mwud_wrapper_factory(tval), ratio, N_eval, Nsim_eval, "with_drift",    alpha_sig)
    
    type1_power_results <- rbind(type1_power_results, data.frame(
      Design="MWUD", TuningParam=tval, Ratio=ratio_label,
      TypeI_without_drift=round(out_wo["type1"], 4),
      Power_without_drift=round(out_wo["power"], 4),
      TypeI_with_drift=round(out_w["type1"], 4),
      Power_with_drift=round(out_w["power"], 4)
    ))
    
    # DLUD
    out_wo <- evaluate_design_simple(dlud_wrapper_factory(tval), ratio, N_eval, Nsim_eval, "without_drift", alpha_sig)
    out_w  <- evaluate_design_simple(dlud_wrapper_factory(tval), ratio, N_eval, Nsim_eval, "with_drift",    alpha_sig)
    
    type1_power_results <- rbind(type1_power_results, data.frame(
      Design="DLUD", TuningParam=tval, Ratio=ratio_label,
      TypeI_without_drift=round(out_wo["type1"], 4),
      Power_without_drift=round(out_wo["power"], 4),
      TypeI_with_drift=round(out_w["type1"], 4),
      Power_with_drift=round(out_w["power"], 4)
    ))
    
    # BSD
    out_wo <- evaluate_design_simple(bsd_wrapper_factory(tval), ratio, N_eval, Nsim_eval, "without_drift", alpha_sig)
    out_w  <- evaluate_design_simple(bsd_wrapper_factory(tval), ratio, N_eval, Nsim_eval, "with_drift",    alpha_sig)
    
    type1_power_results <- rbind(type1_power_results, data.frame(
      Design="BSD", TuningParam=tval, Ratio=ratio_label,
      TypeI_without_drift=round(out_wo["type1"], 4),
      Power_without_drift=round(out_wo["power"], 4),
      TypeI_with_drift=round(out_w["type1"], 4),
      Power_with_drift=round(out_w["power"], 4)
    ))
  }
}

rownames(type1_power_results) <- NULL
print(type1_power_results)

write.csv(type1_power_results, "Type1_Power_BSD_MWUD_DLUD_tuning_2_6.csv", row.names = FALSE)
type1_power_results%>% type1_power_results[""]