library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

theme_set(theme_bw(base_size = 12))

# --------------------------------------------------
# Helper
# --------------------------------------------------

parse_ratio <- function(r) {
  as.numeric(strsplit(r, ":", fixed = TRUE)[[1]])
}

# --------------------------------------------------
# 1. Imbalance across designs and ratios
# --------------------------------------------------

plot_imbalance_across_designs <- function(results) {
  
  res_imb <- results %>%
    filter(Design %in% c("BSD", "MWUD", "DLUD")) %>%
    mutate(
      Design = factor(Design, levels = c("BSD", "DLUD", "MWUD"))
    )
  
  ggplot(res_imb, aes(x = Ratio, y = Imbalance, fill = Design)) +
    geom_boxplot(alpha = 0.75, position = position_dodge(width = 0.8)) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Imbalance across Designs and Ratios",
      x = "Allocation ratio",
      y = "Imbalance",
      fill = "Design"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "plain")
    )
}

# --------------------------------------------------
# 2. Overall performance G
# --------------------------------------------------

plot_G_across_designs <- function(results, chosen_param = 2) {
  
  agg_G <- results %>%
    filter(TuningParam == chosen_param) %>%
    group_by(Ratio, Design) %>%
    summarise(
      G_mean = mean(G, na.rm = TRUE),
      G_sd   = sd(G, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(agg_G, aes(x = Design, y = G_mean, fill = Design)) +
    geom_col(width = 0.7) +
    geom_errorbar(
      aes(ymin = G_mean - G_sd, ymax = G_mean + G_sd),
      width = 0.15
    ) +
    facet_wrap(~ Ratio, nrow = 1) +
    labs(
      title = paste0("Overall Performance G across Designs (TuningParam = ", chosen_param, ")"),
      x = "Design",
      y = "G (aggregate index)"
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

# --------------------------------------------------
# 3. UI vs UR trade-off
# --------------------------------------------------

plot_ui_ur_tradeoff <- function(results, ratio0 = "37:21:21:21", Nsim0 = 100) {
  
  plot_df <- results %>%
    filter(Ratio == ratio0, Nsim == Nsim0) %>%
    mutate(
      TuningLabel = ifelse(
        Design %in% c("CRD", "MaxEnt"),
        "ref",
        as.character(TuningParam)
      ),
      TuningLabel = factor(TuningLabel, levels = c("ref", "2", "4", "6", "8")),
      Design = factor(Design, levels = c("BSD", "CRD", "DLUD", "MaxEnt", "MWUD"))
    )
  
  ggplot(
    plot_df,
    aes(x = UI, y = UR, color = Design, shape = TuningLabel)
  ) +
    geom_point(size = 3.5, alpha = 0.9) +
    labs(
      title = paste0(
        "Trade-off between Balance (UI) and Randomness (UR) (Ratio = ", ratio0, ")"
      ),
      x = "UI (0 = MaxEnt, 1 = CRD)",
      y = "UR (0 = CRD, 1 = MaxEnt)",
      color = "Design",
      shape = "Tuning"
    ) +
    scale_color_brewer(palette = "Set2") +
    scale_shape_manual(
      values = c(
        "ref" = 4,
        "2"   = 16,
        "4"   = 17,
        "6"   = 15,
        "8"   = 18
      )
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "plain")
    )
}

# --------------------------------------------------
# 4. BSD ARP plot
# --------------------------------------------------

plot_bsd_arp <- function(
    K = 3,
    N = 200,
    Nsim = 10000,
    ratio = c(1, 1, 1),
    groups = as.character(1:K),
    mti = 2,
    seed = 123
) {
  
  p0 <- 1 / K
  
  set.seed(seed)
  
  bsd_par <- bsdPar(
    N      = N,
    mti    = mti,
    K      = K,
    ratio  = ratio,
    groups = groups
  )
  
  bsd_seq <- genSeq(bsd_par, r = Nsim, seed = seed)
  M <- bsd_seq@M
  
  arp_df <- expand.grid(
    step = 1:N,
    treatment = 1:K
  ) %>%
    arrange(treatment, step) %>%
    mutate(
      pi_hat = purrr::map2_dbl(step, treatment, ~ mean(M[, .x] == .y)),
      treatment = factor(treatment, levels = 1:K, labels = groups)
    )
  
  se <- sqrt(p0 * (1 - p0) / Nsim)
  
  ci_df <- data.frame(
    step = 1:N,
    lower = p0 - 1.96 * se,
    upper = p0 + 1.96 * se,
    target = p0
  )
  
  ggplot() +
    geom_ribbon(
      data = ci_df,
      aes(x = step, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      alpha = 0.15,
      fill = "grey60"
    ) +
    geom_line(
      data = arp_df,
      aes(x = step, y = pi_hat, color = treatment),
      linewidth = 1
    ) +
    geom_hline(
      yintercept = p0,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    ) +
    labs(
      title = "Empirical unconditional allocation probabilities under BSD",
      x = "Allocation step",
      y = expression(hat(pi)[ij]),
      color = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right"
    )
}

# --------------------------------------------------
# 5. Imbalance over time
# --------------------------------------------------

plot_imbalance_over_time <- function(
    assignments_df,
    ratio0 = "4:3:2:1",
    Nsim0 = 10,
    tuning0 = 6
) {
  
  scen <- assignments_df %>%
    filter(
      Ratio == ratio0,
      Nsim == Nsim0,
      TuningParam == tuning0,
      Design %in% c("BSD", "DLUD", "MWUD")
    ) %>%
    mutate(
      Design = factor(Design, levels = c("BSD", "DLUD", "MWUD"))
    )
  
  ratio_vec <- parse_ratio(ratio0)
  rho <- ratio_vec / sum(ratio_vec)
  groups <- LETTERS[seq_along(rho)]
  
  scen <- scen %>%
    mutate(
      treatment = factor(treatment, levels = groups)
    )
  
  wide_cum <- scen %>%
    mutate(one = 1L) %>%
    pivot_wider(
      id_cols = c(Design, simulation, patient),
      names_from = treatment,
      values_from = one,
      values_fill = 0L
    ) %>%
    arrange(Design, simulation, patient) %>%
    group_by(Design, simulation) %>%
    mutate(across(all_of(groups), cumsum)) %>%
    ungroup()
  
  imb_over_time <- wide_cum %>%
    pivot_longer(
      cols = all_of(groups),
      names_to = "Arm",
      values_to = "N_im"
    ) %>%
    mutate(
      m = patient,
      arm_index = match(Arm, groups),
      expected = m * rho[arm_index],
      sq_diff = (N_im - expected)^2
    ) %>%
    group_by(Design, simulation, m) %>%
    summarise(
      Imb = sqrt(sum(sq_diff)),
      .groups = "drop"
    )
  
  ggplot(
    imb_over_time,
    aes(x = m, y = Imb, color = Design, fill = Design)
  ) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.1) +
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.15, color = NA) +
    scale_color_manual(values = c(
      "BSD"  = "#66C2A5",
      "DLUD" = "#FC8D62",
      "MWUD" = "#8DA0CB"
    )) +
    scale_fill_manual(values = c(
      "BSD"  = "#66C2A5",
      "DLUD" = "#FC8D62",
      "MWUD" = "#8DA0CB"
    )) +
    labs(
      title = "Imbalance over time (Imb)",
      subtitle = paste0(
        "Ratio = ", ratio0,
        " | Nsim = ", Nsim0,
        " | Tuning = ", tuning0
      ),
      x = "Allocation step (m)",
      y = expression(Imb(m)),
      color = "Design",
      fill = "Design"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "plain")
    )
}

# --------------------------------------------------
# 6. AFI curves
# --------------------------------------------------

plot_afi_curves <- function(
    ratio_vec = c(4, 3, 2, 1),
    tuning_vals = c(2, 8),
    N_total = 200,
    Nsim_plot = 100
) {
  
  get_curve <- function(design, tuning, ratio_vec, N, Nsim_scalar) {
    
    alloc_fun <- switch(
      design,
      BSD  = bsd_wrapper_factory(mti   = tuning),
      MWUD = mwud_wrapper_factory(alpha = tuning),
      DLUD = dlud_wrapper_factory(a     = tuning)
    )
    
    K <- length(ratio_vec)
    groups <- LETTERS[1:K]
    
    p <- list(
      N = N,
      K = K,
      ratio = ratio_vec,
      groups = groups
    )
    
    out <- evaluate_assignments(
      allocation_function = alloc_fun,
      params  = p,
      Nsim    = Nsim_scalar,
      ratio   = ratio_vec,
      groups  = groups,
      method  = design,
      plot    = FALSE
    )
    
    data.frame(
      Design = design,
      Tuning = paste0("t=", tuning),
      step   = seq_along(out$AFI_by_step),
      AFI    = out$AFI_by_step
    )
  }
  
  afi_6lines <- bind_rows(
    lapply(c("BSD", "MWUD", "DLUD"), function(design) {
      bind_rows(lapply(tuning_vals, function(t) {
        get_curve(design, t, ratio_vec, N_total, Nsim_plot)
      }))
    })
  )
  
  ggplot(
    afi_6lines,
    aes(
      x = step,
      y = AFI,
      linetype = Tuning,
      colour = Design,
      group = interaction(Design, Tuning)
    )
  ) +
    geom_line(linewidth = 1) +
    labs(
      title = "AFI by Allocation Step",
      subtitle = paste0(
        "Ratio = ", paste(ratio_vec, collapse = ":"),
        " | Tunings = {", paste(tuning_vals, collapse = ", "), "}",
        " | Nsim = ", Nsim_plot
      ),
      x = "Allocation step",
      y = "AFI(m)",
      linetype = "Tuning",
      colour = "Design"
    ) +
    theme_minimal()
}