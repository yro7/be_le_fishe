# Illustrations des différents modes et 
# vitesses de convergence de l'EMV pour processus de Weibull
# (cf Cocozza-Thivent, Propositions 3.7, 3.11 et 3.12)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Chargement des fonctions du modèle (simulation, EMV, IC, test) depuis power_law.R
if (!exists("mle_power_law")) {
  src_dir <- if (file.exists("src/power_law.R")) {
    "src"
  } else if (file.exists("power_law.R")) {
    "."
  } else {
    file.path(rprojroot::find_rstudio_root_file(), "src")
  }
  source(file.path(src_dir, "power_law.R"))
}

# 2. Fig 1 : Convergence presque-sûrement (Proposition 3.7)

#' Simule le comportement de l'EMV le long d'une même trajectoire quand t croît
#'
#' Proposition 3.7 : Lorsque t -> +∞, alpha_hat(t) -> alpha et beta_hat(t) -> beta p.s.
#'
#' @param alpha_true Vraie valeur de alpha
#' @param beta_true  Vraie valeur de beta
#' @param T_max      Horizon maximal pour la trajectoire
#' @param n_points   Nombre de points d'évaluation temporelle
#' @param seed       Graine aléatoire
#' @return data.frame contenant t, n, beta_hat, alpha_hat
sim_trajectory_as_convergence <- function(alpha_true, beta_true, T_max, n_points = 400, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  w_all <- sim_power_law_process(alpha_true, beta_true, T_max)
  while (length(w_all) < 10) {
    w_all <- sim_power_law_process(alpha_true, beta_true, T_max)
  }
  
  t_min <- max(w_all[3], T_max * 0.01)
  t_grid <- seq(t_min, T_max, length.out = n_points)
  
  res <- vector("list", length(t_grid))
  for (k in seq_along(t_grid)) {
    t_curr <- t_grid[k]
    w_sub <- w_all[w_all <= t_curr]
    est <- mle_power_law(w_sub, t_curr)
    res[[k]] <- data.frame(
      t = t_curr,
      n = est$n,
      beta_hat = est$beta_hat,
      beta_tilde = est$beta_tilde,
      alpha_hat = est$alpha_hat
    )
  }
  
  do.call(rbind, res)
}

#' Génère le graphique de convergence presque-sûre pour plusieurs trajectoires
plot_as_convergence <- function(alpha_true = 1.0, beta_true = 1.5, T_max = 200, n_traj = 3) {
  traj_list <- vector("list", n_traj)
  for (i in seq_len(n_traj)) {
    df_i <- sim_trajectory_as_convergence(alpha_true, beta_true, T_max, seed = 100 + i)
    df_i$traj_id <- paste0("Trajectoire ", i)
    traj_list[[i]] <- df_i
  }
  df_traj <- do.call(rbind, traj_list)
  
  theme_clean <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey30", size = 10),
      legend.position = "bottom"
    )
  
  p_beta <- ggplot(df_traj, aes(x = t, y = beta_hat, color = traj_id)) +
    geom_line(alpha = 0.85, linewidth = 0.7) +
    geom_hline(yintercept = beta_true, linetype = "dashed", color = "black", linewidth = 0.8) +
    labs(
      title = expression("Convergence presque-sûre de " ~ hat(beta)(t)),
      x = "Horizon d'observation t",
      y = expression(hat(beta)(t)),
      color = "Trajectoire"
    ) +
    theme_clean
  
  p_alpha <- ggplot(df_traj, aes(x = t, y = alpha_hat, color = traj_id)) +
    geom_line(alpha = 0.85, linewidth = 0.7) +
    geom_hline(yintercept = alpha_true, linetype = "dashed", color = "black", linewidth = 0.8) +
    labs(
      title = expression("Convergence presque-sûre de " ~ hat(alpha)(t)),
      x = "Horizon d'observation t",
      y = expression(hat(alpha)(t)),
      color = "Trajectoire"
    ) +
    theme_clean
  
  p_combined <- p_beta / p_alpha + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  return(p_combined)
}

# ------------------------------------------------------------------------------
# 3. Simulation Monte-Carlo pour l'étude asymptotique (Prop. 3.11 et 3.12)
# ------------------------------------------------------------------------------

#' Lance une simulation Monte-Carlo sur plusieurs horizons T
#'
#' Statistiques normalisées (Cocozza-Thivent 1997, Prop 3.11 & 3.12) :
#'   Z_beta        = sqrt(N(T)) * (beta_hat / beta - 1)      -> N(0, 1)
#'   Y_beta        = sqrt(N(T)) * (beta / beta_hat - 1)      -> N(0, 1)
#'   Z_alpha       = (sqrt(N(T)) / log(T)) * (alpha_hat / alpha - 1) -> N(0, 1)
#'   Z_alpha_log   = (sqrt(N(T)) / log(T)) * log(alpha_hat / alpha)  -> N(0, 1)
#'   Z_alpha_naive = sqrt(N(T)) * (alpha_hat / alpha - 1)   (variance explose en log^2(T) !)
#'
#' Remarque sur le choix d'alpha_true = 1 :
#'   Dans la preuve de Cocozza-Thivent (ligne 310), le terme de développement contient :
#'     - (1 - log(alpha) / log(t)) * sqrt(Lambda(t)) * (beta / beta_hat - 1)
#'   Lorsque alpha = 1, log(alpha) = 0 exactement, éliminant le biais de vitesse
#'   logarithmique O(1 / log(T)) à distance finie et révélant la variance unitaire.
#'
#' @param alpha_true Vrai alpha (défaut = 1.0)
#' @param beta_true  Vrai beta (défaut = 1.5)
#' @param T_values   Horizons T testés
#' @param M          Nombre de réplications par horizon
#' @param seed       Graine aléatoire
#' @return data.frame structuré
monte_carlo_mle <- function(alpha_true = 1.0, beta_true = 1.5,
                            T_values = c(10, 25, 60, 150),
                            M = 2000, seed = 2026) {
  set.seed(seed)
  
  records <- vector("list", length(T_values) * M)
  idx <- 1
  
  for (T_curr in T_values) {
    mean_N <- (T_curr / alpha_true)^beta_true
    log_T <- log(T_curr)
    
    for (m in seq_len(M)) {
      w <- sim_power_law_process(alpha_true, beta_true, T_curr)
      while (length(w) < 2) {
        w <- sim_power_law_process(alpha_true, beta_true, T_curr)
      }
      
      est <- mle_power_law(w, T_curr)
      n <- est$n
      sqrt_n <- sqrt(n)
      
      z_beta <- sqrt_n * (est$beta_hat / beta_true - 1)
      y_beta <- sqrt_n * (beta_true / est$beta_hat - 1)
      
      z_alpha <- (sqrt_n / log_T) * (est$alpha_hat / alpha_true - 1)
      z_alpha_log <- (sqrt_n / log_T) * log(est$alpha_hat / alpha_true)
      z_alpha_naive <- sqrt_n * (est$alpha_hat / alpha_true - 1)
      
      records[[idx]] <- list(
        T_val = T_curr,
        mean_N = round(mean_N, 1),
        T_label = paste0("T = ", T_curr, " (E[N] ~ ", round(mean_N), ")"),
        n = n,
        beta_hat = est$beta_hat,
        beta_tilde = est$beta_tilde,
        alpha_hat = est$alpha_hat,
        z_beta = z_beta,
        y_beta = y_beta,
        z_alpha = z_alpha,
        z_alpha_log = z_alpha_log,
        z_alpha_naive = z_alpha_naive
      )
      idx <- idx + 1
    }
  }
  
  df_res <- do.call(rbind.data.frame, records)
  df_res$T_label <- factor(df_res$T_label, levels = unique(df_res$T_label))
  return(df_res)
}

# ------------------------------------------------------------------------------
# 4. Illustration 2 : Normalité Asymptotique (Proposition 3.11 & 3.12)
# ------------------------------------------------------------------------------

#' Graphique de normalité asymptotique pour beta_hat
plot_normality_beta <- function(mc_df) {
  theme_clean <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey30", size = 10),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  p_hist <- ggplot(mc_df, aes(x = z_beta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 35, fill = "#3182bd", color = "white", alpha = 0.65) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red3", linewidth = 0.9) +
    facet_wrap(~ T_label, scales = "free_y", nrow = 2) +
    xlim(-4, 4) +
    labs(
      title = expression("Densité empirique de " ~ Z[beta] == sqrt(N(t)) * (frac(hat(beta), beta) - 1)),
      x = expression(Z[beta]),
      y = "Densité"
    ) +
    theme_clean
  
  p_qq <- ggplot(mc_df, aes(sample = z_beta)) +
    stat_qq(color = "#3182bd", alpha = 0.35, size = 0.8) +
    stat_qq_line(color = "red3", linewidth = 0.8) +
    facet_wrap(~ T_label, nrow = 2) +
    labs(
      title = expression("QQ-plot gaussien de " ~ Z[beta]),
      x = "Quantiles théoriques N(0, 1)",
      y = "Quantiles empiriques"
    ) +
    theme_clean
  
  p_combined <- p_hist / p_qq
  return(p_combined)
}

#' Graphique de normalité asymptotique pour alpha_hat et effet du facteur log(T)
plot_normality_alpha <- function(mc_df) {
  theme_clean <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey30", size = 10),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  p_hist_log <- ggplot(mc_df, aes(x = z_alpha)) +
    geom_histogram(aes(y = after_stat(density)), bins = 35, fill = "#31a354", color = "white", alpha = 0.65) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red3", linewidth = 0.9) +
    facet_wrap(~ T_label, scales = "free_y", nrow = 2) +
    xlim(-4, 4) +
    labs(
      title = expression("Densité empirique de " ~ Z[alpha] == frac(sqrt(N(t)), log(t)) * (frac(hat(alpha), alpha) - 1)),
      x = expression(Z[alpha]),
      y = "Densité"
    ) +
    theme_clean
  
  var_stats <- aggregate(
    cbind(var_correct = mc_df$z_alpha, var_naive = mc_df$z_alpha_naive),
    by = list(T_val = mc_df$T_val, mean_N = mc_df$mean_N),
    FUN = var
  )
  
  p_var <- ggplot(var_stats, aes(x = T_val)) +
    geom_line(aes(y = var_correct, color = "Avec diviseur log(T)"), linewidth = 1) +
    geom_point(aes(y = var_correct, color = "Avec diviseur log(T)"), size = 3) +
    geom_line(aes(y = var_naive, color = "Sans diviseur log(T)"), linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = var_naive, color = "Sans diviseur log(T)"), size = 3) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black", linewidth = 0.7) +
    scale_color_manual(values = c("Avec diviseur log(T)" = "#31a354",
                                  "Sans diviseur log(T)" = "firebrick")) +
    labs(
      title = "Variance empirique selon l'horizon d'observation T",
      x = "Horizon d'observation T",
      y = "Variance empirique",
      color = "Normalisation"
    ) +
    theme_clean +
    theme(legend.position = "bottom")
  
  p_combined <- p_hist_log / p_var
  return(p_combined)
}

# 5. figure 3 : Dégénérescence asymptotique conjointe (Prop 3.11 et 3.12)

plot_joint_degeneracy <- function(mc_df) {
  theme_clean <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(color = "grey30", size = 9),
      legend.position = "none"
    )
  
  max_T_label <- levels(mc_df$T_label)[length(levels(mc_df$T_label))]
  df_sub <- mc_df[mc_df$T_label == max_T_label, ]
  
  cor_neg <- cor(df_sub$z_alpha_log, df_sub$y_beta)
  cor_pos <- cor(df_sub$z_alpha, df_sub$z_beta)
  
  # 1. Cas Prop 3.11 3.a : (X, -X)
  p_neg <- ggplot(df_sub, aes(x = z_alpha_log, y = y_beta)) +
    geom_point(alpha = 0.35, color = "#756bb1", size = 1.2) +
    geom_abline(intercept = 0, slope = -1, color = "black", linewidth = 1, linetype = "solid") +
    xlim(-3.5, 3.5) + ylim(-3.5, 3.5) +
    labs(
      title = expression("Forme (3.a) : " ~ paste("(", Z[alpha]^log, ", ", Y[beta], ")")),
      x = expression(Z[alpha]^log == frac(sqrt(N), log(T)) * log(hat(alpha)/alpha)),
      y = expression(Y[beta] == sqrt(N) * (frac(beta, hat(beta)) - 1))
    ) +
    theme_clean
  
  # 2. Cas Prop 3.11 3.b : (X, X)
  p_pos <- ggplot(df_sub, aes(x = z_alpha, y = z_beta)) +
    geom_point(alpha = 0.35, color = "#e6550d", size = 1.2) +
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1, linetype = "solid") +
    xlim(-3.5, 3.5) + ylim(-3.5, 3.5) +
    labs(
      title = expression("Forme (3.b) : " ~ paste("(", Z[alpha], ", ", Z[beta], ")")),
      x = expression(Z[alpha] == frac(sqrt(N), log(T)) * (frac(hat(alpha), alpha) - 1)),
      y = expression(Z[beta] == sqrt(N) * (frac(hat(beta), beta) - 1))
    ) +
    theme_clean
  
  p_combined <- p_neg + p_pos +
    plot_annotation(
      title = "Distribution limite conjointe du couple d'estimateurs",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13)
      )
    )
  
  return(p_combined)
}

# 6. fig 4 : Biais à Taille Finie de beta_hat vs beta_tilde

plot_finite_sample_bias <- function(mc_df, beta_true = 1.5) {
  theme_clean <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
  
  bias_summary <- aggregate(
    cbind(
      mean_beta_hat = mc_df$beta_hat,
      mean_beta_tilde = mc_df$beta_tilde,
      mean_N = mc_df$n
    ),
    by = list(T_val = mc_df$T_val, T_label = mc_df$T_label),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  bias_summary$bias_hat <- (bias_summary$mean_beta_hat - beta_true) / beta_true
  bias_summary$bias_tilde <- (bias_summary$mean_beta_tilde - beta_true) / beta_true
  bias_summary$theo_bias_hat <- 1 / (bias_summary$mean_N - 1)
  
  p_bias <- ggplot(bias_summary, aes(x = mean_N)) +
    geom_line(aes(y = bias_hat * 100, color = "EMV brut hat(beta)"), linewidth = 1) +
    geom_point(aes(y = bias_hat * 100, color = "EMV brut hat(beta)"), size = 3) +
    geom_line(aes(y = theo_bias_hat * 100, color = "Courbe théorique 1 / (n - 1)"), linetype = "dotted", linewidth = 1) +
    geom_line(aes(y = bias_tilde * 100, color = "Estimateur débiaisé tilde(beta)"), linewidth = 1) +
    geom_point(aes(y = bias_tilde * 100, color = "Estimateur débiaisé tilde(beta)"), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = c(
      "EMV brut hat(beta)" = "firebrick",
      "Courbe théorique 1 / (n - 1)" = "black",
      "Estimateur débiaisé tilde(beta)" = "#31a354"
    )) +
    labs(
      title = expression("Biais relatif moyen en fonction de " ~ E(N(T))),
      x = expression("Nombre moyen d'événements observés " ~ E(N(T))),
      y = "Biais relatif moyen (%)",
      color = "Estimateur"
    ) +
    theme_clean
  
  return(p_bias)
}

run_all_illustrations <- function(save_figs = TRUE, output_dir = "figures") {
  cat("Illustrations des Convergences de l'EMV\n")
  
  if (save_figs && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  alpha_true <- 1.0
  beta_true <- 1.5
  
  p_as <- plot_as_convergence(alpha_true = alpha_true, beta_true = beta_true, T_max = 200, n_traj = 3)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig01_convergence_ps.png"), p_as, width = 8.5, height = 7, dpi = 200)
  }
  
  T_vals <- c(10, 25, 60, 150)
  mc_data <- monte_carlo_mle(alpha_true = alpha_true, beta_true = beta_true, T_values = T_vals, M = 2000, seed = 2026)
  
  p_beta_norm <- plot_normality_beta(mc_data)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig02_normalite_beta.png"), p_beta_norm, width = 8.5, height = 7, dpi = 200)
  }
  
  p_alpha_norm <- plot_normality_alpha(mc_data)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig03_normalite_alpha_vitesse_log.png"), p_alpha_norm, width = 8.5, height = 7, dpi = 200)
  }
  
  p_joint <- plot_joint_degeneracy(mc_data)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig04_degenerescence_conjointe.png"), p_joint, width = 9.5, height = 5, dpi = 200)
  }
  
  p_bias <- plot_finite_sample_bias(mc_data, beta_true = beta_true)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig05_biais_beta_taille_finie.png"), p_bias, width = 8, height = 4.5, dpi = 200)
  }
  
  if (save_figs) {
    cat("Figures sauvegardées dans :", normalizePath(output_dir), "\n")
  }
  
  invisible(list(
    mc_data = mc_data,
    plots = list(
      p_as = p_as,
      p_beta_norm = p_beta_norm,
      p_alpha_norm = p_alpha_norm,
      p_joint = p_joint,
      p_bias = p_bias
    )
  ))
}

# pour CLI
if (sys.nframe() == 0) {
  run_all_illustrations(save_figs = TRUE)
}