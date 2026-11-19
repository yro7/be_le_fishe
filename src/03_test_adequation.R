# Partie tests d'adéquation

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Glue config
if (!exists("CONFIG") || !exists("mle_power_law")) {
  src_dir <- if (file.exists("src/power_law.R")) {
    "src"
  } else if (file.exists("power_law.R")) {
    "."
  } else {
    file.path(rprojroot::find_rstudio_root_file(), "src")
  }
  if (!exists("CONFIG")) source(file.path(src_dir, "config.R"))
  if (!exists("mle_power_law")) source(file.path(src_dir, "power_law.R"))
}

# 2. Simulation de processus différents (qui ne suivent pas une loi de weibull)

#' Processus de Poisson à intensité log-linéaire (Cox & Lewis, 1966)
#'
#' lambda(t) = c * exp(delta * t / T) sur [0, T], c calibré pour E[N(T)] = mean_N.
#' La forme de l'intensité sur [0, T] ne dépend que de delta.
sim_loglinear_process <- function(delta, mean_N, T_end) {
  n <- rpois(1, mean_N)
  if (n == 0) return(numeric(0))
  v <- runif(n)
  u <- if (abs(delta) < 1e-12) v else log1p(v * expm1(delta)) / delta
  sort(T_end * u)
}

#' Processus de Poisson à rupture d'intensité (ex. révision majeure)
#'
#' Intensité c sur [0, tau * T) puis rho * c sur [tau * T, T].
sim_changepoint_process <- function(rho, mean_N, T_end, tau = 0.5) {
  n <- rpois(1, mean_N)
  if (n == 0) return(numeric(0))
  p_first <- tau / (tau + rho * (1 - tau))
  in_first <- runif(n) < p_first
  u <- ifelse(in_first, runif(n, 0, tau), runif(n, tau, 1))
  sort(T_end * u)
}

#' Processus de renouvellement à inter-arrivées de Weibull (à ne pas confondre avec le processus de Weibull)
#'
#' Inter-arrivées i.i.d. Weibull(shape), d'espérance 1 => E[N(T)] ~ T.
sim_weibull_renewal <- function(shape, T_end) {
  scale <- 1 / gamma(1 + 1 / shape)
  w <- numeric(0)
  t_curr <- 0
  while (t_curr <= T_end) {
    x <- cumsum(rweibull(ceiling(1.5 * T_end) + 10, shape = shape, scale = scale))
    w <- c(w, t_curr + x)
    t_curr <- t_curr + x[length(x)]
  }
  w[w <= T_end]
}

# ------------------------------------------------------------------------------
# 3. Études Monte-Carlo : niveau et puissance
# ------------------------------------------------------------------------------

#' Fonction générique : taux de rejet d'un générateur de trajectoires
#'
#' @param sim_fun Fonction sans argument renvoyant un vecteur d'instants sur [0, T]
#' @param T_end   Horizon d'observation
#' @param M       Nombre de réplications
#' @param B       Taille de la table de calibration
#' @param risk    Risque de première espèce du test
#' @param min_n   Nombre minimal d'événements (rejet des trajectoires trop courtes)
#' @return Liste : résumé (taux de rejet, n moyen) et p-valeurs brutes
rejection_rate <- function(sim_fun, T_end, M, B = CONFIG$gof$B, risk = 1 - CONFIG$conf_level,
                           min_n = CONFIG$gof$min_n) {
  p_cal <- numeric(M)
  p_naive <- numeric(M)
  n_obs <- integer(M)
  for (r in seq_len(M)) {
    w <- sim_fun()
    while (length(w) < min_n) w <- sim_fun()
    res <- gof_power_law(w, T_end = T_end, B = B, risk = risk)
    p_cal[r] <- res$p_value
    p_naive[r] <- res$p_value_naive
    n_obs[r] <- res$m
  }
  list(
    summary = data.frame(
      mean_n = mean(n_obs),
      reject_cal = mean(p_cal <= risk),
      reject_naive = mean(p_naive <= risk)
    ),
    p_cal = p_cal,
    p_naive = p_naive
  )
}

#' Niveau empirique sous H0 (processus de Weibull) pour plusieurs beta
monte_carlo_gof_level <- function(alpha_true = CONFIG$model$alpha_true,
                                  beta_values = CONFIG$gof$level_beta_values,
                                  mean_N_values = CONFIG$gof$mean_N_values,
                                  M = CONFIG$gof$M_level, B = CONFIG$gof$B,
                                  risk = 1 - CONFIG$conf_level,
                                  seed = CONFIG$gof$seed_level) {
  set.seed(seed)
  summary_records <- list()
  pval_records <- list()
  for (beta_true in beta_values) {
    for (mean_N in mean_N_values) {
      # Horizon T tel que E[N(T)] = (T / alpha)^beta = mean_N
      T_curr <- alpha_true * mean_N^(1 / beta_true)
      rr <- rejection_rate(
        function() sim_power_law_process(alpha_true, beta_true, T_curr),
        T_end = T_curr, M = M, B = B, risk = risk
      )
      summary_records[[length(summary_records) + 1]] <- cbind(
        data.frame(beta = beta_true, mean_N = mean_N, T_val = T_curr), rr$summary
      )
      pval_records[[length(pval_records) + 1]] <- data.frame(
        beta = beta_true, mean_N = mean_N,
        method = rep(c("Calibré (Lilliefors MC)", "Naïf (loi de Kolmogorov)"), each = M),
        p_value = c(rr$p_cal, rr$p_naive)
      )
    }
  }
  list(
    summary = do.call(rbind, summary_records),
    pvalues = do.call(rbind, pval_records),
    M = M, risk = risk
  )
}

#' Puissance contre trois familles d'alternatives
monte_carlo_gof_power <- function(mean_N_values = CONFIG$gof$mean_N_values,
                                  delta_values = CONFIG$gof$delta_values,
                                  rho_values = CONFIG$gof$rho_values,
                                  shape_values = CONFIG$gof$shape_values,
                                  M = CONFIG$gof$M_power, B = CONFIG$gof$B,
                                  risk = 1 - CONFIG$conf_level,
                                  seed = CONFIG$gof$seed_power) {
  set.seed(seed)
  T_fix <- CONFIG$gof$T_alt  # sans effet pour les alternatives normalisées (log-linéaire, rupture)
  records <- list()
  add <- function(family, param_label, mean_N, rr) {
    records[[length(records) + 1]] <<- cbind(
      data.frame(family = family, param = param_label, mean_N = mean_N), rr$summary
    )
  }
  for (mean_N in mean_N_values) {
    for (delta in delta_values) {
      rr <- rejection_rate(function() sim_loglinear_process(delta, mean_N, T_fix),
                           T_end = T_fix, M = M, B = B, risk = risk)
      add("Intensité log-linéaire", paste0("delta = ", delta), mean_N, rr)
    }
    for (rho in rho_values) {
      rr <- rejection_rate(function() sim_changepoint_process(rho, mean_N, T_fix),
                           T_end = T_fix, M = M, B = B, risk = risk)
      add("Rupture d'intensité en T/2", paste0("rho = ", rho), mean_N, rr)
    }
    for (k in shape_values) {
      # Inter-arrivées d'espérance 1 : T = mean_N donne E[N(T)] ~ mean_N
      rr <- rejection_rate(function() sim_weibull_renewal(k, mean_N),
                           T_end = mean_N, M = M, B = B, risk = risk)
      add("Renouvellement de Weibull", paste0("k = ", k), mean_N, rr)
    }
  }
  df <- do.call(rbind, records)
  df$family <- factor(df$family, levels = c("Intensité log-linéaire",
                                            "Rupture d'intensité en T/2",
                                            "Renouvellement de Weibull"))
  list(summary = df, M = M, risk = risk)
}

# ------------------------------------------------------------------------------
# 4. Visualisations Graphiques (sans interprétation textuelle)
# ------------------------------------------------------------------------------

theme_gof <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey30", size = 10),
      legend.position = "bottom"
    )
}

#' Niveau empirique (calibré vs naïf) et histogrammes des p-valeurs sous H0
plot_gof_level <- function(level_res, mean_N_hist = CONFIG$gof$mean_N_hist) {
  df <- level_res$summary
  risk <- level_res$risk
  M <- level_res$M
  mc_band <- 1.96 * sqrt(risk * (1 - risk) / M)

  df_long <- rbind(
    data.frame(df[, c("beta", "mean_n")], method = "Calibré (Lilliefors MC)", reject = df$reject_cal),
    data.frame(df[, c("beta", "mean_n")], method = "Naïf (loi de Kolmogorov)", reject = df$reject_naive)
  )
  df_long$beta_lab <- factor(paste0("beta = ", df_long$beta))

  p_level <- ggplot(df_long, aes(x = mean_n, y = reject * 100, color = method,
                                 linetype = beta_lab, shape = beta_lab)) +
    annotate("rect", xmin = min(df$mean_n), xmax = max(df$mean_n),
             ymin = (risk - mc_band) * 100, ymax = (risk + mc_band) * 100,
             fill = "grey80", alpha = 0.5) +
    geom_hline(yintercept = risk * 100, linetype = "dotted", linewidth = 0.8) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.3) +
    scale_x_log10() +
    scale_color_manual(values = c("Calibré (Lilliefors MC)" = "#31a354",
                                  "Naïf (loi de Kolmogorov)" = "firebrick")) +
    labs(
      title = paste0("Niveau empirique du test sous H0 (nominal ", 100 * risk, "%)"),
      subtitle = paste0("M = ", M, " réplications par point ; bande grise : fluctuation Monte-Carlo à 95%"),
      x = expression("Nombre moyen d'événements observés " ~ bar(n) ~ " (échelle log)"),
      y = "Taux de rejet (%)",
      color = "Calibration", linetype = expression(beta), shape = expression(beta)
    ) +
    theme_gof()

  df_p <- level_res$pvalues[level_res$pvalues$mean_N == mean_N_hist, ]
  p_hist <- ggplot(df_p, aes(x = p_value, fill = method)) +
    geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, by = 0.05),
                   color = "white", boundary = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~ method) +
    scale_fill_manual(values = c("Calibré (Lilliefors MC)" = "#31a354",
                                 "Naïf (loi de Kolmogorov)" = "firebrick"), guide = "none") +
    labs(
      title = paste0("Distribution des p-valeurs sous H0 (E[N(T)] = ", mean_N_hist,
                     ", tous beta confondus)"),
      x = "p-valeur", y = "Densité"
    ) +
    theme_gof()

  p_level / p_hist + plot_layout(heights = c(1.3, 1))
}

#' Puissance du test calibré contre les trois familles d'alternatives
plot_gof_power <- function(power_res) {
  df <- power_res$summary
  risk <- power_res$risk
  ggplot(df, aes(x = mean_n, y = reject_cal * 100, color = param)) +
    geom_hline(yintercept = risk * 100, linetype = "dotted", linewidth = 0.8) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_wrap(~ family, ncol = 3) +
    scale_x_log10() +
    ylim(0, 100) +
    labs(
      title = paste0("Puissance empirique du test calibré (niveau ", 100 * risk, "%)"),
      subtitle = paste0("M = ", power_res$M, " réplications par point"),
      x = expression("Nombre moyen d'événements observés " ~ bar(n) ~ " (échelle log)"),
      y = "Taux de rejet (%)",
      color = "Paramètre"
    ) +
    theme_gof() +
    theme(legend.position = "right")
}

#' Diagramme des fractiles exponentiel (Møller, 1976, Fig. 1-4) avec bande KS
#'
#' Points : (quantile Exp(1) théorique, beta_hat * Y_(i)) ; sous H0 ils sont
#' alignés sur la première bissectrice. La bande correspond à la valeur critique
#' calibrée du test au risque risk (inversion de la bande de Kolmogorov).
#'
#' @param w     Instants de sauts
#' @param T_end Horizon (NULL = troncature par le n-ième défaut)
#' @param title Titre du graphique
plot_gof_fractile <- function(w, T_end = NULL, title = NULL, B = CONFIG$gof$B,
                              risk = 1 - CONFIG$conf_level) {
  y <- gof_transform(w, T_end)
  res <- gof_power_law(w, T_end = T_end, B = B, risk = risk)
  m <- res$m
  i <- seq_len(m)
  p_i <- (i - 0.5) / m
  d <- res$crit_value
  q_exp <- function(p) -log(1 - pmin(pmax(p, 0), 1 - 1e-6))
  df <- data.frame(
    theo = q_exp(p_i),
    obs = res$beta_hat * y,
    # Bande KS : |H - F| <= d  <=>  F^{-1}(H - d) <= x <= F^{-1}(H + d)
    low = q_exp(i / m - d),
    up = q_exp((i - 1) / m + d)
  )
  if (is.null(title)) title <- "Diagramme des fractiles exponentiel"
  ggplot(df, aes(x = theo)) +
    geom_ribbon(aes(ymin = low, ymax = up), fill = "#9ecae1", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(aes(y = obs), size = 1.8) +
    labs(
      title = title,
      subtitle = sprintf("n = %d ; beta_hat = %.3f ; D = %.3f ; p-valeur = %.3f",
                         m, res$beta_hat, res$D, res$p_value),
      x = "Quantile théorique Exp(1)",
      y = expression(hat(beta) %*% Y[(i)])
    ) +
    coord_cartesian(ylim = c(0, 1.15 * max(df$obs, df$theo))) +
    theme_gof()
}

#' Illustration : un échantillon sous H0 et deux sous alternatives
plot_gof_examples <- function(ex = CONFIG$gof$examples, T_alt = CONFIG$gof$T_alt,
                              alpha = CONFIG$model$alpha_true) {
  set.seed(ex$seed)
  T_h0 <- alpha * ex$mean_N^(1 / ex$beta)  # E[N(T)] = (T / alpha)^beta = mean_N
  w_h0 <- sim_power_law_process(alpha = alpha, beta = ex$beta, T_end = T_h0)
  w_ll <- sim_loglinear_process(delta = ex$delta, mean_N = ex$mean_N, T_end = T_alt)
  w_cp <- sim_changepoint_process(rho = ex$rho, mean_N = ex$mean_N, T_end = T_alt)
  p1 <- plot_gof_fractile(w_h0, T_end = T_h0, title = paste0("H0 : processus de Weibull (beta = ", ex$beta, ")"))
  p2 <- plot_gof_fractile(w_ll, T_end = T_alt, title = paste0("Log-linéaire (delta = ", ex$delta, ")"))
  p3 <- plot_gof_fractile(w_cp, T_end = T_alt, title = paste0("Rupture (rho = ", ex$rho, ")"))
  p1 + p2 + p3
}

# ------------------------------------------------------------------------------
# 5. Bloc Principal d'Exécution et Sauvegarde des Figures
# ------------------------------------------------------------------------------

run_all_gof_validation <- function(save_figs = TRUE, output_dir = "figures",
                                   M_level = CONFIG$gof$M_level,
                                   M_power = CONFIG$gof$M_power,
                                   B = CONFIG$gof$B) {
  cat("====================================================================\n")
  cat("BE Poisson : Validation pratique du test d'adéquation (Møller, 1976)\n")
  cat("====================================================================\n\n")

  if (save_figs && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cat("1. Niveau empirique sous H0 (M =", M_level, "réplications par point)...\n")
  level_res <- monte_carlo_gof_level(M = M_level, B = B)
  p_level <- plot_gof_level(level_res)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig09_niveau_test_adequation.png"), p_level, width = 9, height = 8, dpi = 200)
  }

  cat("2. Puissance contre les alternatives (M =", M_power, "réplications par point)...\n")
  power_res <- monte_carlo_gof_power(M = M_power, B = B)
  p_power <- plot_gof_power(power_res)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig10_puissance_test_adequation.png"), p_power, width = 11, height = 4.5, dpi = 200)
  }

  cat("3. Diagrammes des fractiles illustratifs...\n")
  p_examples <- plot_gof_examples()
  if (save_figs) {
    ggsave(file.path(output_dir, "fig11_fractiles_exemples.png"), p_examples, width = 12, height = 4.5, dpi = 200)
  }

  cat("\nValidation numérique du test d'adéquation terminée avec succès !\n")
  if (save_figs) {
    cat("Figures sauvegardées dans :", normalizePath(output_dir), "\n")
  }

  invisible(list(
    level_res = level_res,
    power_res = power_res,
    plots = list(
      p_level = p_level,
      p_power = p_power,
      p_examples = p_examples
    )
  ))
}

# Exécution autonome si appelé en CLI ou source
if (sys.nframe() == 0) {
  run_all_gof_validation(save_figs = TRUE)
}
