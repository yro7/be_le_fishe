# ==============================================================================
# Fichier : src/02_intervalles_confiance.R
# Projet  : BE Processus Ponctuels & Fiabilité — INSA Toulouse, 5ModIA
# Auteur  : Marin Décanini
#
# Objet   : Validation pratique des Intervalles de Confiance pour le Processus
#           Loi de Puissance (Power Law / Weibull Process) à horizon T fixé
#           Références : Cocozza-Thivent (1997), Corollaires 3.6 et 3.13
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Chargement des fonctions de simulation et d'estimation depuis 01_convergences_emv.R
if (!exists("sim_power_law_process") || !exists("mle_power_law")) {
  src_file <- if (file.exists("src/01_convergences_emv.R")) {
    "src/01_convergences_emv.R"
  } else if (file.exists("01_convergences_emv.R")) {
    "01_convergences_emv.R"
  } else {
    file.path(rprojroot::find_rstudio_root_file(), "src", "01_convergences_emv.R")
  }
  source(src_file)
}

# ------------------------------------------------------------------------------
# 1. Calcul des Intervalles de Confiance (Cocozza-Thivent, 1997)
# ------------------------------------------------------------------------------

#' Calcule les intervalles de confiance exacts et asymptotiques pour alpha et beta
#'
#' @param n         Nombre d'événements observés N(T)
#' @param beta_hat  EMV de beta
#' @param alpha_hat EMV de alpha
#' @param T_end     Horizon d'observation T
#' @param level     Niveau de confiance nominal (défaut = 0.95)
#' @return Liste des bornes [lower, upper] pour chaque méthode
confint_power_law <- function(n, beta_hat, alpha_hat, T_end, level = 0.95) {
  if (is.na(n) || n < 1 || is.na(beta_hat) || is.na(alpha_hat)) {
    return(list(
      beta_exact = c(lower = NA_real_, upper = NA_real_),
      beta_asymp = c(lower = NA_real_, upper = NA_real_),
      alpha_log   = c(lower = NA_real_, upper = NA_real_),
      alpha_lin   = c(lower = NA_real_, upper = NA_real_)
    ))
  }
  
  alpha_risk <- 1 - level
  u_val <- qnorm(1 - alpha_risk / 2)
  sqrt_n <- sqrt(n)
  log_T <- log(T_end)
  
  # 1. Intervalle de confiance EXACT pour beta (Corollaire 3.6)
  # Conditionnellement à N(T) = n, 2*n*beta / beta_hat ~ Chi2(2*n)
  chi2_lower <- qchisq(alpha_risk / 2, df = 2 * n)
  chi2_upper <- qchisq(1 - alpha_risk / 2, df = 2 * n)
  beta_exact <- c(
    lower = (beta_hat / (2 * n)) * chi2_lower,
    upper = (beta_hat / (2 * n)) * chi2_upper
  )
  
  # 2. Intervalle de confiance ASYMPTOTIQUE pour beta (Corollaire 3.13-2c)
  # Basé sur la normalité de sqrt(n) * (beta_hat / beta - 1)
  beta_asymp <- c(
    lower = max(0, beta_hat * (1 - u_val / sqrt_n)),
    upper = beta_hat * (1 + u_val / sqrt_n)
  )
  
  # 3. Intervalle de confiance ASYMPTOTIQUE pour alpha : forme log (Corollaire 3.13-1c)
  # Basé sur la normalité de (sqrt(n) / log(T)) * log(alpha_hat / alpha)
  scale_factor_log <- (log_T / sqrt_n) * u_val
  alpha_log <- c(
    lower = alpha_hat * exp(-scale_factor_log),
    upper = alpha_hat * exp(scale_factor_log)
  )
  
  # 4. Intervalle de confiance ASYMPTOTIQUE pour alpha : forme linéaire (Corollaire 3.13-1f)
  # Basé sur (sqrt(n) / log(T)) * (alpha_hat / alpha - 1)
  scale_factor_lin <- (log_T / sqrt_n) * u_val
  alpha_lin <- c(
    lower = alpha_hat - scale_factor_lin,
    upper = alpha_hat + scale_factor_lin
  )
  
  return(list(
    beta_exact = beta_exact,
    beta_asymp = beta_asymp,
    alpha_log = alpha_log,
    alpha_lin = alpha_lin
  ))
}

# ------------------------------------------------------------------------------
# 2. Étude Monte-Carlo du Taux de Couverture et de la Largeur
# ------------------------------------------------------------------------------

#' Évalue le taux de couverture et la largeur moyenne sur une grille d'horizons T
#'
#' @param alpha_true Valeur théorique d'alpha
#' @param beta_true  Valeur théorique de beta
#' @param T_values   Vecteur d'horizons d'observation
#' @param M          Nombre de réplications Monte-Carlo par horizon
#' @param level      Niveau de confiance nominal
#' @param seed       Graine aléatoire
#' @return Liste contenant le résumé agrégé et les données brutes des réplications
monte_carlo_ci_study <- function(alpha_true = 1.0, beta_true = 1.5,
                                T_values = c(3, 5, 10, 20, 40, 80, 150),
                                M = 2000, level = 0.95, seed = 2026) {
  set.seed(seed)
  
  summary_records <- vector("list", length(T_values))
  sample_records <- vector("list", length(T_values))
  
  for (idx in seq_along(T_values)) {
    T_curr <- T_values[idx]
    mean_N <- (T_curr / alpha_true)^beta_true
    
    covered_beta_exact <- logical(M)
    covered_beta_asymp <- logical(M)
    covered_alpha_log   <- logical(M)
    covered_alpha_lin   <- logical(M)
    
    width_beta_exact <- numeric(M)
    width_beta_asymp <- numeric(M)
    width_alpha_log   <- numeric(M)
    width_alpha_lin   <- numeric(M)
    
    min_lower_alpha_lin <- Inf
    n_negative_lower_alpha_lin <- 0
    
    # Échantillon de quelques réplications pour affichage des segments
    rep_sample_list <- vector("list", min(M, 60))
    
    for (m in seq_len(M)) {
      w <- sim_power_law_process(alpha_true, beta_true, T_curr)
      while (length(w) < 2) {
        w <- sim_power_law_process(alpha_true, beta_true, T_curr)
      }
      
      est <- mle_power_law(w, T_curr)
      ci <- confint_power_law(est$n, est$beta_hat, est$alpha_hat, T_curr, level = level)
      
      # Couverture pour beta
      c_be <- (beta_true >= ci$beta_exact["lower"]) && (beta_true <= ci$beta_exact["upper"])
      c_ba <- (beta_true >= ci$beta_asymp["lower"]) && (beta_true <= ci$beta_asymp["upper"])
      
      # Couverture pour alpha
      c_al <- (alpha_true >= ci$alpha_log["lower"]) && (alpha_true <= ci$alpha_log["upper"])
      c_an <- (alpha_true >= ci$alpha_lin["lower"]) && (alpha_true <= ci$alpha_lin["upper"])
      
      covered_beta_exact[m] <- c_be
      covered_beta_asymp[m] <- c_ba
      covered_alpha_log[m]   <- c_al
      covered_alpha_lin[m]   <- c_an
      
      width_beta_exact[m] <- ci$beta_exact["upper"] - ci$beta_exact["lower"]
      width_beta_asymp[m] <- ci$beta_asymp["upper"] - ci$beta_asymp["lower"]
      width_alpha_log[m]   <- ci$alpha_log["upper"] - ci$alpha_log["lower"]
      width_alpha_lin[m]   <- ci$alpha_lin["upper"] - ci$alpha_lin["lower"]
      
      if (ci$alpha_lin["lower"] < min_lower_alpha_lin) {
        min_lower_alpha_lin <- ci$alpha_lin["lower"]
      }
      if (ci$alpha_lin["lower"] < 0) {
        n_negative_lower_alpha_lin <- n_negative_lower_alpha_lin + 1
      }
      
      if (m <= 60) {
        rep_sample_list[[m]] <- data.frame(
          rep = m,
          T_val = T_curr,
          n = est$n,
          beta_hat = est$beta_hat,
          be_low = ci$beta_exact["lower"],
          be_up = ci$beta_exact["upper"],
          be_cov = c_be,
          ba_low = ci$beta_asymp["lower"],
          ba_up = ci$beta_asymp["upper"],
          ba_cov = c_ba
        )
      }
    }
    
    summary_records[[idx]] <- data.frame(
      T_val = T_curr,
      mean_N = mean_N,
      coverage_beta_exact = mean(covered_beta_exact),
      coverage_beta_asymp = mean(covered_beta_asymp),
      coverage_alpha_log   = mean(covered_alpha_log),
      coverage_alpha_lin   = mean(covered_alpha_lin),
      width_beta_exact = mean(width_beta_exact),
      width_beta_asymp = mean(width_beta_asymp),
      width_alpha_log   = mean(width_alpha_log),
      width_alpha_lin   = mean(width_alpha_lin),
      prop_negative_alpha_lin = n_negative_lower_alpha_lin / M
    )
    
    sample_records[[idx]] <- do.call(rbind, rep_sample_list)
  }
  
  df_summary <- do.call(rbind, summary_records)
  df_samples <- do.call(rbind, sample_records)
  
  return(list(
    summary = df_summary,
    samples = df_samples,
    alpha_true = alpha_true,
    beta_true = beta_true,
    level = level
  ))
}

# ------------------------------------------------------------------------------
# 3. Visualisations Graphiques (sans interprétation textuelle)
# ------------------------------------------------------------------------------

#' Graphique de couverture et de largeur pour beta (Exact Chi-2 vs Asymptotique)
plot_ci_comparison_beta <- function(ci_res) {
  df <- ci_res$summary
  level <- ci_res$level
  
  theme_clean <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
  
  # 1. Taux de couverture empirique
  p_cov <- ggplot(df, aes(x = mean_N)) +
    geom_line(aes(y = coverage_beta_exact * 100, color = "Exact (Chi-deux)"), linewidth = 1) +
    geom_point(aes(y = coverage_beta_exact * 100, color = "Exact (Chi-deux)"), size = 2.5) +
    geom_line(aes(y = coverage_beta_asymp * 100, color = "Asymptotique (Gaussien)"), linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = coverage_beta_asymp * 100, color = "Asymptotique (Gaussien)"), size = 2.5) +
    geom_hline(yintercept = level * 100, linetype = "dotted", color = "black", linewidth = 0.8) +
    scale_color_manual(values = c("Exact (Chi-deux)" = "#31a354",
                                  "Asymptotique (Gaussien)" = "firebrick")) +
    ylim(75, 100) +
    labs(
      title = expression("Taux de couverture empirique pour " ~ beta ~ " (nominal 95%)"),
      x = expression("Nombre moyen d'événements " ~ E(N(T))),
      y = "Couverture empirique (%)",
      color = "Méthode"
    ) +
    theme_clean
  
  # 2. Largeur moyenne des intervalles
  p_width <- ggplot(df, aes(x = mean_N)) +
    geom_line(aes(y = width_beta_exact, color = "Exact (Chi-deux)"), linewidth = 1) +
    geom_point(aes(y = width_beta_exact, color = "Exact (Chi-deux)"), size = 2.5) +
    geom_line(aes(y = width_beta_asymp, color = "Asymptotique (Gaussien)"), linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = width_beta_asymp, color = "Asymptotique (Gaussien)"), size = 2.5) +
    scale_color_manual(values = c("Exact (Chi-deux)" = "#31a354",
                                  "Asymptotique (Gaussien)" = "firebrick")) +
    labs(
      title = expression("Largeur moyenne de l'intervalle de confiance pour " ~ beta),
      x = expression("Nombre moyen d'événements " ~ E(N(T))),
      y = "Largeur moyenne",
      color = "Méthode"
    ) +
    theme_clean
  
  p_combined <- p_cov / p_width + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  return(p_combined)
}

#' Graphique de couverture et de largeur pour alpha (Forme log vs Linéaire)
plot_ci_comparison_alpha <- function(ci_res) {
  df <- ci_res$summary
  level <- ci_res$level
  
  theme_clean <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
  
  # 1. Taux de couverture empirique
  p_cov <- ggplot(df, aes(x = mean_N)) +
    geom_line(aes(y = coverage_alpha_log * 100, color = "Forme log-exponentielle"), linewidth = 1) +
    geom_point(aes(y = coverage_alpha_log * 100, color = "Forme log-exponentielle"), size = 2.5) +
    geom_line(aes(y = coverage_alpha_lin * 100, color = "Forme linéaire"), linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = coverage_alpha_lin * 100, color = "Forme linéaire"), size = 2.5) +
    geom_hline(yintercept = level * 100, linetype = "dotted", color = "black", linewidth = 0.8) +
    scale_color_manual(values = c("Forme log-exponentielle" = "#2b8cbe",
                                  "Forme linéaire" = "#e6550d")) +
    ylim(80, 100) +
    labs(
      title = expression("Taux de couverture empirique pour " ~ alpha ~ " (nominal 95%)"),
      x = expression("Nombre moyen d'événements " ~ E(N(T))),
      y = "Couverture empirique (%)",
      color = "Formulation"
    ) +
    theme_clean
  
  # 2. Largeur moyenne des intervalles
  p_width <- ggplot(df, aes(x = mean_N)) +
    geom_line(aes(y = width_alpha_log, color = "Forme log-exponentielle"), linewidth = 1) +
    geom_point(aes(y = width_alpha_log, color = "Forme log-exponentielle"), size = 2.5) +
    geom_line(aes(y = width_alpha_lin, color = "Forme linéaire"), linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = width_alpha_lin, color = "Forme linéaire"), size = 2.5) +
    scale_color_manual(values = c("Forme log-exponentielle" = "#2b8cbe",
                                  "Forme linéaire" = "#e6550d")) +
    labs(
      title = expression("Largeur moyenne de l'intervalle de confiance pour " ~ alpha),
      x = expression("Nombre moyen d'événements " ~ E(N(T))),
      y = "Largeur moyenne",
      color = "Formulation"
    ) +
    theme_clean
  
  p_combined <- p_cov / p_width + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  return(p_combined)
}

#' Visualisation des segments d'intervalles de confiance sur un sous-échantillon
plot_ci_segments_beta <- function(ci_res, T_target = 5) {
  df_sample <- ci_res$samples[ci_res$samples$T_val == T_target, ]
  if (nrow(df_sample) == 0) {
    df_sample <- ci_res$samples[ci_res$samples$T_val == unique(ci_res$samples$T_val)[1], ]
  }
  
  df_sample <- df_sample[seq_len(min(45, nrow(df_sample))), ]
  beta_true <- ci_res$beta_true
  
  theme_clean <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )
  
  # Segments Exacts
  p_exact <- ggplot(df_sample, aes(x = rep)) +
    geom_linerange(aes(ymin = be_low, ymax = be_up, color = be_cov), linewidth = 0.8) +
    geom_point(aes(y = beta_hat), size = 1.2, color = "grey20") +
    geom_hline(yintercept = beta_true, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "#31a354", "FALSE" = "firebrick"),
                       labels = c("TRUE" = "Couvert", "FALSE" = "Non couvert")) +
    labs(
      title = expression("IC Exacts (Chi-deux) sur 45 réplications (" ~ E(N) ~ " ~ 11)"),
      x = "Indice de réplication",
      y = expression(beta),
      color = "Statut"
    ) +
    theme_clean
  
  # Segments Asymptotiques
  p_asymp <- ggplot(df_sample, aes(x = rep)) +
    geom_linerange(aes(ymin = ba_low, ymax = ba_up, color = ba_cov), linewidth = 0.8) +
    geom_point(aes(y = beta_hat), size = 1.2, color = "grey20") +
    geom_hline(yintercept = beta_true, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "#31a354", "FALSE" = "firebrick"),
                       labels = c("TRUE" = "Couvert", "FALSE" = "Non couvert")) +
    labs(
      title = expression("IC Asymptotiques (Gaussiens) sur 45 réplications (" ~ E(N) ~ " ~ 11)"),
      x = "Indice de réplication",
      y = expression(beta),
      color = "Statut"
    ) +
    theme_clean
  
  p_combined <- p_exact + p_asymp + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  return(p_combined)
}

# ------------------------------------------------------------------------------
# 4. Bloc Principal d'Exécution et Sauvegarde des Figures
# ------------------------------------------------------------------------------

run_all_ci_validation <- function(save_figs = TRUE, output_dir = "figures") {
  cat("====================================================================\n")
  cat("BE Poisson : Validation pratique des IC (Cocozza-Thivent, 1997)\n")
  cat("====================================================================\n\n")
  
  if (save_figs && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  alpha_true <- 1.0
  beta_true <- 1.5
  T_vals <- c(3, 5, 10, 20, 40, 80, 150)
  
  cat("1. Simulation Monte-Carlo des IC (M = 2000 réplications sur 7 horizons T)...\n")
  ci_res <- monte_carlo_ci_study(
    alpha_true = alpha_true,
    beta_true = beta_true,
    T_values = T_vals,
    M = 2000,
    level = 0.95,
    seed = 2026
  )
  
  cat("2. Génération de la comparaison de couverture et largeur pour beta...\n")
  p_beta <- plot_ci_comparison_beta(ci_res)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig06_couverture_largeur_beta.png"), p_beta, width = 8.5, height = 7, dpi = 200)
  }
  
  cat("3. Génération de la comparaison de couverture et largeur pour alpha...\n")
  p_alpha <- plot_ci_comparison_alpha(ci_res)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig07_couverture_largeur_alpha.png"), p_alpha, width = 8.5, height = 7, dpi = 200)
  }
  
  cat("4. Génération de la visualisation des segments d'IC pour beta...\n")
  p_segments <- plot_ci_segments_beta(ci_res, T_target = 5)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig08_segments_ic_beta.png"), p_segments, width = 9.5, height = 4.5, dpi = 200)
  }
  
  cat("\nValidation numérique des IC terminée avec succès !\n")
  if (save_figs) {
    cat("Figures sauvegardées dans :", normalizePath(output_dir), "\n")
  }
  
  invisible(list(
    ci_res = ci_res,
    plots = list(
      p_beta = p_beta,
      p_alpha = p_alpha,
      p_segments = p_segments
    )
  ))
}

# Exécution autonome si appelé en CLI ou source
if (sys.nframe() == 0) {
  run_all_ci_validation(save_figs = TRUE)
}

