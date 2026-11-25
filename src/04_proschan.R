# Partie données réelles : pannes de climatisation des Boeing 720 (Proschan, 1963)
#
# Unité d'analyse : le segment (avion x avant/après révision), horloge remise à 0.
# Horizon : T = T_n, instant de la dernière panne du segment.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Glue config
if (!exists("CONFIG") || !exists("mle_power_law") || !exists("plot_gof_fractile")) {
  src_dir <- if (file.exists("src/power_law.R")) {
    "src"
  } else if (file.exists("power_law.R")) {
    "."
  } else {
    file.path(rprojroot::find_rstudio_root_file(), "src")
  }
  if (!exists("CONFIG")) source(file.path(src_dir, "config.R"))
  if (!exists("mle_power_law")) source(file.path(src_dir, "power_law.R"))
  if (!exists("plot_gof_fractile")) source(file.path(src_dir, "03_test_adequation.R"))
}

# 1. Données

#' Charge les segments (avion x avant/après révision)
#'
#' Proschan omet l'intervalle qui contient la révision : arrival_time (cumul sur toute
#' la vie de l'avion) est donc faux après elle. On recumule interval dans chaque segment.
#'
#' @param file Chemin du CSV
#' @return Liste nommée par segment ("7908", "7908 MO", ... MO = après révision,
#'         comme Proschan) de listes (plane, phase, w)
load_segments <- function(file = CONFIG$proschan$file) {
  df <- read.csv(file)
  df$after_overhaul <- as.logical(df$after_overhaul)
  df <- df[order(df$plane, df$failure_index), ]
  df$segment <- paste0(df$plane, ifelse(df$after_overhaul, " MO", ""))

  lapply(split(df, factor(df$segment, levels = unique(df$segment))), function(d) {
    list(
      plane = d$plane[1],
      phase = if (d$after_overhaul[1]) "après" else "avant",
      w = cumsum(d$interval)
    )
  })
}

# 2. Analyse d'un segment

#' Test d'adéquation puis, si le modèle est retenu, estimation et IC
#'
#' Statut : "non testable" (trop peu de pannes), "rejeté" (p-valeur <= risk) ou "retenu".
#' Les colonnes d'estimation valent NA si le modèle n'est pas retenu.
#'
#' @param w     Instants de pannes du segment
#' @param level Niveau de confiance des IC
#' @param risk  Risque de première espèce du test
#' @return data.frame à une ligne
analyze_segment <- function(w, level = CONFIG$conf_level, risk = 1 - CONFIG$conf_level) {
  T_end <- max(w)
  gof <- gof_power_law(w, risk = risk)  # sans T_end : T = T_n, la dernière panne est retirée
  status <- if (is.na(gof$p_value)) {
    "non testable"
  } else if (gof$p_value <= risk) {
    "rejeté"
  } else {
    "retenu"
  }

  # Estimation seulement si le modèle est retenu
  est <- list(beta_hat = NA_real_, alpha_hat = NA_real_)
  if (status == "retenu") est <- mle_power_law(w, T_end)
  # confint_* renvoient des NA si les estimateurs sont NA
  ci_beta <- confint_power_law(length(w), est$beta_hat, est$alpha_hat, T_end, level)$beta_exact
  ci_alpha <- confint_alpha_exact(length(w), est$beta_hat, est$alpha_hat, level)

  # Position de 1 par rapport à l'IC de beta
  trend <- if (status != "retenu") {
    NA_character_
  } else if (ci_beta["lower"] > 1) {
    "usure"
  } else if (ci_beta["upper"] < 1) {
    "amélioration"
  } else {
    "hPP possible"
  }

  data.frame(
    n = length(w),
    T_n = T_end,
    D = gof$D,
    p_value = gof$p_value,
    status = status,
    beta_hat = est$beta_hat,
    beta_low = unname(ci_beta["lower"]),
    beta_up = unname(ci_beta["upper"]),
    alpha_hat = est$alpha_hat,
    alpha_low = unname(ci_alpha["lower"]),
    alpha_up = unname(ci_alpha["upper"]),
    trend = trend
  )
}

#' Analyse tous les segments : une ligne par segment
analyze_proschan <- function(segments = load_segments(), level = CONFIG$conf_level,
                             risk = 1 - CONFIG$conf_level, seed = CONFIG$proschan$seed) {
  set.seed(seed)
  rows <- lapply(names(segments), function(id) {
    s <- segments[[id]]
    cbind(data.frame(segment = id, plane = s$plane, phase = s$phase),
          analyze_segment(s$w, level = level, risk = risk))
  })
  do.call(rbind, rows)
}

# 3. graphs

phase_colors <- c("avant" = "#2b8cbe", "après" = "#e6550d")

#' beta_hat et IC par segment, avec la référence beta = 1 (hPP)
plot_proschan_beta <- function(res, level = CONFIG$conf_level) {
  res$label <- factor(sprintf("%s (n = %d)", res$segment, res$n),
                      levels = rev(sprintf("%s (n = %d)", res$segment, res$n)))
  kept <- res[res$status == "retenu", ]
  dropped <- res[res$status != "retenu", ]

  ggplot(kept, aes(x = beta_hat, y = label, color = phase)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_segment(aes(x = beta_low, xend = beta_up, yend = label), linewidth = 0.9) +
    geom_point(size = 2.4) +
    geom_text(data = dropped, aes(x = 1, label = status), color = "grey40",
              size = 3.2, fontface = "italic", hjust = -0.1) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    scale_y_discrete(drop = FALSE) +
    scale_color_manual(values = phase_colors) +
    labs(
      title = expression("Estimation de " ~ beta ~ " par segment"),
      subtitle = sprintf("IC « exact » à %g%% avec T = T_n ; MO : après révision majeure", 100 * level),
      x = expression(hat(beta) ~ " (échelle log)"),
      y = NULL,
      color = "Révision"
    ) +
    theme_gof()
}

#' Nombre de pannes observé et moyenne cumulée estimée, pour chaque segment retenu
plot_proschan_counts <- function(segments, res) {
  kept <- res[res$status == "retenu", ]
  make_label <- function(i) sprintf("%s : beta_hat = %.2f", kept$segment[i], kept$beta_hat[i])
  labels <- vapply(seq_len(nrow(kept)), make_label, character(1))

  obs <- do.call(rbind, lapply(seq_len(nrow(kept)), function(i) {
    w <- segments[[kept$segment[i]]]$w
    data.frame(label = labels[i], t = c(0, w), N = 0:length(w))
  }))
  fit <- do.call(rbind, lapply(seq_len(nrow(kept)), function(i) {
    t <- seq(0, kept$T_n[i], length.out = 200)
    data.frame(label = labels[i], t = t, Lambda = (t / kept$alpha_hat[i])^kept$beta_hat[i])
  }))
  obs$label <- factor(obs$label, levels = labels)
  fit$label <- factor(fit$label, levels = labels)

  ggplot() +
    geom_step(data = obs, aes(x = t, y = N, color = "Observé"), linewidth = 0.7) +
    geom_line(data = fit, aes(x = t, y = Lambda, color = "Estimé"), linewidth = 0.8) +
    facet_wrap(~ label, scales = "free", ncol = 4) +
    scale_color_manual(
      values = c("Observé" = "grey25", "Estimé" = "firebrick"),
      labels = c("Observé" = expression(N[t] ~ " observé"),
                 "Estimé" = expression(hat(Lambda)(t) ~ " estimée"))
    ) +
    labs(
      title = "Nombre cumulé de pannes : observé et modèle ajusté",
      subtitle = "Segments retenus ; la courbe passe par (T_n, n) par construction",
      x = "Temps depuis le début du segment (h)", y = "Pannes cumulées", color = NULL
    ) +
    theme_gof()
}

#' Diagrammes des fractiles de quelques segments
plot_proschan_fractiles <- function(segments, res, examples = CONFIG$proschan$examples) {
  plots <- lapply(examples, function(id) {
    r <- res[res$segment == id, ]
    plot_gof_fractile(segments[[id]]$w, title = paste0("Segment ", id, " (", r$status, ")")) +
      labs(subtitle = sprintf("n = %d ; D = %.3f ; p-valeur = %.3f", r$n, r$D, r$p_value))
  })
  wrap_plots(plots, nrow = 1)
}

# 4. main

run_proschan <- function(save_figs = TRUE, output_dir = "figures") {
  cat("BE Poisson : données de Proschan (1963)\n")

  if (save_figs && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  segments <- load_segments()
  res <- analyze_proschan(segments)
  print(res, digits = 3, row.names = FALSE)

  p_beta <- plot_proschan_beta(res)
  p_counts <- plot_proschan_counts(segments, res)
  p_fractiles <- plot_proschan_fractiles(segments, res)
  if (save_figs) {
    ggsave(file.path(output_dir, "fig12_proschan_beta.png"), p_beta, width = 8.5, height = 6, dpi = 200)
    ggsave(file.path(output_dir, "fig13_proschan_ajustement.png"), p_counts, width = 12, height = 9, dpi = 200)
    ggsave(file.path(output_dir, "fig14_proschan_fractiles.png"), p_fractiles, width = 12, height = 4.5, dpi = 200)
    cat("Figures sauvegardées dans :", normalizePath(output_dir), "\n")
  }

  invisible(list(
    segments = segments,
    results = res,
    plots = list(p_beta = p_beta, p_counts = p_counts, p_fractiles = p_fractiles)
  ))
}

# Exécution autonome si appelé en CLI ou source
if (sys.nframe() == 0) {
  run_proschan(save_figs = TRUE)
}
