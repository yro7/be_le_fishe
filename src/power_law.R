# ==============================================================================
# Fichier : src/power_law.R
# Projet  : BE Processus Ponctuels & Fiabilité — INSA Toulouse, 5ModIA
#
# Objet   : Fonctions du modèle "processus de loi de puissance" (iPP d'intensité
#           lambda(s) = (beta / alpha) * (s / alpha)^(beta - 1)) observé sur [0, T] :
#           simulation, EMV, intervalles de confiance et test d'adéquation.
#           Chaque fonction renvoie à son résultat dans la partie théorique
#           (src/rapport/00_partie_theorique.tex).
#           Ce fichier ne contient que des formules : aucune étude Monte-Carlo.
#           Les valeurs par défaut des arguments sont lues dans CONFIG
#           (src/config.R), qui doit être chargé avant tout appel.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Simulation et EMV (Remarque "Simulation", Proposition "Vraisemblance et EMV")
# ------------------------------------------------------------------------------

# Simule les temps d'arrivée d'un Processus Loi de Puissance sur [0, T]
#
# Intensité : lambda(s) = (beta / alpha) * (s / alpha)^(beta - 1)
# Moyenne cumulée : Lambda(t) = (t / alpha)^beta
#
#' @param alpha Paramètre d'échelle (alpha > 0)
#' @param beta  Paramètre de forme (beta > 0)
#' @param T_end Horizon d'observation (T_end > 0)
#' @return Vecteur numérique des instants de sauts W_1 < ... < W_n sur [0, T_end]
sim_power_law_process <- function(alpha, beta, T_end) {
  if (alpha <= 0 || beta <= 0 || T_end <= 0) {
    stop("alpha, beta et T_end doivent être strictement positifs.")
  }
  
  # 1. Nombre total d'événements N(T) ~ Poisson( (T / alpha)^beta )
  mean_N <- (T_end / alpha)^beta
  n <- rpois(1, lambda = mean_N)
  
  if (n == 0) {
    return(numeric(0))
  }
  
  # 2. Conditionnellement à N(T) = n, (W_i / T)^beta sont des statistiques
  #    d'ordre de variables uniformes U(0, 1) i.i.d.
  #    => W_i = T * U_(i)^(1 / beta)
  u <- sort(runif(n, min = 0, max = 1))
  w <- T_end * (u^(1 / beta))
  return(w)
}

# Calcule les Estimateurs du Maximum de Vraisemblance (EMV) à horizon T
#
# Formules issues de Cocozza-Thivent (1997, section 3.2.3) :
#   1 / beta_hat = log(T) - (1 / n) * sum(log(W_i)) = (1 / n) * sum(log(T / W_i))
#   log(alpha_hat) = log(T) - (1 / beta_hat) * log(n) => alpha_hat = T / (n^(1 / beta_hat))
#   Estimateur débiaisé : beta_tilde = ((n - 1) / n) * beta_hat
#
#' @param w     Vecteur ordonné des instants de sauts observés dans [0, T]
#' @param T_end Horizon d'observation T
#' @return Liste contenant n, beta_hat, beta_tilde, alpha_hat (ou NA si n == 0)
mle_power_law <- function(w, T_end) {
  n <- length(w)
  if (n == 0) {
    return(list(n = 0, beta_hat = NA_real_, beta_tilde = NA_real_, alpha_hat = NA_real_))
  }
  
  sum_log_ratio <- sum(log(T_end / w))
  
  beta_hat <- n / sum_log_ratio
  beta_tilde <- if (n > 1) ((n - 1) / n) * beta_hat else NA_real_
  alpha_hat <- T_end / (n^(1 / beta_hat))
  
  return(list(
    n = n,
    beta_hat = beta_hat,
    beta_tilde = beta_tilde,
    alpha_hat = alpha_hat
  ))
}

# ------------------------------------------------------------------------------
# 2. Intervalles de confiance (Corollaires "IC exacts" et "IC asymptotiques")
# ------------------------------------------------------------------------------

#' Calcule les intervalles de confiance exacts et asymptotiques pour alpha et beta
#'
#' @param n         Nombre d'événements observés N(T)
#' @param beta_hat  EMV de beta
#' @param alpha_hat EMV de alpha
#' @param T_end     Horizon d'observation T
#' @param level     Niveau de confiance nominal (défaut : CONFIG$conf_level)
#' @return Liste des bornes [lower, upper] pour chaque méthode
confint_power_law <- function(n, beta_hat, alpha_hat, T_end, level = CONFIG$conf_level) {
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
  # Basé sur (sqrt(n) / log(T)) * (alpha_hat / alpha - 1) -> N(0, 1), d'où par Slutsky
  # alpha dans alpha_hat * (1 -/+ u * log(T) / sqrt(n)).
  # Attention : [CT97] écrit alpha_hat -/+ u * log(T) / sqrt(n), sans le facteur alpha_hat
  # (borne non homogène à un temps, correcte seulement si alpha = 1).
  scale_factor_lin <- (log_T / sqrt_n) * u_val
  alpha_lin <- c(
    lower = alpha_hat * (1 - scale_factor_lin),
    upper = alpha_hat * (1 + scale_factor_lin)
  )
  
  return(list(
    beta_exact = beta_exact,
    beta_asymp = beta_asymp,
    alpha_log = alpha_log,
    alpha_lin = alpha_lin
  ))
}

# ------------------------------------------------------------------------------
# 3. Test d'adéquation (Møller, 1976 ; Lilliefors, 1969)
# ------------------------------------------------------------------------------

#' Transforme les instants de défaillance en échantillon exponentiel (Théorème 1)
#'
#' @param w     Vecteur des instants de sauts W_1 < ... < W_n
#' @param T_end Horizon d'observation T. Si NULL : troncature par le n-ième
#'              défaut (Møller), on normalise par W_n et on retire ce point.
#' @return Vecteur trié des Y_i, i.i.d. Exp(beta) sous H0
gof_transform <- function(w, T_end = NULL) {
  w <- sort(w)
  if (is.null(T_end)) {
    n <- length(w)
    y <- log(w[n] / w[-n])
  } else {
    y <- log(T_end / w)
  }
  sort(y)
}

#' Statistique de Kolmogorov-Smirnov avec paramètre estimé (Lilliefors)
#'
#' @param y Échantillon (supposé exponentiel)
#' @return Liste : D (statistique) et beta_hat (EMV du taux exponentiel)
ks_exp_stat <- function(y) {
  y <- sort(y)
  m <- length(y)
  beta_hat <- m / sum(y)
  F_hat <- 1 - exp(-beta_hat * y)
  i <- seq_len(m)
  D <- max(i / m - F_hat, F_hat - (i - 1) / m)
  list(D = D, beta_hat = beta_hat)
}

# Cache des lois nulles simulées, indexé par (m, B)
.gof_null_cache <- new.env(parent = emptyenv())

#' Loi nulle de D pour un échantillon exponentiel de taille m (paramètre estimé)
#'
#' La loi est pivotale (indépendante de beta) : on simule sous Exp(1).
#'
#' @param m Taille de l'échantillon exponentiel
#' @param B Nombre de simulations
#' @return Vecteur de B réalisations de D sous H0
gof_null_distribution <- function(m, B = CONFIG$gof$B) {
  key <- paste(m, B, sep = "_")
  if (!is.null(.gof_null_cache[[key]])) {
    return(.gof_null_cache[[key]])
  }
  # Graine dédiée : la table ne dépend pas de l'état du générateur appelant
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(CONFIG$gof$null_seed_base + m)
  y_mat <- matrix(rexp(m * B), nrow = m)
  d_null <- apply(y_mat, 2, function(y) ks_exp_stat(y)$D)
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())
  .gof_null_cache[[key]] <- d_null
  d_null
}

#' Test d'adéquation au Processus Loi de Puissance (Møller, 1976)
#'
#' @param w     Vecteur des instants de sauts observés
#' @param T_end Horizon d'observation T (NULL = troncature par le n-ième défaut)
#' @param B     Nombre de simulations pour la calibration de Lilliefors
#' @param risk  Risque de première espèce du test (défaut : 1 - CONFIG$conf_level)
#' @return Liste : m, beta_hat, D, p-valeur calibrée, p-valeur naïve (loi de
#'         Kolmogorov sans correction), valeur critique calibrée au risque risk
gof_power_law <- function(w, T_end = NULL, B = CONFIG$gof$B, risk = 1 - CONFIG$conf_level) {
  y <- gof_transform(w, T_end)
  m <- length(y)
  if (m < CONFIG$gof$min_m) {
    return(list(m = m, beta_hat = NA_real_, D = NA_real_, p_value = NA_real_,
                p_value_naive = NA_real_, crit_value = NA_real_))
  }
  st <- ks_exp_stat(y)
  d_null <- gof_null_distribution(m, B)

  # p-valeur Monte-Carlo (correction +1 pour garantir un test de niveau exact)
  p_value <- (1 + sum(d_null >= st$D)) / (B + 1)

  # p-valeur naïve : on fait comme si beta était connu (loi de Kolmogorov)
  p_value_naive <- suppressWarnings(
    ks.test(y, "pexp", rate = st$beta_hat, exact = (m < 100))$p.value
  )

  list(
    m = m,
    beta_hat = st$beta_hat,
    D = st$D,
    p_value = p_value,
    p_value_naive = p_value_naive,
    crit_value = unname(quantile(d_null, 1 - risk, type = 1))
  )
}
