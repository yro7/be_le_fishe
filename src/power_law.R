# Centralise les formules des papiers (MLE IC etc). La config et les csts sont dans @config.R.

# 1. Simulation et MLE

# Simule les temps d'arrivée d'un Processus de Weibull sur [0, T]
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
  
  # Changement de temps (Proposition 4.8 du poly) : si S_1 < S_2 < ... sont les temps
  # d'arrivée d'un hPP d'intensité 1, alors les Lambda^{-1}(S_k) sont ceux d'un iPP
  # d'intensité lambda, avec ici Lambda^{-1}(s) = alpha * s^(1 / beta).

  # 1. hPP d'intensité 1 sur [0, Lambda(T)] : S_k = E_1 + ... + E_k, E_j i.i.d. Exp(1).
  #    On tire les E_j par blocs jusqu'à dépasser Lambda(T) (pr éviter un 2eme while), puis on tronque.
  Lambda_T <- (T_end / alpha)^beta
  block_size <- ceiling(Lambda_T + 4 * sqrt(Lambda_T)) + 10  
  s <- numeric(0)
  s_last <- 0
  while (s_last <= Lambda_T) {
    s <- c(s, s_last + cumsum(rexp(block_size, rate = 1)))
    s_last <- s[length(s)]
  }
  s <- s[s <= Lambda_T]  # S_k <= Lambda(T) <=> W_k <= T

  # 2. Retour à l'échelle de temps d'origine : W_k = Lambda^{-1}(S_k)
  w <- alpha * s^(1 / beta)
  return(w)
}

# MLE à Horizon T
#
# Cf Cocozza-Thivent (1997, section 3.2.3) :
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

# 2. IC 

#' Calcule les intervalles de confiance exacts et asymptotiques pour alpha et beta
#'
#' @param n         Nombre d'événements observés N(T)
#' @param beta_hat  MLE de beta
#' @param alpha_hat MLE de alpha
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
  # Conditionnellement à N(T) = n, 2*n*beta / beta_hat ~ χ²(2*n)
  chi2_lower <- qchisq(alpha_risk / 2, df = 2 * n)
  chi2_upper <- qchisq(1 - alpha_risk / 2, df = 2 * n)
  beta_exact <- c(
    lower = (beta_hat / (2 * n)) * chi2_lower,
    upper = (beta_hat / (2 * n)) * chi2_upper
  )
  
  # 2. Intervalle de confiance ASYMPTOTIQUE pour β
  # Basé sur la normalité de sqrt(n) * (beta_hat / beta - 1)
  beta_asymp <- c(
    lower = max(0, beta_hat * (1 - u_val / sqrt_n)),
    upper = beta_hat * (1 + u_val / sqrt_n)
  )
  
  # 3. Intervalle de confiance ASYMPTOTIQUE pour α : forme log
  # Basé sur la normalité de (sqrt(n) / log(T)) * log(alpha_hat / alpha)
  scale_factor_log <- (log_T / sqrt_n) * u_val
  alpha_log <- c(
    lower = alpha_hat * exp(-scale_factor_log),
    upper = alpha_hat * exp(scale_factor_log)
  )
  
  # 4. Intervalle de confiance ASYMPTOTIQUE pour alpha : forme linéaire
  # Basé sur (sqrt(n) / log(T)) * (alpha_hat / alpha - 1) -> N(0, 1), d'où par Slutsky
  # alpha dans alpha_hat * (1 -/+ u * log(T) / sqrt(n)).
  # NOTE Coquille dans papier CT97 : écrit alpha_hat -/+ u * log(T) / sqrt(n), sans le facteur alpha_hat
  # (correct seulement si alpha = 1).
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

# 3. Test d'adéquation

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

#' Test d'adéquation au Processus de Weibull (Møller, 1976)
#'
#' Les Y_i sont testés par KScorrect::LcKS : KS contre Exp(beta) avec beta estimé,
#' p-valeur calibrée par Monte-Carlo (Lilliefors). Sous H0 la loi de D ne dépend
#' que de m : cache = TRUE (fork du package) réutilise D.sim d'un appel à l'autre.
#'
#' @param w     Vecteur des instants de sauts observés
#' @param T_end Horizon d'observation T (NULL = troncature par le n-ième défaut)
#' @param B     Nombre de simulations de LcKS (nreps)
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
  beta_hat <- 1 / mean(y)  # MLE du taux exponentiel (non renvoyé par LcKS)
  lc <- KScorrect::LcKS(y, cdf = "pexp", nreps = B, cache = TRUE)

  # p-valeur naïve : on fait comme on connaissait β 
  p_value_naive <- suppressWarnings(
    ks.test(y, "pexp", rate = beta_hat, exact = (m < 100))$p.value
  )

  list(
    m = m,
    beta_hat = beta_hat,
    D = lc$D.obs,
    p_value = lc$p.value,
    p_value_naive = p_value_naive,
    crit_value = unname(quantile(lc$D.sim, 1 - risk, type = 1))
  )
}
