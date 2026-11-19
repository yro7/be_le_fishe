# ==============================================================================
# Fichier : src/config.R
# Projet  : BE Processus Ponctuels & Fiabilité — INSA Toulouse, 5ModIA
#
# Objet   : Paramètres de TOUTES les études numériques du projet, au même endroit
#           pour pouvoir les relire et les modifier sans fouiller le code.
#           Les formules (EMV, IC, test) sont dans src/power_law.R ; les études
#           (01, 02, 03) lisent leurs valeurs par défaut dans CONFIG.
#
# Conventions :
#   - conf_level est le niveau de confiance des IC ; le risque des tests vaut
#     toujours 1 - conf_level (jamais de second paramètre "level" ambigu).
#   - Chaque graine est nommée : changer une étude ne modifie pas les autres.
#   - Modifier ce fichier invalide le cache knitr des chunks de simulation
#     (cache.extra dans rapport.Rmd).
# ==============================================================================

CONFIG <- list(

  # --- Vraies valeurs du processus de loi de puissance simulé ----------------
  model = list(
    alpha_true = 1.0,   # échelle (alpha = 1 : voir la remarque de monte_carlo_mle)
    beta_true  = 1.5    # forme (beta > 1 : usure)
  ),

  # --- Niveau commun des IC ; risque des tests = 1 - conf_level ---------------
  conf_level = 0.95,

  # --- 01 : convergences de l'EMV (Prop. 3.7, 3.11, 3.12 de Cocozza-Thivent) --
  mle = list(
    # Convergence presque sûre le long de trajectoires
    as_T_max      = 200,   # horizon maximal d'une trajectoire
    as_n_traj     = 3,     # nombre de trajectoires tracées
    as_n_points   = 400,   # points de la grille en t
    as_min_events = 10,    # nombre minimal d'événements sur [0, T_max]
    as_seed_base  = 100,   # graine de la trajectoire i : as_seed_base + i
    # Normalité asymptotique et biais (Monte-Carlo)
    T_values = c(10, 25, 60, 150),
    M        = 2000,       # réplications par horizon
    min_n    = 2,          # trajectoires avec moins d'événements : retirées
    seed     = 2026
  ),

  # --- 02 : couverture des intervalles de confiance ---------------------------
  ci = list(
    T_values   = c(3, 5, 10, 20, 40, 80, 150),
    M          = 2000,     # réplications par horizon
    min_n      = 2,
    seed       = 2026,
    n_kept     = 60,       # réplications conservées par horizon pour les segments
    T_segments = 5,        # horizon de la figure des segments
    n_segments = 45        # segments affichés
  ),

  # --- 03 : niveau et puissance du test d'adéquation (Møller, 1976) -----------
  gof = list(
    B              = 4000, # taille de la table de Lilliefors (loi nulle de D)
    null_seed_base = 1969, # graine de la table de taille m : null_seed_base + m
    min_m          = 3,    # taille minimale de l'échantillon exponentiel testé
    min_n          = 3,    # trajectoires avec moins d'événements : retirées
    mean_N_values  = c(5, 10, 20, 40, 80, 150),
    # Niveau sous H0
    level_beta_values = c(0.5, 1.0, 2.0),
    M_level           = 2000,
    seed_level        = 2026,
    mean_N_hist       = 20,  # E[N(T)] de l'histogramme des p-valeurs
    # Puissance contre les alternatives hors modèle
    M_power      = 1000,
    seed_power   = 2027,
    T_alt        = 100,      # horizon des alternatives normalisées
    delta_values = c(-3, 2, 4),       # intensité log-linéaire
    rho_values   = c(0.25, 4),        # rupture d'intensité en T/2
    shape_values = c(0.5, 2, 3),      # renouvellement de Weibull
    # Diagrammes des fractiles d'exemple
    examples = list(seed = 7, beta = 1.5, mean_N = 30, delta = -3, rho = 0.25)
  )
)
