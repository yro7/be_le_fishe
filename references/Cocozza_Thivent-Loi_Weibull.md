<!-- page 1 -->

Ce résultat fournit une approximation exponentielle de la fiabilité qui est pessimiste. Ce type de résultat est fort utile en pratique (du moins lorsqu'on sait calculer assez facilement le paramètre, ici $\lambda(\infty)$, de l'approximation exponentielle). Nous y reviendrons dans le chapitre 7 (paragraphe 7.3.4) et dans le chapitre 9 (paragraphes 9.3 et 9.4).

Des développements sur d'autres approximations exponentielles de la fiabilité et d'autres notions de vieillissement ont été réalisés par l'école probabiliste russe. Le lecteur intéressé est invité à se reporter au livre de J.L. Bon [17], chapitres VII et XIII.

### 1.2.4 Deux familles de lois classiques en fiabilité

#### La loi gamma

Soit $\alpha$ et $\beta$ deux réels strictement positifs. La **loi gamma** de paramètres $(\alpha, \beta)$ est la loi de densité

$$
f(x) = \frac{1}{\Gamma(\alpha)\beta^\alpha} x^{\alpha-1} e^{-\frac{x}{\beta}}
$$

par rapport à la mesure de Lebesgue sur $\mathbb{R}_+$, la fonction $\Gamma$ étant définie pour $\alpha > 0$ par :

$$
\Gamma(\alpha) = \int_0^{+\infty} x^{\alpha-1} e^{-x} \, dx.
$$

Le paramètre $\alpha$ est le paramètre de forme et $\beta$ le paramètre d'échelle. Nous résumons les principales propriétés de la loi gamma dans la proposition suivante.

**Proposition 1.11 :**

1. (a) Pour $0 < \alpha < 1$, la loi gamma de paramètres $(\alpha, \beta)$ est DFR et son taux de hasard $h$ vérifie :

$$
\lim_{x \to 0} h(x) = +\infty, \quad \lim_{x \to +\infty} h(x) = \frac{1}{\beta}.
$$

(b) Pour $\alpha > 1$, la loi gamma de paramètres $(\alpha, \beta)$ est IFR et son taux de hasard $h$ vérifie :

$$
\lim_{x \to 0} h(x) = 0, \quad \lim_{x \to +\infty} h(x) = \frac{1}{\beta}.
$$

2. Soit $T$ une variable aléatoire de loi gamma de paramètres $(\alpha, \beta)$, alors :

$$
\mathbb{E}(T) = \alpha\beta, \quad \mathrm{var}(T) = \alpha\beta^2
$$

$$
\mathbb{E}(e^{-sT}) = \frac{1}{(1+\beta s)^\alpha}, \quad \mathbb{E}(e^{iuT}) = \frac{1}{(1-i\beta u)^\alpha}.
$$

3. Soit $T_1$ et $T_2$ deux variables aléatoires indépendantes de loi gamma de paramètres respectifs $(\alpha_1, \beta)$ et $(\alpha_2, \beta)$. Alors $T_1 + T_2$ est de loi gamma de paramètres $(\alpha_1 + \alpha_2, \beta)$.

La famille des lois gamma contient des lois connues. La loi gamma de paramètres

- $\alpha = 1$ et $\beta$ est la **loi exponentielle** de paramètre $\frac{1}{\beta}$,
- $\alpha = n$ et $\beta = \frac{1}{\lambda}$ est la **loi d'Erlang** d'ordre $n$ et de paramètre $\lambda$ (loi de la somme de $n$ variables aléatoires indépendantes de même loi exponentielle de paramètre $\lambda$),
- $\alpha = \frac{n}{2}$ et $\beta = 2$ est la **loi du $\chi^2$** à $n$ degrés de liberté (loi de la somme des carrés de $n$ variables aléatoires indépendantes de même loi gaussienne centrée de variance 1).

Les deux derniers résulats sont des conséquences de la partie 3 de la proposition 1.11.

#### La loi de Weibull

Donnons-nous trois paramètres $\alpha$, $\beta$ et $\theta$, les deux premiers étant strictement positifs, le troisième étant positif ou nul. La loi de densité :

$$
f(t) = \frac{\beta}{\alpha}\left(\frac{t-\theta}{\alpha}\right)^{\beta-1}\exp\left\{-\left(\frac{t-\theta}{\alpha}\right)^\beta\right\}
$$

par rapport à la mesure de Lebesgue sur $[\theta, +\infty[$ est la **loi de Weibull** de paramètres $(\alpha, \beta, \theta)$.

Le paramètre $\alpha$ est le paramètre d'échelle, $\beta$ est le paramètre de forme et $\theta$ le paramètre de translation.

Remarquons que le cas $\beta = 1$ et $\theta = 0$ correspond à la loi exponentielle de paramètre $1/\alpha$.

Nous ne considèrerons ici que le cas $\theta = 0$ et nous parlerons alors de la loi de Weibull de paramètres $(\alpha, \beta)$.

**Proposition 1.12** *Soit $T$ une variable aléatoire de loi de Weibull de paramètres $(\alpha, \beta)$. Notons $h$ son taux de hasard.*

1. *Alors :*

$$
h(t) = \frac{\beta}{\alpha}\left(\frac{t}{\alpha}\right)^{\beta-1}, \quad \mathbb{P}(T > t) = \exp\left\{-\left(\frac{t}{\alpha}\right)^\beta\right\},
$$

$$
\mathbb{E}(T) = \alpha\Gamma\left(1+\frac{1}{\beta}\right), \quad \mathrm{var}(T) = \alpha^2\left[\Gamma\left(1+\frac{2}{\beta}\right) - \Gamma^2\left(1+\frac{1}{\beta}\right)\right].
$$

2. *La variable aléatoire $\left(\frac{T}{\alpha}\right)^\beta$ est de loi exponentielle de paramètre 1.*

<!-- page 2 -->

[FIGURE page 2: Graphique comparatif de la densité de la loi de Weibull de paramètres $\theta = 0$, $\alpha = 1$ et $\beta \in \{0.5, 1, 2, 4\}$ sur l'intervalle $[0, 2]$. Pour $\beta=0.5$, la densité décroît depuis $+\infty$ ; pour $\beta=1$, elle correspond à une décroissance exponentielle standard issue de 1 ; pour $\beta=2$ et $\beta=4$, elle présente une forme en cloche asymétrique unimodale.]

**Densité de la loi de Weibull de paramètres $\theta = 0$, $\alpha = 1$ et $\beta$.**

[FIGURE page 2: Graphique comparatif du taux de hasard de la loi de Weibull de paramètres $\theta = 0$, $\alpha = 1$ et $\beta \in \{0.5, 1, 2, 4\}$ sur l'intervalle $[0, 2]$. Pour $\beta=0.5$, $h(t)$ est strictement décroissant (DFR) divergeant en 0 ; pour $\beta=1$, $h(t)=1$ est constant ; pour $\beta=2$, $h(t)=2t$ est une droite croissante (IFR) ; pour $\beta=4$, $h(t)=4t^3$ est strictement convexe croissant (IFR).]

**Taux de hasard de la loi de Weibull de paramètres $\theta = 0$, $\alpha = 1$ et $\beta$.**

La loi de Weibull est très utilisée en fiabilité. La première raison est que cette loi est apparue "expérimentalement" lors d'études sur le taux de défaillance de matériels. Notons $\lambda$ le taux de défaillance du matériel étudié et $\Lambda$ son taux de défaillance cumulé : $\Lambda(t) = \int_0^t \lambda(s)\,ds$. J. Duane a constaté que $\log \Lambda(t)$ était approximativement une fonction linéaire de $\log t$ (voir à ce propos l'exercice 3.1).

On peut également considérer que le théorème des valeurs extrêmes fournit une explication mathématique au fait que cette loi se rencontre "dans la nature". Il peut s'énoncer par exemple sous la forme suivante ([11], chapitre 8, théorème 3.3) :

**Théorème 1.13 (extrait du théorème des valeurs extrêmes)** *Considérons des variables aléatoires $(X_k)_{k \ge 1}$ indépendantes de même loi de fonction de répartition $F$. Supposons qu'il existe :*

*a. un réel $x_0$ pour lequel $F(x_0) = 0$ et $F(x) > 0$ pour tout $x > x_0$,*

*b. un réel $\beta > 0$ tel que, pour tout $x > 0$ :*

$$
\lim_{t \downarrow 0} \frac{F(xt+x_0)}{F(t+x_0)} = x^\beta.
$$

*Alors il existe des suites de réels $a_n$ et $b_n$ ($a_n > 0$) telles que la suite de variables aléatoires :*

$$
Z_n = \min_{1 \le k \le n} (a_n X_k + b_n)
$$

*converge en loi, lorsque $n$ tend vers l'infini, vers une variable aléatoire de loi de Weibull de paramètre de forme $\beta$.*

La variable $Z_n$ peut représenter la durée de vie d'un matériel décomposé en $n$ éléments, le matériel étant défaillant dès que l'un des éléments le constituant est défaillant (imaginer par exemple un câble découpé "virtuellement" en $n$ tronçons).

**Exemple 1.14 :** Si les $X_k$ sont de loi gamma de paramètres de forme $\alpha$, alors $Z_n$ converge en loi vers une variable aléatoire de loi de Weibull de paramètre de forme $\alpha$. En effet, $x_0 = 0$ et en appliquant la formule généralisée des accroissements finis (règle de l'hôpital), on obtient :

$$
\lim_{t \downarrow 0} \frac{F(xt)}{F(t)} = \lim_{t \downarrow 0} \frac{F(xt)-F(0)}{F(t)-F(0)} = \lim_{t \downarrow 0} \frac{x f(xt)}{f(t)} = x^\alpha.
$$
