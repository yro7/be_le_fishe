<!-- page 1 -->

1. $[0, \, \frac{1}{2}\chi_\gamma^2(2n+2)]$ est un intervalle de confiance pour $\Lambda(t)$ de niveau $\gamma$,

2. $[\exp(-\frac{1}{2}\chi_\gamma^2(2n+2)), \, 1]$ est un intervalle de confiance pour $R(t)$ de niveau $\gamma$.

En généralisant au cas non paramétrique la notion d'exhaustivité, on peut voir que $N(t)$ est exhaustif dans le cadre d'une estimation non paramétrique de $\Lambda(t)$ à partir de l'observation du processus de Poisson sur $[0, t]$. En effet, nous pouvons "paramétrer" la fonction $\lambda$ sur $[0, t]$ par son intégrale $\Lambda(t)$ et sa "renormalisée" $\tilde{\lambda}$ :

$$
\tilde{\lambda}(s) = \frac{\lambda(s)}{\Lambda(t)}.
$$

Or la loi de l'observation sachant $\{N(t) = n\}$ a pour densité (proposition 2.27) la fonction :

$$
\frac{n!}{\Lambda(t)^n}\lambda(t_1)\cdots\lambda(t_n) = n! \, \tilde{\lambda}(t_1)\cdots\tilde{\lambda}(t_n),
$$

qui ne contient pas le paramètre $\Lambda(t)$.

### 3.2.3 Cas d'un processus de Weibull

Nous regardons maintenant un processus de Weibull (*power law process* dans la terminologie américaine) de paramètres $\alpha$ et $\beta$ inconnus, c'est-à-dire un processus de Poisson d'intensité :

$$
\lambda(s) = \frac{\beta}{\alpha}\left(\frac{s}{\alpha}\right)^{\beta-1}.
$$

Il modélise les instants successifs de défaillance d'un matériel dont le taux de défaillance est de type Weibull et qui subit de "petites réparations".

La log-vraisemblance des paramètres $\alpha$ et $\beta$ s'écrit, d'après la proposition 2.27 :

$$
L(\alpha, \beta) = (\log \beta - \beta \log \alpha)N(t) + (\beta - 1)\sum_{i=1}^{N(t)}\log T_i - \left(\frac{t}{\alpha}\right)^\beta. \tag{3.1}
$$

Les estimateurs du maximum de vraisemblance n'existent ici que si $N(t)$ est non nul. Dans ce cas, les estimateurs du maximum de vraisemblance de $\alpha$ et $\beta$ sont respectivement $\hat{\alpha}$ et $\hat{\beta}$, avec :

$$
\frac{1}{\hat{\beta}} = \log t - \frac{1}{N(t)}\sum_{i=1}^{N(t)}\log T_i,
$$

$$
\log \hat{\alpha} = \log t - \frac{1}{\hat{\beta}}\log N(t).
$$

Cherchons des intervalles de confiance pour les paramètres $\alpha$ et $\beta$. Il faut essayer de mettre en évidence des variables aléatoires (fonctions des paramètres) dont la loi (ou la loi conditionnelle en une variable observable) ne dépend pas des paramètres à estimer. Nous allons le faire pour $\beta$ à l'aide de la variable aléatoire $\frac{\beta}{\hat{\beta}}$. Par contre nous ne disposons pas de telle variable pour $\alpha$ et donc nous résoudrons seulement plus loin, dans un cadre asymptotique, le problème de la construction d'un intervalle de confiance pour $\alpha$.

Nous utiliserons fréquemment le lemme suivant :

**Lemme 3.4** *Soit $X$ une variable aléatoire de loi gamma de paramètres $(n, \beta)$, alors $\frac{2X}{\beta}$ est de loi gamma de paramètres $n$ et $2$, c'est-à-dire de loi du $\chi^2$ à $2n$ degrés de liberté.*

**Proposition 3.5** *Dans le cas de l'observation d'un processus de Weibull sur $[0, t]$, la loi de $\frac{\beta}{\hat{\beta}}$ sachant $\{N(t) = n\}$ ($n \neq 0$) est la loi gamma de paramètres $(n, \frac{1}{n})$, par conséquent la loi de $\frac{2n\beta}{\hat{\beta}}$ sachant $\{N(t) = n\}$ est la loi du $\chi^2$ à $2n$ degrés de liberté.*

*Démonstration :* D'après la proposition 2.27, la loi de $\frac{\beta}{\hat{\beta}}$ sachant $\{N(t) = n\}$ est la loi de :

$$
Z_n = \beta \log t - \frac{\beta}{n}\sum_{i=1}^n \log X_i,
$$

les variables aléatoires $X_i$ ($1 \leq i \leq n$) étant indépendantes, de même loi de densité

$$
\frac{1}{\Lambda(t)}\frac{\beta}{\alpha^\beta}x^{\beta-1}\mathbb{1}_{[0, t]}(x) = \frac{\beta}{t^\beta}x^{\beta-1}\mathbb{1}_{[0, t]}(x)
$$

par rapport à la mesure de Lebesgue. Un calcul sans difficulté montre alors que la fonction caractéristique de $Z_n$ est égale à :

$$
\mathbb{E}(e^{iuZ_n}) = \frac{1}{\left(1 - \frac{iu}{n}\right)^n},
$$

d'où le résultat en utilisant la proposition 1.11. ♣

**Corollaire 3.6** *Dans le cas de l'observation d'un processus de Weibull sur $[0, t]$, si nous observons $N(t) = n$ ($n \neq 0$), alors chacun des intervalles $I(n)$ suivants :*
- $[0, \, \frac{\hat{\beta}}{2n}\chi_\gamma^2(2n)]$,
- $[\frac{\hat{\beta}}{2n}\chi_{1-\gamma}^2(2n), \, +\infty[$,
- $[\frac{\hat{\beta}}{2n}\chi_{(1-\gamma)/2}^2(2n), \, \frac{\hat{\beta}}{2n}\chi_{(1+\gamma)/2}^2(2n)]$,

<!-- page 2 -->

*est un intervalle de confiance de $\beta$ de niveau $\gamma$, c'est-à-dire que :*

$$
\mathbb{P}(\beta \in I(n) / N(t) = n) = \gamma.
$$

*Par suite $\mathbb{P}(\beta \in I(N(t)) / N(t) \geq 1) = \gamma$.*

Si nous nous intéressons à une estimation de la fiabilité à l'instant $t$, (estimation ponctuelle ou par intervalle de confiance), nous cherchons à estimer $\Lambda(t) = \left(\frac{t}{\alpha}\right)^\beta$. Il est alors naturel de prendre comme estimateur $\left(\frac{t}{\hat{\alpha}}\right)^{\hat{\beta}}$. Lorsque $N(t) \neq 0$,

$$
\left(\frac{t}{\hat{\alpha}}\right)^{\hat{\beta}} = N(t).
$$

Nous sommes donc conduits à prendre $N(t)$ comme estimateur de $\Lambda(t)$ (que $N(t)$ soit nul ou non) : c'est l'estimateur que nous avons proposé dans le cas non paramétrique et la proposition 3.3 s'applique.

Si nous souhaitons une estimation de la fiabilité à un instant $t' \neq t$, ou de façon équivalente une estimation de :

$$
\Lambda(t') = \left(\frac{t'}{t}\right)^\beta \Lambda(t),
$$

nous pouvons prendre :

$$
\hat{\Lambda}_1(t') = \left(\frac{t'}{t}\right)^{\hat{\beta}}N(t).
$$

Si $t' < t$, nous pouvons également prendre $\hat{\Lambda}_2(t') = N(t')$. Pour des développements sur ce thème, le lecteur se reportera à l'exercice 3.2.

Intéressons-nous maintenant au comportement asymptotique, lorsque $t$ tend vers l'infini, des estimateurs $\hat{\alpha}$ et $\hat{\beta}$. Pour ne pas alourdir les notations, nous ne faisons pas apparaître la dépendance en $t$ de ces estimateurs.

Puisque, d'après la proposition 2.20, $N(t)$ tend presque-sûrement vers l'infini lorsque $t$ tend vers l'infini, le problème de la non-définition des estimateurs $\hat{\alpha}$ et $\hat{\beta}$ pour $N(t)$ nul ne pose plus de problème. Nous omettrons donc dans les démonstrations de faire état du cas $N(t) = 0$.

Commençons par la consistance. Notre démonstration s'appuie sur la proposition F.40 qui donne une condition suffisante de convergence vers 0 d'une martingale. Pour avoir des précisions sur la terminologie employée, se reporter à l'appendice F.

Posons :

$$
\Lambda(t) = \int_0^t \lambda(s) \, ds = \left(\frac{t}{\alpha}\right)^\beta,
$$

et :

$$
M(t) = N(t) - \Lambda(t).
$$

**Proposition 3.7** *Lorsque $t$ tend vers l'infini, alors :*
- $\hat{\alpha}$ *converge presque-sûrement vers $\alpha$,*
- $\hat{\beta}$ *converge presque-sûrement vers $\beta$.*

*Démonstration :* Commençons par montrer la consistance de $\hat{\beta}$. Un petit calcul conduit à :

$$
\frac{1}{\hat{\beta}} = \frac{\Lambda(t)}{N(t)} \left( \frac{\log t}{\Lambda(t)}M(t) - \frac{1}{\Lambda(t)}\int_0^t \log s \, dM(s) + \frac{1}{\beta} \right),
$$

car :

$$
\int_0^t \log s \, d\Lambda(s) = \Lambda(t)\left(\log t - \frac{1}{\beta}\right).
$$

D'après la proposition 2.20, $\frac{N(t)}{\Lambda(t)}$ converge presque-sûrement vers 1 lorsque $t$ tend vers l'infini.

D'autre part, les propositions 2.38 et F.47 entraînent que $M$ est une martingale de processus croissant $\Lambda$ et que $\int_0^\cdot \log s \, dM(s)$ est une martingale de processus croissant $A$ donné par :

$$
\begin{aligned}
A(t) &= \int_0^t (\log s)^2 \, d\Lambda(s) \\
&= \Lambda(t)\left[ \log^2 t - \frac{2}{\beta}\log t + \frac{2}{\beta^2} \right] \\
&\approx \Lambda(t)\log^2 t.
\end{aligned}
$$

Soit $0 < \varepsilon < \frac{1}{2}$, la proposition F.40 entraîne que

$$
\frac{\log t}{\Lambda(t)}M(t) = \frac{\log t}{\Lambda(t)^{1/2-\varepsilon}} \frac{M(t)}{\Lambda(t)^{1/2+\varepsilon}},
$$

et que

$$
\frac{1}{\Lambda(t)}\int_0^t \log s \, dM(s) = \frac{A(t)^{1/2+\varepsilon}}{\Lambda(t)} \frac{\int_0^t \log s \, dM(s)}{A(t)^{1/2+\varepsilon}}
$$

convergent presque-sûrement vers 0 lorsque $t$ tend vers l'infini, ce qui achève de prouver que $1/\hat{\beta}$ converge presque-sûrement vers $1/\beta$.

Montrons maintenant la convergence de $\log \hat{\alpha}$ vers $\log \alpha$. Nous avons :

$$
\begin{aligned}
\log \hat{\alpha} &= \log t - \frac{1}{\hat{\beta}}\log N(t) \\
&= \log t - \frac{1}{\hat{\beta}}\log \frac{N(t)}{\Lambda(t)} - \frac{1}{\hat{\beta}}\log \Lambda(t) \\
&= \frac{\beta}{\hat{\beta}}\log \alpha - \frac{1}{\hat{\beta}}\log \frac{N(t)}{\Lambda(t)} - \left(\frac{1}{\hat{\beta}} - \frac{1}{\beta}\right)\beta \log t.
\end{aligned}
$$

<!-- page 3 -->

Nous savons déjà que $\hat{\beta}$ et $\frac{N(t)}{\Lambda(t)}$ convergent presque-sûrement respectivement vers $\beta$ et 1, il reste donc à montrer la convergence vers 0 de $\left(\frac{1}{\hat{\beta}} - \frac{1}{\beta}\right)\beta \log t$.

Or :

$$
\left(\frac{1}{\hat{\beta}} - \frac{1}{\beta}\right)\log t = \frac{\Lambda(t)}{N(t)} \left[ \frac{(\log t - \frac{1}{\beta})\log t}{\Lambda(t)^{1/2-\varepsilon}} \frac{M(t)}{\Lambda(t)^{1/2+\varepsilon}} - \frac{A(t)^{1/2+\varepsilon}\log t}{\Lambda(t)} \frac{\int_0^t \log s \, dM(s)}{A(t)^{1/2+\varepsilon}} \right]
$$

et nous appliquons à nouveau la proposition F.40. ♣

Regardons maintenant la normalité asymptotique des estimateurs.

Nous utiliserons les résultats suivants dont les principes de démonstrations figurent dans [35] (chapitres 0 et 3).

**Proposition 3.8** *L'ensemble $T$ désigne soit $\mathbb{R}_+$ soit $\mathbb{N}$. Considérons deux familles $(X(t))_{t \in T}$ et $(Y(t))_{t \in T}$ de variables aléatoires à valeurs respectivement dans $\mathbb{R}^p$ et $\mathbb{R}^q$ et $f$ une fonction continue définie sur $\mathbb{R}^{p+q}$ et à valeurs dans $\mathbb{R}^r$. Supposons que lorsque $t$ tend vers l'infini, $X(t)$ converge en loi vers $X$, alors :*
1. *Si $Y(t)$ converge en probabilité vers une constante $a$, alors $f(X(t), Y(t))$ converge en loi vers $f(X, a)$,*
2. *si $p = q = 1$ et si $Y(t)$ converge en probabilité vers 0, alors $X(t) Y(t)$ converge en probabilité vers 0.*

**Proposition 3.9** *L'ensemble $T$ désigne soit $\mathbb{R}_+$ soit $\mathbb{N}$. Soit $(X(t))_{t \in T}$ une famille de variables aléatoires à valeurs dans $\mathbb{R}^k$, $c(t)$ une famille de réels tendant vers l'infini et $m \in \mathbb{R}^k$. On suppose que :*

$$
X(t) \xrightarrow[t \to +\infty]{\mathbb{P}} m,
$$

$$
c(t)(X(t) - m) \xrightarrow[t \to +\infty]{\mathcal{L}} \mathcal{N}(0, \Gamma),
$$

*où $\mathcal{N}(0, \Gamma)$ désigne une variable aléatoire de loi gaussienne centrée de dispersion $\Gamma$.*

*Soit $f$ une fonction définie sur un voisinage $U$ de $m$, à valeurs dans $\mathbb{R}^p$, deux fois différentiable. On suppose que ses dérivées secondes sont bornées sur $U$ et on note $Jf$ la matrice jacobienne de $f$.*

*Alors :*

$$
c(t)(f(X(t)) - f(m)) \xrightarrow[t \to +\infty]{\mathcal{L}} \mathcal{N}(0, K),
$$

*où $\mathcal{N}(0, K)$ désigne une variable aléatoire de loi gaussienne centrée de dispersion :*

$$
K = Jf(m)\Gamma Jf(m)^T.
$$

Lorsque les différentes composantes du vecteur $X_t$ ne convergent pas à la même vitesse, on peut adapter la démonstration de la proposition 3.9 et on obtient :

**Proposition 3.10** *L'ensemble $T$ désigne soit $\mathbb{R}_+$ soit $\mathbb{N}$. Soit $(X(t))_{t \in T}$ une famille de variables aléatoires à valeurs dans $\mathbb{R}^k$, $c(t)$ une famille d'éléments de $\mathbb{R}^k$ tendant vers l'infini et $m \in \mathbb{R}^k$. Notons respectivement $X_i(t)$, $c_i(t)$ et $m_i$ la $i^{\text{ème}}$ composante des vecteurs $X(t)$, $c(t)$ et $m$. On suppose que :*

$$
X(t) \xrightarrow[t \to +\infty]{\mathbb{P}} m,
$$

$$
(c_1(t)(X_1(t) - m_1), \ldots, c_k(t)(X_k(t) - m_k)) \xrightarrow[t \to +\infty]{\mathcal{L}} \mathcal{N}(0, \Gamma),
$$

*où $\mathcal{N}(0, \Gamma)$ désigne une variable aléatoire de loi gaussienne centrée de dispersion $\Gamma$.*

*Soit $f_i$ ($1 \leq i \leq k$) des fonctions à valeurs réelles définies respectivement sur des voisinages de $m_i$, deux fois dérivables et avec des dérivées secondes bornées sur ces voisinages.*

*Alors :*

$$
(c_1(t)(f_1(X_1(t)) - f_1(m_1)), \ldots, c_k(t)(f_k(X_k(t)) - f_k(m_k))) \xrightarrow[t \to +\infty]{\mathcal{L}} \mathcal{N}(0, K),
$$

*avec $K(i, j) = f_i'(m_i)f_j'(m_j)K(i, j)$ [sic], c'est-à-dire $K = Jf(m)\Gamma Jf(m)^T$ où la matrice $Jf$ est la matrice jacobienne de la fonction $f$, de $\mathbb{R}^k$ dans $\mathbb{R}^k$, définie par $f(m_1, \ldots, m_k) = (f_1(m_1), \ldots, f_k(m_k))$.*

**Proposition 3.11** *Lorsque $t$ tend vers l'infini, alors :*

*1.a $\frac{\sqrt{\Lambda(t)}}{\log t}(\log \hat{\alpha} - \log \alpha)$ converge en loi vers une loi gaussienne centrée de variance 1,*

*1.b $\frac{\sqrt{\Lambda(t)}}{\log t}(\hat{\alpha} - \alpha)$ converge en loi vers une loi gaussienne centrée de variance $\alpha^2$,*

*2.a $\sqrt{\Lambda(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right)$ converge en loi vers une loi gaussienne centrée de variance 1,*

*2.b $\sqrt{\Lambda(t)}(\hat{\beta} - \beta)$ converge en loi vers une loi gaussienne centrée de variance $\beta^2$,*

<!-- page 4 -->

*3.a le vecteur bidimensionnel $\left(\frac{\sqrt{\Lambda(t)}}{\log t}(\log \hat{\alpha} - \log \alpha), \, \sqrt{\Lambda(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right)\right)$ converge en loi vers $(X, -X)$ où $X$ est une variable aléatoire gaussienne centrée réduite, c'est-à-dire vers la loi gaussienne bidimensionnelle centrée de matrice de dispersion $\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$,*

*3.b le vecteur bidimensionnel $\left(\frac{\sqrt{\Lambda(t)}}{\log t}(\hat{\alpha} - \alpha), \, \sqrt{\Lambda(t)}(\hat{\beta} - \beta)\right)$ converge en loi vers la loi gaussienne bidimensionnelle centrée de matrice de dispersion $\begin{pmatrix} \alpha^2 & \alpha\beta \\ \alpha\beta & \beta^2 \end{pmatrix}$.*

*Démonstration :* Commençons par démontrer les résultats concernant $\hat{\beta}$. D'après les propositions 2.20 et 3.8, pour obtenir 2.a il suffit de démontrer la convergence en loi de

$$
Z_t = \frac{N(t)}{\Lambda(t)}\sqrt{\Lambda(t)}\left[\frac{\beta}{\hat{\beta}} - 1\right]
$$

vers une gaussienne centrée de variance 1. Pour cela calculons la fonction caractéristique de $Z_t$. En utilisant la proposition 3.5 et la proposition 1.11 nous obtenons :

$$
\begin{aligned}
\mathbb{E}(e^{iuZ_t} / N(t) = n) &= e^{-iun/\sqrt{\Lambda(t)}}\mathbb{E}\left(e^{iun\beta / (\hat{\beta}\sqrt{\Lambda(t)})} \middle/ N(t) = n\right) \\
&= \frac{e^{-iun/\sqrt{\Lambda(t)}}}{\left(1 - \frac{iu}{\sqrt{\Lambda(t)}}\right)^n}.
\end{aligned}
$$

Par suite :

$$
\begin{aligned}
\mathbb{E}(e^{iuZ_t}) &= \sum_{n \geq 0} \frac{e^{-iun/\sqrt{\Lambda(t)}}}{\left(1 - \frac{iu}{\sqrt{\Lambda(t)}}\right)^n} e^{-\Lambda(t)}\frac{\Lambda(t)^n}{n!} \\
&= \exp\left\{ -\Lambda(t) + \frac{\Lambda(t)}{1 - \frac{iu}{\sqrt{\Lambda(t)}}} e^{-iu/\sqrt{\Lambda(t)}} \right\} \\
&= \exp\{\Phi(t)\},
\end{aligned}
$$

avec :

$$
\begin{aligned}
\Phi(t) &= \frac{\Lambda(t)}{1 - \frac{iu}{\sqrt{\Lambda(t)}}} \left[ e^{-iu/\sqrt{\Lambda(t)}} - 1 + \frac{iu}{\sqrt{\Lambda(t)}} \right] \\
&\simeq -\frac{u^2}{2},
\end{aligned}
$$

ce qui achève la démonstration de 2.a.

Le résultat 2.b découle de 2.a et de la proposition 3.9.

Pour 1.a, nous écrivons :

$$
\begin{aligned}
\frac{\sqrt{\Lambda(t)}}{\log t}(\log \hat{\alpha} - \log \alpha) &= \frac{\log \alpha}{\log t}\sqrt{\Lambda(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right) - \frac{1}{\hat{\beta}\log t}\sqrt{\Lambda(t)}\log\frac{N(t)}{\Lambda(t)} - \sqrt{\Lambda(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right) \\
&= H(t) - \sqrt{\Lambda(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right).
\end{aligned}
$$

D'après les propositions 2.20 et 3.9, $\sqrt{\Lambda(t)}\log\frac{N(t)}{\Lambda(t)}$ converge en loi. En appliquant les résultats sur la convergence de $\hat{\beta}$ et la proposition 3.8, nous voyons que $H(t)$ tend en probabilité vers 0 lorsque $t$ tend vers l'infini, puis que $\frac{\sqrt{\Lambda(t)}}{\log t}(\log \hat{\alpha} - \log \alpha)$ converge en loi vers $-X$, la variable $X$ étant gaussienne centrée de variance 1. La loi de $-X$ étant égale à la loi de $X$, l'assertion 1.a est établie.

L'assertion 1.b résulte encore une fois de la proposition 3.9.

Pour démontrer 3.a, posons :

$$
X(t) = \left( \frac{\sqrt{\Lambda(t)}}{\log t}(\log \hat{\alpha} - \log \alpha), \, \sqrt{\Lambda(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right) \right),
$$

$$
Y(t) = \left( -\sqrt{\Lambda(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right), \, \sqrt{\Lambda(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right) \right),
$$

alors :

$$
X(t) = (H(t), 0) + Y(t).
$$

La proposition 3.8 montre qu'il suffit de prouver la convergence en loi de $Y(t)$, ce qui est immédiat en utilisant 2.a.

Le résultat 3.b se déduit de 3.a par la proposition 3.10. ♣

Les vitesses de convergence des estimateurs dépendant des paramètres, les résultats ci-dessus ne permettent pas de construire des intervalles de confiance asymptotiques. Nous remédions à cette difficulté dans la proposition ci-dessous.

**Proposition 3.12** *Soit $X$ une variable aléatoire gaussienne centrée de variance 1. Lorsque $t$ tend vers l'infini, les variables aléatoires suivantes :*
- $\frac{\sqrt{N(t)}}{\log t}\log\frac{\hat{\alpha}}{\alpha}$,

<!-- page 5 -->

- $\frac{\sqrt{N(t)}}{\log t}\left(\frac{\hat{\alpha}}{\alpha} - 1\right)$,
- $\sqrt{N(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right)$,
- $\sqrt{N(t)}\left(\frac{\hat{\beta}}{\beta} - 1\right)$

convergent en loi vers $X$ et $\left( \frac{\sqrt{N(t)}}{\log t}(\log \hat{\alpha} - \log \alpha), \, \sqrt{N(t)}\left(\frac{\beta}{\hat{\beta}} - 1\right) \right)$ converge en loi vers $(X, -X)$.

La démonstration de cette proposition résulte de la proposition précédente et des propositions 3.8 et 2.20.

Lorsque $t$ est grand, la construction d'intervalles de confiance ne pose donc pas plus de problème pour $\alpha$ que pour $\beta$. Soit $n_\gamma$ le quantile d'ordre $\gamma$ de la loi normale centrée réduite ($\mathbb{P}(X \leq n_\gamma) = \gamma$, $X$ étant de loi gaussienne de moyenne 0 et de variance 1). Nous obtenons :

**Corollaire 3.13** *Lorsqu'on observe un processus de Weibull sur $[0, t]$, les intervalles et pavés suivants sont des régions de confiance asymptotiques (lorsque $t$ tend vers l'infini) de niveau $\gamma$ :*

*1. pour $\alpha$ :*
$$
\text{(a) } [0, \, \hat{\alpha}\exp(-\frac{\log t}{\sqrt{N(t)}}n_\gamma)],
$$
$$
\text{(b) } [\hat{\alpha}\exp(-\frac{\log t}{\sqrt{N(t)}}n_\gamma), \, +\infty[,
$$
$$
\text{(c) } [\hat{\alpha}\exp(-\frac{\log t}{\sqrt{N(t)}}n_{(1+\gamma)/2}), \, \hat{\alpha}\exp(\frac{\log t}{\sqrt{N(t)}}n_{(1+\gamma)/2})],
$$
$$
\text{(d) } [0, \, \hat{\alpha} + \frac{\log t}{\sqrt{N(t)}}n_\gamma],
$$
$$
\text{(e) } [\hat{\alpha} - \frac{\log t}{\sqrt{N(t)}}n_\gamma, \, +\infty[,
$$
$$
\text{(f) } [\hat{\alpha} - \frac{\log t}{\sqrt{N(t)}}n_{(1+\gamma)/2}, \, \hat{\alpha} + \frac{\log t}{\sqrt{N(t)}}n_{(1+\gamma)/2}],
$$

*2. pour $\beta$ :*
$$
\text{(a) } [0, \, \hat{\beta}(1 + \frac{n_\gamma}{\sqrt{N(t)}})],
$$
$$
\text{(b) } [\hat{\beta}(1 - \frac{n_\gamma}{\sqrt{N(t)}}), \, +\infty[,
$$
$$
\text{(c) } [\hat{\beta}(1 - \frac{n_{(1+\gamma)/2}}{\sqrt{N(t)}}), \, \hat{\beta}(1 + \frac{n_{(1+\gamma)/2}}{\sqrt{N(t)}})],
$$
$$
\text{(d) } [0, \, \hat{\beta} + \frac{n_\gamma}{\sqrt{N(t)}}],
$$
$$
\text{(e) } [\hat{\beta} - \frac{n_\gamma}{\sqrt{N(t)}}, \, +\infty[,
$$
$$
\text{(f) } [\hat{\beta} - \frac{n_{(1+\gamma)/2}}{\sqrt{N(t)}}, \, \hat{\beta} + \frac{n_{(1+\gamma)/2}}{\sqrt{N(t)}}],
$$

*3. pour $(\alpha, \beta)$ :*
$$
\text{(a) } [0, \, \hat{\alpha}\exp(-\frac{\log t}{\sqrt{N(t)}}n_\gamma)] \times [0, \, \hat{\beta}(1 + \frac{n_\gamma}{\sqrt{N(t)}})],
$$
$$
\text{(b) } [\hat{\alpha}\exp(-\frac{\log t}{\sqrt{N(t)}}n_\gamma), \, +\infty[ \times [\hat{\beta}(1 - \frac{n_\gamma}{\sqrt{N(t)}}), \, +\infty[,
$$
$$
\begin{aligned}
\text{(c) } &[\hat{\alpha}\exp(-\frac{\log t}{\sqrt{N(t)}}n_{(1+\gamma)/2}), \, \hat{\alpha}\exp(\frac{\log t}{\sqrt{N(t)}}n_{(1+\gamma)/2})] \\
&\times [\hat{\beta}(1 - \frac{n_{(1+\gamma)/2}}{\sqrt{N(t)}}), \, \hat{\beta}(1 + \frac{n_{(1+\gamma)/2}}{\sqrt{N(t)}})].
\end{aligned}
$$

**Remarque 3.14 :** Jusqu'à maintenant nous avons supposé que les paramètres $\alpha$ et $\beta$ étaient inconnus. Lorsque $\alpha$ est connu, l'estimateur du maximum de vraisemblance de $\beta$, est solution de l'équation :

$$
\left(\frac{1}{\beta} - \log \alpha\right)N(t) + \sum_{i=1}^{N(t)}\log T_i - \left(\frac{t}{\alpha}\right)^\beta \log \frac{t}{\alpha} = 0,
$$

il n'a donc pas de forme explicite. Nous pouvons cependant travailler avec l'estimateur $\hat{\beta}$ donné ci-dessus, bien qu'il ne soit plus l'estimateur du maximum de vraisemblance. Ses propriétés, que nous avons données dans les propositions précédentes, et les résultats qui en découlent restent valables.

### 3.2.4 Cas de $p$ processus de Weibull indépendants

Si nous observons non pas un processus de Weibull sur $[0, t]$, mais $p$ processus de Weibull indépendants sur $[0, t]$, l'observation globale est celle de la superposition de $p$ processus de Poisson indépendants de même intensité $\lambda(u) = \frac{\beta}{\alpha}\left(\frac{u}{\alpha}\right)^{\beta-1}$, c'est donc, d'après la proposition 2.47, un processus de Poisson d'intensité $p\lambda(t)$.

Notons $(\tilde{T}_n)_{n \geq 1}$ le processus global et $\tilde{N}$ sa fonction de comptage. Un calcul facile montre que les estimateurs du maximum de vraisemblance $\hat{\alpha}$ et $\hat{\beta}$ de $\alpha$ et $\beta$ sont donnés par :

<!-- page 6 -->

[NOTE : Dans le document source scanné, cette page saute des pages 72-73 aux pages 78-79 du livre (Section 3.3.2).]

**Proposition 3.21** *Dans le cas de l'observation des $n$ premiers points d'un processus de Weibull, nous avons :*

$$
\mathbb{P}\left(\hat{\beta}\log\frac{\hat{\alpha}}{\alpha} \leq x\right) = \int_0^{+\infty} F_{2n}\left(2e^{(x+\log n)t}\right)\frac{n^{n-1}}{(n-2)!}t^{n-2}e^{-tn} \, dt,
$$

*où $F_{2n}$ est la fonction de répartition de loi du $\chi^2$ à $2n$ degrés de liberté.*

*Démonstration :* Posons :

$$
Y = \hat{\beta}\log\frac{\hat{\alpha}}{\alpha} + \log n = \frac{\hat{\beta}}{\beta}\log U_n.
$$

La loi de $U_n$ est la loi gamma de paramètres $(n, 1)$, par conséquent (d'après le lemme 3.4) la loi de $X = 2U_n$ est la loi du $\chi^2$ à $2n$ degrés de liberté et la proposition 3.19 donne :

$$
\begin{aligned}
\mathbb{P}(Y \leq y) &= \mathbb{P}\left(2U_n \leq 2e^{\beta y / \hat{\beta}}\right) \\
&= \int_0^{+\infty} \mathbb{P}(X \leq 2e^{ty})g(t) \, dt
\end{aligned}
$$

où $g$ est la densité de la loi gamma de paramètres $(n-1, \frac{1}{n})$. ♣

Si $n$ est grand, il est plus simple d'utiliser le comportement asymptotique des estimateurs pour obtenir les intervalles de confiance.

**Proposition 3.22** *Lorsque $n$ tend vers l'infini, alors :*

*1.a $\hat{\alpha}$ converge presque-sûrement vers $\alpha$,*

*1.b $\frac{\beta\sqrt{n}}{\log n}(\log\hat{\alpha} - \log\alpha)$ converge en loi vers une loi gaussienne centrée de variance 1,*

*2.a $\hat{\beta}$ converge presque-sûrement vers $\beta$,*

*2.b $\sqrt{n}\left(\frac{\beta}{\hat{\beta}} - 1\right)$ converge en loi vers une loi gaussienne centrée de variance 1,*

*3. le vecteur $\left(\frac{\beta\sqrt{n}}{\log n}(\log\hat{\alpha} - \log\alpha), \, \sqrt{n}\left(\frac{\beta}{\hat{\beta}} - 1\right)\right)$ converge en loi vers la loi gaussienne bidimensionnelle centrée de dispersion $\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$, et le vecteur $\left(\frac{\beta\sqrt{n}}{\log n}(\hat{\alpha} - \alpha), \, \sqrt{n}(\hat{\beta} - \beta)\right)$ converge en loi vers la loi gaussienne bidimensionnelle centrée de dispersion $\begin{pmatrix} \alpha^2 & \alpha\beta \\ \alpha\beta & 1\beta^2 \end{pmatrix}$.*

*Démonstration :* Posons :

$$
M(t) = \sum_{i=1}^n \mathbb{1}_{\{U_i \leq t\}} - t, \quad \tilde{M}(t) = \int_0^t \log s \, dM(s).
$$

D'après les propositions 2.38 et F.47, les processus $M$ et $\tilde{M}$ sont des martingales de processus croissants respectifs $t$ et $A(t) = \int_0^t \log^2 s \, ds \simeq t \log^2 t$. Nous en déduisons (proposition F.40) que pour tout $\varepsilon > 0$, $\frac{M(t)}{t^{1/2+\varepsilon}}$ et $\frac{\tilde{M}(t)}{(t\log^2 t)^{1/2+\varepsilon}}$ convergent presque-sûrement vers 0 lorsque $t$ tend vers l'infini. Posons :

$$
X_n = M(U_n) = n - U_n,
$$

$$
Y_n = \tilde{M}(U_n) = \sum_{i=1}^n \log U_i - \int_0^{U_n} \log s \, ds = \sum_{i=1}^n \log U_i - U_n(\log U_n - 1).
$$

Puisque $U_n$ tend presque-sûrement vers l'infini lorsque $n$ tend vers l'infini, il s'ensuit que $\frac{X_n}{U_n^{1/2+\varepsilon}}$ et $\frac{Y_n}{(U_n \log^2 U_n)^{1/2+\varepsilon}}$ tendent presque-sûrement vers 0 lorsque $n$ tend vers l'infini.

En prenant $0 < \varepsilon < \frac{1}{4}$, nous en déduisons que :

$$
\left(\frac{1}{\hat{\beta}} - \frac{1}{\beta}\right)\beta \log n = \frac{X_n}{U_n^{1/2+2\varepsilon}}\left(\frac{U_n}{n}\right)^{1/2+2\varepsilon}\frac{\log n}{n^{1/2-2\varepsilon}U_n^\varepsilon} - \frac{Y_n}{U_n^{1/2+\varepsilon}}\left(\frac{U_n}{n}\right)^{1/2+\varepsilon}\frac{\log n}{n^{1/2-\varepsilon}}
$$

converge presque-sûrement vers 0 lorsque $n$ tend vers l'infini car $U_n/n$ tend presque-sûrement vers 1.

Ceci montre d'une part que $\frac{1}{\hat{\beta}}$ converge presque-sûrement vers $\frac{1}{\beta}$ et d'autre part que :

$$
\log \hat{\alpha} = \log \alpha + \frac{1}{\beta}\log\frac{U_n}{n} - \left(\frac{1}{\hat{\beta}} - \frac{1}{\beta}\right)\log n
$$

converge presque-sûrement vers $\log \alpha$.

Le schéma de la démonstration de la normalité asymptotique est donné dans l'exercice 3.6. ♣

Pour le calcul d'intervalles de confiance asymptotiques pour $\alpha$ il est préférable d'utiliser la proposition 3.8 pour obtenir une formulation ne faisant pas apparaître le paramètre $\beta$, qu'on remplace par son estimation $\hat{\beta}$, ce qui donne par exemple :

**Proposition 3.23** *Lorsque $n$ tend vers l'infini, les variables aléatoires suivantes :*
- $\frac{\hat{\beta}\sqrt{n}}{\log n}\log\frac{\hat{\alpha}}{\alpha}$,
