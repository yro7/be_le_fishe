<!-- page 1 -->

### 3.3.2 Cas d'un processus de Weibull

Nous observons les $n$ premiers points $T_1, \ldots, T_n$ d'un processus de Weibull de paramètres $\alpha$ et $\beta$.

D'après la proposition 2.25, la log-vraisemblance est :

$$
L(\alpha, \beta) = n(\log \beta - \beta \log \alpha) + (\beta - 1)\sum_{i=1}^n \log T_i - \left(\frac{T_n}{\alpha}\right)^\beta. \tag{3.2}
$$

Les estimateurs du maximum de vraisemblance $\hat{\alpha}$ et $\hat{\beta}$ de $\alpha$ et $\beta$ sont donnés par :

$$
\frac{1}{\hat{\beta}} = \log T_n - \frac{1}{n}\sum_{i=1}^n \log T_i,
$$

$$
\log \hat{\alpha} = \log T_n - \frac{1}{\hat{\beta}}\log n.
$$

Cherchons des intervalles de confiance pour $\alpha$ et $\beta$. Posons :

$$
U_i = \Lambda(T_i) = \left(\frac{T_i}{\alpha}\right)^\beta,
$$

ce qui équivaut à :

$$
\log T_i = \log \alpha + \frac{1}{\beta}\log U_i.
$$

La proposition 2.28 montre que $U_1, \ldots, U_n$ sont les $n$ premiers points d'un processus de Poisson homogène de paramètre 1. La relation :

$$
\frac{\beta}{\hat{\beta}} = \log U_n - \frac{1}{n}\sum_{i=1}^n \log U_i
$$

entraîne donc que la loi de $\beta/\hat{\beta}$ ne dépend pas des paramètres $\alpha$ et $\beta$. Plus précisément :

**Proposition 3.19** *Les variables aléatoires $U_n$ et $\frac{\beta}{\hat{\beta}} = \log U_n - \frac{1}{n}\sum_{i=1}^n \log U_i$ sont indépendantes, de lois gamma, de paramètres respectifs $(n, 1)$ et $(n-1, \frac{1}{n})$.*

*Démonstration :* On vérifie facilement que la loi conditionnelle de $(U_1, \ldots, U_{n-1})$ sachant $\{U_n = t\}$ est la loi de la statistique d'ordre de $n-1$ variables aléatoires de loi uniforme sur $[0, t]$. Par suite :

$$
\mathbb{E}\left(e^{iu\beta / \hat{\beta}} \middle/ U_n = t\right) = \mathbb{E}\left(\exp\left\{ iu\left(1 - \frac{1}{n}\right)\log t - \frac{iu}{n}\sum_{i=1}^{n-1} \log X_i \right\}\right),
$$

les variables aléatoires $X_i$ étant indépendantes et de loi uniforme sur $[0, t]$. Notons $\varphi$ la fonction caractéristique de $\log X_i$. Un petit calcul donne :

$$
\varphi(u) = \frac{e^{iu \log t}}{1 + iu},
$$

et :

$$
\mathbb{E}\left(e^{iu\beta / \hat{\beta}} \middle/ U_n = t\right) = \frac{1}{\left(1 - \frac{iu}{n}\right)^{n-1}}.
$$

La loi conditionnelle de $\beta/\hat{\beta}$ sachant $\{U_n = t\}$ est indépendante de $t$, donc les variables aléatoires $\beta/\hat{\beta}$ et $U_n$ sont indépendantes et la loi de $\beta/\hat{\beta}$ est la loi gamma de paramètres $(n-1, \frac{1}{n})$. ♣

En utilisant le lemme 3.4, nous obtenons :

**Corollaire 3.20** *Dans le cas de l'observation des $n$ premiers points d'un processus de Weibull, la loi de $2n\beta/\hat{\beta}$ est la loi du $\chi^2$ à $2n - 2$ degrés de liberté, et les intervalles suivants :*
- $]0, \, \frac{\hat{\beta}}{2n}\chi_\gamma^2(2n-2)]$,
- $[\frac{\hat{\beta}}{2n}\chi_{1-\gamma}^2(2n-2), \, +\infty[$,
- $[\frac{\hat{\beta}}{2n}\chi_{(1-\gamma)/2}^2(2n-2), \, \frac{\hat{\beta}}{2n}\chi_{(1+\gamma)/2}^2(2n-2)]$,

*sont des intervalles de confiance de $\beta$ de niveau $\gamma$.*

La construction d'intervalles de confiance pour $\alpha$ est plus difficile. Nous n'en obtenons pas à partir des tables existantes, car la loi de

$$
\hat{\beta}\log\frac{\hat{\alpha}}{\alpha} = \frac{\log U_n}{\log U_n - \frac{1}{n}\sum_{i=1}^n \log U_i} - \log n
$$

est bien indépendante de $\alpha$ mais ce n'est pas une loi tabulée. Une formule de calcul de sa fonction de répartition fait l'objet de la proposition 3.21 ci-dessous. Si on en déduit numériquement ses fractiles d'ordre $\gamma$, soit $s_\gamma(n)$, alors chacun des intervalles :
- $[\hat{\alpha}\exp\left(-\frac{s_\gamma(n)}{\hat{\beta}}\right), \, +\infty[$,
- $]0, \, \hat{\alpha}\exp\left(-\frac{s_{1-\gamma}(n)}{\hat{\beta}}\right)]$,
- $[\hat{\alpha}\exp\left(-\frac{s_{(1+\gamma)/2}(n)}{\hat{\beta}}\right), \, \hat{\alpha}\exp\left(-\frac{s_{(1-\gamma)/2}(n)}{\hat{\beta}}\right)]$

est un intervalle de confiance pour $\alpha$ de niveau $\gamma$.

<!-- page 2 -->

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

<!-- page 3 -->

- $\frac{\hat{\beta}\sqrt{n}}{\log n}\left(\frac{\hat{\alpha}}{\alpha} - 1\right)$,

*convergent en loi vers une loi gaussienne centrée de variance 1.*

La remarque 3.14 est également valable ici.

**Remarque 3.24 :** Les formules (3.1) et (3.2) sont similaires. On passe de la première à la deuxième en remplaçant $t$ par $T_n$ et $N(t)$ par $n = N(T_n)$. Il est donc normal que les estimateurs $\hat{\alpha}$ et $\hat{\beta}$ dans le cas de l'observation des $n$ premiers points s'obtiennent à partir de ceux correspondant à une observation sur $[0, t]$ par les mêmes transformations. Cependant ces deux familles d'estimateurs n'ont pas les mêmes propriétés car, dans les transformations effectuées, le $t$ qui était constant devient un $T_n$ aléatoire et le $N(t)$ qui était aléatoire devient un $n$ qui est constant. Cela se voit dans la différence entre les corollaires 3.6 et 3.20. Par contre pour les vitesses de convergence les transformations fonctionnent. En effet, si dans la proposition 3.12, on remplace $t$ par $T_n$ (et donc $N(t)$ par $n$), on obtient bien les vitesses de convergence de la proposition 3.22 car $\log t$ devient $\log T_n$ et comme $N(t)/(t/\alpha)^\beta$ tend presque-sûrement vers 1 lorsque $t$ tend vers l'infini (proposition 2.20) il en est de même pour $N(T_n)/(T_n/\alpha)^\beta = n/(T_n/\alpha)^\beta$ ce qui entraîne que $\log T_n$ est équivalent à $(\log n)/\beta$.

## 3.4 Tests

### 3.4.1 Un test d'homogénéité non paramétrique

Nous supposons observer le processus de Poisson sur $[0, t]$ et nous souhaitons tester l'hypothèse $H_0$ : "le processus de Poisson observé est homogène".

Divisons l'intervalle $[0, t]$ en $d$ intervalles de même longueur notés $I_1, \ldots, I_d$. Soit $N_k$ ($1 \leq k \leq d$), le nombre d'observations appartenant au $k^{\text{ème}}$ intervalle :

$$
N_k = \sum_{i=1}^{N(t)}\mathbb{1}_{\{T_i \in I_k\}}.
$$

D'après la proposition 2.3, sous $H_0$ la loi de $(N_1, \ldots, N_d)$ sachant $\{N(t) = n\}$ est la loi de $(\sum_{i=1}^n \mathbb{1}_{\{X_i \in I_1\}}, \ldots, \sum_{i=1}^n \mathbb{1}_{\{X_i \in I_d\}})$, les variables aléatoires $X_i$ étant indépendantes et de loi uniforme sur $[0, t]$. Par conséquent, d'après le théorème classique de convergence vers la loi du $\chi^2$, sous $H_0$ quand $n$ tend vers l'infini la loi, conditionnellement à $\{N(t) = n\}$, de la variable aléatoire

$$
Z_n = \frac{d}{n}\sum_{k=1}^d \left(N_k - \frac{n}{d}\right)^2
$$

tend vers la loi du $\chi^2$ à $d-1$ degrés de liberté.

Donc, en pratique, pour construire un test de niveau approximativement égal à $\gamma$, si l'on observe $\{N(t) = n\}$ avec $n$ grand, on prend comme région critique :

$$
D = \{Z_n > \chi_{1-\gamma}^2(d-1)\}
$$

(rappelons que $\chi_{1-\gamma}^2(d-1)$ est le quantile d'ordre $1-\gamma$ de la loi du $\chi^2$ à $d-1$ degrés de liberté).

### 3.4.2 Un test paramétrique

Nous supposons maintenant que le processus observé est un processus de Weibull de paramètres $\alpha$ et $\beta$. Nous voulons tester $H_0$ "$\beta \leq \beta_0$" contre $H_1$ "$\beta > \beta_0$".

Nous allons nous limiter à la recherche de tests semblables (*similar* en anglais, terminologie de [71] chapitre 4, paragraphe 3) de niveau $\gamma$, c'est-à-dire tels que, pour toutes valeurs des paramètres $(\alpha, \beta)$ situées sur la frontière entre $H_0$ et $H_1$ (c'est-à-dire telles que $\beta = \beta_0$) la probabilité de rejet soit égale à $\gamma$.

D'après la proposition 2.27, la loi conditionnelle de l'observation $(T_1, \ldots, T_n)$ sachant $\{N(t) = n\}$ a pour densité

$$
f_\beta(t_1, \ldots, t_n) = \frac{n!\beta^n}{t^{\beta n}}\prod_{i=1}^n t_i^{\beta-1}
$$

par rapport à la mesure de Lebesgue sur $[0, t]^n$. Par suite, pour $\beta < \beta'$, la fonction

$$
\frac{f_{\beta'}(T_1, \ldots, T_n)}{f_\beta(T_1, \ldots, T_n)} = \left(\frac{\beta'}{\beta}\right)^n t^{(\beta - \beta')n}\exp\left( (\beta' - \beta)\sum_{i=1}^n \log T_i \right) \tag{3.3}
$$

est une fonction croissante de $\sum_{i=1}^n \log T_i$. Nous en déduisons le test suivant :

**Proposition 3.25** *Nous observons un processus de Weibull sur $[0, t]$.*

*1. Pour tester [sic: tester tester] $H_0$ "$\beta \leq \beta_0$" contre $H_1$ "$\beta > \beta_0$", le test de région critique :*

$$
D = \left\{ \sum_{i=1}^{N(t)}\log T_i \geq N(t)\log t - \frac{1}{2\beta_0}\chi_\gamma^2(2N(t)) \right\}
$$

*est de niveau $\gamma$ et est uniformément le plus puissant à son niveau.*

*2. Pour tester [sic: tester tester] $H_0$ "$\beta \geq \beta_0$" contre $H_1$ "$\beta < \beta_0$", le test de région critique :*

$$
D = \left\{ \sum_{i=1}^{N(t)}\log T_i \leq N(t)\log t - \frac{1}{2\beta_0}\chi_{1-\gamma}^2(2N(t)) \right\}
$$

*est de niveau $\gamma$ et est uniformément le plus puissant à son niveau.*
