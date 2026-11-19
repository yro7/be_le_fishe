<!-- page 1 -->

# The Rasch-Weibull Process

**Author(s):** Svend Kreiner Møller  
**Source:** *Scandinavian Journal of Statistics*, 1976, Vol. 3, No. 3 (1976), pp. 107-115  
**Published by:** Wiley on behalf of Board of the Foundation of the Scandinavian Journal of Statistics  
**Stable URL:** http://www.jstor.com/stable/4615622  

<!-- page 2 -->

*Scand J Statist* 3: 107-115, 1976

# The Rasch-Weibull Process

**SVEND KREINER MØLLER**  
*The Danish Institute for Educational Research*  
Received June 1974, revised December 1975

**ABSTRACT.** The Rasch-Weibull process—a non-homogeneous Poisson process with intensity $\lambda \gamma t^{\gamma-1}$—is introduced, and it is shown that the process has a number of qualities which make a statistical analysis of data which can be described by the Rasch-Weibull process relatively simple. Two examples are included. This process must not be confused with the Weibull Renewal process (see Lomnicki, 1966; Smith & Leadbetter, 1963).

**Key words:** Point process, non-homogeneous Poisson process, intensity function, Rasch-Weibull, characterization

## I. Introduction

The intent of this article is to present a non-homogeneous Poisson process known as the Rasch-Weibull process.

This model was originally described in a comment by Rasch (1970), who suggested that it could be used in connection with the analysis of oral reading speed measured in reading tests then being conducted by the Danish Institute for Educational Research. This project is described in a series of duplicated reports; especially Kreiner Møller (1970).

In the following, the R-W process will be defined and the characteristic features of the process within the class of non-homogeneous Poisson processes will be discussed. Two examples of statistical analysis with the R-W process will be presented: First an analysis of data from Proschan (1963), as presented by Cox & Lewis (1966), and then in connection with data from the above-mentioned tests of oral reading speed.

## II. The model

We will use the same notation as Parzen (1962), chapter 4.

It is assumed that we are working with phenomena which can be described as events occurring randomly in time, in a manner which is consistent with a non-homogeneous Poisson process $\{N(t), t > 0\}$ with strictly monotonally increasing mean-value function $m(t)$ and, if the mean-value function is differentiable, intensity $
u(t) = \partial m(t)/\partial t$.

The corresponding arrival time process is $\{W_n, n = 1, 2, \dots\}$, where $W_n$ represents the waiting time to the occurrence of the $n$th event.

**Definition.** The Rasch-Weibull process with parameters $(\lambda, \gamma)$, $\lambda > 0$ and $\gamma > 0$, is to be understood as the non-homogeneous Poisson process with mean-value function $m(t) = \lambda t^\gamma$ and intensity $
u(t) = \lambda \gamma t^{\gamma-1}$.

In the special case when $\gamma = 1$, the R-W process reduces to an homogeneous Poisson process.

When $\gamma < 1$, the intensity is decreasing, so that as the observation period increases, the intervals between events increase in length.

When $\gamma > 1$, the intensity is increasing, so that the intervals between events on the average become shorter and shorter.

Because of the monotonal character of the intensity, the process can be expected to be relevant in situations with increasing or decreasing "failure rates" as in Barlow & Proschan (1967) and Hollander & Proschan (1972). Notice however that the change in "failure rate" must occur in time, not in correlation to serial number.

Additional important properties of the R-W process will be described below.

### Characterization

A necessary and sufficient condition for a non-homogeneous Poisson process with strictly monotonically increasing mean-value function $m(t)$ to be a Rasch-Weibull process is:

For every $n \ge 2$, the stochastic vector
$$(W_1/W_n, W_2/W_n, \dots, W_{n-1}/W_n)$$
is stochastically independent of $W_n$.

*Proof.* It is known that a non-homogeneous Poisson process can be transformed to a homogeneous Poisson process (Parzen, 1962, p. 126) by transforming the time scale with the aid of the mean-value function $s = m(t)$, since the stochastic process

<!-- page 3 -->

$\{M(s), s \ge 0\}$, defined by $M(s) = N(m^{-1}(s))$, is a homogeneous Poisson process with intensity 1.

The relation between the arrival-time process $\{U_n, n = 1, 2, \dots\}$, which corresponds to $\{M(s), s \ge 0\}$, and the arrival-time process $\{W_n, n = 1, 2, \dots\}$ from the non-homogeneous Poisson process is given by $U_n = m(W_n)$. Therefore, for all $n \ge 2$, the conditional distribution of $(m(W_1), \dots, m(W_{n-1}))$, given that $W_n = w_n$, will be the same as the distribution of $n - 1$ ordered, independent observations from a rectangular distribution over the interval $(0, m(w_n)]$ as shown by Parzen (1962, p. 143).

An immediate result is that the conditional distribution of $(W_1/W_n, \dots, W_{n-1}/W_n)$, given $W_n = w_n$, is the same as the distribution of $n - 1$ ordered, independent observations from the distribution
$$F(x) = rac{m(w_n x)}{m(w_n)} \quad 	ext{for } x \in (0, 1]. 	ag{1}$$

Thus, $(W_1/W_n, \dots, W_{n-1}/W_n)$ is stochastically independent of $W_n$ precisely in the case where the function $m(w_n x)/m(w_n)$ is independent of $w_n$ or, in other words, where $m$ satisfies the functional equation $m(wx)m(1) = m(w)m(x)$.

The positive, continuous, non-decreasing solutions of this functional equation, under the condition that $m(x)$ need only be defined for the set of positive real numbers, are given by $m(x) = \lambda x^\gamma$, where $\lambda$ and $\gamma$ are positive constants (Aczél, 1966, p. 41).

This completes the proof.

## III. Statistical analysis with the Rasch-Weibull process

A point process will most often be observed in one of the following ways:

- **A.** For a given $n$, $W_1, \dots, W_n$ have been observed.
- **B.** In a given interval, $[0, T]$, the number of events $N(T) = n$ and the arrival time $W_1, \dots, W_n$ of these events have been observed.
- **C.** At given points in time $t_1, \dots, t_r$ the number of events $N(t_1), \dots, N(t_r)$ have been recorded.

In all of these three forms of data, the situation is such that we are dealing with statistics which are ancillary for the trend parameter $\gamma$, so that an inference about this parameter should be based upon the conditional distribution of the observations given these ancillary statistics.

The primary aim of the following discussion will be to illustrate questions relating to material in the forms described above as (A) and (C). Only estimation and test of fit of the model will be discussed, but the results should be immediately applicable to hypothesis testing as well. Finally, emphasis is given to inference about the trend parameter, but estimation of the intensity parameter $\lambda$ will also be described.

### Data Form A

It is assumed that $n$ arrival times $W_1, \dots, W_n$ have been observed.

Since the variables $U_1 = \lambda W_1^\gamma, \dots, U_n = \lambda W_n^\gamma$ have the same distributions as the first $n$ arrival times from a homogeneous Poisson process with intensity 1, it follows that $W_n$ conforms to a generalized gamma distribution (Stacy, 1962) with density
$$rac{\lambda^n \gamma w_n^{n\gamma - 1} \exp(-\lambda w_n^\gamma)}{(n-1)!} \quad 	ext{for } w_n > 0$$
and secondly $W_1, \dots, W_n$ has the simultaneous density
$$\gamma^n \left(\prod_{i=1}^n w_i^{\gamma-1}ight) \lambda^n \exp(-\lambda w_n^\gamma) \quad 	ext{for } 0 < w_1 < w_2 < \dots < w_n.$$

For each value of $\gamma$ and for each value of $w_n$, it is possible to determine a $\lambda$ so that $w_n$ is mode in the generalized gamma distribution having the above-mentioned density. Thus, $W_n$ is M-ancillary with respect to $\gamma$ (Barndorff-Nielsen, 1973) which provides an argument for studying inference about the parameter $\gamma$ in the conditional distribution of $W_1, \dots, W_{n-1}$ given $W_n$. This is most easily done by considering the variables $X_1 = -\ln(W_{n-1}/W_n), X_2 = -\ln(W_{n-2}/W_n), \dots, X_{n-1} = -\ln(W_1/W_n)$ given $W_n$, which by virtue of the characterization above, are stochastically independent of $W_n$.

**Theorem 1.** *For all $n \ge 2$, the stochastic vector*
$$\left(-\ln(W_{n-1}/W_n), -\ln(W_{n-2}/W_n), \dots, -\ln(W_1/W_n)ight)$$
*is distributed as $n-1$ ordered, independent observations from an exponential distribution with the scale parameter $\gamma^{-1}$ if and only if the recorded events are from a Rasch-Weibull process with trend parameter $\gamma$.*

*Proof.* The theorem is the immediate consequence of the characterization, where it may be noted from (1) that $(W_1/W_n, \dots, W_{n-1}/W_n)$ is distributed as $n-1$ ordered observations from the distribution having the distribution function
$$F(x) = x^\gamma, \quad x \in (0, 1]. 	ag{2}$$

Consequently, conditional inference on $\gamma$ for given $W_n$ amounts to studying ordered observations from the exponential distribution, which is exhaustively treated in the literature (e.g. Sarhan & Greenberg (1962), chapter 11). It can be shown that the sum of the components of the stochastic vector is sufficient for $\gamma$, and that the sum follows a gamma distribution with scale parameter $\gamma^{-1}$ and shape parameter $n - 1$. Exact statistical analysis about the parameter $\gamma$ is therefore possible.

<!-- page 4 -->

The maximum-likelihood estimate of $\gamma$ is then given by
$$\hat{\gamma} = rac{n-1}{\sum_{i=1}^{n-1} (-\ln(W_{n-i}/W_n))}.$$

In order to achieve an unbiased estimate of $\gamma$, it is necessary to multiply $\hat{\gamma}$ by $rac{n-2}{n-1}$, while multiplication by the factor $rac{n-3}{n-1}$ provides an estimate of $\gamma$ with minimum mean square error.

According to the characterization and Theorem 1, a necessary and sufficient condition for the validity of the R-W process is that $(-\ln(W_{n-1}/W_n), \dots, -\ln(W_1/W_n))$ follows the same distribution as a set of $n-1$ ordered statistics from an exponential distribution with (unknown) scale parameter $\gamma^{-1}$. This provides a test of fit of the R-W process. The test may be made by graphic methods, where the empirical distribution is plotted against the theoretical, or by numerical methods such as a Kolmogorov-Smirnov test (Lilliefors, 1969).

In an analysis of conditions about the intensity parameter $\lambda$, an estimation of the trend parameter is necessary, so that the analysis is not so simple as in the previous case. The maximum-likelihood estimate is given by $\hat{\lambda} = n/W_n^{\hat{\gamma}}$ and, while it is quite possible to determine the simultaneous density of $(\hat{\gamma}, \hat{\lambda})$, this is so complicated as to prevent any immediate gains in information about exact qualities of $\hat{\lambda}$. An evaluation of $\hat{\lambda}$ must therefore be based upon large-sample properties.

### Data Form C

At given points in time $t_1, \dots, t_r$ the corresponding number of events $N(t_1), \dots, N(t_r)$ have been observed.

The joint distribution of these variables is given by the probabilities
$$P(N(t_1) = n_1, N(t_2)-N(t_1) = n_2, \dots, N(t_r)-N(t_{r-1}) = n_r)$$
$$= \prod_{i=1}^r rac{\lambda^{n_i} (t_i^\gamma - t_{i-1}^\gamma)^{n_i} \exp(-\lambda(t_i^\gamma - t_{i-1}^\gamma))}{n_i!}$$
$$= rac{\lambda^n}{n!} \exp(-\lambda t_r^\gamma) inom{n}{n_1 \dots n_r} \prod_{i=1}^r (t_i^\gamma - t_{i-1}^\gamma)^{n_i}$$
where $n = \sum_{i=1}^r n_i$ and $t_0 = 0$.

$N(t_r)$ is Poisson distributed with the mean value $\lambda t_r^\gamma$ and thus $N(t_r)$ is M-ancillary for $\gamma$. Inference about the trend parameter $\gamma$ should be drawn, in this case, from the conditional distribution of $N(t_1), N(t_2)-N(t_1), \dots, N(t_r)-N(t_{r-1})$ given $N(t_r)$. This distribution is given by the probabilities from the multinomial distribution
$$P(N(t_1) = n_1, \dots, N(t_r)-N(t_{r-1}) = n_r \mid N(t_r) = n) = inom{n}{n_1 \dots n_r} \prod_{i=1}^r \left(rac{t_i^\gamma - t_{i-1}^\gamma}{t_r^\gamma}ight)^{n_i}. 	ag{3*}$$

The maximum-likelihood estimate of $\gamma$ is in this case given implicitly, but it can naturally be determined with some difficulty. Asymptotically (for $n 	o \infty$) this estimate will be optimal, but no simple sufficient statistic for $\gamma$ is known, so there is no surety that the maximum-likelihood estimate has small-sample qualities which recommend its use. With this in mind, we will now concentrate on an alternative method of inference concerning $\gamma$, based on an explicit, asymptotically fully efficient estimator, which invites a simple test of fit of the model.

Note first that for $i = 1, \dots, r-1$ the definition of the non-homogeneous Poisson process has the direct result that $N(t_i)$ and $N(t_{i+1}) - N(t_i)$ are stochastically independent Poisson-distributed variables with the mean values $\lambda t_i^\gamma$ and $\lambda(t_{i+1}^\gamma - t_i^\gamma)$ respectively.

As a consequence, given $N(t_{i+1}) = n_{i+1}$, the conditional distribution of $N(t_i)$ is a binomial distribution with the probability $(t_i/t_{i+1})^\gamma$ and $n_{i+1}$ observations, which motivates the study of the $r-1$ estimates defined by
$$rac{N(t_i)}{N(t_{i+1})} pprox \left(rac{t_i}{t_{i+1}}ight)^\gamma \implies \hat{\gamma}_i = rac{\ln(N(t_i)/N(t_{i+1}))}{\ln(t_i/t_{i+1})}.$$

A test of fit of the R-W process is now available in the form of a comparison of the $r-1$ estimates. It seems that a natural choice for the final estimate of $\gamma$ is a weighted mean of the $\hat{\gamma}_i$ values.

Asymptotic properties of the conditional distribution of $(\hat{\gamma}_1, \dots, \hat{\gamma}_{r-1})$ given $N(t_r) = n$ are given in the following theorem.

**Theorem 2.** *Consider the column vector $V$, defined by $V' = (\sqrt{n}(\hat{\gamma}_1 - \gamma), \dots, \sqrt{n}(\hat{\gamma}_{r-1} - \gamma))$. For $n 	o \infty$, the asymptotic distribution of $V$ is $(r-1)$-variate normal with mean value 0 and dispersion matrix $D = (d_{ij})$, where $d_{ij} = 0$ for $i 
e j$ and*
$$d_{ii} = rac{(t_r/t_i)^\gamma - (t_r/t_{i+1})^\gamma}{[\ln(t_i/t_{i+1})]^2}.$$
*Thus the $r-1$ estimates are asymptotically independent.*

*Proof.* Theorem 2 is a consequence of the asymptotic normality of the multinomial distribution $(N(t_1), N(t_2)-N(t_1), \dots, N(t_r)-N(t_{r-1}) \mid N(t_r) = n)$ and a result of the transformation of the asymptotic

<!-- page 5 -->

normally-distributed variables (Rao, 1965, p. 322, (iii)).

An obvious choice is to use a weighted average of the $\hat{\gamma}_i$ values as an estimate of $\gamma$. The coefficients which minimize the asymptotic variance depend upon $D$, however, and thereby also upon the unknown $\gamma$. If it is suspected that the value of $\gamma$ is near 1, it will be reasonable to use coefficients which are optimal when $\gamma = 1$, and this gives the following estimate of $\gamma$:
$$	ilde{\gamma} = \sum_{i=1}^{r-1} k_i \hat{\gamma}_i, \quad 	ext{where } k_i = D_{ii}^{-1}(1) \Big/ \sum_{j=1}^{r-1} D_{jj}^{-1}(1).$$

The test of fit of the R-W process can be made by comparing the $\hat{\gamma}_i$ values with $	ilde{\gamma}$. This can be done, for example, with the aid of
$$H = n \sum_{i=1}^{r-1} D_{ii}^{-1}(	ilde{\gamma}) (\hat{\gamma}_i - 	ilde{\gamma})^2$$
which is asymptotically $\chi^2$-distributed with $r-1$ degrees of freedom. (See Rao, 1965, p. 323, where a similar result is proved.)

If the model is accepted, it is possible to improve the estimate by using weights determined from the preceding estimate. This gives a new weighted average

#### Table 1. The asymptotic efficiency of the estimate $	ilde{\gamma}$, for different values of $r$ and $\gamma$

| $r$ \ $\gamma$: | 0.50 | 0.70 | 0.90 | 0.95 | 1.00 | 1.05 | 1.10 | 1.20 | 1.50 | 2.00 |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0.965 | 0.987 | 0.999 | 1.000 | 1.000 | 1.000 | 0.999 | 0.994 | 0.964 | 0.862 |
| 6 | 0.953 | 0.983 | 0.998 | 1.000 | 1.000 | 1.000 | 0.998 | 0.992 | 0.952 | 0.820 |
| 7 | 0.943 | 0.979 | 0.998 | 0.999 | 1.000 | 0.999 | 0.998 | 0.990 | 0.942 | 0.787 |
| 8 | 0.935 | 0.976 | 0.997 | 0.999 | 1.000 | 0.999 | 0.997 | 0.989 | 0.933 | 0.759 |
| 9 | 0.927 | 0.973 | 0.997 | 0.999 | 1.000 | 0.999 | 0.997 | 0.988 | 0.926 | 0.736 |
| 10 | 0.921 | 0.971 | 0.997 | 0.999 | 1.000 | 0.999 | 0.997 | 0.987 | 0.920 | 0.716 |
| 15 | 0.898 | 0.962 | 0.996 | 0.999 | 1.000 | 0.999 | 0.996 | 0.983 | 0.897 | 0.646 |
| 20 | 0.883 | 0.956 | 0.995 | 0.999 | 1.000 | 0.999 | 0.995 | 0.980 | 0.881 | 0.602 |
| 30 | 0.863 | 0.948 | 0.994 | 0.999 | 1.000 | 0.999 | 0.994 | 0.977 | 0.862 | 0.547 |
| 100 | 0.818 | 0.931 | 0.992 | 0.998 | 1.000 | 0.998 | 0.992 | 0.969 | 0.817 | 0.423 |

#### Table 2. Numbers of operating hours between successive failures of airconditioning equipment in 13 aircraft

| Aircraft 1 | Aircraft 2 | Aircraft 3 | Aircraft 4 | Aircraft 5 | Aircraft 6 | Aircraft 7 | Aircraft 8 | Aircraft 9 | Aircraft 10 | Aircraft 11 | Aircraft 12 | Aircraft 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 194 | 413 | 90 | 74 | 55 | 23 | 97 | 50 | 359 | 50 | 130 | 487 | 102 |
| 15 | 14 | 10 | 57 | 320 | 261 | 51 | 44 | 9 | 254 | 493 | 18 | 209 |
| 41 | 58 | 60 | 48 | 56 | 87 | 11 | 102 | 12 | 5 |  | 100 | 14 |
| 29 | 37 | 186 | 29 | 104 | 7 | 4 | 72 | 270 | 283 |  | 7 | 57 |
| 33 | 100 | 61 | 502 | 220 | 120 | 141 | 22 | 603 | 35 |  | 98 | 54 |
| 181 | 65 | 49 | 12 | 239 | 14 | 18 | 39 | 3 | 12 |  | 5 | 32 |
|  | 9 | 14 | 70 | 47 | 62 | 142 | 3 | 104 |  |  | 85 | 67 |
|  | 169 | 24 | 21 | 246 | 47 | 68 | 15 | 2 |  |  | 91 | 59 |
|  | 447 | 56 | 29 | 176 | 225 | 77 | 197 | 438 |  |  | 43 | 134 |
|  | 184 | 20 | 386 | 182 | 71 | 80 | 188 |  |  |  | 230 | 152 |
|  | 36 | 79 | 59 | 33 | 246 | 1 | 79 |  |  |  | 3 | 27 |
|  | 201 | 84 | 27 | 15 | 21 | 16 | 88 |  |  |  | 130 | 14 |
|  | 118 | 44 | 153 | 104 | 42 | 106 | 46 |  |  |  |  | 230 |
|  | 34 | 59 | 26 | 35 | 20 | 206 | 5 |  |  |  |  | 66 |
|  | 31 | 29 | 326 |  | 5 | 82 | 5 |  |  |  |  | 61 |
|  | 18 | 118 |  |  | 12 | 54 | 36 |  |  |  |  | 34 |
|  | 18 | 25 |  |  | 120 | 31 | 22 |  |  |  |  |  |
|  | 67 | 156 |  |  | 11 | 216 | 139 |  |  |  |  |  |
|  | 57 | 310 |  |  | 3 | 46 | 210 |  |  |  |  |  |
|  | 62 | 76 |  |  | 14 | 111 | 97 |  |  |  |  |  |
|  | 7 | 26 |  |  | 71 | 39 | 30 |  |  |  |  |  |
|  | 22 | 44 |  |  | 11 | 63 | 23 |  |  |  |  |  |
|  | 34 | 23 |  |  | 14 | 18 | 13 |  |  |  |  |  |
|  |  | 62 |  |  | 11 | 191 | 14 |  |  |  |  |  |
|  |  | 130 |  |  | 16 | 18 |  |  |  |  |  |  |
|  |  | 208 |  |  |  | 90 |  |  |  |  |  | 163 |
|  |  | 70 |  |  |  | 1 |  |  |  |  |  | 24 |
|  |  | 101 |  |  |  | 16 |  |  |  |  |  |  |
|  |  | 208 |  |  |  | 52 |  |  |  |  |  |  |
|  |  |  |  |  |  | 95 |  |  |  |  |  |  |

<!-- page 6 -->

$$\hat{\gamma}^* = \sum_{i=1}^{r-1} k_i^* \hat{\gamma}_i \quad 	ext{where the weights are} \quad k_i^* = D_{ii}^{-1}(	ilde{\gamma}) \Big/ \sum_{j=1}^{r-1} D_{jj}^{-1}(	ilde{\gamma}).$$

Both estimates are asymptotically normal-distributed, for $n 	o \infty$, $	ilde{\gamma}$ with variance
$$rac{1}{n} \left(\sum_{i=1}^{r-1} D_{ii}^{-1}(1)ight)^{-2} \sum_{i=1}^{r-1} D_{ii}^{-2}(1) D_{ii}(\gamma)$$
and, as $D_{ii}(	ilde{\gamma}) 	o D_{ii}(\gamma)$, for $\hat{\gamma}^*$ with variance
$$rac{1}{n \sum_{i=1}^{r-1} D_{ii}^{-1}(\gamma)}.$$

By twice differentiating the probability (3*) with respect to $\gamma$ and then making use of
$$E[N(t_i) - N(t_{i-1}) \mid N(t_r) = n] = n rac{t_i^\gamma - t_{i-1}^\gamma}{t_r^\gamma}$$
it can be seen that the information in the conditional distribution of $N(t_1), \dots, N(t_{r-1})$ given $N(t_r)$ is equal to
$$n \sum_{i=1}^{r-1} D_{ii}^{-1}(\gamma).$$

The estimate $\hat{\gamma}^*$ is thus asymptotically fully efficient. The extent of loss of efficiency in the estimate $	ilde{\gamma}$ depends upon $r$, $\gamma$ and the values of $t_1, \dots, t_r$. This loss has been computed under the condition that $t_i = i	au$ for $i = 1, \dots, r$ and the results are presented in Table 1. It can be seen that the loss is negligible for $\gamma$ near 1 and moderate $r$.

#### Table 3. Values of the stochastic vectors corresponding to the data in Table 2

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.458 | 0.016 | 0.090 | 0.197 | 0.019 | 0.055 | 0.012 | 0.009 | 0.279 | 0.019 | 1.567 | 0.106 | 0.026 |
| 0.569 | 0.026 | 0.136 | 0.215 | 0.079 | 0.086 | 0.094 | 0.018 | 0.280 | 0.076 |  | 0.108 | 0.075 |
| 0.679 | 0.029 | 0.170 | 0.325 | 0.088 | 0.096 | 0.104 | 0.033 | 0.360 | 0.727 |  | 0.328 | 0.131 |
| 0.858 | 0.058 | 0.278 | 0.346 | 0.108 | 0.096 | 0.212 | 0.053 | 0.362 | 0.743 |  | 0.375 | 0.354 |
| 0.933 | 0.086 | 0.351 | 0.393 | 0.225 | 0.153 | 0.223 | 0.122 | 1.019 | 2.548 |  | 0.483 | 0.369 |
|  | 0.120 | 0.388 | 0.770 | 0.353 | 0.164 | 0.261 | 0.290 | 1.555 |  |  | 0.596 | 0.399 |
|  | 0.129 | 0.402 | 0.805 | 0.565 | 0.171 | 0.286 | 0.418 | 1.587 |  |  | 0.603 | 0.589 |
|  | 0.139 | 0.430 | 0.831 | 0.611 | 0.180 | 0.360 | 0.440 | 1.612 |  |  | 0.751 | 0.792 |
|  | 0.155 | 0.446 | 0.924 | 0.886 | 0.188 | 0.392 | 0.470 |  |  |  | 0.763 | 0.897 |
|  | 0.173 | 0.497 | 0.941 | 1.231 | 0.237 | 0.560 | 0.482 |  |  |  | 0.943 | 1.031 |
|  | 0.239 | 0.733 | 2.169 | 1.447 | 0.247 | 0.586 | 0.488 |  |  |  | 0.980 | 1.102 |
|  | 0.362 | 0.877 | 2.319 | 1.586 | 0.249 | 0.634 | 0.538 |  |  |  | 1.234 |  |
|  | 0.386 | 0.902 | 2.631 | 3.506 | 0.257 | 0.712 | 0.641 |  |  |  | 1.395 |  |
|  | 0.517 | 1.030 | 3.202 |  | 0.348 | 0.938 | 0.743 |  |  |  | 1.440 |  |
|  | 0.934 | 1.064 |  |  | 0.357 | 1.078 | 1.040 |  |  |  |  | 2.554 |
|  | 1.151 | 1.137 |  |  | 0.361 | 1.101 | 1.490 |  |  |  |  |  |
|  | 1.164 | 1.195 |  |  | 0.377 | 1.102 | 1.534 |  |  |  |  |  |
|  | 1.264 | 1.317 |  |  | 0.412 | 1.225 | 1.543 |  |  |  |  |  |
|  | 1.439 | 1.447 |  |  | 0.430 | 1.361 | 1.669 |  |  |  |  |  |
|  | 1.513 | 1.482 |  |  | 0.668 | 1.497 | 1.748 |  |  |  |  |  |
|  | 1.640 | 1.590 |  |  | 0.748 | 1.863 | 2.061 |  |  |  |  |  |
|  | 1.673 | 1.640 |  |  | 1.058 | 1.920 | 2.796 |  |  |  |  |  |
|  |  | 1.670 |  |  | 1.136 | 2.543 | 3.427 |  |  |  |  |  |
|  |  | 1.784 |  |  | 1.251 |  | 2.568 |  |  |  |  |  |
|  |  | 1.946 |  |  | 1.278 |  | 2.640 |  |  |  |  |  |
|  |  | 2.717 |  |  | 1.554 |  | 3.063 |  |  |  |  |  |
|  |  | 3.187 |  |  | 1.573 |  |  |  |  |  |  |  |
|  |  | 3.293 |  |  | 1.840 |  |  |  |  |  |  |  |
|  |  | 4.353 |  |  |  |  |  |  |  |  |  |  |
| **$\hat{\gamma}$:** | 0.699 | 0.601 | 1.150 | 1.148 | 0.823 | 0.687 | 1.051 | 0.959 | 0.882 | 0.823 | 1.567 | 0.549 | 0.826 |

## Examples

### Failure of air-conditioning equipment in aircraft

The first example uses data taken from "Numbers of operating hours between successive failures of air-conditioning equipment in 13 Boeing 720 aircraft". These data from Cox & Lewis (1966) are reprinted in Table 2.

The data were originally published by Proschan (1963) and it is not completely clear how the observations were made or, especially, how it was decided to conclude the observations. It is assumed,

<!-- page 7 -->

however, in Cox & Lewis, that the data have form (A); we will do likewise here.

Taking our origin in Theorem 1, each aircraft was individually analysed to determine whether or not the R-W process can be used as a satisfactory description of the data.

Table 3 presents the values of the stochastic vectors
$$(-\ln(W_{n-1}/W_n), -\ln(W_{n-2}/W_n), \dots, -\ln(W_1/W_n))$$
derived from the corresponding columns of Table 2. These values can be seen in relation to $n-1$ ordered observations from an exponential distribution by use of fractile diagrams where the $n-1$ points defined by
$$(-\ln(W_{n-m}/W_n), -\ln((n-m-1/2)/(n-1)))$$
are plotted for $m = 1, \dots, n-1$. These points should approximate a line with slope $eta^{-1} = \gamma$. The value of $eta$ is unknown, but it can be estimated by
$$\hat{eta} = rac{1}{n-1} \sum_{m=1}^{n-1} (-\ln(W_{n-m}/W_n)).$$

The line having the slope $eta^{-1} = \hat{\gamma}$ has been computed and plotted in the fractile diagrams.

The graphical control can be supplemented by a Kolmogorov-Smirnov test with the estimate of the scale parameter inserted, i.e. with value of

[FIGURE page 7: Fig. 1. Fractile diagram for airplane no. 3. Estimated slope: $\hat{eta} = 0.8697$.]

[FIGURE page 7: Fig. 2. Fractile diagram for airplane no. 6. Estimated slope: $\hat{eta} = 1.456$.]

[FIGURE page 7: Fig. 3. Fractile diagram for airplane no. 7. Estimated slope: $\hat{eta} = 0.9512$.]

[FIGURE page 7: Fig. 4. Fractile diagram for airplane no. 8. Estimated slope: $\hat{eta} = 1.043$.]

<!-- page 8 -->

$$D = \sup_{0 < x < \infty} |H(x) - P(x)|$$
where $H$ is the empirical distribution function
$$H(x) = egin{cases}
0 & 	ext{for } x < -\ln(W_{n-1}/W_n) \
m/(n-1) & 	ext{for } -\ln(W_{n-m+1}/W_n) \le x < -\ln(W_{n-m}/W_n) \
1 & 	ext{for } -\ln(W_1/W_n) \le x
\end{cases}$$
and $P$ is given by $P(x) = 1 - \exp(-\hat{\gamma} x)$. The critical values of $D$ at the levels of 5% and 1% are given by Lilliefors (1969).

Fractile diagrams for the cases of the four aircraft with the greatest numbers of observed defects are reproduced as Figs. 1 to 4, where the confidence limits corresponding to the 5% level of the K-S test are also plotted. All of the points in these four plots can be seen to be within the 5% confidence levels and the control values are also tabulated together with the corresponding 5% critical values in Table 4. In no case does a plot show a systematic variation from the expected values: not in any of these four figures nor, for that matter, in the fractile diagrams for the other nine aircraft. It thus seems reasonable to use the R-W process to describe the failures of the air-conditioning systems in the 13 aircraft.

Analysis could now be continued with likelihood-ratio tests of the hypotheses $\gamma_1 = \gamma_2 = \dots = \gamma_{13} = \gamma$ and $\gamma = 1$. Both of these hypotheses are acceptable, which is not surprising when it is recalled that Cox & Lewis (1966) used a non-homogeneous Poisson process in their analysis and concluded that these phenomena could be described with a homogeneous Poisson process.

#### Table 4. Results of the Kolmogorov-Smirnov tests for each of the 13 aircraft

| Aircraft no. | Number of observations | Test values | 5% critical values |
|---|---|---|---|
| 1 | 5 | 0.480 | 0.442* |
| 2 | 22 | 0.204 | 0.226 |
| 3 | 28 | 0.141 | 0.200 |
| 4 | 14 | 0.158 | 0.277 |
| 5 | 13 | 0.185 | 0.287 |
| 6 | 29 | 0.190 | 0.197 |
| 7 | 26 | 0.103 | 0.208 |
| 8 | 23 | 0.136 | 0.221 |
| 9 | 8 | 0.271 | 0.360 |
| 10 | 5 | 0.311 | 0.442 |
| 11 | 1 | — | — |
| 12 | 11 | 0.269 | 0.311 |
| 13 | 15 | 0.150 | 0.269 |

*\* Significant at 5% level. \*\* Significant at 1% level.*

### Oral reading

Data for this example is taken from tape recordings of 25 pupils reading aloud a series of reading tests

#### Table 5. Words read in ten-second intervals using Oral Reading Test T3

| No. | Numbers of ten-second periods | Number of words read in these periods (in chronological order) | $	ilde{\gamma}$ value | $\chi^2$-test value | Degrees of freedom |
|---|---|---|---|---|---|
| 1 | 14 | 21, 13, 8, 11, 17, 14, 13, 11, 17, 10, 9, 14, 14, 14 | 0.916 | 12.64 | 12 |
| 2 | 18 | 16, 7, 13, 7, 9, 18, 8, 14, 15, 9, 20, 3, 6, 8, 11, 13, 13 | 0.933 | 33.82** | 16 |
| 3 | 25 | 21, 12, 5, 2, 8, 3, 8, 9, 8, 6, 7, 15, 13, 10, 3, 8, 2, 11, 1, 5, 2, 5, 11, 14, 10 | 0.844 | 76.36*** | 23 |
| 4 | 10 | 19, 19, 11, 14, 23, 29, 17, 16, 21, 22 | 1.059 | 12.33 | 8 |
| 5 | 16 | 19, 17, 9, 5, 10, 13, 14, 12, 14, 13, 9, 14, 9, 8, 10, 16 | 0.879 | 24.85* | 14 |
| 6 | 8 | 32, 23, 28, 24, 18, 22, 29, 24 | 0.913 | 4.17 | 6 |
| 7 | 8 | 32, 22, 26, 28, 24, 26, 18, 24 | 0.901 | 3.26 | 6 |
| 8 | 12 | 22, 22, 19, 6, 23, 19, 19, 10, 18, 7, 11, 21 | 0.865 | 34.98*** | 10 |
| 9 | 21 | 7, 9, 9, 10, 6, 6, 4, 14, 17, 15, 14, 4, 10, 6, 12, 12, 3, 1, 11, 13, 13 | 1.068 | 66.50*** | 19 |
| 10 | 18 | 17, 12, 9, 10, 3, 12, 12, 16, 7, 11, 12, 13, 9, 10, 10, 7, 6, 15 | 0.910 | 35.27*** | 16 |
| 11 | 13 | 21, 17, 8, 14, 19, 18, 16, 10, 16, 16, 9, 16, 16 | 0.919 | 17.65 | 11 |
| 12 | 13 | 22, 16, 15, 17, 20, 18, 18, 10, 16, 9, 9, 7, 19 | 0.839 | 14.44 | 11 |
| 13 | 8 | 21, 16, 24, 33, 22, 21, 22, 26 | 1.027 | 3.53 | 6 |
| 14 | 9 | 21, 20, 14, 28, 29, 23, 24, 18, 23 | 1.049 | 7.70 | 7 |
| 15 | 6 | 29, 26, 27, 31, 34, 36 | 0.965 | 1.75 | 4 |
| 16 | 7 | 33, 26, 33, 26, 22, 31, 25 | 0.925 | 2.77 | 5 |
| 17 | 7 | 30, 19, 24, 24, 29, 27, 21 | 0.958 | 4.57 | 5 |
| 18 | 6 | 32, 29, 30, 34, 31, 31, 41 | 1.059 | 2.16 | 4 |
| 19 | 9 | 28, 18, 19, 17, 25, 22, 16, 23, 20 | 0.926 | 4.86 | 7 |
| 20 | 7 | 29, 22, 26, 28, 22, 22, 29 | 0.965 | 2.56 | 5 |
| 21 | 7 | 33, 22, 23, 31, 26, 22, 25 | 0.918 | 3.09 | 5 |
| 22 | 8 | 26, 24, 31, 26, 26, 18, 18, 21 | 0.911 | 4.08 | 6 |

*\* Significant at 5% level. \*\* Significant at 1% level. \*\*\* Significant at 0.5% level.*

<!-- page 9 -->

consisting of prose Danish-language texts of from 75 to 200 words. These tapes are used to determine the number of words read within each successive ten-second interval of recitation.

Table 5 contains a tabulation of word counts made with test T3, which has been selected from the seven possibilities because it is the test read by the greatest number of children.

Due to the fact that reading was stopped at a previously-determined number of words, the data is not of Type C, but of a mixture of types A and C. It is, however, a direct consequence of the characterization and (2) that the conditional distribution of $N(10), N(20) - N(10), \dots, N(10R) - N(10(R - 1))$, given $N(10R) = n$—where $R$ is the number of ten-second intervals completed when reading was stopped—is described by the multinomial distribution (3*) with $t_i = 10i$, i.e. the same as if $R$ had already been determined before the observations.

The suitability of the R-W process when used to describe the recitation-test results can therefore be controlled with a $\chi^2$ test based upon Theorem 2. The results of the estimation and the $\chi^2$ tests are also recorded in Table 5.

Of the 22 computed $\chi^2$ values, it can be seen that four are significant at the 0.5% level, one is significant at the 1% level, and one at the 5% level. Thus, it seems that the R-W process is not an acceptable model for describing the oral reading performances.

## IV. Discussion

It now seems natural to make a comparison between the Rasch-Weibull process and the non-homogeneous Poisson process with exponential linear trend as used by Cox & Lewis (1966), because either process can be used to describe situations where there is a monotonally increasing or decreasing intensity.

At the present time, most non-homogeneous processes present insurmountable problems as soon as an attempt is made to actually use the process to describe data, so that when such an attempt is to be made, it is desirable that the process used be as simple as possible. Either of the above-mentioned processes is a suitable possibility, as both are formally very simple, and both can be statistically analysed. Which of them to choose in a given circumstance must naturally depend upon the statement of the problem, but the R-W process possesses a number of obvious advantages which can reasonably be pointed out here.

Firstly, from the mathematical-statistical point of view, it is clearly an advantage to make an exact analysis without the use of asymptotic results. Whenever there is a complete set of waiting times, this can be accomplished with the R-W process and not with the other. When the record is incomplete, there is no principal difference between the two processes, as either may be used to make the analysis.

Secondly, and of greater importance in some cases, it is possible to give a clear characterization based upon the R-W process, which to a certain extent permits an evaluation of the relevancy of the model to the specific problem. A corresponding characterization is not possible in the case of the other process, so that it is difficult to evaluate whether or not the process does, in fact, function in a given situation.

It should be noted, however, that the non-homogeneous Poisson process with exponential linear trend does occur naturally in the theory of the pure birth process (e.g. Keiding, 1974).

In all fairness, a significant disadvantage in the use of the R-W process should be mentioned. If the characterization of the model is required together with a relatively simple statistical analysis, it is necessary that the temporal zero-point of the model be coincident with the starting point of the process or, at the least, that the point be known. In cases where this is not so, it seems in principle that the other process provides a more suitable basis for statistical analysis. Another possibility in this situation is to use the R-W process with the starting point as an unknown parameter, but the analysis becomes significantly more complicated.

## References

- Aczél, J. (1966). *Lectures on functional equations and their applications*. Academic Press, New York.
- Barlow, R. E. & Proschan, F. (1967). Exponential life test procedures when the distribution has monotone failure rate. *J. Am. Stat. Ass.* 62, 548-560.
- Barndorff-Nielsen, O. (1973). *Exponential families and conditioning*. Sc.D. thesis, University of Copenhagen.
- Cox, D. R. & Lewis, P. A. W. (1966). *The statistical analysis of series of events*. Methuen, London.
- Hollander, M. & Proschan, F. (1972). Testing whether new is better than used. *Ann. Math. Statist.* 43, 1136-1146.
- Keiding, N. (1974). Estimation of the birth process. *Biometrika* 61, 71-80.
- Kreiner Møller, S. (1970). *Statistical analysis of reading rates in oral reading* (in Danish). The Danish Institute for Educational Research, Copenhagen.
- Lilliefors, H. W. (1969). On the Kolmogorov-Smirnov test for the exponential distribution with mean unknown. *J. Am. Stat. Ass.* 64, 387-390.
- Lomnicki, Z. A. (1966). A note on the Weibull renewal process. *Biometrika* 53, 375-381.
- Parzen, E. (1962). *Stochastic processes*. Holden Day, San Francisco.
- Proschan, F. (1963). Theoretical explanation of observed decreasing failure rate. *Technometrics* 5, 375-383.
- Rao, C. Radhakrishna (1965). *Linear statistical inference and its applications*. Wiley, New York.

<!-- page 10 -->

- Rasch, G. (1970). Comments on report No. 12 par D. A. Srour (mimeographed).
- Sarhan, A. E. & Greenberg, B. G. (1962). Exponential distribution. In *Contributions to order statistics* (ed. A. E. Sarhan & B. G. Greenberg), pp. 343-382. Wiley, New York.
- Smith, W. L. & Leadbetter, M. R. (1963). On the renewal function for the Weibull distribution. *Technometrics* 5, 393-396.
- Stacy, E. W. (1962). A generalization of the gamma distribution. *Ann. Math. Statist.* 33, 1187-1192.

**Svend Kreiner Møller**  
The Danish Institute for Educational Research  
28, Hermodsgade  
DK-2200 Copenhagen N  
Denmark
