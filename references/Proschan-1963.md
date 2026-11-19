<!-- page 1 -->

# Theoretical Explanation of Observed Decreasing Failure Rate

**Author(s):** Frank Proschan  
**Source:** *Technometrics*, Feb., 2000, Vol. 42, No. 1, Special 40th Anniversary Issue (Feb., 2000), pp. 7-11  
*(Originally published in Technometrics, Aug., 1963, Vol. 5, No. 3, pp. 375-383)*  
**Published by:** Taylor & Francis, Ltd. on behalf of American Statistical Association and American Society for Quality  
**Stable URL:** https://www.jstor.org/stable/1271427  

<!-- page 2 -->

# Theoretical Explanation of Observed Decreasing Failure Rate

**Frank PROSCHAN\***  
*Boeing Scientific Research Laboratories, Seattle, Washington*  
*(8315 Stedmall Terrace, Gaithersburg, MD 20886)*  

> Pooled data on the times of successive failures of the air conditioning system of a fleet of jet airplanes seemed to indicate that the life distribution had a decreasing failure rate. More refined analysis showed that the failure distribution for each airplane separately was exponential, but with a different failure rate. Using the theorem that a mixture of distributions each having a non-increasing failure rate (e.g., a mixture of exponential distributions) itself has a non-increasing failure rate, the apparent decreasing failure rate of the pooled air-conditioning life distribution was satisfactorily explained. The present study has implications in other areas where an observed decreasing failure rate may well be the result of mixing exponential distributions having different parameters.

## 1. INTRODUCTION

It is easy to accept the assumption that the life of a unit has an increasing failure rate (IFR). This assumption is readily explained by wearout, so that the older the unit, the greater the chance that failure will occur in the next unit of time. By the failure rate $r(t)$ of a distribution $F$ with density $f$ corresponding to a non-negative random variable we mean $r(t) = f(t)/\bar{F}(t)$, where $\bar{F}(t) = 1 - F(t)$ throughout.

It is more difficult to explain why the life of a unit would have a decreasing failure rate (DFR). Ordinarily a decreasing failure rate would seem to correspond to some physical mechanism of improvement with age, so that the longer the unit survives the less the chance of failure in the next unit of time. However, there are some situations in which no process of improvement seems to be involved, and yet a decreasing failure rate is observed. (Throughout, we use "decreasing" for "non-increasing" and "increasing" for "non-decreasing.")

In this paper we show theoretically why a life distribution would under certain circumstances exhibit a decreasing failure rate. We then apply a new test for monotonicity of failure rate to some observed life test data to verify that the life distribution does exhibit the decreasing failure rate expected under the circumstances.

The data arose in the following way. Records were kept for the time of successive failures of the air conditioning system of each member of a fleet of Boeing 720 jet airplanes. The intervals between successive failures are listed in order of occurrence in Table 1. After roughly 2000 hours of service the planes received major overhaul; the failure interval containing major overhaul is omitted from the listing since the length of that failure interval may have been affected by major overhaul.

The purpose of the data collection was, of course, to obtain information as to the distribution of failure intervals for the air conditioning system. Obvious applications of such information would be in predicting reliability, scheduling maintenance, and providing spare parts.

## 2. EVIDENCE AGAINST EXPONENTIAL FIT OF POOLED DATA

Since the original aim was to make predictions and decisions for the entire fleet of 720 jets (rather than for individual planes), in the original preliminary analysis:
1. failure intervals for the different jets were pooled together and a single empirical fraction surviving histogram was plotted as a function of time (see Fig. 1).
2. the mean $\bar{x}$ of all the failure intervals was computed, yielding $\bar{x} = 93.14$ hours.
3. An exponential survival probability $\bar{F}(t)$ was fitted to the pooled failure intervals:
$$\bar{F}(t) = e^{-t/93.14} \tag{1}$$
4. Finally the Kolmogorov-Smirnov test was performed to determine whether the exponential distribution offered a good fit of the pooled failure intervals, Kendall and Stuart, 1961, pp. 452-461. Under the hypothesis of exponential fit, the probability is approximately .95 that the empirical fraction surviving does not deviate from (1) by more than $1.3581 n^{-1/2}$ for large $n$; in the present case the upper and lower limits for the fraction surviving are obtained by adding and subtracting $1.3581(213)^{-1/2} = .0931$ from (1), as displayed in Fig. 1. Note that the empirical fraction surviving is completely contained between the upper and lower limits. (It should be pointed out that in the Kolmogorov-Smirnov test no provision is made for estimating a parameter from the data, i.e., the test compares an empirical distribution with a known (not estimated) distribution. However, since the preliminary analysis is based on 213 observations, the discrepancy between the estimate and the true value of the exponential mean should be quite small.)

Actually, however, the exponential distribution is not a good fit of the pooled failure intervals even though the Kolmogorov-Smirnov test does not reject the hypothesis of exponential fit. This is shown by the fact that the empirical

<!-- page 3 -->

fraction surviving lies consistently below the corresponding exponential for $0 < t < 150$ and then lies consistently above the exponential for $t > 150$. The fact that the survival probability curve seems to cross the exponential once only, and from below, suggests strongly the possibility that the failure distribution for the pooled failure intervals has a decreasing failure rate, as supported by

### Table 1. Intervals Between Failures

| Plane: 7907 | 7908 | 7909 | 7910 | 7911 | 7912 | 7913 | 7914 | 7915 | 7916 | 7917 | 8044 | 8045 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 194 | 413 | 90 | 74 | 55 | 23 | 97 | 50 | 359 | 50 | 130 | 487 | 102 |
| 15 | 14 | 10 | 57 | 320 | 261 | 51 | 44 | 9 | 254 | 493 | 18 | 209 |
| 41 | 58 | 60 | 48 | 56 | 87 | 11 | 102 | 12 | 5 | | 100 | 14 |
| 29 | 37 | 186 | 29 | 104 | 7 | 4 | 72 | 270 | 283 | | 7 | 57 |
| 33 | 100 | 61 | 502 | 220 | 120 | 141 | 22 | 603 | 35 | | 98 | 54 |
| 181 | 65 | 49 | 12 | 239 | 14 | 18 | 39 | 3 | 12 | | 5 | 32 |
| | 9 | 14 | 70 | 47 | 62 | 142 | 3 | 104 | | | 85 | 67 |
| | 169 | 24 | 21 | 246 | 47 | 68 | 15 | 2 | | | 91 | 59 |
| | 447 | 56 | 29 | 176 | 225 | 77 | 197 | 438 | | | 43 | 134 |
| | 184 | 20 | 386 | 182 | 71 | 80 | 188 | | | | 230 | 152 |
| | 36 | 79 | 59 | 33 | 246 | 1 | 79 | | | | 3 | 27 |
| | 201 | 84 | 27 | 15 | 21 | 16 | 88 | | | | 130 | 14 |
| | 118 | 44 | ** | 104 | 42 | 106 | 46 | | | | | 230 |
| | ** | ** | 15 | 35 | 20 | 206 | 5 | | | | | 66 |
| | 34 | 59 | 153 | ** | 5 | 82 | 5 | | | | | 61 |
| | 31 | 29 | 26 | 12 | 12 | 54 | 36 | | | | | 34 |
| | 18 | 118 | 326 | 25 | 120 | 31 | 22 | | | | | |
| | 18 | 25 | | | 11 | 216 | 139 | | | | | |
| | 67 | 156 | | | 3 | 46 | 210 | | | | | |
| | 57 | 310 | | | 14 | 111 | 97 | | | | | |
| | 62 | 76 | | | 71 | 39 | 30 | | | | | |
| | 7 | 26 | | | 11 | 63 | 23 | | | | | |
| | 22 | 44 | | | 14 | 18 | 13 | | | | | |
| | 34 | 23 | | | 11 | 191 | 14 | | | | | |
| | | 62 | | | 16 | 18 | | | | | | |
| | | 130 | | | | 90 | | | | | | 163 |
| | | 208 | | | | 1 | | | | | | 24 |
| | | 70 | | | | 16 | | | | | | |
| | | 101 | | | | 52 | | | | | | |
| | | 208 | | | | 95 | | | | | | |

*\*\* Major overhaul.*

[FIGURE page 3: Figure 1. Empirical Fraction Surviving Histogram vs. Fitted Exponential Survival Probability Curve and 95% Kolmogorov-Smirnov Limits.]

**Theorem 1.** *Let $F$ be a distribution with density $f$, satisfying (a) $F(t) = 0$ for $t \le 0$, (b) the failure rate $r(t) = f(t)/\bar{F}(t)$ is decreasing for $t > 0$, (c) $\bar{F}(\xi_p) = p$; i.e., $\xi_p$ is the $p$th percentile. Then*
$$\bar{F}(t) \begin{cases} < e^{-at} & \text{for } t < \xi_p \\ > e^{-at} & \text{for } t > \xi_p \end{cases}$$

<!-- page 4 -->

where $a = -\log(1 - p)/\xi_p$. (See Barlow, Marshall, Proschan (1963).)

*Proof.* The result follows immediately from the following lemmas:

**Lemma 1.** *Under assumptions (a) and (b) of the theorem $\log \bar{F}(t)$ is convex for $0 < t < \infty$.*

*Proof.* The derivative of $\log \bar{F}(t)$ is $-f(t)/\bar{F}(t)$, an increasing function for $0 < t < \infty$ by hypothesis. $\parallel$

**Lemma 2.** *Under assumptions (a) and (b), $\{\bar{F}(t)\}^{1/t}$ is increasing in $t > 0$.*

*Proof.* By Lemma 1
$$\frac{\log \bar{F}(t) - \log \bar{F}(0)}{t - 0}$$
is increasing in $t > 0$. Hence $\{\bar{F}(t)\}^{1/t}$ is increasing in $t > 0$. $\parallel$

In Section 5 we shall see why the hypothesis of a decreasing failure rate for the distribution of pooled failure intervals is a natural one to consider on theoretical grounds. Then we shall perform a statistical test on the observed data to decide whether the distribution of pooled failure intervals has a constant failure rate (exponential distribution) or a non-constant decreasing failure rate. Before we analyze the pooled data for the entire fleet, it would be wise, however, to obtain a clearer picture concerning the behavior of successive failure intervals for individual airplanes.

## 3. DO SUCCESSIVE INTERVALS BETWEEN FAILURES SHOW A TREND?

For a given airplane it might well be that successive failure intervals would tend to get either longer or shorter with increased operation of the plane. A trend toward longer intervals, if established, could be the result of greater experience, debugging, or elimination of faulty parts, whereas a trend toward shorter intervals could be the result of wearout, aging, or poor maintenance.

To determine whether successive intervals show a trend, a nonparametric test against trend developed by Mann (1945) is used. Let $X_1, X_2, \dots, X_n$ represent the successive intervals between failures for a given plane. Let $T_n$ represent the number of distinct comparisons between a pair of intervals for which the earlier interval is less than the later interval (the number of inequalities $X_i < X_j, i < j$). If two intervals are actually equal, a contribution of $\frac{1}{2}$ is made toward $T_n$. Thus for plane 7907, $T_6 = 8$ since among the 15 comparisons possible, in 8 cases an earlier interval is less than a later interval.

Table 2 gives the detailed results. Note that the intervals between failures preceding major overhaul are treated separately from the intervals following major overhaul (the latter intervals are indicated by "MO" after the plane number). The second column gives $n$, the number of intervals in question. The fourth column gives the expected value $E T_n$ of $T_n$ under $H_0$, the hypothesis of no trend, shown by Mann (1945) to be $n(n - 1)/4$. The fifth column gives the variance $\sigma^2(T_n)$ of $T_n$ under $H_0$, shown by Mann to be $(2n + 5)(n - 1)n/72$. Mann also shows that asymptotically under $H_0$, $T_n$ follows a normal distribution. Thus the standardized deviate $(T_n + \frac{1}{2} - E T_n)/\sigma(T_n)$ shown in column 6 is approximately normally distributed with 0 mean, unit variance, under $H_0$. (The addition of "$\frac{1}{2}$" is required because of the discrete nature of the data.) The probability of getting a value of $T_n$ larger than the observed value under $H_0$ is shown in column 7. For $n \ge 10$ the normal approximation is used, while for $n < 10$ the exact probability is tabulated in Mann (1945). Mann shows that the normal approximation is quite good at $n = 10$. Finally in column 8 the logarithm to base 10 of this probability is listed.

To combine the outcomes of the individual tests of significance corresponding to each airplane, the Fisher (1958)

### Table 2. Mann Test Against Trend For Successive Intervals Between Failures

| Plane | $n$ | $T_n$ | $E T_n$ | $\sigma^2(T_n)$ | $\frac{T_n + \frac{1}{2} - E T_n}{\sigma(T_n)}$ | $P[T_n \ge \text{obs. } T_n]$ | $\log_{10} P$ |
|---|---|---|---|---|---|---|---|
| 7907 | 6 | 8 | 7.5 | 7.08 | | .5000 | .69897-1 |
| 7908 | 13 | 47 | 39 | 67.17 | 1.04 | .1492 | .17377-1 |
| 7908MO | 10 | 22 | 22.5 | 33.75 | | .5000 | .69897-1 |
| 7909 | 24 | 142.5 | 138 | 406.33 | .22 | .4129 | .61584-1 |
| 7909MO | 5 | 5.5 | 5 | 4.17 | | .3250 | .51188-1 |
| 7910 | 12 | 27.5 | 33 | 53.17 | -.69 | .7549 | .87789-1 |
| 7910MO | 3 | 2 | 1.5 | .92 | | .1670 | .22272-1 |
| 7911 | 11 | 26 | 27.5 | 41.25 | -.16 | .5636 | .75097-1 |
| 7911MO | 3 | 2 | 1.5 | .92 | | .1670 | .22272-1 |
| 7912 | 30 | 173.5 | 217.5 | 785.42 | -1.55 | .9394 | .97285-1 |
| 7913 | 27 | 177.5 | 175.5 | 575.25 | .10 | .4602 | .66295-1 |
| 7914 | 24 | 120 | 138 | 406.33 | -.87 | .8078 | .90730-1 |
| 7915 | 9 | 17 | 18 | 23.00 | | .5400 | .73239-1 |
| 7916 | 6 | 6 | 7.5 | 7.08 | | .6400 | .80618-1 |
| 7917 | 2 | 1 | .5 | .25 | | .5000 | .69897-1 |
| 8044 | 12 | 30 | 33 | 53.17 | -.34 | .6331 | .80147-1 |
| 8045 | 16 | 58.5 | 60 | 123.33 | -.090 | .5359 | .72908-1 |
| **Total** | | | | **2888.59** | | | **-5.91508** |

*NOTE: $T_n = \text{Number of inequalities } X_i < X_j, i < j, + \text{ half the number of equalities } X_i = X_j, i < j$. Expected value $E T_n = n(n-1)/4$, variance $\sigma^2(T_n) = (2n+5)(n-1)n/72$. (Note: actual sum of column 5 is 2588.59).*

<!-- page 5 -->

procedure is used. This procedure is based on the fact that under the null hypothesis each $\log_{10} P$ of column 8 after multiplication by -4.605 follows the $\chi^2$ distribution with 2 degrees of freedom, so that $-4.605 \sum \log_{10} P$, the summation taken over the 17 entries, follows the $\chi^2$ distribution with 34 degrees of freedom. The computed value of $-4.605 \sum \log_{10} P = 27.239$ corresponds to a probability of .21 in the left hand tail of the $\chi^2$ distribution, or about .42 for the two tails. Thus no significant evidence of a trend among successive intervals exists.

A corroboration of this conclusion may be obtained as follows. Compute $U = \sum (T_n + \frac{1}{2} - E T_n) / \sqrt{\sum \sigma^2(T_n)}$, the summations being taken over all 17 entries. $U$ should be approximately normally distributed even though some of the entries summed correspond to quite small values of $n$. The resulting $U = -1.00$ corresponds to a probability of .16 in the left hand tail of the normal distribution or .33 for the two tails. Again no significant evidence of a trend among successive intervals exists.

## 4. DO THE INTERVALS BETWEEN FAILURES FOR A SINGLE AIRPLANE FOLLOW THE EXPONENTIAL DISTRIBUTION?

From the analysis of Section 3 we see that it would be appropriate to consider the successive intervals between failures for a single airplane to be governed by a single probability distribution. To test whether the distribution is exponential we construct Table 3. As before column 2 lists $n$, the number of observations for the interval between successive failures. To obtain $V_n$, the entry in column 3, we must carry out the following steps:

1. For a single airplane rank the intervals between failures according to size; let the resulting ranked intervals be $Y_1 \le Y_2 \le \dots \le Y_n$. Thus for plane 7907 the ranked values are 15, 29, 33, 41, 181, 194 (see Table 1).
2. Compute the successive spacings $D_i' = Y_i - Y_{i-1}$, where $Y_0 = 0$ by definition. Thus $D_1' = 15 - 0 = 15$, $D_2' = 29 - 15 = 14$, $D_3' = 4$, $D_4' = 8$, $D_5' = 140$, $D_6' = 13$ for plane 7907.
3. Compute the "normalized" spacings $D_i = (n - i + 1)D_i'$. For plane 7907 $D_1 = 6(15) = 90$, $D_2 = 5(14) = 70$, $D_3 = 4(4) = 16$, $D_4 = 3(8) = 24$, $D_5 = 2(140) = 280$, $D_6 = 1(13) = 13$.
4. For each plane $V_n$ is computed as the number of pairs $i < j$ for which $D_i > D_j$; in case of equality count $\frac{1}{2}$ towards $V_n$. Thus $V_6 = 10$ for plane 7907.

As pointed out in Proschan [and] Pyke, under the hypothesis $H_0$ that the distribution for the interval between successive failures is exponential the expected value $E V_n$ for $V_n$ is $n(n - 1)/4$ while the variance $\sigma^2(V_n)$ for $V_n$ is $(2n + 5)(n - 1)n/72$. $E V_n$ is listed in column 4, $\sigma^2(V_n)$ in column 5 of Table 3. Moreover it is also shown that $V_n$ is asymptotically normally distributed under $H_0$. A natural alternative hypothesis $H_1$ to test against would be that the failure rate is increasing, rather than constant as under $H_0$. Under $H_1$ the successive normalized spacings $D_i$ would tend to get smaller and smaller, so that $V_n$ would tend to be larger under $H_1$ than under $H_0$. Thus the test developed in Proschan, Pyke consists in rejecting $H_0$ in favor of $H_1$ if $V_n$ is larger than a critical value depending on $n$.

### Table 3. To Test Whether Intervals Between Failures for a Single Airplane Have an Exponential Distribution

| Plane | $n$ | $V_n$ | $E V_n$ | $\sigma^2(V_n)$ |
|---|---|---|---|---|
| 7907 | 6 | 10 | 7.5 | 7.08 |
| 7908 | 13 | 38.5 | 39 | 67.17 |
| 7908MO | 10 | 30.5 | 22.5 | 33.75 |
| 7909 | 24 | 142 | 138 | 406.33 |
| 7909MO | 5 | 8 | 5 | 4.17 |
| 7910 | 12 | 35 | 33 | 53.17 |
| 7910MO | 3 | 1 | 1.5 | .92 |
| 7911 | 11 | 34 | 27.5 | 41.25 |
| 7911MO | 3 | 1 | 1.5 | .92 |
| 7912 | 30 | 193 | 217.5 | 785.42 |
| 7913 | 27 | 197 | 175.5 | 575.25 |
| 7914 | 24 | 134.5 | 138 | 406.33 |
| 7915 | 9 | 10.5 | 18 | 23.00 |
| 7916 | 6 | 6 | 7.5 | 7.08 |
| 7917 | 2 | 0 | .5 | .25 |
| 8044 | 12 | 25.5 | 33 | 53.17 |
| 8045 | 16 | 67 | 60 | 123.33 |
| **Total** | | **933.5** | **925.5** | **2588.59** |

*NOTE: $V_n = \text{Number of inequalities } D_i > D_j, i < j, + \text{ half the number of equalities } D_i = D_j, i < j$. $E V_n = n(n-1)/4, \sigma^2(V_n) = (2n+5)(n-1)n/72$.*

Now the amount of data available for each individual airplane is relatively limited, while generally to test whether a distribution follows a certain form requires large amounts of data. To overcome this difficulty we pool the data available for the different airplanes without necessarily assuming that the failure intervals for different planes have identical distributions; rather all we need test under $H_0$ is that each plane has a constant failure rate (equivalent to the assumption of an exponential distribution), the constant being different for the different airplanes. We thus compute
$$\frac{\sum (V_n + \frac{1}{2} - E V_n)}{\sqrt{\sum \sigma^2(V_n)}},$$
the summation being taken over all 17 entries. We obtain
$$\frac{933.5 + 8.5 - 925.5}{\sqrt{2588.59}} = .32.$$

Using the fact that under $H_0$ the ratio is asymptotically normally distributed, we find that the probability of observing a value as large as $.32$ is $.25$, so that no conclusive evidence exists that the intervals between failures for the individual planes have increasing failure rates rather than constant failure rates. Putting it more positively, it seems safe to accept the exponential distribution as describing the failure interval, although to each plane may correspond a different failure rate.

<!-- page 6 -->

## 5. DO THE DIFFERENT PLANES HAVE DIFFERING FAILURE RATES?

If all the planes in question had the same failure rate, then the failure intervals pooled together for the different planes would be governed by a single exponential distribution. (We call this hypothesis $H_0$.) On the other hand if to each plane there corresponded a different failure rate, it would then follow that the distribution of the pooled failure intervals would have a decreasing failure rate. (We call this italicized hypothesis $H_1$.) To validate this inference, we invoke the following theorem from Barlow, Marshall, Proschan (1963):

**Theorem 2.** *If $F_i(t)$ has a decreasing failure rate, $i = 1, 2, \dots, n$, then $G(t) = \sum_{i=1}^n p_i F_i(t)$ has a decreasing failure rate, where each*
$$p_i > 0, \quad \sum_{i=1}^n p_i = 1.$$

*Proof.* First suppose that $F_i(t)$ has a differentiable density $f_i(t), i = 1, 2, \dots, n$. Since the density of any DFR distribution must be a decreasing function, we have by Schwartz' inequality that
$$\sum_{i=1}^n p_i \bar{F}_i(t) \sum_{i=1}^n \{-p_i f_i'(t)\} \ge \left\{\sum_{i=1}^n p_i [-\bar{F}_i(t) f_i'(t)]^{1/2}\right\}^2.$$
Since $f_i(t)/\bar{F}_i(t)$ is decreasing in $t$, we must have
$$\bar{F}_i(t) f_i'(t) \le -\{f_i(t)\}^2.$$
Hence
$$\sum_{i=1}^n p_i \bar{F}_i(t) \sum_{i=1}^n \{-p_i f_i'(t)\} \ge \left\{\sum_{i=1}^n p_i f_i(t)\right\}^2,$$
that is,
$$\bar{G}(t) g'(t) \le -[g(t)]^2,$$
so that $G$ is DFR.

If the $F_i$ do not have differentiable densities, then the same result may be obtained by limiting arguments. $\parallel$

Now a special case of a distribution with decreasing (i.e., non-increasing) failure rate is the exponential distribution. It follows from Theorem 2 that convex combinations of exponentials have decreasing failure rate. Since such a convex combination cannot itself be an exponential, it must have a non-constant decreasing failure rate.

In the present situation, pooling of failure intervals following different exponential distributions for the different planes corresponds approximately to forming convex combinations of these exponential distributions. This would explain the choice of $H_1$. We say that pooling of failure intervals only approximately corresponds to forming convex combinations of exponential distributions for the following reason. The number of failure intervals obtained for the $i$th plane before major overhaul at 2000 hours of service will be approximately $2000 \lambda_i$, where $\lambda_i$ is the failure rate of the exponential life distribution for the $i$th airplane. Thus in Theorem 2 we may take $p_i = 2000 \lambda_i / \sum 2000 \lambda_j$. We can claim that the pooling corresponds only approximately to forming convex combination; since the actual number of failure intervals observed for a given plane is a random variable which depends on the actual lengths of the observed failure intervals. However by the fundamental renewal theorem this random variable should not deviate much from the constant $p_i$ specified above.

To perform the test of $H_0$ vs. $H_1$, we rank all the failure intervals for the different airplanes into one set $Y_1 \le Y_2 \le \dots \le Y_{213}$. We then compute as before $D_i' = Y_i - Y_{i-1}$, $D_i = (n - i + 1)D_i'$ for $i = 1, 2, \dots, 213$. Finally we compute the number $V_{213}$ of distinct pairs in which a $D_i > D_j, i < j$. Under the hypothesis $H_1$ of a decreasing failure rate we would expect smaller values of $V_{213}$ than under the hypothesis $H_0$ of a constant failure rate. (See Proschan and Pyke.) In fact under $H_0$, $V_{213}$ is asymptotically normally distributed with mean $E V_{213} = (213 \times 212)/4 = 11,289$ and standard deviation
$$\sigma(V_{213}) = [431(212)(213)/72]^{1/2} = 519.9.$$

We readily obtain from normal probability tables a probability of only .007 of getting a value less than the observed
$$\frac{V_{213} + \frac{1}{2} - E V_{213}}{\sigma(V_{213})} = -2.46.$$
Thus we must reject $H_0$ in favor of $H_1$. That is, we conclude that the pooled distribution has a decreasing failure rate, as would be expected if the individual airplanes each displayed a different constant failure rate.

## 6. SIMILAR SITUATIONS IN WHICH DFR WOULD ARISE

There seems to be a number of commonly occurring situations in which life distributions with decreasing failure rate would arise in a similar fashion. For example, semiconductors seem to exhibit a DFR life distribution. One reasonable explanation might be that the units within a given production lot exhibit a constant failure rate since semiconductors do not seem to wear out, but that the failure rate varies from lot to lot as a result of inescapable manufacturing variability. On the basis of Theorem 2, we would conclude that the life distribution of units for the various lots combined would have a DFR.

More generally, it has been observed that many types of electronic components have an exponential life distribution (Davis 1952). It would seem reasonable that the distribution formed by pooling failure intervals for a given component operating in different types of equipment would have a DFR.

## REFERENCES

- Barlow, R. E., Marshall, A. W., and Proschan, F. (1963), "Properties of Probability Distributions With Monotone Hazard Rate," *The Annals of Mathematical Statistics*, 34, 375-389.
- Barlow, R. E., and Proschan, F. (in press), *Probabilistic Models in Reliability Theory*, New York: Wiley.
- Davis, D. J. (1952), "The Analysis of Some Failure Data," *Journal of the American Statistical Association*, 47, 113-150.
- Fisher, R. W. (1958), *Statistical Methods for Research Workers* (13th ed.), London: Oliver & Boyd.
- Kendall, M. G., and Stuart, A. (1961), *The Advanced Theory of Statistics* (Vol. 2), New York: Hafner.
- Mann, H. B. (1945), "Nonparametric Tests Against Trend," *Econometrica*, 13, 245-249.
- Proschan, F., and Pyke, R. (in press), "Tests for Increasing Failure Rate," report, Boeing Scientific Research Laboratories.
