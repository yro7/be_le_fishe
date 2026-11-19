<!-- page 1 -->
# Mathematical Modeling and AI — 5th year
## Poisson processes and application to reliability theory and actuarial science
Academic year 2026-2027
**Mélisande ALBERT**
Office 115 GMM
`melisande.albert@insa-toulouse.fr`

<!-- page 2 -->

<!-- page 3 -->
# Contents

- **Introduction** (p. 5)
- **Bibliography** (p. 7)
- **1. Some reminders in probability** (p. 9)
  - I. Poisson random variable (p. 9)
  - II. Exponential random variable (p. 11)
  - III. Memoryless distribution (p. 13)
  - IV. Hazard rate in lifetime models (p. 15)
  - V. On the gamma and Erlang's distributions (p. 18)
- **2. Homogeneous Poisson processes** (p. 21)
  - I. Definitions (p. 21)
    - I.1. Point processes and counting processes (p. 21)
    - I.2. Homogeneous Poisson process on $\mathbb{R}_+$ (p. 25)
  - II. Arrival and interarrival times (p. 29)
    - II.1. First arrival time after a given instant (p. 29)
    - II.2. Exponential interarrival times (p. 29)
    - II.3. Conditional distribution of the arrival times (p. 32)
    - II.4. Construction of a homogeneous Poisson process (p. 35)
  - III. Divisibility of homogeneous Poisson processes (p. 38)
    - III.1. Two states divisibility (p. 38)
    - III.2. Infinite divisibility (p. 39)
- **3. Statistics for homogeneous Poisson processes** (p. 41)
  - I. Introduction (p. 41)
  - II. Fixed window, random number of points (p. 41)
    - II.1. Maximum Likelihood Estimator (p. 42)
    - II.2. Non-asymptotic properties of the MLE (p. 43)
    - II.3. Asymptotic distribution of the MLE (p. 45)
    - II.4. Confidence intervals and statistical tests (p. 47)
  - III. Fixed number of points, random observation time (p. 49)
    - III.1. Maximum Likelihood Estimator (MLE) (p. 49)
    - III.2. Non-asymptotic properties of the MLE (p. 50)
    - III.3. Asymptotic distribution of the MLE (p. 51)
    - III.4. Confidence intervals and statistical tests (p. 52)

<!-- page 4 -->
- **4. Inhomogeneous Poisson processes** (p. 57)
  - I. Definitions (p. 57)
  - II. Construction of an inhomogeneous Poisson process (p. 60)
    - II.1. Time change (p. 60)
    - II.2. Acceptance/Rejection (or Thinning) (p. 64)
  - III. Arrival times (p. 66)
    - III.1. First arrival time after a given instant (p. 66)
    - III.2. (Conditional) distribution of the arrival times (p. 67)
  - IV. Parametric statistics for inhomogeneous Poisson processes (p. 69)
- **A. Mathematical tools and complementary proofs** (p. 71)
  - I. Integration by substitution (p. 71)
  - II. Generalized inverse function (p. 71)
  - III. Complementary proofs (p. 73)
    - III.1. Proof of Lemma 1.13 (p. 73)
    - III.2. Proof of Proposition 4.13 (p. 75)
- **B. Tutorials and Computer Lab** (p. 77)
  - I. Worksheet 1 (p. 78)
  - II. Worksheet 2 (p. 80)
  - III. Worksheet 3 (p. 82)
  - IV. Worksheet 4 (p. 83)
  - V. Computer Lab (p. 85)
  - VI. Last year Exam (p. 88)

<!-- page 5 -->
# Introduction

The Poisson process, named after the French mathematician Siméon Denis Poisson, lies among the most important stochastic processes in probability theory. It is widely used to model random "points" in time and space, such as radioactive emission times, the arrival times of customers at a service center, earthquakes, eye movements or the positions of flaws in a piece of material.

The Poisson process is used as a foundation for building a number of other, more complicated random processes (e.g. the compound Poisson process). The analysis of such processes differs from usual real or vectorial frameworks and has a particular mathematical formalism.

```
Usual vectorial framework (R+)^n case: fixed number of observations
0 -------×---×-×-------------×------------------------------------> t

Point process framework: random number of observations
0 -------------×-----×-----×-×-----------×-----------···----------> t
```

The aim of this lecture is to introduce such processes in time and to study their basic properties from both probabilistic and statistical aspects. Application to reliability theory and to actuarial sciences will be done in some examples.

- **Application to reliability theory.** We consider a system subject to recurrent failures. We assume that when the system fails, it is repaired to a functioning state, and that the repair time is negligible. Different questions arise:
  - How can we model the successive failure times? How are they distributed in time?
  - Are the failure occurrences stable?
  - Can we estimate the rate of occurence failure per unit of time?
  - ...

- **Application to actuarial sciences.** We aim at describing and modeling insurer's vulnerability to ruin. Claims happen at random times, with random sizes.
  - How are distributed the claim times? and the amounts?
  - How can we model the total claim amount or the risk of failure?
  - Can we estimate the ruin probability?
  - How should the insurance company calibrate the premium?
  - ...

<!-- page 6 -->
### Abbreviations

| Abbreviation | Meaning |
| :--- | :--- |
| **a.s.** | almost surely |
| **c.d.f** | cumulative distribution function $F_X : t \in \mathbb{R} \mapsto \mathbb{P}(X \le t)$ |
| **CLT** | Central Limit Theorem |
| **i.i.d.** | independent and identically distributed |
| **LLN** | Law of Large Numbers |
| **m.g.f.** | moment-generating function $L_X : u \mapsto \mathbb{E}[e^{uX}]$ |
| **MLE** | Maximum Likelihood Estimator |
| **r.v.** | random variable |
| **s.t.** | such that |
| **w.r.t** | with respect to |

<!-- page 7 -->
# Bibliography

- **[CT97]** Christiane Cocozza-Thivent. *Processus stochastiques et fiabilité des systèmes*, volume 28. Springer Science & Business Media, 1997.
- **[FF04]** Dominique Foa and Aimé Fuchs. *Processus stochastiques, processus de Poisson, chaînes de Markov et martingales*. Dunod, Paris, 2004.
- **[Mik09]** Thomas Mikosch. *Non-Life Insurance Mathematics: An Introduction with the Poisson Process*. Springer-Verlag Berlin Heidelberg, 2nd edition, 2009.
- **[Ros06]** Sheldon M Ross. *Introduction to Probability Models*. Academic press, 9th edition, 2006.

<!-- page 8 -->

<!-- page 9 -->
# Chapter 1. Some reminders in probability

Several important probability distributions arise naturally from the Poisson process: the Poisson distribution, the exponential distribution, and the gamma distribution. This chapter aims to recall basic definitions and properties in probability theory necessary to define and study Poisson processes.

## I. Poisson random variable

The Poisson distribution is a discrete probability distribution commonly used to model the number of (rare) events occurring in a fixed interval of time or space. For instance, an insurer might want to model the number of claims during a month, or in quality control, one might want to model the number of failures of a manufacturing machine during a week.

> **Definition 1.1.**
> Let $N$ be a real random variable. Then $N$ is said to have a Poisson distribution with parameter $\lambda > 0$, denoted $N \sim \mathcal{P}(\lambda)$, if $N$ is discrete with values in $\mathbb{N}$ such that, for $k \in \mathbb{N}$,
> $$
> \mathbb{P}(N = k) = \frac{\lambda^k}{k!} e^{-\lambda}.
> $$

The parameter $\lambda$ represents the (constant) mean rate of occurence. The Poisson distribution satisfies the following properties.

> **Proposition 1.2.**
> Let $N \sim \mathcal{P}(\lambda)$ with $\lambda > 0$.
> 1. The Poisson distribution is characterized by its moment-generating function (m.g.f.), or Laplace transform:
>    $$
>    L_N : u \in \mathbb{R} \mapsto \mathbb{E}\left[e^{uN}\right] = \exp\left(\lambda(e^u - 1)\right).
>    $$
> 2. Its characteristic function equals:
>    $$
>    \phi_N : u \in \mathbb{R} \mapsto \mathbb{E}\left[e^{iuN}\right] = \exp\left(\lambda(e^{iu} - 1)\right).
>    $$

<!-- page 10 -->
*Proof. In exercise.* Let $N \sim \mathcal{P}(\lambda)$, with $\lambda > 0$. Let $u \le 0$.

1. Compute $\mathbb{E}\left[e^{uN}\right]$.
   **Solution.**
   The m.g.f. of $N$ equals
   $$
   \mathbb{E}\left[e^{uN}\right] = \sum_{k \in \mathbb{N}} e^{uk} \frac{\lambda^k}{k!} e^{-\lambda} = \sum_{k \in \mathbb{N}} \frac{(\lambda e^u)^k}{k!} e^{-\lambda} = \exp(\lambda e^u - \lambda),
   $$
   which ends the proof.

2. Same computation as in 1.

Note that both the moment-generating function and the characteristic function characterize the distribution of a random variable. Hence, if two random variables have the same moment-generating function, or the same characteristic function, then they have the same distribution.

> **Property 1.3.**
> Let $N \sim \mathcal{P}(\lambda)$ with $\lambda > 0$. Then $\mathbb{E}[N] = \text{Var}(N) = \lambda$.

To prove Property 1.3, we may use Lemma 1.4.

> **Lemma 1.4 (Admitted).**
> If $L_X : u \mapsto \mathbb{E}\left[e^{uX}\right]$ denotes the moment-generating function of a real r.v. $X$, then, for all $k$ in $\mathbb{N}^*$, the $k$th-order moment of $X$ is linked to the $k$th-order derivative of $L_X$ as follows:
> $$
> \mathbb{E}\left[X^k\right] = L_X^{(k)}(0).
> $$

*Proof of Property 1.3. In exercise.* Deduce from Proposition 1.2 and Lemma 1.4, $\mathbb{E}[N]$ and $\text{Var}(N)$.
**Solution.**
Recall that
$$
\forall u \in \mathbb{R}, \quad L_N(u) = \exp\left(\lambda(e^u - 1)\right).
$$

- **Expectation:** Since for all $u$ in $\mathbb{R}$,
  $$
  L'_N(u) = \lambda e^u \exp\left(\lambda(e^u - 1)\right) = \lambda \exp(\lambda e^u + u - \lambda),
  $$
  we directly obtain
  $$
  \mathbb{E}[N] = L'_N(0) = \lambda.
  $$

- **Variance:** Since for all $u$ in $\mathbb{R}$,
  $$
  L''_N(u) = \lambda(\lambda e^u + 1)\exp(\lambda e^u + u - \lambda),
  $$
  we directly obtain
  $$
  \mathbb{E}\left[N^2\right] = L''_N(0) = \lambda(\lambda + 1),
  $$
  and thus
  $$
  \text{Var}(N) = \mathbb{E}\left[N^2\right] - \mathbb{E}[N]^2 = \lambda.
  $$

<!-- page 11 -->
## II. Exponential random variable

The exponential distribution is a continuous probability distribution commonly used to model lifetimes, or more generally time elapsed between events.

> **Definition 1.5.**
> Let $X$ be a real random variable. Then $X$ is said to have an exponential distribution with parameter $\lambda > 0$, denoted $X \sim \mathcal{E}(\lambda)$, if $X$ has a density w.r.t. the Lebesgue measure defined on $\mathbb{R}$ by
> $$
> x \mapsto \lambda e^{-\lambda x} \mathbf{1}_{\{x > 0\}}.
> $$

The exponential distribution can be characterized by the following functions.

> **Proposition 1.6.**
> Let $X \sim \mathcal{E}(\lambda)$ with $\lambda > 0$.
> 1. Its moment-generating function (m.g.f.), or Laplace transform, equals for all $u < \lambda$,
>    $$
>    L_X(u) := \mathbb{E}\left[e^{uX}\right] = \frac{\lambda}{\lambda - u}.
>    $$
> 2. Its distribution function equals:
>    $$
>    F : t \in \mathbb{R} \mapsto (1 - e^{-\lambda t})\mathbf{1}_{\{t \ge 0\}}.
>    $$
> 3. Its survival function satisfies:
>    $$
>    \forall t \ge 0, \quad \mathbb{P}(X > t) = e^{-\lambda t}.
>    $$

*Proof. In exercise.*

1. Let $u < \lambda$. Prove that $\mathbb{E}\left[e^{uX}\right] = \lambda \int_0^{+\infty} e^{-(\lambda - u)x}\,dx$. Deduce that $\mathbb{E}\left[e^{uX}\right]$ is well-defined only if $\lambda - u > 0$ and compute the Laplace transform of $X$.
   **Solution.**
   $$
   \mathbb{E}\left[e^{uX}\right] = \int_0^{+\infty} e^{ux} \lambda e^{-\lambda x}\,dx = \lambda \int_0^{+\infty} e^{-(\lambda - u)x}\,dx \stackrel{(*)}{=} \lambda \left[ \frac{-1}{\lambda - u} e^{-(\lambda - u)x} \right]_0^{+\infty} = \frac{\lambda}{\lambda - u}.
   $$
   Note that the integral $(*)$ is well-defined only if $\lambda - u > 0$. Note that in this case, the Laplace transform is positive as expected.

2. Compute $\mathbb{P}(X \le t)$ by distinguishing the cases $t < 0$ and $t \ge 0$.
   **Solution.**
   - Since $X$ takes positive values, if $t < 0$, then $\mathbb{P}(X \le t) = 0$.
   - Let $t \ge 0$. Then
     $$
     \mathbb{P}(X \le t) = \int_0^t \lambda e^{-\lambda x}\,dx = \left[-e^{-\lambda x}\right]_0^t = 1 - e^{-\lambda t}.
     $$

3. Immediate from 2 since $\mathbb{P}(X > t) = 1 - \mathbb{P}(X \le t)$.

<!-- page 12 -->
> **Property 1.7.**
> Let $X \sim \mathcal{E}(\lambda)$ with $\lambda > 0$. Then $\mathbb{E}[X] = 1/\lambda$ and $\text{Var}(X) = 1/\lambda^2$.

*Proof. In exercise, using Lemma 1.4 as in the proof of Property 1.6.*
**Solution.**
Let $X \sim \mathcal{E}(\lambda)$. By Proposition 1.6, its m.g.f. equals for all $u \in ]-\infty; \lambda[$, $L_X(u) = \frac{\lambda}{\lambda - u}$.
Then $L_X$ is differentiable on $]-\infty; \lambda[$, with for all $u \le \lambda$,
$$
L'_X(u) = \frac{\lambda}{(\lambda - u)^2} \quad \text{and} \quad L''_X(u) = \frac{2\lambda}{(\lambda - u)^3}.
$$
Hence, since $0 \in ]-\infty; \lambda[$,
$$
\mathbb{E}[X] = L'_X(0) = \frac{1}{\lambda}, \quad \mathbb{E}\left[X^2\right] = L''_X(0) = \frac{2}{\lambda^2}, \quad \text{and} \quad \text{Var}(X) = \mathbb{E}\left[X^2\right] - \mathbb{E}[X]^2 = \frac{1}{\lambda^2}.
$$

> **Property 1.8.**
> Let $X_1 \sim \mathcal{E}(\lambda_1)$ and $X_2 \sim \mathcal{E}(\lambda_2)$ be two independent exponential random variables. Then
> $$
> \mathbb{P}(X_1 < X_2) = \frac{\lambda_1}{\lambda_1 + \lambda_2}.
> $$

*Proof. In exercise. Hint: $\mathbb{P}(X_1 < X_2) = \mathbb{E}\left[\mathbf{1}_{\{X_1 < X_2\}}\right]$.*
**Solution.**
$$
\begin{aligned}
\mathbb{P}(X_1 < X_2) &= \mathbb{E}\left[\mathbf{1}_{\{X_1 < X_2\}}\right] = \int_{\mathbb{R}^2} \mathbf{1}_{\{x_1 < x_2\}} \lambda_1 e^{-\lambda_1 x_1}\mathbf{1}_{\{x_1 > 0\}} \lambda_2 e^{-\lambda_2 x_2}\mathbf{1}_{\{x_2 > 0\}}\,dx_1 dx_2 \\
&= \int_0^{+\infty} \lambda_1 e^{-\lambda_1 x_1} \left( \int_{x_1}^{+\infty} \lambda_2 e^{-\lambda_2 x_2}\,dx_2 \right) dx_1 \\
&= \int_0^{+\infty} \lambda_1 e^{-\lambda_1 x_1} \left[-e^{-\lambda_2 x_2}\right]_{x_1}^{+\infty}\,dx_1 \\
&= \int_0^{+\infty} \lambda_1 e^{-(\lambda_1 + \lambda_2)x_1}\,dx_1 \\
&= \left[ \frac{-\lambda_1}{\lambda_1 + \lambda_2} e^{-(\lambda_1 + \lambda_2)x_1} \right]_0^{+\infty} \\
&= \frac{\lambda_1}{\lambda_1 + \lambda_2}.
\end{aligned}
$$

**Example:** A system consists of two independent components in series (the system fails if at least one of the components fails). Denote $X_i$ the time (in years) to failure of component $i$ and assume that $X_1 \sim \mathcal{E}(2)$ and $X_2 \sim \mathcal{E}(6)$. In average, component 1 runs 6 months before failure, whereas component 2 runs only 2 months before failure. So in average, component 2 fails before component 1. In fact, according to Property 1.8, the probability that the first component fails before the second one equals $1/4$ (it is coherent).

<!-- page 13 -->
## III. Memoryless distribution

> **Proposition 1.9.**
> Let $(X_i)_{1 \le i \le n}$ be independent exponential r.v. such that for all $1 \le i \le n$, $X_i \sim \mathcal{E}(\lambda_i)$. Then
> $$
> \min_{1 \le i \le n} X_i \sim \mathcal{E}\left(\sum_{i=1}^n \lambda_i\right).
> $$

*Proof. In exercise. Hint: Compute the survival function of $\min_{1 \le i \le n} X_i$.*
**Solution.**
Let $t > 0$. The survival function equals
$$
\begin{aligned}
\mathbb{P}\left(\min_{1 \le i \le n} X_i > t\right) &= \mathbb{P}\left( \bigcap_{1 \le i \le n} \{X_i > t\} \right) \stackrel{(*)}{=} \prod_{1 \le i \le n} \mathbb{P}(X_i > t) \\
&= \prod_{1 \le i \le n} e^{-\lambda_i t} = \exp\left( -\left[ \sum_{i=1}^n \lambda_i \right] t \right),
\end{aligned}
$$
where $(*)$ comes from the independence between the $X_i$'s. We recognize the survival function of an exponential r.v. with parameter $\left(\sum_{i=1}^n \lambda_i\right)$.

**Example:** A system consists of 12 independent components in series with exponential time to failure (in years) $\mathcal{E}(2)$. Individually, the average time to failure of each component is 6 months. However, the time to failure of the system has an exponential distribution $\mathcal{E}(24)$. In particular, its average time to failure is about 15 days.

> **Definition 1.10.**
> A random variable $X$ is said to have a memoryless distribution if for all $t, s > 0$,
> $$
> \mathbb{P}(X > t + s \mid X > t) = \mathbb{P}(X > s). \tag{1.1}
> $$

Note that this definition is well defined if for all $t > 0$, $\mathbb{P}(X > t) > 0$.

**Example:** In reliability theory, if we think of $X$ being the time to failure of a component, then (1.1) states that the probability that the component still runs after $s + t$ hours knowing that it is still working after $t$ hours is the same as the initial probability that it runs for at least $s$ hours. In other words, if the component is still working at time $t$, then the distribution of the remaining amount of time that it works is the same as the original operating distribution; that is, the instrument does not "remember" that is has already been used for a time $t$.

> **Proposition 1.11.**
> The exponential distribution is the only positive continuous memoryless distribution, i.e.,
> 1. If $X \sim \mathcal{E}(\lambda)$, then $X$ has a memoryless distribution.
> 2. If $X$ has a continuous and memoryless distribution, then there exists $\lambda > 0$ such that $X \sim \mathcal{E}(\lambda)$.

<!-- page 14 -->
> **Comment 1.12.** Discrete memoryless distributions are geometric distributions.

The proof of Proposition 1.11 relies on the following Lemma.

> **Lemma 1.13.**
> The solutions on $\mathbb{R}_+$ of the functional equation
> $$
> f(x + y) = f(x)f(y), \quad \forall x, y \in \mathbb{R}_+
> $$
> are either the null function, or $f : x \mapsto e^{ax}$, with $a$ in $\mathbb{R}$.

The proof of Lemma 1.13 and its extension on $\mathbb{R}$ is detailed in Appendix A.III.1.

*Proof of Proposition 1.11. In exercise.*

1. Let $X \sim \mathcal{E}(\lambda)$. Let $t, s \ge 0$ and compute $\mathbb{P}(X > t + s \mid X > t)$.
   **Solution.**
   $$
   \mathbb{P}(X > t + s \mid X > t) = \frac{\mathbb{P}(X > t + s)}{\mathbb{P}(X > t)} = \frac{\exp(-\lambda(t + s))}{\exp(-\lambda t)} = \exp(-\lambda s) = \mathbb{P}(X > s).
   $$

2. Let $X$ be a positive continuous r.v. with a memoryless distribution. Denote $G : t \mapsto \mathbb{P}(X > t)$ its survival function.
   - **a)** Justify that $G$ is a continuous function such that $G(0) = 1$.
     **Solution.**
     Since $X$ is a continuous r.v., then its cumulative distribution function, and thus its survival function, are continuous functions. Moreover, $G(0) = \mathbb{P}(X > 0) = 1$ since $X$ is positive.

   - **b)** Let $s, t \ge 0$. Prove that $G(t + s) = G(s)G(t)$.
     **Solution.**
     By definition of the memoryless property,
     $$
     \begin{aligned}
     G(t + s) &= \mathbb{P}(X > t + s) = \mathbb{P}(\{X > t + s\} \cap \{X > t\}) \\
     &= \mathbb{P}(X > t + s \mid X > t) \mathbb{P}(X > t) \\
     &= \mathbb{P}(X > s) \mathbb{P}(X > t) = G(s)G(t).
     \end{aligned}
     $$

   - **c)** Deduce from Lemma 1.13, that for all $t \ge 0$, $G(t) = e^{-\lambda t}$ where $\lambda$ is a positive constant.
     **Solution.**
     Since $G$ is a continuous function satisfying the functional equality of Lemma 1.13, and since $G$ is not the null function as $G(0) = 1$, then there exists $a$ in $\mathbb{R}$ such that for all nonnegative $t$,
     $$
     G(t) = \exp(at).
     $$
     Moreover, since $\lim_{t \to +\infty} G(t) = 0$, then $a < 0$. In particular, setting $\lambda = -a > 0$ leads to $G(t) = e^{-\lambda t}$.

   - **d)** Conclude.
     **Solution.**
     We recognize the survival function of an exponential distribution $\mathcal{E}(\lambda)$.

<!-- page 15 -->
## IV. Hazard rate in lifetime models

In this section, we only consider non-negative continuous real-valued random variables, that represent lifetimes. For instance, it may be the time before death or relapse of a patient in a medical study, the time before failure of a component in reliability theory, the time before a new claim in actuarial science, or the time before an aircraft loses performance in aeronautics.

The hazard rate represents the instantaneous probability of occurence of the event. It is mathematically defined as follows.

> **Definition 1.14.**
> Let $X$ be a non-negative real random variable with density function $f$ w.r.t. the Lebesgue measure. The hazard rate function $h$ is defined by
> $$
> h(t) = \begin{cases}
> \frac{f(t)}{\mathbb{P}(X > t)} & \text{if } \mathbb{P}(X > t) > 0, \\
> 0 & \text{otherwise.}
> \end{cases}
> $$

**Remark:** Notice that, heuristically,
$$
h(t)\,dt = \frac{f(t)\,dt}{\mathbb{P}(X > t)} \approx \frac{\mathbb{P}(X \in ]t, t + dt])}{\mathbb{P}(X > t)} = \mathbb{P}(X \in ]t, t + dt] \mid X > t).
$$
For instance, in reliability theory, when $X$ denotes the (random) time to failure of a component, the hazard rate corresponds to the probability of a failure occurring just after time $t$ given that the device is still working at time $t$. The survival function $R : t \mapsto \mathbb{P}(X > t)$ is called reliability function, and the hazard rate is called failure rate and is usually denoted $\lambda$. In particular, if the failure rate $h$ is decreasing, it means that the instantaneous probability of failure decreases with the age; the component state improves. On the contrary, if the failure rate $h$ is increasing, it means that the instantaneous probability of failure increases with the age; the component deteriorates. In the reliability theory, the failure rate is usually shaped as a bathtub curve (see below).

```
Hazard rate h(t)
 ^
 | \                                    /
 |  \                                  /
 |   \________________________________/
 |    Running-in    Constant rate   Wearing-out
 +---------------------------------------------> Time t
```
*(Bathtub curve of failure rate over time)*

The survival function and the hazard rate are linked as below.

> **Proposition 1.15.**
> Let $X$ be a continuous non-negative r.v. with positive survival function (i.e. $X$ is not upper bounded). Then, the hazard rate of $X$ is $h$, if and only if, for all $t > 0$,
> $$
> \mathbb{P}(X > t) = \exp\left( -\int_0^t h(x)\,dx \right).
> $$
> In particular, we deduce that the hazard rate characterizes the distribution of a r.v.

<!-- page 16 -->
*Proof. In exercise. Denote $G : t \mapsto \mathbb{P}(X > t)$ the survival function of $X$.*

$\Rightarrow$ Denote $h$ the hazard rate of $X$.

- **a)** Prove that for all $t > 0$, $\mathbb{P}(X \le t) = \int_0^t h(s)\mathbb{P}(X > s)\,ds$, and deduce that
  $$
  G(t) = 1 - \int_0^t h(s)G(s)\,ds. \tag{1.2}
  $$
  **Solution.**
  Since the hazard rate of $X$ is $h$, then its density equals for all $s > 0$, $f(s) = h(s)\mathbb{P}(X > s)$. In particular, its c.d.f equals
  $$
  \mathbb{P}(X \le t) = \int_0^t f(s)\,ds = \int_0^t h(s)\mathbb{P}(X > s)\,ds,
  $$
  that is $1 - G(t) = \int_0^t h(s)G(s)\,ds$.

- **b)** Deduce that $G$ is differentiable and satisfies
  $$
  \begin{cases}
  G(0) = g_0 \\
  G'(t) = -h(t)G(t),
  \end{cases}
  $$
  where $g_0$ is a constant to determine.
  **Solution.**
  Since $G(0) = \mathbb{P}(X > 0) = 1$, $g_0 = 1$. By the Fundamental Theorem of Calculus, since $h$ and $G$ are integrable, $t \mapsto \int_0^t h(s)G(s)\,ds$ is differentiable, hence $G$ is differentiable on $\mathbb{R}_+^*$ and $G'(t) = -h(t)G(t)$.

- **c)** Solve this differential equation and conclude.
  **Solution.**
  The unique solution of $G'(t) = -h(t)G(t)$ with $G(0) = 1$ is
  $$
  G(t) = G(0)\exp\left(-\int_0^t h(s)\,ds\right) = \exp\left(-\int_0^t h(s)\,ds\right).
  $$

$\Leftarrow$ Assume that for all $t > 0$, $\mathbb{P}(X > t) = \exp\left(-\int_0^t \psi(x)\,dx\right)$, where $\psi$ is a positive continuous function.

- **i)** Compute the density function of $X$.
  **Solution.**
  The c.d.f. is $F(t) = 1 - G(t) = 1 - \exp\left(-\int_0^t \psi(x)\,dx\right)$. Differentiating w.r.t. $t$ gives the density:
  $$
  f(t) = F'(t) = \psi(t)\exp\left(-\int_0^t \psi(x)\,dx\right) = \psi(t)G(t).
  $$

<!-- page 17 -->
- **ii)** Deduce the hazard rate of $X$.
  **Solution.**
  Hence, by definition, since for all $t$, $G(t) = \mathbb{P}(X > t)$ is positive, the hazard rate of $X$ equals
  $$
  h(t) = \frac{f(t)}{G(t)} = \psi(t),
  $$
  which ends the proof.

> **Property 1.16.**
> Let $\lambda > 0$. The random variable $X \sim \mathcal{E}(\lambda)$ if and only if its hazard rate function is constant: $h(t) = \lambda$ for all $t > 0$.

*Proof. In exercise.*

$\Rightarrow$ Assume $X \sim \mathcal{E}(\lambda)$ and compute the hazard rate function of $X$.
**Solution.**
Recall the density function $f$ of $X$ equals for all $t > 0$, $f(t) = \lambda e^{-\lambda t}$, and the survival function is $\mathbb{P}(X > t) = e^{-\lambda t}$. Hence the hazard rate function of $X$ equals
$$
h(t) = \frac{f(t)}{\mathbb{P}(X > t)} = \lambda.
$$

$\Leftarrow$ Assume the hazard rate function of $X$ is constant equal to $\lambda$. Prove that $X \sim \mathcal{E}(\lambda)$ using Proposition 1.15 without any computation.
**Solution.**
By Proposition 1.15, the hazard rate characterizes the distribution. Hence a r.v. with constant hazard rate function equal to $\lambda$ has the same hazard rate function, and thus the same distribution than an exponential r.v. (by $\Rightarrow$), that is $\mathcal{E}(\lambda)$.

<!-- page 18 -->
## V. On the gamma and Erlang's distributions

The Erlang distribution generalizes the exponential distribution. It is commonly used to model the waiting time between more than one events. It was developed by A. K. Erlang to examine the number of telephone calls which might be made at the same time to the operators of the switching stations. More precisely, the Erlang distribution is the distribution of the sum of $n$ i.i.d. r.v. with exponential distribution. The gamma distribution then generalizes the Erlang distribution by allowing $n$ to be any positive real number.

> **Definition 1.17.**
> - **Gamma distribution.** Let $\lambda > 0$ and $\alpha > 0$. A non-negative real-valued random variable $X$ with density (w.r.t. the Lebesgue measure)
>   $$
>   x \mapsto \frac{(\lambda x)^{\alpha - 1}}{\Gamma(\alpha)} \lambda e^{-\lambda x} \mathbf{1}_{\{x > 0\}},
>   $$
>   where $\Gamma(\alpha) = \int_0^{+\infty} e^{-t} t^{\alpha - 1}\,dt$ is the gamma function, is said to have a gamma distribution with parameter $(\alpha, \lambda)$. We denote $X \sim \Gamma(\alpha, \lambda)$.
> - **Erlang distribution.** If $\alpha = n \in \mathbb{N}^*$, then the gamma distribution $\Gamma(n, \lambda)$ is called Erlang's distribution. In particular, since $\Gamma(n) = (n-1)!$, the density of a r.v. $X \sim \Gamma(n, \lambda)$ equals
>   $$
>   x \mapsto \frac{\lambda^n x^{n-1}}{(n - 1)!} e^{-\lambda x} \mathbf{1}_{\{x > 0\}}.
>   $$

The gamma distribution can be characterized by its moment-generating function.

> **Proposition 1.18.**
> Let $\alpha, \lambda > 0$ and consider $X \sim \Gamma(\alpha, \lambda)$.
> 1. The moment-generating function, or Laplace transform of $X$ equals for all $u < \lambda$,
>    $$
>    L_X(u) := \mathbb{E}\left[e^{uX}\right] = \left( \frac{\lambda}{\lambda - u} \right)^\alpha.
>    $$
> 2. We deduce that $\mathbb{E}[X] = \alpha/\lambda$ and $\text{Var}(X) = \alpha/\lambda^2$.

*Proof. In exercise (Worksheet 1, Exercise 2).*

> **Proposition 1.19.**
> Links with the exponential and the Chi-square distributions.
> - **i)** If $X \sim \Gamma(\alpha, \lambda)$, then, for all $\mu > 0$, $\mu X \sim \Gamma(\alpha, \lambda/\mu)$.
> - **ii)** $\Gamma(1, \lambda) \sim \mathcal{E}(\lambda)$.
> - **iii)** If $X_1, \dots, X_n$ are i.i.d. with distribution $\mathcal{E}(\lambda)$, then $\sum_{i=1}^n X_i \sim \Gamma(n, \lambda)$.
> - **iv)** $\Gamma(n, 1/2) \stackrel{(d)}{=} \chi^2(2n)$.

<!-- page 19 -->
*Proof. In exercise.*

- **i)** Use the m.g.f. to prove this point.
  **Solution.**
  For all $u < \lambda/\mu$, the m.g.f. of $\mu X$ equals
  $$
  \mathbb{E}\left[e^{u(\mu X)}\right] = \mathbb{E}\left[e^{(\mu u)X}\right] = \left( \frac{\lambda}{\lambda - \mu u} \right)^\alpha = \left( \frac{\lambda/\mu}{\frac{\lambda}{\mu} - u} \right)^\alpha.
  $$
  We recognize the m.g.f. of a $\Gamma(\alpha, \lambda/\mu)$ distribution.

- **ii)** By taking $\alpha = 1$ in the definition, since $\Gamma(1) = 0! = 1$, we recognize the density of an exponential r.v. with parameter $\lambda$.

- **iii)** Worksheet 1, Exercise 2.

- **iv)** Immediate since the density of a $\chi^2(d)$ is
  $$
  x \mapsto \frac{1}{2^{d/2}\Gamma(d/2)} x^{d/2 - 1} e^{-x/2} \mathbf{1}_{\{x > 0\}}.
  $$

<!-- page 20 -->

<!-- page 21 -->
# Chapter 2. Homogeneous Poisson processes

Homogeneous Poisson processes may provide very simple models in Reliability theory. For instance, if at each failure a system is repaired as new and placed in service again, then the failures can be modeled by an homogeneous Poisson process, provided the repair times can be neglected.

## I. Definitions

Before defining homogeneous Poisson processes, let us now introduce the notion of point processes on $\mathbb{R}_+$. We observe the occurrences times of a phenomenon, such as failures in reliability theory, or claims in actuarial sciences.

### I.1. Point processes and counting processes

There are several ways of defining a point process depending on which quantity is considered.

#### A) Counting process

The first perspective consists in counting the number of events occurring in any time interval (or time window).

> **Definition 2.1 (Counting process).**
> A stochastic process $N = (N_t)_{t \in \mathbb{R}_+}$ is said to be a counting process if, for all $t \ge 0$, $N_t$ represents the total number of "events" that occur by time $t$.

> **Example 2.2.**
> Here are some examples of counting processes.
> - **a)** In reliability theory, one may denote $N_t$ the number of failures of a reparable component on or before time $t$. An event thus corresponds to a failure.
> - **b)** In actuarial science, $N_t$ may represent the number of accidents an insurance company has to cover by time $t$. Here an event is a claim.
> - **c)** If we say an event occurs whenever a child is born, then $N_t$ equals the total number of people who were born by time $t$.
> - **d)** Let $N_t$ count the number of goals scored by the French team during the Football World Cup final after $t$ minutes of game. An event of this process occurs when a player scores a goal.

<!-- page 22 -->
> **Property 2.3.**
> A counting process satisfies the following properties.
> - **(i)** For all $t \ge 0$, $N_t \ge 0$.
> - **(ii)** $t \mapsto N_t$ is non-decreasing, that is for all $s \le t$, $N_s \le N_t$.
> - **(iii)** $t \mapsto N_t$ is piecewise constant, with values in $\mathbb{N}$.
> - **(iv)** $(N_t)_{t \in \mathbb{R}_+}$ is càdlàg, that is right-continuous and left-limited (*continue à droite, limitée à gauche* in French).
> - **(v)** For all $0 \le s < t$, $N_t - N_s$ counts the number of events that occur in the interval $(s, t]$.

In this lecture, we restrict attention to counting processes $N$ whose realizations are locally finite (i.e. the number of events in any bounded domain in $\mathbb{R}_+$ is finite). In particular, this is true as soon as for all $t \ge 0$, $N_t$ is almost surely finite.

#### B) Point process

The second perspective consists in considering the times each "event" occurs.

> **Definition 2.4 (Point process).**
> A point process on $\mathbb{R}_+$ is a random countable subset of $\mathbb{R}_+$, each point representing the time occurence of an "event". It may be formally defined as the intersection between $\mathbb{R}_+$ and the non-decreasing sequence of times the event occurs, also called arrival times, and denoted
> $$
> 0 \le T_1 \le T_2 \le \dots \le T_n \le \dots
> $$
> where the $T_n$ belong to $\mathbb{R}_+ \cup \{+\infty\}$.

The knowledge of the counting process $N$ is equivalent to the knowledge of the arrival times $0 \le T_1 \le T_2 \le \dots \le T_n \le \dots$

```
count Nt
  ^
5 |                                              •-------
4 |                                      •-------
3 |                             •--------
2 |                    •--------
1 |           •--------
0 +-----------•--------•--------•--------•-------•-------> time t
              T1       T2       T3       T4      T5
```
*(Staircase trajectory of counting process $N_t$ jumping by 1 at each arrival time $T_n$)*

Indeed, on the one hand, if we know the counting process, it is possible to recover the arrival times as
$$
T_n = \inf\{t \ge 0 \; ; \; N_t \ge n\}, \tag{2.1}
$$
with the convention that $\inf\{\emptyset\} = +\infty$.

<!-- page 23 -->
On the other hand, if we know the point process, it is possible to define the corresponding counting process for all $t \ge 0$ by
$$
N_t = \sum_{n \ge 1} \mathbf{1}_{\{T_n \le t\}}. \tag{2.2}
$$
Notice that this sum converges as soon as the point process does not have any accumulation point, that is the process is locally finite.

> **Proposition 2.5 (Link between counting and point processes).**
> Let $N$ be a counting process and denote $(T_n)_{n \in \mathbb{N}^*}$ the corresponding arrival times. Then, for all $n$ in $\mathbb{N}^*$ and all $t \ge 0$,
> - **(i)** $\{T_n \le t\} = \{N_t \ge n\}$.
> - **(ii)** $\{T_n > t\} = \{N_t < n\}$.
> - **(iii)** $\{N_t = n\} = \{T_n \le t < T_{n+1}\}$.
>
> More generally, for all $n$ in $\mathbb{N}^*$ and all $(t_1, t_2, \dots, t_n)$ in $(\mathbb{R}_+)^n$,
> $$
> \{T_1 \le t_1, T_2 \le t_2, \dots, T_n \le t_n\} = \{N_{t_1} \ge 1, N_{t_2} \ge 2, \dots, N_{t_n} \ge n\}.
> $$
> It is thus equivalent to define the distribution of the point process by the one of $(T_n)_{n \ge 1}$ or the one of $(N_t)_{t \in \mathbb{R}_+}$.

From now on, unless specified, $N$ refers to both the counting process $(N_t)_{t \in \mathbb{R}_+}$ and its associated point process $(T_n)_{n \in \mathbb{N}^*}$, and we do not distinguish between "counting" and "point" processes.

#### C) Some basic concepts

In this lecture, we only consider simple point processes, for which several events cannot occur at the same time.

- On the one hand, from a "point process" perspective, this means that the random variables $(T_n)$ are a.s. pairwise distinct on $\mathbb{R}_+$. More precisely, for all $n$ in $\mathbb{N}^*$ and almost all $\omega$ in $\Omega$,
  $$
  T_n(\omega) < +\infty \implies T_n(\omega) < T_{n+1}(\omega).
  $$
  Essentially, it means that for almost all $\omega \in \Omega$, the sequence $(T_n(\omega))_{n \in \mathbb{N}^*}$ is increasing as long as the $T_n(\omega)$ are finite.
- On the other hand, from a "counting process" perspective, it means that all the jumps of the counting process $(N_t)_{t \in \mathbb{R}_+}$ are equal to 1.

> **Definition 2.6 (Regular point process).**
> A point process $N$ is said to be regular if for all $t \ge 0$,
> $$
> \lim_{h \to 0} \frac{\mathbb{P}(N_{t+h} - N_t \ge 2)}{h} = 0.
> $$

For instance in reliability theory, it means that the system will not experience two or more failures simultaneously. In particular, a regular point process is simple.

<!-- page 24 -->
> **Definition 2.7 (Mean function and rate of a counting process).**
> Let $N$ be a counting process.
> - The **mean function** $m$ of $N$ is defined by
>   $$
>   m : t \in \mathbb{R}_+ \mapsto m(t) = \mathbb{E}[N_t].
>   $$
> - The **rate** $w$ of $N$ is defined by
>   $$
>   w : t \in \mathbb{R}_+ \mapsto w(t) = m'(t) = \lim_{h \to 0} \frac{\mathbb{E}[N_{t+h} - N_t]}{h}.
>   $$

> **Comment 2.8.**
> - **i.** Note that the mean function and the rate of a counting process are deterministic.
> - **ii.** Moreover, even if the counting process is piecewise constant (it has "jumps"), its mean function is not necessarily discontinuous.
> - **iii.** If $N$ is a regular point process, the probability of two or more events in $]t, t + h]$ is negligible when $h$ is small. Thus, we may assume that, for small $h$, $N_{t+h} - N_t$ takes values in $\{0, 1\}$, and that
>   $$
>   \mathbb{E}[N_{t+h} - N_t] \approx \mathbb{P}(N_{t+h} - N_t = 1).
>   $$
>   Hence, the rate of the process is approximately equal to
>   $$
>   w(t) \approx \frac{\mathbb{P}(N_{t+h} - N_t = 1)}{h}.
>   $$
>   Hence, we can think of $w(t)\,dt$ as the probability of an occurence "just after time $t$".
> - **iv.** In reliability theory, the rate is called **Rate of Occurence of Failure (RoCoF)**. It represents the mean number of failures per unit of time.
> - **v.** The rate can be generalized to the conditional rate, in the case where we work conditionally on all the history of the process. Unlike the rate, the conditional rate is usually stochastic. In reliability theory, it is called **Failure intensity**.

#### D) Increments of a point process

> **Definition 2.9 (Stationary increments).**
> A point process $N$ is said to have stationary increments if the distribution of the number of points in any interval of time depends only on the length of the interval. In other words, for any $s \ge 0$, the number of points in the interval $(t, t + s]$, that is $N_{t+s} - N_t$ has the same distribution for all $t$.

> **Example (2.2 continued).**
> Recall the examples enumerated above.
> - **a)** In reliability theory, the stationarity seems reasonable if we consider that there are no running-in period or aging effect, and if after a repair, the component is "as new" (this case is often referred to as perfect repair).
> - **b)** In actuarial science, it is a reasonable assumption if, for instance, there are no periods in the year where the accidents seem more likely to happen.

<!-- page 25 -->
> - **c)** If we believe that the earth's population is basically constant (a belief not held at present by most scientists), then the assumption of stationary increments might be reasonable.
> - **d)** Not reasonable if we take into account the tiredness of the players.

> **Definition 2.10 (Independent increments).**
> A point process $N$ is said to have independent increments if the number of points belonging to disjoint intervals are independent. In other words, for all $n$ in $\mathbb{N}^*$, for all $0 < t_1 < t_2 < \dots < t_n$, the random variables $N_{t_1} - N_0, N_{t_2} - N_{t_1}, \dots, N_{t_n} - N_{t_{n-1}}$ are independent.

> **Example (2.2 continued).**
> Recall the examples enumerated above.
> - **a)** In reliability theory, it seems reasonable to assume that the increments are independent if we can assume that a failure does not have repercussions (e.g. perfect repair case).
> - **b)** It depends on the kind of insurance: in life insurance, the increments can be assumed as independent (the number of deaths at time $t$ does not impact the other customers). Yet, in non-life insurance, the independence assumption for the increments is less realistic. For instance, if an earthquake happens, then there will probably be aftershocks and thus other damages.
> - **c)** No. If, for instance, $N_s$ is very large, then it is probable that there are many people alive at time $s$; this would lead us to believe that the number of new births between time $s$ and time $s + t$ would also tend to be large. Hence it does not seem reasonable to assume that $N_s$ is independent of $N_{s+t} - N_s$.
> - **d)** Could be justified if we assume the players' chances of scoring a goal do not depend on how the game has been going (unlikely).

### I.2. Homogeneous Poisson process on $\mathbb{R}_+$

The most elementary point process is the homogeneous Poisson point process which is defined as follows.

> **Definition 2.11.**
> A point process $N$ is said to be a **homogeneous Poisson process** with rate $\lambda > 0$ if:
> 1. $N_0 = 0$.
> 2. The process has independent increments.
> 3. The number of points in any interval of length $t > 0$ has a Poisson distribution with parameter $\lambda t$, that is, for all $s \ge 0$,
>    $$
>    N_{s+t} - N_s \sim \mathcal{P}(\lambda t).
>    $$
> The parameter $\lambda$ is also called the **intensity** of the Poisson point process.

<!-- page 26 -->
> **Property 2.12.**
> Let $N$ be a homogeneous Poisson process with intensity $\lambda > 0$.
> 1. $N$ has stationary increments.
> 2. The mean function of $N$ equals
>    $$
>    \forall t \ge 0, \quad m(t) = \mathbb{E}[N_t] = \lambda t.
>    $$
>    In particular, the rate (see Definition 2.7) of the Poisson process $N$ is constant (equal, as expected, to $\lambda$):
>    $$
>    \forall t \ge 0, \quad w(t) = \lambda.
>    $$
> 3. One has as $t \to 0$,
>    - **(a)** $\mathbb{P}(N_t = 0) = e^{-\lambda t} = 1 - \lambda t + o(t)$.
>    - **(b)** $\mathbb{P}(N_t = 1) = \lambda t e^{-\lambda t} = \lambda t + o(t)$.
>    - **(c)** $\mathbb{P}(N_t \ge 2) = \sum_{k=2}^{+\infty} \frac{(\lambda t)^k}{k!} e^{-\lambda t} = o(t)$.

*Proof.* Immediate by Definition 2.11, and since $N_t \sim \mathcal{P}(\lambda t)$.

> **Comment 2.13.** In particular, $\forall t \ge 0$,
> $$
> \frac{\mathbb{P}(N_{t+h} - N_t \ge 2)}{h} = \frac{\mathbb{P}(N_h \ge 2)}{h} = o(1) \xrightarrow[h \to 0]{} 0.
> $$
> Hence a homogeneous Poisson process is regular (and thus simple).

> **Definition 2.14.**
> The point process $N$ is said to be a **homogeneous Poisson process** with rate (or with intensity) $\lambda > 0$ if:
> - **(A)** $N_0 = 0$.
> - **(B)** The process has stationary and independent increments.
> - **(C)** As $t \to 0$, $t > 0$,
>   - **(a)** $\mathbb{P}(N_t = 1) = \lambda t + o(t)$.
>   - **(b)** $\mathbb{P}(N_t \ge 2) = o(t)$.

> **Proposition 2.15.**
> Definition 2.11 and Definition 2.14 are equivalent.

*Proof. In exercise.*

$\Rightarrow$ Done in Property 2.12.

$\Leftarrow$ Assume Definition 2.14. Then, the points (1) and (2) in Definition 2.11 are immediate from (A) and (B) in Definition 2.14.
Remains to prove (3). To do so, let us compute the moment-generating function of $N_t$ for all $t > 0$.

<!-- page 27 -->
- **a)** Let $h > 0$. Using Assumption (C), prove that for all $u \le 0$,
  $$
  \mathbb{E}\left[e^{uN_h}\right] = 1 + \lambda h(e^u - 1) + o(h).
  $$
  **Solution.**
  From (C), we obtain
  $$
  \mathbb{P}(N_h = 0) = 1 - [\mathbb{P}(N_h = 1) + \mathbb{P}(N_h \ge 2)] = 1 - \lambda h + o(h).
  $$
  Hence, since $u \le 0$,
  $$
  \begin{aligned}
  \mathbb{E}\left[e^{uN_h}\right] &= \underbrace{\mathbb{P}(N_h = 0)}_{1 - \lambda h + o(h)} + e^u \underbrace{\mathbb{P}(N_h = 1)}_{\lambda h + o(h)} + \sum_{k \ge 2} e^{uk} \mathbb{P}(N_h = k) \\
  &\le 1 - \lambda h + o(h) + e^u(\lambda h + o(h)) + \underbrace{\sum_{k \ge 2} \mathbb{P}(N_h = k)}_{=\mathbb{P}(N_h \ge 2) = o(h)} \\
  &= 1 + \lambda h(e^u - 1) + o(h).
  \end{aligned}
  $$

- **b)** Fix $u \le 0$ and denote $g : t \mapsto \mathbb{E}\left[e^{uN_t}\right]$. We now aim at proving that for any small $h$,
  $$
  g(t + h) = g(t) [1 + \lambda h(e^u - 1) + o(h)]. \tag{2.3}
  $$
  - **i)** Let $h > 0$. Using both properties in Assumption (B), prove that
    $$
    g(t + h) = \mathbb{E}\left[e^{uN_h}\right] g(t),
    $$
    and deduce (2.3). *Hint: $N_{t+h} = (N_{t+h} - N_t) + N_t$.*
    **Solution.**
    By stationarity of the increments, $(N_{t+h} - N_t)$ has the same distribution as $N_h$ $(\star)$, and by independence of the increments, since $(0; t] \cap (t; t + h] = \emptyset$,
    $$
    (N_{t+h} - N_t) \perp\!\!\!\perp N_t - N_0 = N_t \quad (\dagger).
    $$
    Therefore,
    $$
    \begin{aligned}
    g(t + h) = \mathbb{E}\left[e^{uN_{t+h}}\right] &= \mathbb{E}\left[e^{u(N_{t+h} - N_t)} e^{uN_t}\right] \\
    &\stackrel{(\dagger)}{=} \mathbb{E}\left[e^{u(N_{t+h} - N_t)}\right] \mathbb{E}\left[e^{uN_t}\right] \stackrel{(\star)}{=} \mathbb{E}\left[e^{uN_h}\right] g(t).
    \end{aligned}
    $$
    Equation (2.3) is then immediate from a).

  - **ii)** Let $h < 0$. In the same way, prove that
    $$
    g(t) = \mathbb{E}\left[e^{uN_{-h}}\right] g(t + h),
    $$
    and deduce (2.3).

<!-- page 28 -->
    **Solution.**
    As in the previous case, since $N_t = (N_t - N_{t+h}) + N_{t+h}$, by stationarity and independence (since $(0; t + h] \cap (t + h; t] = \emptyset$) of the increments, we obtain
    $$
    \begin{aligned}
    g(t) = \mathbb{E}\left[e^{uN_t}\right] &= \mathbb{E}\left[e^{u(N_t - N_{t+h})} e^{uN_{t+h}}\right] \\
    &= \mathbb{E}\left[e^{u(N_t - N_{t+h})}\right] \mathbb{E}\left[e^{uN_{t+h}}\right] \quad \text{(independence)} \\
    &= \mathbb{E}\left[e^{uN_{-h}}\right] g(t + h) \quad \text{(stationarity)}.
    \end{aligned}
    $$
    We deduce from a) that $g(t) = g(t + h) [1 - \lambda h(e^u - 1) + o(h)]$, and thus
    $$
    g(t + h) = \frac{g(t)}{1 - \lambda h(e^u - 1) + o(h)} = g(t) [1 + \lambda h(e^u - 1) + o(h)].
    $$

- **c)**
  - **i.** Deduce from (2.3) that $g$ is differentiable at point $t$, and satisfies the homogeneous linear differential equation with initial condition $g(0)$ (to determine)
    $$
    \begin{cases}
    g(0) = 1 \\
    g'(t) = g(t)\lambda(e^u - 1).
    \end{cases} \tag{2.4}
    $$
    **Solution.**
    From (2.3), we directly deduce that
    $$
    \frac{g(t + h) - g(t)}{h} = g(t) [\lambda(e^u - 1) + o(1)]
    $$
    has a finite limit when $h \to 0$, which equals
    $$
    g'(t) = \lim_{h \to 0} \frac{g(t + h) - g(t)}{h} = g(t) [\lambda(e^u - 1)].
    $$
    Moreover, $g(0) = \mathbb{E}[e^0] = 1$.

  - **ii.** Solve the differential equation and deduce that
    $$
    g(t) = e^{\lambda t(e^u - 1)}.
    $$
    **Solution.**
    Hence $g$ is the unique solution of the homogeneous linear differential equation with initial condition (2.4), that is
    $$
    g(t) = g(0) e^{\lambda(e^u - 1)t} = e^{\lambda t(e^u - 1)}.
    $$

- **d)** Deduce that $N_t \sim \mathcal{P}(\lambda t)$, and then (3).
  **Solution.**
  We proved that for any $t > 0$ and any $u \le 0$,
  $$
  \mathbb{E}\left[e^{uN_t}\right] = e^{\lambda t(e^u - 1)}.
  $$
  We recognize the m.g.f. of a $\mathcal{P}(\lambda t)$ distributed r.v. (c.f. Proposition 1.2). Finally, by the stationarity assumption, for all $s \ge 0$, $N_{s+t} - N_s$ and $N_t$ have the same distribution, that is $\mathcal{P}(\lambda t)$.

<!-- page 29 -->
## II. Arrival and interarrival times

In all this section, $N$ denotes a point process, and $(T_n)_{n \in \mathbb{N}^*}$ is the sequence of the corresponding arrival times (see Equations (2.1) and (2.2) for more details).

### II.1. First arrival time after a given instant

> **Property 2.16.**
> Let $N$ be a homogeneous Poisson process with rate $\lambda > 0$ and denote $(T_n)_{n \in \mathbb{N}^*}$ the corresponding arrival times.
> - Then, the first arrival time satisfies $T_1 \sim \mathcal{E}(\lambda)$.
> - Moreover, for all fixed $t > 0$, one has
>   $$
>   T_{N_t} \le t < T_{N_t + 1}.
>   $$
>   In particular, the first arrival time after a given instant $t$ satisfies
>   $$
>   T_{N_t + 1} - t \sim \mathcal{E}(\lambda).
>   $$

```
---------|-------------------|-------------------|-------------------> t
       TNt                   t                 TNt+1
```

*Proof. In exercise (c.f. Worksheet 2, Exercise 2).*

> **Example (The inspection paradox or "le paradoxe de l'autobus").**
> You are waiting for a bus to arrive. It is written on your timetable that there is a bus every 10 minutes.
> - If the arrival times of the bus were deterministic, and a bus arrived every 10 minutes, then the time you wait at the bus stop is uniformly distributed on $[0, 10]$ and, in expectation you wait for $10/2 = 5$ minutes.
> - However, due to traffic, the arrival times of the buses are random. You assume they can be modeled by a Poisson process with rate $1/10$ (since in average, there is 1 bus in 10 min). Then, your waiting time equals $T_{N_t + 1} - t$ which has an exponential distribution with parameter $1/10$. Hence, in average, you wait for 10 minutes (instead of the expected 5 mins).

### II.2. Exponential interarrival times

> **Definition 2.17.**
> The **interarrival times** (also called waiting times) of a point process are defined by
> $$
> \begin{cases}
> W_1 = T_1, \\
> W_i = T_i - T_{i-1}, \quad \forall i \ge 2.
> \end{cases}
> $$
> The $i$th waiting time corresponds to the delay between the $(i - 1)$th and the $i$th occurrences. In particular, for all $n$ in $\mathbb{N}^*$,
> $$
> T_n = \sum_{i=1}^n W_i.
> $$

```
                   Wi
              |----------|
--------------|----------|-------------------> t
            Ti-1         Ti
```

<!-- page 30 -->
> **Proposition 2.18.**
> Let $N$ be a homogeneous Poisson process with rate $\lambda$.
> 1. Then the corresponding interarrival times $(W_i)_{i \ge 1}$ are independent and identically distributed exponential random variables with parameter $\lambda$.
> 2. It follows that for all $n$ in $\mathbb{N}^*$, the $n$th arrival time $T_n$ has a gamma (or Erlang) distribution with parameters $n$ and $\lambda$.

Notice that, heuristically, the assumptions of stationarity and independent increments mean that, at any point in time, the process probabilistically restarts itself, independently of the past. That is, the process from any time $t$ is independent of all that has previously occurred, and also has the same distribution as the original process. In other words, the process has no memory, and hence, exponential interarrival times are to be expected (c.f. Section III of Chapter 1).

This result relies on the following Lemma.

> **Lemma 2.19.**
> Let $N$ be a homogeneous Poisson process with rate $\lambda$, and denote $T_1 < T_2 < \dots < T_n < \dots$ the corresponding arrival times. Then, for all $n$ in $\mathbb{N}^*$, $(T_1, \dots, T_n)$ has the following density w.r.t. the Lebesgue measure on $\mathbb{R}^n$:
> $$
> (t_1, \dots, t_n) \mapsto \lambda^n e^{-\lambda t_n} \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}.
> $$

*Idea of proof of Lemma 2.19.*

- **(Reminder) For a real-valued continuous random variable:** Let $X$ be a random variable with density $f$ and c.d.F. $F$. Then, its density
  $$
  f(x) = F'(x) = \lim_{h \to 0} \frac{F(x + h) - F(x)}{h} = \lim_{h \to 0} \frac{\mathbb{P}(X \in (x; x + h])}{h}.
  $$
- **In the multidimensional case:** Similarly, consider an increasing sequence $0 < t_1 < \dots < t_n$, and introduce for $h = (h_1, \dots, h_n)$ small enough such that all $(t_i, t_i + h_i]$ are disjoint:
  $$
  A_h = \{T_1 \in (t_1; t_1 + h_1], \dots, T_n \in (t_n; t_n + h_n]\}.
  $$

```
   T1            T2                       Tn
---|---|---------|---|--------------------|---|-------------> t
   t1  t1+h1     t2  t2+h2                tn  tn+hn
```

Then,
$$
f(t_1, \dots, t_n) = \lim_{\forall i, \, h_i \to 0} \frac{\mathbb{P}(A_h)}{h_1 \dots h_n}.
$$
Yet, the event $A_h$ is satisfied if and only if there are:
- 0 points in $(0, t_1]$ and 1 point in $(t_1, t_1 + h_1]$,
- 0 points in $(t_1 + h_1, t_2]$ and 1 point in $(t_2, t_2 + h_2]$,
- ...
- 0 points in $(t_{n-1} + h_{n-1}, t_n]$ and 1 point in $(t_n, t_n + h_n]$.

<!-- page 31 -->
Thus, since for $h_1, \dots, h_n$ small enough, the intervals appearing in $A_h$ are disjoint, by independence of the increments:
$$
\begin{aligned}
\mathbb{P}(A_h) &= \mathbb{P}(N_{t_1} = 0) \times \mathbb{P}(N_{t_1 + h_1} - N_{t_1} = 1) \times \mathbb{P}(N_{t_2} - N_{t_1 + h_1} = 0) \times \mathbb{P}(N_{t_2 + h_2} - N_{t_2} = 1) \\
&\quad \times \dots \times \mathbb{P}(N_{t_n} - N_{t_{n-1} + h_{n-1}} = 0) \times \mathbb{P}(N_{t_n + h_n} - N_{t_n} = 1) \\
&= e^{-\lambda t_1} \times [\lambda h_1 e^{-\lambda h_1}] \times e^{-\lambda(t_2 - [t_1 + h_1])} \times [\lambda h_2 e^{-\lambda h_2}] \\
&\quad \times \dots \times e^{-\lambda(t_n - [t_{n-1} + h_{n-1}])} \times [\lambda h_n e^{-\lambda h_n}] \\
&= \lambda^n h_1 \dots h_n e^{-\lambda(t_n + h_n)}.
\end{aligned}
$$
Hence,
$$
f(t_1, \dots, t_n) = \lim_{\forall i, \, h_i \to 0} \lambda^n e^{-\lambda(t_n + h_n)} = \lambda^n e^{-\lambda t_n}.
$$
Note that, if we do not have $0 < t_1 < \dots < t_n$, then $\mathbb{P}(A_h) = 0$ since $0 < T_1 < \dots < T_n$ a.s.

We can now prove Proposition 2.18.

*Proof. In exercise.*

1. Let $N$ be a Poisson process with rate $\lambda$, and consider $(W_i)_{i \ge 1}$ the interarrival times. Let us prove that they are i.i.d. with distribution $\mathcal{E}(\lambda)$.
   - **(a)** Fix $n \in \mathbb{N}^*$. Consider $g : \mathbb{R}^n \to \mathbb{R}^n$ a measurable bounded function. Prove that
     $$
     \mathbb{E}[g(W_1, \dots, W_n)] = \int_{\mathbb{R}^n} g(t_1, t_2 - t_1, \dots, t_n - t_{n-1}) \lambda^n e^{-\lambda t_n} \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}\,dt_1 \dots dt_n.
     $$
     **Solution.**
     Consider the arrival times $T_i = \sum_{j=1}^i W_j$. Knowing the density of $(T_1, \dots, T_n)$ given in Lemma 2.19 allows us to write
     $$
     \begin{aligned}
     \mathbb{E}[g(W_1, \dots, W_n)] &= \mathbb{E}[g(T_1, T_2 - T_1, \dots, T_n - T_{n-1})] \\
     &= \int_{\mathbb{R}^n} g(t_1, t_2 - t_1, \dots, t_n - t_{n-1}) \lambda^n e^{-\lambda t_n} \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}\,dt_1 \dots dt_n.
     \end{aligned}
     $$

   - **(b)** By a well-chosen change of variables (see Appendix A.I), show that
     $$
     \mathbb{E}[g(W_1, \dots, W_n)] = \int_{\mathbb{R}^n} g(w_1, w_2, \dots, w_n) \prod_{i=1}^n \left[ \lambda e^{-\lambda w_i} \mathbf{1}_{\{w_i > 0\}} \right] dw_1 \dots dw_n.
     $$
     **Solution.**
     Consider the following change of variable:
     $$
     \begin{cases}
     w_1 = t_1 \\
     w_2 = t_2 - t_1 \\
     \vdots \\
     w_n = t_n - t_{n-1}
     \end{cases}
     \iff
     \begin{cases}
     t_1 = w_1 = \phi_1(w) \\
     t_2 = w_1 + w_2 = \phi_2(w) \\
     \vdots \\
     t_n = \sum_{i=1}^n w_i = \phi_n(w)
     \end{cases}
     $$

<!-- page 32 -->
     The function $\phi : w = (w_1, \dots, w_n) \mapsto (w_1, w_1 + w_2, \dots, \sum_{i=1}^n w_i)$ is a one-to-one $C^1$-diffeomorphism with Jacobian matrix
     $$
     J_\phi(w) = \begin{pmatrix}
     \frac{\partial \phi_1}{\partial w_1} & \dots & \frac{\partial \phi_1}{\partial w_n} \\
     \vdots & \ddots & \vdots \\
     \frac{\partial \phi_n}{\partial w_1} & \dots & \frac{\partial \phi_n}{\partial w_n}
     \end{pmatrix} = \begin{pmatrix}
     1 & 0 & \dots & 0 \\
     1 & 1 & \dots & 0 \\
     \vdots & \vdots & \ddots & \vdots \\
     1 & 1 & \dots & 1
     \end{pmatrix}.
     $$
     Hence, for all $w$ in $\mathbb{R}^n$, $\det(J_\phi(w)) = 1 \ne 0$. Moreover,
     $$
     \mathbf{1}_{\{0 < t_1 < \dots < t_n\}} = \mathbf{1}_{\bigcap_{i=1}^n \{w_i > 0\}} = \prod_{i=1}^n \mathbf{1}_{\{w_i > 0\}}.
     $$
     Hence, we obtain the following:
     $$
     \begin{aligned}
     &\int_{\mathbb{R}^n} g(t_1, t_2 - t_1, \dots, t_n - t_{n-1}) \lambda^n e^{-\lambda t_n} \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}\,dt_1 \dots dt_n \\
     &= \int_{\mathbb{R}^n} g(w_1, w_2, \dots, w_n) \lambda^n e^{-\lambda \sum_{i=1}^n w_i} \prod_{i=1}^n \mathbf{1}_{\{w_i > 0\}} |\det(J_\phi(w))|\,dw_1 \dots dw_n \\
     &= \int_{\mathbb{R}^n} g(w_1, w_2, \dots, w_n) \underbrace{\prod_{i=1}^n \left[ \lambda e^{-\lambda w_i} \mathbf{1}_{\{w_i > 0\}} \right]}_{(*)} dw_1 \dots dw_n.
     \end{aligned}
     $$

   - **(c)** Conclude.
     **Solution.**
     We recognize in $(*)$ the density of $n$ i.i.d. r.v. with $\mathcal{E}(\lambda)$ distribution.

2. Immediate from Proposition 1.19, iii).

The reverse of Proposition 2.18 is true and thus provides a way of constructing a homogeneous Poisson process (see Section II.4). Before proving this, we need the notion of order statistics.

### II.3. Conditional distribution of the arrival times

> **Definition 2.20 (Order statistic).**
> Let $Y_1, Y_2, \dots, Y_n$ be $n$ random variables. We say that $Y_{(1)}, Y_{(2)}, \dots, Y_{(n)}$ are the order statistics corresponding to $Y_1, \dots, Y_n$ if $Y_{(k)}$ is the $k$th smallest value among the $Y_i$'s. In other words
> $$
> Y_{(1)} \le Y_{(2)} \le \dots \le Y_{(n)}.
> $$

As an example, if $n = 3$, $Y_1 = 4$, $Y_2 = 5$ and $Y_3 = 1$, then $Y_{(1)} = 1$, $Y_{(2)} = 4$ and $Y_{(3)} = 5$.
Note that even if the $Y_i$ are independent and identically distributed, the order statistics are **not** independent, nor identically distributed.

<!-- page 33 -->
> **Lemma 2.21.**
> If $Y_1, \dots, Y_n$ are i.i.d. continuous random variables with density $f$, then the joint density of the order statistics $(Y_{(1)}, \dots, Y_{(n)})$ is given by
> $$
> (u_1, \dots, u_n) \mapsto n! \left( \prod_{i=1}^n f(u_i) \right) \mathbf{1}_{\{u_1 < \dots < u_n\}}.
> $$

*Proof of Lemma 2.21.* Let $B$ be a Borel set in $\mathbb{R}^n$. Let us denote $\mathcal{S}_n$ the set of permutations of $\{1, \dots, n\}$. Then
$$
\begin{aligned}
\mathbb{P}((Y_{(1)}, \dots, Y_{(n)}) \in B) &= \sum_{\sigma \in \mathcal{S}_n} \mathbb{P}(\{Y_{\sigma(1)} < Y_{\sigma(2)} < \dots < Y_{\sigma(n)}\} \cap \{(Y_{\sigma(1)}, \dots, Y_{\sigma(n)}) \in B\}) \\
&= \sum_{\sigma \in \mathcal{S}_n} \mathbb{E}\left[ \mathbf{1}_{\{(Y_{\sigma(1)}, \dots, Y_{\sigma(n)}) \in B\}} \mathbf{1}_{\{Y_{\sigma(1)} < \dots < Y_{\sigma(n)}\}} \right].
\end{aligned}
$$
Yet, if $\sigma \in \mathcal{S}_n$ is deterministic, then $Y_{\sigma(1)}, \dots, Y_{\sigma(n)}$ are i.i.d. with density $f$. Hence
$$
\mathbb{P}((Y_{(1)}, \dots, Y_{(n)}) \in B) = \sum_{\sigma \in \mathcal{S}_n} \int_{\mathbb{R}^n} \mathbf{1}_{\{(u_1, \dots, u_n) \in B\}} \mathbf{1}_{\{u_1 < \dots < u_n\}} f(u_1) \dots f(u_n)\,du_1 \dots du_n.
$$
Finally, since the cardinality of $\mathcal{S}_n$ equals $n!$:
$$
\mathbb{P}((Y_{(1)}, \dots, Y_{(n)}) \in B) = \int_{\mathbb{R}^n} \mathbf{1}_{\{(u_1, \dots, u_n) \in B\}} \underbrace{n! \left[ \prod_{i=1}^n f(u_i) \right] \mathbf{1}_{\{u_1 < \dots < u_n\}}}_{(*)} du_1 \dots du_n.
$$
This being true for any $B$, we deduce that $(*)$ is the joint density of the order statistics $(Y_{(1)}, \dots, Y_{(n)})$.

In the particular case of uniform random variables, we obtain the following property.

> **Property 2.22.**
> - **i)** If $U_1, \dots, U_n$ are i.i.d. r.v. with uniform distribution on $[0, t]$, then the joint density of the order statistics $U_{(1)}, U_{(2)}, \dots, U_{(n)}$ is given by
>   $$
>   (u_1, u_2, \dots, u_n) \mapsto \frac{n!}{t^n} \mathbf{1}_{\{0 < u_1 < u_2 < \dots < u_n < t\}}.
>   $$
> - **ii)** We deduce the following integral:
>   $$
>   \int_{\mathbb{R}^n} \mathbf{1}_{\{0 < u_1 < u_2 < \dots < u_n < t\}}\,du_1 \dots du_n = \frac{t^n}{n!}.
>   $$

*Proof.* Trivial from Lemma 2.21.

<!-- page 34 -->
The following result is very useful to prove the reverse of Proposition 2.18 and thus provides a way of sampling a homogeneous Poisson process. It can also be used to test for homogeneity of a Poisson process.

> **Proposition 2.23 (Conditional distribution).**
> Let $N$ be a homogeneous Poisson process with rate $\lambda > 0$ and fix $t > 0$. Let $n \in \mathbb{N}^*$. Given that $N_t = n$, the $n$ first arrival times $(T_1, \dots, T_n)$ have the same distribution as the order statistic corresponding to $n$ independent random variables uniformly distributed on the interval $[0, t]$, that is
> $$
> (T_1, T_2, \dots, T_n) \mid \{N_t = n\} \stackrel{(d)}{=} (U_{(1)}, U_{(2)}, \dots, U_{(n)}) \quad \text{where } U_1, \dots, U_n \stackrel{i.i.d.}{\sim} \mathcal{U}([0, t]),
> $$
> where $\stackrel{(d)}{=}$ means that both r.v. have the same distribution.

*Proof. In exercise. Let $n \in \mathbb{N}^*$ and $B$ be a Borel set.*

- **a)** Write the event $\{N_t = n\}$ in terms of the arrival times $(T_i)_i$.
  **Solution.**
  The event $\{N_t = n\}$ means that there are exactly $n$ events by time $t$. This is equivalent to say that the $n$th event occurred by time $t$, and the $(n + 1)$th after time $t$. Hence
  $$
  \{N_t = n\} = \{T_n \le t < T_{n+1}\}.
  $$

- **b)** Using the definition of a Poisson process, and Lemma 2.19, deduce that
  $$
  \mathbb{P}((T_1, \dots, T_n) \in B \mid N_t = n) = \frac{n!}{(\lambda t)^n e^{-\lambda t}} \int_{\mathbb{R}^{n+1}} \mathbf{1}_{\{(t_1, \dots, t_n) \in B\}} \lambda^{n+1} e^{-\lambda t_{n+1}} \mathbf{1}_{\{0 < t_1 < \dots < t_n \le t < t_{n+1}\}}\,dt_1 \dots dt_n dt_{n+1}.
  $$
  **Solution.**
  By definition of the conditional probability:
  $$
  \mathbb{P}((T_1, \dots, T_n) \in B \mid N_t = n) = \frac{\mathbb{P}(\{(T_1, \dots, T_n) \in B\} \cap \{N_t = n\})}{\mathbb{P}(N_t = n)}.
  $$
  Yet, $N_t \sim \mathcal{P}(\lambda t)$, hence $\mathbb{P}(N_t = n) = \frac{(\lambda t)^n}{n!} e^{-\lambda t}$.
  Moreover, $\{N_t = n\} = \{T_n \le t < T_{n+1}\}$. Thus,
  $$
  \mathbb{P}((T_1, \dots, T_n) \in B \mid N_t = n) = \frac{n!}{(\lambda t)^n e^{-\lambda t}} \mathbb{P}(\{(T_1, \dots, T_n) \in B\} \cap \{T_n \le t < T_{n+1}\}).
  $$
  Moreover, Lemma 2.19 provides the density of $(T_1, \dots, T_n, T_{n+1})$, directly leading to the desired formula.

<!-- page 35 -->
- **c)** Deduce that
  $$
  \mathbb{P}((T_1, \dots, T_n) \in B \mid N_t = n) = \int_{\mathbb{R}^n} \mathbf{1}_{\{(t_1, \dots, t_n) \in B\}} \frac{n!}{t^n} \mathbf{1}_{\{0 < t_1 < \dots < t_n < t\}}\,dt_1 \dots dt_n.
  $$
  **Solution.**
  By simplification of the $\lambda$ powers, and by the Fubini-Tonelli theorem, we can integrate first w.r.t $t_{n+1}$, and obtain
  $$
  \begin{aligned}
  \mathbb{P}((T_1, \dots, T_n) \in B \mid N_t = n) &= \frac{n!}{t^n e^{-\lambda t}} \int_{\mathbb{R}^n} \mathbf{1}_{\{(t_1, \dots, t_n) \in B\}} \underbrace{\left( \int_t^{+\infty} \lambda e^{-\lambda t_{n+1}}\,dt_{n+1} \right)}_{e^{-\lambda t}} \mathbf{1}_{\{0 < t_1 < \dots < t_n \le t\}}\,dt_1 \dots dt_n \\
  &= \int_{\mathbb{R}^n} \mathbf{1}_{\{(t_1, \dots, t_n) \in B\}} \frac{n!}{t^n} \mathbf{1}_{\{0 < t_1 < \dots < t_n < t\}}\,dt_1 \dots dt_n.
  \end{aligned}
  $$

- **d)** Conclude.
  **Solution.**
  We recognize the density of the order statistic of $n$ i.i.d. uniformly distributed r.v. on $[0, t]$ given in Property 2.22.

### II.4. Construction of a homogeneous Poisson process

Let us now prove the reverse of Proposition 2.18, which provides an alternative definition of a homogeneous Poisson process.

> **Theorem 2.24.**
> Let $(W_i)_{i \ge 1}$ be a sequence of i.i.d. random variables with distribution $\mathcal{E}(\lambda)$. Then the point process with arrival times defined for all $n \ge 1$ by $T_n = \sum_{i=1}^n W_i$, is a homogeneous Poisson process with rate $\lambda$.

This result relies on the following Lemmas.

> **Lemma 2.25.**
> Under the assumptions of Theorem 2.24, for all $t$ in $\mathbb{R}_+$ and all $n$ in $\mathbb{N}$,
> $$
> \mathbb{P}(N_t = n) = \frac{(\lambda t)^n}{n!} e^{-\lambda t}.
> $$
> In particular, for all $t > 0$, $N_t \sim \mathcal{P}(\lambda t)$.

*Proof.* Recall that $\{N_t = n\} = \{T_n \le t < T_{n+1}\}$, where
$$
\begin{cases}
T_n = \sum_{i=1}^n W_i \sim \Gamma(n, \lambda), \\
T_{n+1} = T_n + W_{n+1} \quad \text{with } W_{n+1} \perp\!\!\!\perp T_n, \text{ and } W_{n+1} \sim \mathcal{E}(\lambda).
\end{cases}
$$

<!-- page 36 -->
1. Justify that
   $$
   \mathbb{P}(N_t = n) = \int_{\mathbb{R}^2} \mathbf{1}_{\{s \le t\}} \mathbf{1}_{\{s + w > t\}} \frac{(\lambda s)^{n-1}}{(n - 1)!} \lambda e^{-\lambda s} \mathbf{1}_{\{s > 0\}} \lambda e^{-\lambda w} \mathbf{1}_{\{w > 0\}}\,ds\,dw.
   $$
   **Solution.**
   From Definition 1.17, we have
   $$
   \begin{aligned}
   \mathbb{P}(N_t = n) &= \mathbb{P}(\{T_n \le t\} \cap \{T_n + W_{n+1} > t\}) \\
   &= \int_{\mathbb{R}^2} \mathbf{1}_{\{s \le t\}} \mathbf{1}_{\{s + w > t\}} \underbrace{\frac{(\lambda s)^{n-1}}{(n - 1)!} \lambda e^{-\lambda s} \mathbf{1}_{\{s > 0\}}}_{\text{density of } T_n} \lambda e^{-\lambda w} \mathbf{1}_{\{w > 0\}}\,ds\,dw.
   \end{aligned}
   $$

2. Compute $\mathbb{P}(N_t = n)$.
   **Solution.**
   Hence by the Fubini-Tonelli theorem:
   $$
   \begin{aligned}
   \mathbb{P}(N_t = n) &= \int_0^t \frac{\lambda^n s^{n-1}}{(n - 1)!} e^{-\lambda s} \underbrace{\left( \int_{t-s}^{+\infty} \lambda e^{-\lambda w}\,dw \right)}_{e^{-\lambda(t-s)}} ds \\
   &= \frac{\lambda^n}{(n - 1)!} e^{-\lambda t} \underbrace{\int_0^t s^{n-1}\,ds}_{t^n / n} = \frac{(\lambda t)^n}{n!} e^{-\lambda t}.
   \end{aligned}
   $$

In order to prove the stationarity and the independence of the increments, we also need the joint distribution of the arrival times.

> **Lemma 2.26.**
> Under the assumptions of Theorem 2.24, the arrival times $(T_1, \dots, T_n)$ have the following density w.r.t. the Lebesgue measure:
> $$
> (t_1, \dots, t_n) \mapsto \lambda^n e^{-\lambda t_n} \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}.
> $$

Note that we recover the same density as in Lemma 2.19 under the homogeneous Poisson process assumption (phew!).

*Proof.* "Same" change of variables as in the proof of Proposition 2.18.

We may now prove Theorem 2.24.

<!-- page 37 -->
*Proof of Theorem 2.24.* Notice that Lemmas 2.25 and 2.26 jointly imply that the conditional distribution of $(T_1, \dots, T_n) \mid N_t = n$ is the same as the joint distribution of the order statistic $U_{(1)} < \dots < U_{(n)}$ of i.i.d. r.v. $(U_1, \dots, U_n)$ with distribution $\mathcal{U}([0, t])$ (exactly same proof as Proposition 2.23).

Let $k \in \mathbb{N}^*$ and $0 = t_0 < t_1 < \dots < t_k = t$ be a subdivision of $[0, t]$. Conditionally on $\{N_t = n\}$, since $(T_1, \dots, T_n) \stackrel{(d)}{=} (U_{(1)}, \dots, U_{(n)})$:
$$
N_{t_j} - N_{t_{j-1}} = \sum_{i=1}^n \mathbf{1}_{\{t_{j-1} < T_i \le t_j\}} \stackrel{(d)}{=} \sum_{i=1}^n \mathbf{1}_{\{t_{j-1} < U_{(i)} \le t_j\}} = \sum_{i=1}^n \mathbf{1}_{\{t_{j-1} < U_i \le t_j\}}.
$$
Then $(N_{t_1}, N_{t_2} - N_{t_1}, \dots, N_{t_k} - N_{t_{k-1}})$ has a multinomial distribution with parameters $(n, (p_j)_{1 \le j \le k})$ where
$$
p_j = \mathbb{P}(t_{j-1} < U_1 \le t_j) = \frac{t_j - t_{j-1}}{t}, \quad \text{since } U_1 \sim \mathcal{U}([0, t]).
$$
Hence, for all integers $n_1, \dots, n_k$ such that $\sum_{j=1}^k n_j = n$:
$$
\begin{aligned}
\mathbb{P}(N_{t_1} = n_1, N_{t_2} - N_{t_1} = n_2, \dots, N_{t_k} - N_{t_{k-1}} = n_k \mid N_t = n) &= \frac{n!}{n_1! \dots n_k!} p_1^{n_1} \dots p_k^{n_k} \\
&= \frac{n!}{t^n} \prod_{j=1}^k \frac{(t_j - t_{j-1})^{n_j}}{n_j!}.
\end{aligned}
$$
We thus obtain the unconditional probability:
$$
\begin{aligned}
&\mathbb{P}(N_{t_1} = n_1, N_{t_2} - N_{t_1} = n_2, \dots, N_{t_k} - N_{t_{k-1}} = n_k) \\
&= \mathbb{P}(N_{t_1} = n_1, N_{t_2} - N_{t_1} = n_2, \dots, N_{t_k} - N_{t_{k-1}} = n_k \mid N_t = n) \mathbb{P}(N_t = n) \\
&= \left( \frac{n!}{t^n} \prod_{j=1}^k \frac{(t_j - t_{j-1})^{n_j}}{n_j!} \right) \times \frac{(\lambda t)^n e^{-\lambda t}}{n!} \\
&= \prod_{j=1}^k \frac{[\lambda(t_j - t_{j-1})]^{n_j}}{n_j!} e^{-\lambda(t_j - t_{j-1})},
\end{aligned}
$$
since $\lambda^n = \lambda^{\sum_{j=1}^k n_j} = \prod_{j=1}^k \lambda^{n_j}$ and $e^{-\lambda t} = e^{-\lambda \sum_{j=1}^k (t_j - t_{j-1})} = \prod_{j=1}^k e^{-\lambda(t_j - t_{j-1})}$.

We recognize the joint distribution of $k$ independent r.v. with respective distribution $\mathcal{P}(\lambda(t_j - t_{j-1}))$. This directly implies points (2) (independence of the increments) and (3) (the number of points in an interval $(t_{j-1}, t_j]$ has a Poisson distribution with parameter $\lambda$ times the length $(t_j - t_{j-1})$ of the interval) in Definition 2.11.

A possible generalization of Poisson processes is to consider a point process for which the times between successive events are independent and identically distributed (but not necessarily exponentially distributed). Such a point process is called **renewal process**.

<!-- page 38 -->
## III. Divisibility of homogeneous Poisson processes

### III.1. Two states divisibility

> **Proposition 2.27.**
> Consider a homogeneous Poisson process $N$ with rate $\lambda > 0$ and assume that each time an event occurs it is classified as either a type I with probability $p$ or a type II event with probability $(1 - p)$ independently of all other events. Let $N^I$ and $N^{II}$ respectively denote the point processes corresponding to type I and type II events. Then:
> 1. For all $t > 0$, the total number of events occurring in $[0, t)$ equals $N_t = N_t^I + N_t^{II}$.
> 2. $N^I$ and $N^{II}$ are both Poisson processes having respective rates $p\lambda$ and $(1 - p)\lambda$.
> 3. Furthermore, the two processes $N^I$ and $N^{II}$ are independent.

This result relies on the following key point. The classification of each point being independent of when it occurs, we can say that the conditional distribution of $N_t^I$ given that $N_t = n$ has a binomial distribution with parameters $(n, p)$, that is
$$
N_t^I \mid \{N_t = n\} \sim \mathcal{B}(n, p). \tag{2.5}
$$

*Proof.* Point 1 is trivial. Let us prove 2 and 3.

- **The independence of the increments of $N^I$ and $N^{II}$** comes from the fact that $N$ is a Poisson process. Hence the number of points in disjoint intervals are independent. Moreover, the type is decided independently of the point process.
- **The stationarity of the increments of $N^I$ and $N^{II}$** also comes from the fact that $N$ is a Poisson process. Indeed, the number of points of $N$ in two intervals with same length has the same distribution. Moreover, once again, the classification does not depend on the times the event occurs.
- **Poisson distribution and $N^I \perp\!\!\!\perp N^{II}$:** Let $t > 0$ and $k, l \in \mathbb{N}$. Then,
  $$
  \begin{aligned}
  \mathbb{P}(N_t^I = k, N_t^{II} = l) &= \mathbb{P}(N_t^I = k, N_t = k + l) \\
  &= \mathbb{P}(N_t^I = k \mid N_t = k + l) \mathbb{P}(N_t = k + l).
  \end{aligned}
  $$
  Yet, from the key point stated above, the conditional distribution of $N_t^I$ given $\{N_t = k + l\}$ is a binomial $\mathcal{B}(k + l, p)$. Moreover, $N_t \sim \mathcal{P}(\lambda t)$. Hence
  $$
  \begin{aligned}
  \mathbb{P}(N_t^I = k, N_t^{II} = l) &= \frac{(k + l)!}{k!\,l!} p^k (1 - p)^l \times \frac{(\lambda t)^{k+l}}{(k + l)!} e^{-\lambda t} \\
  &= \underbrace{\frac{(p\lambda t)^k}{k!} e^{-p\lambda t}}_{\mathbb{P}(N_t^I = k)} \times \underbrace{\frac{((1 - p)\lambda t)^l}{l!} e^{-(1-p)\lambda t}}_{\mathbb{P}(N_t^{II} = l)}.
  \end{aligned}
  $$
  We deduce that:
  $$
  \begin{cases}
  N_t^I \sim \mathcal{P}(p\lambda t) \\
  N_t^{II} \sim \mathcal{P}((1 - p)\lambda t) \\
  N_t^I \perp\!\!\!\perp N_t^{II}
  \end{cases}
  $$

<!-- page 39 -->
> **Example 2.28.**
> Clients of an insurance company call at a Poisson process rate of ten per week. In 10% of the cases, it is to sign up a new policy, and in the remaining 90%, to report a claim. What is the probability that there will be no new customers in February?
> **Solution.**
> The number of new customers during February follows a Poisson distribution with parameter $4 \text{ (weeks)} \times 10\,(\lambda) \times 1/10\,(p) = 4$, so the answer is $e^{-4}$.

The reverse of Proposition 2.27 is true.

> **Proposition 2.29.**
> If $N^I$ and $N^{II}$ are two independent Poisson processes with respective rates $\lambda^I$ and $\lambda^{II}$, then the counting process $N = (N_t)_{t \ge 0}$ defined for all $t$ in $\mathbb{R}_+$ by
> $$
> N_t = N_t^I + N_t^{II}
> $$
> is also a Poisson process with rate $\lambda^I + \lambda^{II}$.

*Proof. Facultative homework.*

### III.2. Infinite divisibility

The following generalization of Proposition 2.27 can be easily proved by mathematical induction.

> **Proposition 2.30.**
> Let $N$ be a homogeneous Poisson process with intensity $\lambda$. Assume each time an event occurs is associated to one of $M$ types in $\{1, \dots, M\}$ with probability $1/M$ independently of all other events. For all $m$ in $\{1, \dots, M\}$, denote $N^m$ the point process associated to events of type $m$. Then $N^1, \dots, N^M$ are $M$ i.i.d. homogeneous Poisson processes with intensity $\lambda/M$.

> **Definition 2.31.**
> A random variable $X$ is said to be **infinitely divisible** if for all $M$ in $\mathbb{N}^*$, there exists $M$ i.i.d. random variables $Y_1, \dots, Y_M$ such that $X$ has the same distribution as $Y_1 + \dots + Y_M$.

> **Property 2.32.**
> Random variables with distribution $\mathcal{P}(\lambda)$, $\mathcal{N}(\mu, \sigma^2)$, $\mathcal{E}(\lambda)$ or $\Gamma(n, \lambda)$ are infinitely divisible.

> **Property 2.33.**
> A homogeneous Poisson process is infinitely divisible.

<!-- page 40 -->

<!-- page 41 -->
# Chapter 3. Statistics for homogeneous Poisson processes

In this part, we do some statistical inference for homogeneous Poisson processes. We observe a random point process $N = (N_t)_{t \in \mathbb{R}_+}$ with arrival times $(T_n)_{n \in \mathbb{N}^*}$.

We assume that $N$ is a homogeneous Poisson process with (unknown) intensity $\lambda > 0$, and our goal is to estimate $\lambda$ and to construct statistical tests on $\lambda$.

As for the simulation of a point process, there are two distinct observation frameworks:
- **Either the observation interval $[0, T]$ is fixed**, and we observe the process on a fixed interval, that is we observe $(N_t)_{t \in [0, T]}$.

```
    Nt
     ^
     |                                                         •
     |                                               •
     |                                         •
     |                                   •
     |                             •
     |                       •
   0 +-----------•-----------•-----•-----•-------------•-------+-------> t
     0           T1          T2    T3    T4           TNT      T
```
In this case, the number of events occurring in this interval $N_T$ is random.

- **Either $n$ is fixed**, and we observe the process up to the $n$th event, that is we observe $(T_1, \dots, T_n)$.

```
   0 +-----------•-----------•-----•-----•-------------···-----•-------> t
     0           T1          T2    T3    T4                    Tn
```
In this case, the observation interval/window $[0, T_n]$ is random.

## II. Fixed window, random number of points

In this section $T$ is fixed, and we observe a homogeneous Poisson process $N$ with (unknown) intensity $\lambda > 0$ on the interval $[0, T]$. In particular, the number of "events" that have occurred by time $T$, that is $N_T$, is random.

<!-- page 42 -->
### II.1. Maximum Likelihood Estimator

**Heuristic:** Recall that, if $X$ is a random variable with density $f_\theta$ and c.d.f. $F_\theta$. Then, its likelihood in $\theta$ equals
$$
L(X, \theta) = f_\theta(X),
$$
where
$$
f_\theta(x) = F'_\theta(x) = \lim_{h \to 0} \frac{F_\theta(x + h) - F_\theta(x)}{h} = \lim_{h \to 0} \frac{\mathbb{P}_\theta(X \in (x; x + h])}{h}.
$$

**In the point process case:** Recall that in the "fixed-$T$" case, $N_T = \sum_{i=1}^{+\infty} \mathbf{1}_{\{T_i \le T\}}$ is random. Similarly to the proof of Lemma 2.19, consider an increasing sequence $(t_i)_{i \in \mathbb{N}^*}$ in $\mathbb{R}_+$, and denote
$$
k_t = \sum_{i=1}^{+\infty} \mathbf{1}_{\{t_i \le T\}} \quad \text{s.t.} \quad 0 =: t_0 < t_1 < \dots < t_{k_t} < T < t_{k_t + 1} < \dots.
$$
Hence, $k_t$ counts the number of $t_i$'s that belong to $[0, T]$.
Introduce for $h = (h_1, \dots, h_n)$ small enough such that all $(t_i, t_i + h_i]$ are disjoint (and $t_{k_t} + h_{k_t} \le T$):
$$
B_h = \{1 \text{ point in each } (t_i, t_i + h_i] \text{ and } 0 \text{ points elsewhere in } [0, T]\}.
$$

```
   T1              T2                         Tkt
---|---|-----------|---|----------------------|---|----------]----------> t
   t1  t1+h1       t2  t2+h2                  tkt tkt+hkt    T
```

The main difference with $A_h$ in Lemma 2.19 is that now $n$ is not fixed anymore, and we need to consider $k_t$ which is the number of $(t_i)$ which belong to the observation interval $[0, T]$.

Heuristically, as in the proof of Lemma 2.19,
$$
f((t_i)_{i \in \mathbb{N}^*} \cap [0, T]) = \lim_{\forall i, \, h_i \to 0} \frac{\mathbb{P}(B_h)}{h_1 \dots h_{k_t}}.
$$
Since the intervals appearing in $B_h$ are disjoint, by independence of the increments:
$$
\begin{aligned}
\mathbb{P}(B_h) &= \left( \prod_{i=1}^{k_t} \mathbb{P}(1 \text{ point in } (t_i, t_i + h_i]) \right) \times \mathbb{P}(0 \text{ points elsewhere}) \\
&\stackrel{(*)}{=} \left( \prod_{i=1}^{k_t} \lambda h_i e^{-\lambda h_i} \right) \times e^{-\lambda(T - \sum_{i=1}^{k_t} h_i)} \\
&= \lambda^{k_t} h_1 \dots h_{k_t} e^{-\lambda T},
\end{aligned}
$$
with $(*)$ because "elsewhere" is of length $T - \sum_{i=1}^{k_t} h_i$.

Finally,
$$
f((t_i)_{i \in \mathbb{N}^*} \cap [0, T]) = \lim_{\forall i, \, h_i \to 0} \lambda^{k_t} e^{-\lambda T} = \lambda^{k_t} e^{-\lambda T}.
$$

<!-- page 43 -->
> **Proposition 3.1.**
> Let $N$ be a homogeneous Poisson process with (unknown) intensity $\lambda$. Consider the number $N_T$ of "events" that have occurred in $[0, T]$, and denote $(T_1, \dots, T_{N_T})$ the corresponding arrival times. Then likelihood of the trajectory of $N$ is equal to
> $$
> L((N_t)_{t \in [0, T]}; \lambda) = \lambda^{N_T} e^{-\lambda T}.
> $$
> Hence, the Maximum Likelihood Estimator (MLE) of $\lambda$ is
> $$
> \hat{\lambda}_T = \frac{N_T}{T}.
> $$

Notice that the likelihood can be expressed as follows:
$$
L((N_t)_{t \in [0, T]}; \lambda) = N_T! \left( \prod_{i=1}^{N_T} \frac{1}{T} \mathbf{1}_{\{T_i \in [0, T]\}} \right) \mathbf{1}_{\{T_1 < \dots < T_{N_T}\}} \times \frac{(\lambda T)^{N_T} e^{-\lambda T}}{N_T!}.
$$
We recognize the product of:
- the (conditional) likelihood of $(T_1, \dots, T_{N_T}) \mid N_T$, that is, according to Proposition 2.23, the likelihood of the order statistic corresponding to $N_T$ independent random variables uniformly distributed on the interval $[0, T]$ given $N_T$,
- and the likelihood of $N_T$ which has a Poisson distribution $\mathcal{P}(\lambda T)$.

*Proof. In exercise. Compute the MLE of $\lambda$.*
The expression of the likelihood is admitted here (c.f. heuristic above).
**Solution.**
We deduce that the log-likelihood equals
$$
\ell((N_t)_{t \in [0, T]}; \lambda) = N_T \ln(\lambda) - \lambda T.
$$
Hence
$$
\frac{\partial \ell}{\partial \lambda}((N_t)_{t \in [0, T]}; \lambda) = \frac{N_T}{\lambda} - T = 0 \iff \lambda = \frac{N_T}{T}.
$$
Moreover, $N_T / T$ is a maximum since the log-likelihood function is concave. Indeed,
$$
\frac{\partial^2 \ell}{\partial \lambda^2}((N_t)_{t \in [0, T]}; \lambda) = -\frac{N_T}{\lambda^2} \le 0 \quad \text{a.s.}
$$

### II.2. Non-asymptotic properties of the MLE

> **Proposition 3.2.**
> Let $(N_t)_{t \in [0, T]}$ be a homogeneous Poisson process with rate $\lambda > 0$ observed on $[0, T]$ and consider the MLE $\hat{\lambda}_T = N_T / T$. Since $N_T \sim \mathcal{P}(\lambda T)$,
> $$
> \mathbb{E}[\hat{\lambda}_T] = \lambda \quad \text{and} \quad \text{Var}(\hat{\lambda}_T) = \frac{\lambda}{T}.
> $$
> We deduce that $\hat{\lambda}_T$ is an unbiased and efficient estimator of $\lambda$.

<!-- page 44 -->
*Proof. In exercise.*
**Solution.**
Since $N_T \sim \mathcal{P}(\lambda T)$, $\mathbb{E}[N_T] = \lambda T$ and $\text{Var}(N_T) = \lambda T$. Hence, on the one hand, by linearity of the expectation,
$$
\mathbb{E}[\hat{\lambda}_T] = \frac{\mathbb{E}[N_T]}{T} = \lambda.
$$
We deduce that $\hat{\lambda}_T$ is an unbiased estimator of $\lambda$. On the other hand, by property of the variance,
$$
\text{Var}(\hat{\lambda}_T) = \frac{\text{Var}(N_T)}{T^2} = \frac{\lambda T}{T^2} = \frac{\lambda}{T}.
$$
Moreover, one can compute the Fisher information:
$$
I_T(\lambda) = \mathbb{E}_\lambda\left[ -\frac{\partial^2}{\partial \lambda^2} \ell((N_t)_{t \in [0, T]}; \lambda) \right] = \mathbb{E}_\lambda\left[ \frac{N_T}{\lambda^2} \right] = \frac{T}{\lambda}.
$$
Thus the Cramér-Rao bound equals
$$
B_T(\lambda) = \frac{1}{I_T(\lambda)} = \frac{\lambda}{T} = \text{Var}(\hat{\lambda}_T).
$$
We deduce that $\hat{\lambda}_T$ is efficient.

Moreover, to construct statistical tests with uniform asymptotic level, we need the following lemma.

> **Lemma 3.3.**
> Let $N$ be a homogeneous Poisson process with (unknown) intensity $\lambda > 0$ observed on the interval $[0, T]$ and consider the MLE $\hat{\lambda}_T = N_T / T$ of $\lambda$.
> Then, for all (fixed) $x \ge 0$,
> $$
> \lambda \mapsto \mathbb{P}_\lambda(\hat{\lambda}_T \ge x)
> $$
> is a non-decreasing function on $\mathbb{R}_+^*$.

*Proof. In exercise. Fix $x \ge 0$ and denote $g : \lambda > 0 \mapsto \mathbb{P}_\lambda(\hat{\lambda}_T \ge x)$.*

1. Assume $x = 0$. Compute $g(\lambda)$ for all $\lambda > 0$ and conclude in this case.
   **Solution.**
   If $x = 0$, $g$ is constant equal to 1, since for all $\lambda > 0$, $\mathbb{P}(\hat{\lambda}_T \ge 0) = \mathbb{P}(N_T \ge 0) = 1$.

2. Now assume $x > 0$. Justify that for all $\lambda > 0$,
   $$
   g(\lambda) = \sum_{k=\lceil xT \rceil}^{+\infty} \frac{(\lambda T)^k}{k!} e^{-\lambda T},
   $$
   where $\lceil \cdot \rceil$ denotes the ceiling function.
   **Solution.**
   Note that
   $$
   g(\lambda) = \mathbb{P}(\hat{\lambda}_T \ge x) = \mathbb{P}(N_T \ge xT) = \mathbb{P}(N_T \ge \lceil xT \rceil),
   $$
   since $N_T$ is an integer. The expression of $g(\lambda)$ is then immediate since $N_T \sim \mathcal{P}(\lambda T)$.

<!-- page 45 -->
3. Using the method of differences, that is telescoping series, deduce that
   $$
   g'(\lambda) = \frac{(\lambda T)^{\lceil xT \rceil - 1}}{(\lceil xT \rceil - 1)!} T e^{-\lambda T}.
   $$
   **Solution.**
   We obtain:
   $$
   \begin{aligned}
   g'(\lambda) &= \left[ \sum_{k=\lceil xT \rceil}^{+\infty} \frac{k \lambda^{k-1} T^k}{k!} e^{-\lambda T} \right] - \left[ \sum_{k=\lceil xT \rceil}^{+\infty} \frac{(\lambda T)^k}{k!} T e^{-\lambda T} \right] \\
   &= \left( \sum_{k=\lceil xT \rceil}^{+\infty} \frac{(\lambda T)^{k-1}}{(k-1)!} - \sum_{k=\lceil xT \rceil}^{+\infty} \frac{(\lambda T)^k}{k!} \right) T e^{-\lambda T} \\
   &\stackrel{(*)}{=} \left( \sum_{k=\lceil xT \rceil - 1}^{+\infty} \frac{(\lambda T)^k}{k!} - \sum_{k=\lceil xT \rceil}^{+\infty} \frac{(\lambda T)^k}{k!} \right) T e^{-\lambda T} \\
   &\stackrel{(\dagger)}{=} \frac{(\lambda T)^{\lceil xT \rceil - 1}}{(\lceil xT \rceil - 1)!} T e^{-\lambda T},
   \end{aligned}
   $$
   with $(*)$ by change of variable $k \leftarrow k - 1$ (since $\lceil xT \rceil \ge 1$), and $(\dagger)$ by the difference method.

4. Conclude.
   **Solution.**
   We deduce that $g'(\lambda) > 0$ for all $\lambda > 0$, that is $g$ is increasing when $x > 0$.

### II.3. Asymptotic distribution of the MLE

> **Theorem 3.4.**
> Let $(N_t)_{t \in [0, T]}$ be a Poisson process with rate $\lambda > 0$ observed on $[0, T]$ and consider the MLE $\hat{\lambda}_T = N_T / T$.
> 1. **LLN-type result:** $\hat{\lambda}_T$ is a consistent estimator of $\lambda$, that is
>    $$
>    \frac{N_T}{T} \xrightarrow[T \to +\infty]{\mathbb{P}} \lambda.
>    $$
> 2. **CLT-type result:** Moreover, it is asymptotically gaussian, that is
>    $$
>    \sqrt{T}\left( \frac{N_T}{T} - \lambda \right) \xrightarrow[T \to +\infty]{\mathcal{L}} \mathcal{N}(0, \lambda).
>    $$

Note that it is also possible to prove that $N_T \xrightarrow[T \to +\infty]{\text{a.s.}} +\infty$ and deduce the "strong" version of the LLN-type result:
$$
\frac{N_T}{T} \xrightarrow[T \to +\infty]{\text{a.s.}} \lambda.
$$
This result is based on the following inequality: $\frac{T_{N_T}}{N_T} \le \frac{T}{N_T} < \frac{T_{N_T + 1}}{N_T}$.

*Proof. In exercise.*

<!-- page 46 -->
1. Note that we cannot apply the classical Law of Large Numbers since the index in the limit is not a countable sequence. Apply Chebychev's inequality to prove the consistency of $\hat{\lambda}_T$.
   **Solution.**
   Since $\mathbb{E}[\hat{\lambda}_T] = \lambda$ and $\text{Var}(\hat{\lambda}_T) = \lambda/T < +\infty$, for all $\varepsilon > 0$:
   $$
   \mathbb{P}(|\hat{\lambda}_T - \lambda| \ge \varepsilon) \le \frac{\text{Var}(\hat{\lambda}_T)}{\varepsilon^2} = \frac{\lambda}{T\varepsilon^2} \xrightarrow[T \to +\infty]{} 0.
   $$

2. Using Proposition 1.2, prove that the characteristic function of $Z_T := \sqrt{T}(N_T / T - \lambda)$ equals for all $u \in \mathbb{R}$:
   $$
   \mathbb{E}\left[e^{iuZ_T}\right] = \exp\left( \lambda T \left[ e^{iu/\sqrt{T}} - 1 - \frac{iu}{\sqrt{T}} \right] \right).
   $$
   **Solution.**
   Let $u \in \mathbb{R}$.
   $$
   \mathbb{E}\left[e^{iuZ_T}\right] = \mathbb{E}\left[ e^{iu\sqrt{T}(N_T / T - \lambda)} \right] = \mathbb{E}\left[ e^{i\frac{u}{\sqrt{T}} N_T} \right] e^{-iu\lambda\sqrt{T}}.
   $$
   Yet $N_T \sim \mathcal{P}(\lambda T)$. Hence, by Proposition 1.2, its characteristic function computed in $u/\sqrt{T}$ equals:
   $$
   \mathbb{E}\left[ e^{i\frac{u}{\sqrt{T}} N_T} \right] = \exp\left( \lambda T \left[ e^{iu/\sqrt{T}} - 1 \right] \right).
   $$
   Hence, we directly obtain:
   $$
   \mathbb{E}\left[e^{iuZ_T}\right] = \exp\left( \lambda T \left[ e^{iu/\sqrt{T}} - 1 - \frac{iu}{\sqrt{T}} \right] \right).
   $$

3. Using Taylor's expansion, prove that
   $$
   \lim_{T \to +\infty} \mathbb{E}\left[e^{iuZ_T}\right] = e^{-\frac{\lambda u^2}{2}}.
   $$
   **Solution.**
   By Taylor's expansion, as $x \to 0$:
   $$
   \exp(x) = 1 + x + \frac{x^2}{2} + o(x^2).
   $$
   Hence as $T \to +\infty$:
   $$
   \mathbb{E}\left[e^{iuZ_T}\right] = \exp\left( \lambda T \left[ \frac{-u^2}{2T} + o\left(\frac{1}{T}\right) \right] \right) = \exp\left( -\frac{\lambda u^2}{2} + o(1) \right) \xrightarrow[T \to +\infty]{} \exp\left(-\frac{\lambda u^2}{2}\right).
   $$

4. Deduce the asymptotic distribution of $Z_T$. *Hint: Recall that the characteristic function of a gaussian r.v. $Z$ with distribution $\mathcal{N}(m, \sigma^2)$ equals for all $u \in \mathbb{R}$, $\mathbb{E}[e^{iuZ}] = \exp(imu - \frac{u^2\sigma^2}{2})$.*
   **Solution.**
   We recognize the characteristic function of a $\mathcal{N}(0, \lambda)$, which ends the proof of Theorem 3.4.

<!-- page 47 -->
### II.4. Confidence intervals and statistical tests

Fix $T > 0$ and let $\alpha \in (0, 1)$. We observe a homogeneous Poisson process with unknown intensity $\lambda > 0$ on $[0, T]$ and we aim at constructing confidence intervals for $\lambda$ with confidence level $1 - \alpha$, or at testing, for instance,
$$
H_0 : \lambda \le \lambda_0 \quad \text{against} \quad H_1 : \lambda > \lambda_0, \quad \text{at level } \alpha.
$$

From Theorem 3.4 and Slutsky's Lemma, for all $\alpha \in (0, 1)$, we construct an asymptotic confidence interval for $\lambda$ with confidence level $1 - \alpha$, that is
$$
\left[ \frac{N_T}{T} - \frac{\sqrt{N_T}}{T} z_{1-\alpha/2} \; ; \; \frac{N_T}{T} + \frac{\sqrt{N_T}}{T} z_{1-\alpha/2} \right],
$$
where $z_{1-\alpha/2}$ denotes the $(1 - \alpha/2)$-quantile of a standard gaussian distribution $\mathcal{N}(0, 1)$.

*Proof. In exercise.*
**Solution.**
By Slutsky's Lemma:
$$
\frac{T}{\sqrt{N_T}} \left( \frac{N_T}{T} - \lambda \right) = \underbrace{\sqrt{\frac{\lambda}{N_T / T}}}_{(a)} \times \underbrace{\sqrt{\frac{T}{\lambda}} \left( \frac{N_T}{T} - \lambda \right)}_{(b)} \xrightarrow[T \to +\infty]{\mathcal{L}} \mathcal{N}(0, 1),
$$
Since by the LLN-type result in Theorem 3.4 and by continuity of $x \mapsto \sqrt{\lambda/x}$ on $\mathbb{R}_+^*$:
$$
(a) \xrightarrow[T \to +\infty]{\mathbb{P}} 1,
$$
and by the CLT-type result in Theorem 3.4:
$$
(b) \xrightarrow[T \to +\infty]{\mathcal{L}} \mathcal{N}(0, 1).
$$
We deduce that
$$
\begin{aligned}
\mathbb{P}\left( \frac{N_T}{T} - \frac{\sqrt{N_T}}{T} z_{1-\alpha/2} \le \lambda \le \frac{N_T}{T} + \frac{\sqrt{N_T}}{T} z_{1-\alpha/2} \right) &= \mathbb{P}\left( -z_{1-\alpha/2} \le \frac{T}{\sqrt{N_T}}\left(\frac{N_T}{T} - \lambda\right) \le z_{1-\alpha/2} \right) \\
&\xrightarrow[T \to +\infty]{} \mathbb{P}(-z_{1-\alpha/2} \le Z \le z_{1-\alpha/2}) = 1 - \alpha,
\end{aligned}
$$
where $Z$ denotes a standard gaussian r.v.

<!-- page 48 -->
Moreover, one can construct different tests of
$$
H_0 : \lambda \le \lambda_0 \quad \text{against} \quad H_1 : \lambda > \lambda_0, \quad \text{at level } \alpha.
$$

- **Test 1:** Reject $H_0$ if $\hat{\lambda}_T \ge \lambda_0 + \sqrt{\frac{\lambda_0}{T}} z_{1-\alpha}$, with
  $$
  \begin{cases}
  \text{test statistic: } Z_T^{(1)} = \sqrt{T} \frac{\hat{\lambda}_T - \lambda_0}{\sqrt{\lambda_0}}, \\
  \text{p-value: } 1 - \Phi(Z_T^{(1),\text{obs}}),
  \end{cases}
  $$

- **Test 2:** Reject $H_0$ if $\hat{\lambda}_T \ge \lambda_0 + \sqrt{\frac{\hat{\lambda}_T}{T}} z_{1-\alpha}$, with
  $$
  \begin{cases}
  \text{test statistic: } Z_T^{(2)} = \sqrt{T} \frac{\hat{\lambda}_T - \lambda_0}{\sqrt{\hat{\lambda}_T}}, \\
  \text{p-value: } 1 - \Phi(Z_T^{(2),\text{obs}}),
  \end{cases}
  $$

where $z_{1-\alpha}$ and $\Phi$ respectively denote the $(1 - \alpha)$-quantile and the c.d.f. of a standard gaussian distribution $\mathcal{N}(0, 1)$.

One can prove that both Test 1 and Test 2 are (uniformly) asymptotically of level $\alpha$, that is
$$
\limsup_{T \to +\infty} \sup_{\lambda \le \lambda_0} \left\{ \mathbb{P}_\lambda(\text{reject } H_0) \right\} \le \alpha.
$$
This result relies on Lemma 3.3.

*Proof. Facultative homework.*

<!-- page 49 -->
## III. Fixed number of points, random observation time

In this section $n$ is fixed, and we observe a homogeneous Poisson process $N$ with (unknown) intensity $\lambda > 0$ up to the $n$th arrival times. In particular, we observe the $n$ first arrival times $(T_1, \dots, T_n)$ or equivalently the $n$ first interarrival times $(W_1, \dots, W_n)$. Then, according to Proposition 2.18, the $(W_i)_{i \ge 1}$ are i.i.d. with distribution $\mathcal{E}(\lambda)$ which makes the study easier. In particular, $T_n$ is $\Gamma(n, \lambda)$-distributed.

### III.1. Maximum Likelihood Estimator (MLE)

> **Proposition 3.5.**
> Let $N$ be a homogeneous Poisson process with (unknown) intensity $\lambda > 0$. Consider the $n$ first arrival times $(T_1, \dots, T_n)$ of $N$. Then likelihood of the trajectory of $N$ is equal to
> $$
> L((T_1, \dots, T_n); \lambda) = \lambda^n e^{-\lambda T_n}.
> $$
> Hence, the Maximum Likelihood Estimator (MLE) of $\lambda$ is
> $$
> \hat{\lambda}_n = \frac{n}{T_n}.
> $$

*Proof. In exercise.*

1. Recall the joint density of $(T_1, \dots, T_n)$, and deduce its likelihood.
   **Solution.**
   According to Lemma 2.19, we know that the density of $(T_1, \dots, T_n)$ is
   $$
   (t_1, \dots, t_n) \mapsto \lambda^n e^{-\lambda t_n} \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}.
   $$
   Hence its likelihood equals
   $$
   L((T_1, \dots, T_n); \lambda) = \lambda^n e^{-\lambda T_n} \mathbf{1}_{\{0 < T_1 < \dots < T_n\}},
   $$
   with $\mathbf{1}_{\{0 < T_1 < \dots < T_n\}} = 1$ a.s.

2. Compute the MLE of $\lambda$.
   **Solution.**
   We deduce that the log-likelihood equals
   $$
   \ell((T_1, \dots, T_n); \lambda) = n \ln(\lambda) - \lambda T_n.
   $$
   Hence
   $$
   \frac{\partial \ell}{\partial \lambda}((T_1, \dots, T_n); \lambda) = \frac{n}{\lambda} - T_n = 0 \iff \lambda = \frac{n}{T_n}.
   $$
   Moreover, $n / T_n$ is a maximum since the log-likelihood function is concave. Indeed,
   $$
   \frac{\partial^2 \ell}{\partial \lambda^2}((T_1, \dots, T_n); \lambda) = -\frac{n}{\lambda^2} < 0.
   $$

<!-- page 50 -->
Note that we could recover this result using that knowing the point process is also equivalent to knowing the interarrival times which are i.i.d. with $\mathcal{E}(\lambda)$ distribution. Indeed,
$$
L((W_1, \dots, W_n); \lambda) = \prod_{i=1}^n \lambda e^{-\lambda W_i} \mathbf{1}_{\{W_i > 0\}} = \lambda^n e^{-\lambda \sum_{i=1}^n W_i} = \lambda^n e^{-\lambda T_n}.
$$

### III.2. Non-asymptotic properties of the MLE

> **Lemma 3.6.**
> Let $N$ be a homogeneous Poisson process with intensity $\lambda > 0$ and $n \ge 2$ a fixed integer. Then, the MLE of $\lambda$ satisfies
> $$
> \mathbb{E}[\hat{\lambda}_n] = \frac{n}{n - 1} \lambda.
> $$

In particular, from Theorem 3.9, we can deduce that in the case "fixed $n$", the Maximum Likelihood Estimator $\hat{\lambda}_n$ of $\lambda$ is biased, and asymptotically unbiased. However, we can consider
$$
\tilde{\lambda} = \frac{n - 1}{T_n},
$$
and prove that $\tilde{\lambda}$ is an unbiased estimator of $\lambda$, but not efficient.

*Proof. Worksheet 3, exercise 1.*

Moreover, to construct statistical tests with uniform asymptotic level, we need the following lemma.

> **Lemma 3.7.**
> Let $N$ be a homogeneous Poisson process with (unknown) intensity $\lambda > 0$ observed up to the $n$th first arrival time $T_n$ and consider the MLE $\hat{\lambda}_n = n/T_n$ of $\lambda$.
> Then, for all (fixed) $t > 0$,
> $$
> \lambda \mapsto \mathbb{P}_\lambda(\hat{\lambda}_n \ge t)
> $$
> is a non-decreasing function.

*Proof. In exercise. Consider $\lambda \le \mu$ and prove that $\mathbb{P}_\lambda(\hat{\lambda}_n \ge t) \le \mathbb{P}_\mu(\hat{\lambda}_n \ge t)$. Hint: Use Proposition 1.19.*
**Solution.**
Note that $\hat{\lambda}_n = n/T_n$ with $T_n \sim \Gamma(n, \lambda)$. Hence, by Proposition 1.19,
$$
2\lambda T_n \sim \Gamma(n, 1/2) \stackrel{(d)}{=} \chi^2(2n).
$$
Therefore, for $\lambda \le \mu$:
$$
\begin{aligned}
\mathbb{P}_\lambda(\hat{\lambda}_n \ge t) &= \mathbb{P}_\lambda\left( T_n \le \frac{n}{t} \right) = \mathbb{P}_\lambda\left( 2\lambda T_n \le \frac{2\lambda n}{t} \right) \\
&= \mathbb{P}\left( X \le \frac{2\lambda n}{t} \right) \quad \text{where } X \sim \chi^2(2n) \\
&\stackrel{(*)}{\le} \mathbb{P}\left( X \le \frac{2\mu n}{t} \right) \quad \text{since } \lambda \le \mu \\
&= \mathbb{P}_\mu\left( 2\mu T_n \le \frac{2\mu n}{t} \right) = \mathbb{P}_\mu(\hat{\lambda}_n \ge t),
\end{aligned}
$$
where $(*)$ holds from the fact that a c.d.f. is (always) a non-decreasing function.

<!-- page 51 -->
### III.3. Asymptotic distribution of the MLE

Before studying the asymptotic properties of the MLE, recall the univariate delta method.

> **Lemma 3.8 (Univariate delta method).**
> Let $(X_n)_{n \in \mathbb{N}}$ be a sequence of real-valued random variables, $Y$ be a real-valued random variable and $\theta$ be in $\mathbb{R}$ such that
> $$
> r_n (X_n - \theta) \xrightarrow[n \to +\infty]{\mathcal{L}} Y,
> $$
> where $(r_n)$ is a sequence of non-negative real numbers that tends toward infinity. Let $f : \mathbb{R} \to \mathbb{R}$ be a differentiable function at $\theta$. Then,
> $$
> f(X_n) \xrightarrow[n \to +\infty]{\mathbb{P}} f(\theta) \quad \text{and} \quad r_n (f(X_n) - f(\theta)) \xrightarrow[n \to +\infty]{\mathcal{L}} f'(\theta)Y.
> $$

This allows us to deduce the asymptotic behavior of the MLE.

> **Theorem 3.9.**
> Let $N$ be a homogeneous Poisson process with (unknown) intensity $\lambda > 0$ observed up to the $n$th first arrival time $T_n$ and consider the MLE $\hat{\lambda}_n = n/T_n$ of $\lambda$.
> 1. **LLN-type result:** $\hat{\lambda}_n$ is a consistent estimator of $\lambda$, that is
>    $$
>    \frac{n}{T_n} \xrightarrow[n \to +\infty]{\mathbb{P}} \lambda.
>    $$
> 2. **CLT-type result:** Moreover, it is asymptotically gaussian, that is
>    $$
>    \sqrt{n}\left( \frac{n}{T_n} - \lambda \right) \xrightarrow[n \to +\infty]{\mathcal{L}} \mathcal{N}(0, \lambda^2).
>    $$

> **Comment 3.10.** Note that here, the asymptotic variance equals $\lambda^2$.

*Proof. In exercise.*

1. Apply the CLT to the interarrival times $(W_i)_{i \ge 1}$, and prove that
   $$
   \sqrt{n}\left( \frac{1}{\hat{\lambda}_n} - \frac{1}{\lambda} \right) \xrightarrow[n \to +\infty]{\mathcal{L}} \mathcal{N}\left(0, \frac{1}{\lambda^2}\right). \tag{3.1}
   $$
   **Solution.**
   Recall that $\hat{\lambda}_n = n / T_n$ where
   $$
   \frac{T_n}{n} = \frac{1}{n} \sum_{i=1}^n W_i = \bar{W}_n.
   $$
   Applying the CLT to the sequence $(W_i)_{i \ge 1}$ which are i.i.d. exponential r.v. such that
   $$
   \mathbb{E}[|W_1|] = \mathbb{E}[W_1] = 1/\lambda < +\infty \quad \text{and} \quad \text{Var}(W_1) = 1/\lambda^2 < +\infty,
   $$
   leads to
   $$
   \sqrt{n} \frac{\frac{1}{\hat{\lambda}_n} - \frac{1}{\lambda}}{\sqrt{1/\lambda^2}} = \sqrt{n} \frac{\bar{W}_n - \frac{1}{\lambda}}{\sqrt{1/\lambda^2}} = \sqrt{n} \frac{\bar{W}_n - \mathbb{E}[W_1]}{\sqrt{\text{Var}(W_1)}} \xrightarrow[n \to +\infty]{\mathcal{L}} \mathcal{N}(0, 1).
   $$
   We deduce (3.1) by multiplying by $\sqrt{1/\lambda^2}$ (which multiplies the variance by $(1/\lambda^2)$).

<!-- page 52 -->
2. Deduce Theorem 3.9 from the delta method.
   **Solution.**
   Let us apply the delta method to the differentiable function
   $$
   f : x \in \mathbb{R}_+^* \mapsto \frac{1}{x}, \quad \text{with, by (3.1),} \quad
   \begin{cases}
   r_n = \sqrt{n}, \\
   X_n = 1/\hat{\lambda}_n, \\
   \theta = 1/\lambda, \\
   Y \sim \mathcal{N}(0, 1/\lambda^2).
   \end{cases}
   $$
   - **(a)** Since $f(X_n) = \hat{\lambda}_n$ and $f(\theta) = \lambda$, we first obtain the LLN-type result:
     $$
     \hat{\lambda}_n \xrightarrow[n \to +\infty]{\mathbb{P}} \lambda.
     $$
   - **(b)** Moreover, $\forall x \in \mathbb{R}_+^*$, $f'(x) = -1/x^2$. Hence,
     $$
     f'(\theta)Y = -\lambda^2 Y \sim \mathcal{N}\left(0, (-\lambda^2)^2 \times [1/\lambda^2]\right) = \mathcal{N}(0, \lambda^2),
     $$
     and we obtain the CLT-type result:
     $$
     \sqrt{n}(\hat{\lambda}_n - \lambda) \xrightarrow[n \to +\infty]{\mathcal{L}} \mathcal{N}(0, \lambda^2).
     $$

### III.4. Confidence intervals and statistical tests

Fix $n$ in $\mathbb{N}^*$ and let $\alpha$ in $(0, 1)$. We observe the $n$ first arrival times of a homogeneous Poisson process with intensity $\lambda > 0$ and we aim at constructing confidence intervals for $\lambda$ with confidence level $1 - \alpha$, or at testing, for instance,
$$
H_0 : \lambda \le \lambda_0 \quad \text{against} \quad H_1 : \lambda > \lambda_0, \quad \text{at level } \alpha.
$$

#### A) Non-asymptotic statistical inference

Based on Proposition 1.19, one may construct a (non-asymptotic) confidence interval for $\lambda$ with confidence level $1 - \alpha$, that is
$$
\left[ \frac{x_{2n, \alpha/2}}{2T_n} \; ; \; \frac{x_{2n, 1-\alpha/2}}{2T_n} \right],
$$
where $x_{d, \eta}$ denotes the $\eta$-quantile of a $\chi^2(d)$ distribution.

*Proof. In exercise.*
**Solution.**
By Proposition 1.19, we know that $T_n \sim \Gamma(n, \lambda)$, thus, $2\lambda T_n \sim \Gamma(n, 1/2) \stackrel{(d)}{=} \chi^2(2n)$. In particular:
$$
\mathbb{P}\left( \frac{x_{2n, \alpha/2}}{2T_n} \le \lambda \le \frac{x_{2n, 1-\alpha/2}}{2T_n} \right) = \mathbb{P}(x_{2n, \alpha/2} \le 2\lambda T_n \le x_{2n, 1-\alpha/2}) = 1 - \alpha.
$$

<!-- page 53 -->
Equivalently, one may prove that the test which rejects $H_0$ if $2\lambda_0 T_n \le x_{2n, \alpha}$ is of (non-asymptotic) level $\alpha$, that is
$$
\sup_{\lambda \le \lambda_0} \left\{ \mathbb{P}_\lambda(2\lambda_0 T_n \le x_{2n, \alpha}) \right\} \le \alpha.
$$
The test statistic of this test is
$$
X_n := 2\lambda_0 T_n.
$$
Moreover, if $\lambda = \lambda_0$, $X_n \sim \chi^2(2n)$. Hence, the p-value of this test equals
$$
\mathbb{P}_{\lambda_0}(X_n \le X_n^{\text{obs}}) = F_{2n}(X_n^{\text{obs}}),
$$
where $X_n^{\text{obs}}$ is the observed value of the test statistic (numerical application) and $F_{2n}$ denotes the c.d.f. of a $\chi^2(2n)$ distribution.

*Proof. Facultative homework.*

#### B) Asymptotic statistical inference

From Theorem 3.9, one can construct an asymptotic confidence interval that is
$$
\left[ \frac{n}{T_n \left(1 + \frac{z_{1-\alpha/2}}{\sqrt{n}}\right)} \; ; \; \frac{n}{T_n \left(1 - \frac{z_{1-\alpha/2}}{\sqrt{n}}\right)} \right],
$$
where $z_{1-\alpha/2}$ denotes the $(1 - \alpha/2)$-quantile of a standard gaussian distribution $\mathcal{N}(0, 1)$.

Note that this confidence interval requires that $1 - z_{1-\alpha/2}/\sqrt{n} > 0$, that is $n > z_{1-\alpha/2}^2$, and does not need Slutsky's Lemma.

*Proof. In exercise.*
**Solution.**
From Theorem 3.9, we know that
$$
\sqrt{n}\left( \frac{n}{T_n} - \lambda \right) \xrightarrow[n \to +\infty]{\mathcal{L}} \mathcal{N}(0, \lambda^2).
$$
Hence,
$$
\sqrt{n}\left( \frac{n}{\lambda T_n} - 1 \right) \xrightarrow[n \to +\infty]{\mathcal{L}} \mathcal{N}(0, 1).
$$
We deduce that, if $1 - \frac{z_{1-\alpha/2}}{\sqrt{n}} > 0$:
$$
\begin{aligned}
\mathbb{P}\left( \frac{n}{T_n \left(1 + \frac{z_{1-\alpha/2}}{\sqrt{n}}\right)} \le \lambda \le \frac{n}{T_n \left(1 - \frac{z_{1-\alpha/2}}{\sqrt{n}}\right)} \right) &= \mathbb{P}\left( -z_{1-\alpha/2} \le \sqrt{n}\left( \frac{n}{\lambda T_n} - 1 \right) \le z_{1-\alpha/2} \right) \\
&\xrightarrow[n \to +\infty]{} \mathbb{P}(-z_{1-\alpha/2} \le Z \le z_{1-\alpha/2}) = 1 - \alpha,
\end{aligned}
$$
where $Z$ denotes a standard gaussian r.v.

<!-- page 54 -->
Moreover, one can prove that the test which rejects $H_0$ if $\hat{\lambda}_n \ge \lambda_0(1 + z_{1-\alpha}/\sqrt{n})$ is (uniformly) asymptotically of level $\alpha$, that is
$$
\limsup_{n \to +\infty} \sup_{\lambda \le \lambda_0} \mathbb{P}_\lambda\left( \hat{\lambda}_n \ge \lambda_0\left(1 + \frac{z_{1-\alpha}}{\sqrt{n}}\right) \right) \le \alpha.
$$
This result relies on Lemma 3.7.

Note that, by isolating the quantile in the rejection rule, we recover the test statistic, that is
$$
Z_n := \sqrt{n}\left( \frac{\hat{\lambda}_n}{\lambda_0} - 1 \right) = \sqrt{n}\left( \frac{n}{\lambda_0 T_n} - 1 \right).
$$
Moreover, if $\lambda = \lambda_0$, $Z_n \xrightarrow[n \to +\infty]{\mathcal{L}} \mathcal{N}(0, 1)$. Hence the p-value of this test equals
$$
\mathbb{P}_{\lambda_0}(Z_n \ge Z_n^{\text{obs}}) \approx 1 - \Phi(Z_n^{\text{obs}}),
$$
where $Z_n^{\text{obs}}$ denotes the observed value of the test statistic (numerical application), and $\Phi$ denotes the c.d.f. of a standard gaussian distribution $\mathcal{N}(0, 1)$.

*Proof. Facultative homework.*

<!-- page 55 -->
### Take home message

| Estimator | Unbiased | Efficient | Asymptotically unbiased | Consistent | Asymptotically gaussian |
| :---: | :---: | :---: | :---: | :---: | :---: |
| $\hat{\lambda}_T = \frac{N_T}{T}$ | ✓ | ✓ | ✓ | ✓ | $\sqrt{T}(\hat{\lambda}_T - \lambda) \xrightarrow[T \to +\infty]{\mathcal{L}} \mathcal{N}(0, \lambda)$ |
| $\hat{\lambda}_n = \frac{n}{T_n}$ | ✗ | · | ✓ | ✓ | $\sqrt{n}(\hat{\lambda}_n - \lambda) \xrightarrow[n \to +\infty]{\mathcal{L}} \mathcal{N}(0, \lambda^2)$ |
| $\tilde{\lambda}_n = \frac{n-1}{T_n}$ | ✓ | ✗ | ✓ | ✓ | (?) |

<!-- page 56 -->

<!-- page 57 -->
# Chapter 4. Inhomogeneous Poisson processes

In this chapter, we consider inhomogeneous Poisson processes. They generalize homogeneous Poisson processes by allowing the rate of the process to be a varying function of time. In particular, by doing so, we drop the stationary increments property. Inhomogeneous Poisson processes provide suitable models in Reliability theory. For instance, if at each failure a system is repaired to its condition at the time of failure and placed in service again, then the failures are often modeled by an inhomogeneous Poisson process, provided the repair times can be neglected.

## I. Definitions and basic properties

> **Definition 4.1.**
> Let $\lambda : \mathbb{R}_+ \to \mathbb{R}_+$ be a real valued piecewise continuous (or at least locally integrable) function defined on $\mathbb{R}_+$. A point process $(N_t)_{t \in \mathbb{R}_+}$ is said to be an **inhomogeneous Poisson process** with intensity function $\lambda$ if:
> 1. $N_0 = 0$.
> 2. The process has independent increments.
> 3. For all $t > 0$, the number of points in $[0, t]$ has a Poisson distribution with parameter $\int_0^t \lambda(x)\,dx$, that is
>    $$
>    N_t \sim \mathcal{P}\left( \int_0^t \lambda(x)\,dx \right).
>    $$
> The function $\Lambda : t \mapsto \int_0^t \lambda(x)\,dx$ is called **cumulative intensity** of the Poisson process.

> **Example 4.2.**
> The power law process, sometimes called the Weibull process, is very popular in Reliability theory. It is an inhomogeneous Poisson process with intensity given by
> $$
> \lambda(x) = \frac{\beta}{\alpha} \left( \frac{x}{\alpha} \right)^{\beta - 1},
> $$
> where $\alpha$ and $\beta$ are two positive parameters. Compute its cumulative intensity.
> **Solution.**
> The cumulative intensity of a Weibull process equals
> $$
> \Lambda(t) = \int_0^t \frac{\beta x^{\beta - 1}}{\alpha^\beta}\,dx = \left( \frac{t}{\alpha} \right)^\beta.
> $$

<!-- page 58 -->
> **Proposition 4.3.**
> Let $N$ be an inhomogeneous Poisson process with intensity function $\lambda$.
> Denote for all Borel subset $A$ of $\mathbb{R}_+$, $N(A)$ the number of points in $A$. Note that if $(T_i)_{i \in \mathbb{N}^*}$ denote the arrival times, then
> $$
> N(A) = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{T_n \in A\}}.
> $$
> 1. For all bounded Borel set $A$ in $\mathbb{R}_+$, $N(A)$ is a.s. finite.
> 2. For all $s \ge 0$ and all $t > 0$, the random variable $N(]s, s + t]) = N_{s+t} - N_s$ has a Poisson distribution with parameter $\Lambda(s + t) - \Lambda(s) = \int_s^{s+t} \lambda(x)\,dx$, that is
>    $$
>    N_{s+t} - N_s \sim \mathcal{P}\left( \int_s^{s+t} \lambda(x)\,dx \right).
>    $$

*Proof.*

1. Let $A$ be a bounded Borel set in $\mathbb{R}_+$. Hence, $\exists t_A \in \mathbb{R}_+^*$ s.t. $A \subset [0, t_A]$. Therefore,
   $$
   N(A) = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{T_n \in A\}} \le \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{T_n \le t_A\}} = N_{t_A}.
   $$
   Moreover, $N_{t_A} \sim \mathcal{P}(\int_0^{t_A} \lambda(x)\,dx)$. Hence $\mathbb{E}[N_{t_A}] = \int_0^{t_A} \lambda(x)\,dx < +\infty$. We deduce that $N_{t_A}$ is finite a.s., and thus, so is $N(A)$.

2.
   - **(a)** Recall the moment generating function (or m.g.f.), or Laplace transform, of a Poisson r.v., and for all $t > 0$, express the m.g.f. of $N_t$ w.r.t. the cumulative intensity $\Lambda$ of the process.
     **Solution.**
     Let $X \sim \mathcal{P}(\mu)$. Then for all $u < 0$, $\mathbb{E}[e^{uX}] = \exp(\mu(e^u - 1))$.
     Hence, for any $t > 0$, since $N_t \sim \mathcal{P}(\Lambda(t))$, then for all $u < 0$:
     $$
     \mathbb{E}\left[e^{uN_t}\right] = \exp\left(\Lambda(t)[e^u - 1]\right).
     $$

   - **(b)** Let $s \ge 0$ and $t > 0$. Compute the m.g.f. of $N_{s+t} - N_s$, and conclude. *Hint: Use $N_{s+t} = N_{s+t} - N_s + N_s$.*
     **Solution.**
     Note that $N_{s+t} = N_{s+t} - N_s + N_s$, with $(N_{s+t} - N_s) \perp\!\!\!\perp N_s$ since the intervals $[0, s]$ and $(s, s + t]$ are disjoint. Hence
     $$
     \mathbb{E}\left[e^{uN_{s+t}}\right] = \mathbb{E}\left[e^{u(N_{s+t} - N_s)} e^{uN_s}\right] = \mathbb{E}\left[e^{u(N_{s+t} - N_s)}\right] \mathbb{E}\left[e^{uN_s}\right].
     $$
     Therefore, by the previous question:
     $$
     \exp\left( \Lambda(t + s)[e^u - 1] \right) = \mathbb{E}\left[e^{u(N_{s+t} - N_s)}\right] \exp\left( \Lambda(s)[e^u - 1] \right),
     $$
     i.e.
     $$
     \mathbb{E}\left[e^{u(N_{s+t} - N_s)}\right] = \exp\left( [\Lambda(s + t) - \Lambda(s)][e^u - 1] \right) = \exp\left( \int_s^{s+t} \lambda(x)\,dx\,[e^u - 1] \right).
     $$
     We recognize the m.g.f. of a $\mathcal{P}\left(\int_s^{s+t} \lambda(x)\,dx\right)$, which ends the proof.

<!-- page 59 -->
Let us notice several points.

- **a)** If for all $t \ge 0$, $\lambda(t) = \lambda$, that is $\lambda$ is a constant function, $\int_s^{s+t} \lambda\,dx = \lambda t$, and one recovers the definition of a homogeneous Poisson process with rate $\lambda$.
- **b)** Notice that Definition 4.1 characterizes the distribution of the process $N$. Indeed, for any $n$ and for any $0 \le t_1 < t_2 < \dots < t_n$, Definition 4.1 and Proposition 4.3 both provide the distribution of $(N_{t_1}, N_{t_2} - N_{t_1}, \dots, N_{t_n} - N_{t_{n-1}})$. Hence, the distribution of the process $(N_t)_{t \in \mathbb{R}_+}$ is well known.
- **c)** The terminology of point process theory has been criticized for being too varied. In addition to the word *point* often being omitted, the homogeneous Poisson (point) process is also called a stationary Poisson (point) process, as well as uniform Poisson (point) process. The inhomogeneous Poisson (point) process, as well as being called nonhomogeneous, is also referred to as the non-stationary Poisson process.
- **d)** The positive measure $\mu$ defined on any Borel set $A$ by $\mu(A) = \int_A \lambda(x)\,dx$ (or equivalently $d\mu(x) = \lambda(x)\,dx$) is called the **intensity measure** of the Poisson process $N$. This definition of an inhomogeneous Poisson process can be generalized to any non-negative intensity measure on $\mathbb{R}_+$ that is locally bounded (i.e. for all bounded Borel subset $A$ in $\mathbb{R}_+$, $\mu(A) < +\infty$), not necessarily absolutely continuous with respect to the Lebesgue measure by replacing assumption (3) by: for all Borel set $A$ in $\mathbb{R}_+$, $N(A) \sim \mathcal{P}(\mu(A))$.
  However, in this lecture, we restrict attention to the case $d\mu(x) = \lambda(x)\,dx$, that is Poisson processes with intensity functions.

> **Property 4.4.**
> Let $N$ be an inhomogeneous Poisson process with intensity function $\lambda$ on $\mathbb{R}_+$.
> 1. The mean function of $N$ equals
>    $$
>    \forall t \ge 0, \quad m(t) = \mathbb{E}[N_t] = \int_0^t \lambda(x)\,dx = \Lambda(t).
>    $$
>    In particular, the rate (see Definition 2.7) of the Poisson process $N$ equals:
>    $$
>    \forall t \ge 0, \quad w(t) = m'(t) = \lambda(t).
>    $$
> 2. Moreover, for all $t \ge 0$, as $h \to 0$:
>    - **(a)** $\mathbb{P}(N_{t+h} - N_t = 0) = 1 - \lambda(t)h + o(h)$.
>    - **(b)** $\mathbb{P}(N_{t+h} - N_t = 1) = \lambda(t)h + o(h)$.
>    - **(c)** $\mathbb{P}(N_{t+h} - N_t \ge 2) = o(h)$.

*Proof. Facultative exercise.*

1. Immediate.
2. *Hint: Since $\Lambda'(t) = \lambda(t)$, one has $\Lambda(t + h) - \Lambda(t) = \lambda(t)h + o(h)$.*
   **Solution.**
   By Property 4.4, $N_{t+h} - N_t \sim \mathcal{P}(\Lambda(t + h) - \Lambda(t))$. Then:
   - **(a)** $\mathbb{P}(N_{t+h} - N_t = 0) = e^{-[\Lambda(t+h) - \Lambda(t)]} = e^{-\lambda(t)h + o(h)} = 1 - \lambda(t)h + o(h)$.

<!-- page 60 -->
   - **(b)**
     $$
     \begin{aligned}
     \mathbb{P}(N_{t+h} - N_t = 1) &= [\Lambda(t + h) - \Lambda(t)] e^{-[\Lambda(t+h) - \Lambda(t)]} \\
     &= [\lambda(t)h + o(h)] \times [1 - \lambda(t)h + o(h)] \\
     &= \lambda(t)h + o(h).
     \end{aligned}
     $$
   - **(c)** $\mathbb{P}(N_{t+h} - N_t \ge 2) = 1 - [\mathbb{P}(N_{t+h} - N_t = 0) + \mathbb{P}(N_{t+h} - N_t = 1)] = o(h)$.

As in the homogeneous case, one can provide an equivalent definition for inhomogeneous Poisson processes.

> **Definition 4.5.**
> The point process $N$ is said to be an **inhomogeneous Poisson process** with intensity function $\lambda$ defined on $\mathbb{R}_+$ if:
> 1. $N_0 = 0$.
> 2. The process has independent increments.
> 3. For all positive $t$, one has as $h \to 0$:
>    - **(a)** $\mathbb{P}(N_{t+h} - N_t = 1) = \lambda(t)h + o(h)$.
>    - **(b)** $\mathbb{P}(N_{t+h} - N_t \ge 2) = o(h)$.

The proof is the same as in the homogeneous case, except that we obtain a linear differential equation with non-constant coefficients.

Notice that, a Poisson process $N$ with intensity function $\lambda$ is regular since, for all $t > 0$,
$$
\lim_{h \to 0} \frac{\mathbb{P}(N_{t+h} - N_t \ge 2)}{h} = 0.
$$
In particular, it is a simple process, that is no more than 1 event can occur at the same time, and each jump of the counting process equals 1 (almost surely). Notice that it is not the case anymore when we generalize to inhomogeneous Poisson processes with non-diffuse intensity measure (i.e. that contains atoms).

## II. Construction of an inhomogeneous Poisson process

### II.1. Time change

A first way of constructing an inhomogeneous Poisson process with intensity function $\lambda$ on $\mathbb{R}_+$ (and thus proving that such processes exist) is by changing the time scale of a homogeneous Poisson process. Before, let us recall the definition and some properties of the generalized inverse function.

<!-- page 61 -->
#### A) Generalized inverse function

> **Definition 4.6.**
> Let $F : \mathbb{R} \to \mathbb{R}_+$ be a non-decreasing càdlàg non-negative function defined on $\mathbb{R}$. Its generalized inverse function is defined by
> $$
> F^- : \begin{array}{ccc}
> \mathbb{R}_+ & \longrightarrow & \mathbb{R} \cup \{\pm\infty\} \\
> t & \longmapsto & \inf\{x \in \mathbb{R} \; ; \; F(x) \ge t\}
> \end{array},
> $$
> with the conventions $\inf(\emptyset) = +\infty$ and $\inf(\mathbb{R}) = -\infty$.

Note that if $F$ is a c.d.f., then $F^-$ is the corresponding quantile function.

```
       1 |
         |          ------------------
      t3 |........./
      t2 |........•
         |       /
      t1 |....../
       0 +-----+---+--------------------> x
               |   |
             F-(t1) F-(t2)=F-(t3)
```
*(Generalized inverse function $F^-$ illustration)*

> **Property 4.7.**
> The generalized inverse function satisfies the following properties:
> - **i.** $F^-$ is a non-decreasing function.
> - **ii.** $F(x) \ge t \iff x \ge F^-(t)$.
> - **iii.** $F^-$ is left-continuous, that is for all $t_0$ in $\mathbb{R}$, $\lim_{\substack{t \to t_0 \\ t < t_0}} F^-(t) = F^-(t_0)$.
> - **iv.** For all $t \in \mathbb{R}_+$, $F(F^-(t)) \ge t$.
>   Moreover, if $F$ is continuous, for all $t \in \mathbb{R}_+$, $F(F^-(t)) = t$.
> - **v.** For all $x \in \mathbb{R}$, $F^-(F(x)) \le x$.
>   Moreover, if $F$ is increasing, for all $x \in \mathbb{R}$, $F^-(F(x)) = x$.
> - **vi.** In particular, if $F$ is increasing and continuous, then $F$ is bijective with inverse function $F^{-1} = F^-$.

*Proof.* See Appendix A.II.

<!-- page 62 -->
#### B) Back to inhomogeneous Poisson processes

> **Proposition 4.8.**
> Let $\lambda : \mathbb{R}_+ \to \mathbb{R}_+$ be a locally integrable function, and denote
> $$
> \Lambda : t \in \mathbb{R}_+ \mapsto \int_0^t \lambda(x)\,dx.
> $$
> Note that $\Lambda$ is a non-decreasing continuous function.
> 1. Let $\tilde{N}$ be a homogeneous Poisson process with (constant) rate 1, and denote
>    $$
>    0 < S_1 < S_2 < \dots < S_n < \dots
>    $$
>    its arrival times. Define respectively for all $n \in \mathbb{N}^*$, and for all $t \in \mathbb{R}_+$:
>    $$
>    T_n = \Lambda^-(S_n) \quad \text{and} \quad N_t = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{T_n \le t\}},
>    $$
>    where $\Lambda^-$ denotes the generalized inverse function of $\Lambda$. Then $N = (N_t)_{t \in \mathbb{R}_+}$ is an inhomogeneous Poisson process with intensity function $\lambda$.
> 2. Conversely, assume that $\Lambda$ is increasing (and thus a bijective function, since it is always continuous). Let $N$ be an inhomogeneous Poisson process with intensity function $\lambda$ on $\mathbb{R}_+$. Denote
>    $$
>    0 < T_1 < T_2 < \dots < T_n < \dots,
>    $$
>    its arrival times. Define respectively for all $n$ in $\mathbb{N}^*$, and for all $t$ in $\mathbb{R}_+$:
>    $$
>    S_n = \Lambda(T_n) \quad \text{and} \quad \tilde{N}_t = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{S_n \le t\}}.
>    $$
>    Then the point process $\tilde{N} = (\tilde{N}_t)_{t \in \mathbb{R}_+}$ is a homogeneous Poisson process with (constant) rate 1.

*Proof. In exercise.*

1.
   - **(a)** Let $t > 0$. Express $N_t$ in terms of $(S_n)_{n \in \mathbb{N}^*}$, and deduce, by Property 4.7, that $N_t = \tilde{N}_{\Lambda(t)}$.
     **Solution.**
     $$
     N_t = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{T_n \le t\}} = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{\Lambda^-(S_n) \le t\}}.
     $$
     Yet, by point ii. of Property 4.7:
     $$
     \Lambda^-(S_n) \le t \iff S_n \le \Lambda(t).
     $$
     Hence
     $$
     N_t = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{S_n \le \Lambda(t)\}} = \tilde{N}_{\Lambda(t)}.
     $$

<!-- page 63 -->
   - **(b)** Conclude.
     **Solution.**
     - **i.** $N_0 = \tilde{N}_0 = 0$ since $\tilde{N}$ is a homogeneous Poisson process and $\Lambda(0) = \int_0^0 \lambda(x)\,dx = 0$.
     - **ii.** Let $t_1 < t_2 \le s_1 < s_2$. Then, since $\Lambda$ is non-decreasing:
       $$
       \Lambda(t_1) \le \Lambda(t_2) \le \Lambda(s_1) \le \Lambda(s_2).
       $$
       In particular, since $\tilde{N}$ is a homogeneous Poisson process:
       $$
       N_{t_2} - N_{t_1} = \tilde{N}_{\Lambda(t_2)} - \tilde{N}_{\Lambda(t_1)} \perp\!\!\!\perp \tilde{N}_{\Lambda(s_2)} - \tilde{N}_{\Lambda(s_1)} = N_{s_2} - N_{s_1}.
       $$
       Note that this remains true if $\Lambda(t_2) = \Lambda(t_1)$ (or if $\Lambda(s_2) = \Lambda(s_1)$) since in this case, the corresponding increment is constant (equal to zero).
     - **iii.** Finally, since $\tilde{N}$ is a homogeneous Poisson process with rate 1:
       $$
       N_t = \tilde{N}_{\Lambda(t)} \sim \mathcal{P}(1 \times \Lambda(t)).
       $$

2. Assume that $\Lambda$ is increasing (and thus bijective), and denote $\Lambda^{-1}$ its inverse function.
   - **(a)** What is the monotony of $\Lambda^{-1}$? Deduce that for all $u, t$:
     $$
     \Lambda(u) \le t \iff u \le \Lambda^{-1}(t).
     $$
     **Solution.**
     Since $\Lambda^{-1}$ is the inverse function of an increasing function, it is also increasing. We deduce that:
     $\Rightarrow$ if $\Lambda(u) \le t$, then $u = \Lambda^{-1}(\Lambda(u)) \le \Lambda^{-1}(t)$,
     $\Leftarrow$ if $u \le \Lambda^{-1}(t)$, then $\Lambda(u) \le \Lambda(\Lambda^{-1}(t)) = t$.

   - **(b)** Deduce that for all $t > 0$, $\tilde{N}_t = N_{\Lambda^{-1}(t)}$.
     **Solution.**
     Let $t > 0$. Then,
     $$
     \tilde{N}_t = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{S_n \le t\}} = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{\Lambda(T_n) \le t\}} = \sum_{n \in \mathbb{N}^*} \mathbf{1}_{\{T_n \le \Lambda^{-1}(t)\}} = N_{\Lambda^{-1}(t)}.
     $$

   - **(c)** Conclude.
     **Solution.**
     - **i.** Since $\Lambda(0) = 0$, then $\Lambda^{-1}(0) = 0$ and $\tilde{N}_0 = N_0 = 0$.
     - **ii.** $\tilde{N}$ has independent increments since $N$ does.

<!-- page 64 -->
     - **iii.** By Proposition 4.3:
       $$
       \tilde{N}_{s+t} - \tilde{N}_s = N_{\Lambda^{-1}(t+s)} - N_{\Lambda^{-1}(s)} \sim \mathcal{P}\left( \int_{\Lambda^{-1}(s)}^{\Lambda^{-1}(s+t)} \lambda(x)\,dx \right),
       $$
       with
       $$
       \int_{\Lambda^{-1}(s)}^{\Lambda^{-1}(s+t)} \lambda(x)\,dx = \Lambda(\Lambda^{-1}(s + t)) - \Lambda(\Lambda^{-1}(s)) = s + t - s = t.
       $$
       Hence,
       $$
       \tilde{N}_{s+t} - \tilde{N}_s \sim \mathcal{P}(1 \times t).
       $$
     We deduce from all three points above that $\tilde{N}$ is a homogeneous Poisson process with (constant) rate 1.

### II.2. Acceptance/Rejection (or Thinning)

The following result leads to one of the most popular methods for generating an inhomogeneous Poisson process, sometimes called the thinning method. It is the "point process" analogue of the acceptance/rejection method classically used to simulate r.v.
It generalizes Proposition 2.27 which claims that if you randomly classify each point of a homogeneous Poisson process with rate $\lambda$ into two categories (with probability $p$ and $1 - p$), the two corresponding point processes are also homogeneous Poisson processes with corresponding rates $p\lambda$ and $(1 - p)\lambda$. In the inhomogeneous case, we allow $p$ to be a varying function of time.

> **Theorem 4.9.**
> Let $\lambda : \mathbb{R}_+ \to \mathbb{R}_+^*$ be a locally integrable non-negative function. Assume that $\lambda$ is upper-bounded by $M \in \mathbb{R}_+^*$ on a subset $[0, t)$ of $\mathbb{R}_+$ (with $t$ in $\mathbb{R}_+^* \cup \{+\infty\}$), that is
> $$
> \forall x \in [0, t), \quad \lambda(x) \le M.
> $$
> Let $N$ be a homogeneous Poisson process with rate $M$ and classify each point $u$ of $N$ as type I with probability $p(u) := \lambda(u)/M$ and type II otherwise. Then, the resulting point processes $N^I$ and $N^{II}$ of type I and type II events satisfy:
> - **a)** $N^I$ is an inhomogeneous Poisson process with intensity function $\lambda$.
>   $N^{II}$ is an inhomogeneous Poisson process with intensity function $M - \lambda$.
> - **b)** $N^I \perp\!\!\!\perp N^{II}$.

*Idea of proof.*

- **a)** The independence of the increments comes from the independence of the increments of the homogeneous Poisson process $N$, and the independence between the tags.
- **b)** Let us now prove the Poisson distribution, and that $N^I \perp\!\!\!\perp N^{II}$.
  - **Modeling:** We know (c.f. Proposition 2.23) that conditionally on $\{N_t = n\}$, $T_1, \dots, T_n$ behaves as the order statistic $U_{(1)} \le \dots \le U_{(n)}$ associated to $n$ i.i.d. r.v. $U_1, \dots, U_n$ with distribution $\mathcal{U}([0, t])$. The classifying process can be seen as follows: Let $Y_1, \dots, Y_n$ be $n$ i.i.d. r.v. with distribution $\mathcal{U}([0, M])$, independent of $(U_1, \dots, U_n)$.
    $$
    \text{Classify } U_i \text{ as } \begin{cases}
    \text{type I} & \text{if } Y_i \le \lambda(U_i), \\
    \text{type II} & \text{otherwise.}
    \end{cases}
    $$

<!-- page 65 -->
    Then,
    $$
    \begin{aligned}
    \mathbb{P}(U_i \text{ of type I} \mid U_i = u) &= \mathbb{P}(Y_i \le \lambda(U_i) \mid U_i = u) \\
    &= \mathbb{P}(Y_i \le \lambda(u) \mid U_i = u) \\
    &= \mathbb{P}(Y_i \le \lambda(u)), \quad \text{as } Y_i \perp\!\!\!\perp U_i \\
    &= \frac{\lambda(u)}{M} = p(u), \quad \text{as } Y_i \sim \mathcal{U}([0, M]).
    \end{aligned}
    $$

  - **Towards a binomial distribution:** Conditionally on $\{N_t = n\}$, we thus obtain
    $$
    N_t^I = \sum_{i=1}^n \mathbf{1}_{\{Y_i \le \lambda(U_i)\}},
    $$
    where the $\mathbf{1}_{\{Y_i \le \lambda(U_i)\}}$ are independent (by independence of the couples $(U_i, Y_i)_i$) Bernoulli r.v. with parameter
    $$
    \tilde{p} := \mathbb{P}(Y_i \le \lambda(U_i)) = \int_0^{+\infty} \underbrace{\mathbb{P}(Y_i \le \lambda(U_i) \mid U_i = u)}_{p(u)} \underbrace{f_{U_i}(u)}_{\frac{1}{t}\mathbf{1}_{\{u \in [0, t]\}}}\,du = \frac{1}{t} \int_0^t p(u)\,du.
    $$
    We deduce that the conditional distribution of $N_t^I$ given $\{N_t = n\}$, is a binomial distribution with parameters $n$ and $\tilde{p}$, that is
    $$
    \mathcal{L}(N_t^I \mid N_t = n) = \mathcal{B}(n, \tilde{p}).
    $$

  - **Unconditioning:** Let $k, l \in \mathbb{N}$. Then
    $$
    \begin{aligned}
    \mathbb{P}(N_t^I = k, N_t^{II} = l) &= \mathbb{P}(N_t^I = k, N_t = k + l) \\
    &= \mathbb{P}(N_t^I = k \mid N_t = k + l) \mathbb{P}(N_t = k + l) \\
    &= \binom{k + l}{k} \tilde{p}^k (1 - \tilde{p})^l \times \frac{(Mt)^{k+l}}{(k + l)!} e^{-Mt(\tilde{p} + 1 - \tilde{p})} \\
    &= \frac{(Mt\tilde{p})^k}{k!} e^{-Mt\tilde{p}} \times \frac{(Mt[1 - \tilde{p}])^l}{l!} e^{-Mt(1 - \tilde{p})}.
    \end{aligned}
    $$
    We deduce that
    $$
    \mathbb{P}(N_t^I = k) = \sum_{l=0}^{+\infty} \mathbb{P}(N_t^I = k, N_t^{II} = l) = \frac{(Mt\tilde{p})^k}{k!} e^{-Mt\tilde{p}},
    $$
    thus
    $$
    N_t^I \sim \mathcal{P}(Mt\tilde{p}), \quad \text{with} \quad Mt\tilde{p} = Mt \frac{1}{t} \int_0^t \frac{\lambda(u)}{M}\,du = \int_0^t \lambda(u)\,du.
    $$
    In the same way, we prove that $N_t^{II} \sim \mathcal{P}(Mt[1 - \tilde{p}])$, with
    $$
    Mt[1 - \tilde{p}] = Mt - \int_0^t \lambda(u)\,du = \int_0^t [M - \lambda(u)]\,du.
    $$
    Finally, we also deduce that $N_t^I \perp\!\!\!\perp N_t^{II}$ because of that product form above.

<!-- page 66 -->
```
Algorithm 1 Thinning algorithm
Input: intensity function λ; simulation interval [0, T]; upper-bound M > 0.

1. Generate a homogeneous Poisson process with constant rate M on [0, T];
   Denote t1, ..., tn the observed arrival times on [0, T].
2. Generate independently n i.i.d r.v. with distribution U[0, M];
   Denote y1, ..., yn the observed values.
3. for i in 1 to n, keep ti if yi ≤ λ(ti).
4. Order the kept points.

Output: the set of ordered kept points are the arrival times of an inhomogeneous Poisson process
with intensity function λ on [0, T].
```

```
   M +
     |   u1 -  ◦
     |   u7 -                                       ◦
     |   u4 -                     ◦
     |   u3 -            ◦
     |   u5 -                         ×
     |   u2 -        ×                          ◦      λ(t)
     |   u6 -  ◦     ×   ◦        ◦   ×         ◦     ◦
   0 +---------+-----+---+--------+---+---------+-----+-----> t
               t1    t2  t3       t4  t5        t6    t7

Observed arrival times:
• N   : t1, ..., t7;
• N_I : t2, t5;          (kept: yi ≤ λ(ti))
• N_II: t1, t3, t4, t6, t7.  (rejected: yi > λ(ti))
```
*(Thinning / Acceptance-Rejection method illustration)*

## III. Arrival times

### III.1. First arrival time after a given instant

> **Property 4.10.**
> Let $N$ be a Poisson process with intensity function $\lambda : \mathbb{R}_+ \to \mathbb{R}_+$.
> - **i.** Then the survival function of the first arrival time $T_1$ equals for all $s > 0$:
>   $$
>   \mathbb{P}(T_1 > s) = \mathbb{P}(N_s = 0) = e^{-\int_0^s \lambda(x)\,dx} = e^{-\Lambda(s)}.
>   $$
>   It is not exponentially distributed anymore (except if $\Lambda(s) = \lambda s$, i.e. if the Poisson process is homogeneous). Note that $\lambda$ is the hazard rate of $T_1$.
> - **ii.** Denote $X_t$ the waiting time at time $t$, that is
>   $$
>   X_t = T_{N_t + 1} - t.
>   $$

```
                                          Xt
                                |---------|-----------------|
--------------------------------|---------|-----------------|---------> t
                               TNt        t               TNt+1
```

> Then, for all $s > 0$:
> $$
> \mathbb{P}(X_t > s) = \mathbb{P}(N_{t+s} - N_t = 0) = e^{-\int_t^{s+t} \lambda(x)\,dx} = e^{-\int_0^s \lambda(x+t)\,dx}.
> $$
> In particular, from Proposition 1.15, the hazard rate of $X_t$ is $h_t : x \mapsto \lambda(x + t)$.

<!-- page 67 -->
### III.2. (Conditional) distribution of the arrival times

> **Proposition 4.11.**
> Let $N$ be an inhomogeneous Poisson process with intensity function $\lambda : \mathbb{R}_+ \to \mathbb{R}_+$. Fix $n \in \mathbb{N}^*$ and denote $(T_1, \dots, T_n)$ the $n$ first arrival times. Then $(T_1, \dots, T_n)$ has the following joint density w.r.t. the Lebesgue measure on $\mathbb{R}^n$:
> $$
> (t_1, \dots, t_n) \mapsto e^{-\int_0^{t_n} \lambda(x)\,dx} \left( \prod_{i=1}^n \lambda(t_i) \right) \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}.
> $$

Note that if $\lambda$ is a constant function, we recover the density in the homogeneous case (Lemma 2.19).

*Proof. In exercise (assuming the cumulative intensity $\Lambda$ is increasing).*

1. Assume that $\Lambda$ is increasing and let $g : \mathbb{R}^n \to \mathbb{R}$ be a measurable bounded function. Use the Time change Proposition 4.8 and prove that
   $$
   \mathbb{E}[g(T_1, \dots, T_n)] = \int_{\mathbb{R}^n} g(\Lambda^{-1}(s_1), \dots, \Lambda^{-1}(s_n)) e^{-s_n} \mathbf{1}_{\{0 < s_1 < \dots < s_n\}}\,ds_1 \dots ds_n.
   $$
   **Solution.**
   Define for all $1 \le i \le n$, $S_i = \Lambda(T_i)$. Then, by Proposition 4.8, $S_1, \dots, S_n$ are the $n$ first arrival times of a homogeneous Poisson process with parameter 1. Hence, by Lemma 2.19, the joint density of $(S_1, \dots, S_n)$ equals
   $$
   (s_1, \dots, s_n) \mapsto 1^n e^{-1 \times s_n} \mathbf{1}_{\{0 < s_1 < \dots < s_n\}}.
   $$
   In particular,
   $$
   \begin{aligned}
   \mathbb{E}[g(T_1, \dots, T_n)] &= \mathbb{E}\left[ g(\Lambda^{-1}(S_1), \dots, \Lambda^{-1}(S_n)) \right] \\
   &= \int_{\mathbb{R}^n} g(\Lambda^{-1}(s_1), \dots, \Lambda^{-1}(s_n)) e^{-s_n} \mathbf{1}_{\{0 < s_1 < \dots < s_n\}}\,ds_1 \dots ds_n.
   \end{aligned}
   $$

2. By a well-chosen change of variables, show that
   $$
   \mathbb{E}[g(T_1, \dots, T_n)] = \int_{\mathbb{R}^n} g(t_1, \dots, t_n) e^{-\Lambda(t_n)} \lambda(t_1) \dots \lambda(t_n) \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}\,dt_1 \dots dt_n.
   $$
   **Solution.**
   Consider the diagonal change of variables $s_i = \Lambda(t_i)$ for all $1 \le i \le n$, such that
   $$
   ds_1 \dots ds_n = \lambda(t_1)\,dt_1 \dots \lambda(t_n)\,dt_n,
   $$
   since $\Lambda'(t) = \lambda(t)$ for all $t$. As $\Lambda$ is increasing, $\mathbf{1}_{\{0 < s_1 < \dots < s_n\}} = \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}$, we deduce from the previous question that
   $$
   \mathbb{E}[g(T_1, \dots, T_n)] = \int_{\mathbb{R}^n} g(t_1, \dots, t_n) e^{-\Lambda(t_n)} \lambda(t_1) \dots \lambda(t_n) \mathbf{1}_{\{0 < t_1 < \dots < t_n\}}\,dt_1 \dots dt_n.
   $$

<!-- page 68 -->
3. Conclude.
   **Solution.**
   Since this is true for all bounded measurable functions $g$, we deduce that the joint density of $(T_1, \dots, T_n)$ is
   $$
   (t_1, \dots, t_n) \mapsto e^{-\Lambda(t_n)} \lambda(t_1) \dots \lambda(t_n) \mathbf{1}_{\{0 < t_1 < \dots < t_n\}},
   $$
   which ends the proof.

We can thus deduce the conditional distribution of the arrival times.

> **Proposition 4.12.**
> Let $n \in \mathbb{N}^*$ and $t \in \mathbb{R}_+^*$. Given that $\{N_t = n\}$, the $n$ arrival times $0 < T_1 < \dots < T_n$ in $[0, t]$ of an inhomogeneous Poisson process with intensity function $\lambda : \mathbb{R}_+ \to \mathbb{R}_+$, has conditional density w.r.t. the Lebesgue measure on $\mathbb{R}^n$:
> $$
> (t_1, \dots, t_n) \mapsto \frac{n!}{\left[ \int_0^t \lambda(x)\,dx \right]^n} \left( \prod_{i=1}^n \lambda(t_i) \right) \mathbf{1}_{\{0 < t_1 < \dots < t_n \le t\}}.
> $$
> In particular, conditionally on $\{N_t = n\}$, $(T_1, \dots, T_n)$ behaves as the order statistic associated to $n$ i.i.d. r.v. with common density
> $$
> s \in \mathbb{R} \mapsto \frac{\lambda(s)}{\int_0^t \lambda(x)\,dx} \mathbf{1}_{\{0 < s \le t\}}.
> $$

*Proof.* Similar to the one of Proposition 2.23.

This result leads to a new construction of an inhomogeneous Poisson process.

> **Proposition 4.13.**
> Let $\lambda : \mathbb{R}_+ \to \mathbb{R}_+$ be an intensity function. Let $Y \sim \mathcal{P}\left( \int_0^T \lambda(x)\,dx \right)$.
> - If $Y = 0$, then for all $0 \le t \le T$, $N_t = 0$ (the point process has no points in $[0, T]$).
> - Otherwise, let $S_1, \dots, S_Y$ be i.i.d. r.v. with density w.r.t. the Lebesgue measure on $\mathbb{R}$:
>   $$
>   s \in \mathbb{R} \mapsto \frac{\lambda(s)}{\int_0^T \lambda(x)\,dx} \mathbf{1}_{\{0 < s \le T\}}.
>   $$
>   Sort the values $(T_1, \dots, T_Y) = (S_{(1)}, \dots, S_{(Y)})$ and consider the associated counting process defined for all $0 \le t \le T$:
>   $$
>   N_t = \sum_{n=1}^Y \mathbf{1}_{\{T_n \le t\}} = \sum_{i=1}^Y \mathbf{1}_{\{S_i \le t\}}.
>   $$
> Then $(N_t)_{t \in [0, T]}$ is an inhomogeneous Poisson process with intensity function $\lambda$ on $[0, T]$.

The proof of Proposition 4.13 is detailed in Appendix A.III.2.

<!-- page 69 -->
## IV. Parametric statistics for inhomogeneous Poisson processes

Let us observe an inhomogeneous Poisson process $N$ with intensity function $\lambda_\theta$ which depends on an unknown parameter $\theta \in \Theta$. Then, one can introduce the Maximum Likelihood Estimator (MLE) of $\theta$. To do so, we need to compute the likelihood $L(N; \theta)$, and, if it exists, to compute the MLE of $\theta$:
$$
\hat{\theta}_{MLE} \in \operatorname{argmax}_{\theta \in \Theta} L(N; \theta),
$$
and estimate the intensity function by $\lambda_{\hat{\theta}_{MLE}}$.

As in the homogeneous case, there are two types of observation: either the observation window is fixed (and you observe a random number of events), or $n$ is fixed (and the observation window is random).

> **Proposition 4.14 (Likelihood in the "fixed window case").**
> Let $N$ be an inhomogeneous Poisson process with intensity function $\lambda_\theta : \mathbb{R}_+ \to \mathbb{R}_+$ depending on an unknown parameter $\theta$ in $\Theta$. Consider the number $N_T$ of "events" that have occurred in $[0, T]$, and denote $(T_1, \dots, T_{N_T})$ the corresponding arrival times. Then likelihood of $\theta$ is equal to
> $$
> L((N_t)_{t \in [0, T]}; \theta) = \left( \prod_{i=1}^{N_T} \lambda_\theta(T_i) \right) e^{-\int_0^T \lambda_\theta(x)\,dx}.
> $$

> **Proposition 4.15 (Likelihood in the "fixed number of points case").**
> Let $N$ be an inhomogeneous Poisson process with intensity function $\lambda_\theta : \mathbb{R}_+ \to \mathbb{R}_+$ depending on an unknown parameter $\theta$ in $\Theta$. Fix $n \in \mathbb{N}^*$ and denote $(T_1, \dots, T_n)$ the $n$ first arrival times. Then likelihood of $\theta$ is equal to
> $$
> L((T_i)_{1 \le i \le n}; \theta) = \left( \prod_{i=1}^n \lambda_\theta(T_i) \right) e^{-\int_0^{T_n} \lambda_\theta(x)\,dx}.
> $$

You will have the opportunity to deepen your understanding of statistical inference for inhomogeneous Poisson processes through various mini-projects.

<!-- page 70 -->

<!-- page 71 -->
# Appendix A. Mathematical tools and complementary proofs

## I. Integration by substitution

> **Theorem A.1.**
> Let $U$ and $V$ be two open subsets of $\mathbb{R}^n$, and $\phi : U \to V$ be a $C^1$-diffeomorphism, that is $\phi$ is bijective, of class $C^1$, and with differential invertible at each point of $U$ (i.e. $\det(J_\phi(u)) \ne 0$). Then for any non-negative measurable function $f : V \to \mathbb{R}$:
> $$
> \int_V f(x)\,dx = \int_U f(\phi(u)) |\det(J_\phi(u))|\,du,
> $$
> where $J_\phi(u)$ denotes the Jacobian matrix of $\phi$ at point $u$:
> $$
> J_\phi(u) = \begin{pmatrix}
> \frac{\partial \phi_1}{\partial u_1} & \dots & \frac{\partial \phi_1}{\partial u_n} \\
> \vdots & \ddots & \vdots \\
> \frac{\partial \phi_n}{\partial u_1} & \dots & \frac{\partial \phi_n}{\partial u_n}
> \end{pmatrix}.
> $$

## II. Generalized inverse function

> **Definition A.2.**
> Let $F : \mathbb{R} \to \mathbb{R}_+$ be a non-decreasing càdlàg non-negative function defined on $\mathbb{R}$. Its generalized inverse function is defined by
> $$
> F^- : \begin{array}{ccc}
> \mathbb{R}_+ & \longrightarrow & \mathbb{R} \cup \{\pm\infty\} \\
> t & \longmapsto & \inf\{x \in \mathbb{R} \; ; \; F(x) \ge t\}
> \end{array},
> $$
> with the conventions $\inf(\emptyset) = +\infty$ and $\inf(\mathbb{R}) = -\infty$.
> Note that if $F$ is a c.d.f., then $F^-$ is the corresponding quantile function.

> **Property A.3.**
> The generalized inverse function satisfies the following properties.
> - **i.** $F^-$ is a non-decreasing function.
> - **ii.** $F(x) \ge t \iff x \ge F^-(t)$.

<!-- page 72 -->
> - **iii.** $F^-$ is left-continuous, that is
>   $$
>   \lim_{\substack{t \to t_0 \\ t < t_0}} F^-(t) = F^-(t_0).
>   $$
> - **iv.** For all $t \in \mathbb{R}_+$, $F(F^-(t)) \ge t$.
>
> ```
>        1 |
> F(F-(t)) |........•
>        t |......./
>          |      /
>        0 +-----+------------------------> x
>              F-(t)
> ```
>   Moreover, if $F$ is continuous, for all $t \in \mathbb{R}_+$, $F(F^-(t)) = t$.
> - **v.** For all $x \in \mathbb{R}$, $F^-(F(x)) \le x$.
>
> ```
>        1 |
>          |          --------------------
>     F(x) |........./
>          |        /
>        0 +-------+-----+----------------> x
>               F-(F(x)) x
> ```
>   Moreover, if $F$ is increasing, for all $x \in \mathbb{R}$, $F^-(F(x)) = x$.
> - **vi.** In particular, if $F$ is increasing and continuous, then $F$ is bijective with inverse function $F^{-1} = F^-$.

*Proof.* Denote for all $t \in \mathbb{R}_+$,
$$
A_t = \{x \in \mathbb{R} \; ; \; F(x) \ge t\},
$$
such that $F^-(t) = \inf(A_t)$.

- **i.** If $s \le t$, then, for all $x \in \mathbb{R}$, $F(x) \ge t \implies F(x) \ge s$. Hence, $A_t \subset A_s$ and thus:
  $$
  F^-(t) = \inf(A_t) \ge \inf(A_s) = F^-(s).
  $$
- **ii.**
  $\Rightarrow$ If $F(x) \ge t$, then $x \in A_t$ and thus $F^-(t) = \inf(A_t) \le x$.
  $\Leftarrow$ Assume $x \ge F^-(t)$. Let $\varepsilon > 0$. Then $x + \varepsilon > x \ge F^-(t)$, so $x + \varepsilon > \inf(A_t)$. In particular, there exists $x_\varepsilon \in A_t$ such that $x_\varepsilon \le x + \varepsilon$. Since $F$ is non-decreasing, then $F(x + \varepsilon) \ge F(x_\varepsilon) \ge t$, as $x_\varepsilon \in A_t$.

<!-- page 73 -->
  We proved that, for all $\varepsilon > 0$, $F(x + \varepsilon) \ge t$. Hence, since $F$ is right-continuous, as $\varepsilon \to 0^+$,
  $$
  F(x) = \lim_{\substack{\varepsilon \to 0 \\ \varepsilon > 0}} F(x + \varepsilon) \ge t.
  $$
- **iii.** Admitted here. Note that $F^-$ is not right-continuous as soon as $F$ is not increasing.
- **iv.** Applying ii. to $x = F^-(t)$ directly leads to $F(F^-(t)) = F(x) \ge t$.
  Assume now that $F$ is continuous. Consider $(x_n)_n$ a sequence satisfying:
  $$
  \forall n, \quad x_n < F^-(t) \quad \text{and} \quad \lim_{n \to +\infty} x_n = F^-(t).
  $$
  Then by contraposition of ii., for all $n$, $F(x_n) < t$, and by taking $n \to +\infty$, as $F$ is continuous (and in particular left-continuous):
  $$
  F(F^-(t)) = \lim_{n \to +\infty} F(x_n) \le t.
  $$
  Finally, in this case, both $F(F^-(t)) \le t$ and $F(F^-(t)) \ge t$ imply that $F(F^-(t)) = t$.
- **v.** Applying ii. to $t = F(x)$ leads to $F^-(F(x)) = F^-(t) \le x$.
  Assume now that $F$ is increasing. Assume (by contradiction) that $F^-(F(x)) < x$. Then
  $$
  F(F^-(F(x))) < F(x),
  $$
  which is impossible by the first point applied to $t = F(x)$.
- **vi.** Immediate from the previous point.

## III. Complementary proofs

### III.1. Proof of Lemma 1.13

Let $f : \mathbb{R}_+ \to \mathbb{R}$ be a continuous function satisfying
$$
f(x + y) = f(x)f(y) \quad \text{for all } x, y \in \mathbb{R}_+. \tag{A.1}
$$
Note that setting $x = y = 0$, we obtain $f(0) = f(0)^2$, hence
$$
f(0) = 0 \quad \text{or} \quad f(0) = 1.
$$

- **The case $f(0) = 0$:**
  For any $x$ in $\mathbb{R}_+$, $f(x) = f(x + 0) = f(x)f(0) = 0$, so $f$ is the null function on $\mathbb{R}_+$.

- **The case $f(0) = 1$:**
  We now assume $f(0) = 1$. Let us first prove two useful properties on such solution of (A.1).

<!-- page 74 -->
> **Property A.4.**
> Let $f : \mathbb{R}_+ \to \mathbb{R}$ be a function satisfying (A.1) such that $f(0) = 1$.
> 1. Fix $y$ in $\mathbb{R}_+$. Then, for all $n$ in $\mathbb{N}$,
>    $$
>    f(ny) = [f(y)]^n. \tag{A.2}
>    $$
> 2. **Positiveness:** For all $x$ in $\mathbb{R}_+$, $f(x) > 0$.

*Proof of Property A.4.*

1. **Proof of (A.2) by mathematical induction:**
   - **Base case:** (A.2) holds for $n = 0$ since
     $$
     f(0 \times y) = f(0) = 1 = [f(y)]^0.
     $$
   - **Induction step:** Assume that (A.2) is true for any fixed integer $n$. Then
     $$
     f((n + 1)y) = f(ny + y) = f(ny)f(y) = [f(y)]^n f(y) = [f(y)]^{n+1},
     $$
     that is (A.2) also holds for $n + 1$. This completes the induction step.

2. **Proof of the positiveness of $f$:**
   - **On a neighbourhood of 0:** Since $f$ is continuous at 0, and $f(0) = 1$, there exists $\delta > 0$ such that
     $$
     f(x) > 0, \quad \forall x \in [0, \delta].
     $$
   - **On its complementary:** Now let $x > \delta$. There exists $n$ in $\mathbb{N}^*$ such that $x/n \le \delta$. In particular, $f(x/n) > 0$. Therefore, by Equation (A.2):
     $$
     f(x) = f\left(n \times \frac{x}{n}\right) = \left[ f\left(\frac{x}{n}\right) \right]^n > 0.
     $$
   From both points above, we deduce that $f$ is positive-valued on $\mathbb{R}_+$.

Now, let us prove that there exists $a$ in $\mathbb{R}$ (depending on $f$) such that for all $x$ in $\mathbb{R}_+$,
$$
f(x) = e^{ax}. \tag{A.3}
$$
To do so, we follow the standard steps for solving a functional equation.

**Step 1: For integer arguments**
By Equation (A.2) for $y = 1$, we obtain that
$$
f(n) = [f(1)]^n = e^{an}, \quad \forall n \in \mathbb{N},
$$
where $a = \ln[f(1)]$ is well-defined since, by Property A.4, $f(1) > 0$. This proves (A.3) on $\mathbb{N}^*$.

<!-- page 75 -->
**Step 2: Extension to nonnegative rational numbers**
Let $z = m/n$ be a nonnegative rational number with $m$ in $\mathbb{N}$ and $n$ in $\mathbb{N}^*$. Then, again by (A.2):
$$
f(m) = f\left(n \times \frac{m}{n}\right) = [f(z)]^n.
$$
Since $f$ is positive, we deduce that
$$
f(z) = f(m)^{1/n} = (e^{am})^{1/n} = e^{am/n} = e^{az}.
$$
This proves (A.3) on $\mathbb{Q}_+$.

**Step 3: Extension to nonnegative real numbers (continuity required)**
Let $x$ in $\mathbb{R}_+$. Since $\mathbb{Q}_+$ is a dense subset of $\mathbb{R}_+$, there exists a sequence of nonnegative rational numbers $(z_n)_{n \ge 0}$ that converges to $x$. Hence, by continuity of both $f$, and $z \mapsto e^{az}$ on $\mathbb{R}_+$:
$$
f(x) = \lim_{n \to +\infty} f(z_n) = \lim_{n \to +\infty} e^{az_n} = e^{ax}.
$$
This proves (A.3) on $\mathbb{R}_+$.

Finally, we proved that the solutions of the functional Equation (A.1) on $\mathbb{R}_+$ are either the null function (if $f(0) = 0$) or positive functions $f : x \mapsto e^{ax}$, with $a = \ln[f(1)]$ in $\mathbb{R}$ (if $f(0) = 1$). This ends the proof of Lemma 1.13.

**Optional Step 4: Extension to negative real numbers**
Note that this result extends to functions defined on $\mathbb{R}$. Indeed, assume (A.1) also holds for negative $x$ and $y$ and that $f(0) = 1$. Then, all steps above still hold for nonnegative integers, rational and real numbers. Moreover, for all negative $x$, $1 = f(0) = f(x - x) = f(x)f(-x)$, thus, by (A.3) on $\mathbb{R}_+$:
$$
f(x) = \frac{1}{f(-x)} = \frac{1}{e^{-ax}} = e^{ax}.
$$

### III.2. Proof of Proposition 4.13

- On the one hand, on the event $\{Y = 0\}$, $N_t = 0$.
- On the other hand, on the event $\{Y = n\}$:
  $$
  N_t = \sum_{i=1}^n \mathbf{1}_{\{S_i \le t\}}.
  $$

To prove that $(N_t)_{t \in [0, T]}$ is a Poisson process, we need to prove that for all subdivisions of $[0, T]$:
$$
0 = t_0 < t_1 < \dots < t_{m-1} < t_m = T,
$$
and for all integers $k_1, \dots, k_m$, one has
$$
\mathbb{P}(N_{t_1} = k_1, N_{t_2} - N_{t_1} = k_2, \dots, N_{t_m} - N_{t_{m-1}} = k_m) = \prod_{j=1}^m \frac{c_j^{k_j}}{k_j!} e^{-c_j}, \tag{A.4}
$$
where
$$
c_j = \int_{t_{j-1}}^{t_j} \lambda(x)\,dx = \Lambda(t_j) - \Lambda(t_{j-1}).
$$

<!-- page 76 -->
Yet, since $N_{t_m} = N_T = Y$ (there are $Y$ $S_i$'s that belong to $[0, T]$):
$$
\begin{aligned}
&\mathbb{P}(N_{t_1} = k_1, N_{t_2} - N_{t_1} = k_2, \dots, N_{t_m} - N_{t_{m-1}} = k_m) \\
&= \mathbb{P}\left( N_{t_1} = k_1, N_{t_2} - N_{t_1} = k_2, \dots, N_{t_m} - N_{t_{m-1}} = k_m, \; Y = \sum_{j=1}^m k_j \right) \\
&= \mathbb{P}\left( N_{t_1} = k_1, N_{t_2} - N_{t_1} = k_2, \dots, N_{t_m} - N_{t_{m-1}} = k_m \;\middle|\; Y = \sum_{j=1}^m k_j \right) \mathbb{P}\left( Y = \sum_{j=1}^m k_j \right).
\end{aligned}
$$
On the one hand, for all $j$, for all $i$, by definition of the $S_i$'s:
$$
\mathbb{P}(S_i \in (t_{j-1}, t_j]) = \int_{t_{j-1}}^{t_j} \frac{\lambda(s)}{\Lambda(T)}\,ds = \frac{c_j}{c},
$$
where $c$ denotes
$$
c = \int_0^T \lambda(x)\,dx = c_1 + \dots + c_m = \Lambda(T).
$$
Hence, conditionally on $\{Y = \sum_{j=1}^m k_j\}$, the random variables $\left( \sum_{i=1}^{k_1 + \dots + k_m} \mathbf{1}_{\{S_i \in (t_{j-1}, t_j]\}} \right)_{1 \le j \le m}$, which count the number of $S_i$'s in each interval $(t_{j-1}, t_j]$, have a multinomial distribution with parameters $(k_1 + \dots + k_m)$ and $(c_j / c)_{1 \le j \le m}$.

Therefore,
$$
\begin{aligned}
&\mathbb{P}\left( N_{t_1} = k_1, N_{t_2} - N_{t_1} = k_2, \dots, N_{t_m} - N_{t_{m-1}} = k_m \;\middle|\; Y = \sum_{j=1}^m k_j \right) \\
&= \mathbb{P}\left( \sum_{i=1}^{\sum k_j} \mathbf{1}_{\{S_i \in (t_0, t_1]\}} = k_1, \dots, \sum_{i=1}^{\sum k_j} \mathbf{1}_{\{S_i \in (t_{m-1}, t_m]\}} = k_m \right) \\
&= \frac{(k_1 + k_2 + \dots + k_m)!}{k_1! \dots k_m!} \left(\frac{c_1}{c}\right)^{k_1} \dots \left(\frac{c_m}{c}\right)^{k_m}.
\end{aligned}
$$
On the other hand, $Y \sim \mathcal{P}(\Lambda(T))$, with $\Lambda(T) = c = c_1 + \dots + c_m$. Finally,
$$
\begin{aligned}
&\mathbb{P}(N_{t_1} = k_1, N_{t_2} - N_{t_1} = k_2, \dots, N_{t_m} - N_{t_{m-1}} = k_m) \\
&= \frac{(k_1 + k_2 + \dots + k_m)!}{k_1! \dots k_m!} \left(\frac{c_1}{c}\right)^{k_1} \dots \left(\frac{c_m}{c}\right)^{k_m} \frac{c^{\sum_{j=1}^m k_j}}{(k_1 + \dots + k_m)!} e^{-c} \\
&= \frac{c_1^{k_1}}{k_1!} e^{-c_1} \times \dots \times \frac{c_m^{k_m}}{k_m!} e^{-c_m},
\end{aligned}
$$
which ends the proof of (A.4).

<!-- page 77 -->
# Appendix B. Tutorials and Computer Lab

<!-- page 78 -->
**INSA de Toulouse - 5 ModIA**
**Poisson processes and application**
Academic year 2026-2027
**Worksheet 1**

## Worksheet 1 - Some reminders in probability

### Exercise 1. Memoryless property and Defective cellphones.

**(A) Generalization of the memoryless property.**
Let $X$ be an exponentially distributed r.v. $X \sim \mathcal{E}(\lambda)$, with parameter $\lambda > 0$, and $U$ be a continuous real-valued random variable with density $g$. Assume $X \perp\!\!\!\perp U$. Prove that for all $t, s > 0$,
$$
\mathbb{P}(X + U > t + s \mid X > t) = \mathbb{P}(X + U > s).
$$
*Hint: One may first prove that*
$$
\mathbb{P}(X + U > t + s \mid X > t) = \frac{\mathbb{P}(\{X > t + s - U\} \cap \{s - U > 0\})}{\mathbb{P}(X > t)} + \frac{\mathbb{P}(\{X > t\} \cap \{s - U \le 0\})}{\mathbb{P}(X > t)}.
$$

**(B) Cellphone warranty.** A cellphone manufacturer offers warranty for 8 months with a free replacement if the device is defective. The manufacturer assumes that each unit has a lifetime that can be modeled by an exponential r.v. $X$ with mean equal to 24 months.

1. What is the (theoretical) proportion of defective devices that will be replaced for free?
2. What warranty should the manufacturer offer in order to replace only 10% of the devices for free?
3. Compute the (conditional) distribution function defined by
   $$
   H : t \mapsto \mathbb{P}(X \le t \mid X \ge 8).
   $$
   Deduce that the density function for duration of devices that do not have any failure before 8 months equals
   $$
   h : t \mapsto \lambda e^{-\lambda(t - c)} \mathbf{1}_{\{t \ge c\}},
   $$
   where $c, \lambda$ are parameters to determine. Note that this distribution is called "translated exponential distribution".
4. Charles is a VIP member, he has infinite warranty. He can replace his phone every time it fails. The manufacturer would like to know the probability that 3 phones will not cover 5 years if Charles' first cellphone lasts more than a year. Compute this conditional probability assuming the lifetimes of Charles' cellphones are i.i.d. exponential r.v.
   *Hint: One may use the generalized memoryless property.*

### Exercise 2. Car batteries.

The lifetime (in years) of an automobile battery is described by a r.v. $X$ having an exponential distribution with parameter $\lambda = 1/5$.

1. Determine the expected lifetime of the battery and its standard deviation.
2. Compute the probability that the battery lasts for at least 10 years.
3. If the battery has lasted for 5 years, what is the probability that it will last for at least an additional 5 years?

<!-- page 79 -->
**INSA de Toulouse - 5 ModIA**
**Poisson processes and application**
Academic year 2026-2027
**Worksheet 1**

### Exercise 3. The Gamma distribution.

Let $\alpha, \lambda > 0$ and consider $X \sim \Gamma(\alpha, \lambda)$.

1. Determine the domain and compute the Laplace transform $L_X : u \mapsto \mathbb{E}[e^{uX}]$ of $X$.
2. Deduce $\mathbb{E}[X]$ and $\text{Var}(X)$.
3. Let $\beta > 0$ and consider $Y \sim \Gamma(\beta, \lambda)$ independent of $X$. Show that $X + Y \sim \Gamma(\alpha + \beta, \lambda)$.
4. Deduce the distribution of the sum of $n$ i.i.d. exponential random variables with parameter $\lambda$.

### Exercise 4. The Weibull distribution.

With the exponential distribution, the Weibull distribution is one of the most widely used lifetime distributions in reliability. Let $\alpha, \beta > 0$. A non-negative real-valued random variable $X$ with density (w.r.t. the Lebesgue measure)
$$
x \mapsto \frac{\beta}{\alpha} \left(\frac{x}{\alpha}\right)^{\beta - 1} \exp\left( -\left(\frac{x}{\alpha}\right)^\beta \right) \mathbf{1}_{\{x > 0\}},
$$
is said to have a Weibull distribution with parameter $(\alpha, \beta)$. We denote $X \sim \mathcal{W}(\alpha, \beta)$. Note that if $\beta = 1$, we recover an exponential density with parameter $1/\alpha$.

1. Prove that the survival function of $X$ is $G : t \in \mathbb{R}_+ \mapsto \exp\left( -\left(\frac{t}{\alpha}\right)^\beta \right)$.
2. Deduce that the hazard rate of $X$ equals $h(t) = \frac{\beta}{\alpha} \left(\frac{t}{\alpha}\right)^{\beta - 1}$.
3. Prove that $\left(\frac{X}{\alpha}\right)^\beta$ has an exponential distribution with parameter 1.
4. Show that if $X_1, \dots, X_n$ are $n$ independent r.v. with distribution $\mathcal{W}(\alpha, \beta)$, then
   $$
   \min_{1 \le i \le n} X_i \sim \mathcal{W}\left( \frac{\alpha}{n^{1/\beta}}, \beta \right).
   $$

<!-- page 80 -->
**INSA de Toulouse - 5 ModIA**
**Poisson processes and application**
Academic year 2026-2027
**Worksheet 2**

## Worksheet 2 - Homogeneous Poisson process

### Exercise 1. Defective cellphones (point process approach).

Consider Part (B) "Cellphone warranty" of Exercise 1, Worksheet 1.
A cellphone manufacturer offers warranty for 8 months with a free replacement if the device is defective. The manufacturer assumes that each unit has a lifetime that can be modeled by an exponential variable with mean equal to 24 months.
Recall that Charles is a VIP member, he has infinite warranty. He can replace his phone every time it fails. The manufacturer would like to know the probability that 3 phones will not cover 5 years if Charles' first cellphone lasts more than a year. Compute this conditional probability assuming the number of times Charles replaces his cellphone can be modeled by a homogeneous Poisson process.

### Exercise 2. The inspection paradox (or Le paradoxe de l'autobus in French).

You have an infinite stock of light bulbs with lifetimes which can be modeled by i.i.d. random variables with exponential distribution with parameter $\lambda > 0$. At time $t = 0$, we switch on the light and we change the light bulb as soon as it burns out. The times we change the light bulbs
$$
0 < T_1 < T_2 < \dots < T_n < \dots
$$
thus define a Poisson process $(N_t)_{t \ge 0}$ with rate $\lambda$, where $N_t$ denotes the number of consumed light bulbs by time $t$. We decide to approximate the parameter $\lambda$ by the following "inspection" method:
- Fix a time $t > 0$.
- Observe the lifetime of the "inspected" light bulb at time $t$, that is $T_{N_t + 1} - T_{N_t}$.

We denote $X_t = T_{N_t + 1} - t$ the residual lifetime of the inspected bulb at time $t$, and $Y_t = t - T_{N_t}$ the "age" of the inspected bulb at time $t$.

1. Prove that $X_t \sim \mathcal{E}(\lambda)$.
2. Compute the survival function of $Y_t$. Is it exponentially distributed?
3. Prove that $X_t$ and $Y_t$ are independent random variables.
4. Compute the average lifetime of the inspected bulb at time $t$.
5. Compare it to the average lifetime of the first bulb.

We can see that, for all time $t > 0$, the expected lifetime of the inspected bulb is greater than the average lifetime of the non-inspected bulbs. This phenomenon is often referred to the **inspection paradox**.

<!-- page 81 -->
**INSA de Toulouse - 5 ModIA**
**Poisson processes and application**
Academic year 2026-2027
**Worksheet 2**

### Exercise 3. Small/Large claims.

The time unit considered here is a week. A small insurance company has in average one claim per day. It classifies each claim in two categories:
- small claim if the claim is less than or equal to 5000 €,
- large claim if the claim is greater than 5000 €.

It estimates that large claims occur in 10% of the cases. We model the claim occurrences by a homogeneous Poisson process $N$ with rate $\lambda > 0$, and denote $N^I$ (resp. $N^{II}$) the point process counting the small claims (resp. large claims).

1. Determine the rate $\lambda$.
2. What is the probability of at least 1 large claim in a week?
3. Knowing there has been 5 large claims in 2 weeks, what is the expected number of claims in that interval?
4. There has been 30 accidents in February. What is the probability that 3 of them were large claims?

<!-- page 82 -->
**INSA de Toulouse - 5 ModIA**
**Poisson processes and application**
Academic year 2026-2027
**Worksheet 3**

## Worksheet 3 - Statistics for homogeneous Poisson process

### Exercise 1. Quality control.

A manufacturing plant wants to check the quality of its production. To do so, the person in charge keeps count of the number of defective parts.
He thus denotes $N_t$ the number of defective parts manufactured by time $t$ (for any $t$ in hours) and assumes $(N_t)_{t \in \mathbb{R}_+}$ is a Poisson process with unknown intensity $\lambda > 0$. He believes that the quality of the production is not acceptable if the average number of defective parts per hour is greater than or equal to 5.

1.
   - **(a)** Do you agree with the Poisson process assumption?
   - **(b)** We aim at testing
     $$
     H_0 : \lambda \ge \lambda_0 \quad \text{against} \quad H_1 : \lambda < \lambda_0.
     $$
     Explicit the value of $\lambda_0$ and justify the choice of these testing hypotheses.

2. The person in charge fixes the inspection time to one day and observes $(N_t)_{t \in [0, 24]}$. By the end of the day, he counts 108 defective parts.
   - **(a)** Determine an estimation of the rate $\lambda$.
   - **(b)** According to you, is the production quality satisfactory?

3. Now, the person in charge decides to observe up to the hundredth flaw and notices that it appears after 22 hours and 15 minutes of production.
   - **(a)** What is an estimation of $\lambda$ in this case?
   - **(b)** First, we may consider that 100 points is large enough to work with the asymptotic distribution. Given this new observation, is the production quality satisfactory?
   - **(c)** Second, we do not consider the asymptotic distribution anymore. What is your conclusion in this case?

Let us recall the quantiles of the standard gaussian distribution $\mathcal{N}(0, 1)$ and the $\eta$-quantiles $x_{d, \eta}$ of the $\chi^2(d)$ distribution:

If $Z \sim \mathcal{N}(0, 1)$:
- $\mathbb{P}(Z \le 1.645) = 0.95$
- $\mathbb{P}(Z \le 1.96) = 0.975$

Table of $\chi^2(d)$ quantiles:

| $x_{d, \eta}$ | $\eta = 0.05$ | $\eta = 0.95$ |
| :--- | :---: | :---: |
| $d = 50$ | 34.764 | 67.505 |
| $d = 100$ | 77.929 | 124.342 |
| $d = 200$ | 168.279 | 233.994 |

### Exercise 2. Unbiased estimation for a fixed number of observations.

In this exercise, we assume that $n \ge 2$ is fixed, and we observe a homogeneous Poisson process with rate $\lambda > 0$, up to the $n$th event. Denote $0 < T_1 < T_2 < \dots < T_n$ the corresponding arrival times. We estimate the unknown parameter $\lambda$ by the maximum likelihood estimator (MLE) $\hat{\lambda}_n = n/T_n$.

1. Prove that the MLE is a biased estimator of the rate $\lambda$.
2. Consider $\tilde{\lambda}_n = (n - 1)/T_n$.
   - **(a)** Check that $\tilde{\lambda}_n$ is an unbiased estimator of $\lambda$.
   - **(b)** Is it efficient?

<!-- page 83 -->
**INSA de Toulouse - 5 ModIA**
**Poisson processes and application**
Academic year 2026-2027
**Worksheet 4**

## Worksheet 4 - Inhomogeneous Poisson process

### Exercise 1. Customer's arrival time at a post office.

We count the number of customers who arrive at a post office that opens at 8:00 a.m. and closes at 5:00 p.m.
- From 8:00 a.m. to 11:00 a.m., customers seem to arrive, on average, at a steadily increasing rate that starts with an initial rate of 5 customers per hour at 8:00 a.m., and reaches a maximum of 20 customers per hour at 11:00 a.m.
- From 11:00 a.m. to 1:00 p.m., the (average) rate seems to remain constant at 20 customers per hour.
- Then, the (average) arrival rate drops steadily from 1:00 p.m. until closing time at which time, it has the value of 12 customers per hour.

1. What is a good probability model? Write the intensity.
2. What is the probability that no customer arrives between 8:30 a.m. and 9:30 a.m.?
3. What is the expected number of arrivals between 8:30 a.m. and 9:30 a.m.?

### Exercise 2. Repairable System (from the 2021/2022 Exam).

We study a repairable system for which:
- at each failure, the system is repaired to its condition just before failure,
- the repair times can be neglected.

We model the failure times by an inhomogeneous Poisson process with the following intensity:
$$
\forall x \ge 0, \quad \lambda(x) = \begin{cases}
10 - 4x & \text{if } 0 \le x < 2, \\
2 & \text{if } 2 \le x < 12, \\
x - 10 & \text{if } x \ge 12.
\end{cases}
$$
Denote $N_t$ the number of failures by time $t > 0$ (in months).

1. Recall the definition of an inhomogeneous Poisson process with intensity function $\lambda(\cdot)$.
2. Represent the intensity function and comment its shape.
3. Compute the cumulative intensity function $\Lambda(\cdot)$ of the Poisson process.
   *(Indication: you may distinguish 3 cases and check that $\Lambda(2) = 12$ and $\Lambda(12) = 32$).*
4. Compute the probability that at least one failure occurs during the first month.
5. Compute the expected number of failures between 2 and 12 months.
6. Compute the probability that at least two failures occur by the third month knowing that there was exactly one failure during the first two months.

<!-- page 84 -->
**INSA de Toulouse - 5 ModIA**
**Poisson processes and application**
Academic year 2026-2027
**Worksheet 4**

### Exercise 3. Parametric inhomogeneous Poisson process (fixed $T$).

Let $N$ be an inhomogeneous Poisson process with intensity function
$$
\lambda_\alpha : \begin{array}{ccc}
\mathbb{R}_+^* & \longrightarrow & \mathbb{R}_+ \\
x & \longmapsto & \frac{\alpha}{2\sqrt{x}}
\end{array},
$$
which depends on a positive parameter $\alpha$.

1. Compute the cumulative intensity $\Lambda_\alpha$ of the Poisson process $N$.
2. For given $\alpha$ and $T$, write an algorithm for simulating an inhomogeneous Poisson process with intensity function $\lambda_\alpha$ on $[0, T]$.
3. Let us now assume that $\alpha$ is unknown, and that we observe $N = (N_t)_{t \in [0, T]}$. Denote
   $$
   T_1 < T_2 < \dots < T_{N_T}
   $$
   its corresponding arrival times.
   - **(a)** Compute the maximum likelihood estimator (MLE) of $\alpha$.
   - **(b)** Is the MLE a biased estimator of $\alpha$?
   - **(c)** Is the MLE an efficient estimator of $\alpha$?

### Exercise 4. Parametric inhomogeneous Poisson process (fixed $n$).

Let $N$ be an inhomogeneous Poisson process with intensity function
$$
\lambda_\theta : \begin{array}{ccc}
\mathbb{R}_+^* & \longrightarrow & \mathbb{R}_+ \\
x & \longmapsto & \theta e^{-x}
\end{array},
$$
which depends on a positive parameter $\theta$. Denote $T_1 < T_2 < \dots < T_n$ its $n$ first arrival times.

1. Compute the cumulative intensity $\Lambda_\theta$ of the Poisson process $N$.
2. Compute the maximum likelihood estimator (MLE) of $\theta$.
3. Construct a non-asymptotic confidence interval for $\theta$ with confidence level 95%.
4. For given $\theta$ and $n$, write an algorithm for simulating the $n$ first arrival times of an inhomogeneous Poisson process with intensity function $\lambda_\theta$.

<!-- page 85 -->
**INSA de Toulouse**
**5 ModIA**
**Poisson processes and application 2026-2027**

# Introduction to Poisson processes with R

## Objectives

The aim of this session is to manipulate and illustrate the notions introduced in the lecture on Poisson processes with the R software.

Load the R file `TP-PoissonProcess-etud.Rmd` (to be completed), available on the Moodle page, and fill in the gaps during this session.

## 1. Homogeneous Poisson processes observed on a fixed window

First, we consider the case of a fixed observation window (and thus a random number of points).

### 1.1. Simulation

- **i.** Recall the conditional distribution of the arrival times of a homogeneous Poisson process with rate `lambda` on an interval `[0, Tmax]`, given the number of points in that interval.
- **ii.** Write an R function `simulPPh1` with arguments `lambda` and `Tmax`, which simulates such a process and returns the corresponding arrival times.
- **iii.** For `lambda=2` and `Tmax=10`, simulate a homogeneous Poisson process and plot both the counting process and the arrival times.
  *Indication:* you may understand and use the R function defined in the R file `plot.PP()`.

### 1.2. Maximum likelihood estimator

- **i.** Write an R function `MLE1` which returns the Maximum Likelihood Estimator of a homogeneous Poisson process `PPh` observed on a fixed window `[0, Tmax]`.
- **ii.** Apply it on different simulated data. What do you observe?

### 1.3. Asymptotic behavior of the MLE

In this section, we illustrate the asymptotic behavior[^1], as $T$ goes to $+\infty$, of the MLE $\hat{\lambda}_T$ of the rate $\lambda$ of a homogeneous Poisson process on $[0, T]$.

#### 1.3.1. Strong LLN-type result

First let us illustrate the almost-sure convergence:
$$
\hat{\lambda}_T \xrightarrow[T \to +\infty]{\text{a.s.}} \lambda. \tag{1}
$$
Note that in the lecture, we considered a weak LLN-type result since we only proved the convergence in probability. To illustrate (1):
- **i.** Fix a rate `lambda=2` and a sequence of `Tmax` that tends to $+\infty$, say `Tillustr=1:500`.
- **ii.** For each `Tmax` in `Tillustr`, simulate a homogeneous Poisson process with rate `lambda` on `[0, Tmax]`, and compute the MLE.
- **iii.** Plot the obtained values for the MLE as a function of `Tmax`.
  What do you observe?

[^1]: LLN refers to "Law of Large Numbers" and CLT stands for "Central Limit Theorem".

<!-- page 86 -->
**INSA de Toulouse**
**5 ModIA**
**Poisson processes and application 2026-2027**

#### 1.3.2. CLT-type result

We now illustrate the following result:
$$
\sqrt{T}(\hat{\lambda}_T - \lambda) \xrightarrow[T \to +\infty]{\mathcal{L}} \mathcal{N}(0, \lambda). \tag{2}
$$
To do so, fix `lambda=2` and do the following for different values of `Tmax`.

- **i.** Fix the number of simulations `K=1000` and create a vector of size $K$:
  ```R
  Z = rep(0, K)
  ```
- **ii.** Store in `Z` a $K$-sample with same distribution as $\sqrt{T}(\hat{\lambda}_T - \lambda)$:
  ```R
  for(k in 1:K)
  {
    pph = simulPPh1(lambda, Tmax)
    mle = MLE1(pph, Tmax)
    Z[k] = sqrt(Tmax) * (mle - lambda)
  }
  ```
- **iii. With the density function:** Plot the histogram of the sample `Z` (which approximates the density of the $Z_k$'s), and compare it with the density of the limit distribution $\mathcal{N}(0, \lambda)$:
  ```R
  hist(Z, freq=FALSE, main=paste("Tmax", Tmax, sep="="))
  curve(dnorm(x, mean=0, sd=sqrt(lambda)), col="red", add=TRUE)
  ```
  What do you observe when `Tmax` equals 1, 10, 100 and 500?
- **iv. With the c.d.f.:** Plot the empirical cumulative distribution function of `Z`, and compare it with the cumulative distribution function of a $\mathcal{N}(0, \lambda)$:
  ```R
  plot(ecdf(Z), main=paste("Tmax", Tmax, sep="="))
  curve(pnorm(x, mean=0, sd=sqrt(lambda)), col="red", lwd=2, add=TRUE)
  ```
  What do you observe when `Tmax` equals 1, 10, 100 and 500?

### 1.4. Statistical inference : hypothesis testing

Consider a given rate $\lambda_0 > 0$ to test. Given the observation of a homogeneous Poisson process with (unknown "true") rate $\lambda$ observed on a fixed window $[0, T]$, we aim at testing
$$
H_0 : \lambda = \lambda_0 \quad \text{against} \quad H_1 : \lambda \ne \lambda_0.
$$

- **i.** Construct a test of $H_0$ against $H_1$, and express the corresponding p-value in terms of the standard gaussian c.d.f.
- **ii.** Write an R function `test1`, with arguments the observed homogeneous Poisson process `PPh`, the observation time `Tmax` and the rate `lambda0` to test, and which returns the p-value of the test constructed in the previous question.
- **iii.** To validate this test on simulated data, we estimate for different values of $\lambda$, the probability (or proportion)
  $$
  p(\lambda) = \mathbb{P}_\lambda(\text{reject } H_0).
  $$
  In particular, if $\lambda$ satisfies $H_0$ (i.e. $\lambda = \lambda_0$), then $p(\lambda)$ is the size of the test, and otherwise, it is the power of the test against the alternative $\lambda$.
  The R function `plot.level.power1` in the file `TP-PoissonProcess-etud.R` plots for different values of the rate $\lambda$ (in `TrueLambda`) confidence intervals for the proportion $p(\lambda)$.
  Now, fix `alpha=0.05`, `nsimu=1000` and let `lambda0=2` and `TrueLambda=c(1, 1.5, 1.8, 1.9, 2, 2.1, 2.2, 2.5, 3)`.
  Apply this function for `Tmax = 1, 10, 100` and `500`. Understand and comment the obtained graphs.

<!-- page 87 -->
**INSA de Toulouse**
**5 ModIA**
**Poisson processes and application 2026-2027**

## 2. Inhomogeneous Poisson processes on a fixed window

We now aim at simulating inhomogeneous Poisson processes with given intensity function, on a given interval.

- **i.** Recall the thinning algorithm. In particular, when does it apply?
- **ii.** Write an R function `simulPPi` with arguments the intensity function `lambda_fct`, `Tmax` and an upper bound `M`, which simulates an inhomogeneous Poisson process with intensity function `lambda_fct()` on `[0, Tmax]`, and returns the corresponding arrival times.
- **iii.** Simulate and represent inhomogeneous Poisson processes on $[0, 10]$ with intensity functions:
  - $\lambda_1 : x \mapsto 2x$.
  - $\lambda_2 : x \mapsto \frac{3}{10} x^2$.
  - $\lambda_3 : x \mapsto 100 x e^{-x}$.
  - $\lambda_4 : x \mapsto 10 \times \mathbf{1}_{[0, 7]}(x) + 15 \times \mathbf{1}_{]8, 10]}(x)$.
  - $\lambda_5 : x \mapsto \dots$.

## 3. Homogeneous Poisson processes with fixed number of points

Second, we consider the case of a fixed number of points (and thus a random observation window).

### 3.1. Simulation

- **i.** Recall the distribution of the interarrival times of a homogeneous Poisson process with rate `lambda`.
- **ii.** Write an R function `simulPPh2` with arguments the rate `lambda` and the number of points `n`, which simulates such a process and returns the corresponding arrival times.
- **iii.** As in the "fixed window" case, fix `lambda=2` and `n=20`, simulate a homogeneous Poisson process and plot both the counting process and the arrival times.

### 3.2. Maximum likelihood estimator

- **i.** Write an R function `MLE2` which returns the maximum likelihood estimator of the rate of a homogeneous Poisson process `PPh` observed up to the $n$th point.
- **ii.** Apply it on different simulated data. What do you observe?

### 3.3. Asymptotic behavior of the MLE

As in section 1.3, illustrate the asymptotic behavior of the MLE $\hat{\lambda}_n$ of a homogeneous Poisson process with rate $\lambda$ observed up to the $n$th point as $n$ goes to $+\infty$, that is:

- **i. (Strong) LLN-type result:**
  $$
  \hat{\lambda}_n \xrightarrow[n \to +\infty]{\text{a.s.}} \lambda.
  $$
- **ii. CLT-type result:**
  $$
  \sqrt{n}(\hat{\lambda}_n - \lambda) \xrightarrow[n \to +\infty]{\mathcal{L}} \mathcal{N}(0, \lambda^2).
  $$

### 3.4. Statistical inference : confidence intervals

- **i.** Recall both asymptotic and non-asymptotic confidence intervals for the (unknown) rate `lambda` of a homogeneous Poisson process `PPh` observed up to the $n$th point.
- **ii.** Write an R function with arguments `PPh`, `alpha=0.05` and a boolean `asymptotic`, and which, depending on the value of `asymptotic` (`TRUE` or `FALSE`), returns the corresponding confidence interval for `lambda`.
- **iii.** As for the testing problem, we want to validate the confidence intervals on simulated data.
  To do so, illustrate, in both asymptotic and non-asymptotic cases, the fact that the proportion of times the rate belongs to the obtained confidence interval is larger than the fixed confidence level `1 - alpha`.
  You may fix `lambda=2`, and consider `n=10` or `n=100`.

<!-- page 88 -->
**INSA de Toulouse**
**5 ModIA**
**Academic year 2025-2026**
**Poisson processes and application**

# Written exam (Friday, October 10, 2025)
**Duration: 1 hour 30**

No documents nor calculator are allowed. The exercises are independent. The presentation and the justifications will be taken into account in the assessment. We consider the notation of the lecture notes. You may answer either in French or in English but not both.

---

### Exercise 1. (45 min) Homogeneous Poisson processes.

Insurance claims are made at times $(T_i)_{i \ge 1}$ distributed according to a homogeneous Poisson process with rate $\lambda$. For simplicity, the successive claim amounts are supposed to be constant, denoted by $c$. The total discounted cost of all claims made up to time $t$, denoted $D(t)$, is defined by
$$
D(t) = \sum_{i=1}^{N_t} c e^{-\rho T_i},
$$
where $\rho$ is the discount rate and $N_t$ is the number of claims made by time $t$ (in days, starting in September). We aim at computing the expected value of $D(t)$.

1. Do you agree with the homogeneous Poisson process assumption? Justify your answer.
2. The manager observed that the 50th claim in September arrived on the 11th day of the month. He wants to test that, in average, there are 5 claims per day.
   - **(a)** Construct a non-asymptotic test of
     $$
     H_0 : \lambda = \lambda_0 \quad \text{against} \quad H_1 : \lambda \ne \lambda_0,
     $$
     specifying the value of $\lambda_0$, the test statistic and the critical region.
   - **(b)** Compute the observed test statistic and conclude (in the context) at level $\alpha = 5\%$ using the table below of $\eta$-quantiles $x_{d, \eta}$ of $\chi^2(d)$ distributions.

| $x_{d, \eta}$ | $\eta = 0.025$ | $\eta = 0.05$ | $\eta = 0.95$ | $\eta = 0.975$ |
| :--- | :---: | :---: | :---: | :---: |
| $d = 11$ | 3.86 | 4.57 | 10.34 | 21.92 |
| $d = 22$ | 10.98 | 12.34 | 21.34 | 36.78 |
| $d = 25$ | 13.12 | 14.61 | 37.65 | 40.65 |
| $d = 50$ | 32.36 | 34.76 | 67.50 | 71.42 |
| $d = 100$ | 74.22 | 77.93 | 124.34 | 129.56 |

3. **Computation of the expected total discounted cost.**
   - **(a)** Recall the conditional distribution of $(T_1, \dots, T_n)$ given that $N_t = n$.
   - **(b)** Deduce that
     $$
     \mathbb{E}[D(t) \mid N_t = n] = c n \mathbb{E}[e^{-\rho U}],
     $$
     carefully justifying each step of the calculation, and specifying the distribution of $U$.
   - **(c)** Prove that
     $$
     \mathbb{E}[e^{-\rho U}] = \frac{1 - e^{-\rho t}}{\rho t}.
     $$
   - **(d)** By deconditioning the conditional expectation, prove that
     $$
     \mathbb{E}[D(t)] = c \frac{1 - e^{-\rho t}}{\rho t} \mathbb{E}[N_t],
     $$
     and deduce the expected total discounted cost $\mathbb{E}[D(t)]$.

<!-- page 89 -->
**INSA de Toulouse**
**5 ModIA**
**Academic year 2025-2026**
**Poisson processes and application**

### Exercise 2. (45 min) Parametric inhomogeneous Poisson process.

Let $N_t$ denote the number of failures by time $t$ (in hours). Assume that $N = (N_t)_{t \ge 0}$ is an inhomogeneous Poisson process with intensity function
$$
\lambda_\theta : \begin{array}{ccc}
\mathbb{R}_+^* & \longrightarrow & \mathbb{R}_+ \\
x & \longmapsto & 2\theta x
\end{array},
$$
which depends on a positive parameter $\theta$.

1. Recall the definition of an inhomogeneous Poisson process with intensity function $\lambda_\theta$.
2. Compute the cumulative intensity $\Lambda_\theta$ of the Poisson process $N$.
3. According to you, which plot corresponds to a realization of $N$? Justify.

[FIGURE page 89: Trois trajectoires simulées de processus de comptage (Plot 1, Plot 2, Plot 3) sur $[0, T]$ affichant $N_t$ en fonction de $t$. L'intensité étant $\lambda_\theta(x) = 2\theta x$, celle-ci est strictement croissante avec le temps, ce qui implique que la fréquence des sauts s'accélère au cours du temps (les temps d'inter-arrivée deviennent de plus en plus courts).]

4. Compute the probability that at least five failures will occur during the first three hours, given that there were four failures during the first two hours.
5. Show that $\mathbb{E}[N_{t+1} - N_t] = (2t + 1)\mathbb{E}[N_1]$.
6. For given $\theta$ and $T$, write an algorithm for simulating an inhomogeneous Poisson process with intensity function $\lambda_\theta$ on $[0, T]$.
7. Let us now assume that $\theta$ is unknown, and that we observe $N = (N_t)_{t \in [0, T]}$. Denote
   $$
   T_1 < T_2 < \dots < T_{N_T}
   $$
   its corresponding arrival times.
   - **(a)** Compute the maximum likelihood estimator (MLE) of $\theta$. To do so, one may recall that:
     - in the fixed $T$ case, the likelihood of the process $(N_t)_{t \in [0, T]}$ at $\theta$ equals
       $$
       L((N_t)_{t \in [0, T]}; \theta) = \left( \prod_{i=1}^{N_T} \lambda_\theta(T_i) \right) e^{-\int_0^T \lambda_\theta(x)\,dx},
       $$
     - in the fixed $n$ case, the likelihood of the process $(T_n)_{n \ge 1}$ at $\theta$ equals
       $$
       L((T_i)_{1 \le i \le n}; \theta) = \left( \prod_{i=1}^n \lambda_\theta(T_i) \right) e^{-\int_0^{T_n} \lambda_\theta(x)\,dx}.
       $$
   - **(b)** Is the MLE a biased estimator of $\theta$?
   - **(c)** Is the MLE an efficient estimator of $\theta$?

<!-- page 90 -->
