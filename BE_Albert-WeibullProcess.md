<!-- page 1 -->

INSA de Toulouse, 5ModIA  
**Poisson processes and application**

Academic year 2026-2027  
**Supervisor:** Mélisande ALBERT

# The Power law process and reliability

**Keywords:** Inhomogeneous Poisson process, Weibull distribution, Maximum Likelihood Estimator.

## Context of the project

Inhomogeneous Poisson processes provide nice models in Reliability theory. For instance, if at each failure a system is repaired to its condition at the time of failure and placed in service again, then the failures are often modeled by an inhomogeneous Poisson process, provided the repair times can be neglected. The power law process, sometimes called the Weibull process, is very popular in Reliability theory. It is an inhomogeneous Poisson process with intensity given by

$$
\lambda(s) = \frac{\beta}{\alpha}\left(\frac{s}{\alpha}\right)^{\beta-1},
$$

where $\alpha$ and $\beta$ are two positive parameters.

The main purpose of this project is to study the power law process and statistical inference adapted to this model (maximum likelihood estimator (MLE), confidence intervals, goodness of fit test).

## Work to be performed

### Theoretical results

First, the students are asked to study Section 3.2.3, or Section 3.3.2 of the book [Cocozza-Thivent, 1997] which is dedicated to statistical inference for the power law process, and complete some proofs and calculations. Then, the students may construct a goodness-of-fit test based on Theorem 1 of [Moeller, 1976].

### Numerical results

Second, the students are asked to perform numerical simulations in order to validate, from a practical point of view, the confidence intervals presented in [Cocozza-Thivent, 1997], and the goodness-of-fit test based on [Moeller, 1976]. Illustrations of the different convergences of the MLE are also expected. Finally, [Proschan, 1963] presents the intervals between successive failures of the air conditioning system of each member of a fleet of Boeing 720 jet airplanes. Can we model the counting process by Power law processes? If so, provide confidence intervals for the parameters.

The numerical simulations, will be presented with R Markdown.

## References

[Cocozza-Thivent, 1997] Cocozza-Thivent, C. (1997). *Processus stochastiques et fiabilité des systèmes*, volume 28. Springer Science & Business Media.

[Moeller, 1976] Moeller, S. K. (1976). *The Rasch-Weibull process*. *Scandinavian Journal of Statistics*, pages 107–115.

[Proschan, 1963] Proschan, F. (1963). *Theoretical explanation of observed decreasing failure rate*. *Technometrics*, 5(3):375–383.
