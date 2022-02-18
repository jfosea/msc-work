Summary of A Natural Robustification of the Ordinary Instrumental
Variables Estimator
================

*[Link](https://onlinelibrary.wiley.com/doi/10.1111/biom.12043) to the
paper.*

# 1 Problem

-   Ordinary Least Squares (OLS) estimates are biased and inconsistent
    when one the covariates is correlated with the error terms
    (*endogenous*). This happens when (1) covariate is associated with
    unobserved predictors, (2) is measured with error, and/or (3)
    affected by the response variable.

-   Solution to this problem is to use Ordinary Instrumental Variables
    (OIV) to consistently estimate the parameters of the model. However,
    when used in real datasets, OIV’s are often distorted due to
    outliers in the response variable.

-   There are existing robust methods to obtain IV estimators as listed
    below

    1.  Minimize robust loss function
        [(Amemiya, 1982)](https://www.jstor.org/stable/1912608?seq=1#metadata_info_tab_contents)
    2.  Robustify the estimating equations which define the OIV
        [(Krasker and
        Welsch, 1985)](https://www.jstor.org/stable/1913223?seq=1#metadata_info_tab_contents)
    3.  Orthogonality conditions which give rise to the OIV in the
        Generalized Method of Moments
    4.  Robustify at each of the stages of the OIV since OIV’s can be
        thought of as two-stage least squares estimator (2SLS).
    5.  Robust estimators for simultaneous equation models

# 2 Objective

The authors propose a robust instrumental variable (RIV) estimator based
on the robustification of the OIV’s closed form formula. The authors
replace the sample moments in an OIV with a robust multivariate location
and scatter S-estimators. Under certain regularity conditions, the
advantage of using S-estimators is that they are consistent,
asymptotically normal, affine equivariant, positive definite, have
bounded influence function and can achieve maximal breakdown point
regardless of dimension.

# 3 Methods

The authors robusity OIV by replacing the sample estimators (*μ̂*, *Σ̂*)
with robust multivariate location and scatter S-estimator (*m*, *S*).

*α̂*<sub>*R**I**V*</sub> = *m*<sub>*Y*</sub> − *m*′<sub>*x*</sub>*β̂*<sub>*R**I**V*</sub>

*β̂*<sub>*R**I**V*</sub> = \[*S*<sub>*x**z*</sub>*S*<sub>*z**z*</sub><sup> − 1</sup>*S*<sub>*z**x*</sub>\]<sup> − 1</sup>\[*S*<sub>*x**z*</sub>*S*<sub>*z**z*</sub><sup> − 1</sup>*S*<sub>*z**x*</sub>\]

# 4 Results

# 5 Questions
