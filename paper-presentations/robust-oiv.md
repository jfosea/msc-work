Summary of A Natural Robustification of the Ordinary Instrumental
Variables Estimator
================

*[Link](https://onlinelibrary.wiley.com/doi/10.1111/biom.12043) to the
paper.*

# 1 Problem

-   Ordinary Least Squares (OLS) estimates are biased and inconsistent
    when one of the covariates is correlated with the error terms
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

## 3.1 Estimator

The OIV estimator is given by

*α̂*<sub>*O**I**V*</sub> = *μ̂*<sub>*y*</sub> − *μ̂*<sub>*x*</sub>′*β̂*<sub>*O**I**V*</sub>

*β̂*<sub>*O**I**V*</sub> = \[*Σ̂*<sub>*X**Z*</sub>*Σ̂*<sub>*Z**Z*</sub><sup> − 1</sup>*Σ̂*<sub>*Z**X*</sub>\]<sup> − 1</sup>\[*Σ̂*<sub>*X**Z*</sub>*Σ̂*<sub>*Z**Z*</sub><sup> − 1</sup>*Σ̂*<sub>*Z**X*</sub>\]

where (*μ̂*, *Σ̂*) is the sample mean and sample covariance matrix based
on a sample
{(*x*<sub>*i*</sub>, *z*<sub>*i*</sub>, *y*<sub>*i*</sub>)}<sub>*i* = 1</sub><sup>*n*</sup>
where
(*x*<sub>*i*</sub>, *z*<sub>*i*</sub>, *y*<sub>*i*</sub>) ∈ *R*<sup>*p*</sup> × *R*<sup>*q*</sup> × *R*

(For Jana: *x*<sub>*i*</sub>’s are exposures, *z*<sub>*i*</sub>’s are
instruments, and *y*<sub>*i*</sub>’s are outcomes)

The authors robustify the OIV estimator by replacing the sample
estimators (*μ̂*, *Σ̂*) with robust multivariate location and scatter
S-estimator (*m*, *S*).

*α̂*<sub>*R**I**V*</sub> = *m*<sub>*Y*</sub> − *m*′<sub>*x*</sub>*β̂*<sub>*R**I**V*</sub>

*β̂*<sub>*R**I**V*</sub> = \[*S*<sub>*x**z*</sub>*S*<sub>*z**z*</sub><sup> − 1</sup>*S*<sub>*z**x*</sub>\]<sup> − 1</sup>\[*S*<sub>*x**z*</sub>*S*<sub>*z**z*</sub><sup> − 1</sup>*S*<sub>*z**x*</sub>\]

**Task:** I need to understand what are robust multivariate location and
scatter S-estimators.

## 3.2 Properties of RIV

### (1) Equivariant and Consistent

-   **Task:** I need to understand what is *regression* and *carrier*
    equivariant.

### (2) Bounded and Asymptotically Normal Influence Function

-   **Task:** I need to understand what are influence functions

### (3) High Breakdown Point

-   **Task:** I need to learn what are breakdown points.

# 4 Results

## Simulation

The authors compare performance of RIV with that of OIV and, Two-Stage
Generalized M-estimators (2SGM) (Wagenvoort and Waldmann, 2002),
Weighted Instrumental variables (WIV) (Krasker and Welsh, 1985), and
Multivariate Tau estimator (Tau) (Maronna and Yohia, 1997). They
considered the following regression model

*y*<sub>*i*</sub> = *α* + *x*<sub>1*i*</sub>*β*<sub>1</sub> + *x*<sub>2*i*</sub>*β*<sub>2</sub> + *ϵ*<sub>*i*</sub> for *i* = 1, …, *n*

The authors generated outliers of four types: (1) in the response, *y*,
(2) in the endogenous covariate, *x*<sub>1</sub>, (3) in the exogenous
covariate, *x*<sub>2</sub>, and (4) in the IV, *z*. They use MedSE which
shows estimator’s worst perforamnce obtained across all types of
contamination. Overall, the RIV had the lowest MedSE especially when
contamination percent was 20% or higher.

## Heart Study

In addition, the performance of the RIV was compared to other estimators
when applied to the Framingham Heart Study. The authors focus
specifically on the effect of systolic blood pressure (SBP) on left
atrial dimension (LAD). The other estimators, OIV, 2SGM, WIV, and Tau
did not attain nearly the same level of statistical significance
compared to the RIV. Particularly, the OIV estimand for the SBP is
signficantly difference from the other estimators is because the OIV
gives equal weights to the all observations and so is highly sensitive
to outliers.

# 5 Conclusion

The advantage of uing RIV in datasets with outliers is that it
downweights observations which signficantly depart from the bulk of the
data and helps uncover the true relationships among the variables.

# 5 Questions

-   **Task** Understand the other robust methods to obtain IV estimators
    and list down their pros/cons.
-   **Question** What are common methods to check for validity of IV?
    Are there computational methods already in place?
