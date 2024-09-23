---
title: 'LocalCop: An R package for local likelihood inference for conditional copulas'
params:
  do_calc: FALSE
tags:
  - R
  - C++
  - copula
  - local likelihood
  - covariate effect
authors:
  - name: Elif Fidan Acar
    orcid: 0000-0003-2908-7691
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Martin Lysy
    orcid: 0000-0001-9974-1121
    affiliation: 3
  - name: Alan Kuchinsky
    affiliation: 4
affiliations:
 - name: University of Guelph
   index: 1
 - name: Hospital for Sick Children
   index: 2
 - name: University of Waterloo
   index: 3
 - name: University of Manitoba
   index: 4
citation_author: Acar et. al.
date: 9 September 2024
year: 2024
bibliography: references.bib
output:
  # Note: bookdown cross-links incompatible with joss md processor
  # bookdown::pdf_book:
  #   base_format: rticles::joss_article
  #   toc: false
  #   keep_md: true
  rticles::joss_article:
    keep_md: true
header-includes:
  # to use kableExtra
  - \usepackage{booktabs}
  # - \usepackage{longtable}
  # - \usepackage{array}
  # - \usepackage{multirow}
  # - \usepackage{wrapfig}
  # - \usepackage{float}
  # - \usepackage{colortbl}
  # - \usepackage{pdflscape}
  # - \usepackage{tabu}
  # - \usepackage{threeparttable}
  # - \usepackage{threeparttablex}
  # - \usepackage[normalem]{ulem}
  # - \usepackage{makecell}
  # - \usepackage{xcolor}
csl: apa.csl
journal: JOSS
link-citations: true
---

<!-- 

~~In order to compile this file, you must run `source("joss_article2.R")` first, due to a bug in `rticles::joss_article()`.~~

**Update:** The bug has now been fixed on the development version of [**rticles**](https://github.com/rstudio/rticles).  In order to compile the manuscript, it now suffices to have this version installed: `remotes::install_github('rstudio/rticles', upgrade = TRUE)`.

**Update 2:** Turns out that **bookdown** cross-links (e.g., `\@ref(eq:xyz)` are not processed by `keep_md: true` and therefore incompatible with the joss md processor.  So, the bookdown extensions have been disabled altogether.

-->



<!-- latex macros -->
\newcommand{\R}{{\textsf{R}}}
\newcommand{\cpp}{{\textsf{C++}}}

# Summary

Conditional copulas models allow the dependence structure between multiple response variables to be modelled as a function of covariates.  **LocalCop** [@localcop] is an \R/\cpp package for computationally efficient semiparametric conditional copula modelling using a local likelihood inference framework developed in @ACY2011, @ACY2013 and @ACL2019.

# Statement of Need

There are well-developed \R packages such as **copula** [@copula1; @copula2; @copula3; @copula4] and **VineCopula** [@vinecopula] for fitting copulas in various multivariate data settings. However, these software focus exclusively on unconditional dependence modelling and do not accommodate covariate information. 

Aside from **LocalCop**, \R packages for fitting conditional copulas are **gamCopula** [@gamcopula] and **CondCopulas** [@condcopulas]. **gamCopula** estimates the covariate-dependent copula parameter using spline smoothing.  While this typically has lower variance than the local likelihood estimate provided by **LocalCop**, it also tends to have lower accuracy [@ACL2019].  **CondCopulas** estimates the copula parameter using a semi-parametric maximum-likelihood method based on a kernel-weighted conditional concordance metric.  **LocalCop** also uses kernel weighting, but it uses the full likelihood information of a given copula family rather than just that contained in the concordance metric, and is therefore more statistically efficient.

Local likelihood methods typically involve solving a large number of low-dimensional optimization problems and thus can be computationally intensive.  To address this issue, **LocalCop** implements the local likelihood function in \cpp, using the \R/\cpp package **TMB** [@kristensen.etal16] to efficiently obtain the associated score function using automatic differentiation.  Thus, **LocalCop** is able to solve each optimization problem very quickly using gradient-based algorithms.  It also provides a means of easily parallelizing the optimization across multiple cores, rendering **LocalCop** competitive in terms of speed with other available software for conditional copula estimation.

# Background

For any bivariate response vector $(Y_1, Y_2)$, the conditional joint distribution given a covariate $X$ is given by 
\begin{equation}
F_X(y_1, y_2 \mid x) = C_X (F_{1\mid X} (y_1 \mid x),F_{2\mid X} (y_2 \mid x) \mid x ),
\label{eq:fullmodel}
\end{equation}
where $F_{1\mid X}(y_1 \mid x)$ and $F_{2\mid X}(y_2 \mid x)$ are the conditional marginal distributions of $Y_1$ and $Y_2$ given $X$, and $C_X(u, v \mid x)$ is a conditional copula function.  That is, for given $X = x$, the function $C_X(u, v \mid x)$ is a bivariate CDF with uniform margins. 

The focus of **LocalCop** is on estimating the conditional copula function, which is modelled semi-parametrically as 
\begin{equation}
C_X(u, v \mid x) = \mathcal{C}(u, v\mid \theta(x), \nu),
\label{eq:copmod}
\end{equation}
where $\mathcal{C}(u, v \mid \theta, \nu)$ is a parametric copula family, the copula dependence parameter $\theta \in \Theta$ is an arbitrary function of $X$, and $\nu \in \Upsilon$ is an additional copula parameter present in some models. Since most parametric copula families have a restricted range $\Theta \subsetneq \mathbb{R}$, we describe the data generating model (DGM) in terms of the calibration function $\eta(x)$, such that
\begin{equation}
\theta(x) = g^{-1}(\eta(x)),
\end{equation}
where $g^{-1}: \mathbb{R} \to \Theta$ an inverse-link function which ensures that the copula parameter has the correct range. The choice of $g^{-1}(\eta)$ is not unique and depends on the copula family.

Local likelihood estimation of the conditional copula parameter $\theta(x)$ uses Taylor expansions to approximate the calibration function $\eta(x)$ at an observed covariate value $X = x$ near a fixed point $X = x_0$, i.e.,
$$
\eta(x)\approx \eta(x_0) + \eta^{(1)}(x_0) (x - x_0) + \ldots + \frac{\eta^{(p)}(x_0)}{p!} (x - x_0)^{p}.
$$
One then estimates $\beta_k = \eta^{(k)}(x_0)/k!$ for $k = 0,\ldots,p$ using a kernel-weighted local likelihood function
\begin{equation}
\ell(\boldsymbol{\beta}) = \sum_{i=1}^n \log\left\{ c\left(u_i, v_i \mid g^{-1}( \boldsymbol{x}_{i}^T \boldsymbol{\beta}), \nu \right)\right\} K_h\left(\frac{x_i-x_0}{h}\right),
\label{eq:locallik}
\end{equation}
where $(u_i, v_i, x_i)$ is the data for observation $i$, $\boldsymbol{x}_i = (1, x_i - x_0, (x_i - x_0)^2, \ldots, (x_i - x_0)^p)$, $\boldsymbol{\beta}= (\beta_0, \beta_1, \ldots, \beta_p)$, and $K_h(z)$ is a kernel function with bandwidth parameter $h > 0$.  Having maximized $\ell(\boldsymbol{\beta})$ in \autoref{eq:locallik}, one estimates $\eta(x_0)$ by $\hat \eta(x_0) = \hat \beta_0$.  Usually, a linear fit with $p=1$ suffices to obtain a good estimate in practice.

# Usage

**LocalCop** is available on [CRAN](https://CRAN.R-project.org/package=LocalCop) and [GitHub](https://github.com/mlysy/LocalCop).  The two main package functions are:

- `CondiCopLocFit()`: For estimating the calibration function at a sequence of values $\boldsymbol{x}_0 = (x_{01}, \ldots, x_{0m})$.

- `CondiCopSelect()`: For selecting a copula family and bandwidth parameter using leave-one-out cross-validation (LOO-CV) with subsampling as described in @ACL2019.

In the following example, we illustrate the model selection/tuning and fitting steps for data generated from a Clayton copula with conditional Kendall $\tau$ displayed in \autoref{fig:copcomp}.  The CV metric for each combination of family and bandwidth are displayed in \autoref{fig:select1-plot}.


``` r
library(LocalCop)   # local likelihood estimation
library(VineCopula) # simulate copula data
```


``` r
set.seed(2024)

# simulation setting
family <- 3                    # Clayton Copula
n_obs <- 300                   # number of observations
eta_fun <- function(x) {       # calibration function
  sin(5*pi*x) + cos(8*pi*x^2)
}
```

``` r
# simulate covariate values
x <- sort(runif(n_obs))

# simulate response data
eta_true <- eta_fun(x)                     # calibration parameter eta(x)
par_true <- BiCopEta2Par(family = family,  # copula parameter theta(x)
                         eta = eta_true)
udata <- VineCopula::BiCopSim(n_obs, family = family, par = par_true)
```

``` r
# model selection and tuning
bandset <- c(.02, .05, .1, .2) # set of bandwidth parameters
famset <- c(1, 2, 3, 4, 5)     # set of copula families
kernel <- KernGaus             # kernel function
degree <- 1                    # degree of local polynomial
n_loo <- 100                   # number of LOO-CV observations
                               # (can be much smaller than n_obs)
```

``` r
# calculate cv for each combination of family and bandwidth
cvselect <- CondiCopSelect(u1= udata[,1], u2 = udata[,2],
                           x = x, xind = n_loo,
                           kernel = kernel, degree = degree,
                           family = famset, band = bandset)
```


\begin{figure}
\includegraphics[width=1\linewidth]{paper_files/figure-latex/select1-plot-1} \caption{Cross-validation metric for each combination of family and bandwidth.}\label{fig:select1-plot}
\end{figure}


``` r
# extract the selected family and bandwidth from cvselect
cv_res <- cvselect$cv
i_opt <- which.max(cv_res$cv)
fam_opt <- cv_res[i_opt,]$family
band_opt <- cv_res[i_opt,]$band

# calculate eta(x) on a grid of values
x0 <- seq(0, 1, by = 0.01)
copfit <- CondiCopLocFit(u1 = udata[,1], u2 = udata[,2],
                         x = x, x0 = x0,
                         kernel = kernel, degree = degree,
                         family = fam_opt, band = band_opt)
# convert eta to Kendall tau
tau_loc <- BiCopEta2Tau(copfit$eta, family= fam_opt)
```





``` r
# simulate covariate values
x <- sort(runif(n_obs))

# simulate response data
eta_true <- eta_fun(x)                     # calibration parameter eta(x)
par_true <- BiCopEta2Par(family = family,  # copula parameter theta(x)
                         eta = eta_true)
udata <- VineCopula::BiCopSim(n_obs, family = family, par = par_true)
```

``` r
# model selection and tuning
bandset <- c(.02, .05, .1, .2) # set of bandwidth parameters
famset <- c(1, 2, 3, 4, 5)     # set of copula families
kernel <- KernGaus             # kernel function
degree <- 1                    # degree of local polynomial
n_loo <- 100                   # number of LOO-CV observations
                               # (can be much smaller than n_obs)
```

``` r
# calculate cv for each combination of family and bandwidth
cvselect <- CondiCopSelect(u1= udata[,1], u2 = udata[,2],
                           x = x, xind = n_loo,
                           kernel = kernel, degree = degree,
                           family = famset, band = bandset)
```



``` r
# extract the selected family and bandwidth from cvselect
cv_res <- cvselect$cv
i_opt <- which.max(cv_res$cv)
fam_opt <- cv_res[i_opt,]$family
band_opt <- cv_res[i_opt,]$band

# calculate eta(x) on a grid of values
x0 <- seq(0, 1, by = 0.01)
copfit <- CondiCopLocFit(u1 = udata[,1], u2 = udata[,2],
                         x = x, x0 = x0,
                         kernel = kernel, degree = degree,
                         family = fam_opt, band = band_opt)
# convert eta to Kendall tau
tau_loc <- BiCopEta2Tau(copfit$eta, family= fam_opt)
```



\begin{figure}
\includegraphics[width=1\linewidth]{paper_files/figure-latex/copcomp-1} \caption{True vs estimated conditional Kendall $\tau$ using various methods.}\label{fig:copcomp}
\end{figure}

In \autoref{fig:copcomp}, we compare the true conditional Kendall $\tau$ to estimates using each of the three conditional copula fitting packages **LocalCop**, **gamCopula**, and **CondCopulas**, for sample sizes $n = 300$ and $n = 1000$.  In **gamCopula**, selection of the copula family smoothing splines is done using the generalized CV framework provided by the \R package **mgcv** [@wood17].  In **CondCopulas**, selection of the bandwidth parameter is done using LOO-CV.  In this particular example, the sample size of $n = 300$ is not large enough for **gamCopula** to pick a sufficiently flexible spline basis, and **CondCopulas** picks a large bandwidth which oversmooths the data.  For the larger sample size $n = 1000$, the three methods exhibit similar accuracy.

# Acknowledgements

We acknowledge funding support from the Natural Sciences and Engineering Research Council of Canada Discovery Grants RGPIN-2020-06753 (Acar) and RGPIN-2020-04364 (Lysy).

# References
