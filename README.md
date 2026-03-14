# TCA — Transmission Channel Analysis

[![R-CMD-check](https://github.com/muhammedalkhalaf/TCA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/muhammedalkhalaf/TCA/actions/workflows/R-CMD-check.yaml)

Implements **Transmission Channel Analysis** (TCA) for structural vector autoregressive (SVAR) models, following the methodology of:

> Wegner, E., Lieb, L., and Smeekes, S. (2025). *Transmission Channel Analysis in Dynamic Models.* [arXiv:2405.18987](https://arxiv.org/abs/2405.18987)

TCA decomposes impulse response functions (IRFs) into contributions from distinct transmission channels using a systems form representation and directed acyclic graph (DAG) path analysis.

This repository provides parallel implementations in **R** and **Stata**, both verified against the original [MATLAB toolbox](https://github.com/enweg/tca-matlab-toolbox).

---

## Features

- Systems form construction from VAR coefficient matrices (any lag order)
- Three decomposition modes:
  - **Overlapping**: each channel independently (may overlap)
  - **Exhaustive 3-way**: non-overlapping with inclusive first channel
  - **Exhaustive 4-way**: full inclusion-exclusion (var1-only, var2-only, both, direct)
- Binary additivity validation (total = through + not-through)
- Verified at machine precision against Python and MATLAB reference implementations

---

## R Package

### Installation

```r
# From GitHub
# install.packages("remotes")
remotes::install_github("muhammedalkhalaf/TCA", subdir = "R-package/TCA")
```

### Quick Example

```r
library(TCA)

K <- 4
A1 <- matrix(c(0.7,-0.1,0.05,-0.05, -0.3,0.6,0.10,-0.10,
                -0.2,0.1,0.70,0.05, -0.1,0.2,0.05,0.65), K, K, byrow=TRUE)
Sigma <- matrix(c(1,0.3,0.2,0.1, 0.3,1.5,0.25,0.15,
                   0.2,0.25,0.8,0.1, 0.1,0.15,0.1,0.6), K, K, byrow=TRUE)
Phi0 <- t(chol(Sigma))

sf <- tca_systems_form(Phi0, list(A1), h = 20)
result <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                       intermediates = c(2, 4), K = K, h = 20,
                       order = 1:K, mode = "exhaustive_4way",
                       var_names = c("IntRate","GDP","Inflation","Wages"))
plot_tca(result, target = 3)
```

### Using with the `vars` package

```r
library(vars)
data(Canada)
var_est <- VAR(Canada, p = 2, type = "const")
result <- tca_from_var(var_est, from = "e",
                        intermediates = c("prod", "rw"),
                        h = 20, mode = "exhaustive_4way")
plot_tca(result, target = "U")
```

---

## Stata Package

### Installation

```stata
* From GitHub
net install tca, from("https://raw.githubusercontent.com/muhammedalkhalaf/TCA/main/stata-package/") replace
```

### Quick Example

```stata
matrix A1 = ( 0.7, -0.1,  0.05, -0.05 \ ///
             -0.3,  0.6,  0.10, -0.10 \ ///
             -0.2,  0.1,  0.70,  0.05 \ ///
             -0.1,  0.2,  0.05,  0.65 )
matrix Sigma = ( 1, 0.3, 0.2, 0.1 \ ///
                 0.3, 1.5, 0.25, 0.15 \ ///
                 0.2, 0.25, 0.8, 0.1 \ ///
                 0.1, 0.15, 0.1, 0.6 )
matrix Phi0 = cholesky(Sigma)

tca , phi0(Phi0) ar(A1) horizon(20) from(1) ///
    intermediates(2 4) target(3) mode(exhaustive_4way) ///
    varnames(IntRate GDP Inflation Wages) graph
```

---

## Repository Structure

```
TCA/
├── R-package/
│   └── TCA/
│       ├── DESCRIPTION
│       ├── NAMESPACE
│       ├── R/
│       │   ├── systems_form.R
│       │   ├── transmission.R
│       │   ├── tca_analyze.R
│       │   ├── tca_from_var.R
│       │   └── plot_tca.R
│       ├── tests/testthat/
│       ├── vignettes/
│       └── inst/CITATION
├── stata-package/
│   ├── ado/tca.ado
│   ├── mata/tca.mata
│   ├── help/tca.sthlp
│   ├── test/tca_test.do
│   ├── stata.toc
│   └── tca.pkg
├── .github/workflows/
├── LICENSE
└── README.md
```

---

## Citation

If you use this software, please cite:

```bibtex
@article{wegner2025tca,
  title={Transmission Channel Analysis in Dynamic Models},
  author={Wegner, Emanuel and Lieb, Lenard and Smeekes, Stephan},
  year={2025},
  journal={arXiv preprint arXiv:2405.18987}
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
