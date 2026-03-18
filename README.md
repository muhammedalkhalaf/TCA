# TCA: Transmission Channel Analysis for Stata

Stata implementation of Transmission Channel Analysis (TCA) for Structural VAR models, following Wegner, Lieb & Smeekes (2025).

## Installation

```stata
ssc install tca
```

## Description

TCA decomposes impulse response functions (IRFs) from SVAR models into contributions from distinct causal transmission channels using directed acyclic graph (DAG) path analysis.

## Features

- Channel-specific IRF decomposition
- Overlapping channel handling via inclusion-exclusion principle
- Exhaustive 3-way and 4-way decompositions
- Publication-ready output tables

## R Version

For the R implementation, see [SVARtca](https://github.com/muhammedalkhalaf/SVARtca) (also available on [CRAN](https://cran.r-project.org/package=SVARtca)).

## Reference

Wegner, E., Lieb, L., & Smeekes, S. (2025). Transmission Channel Analysis in Dynamic Models. *arXiv:2405.18987*.

## Author

Muhammad Abdullah Alkhalaf — [ORCID](https://orcid.org/0009-0002-2677-9246)

## License

GPL-3
