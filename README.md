# `SBGLM`: Sparse Bayes Generalized Linear Models

## Installation

``` r
if (!require(devtools)) {
    install.packages('devtools')
}
devtools::install_github('miguelbiron/SBGLM')
```

## Description

The purpose of this package is mainly a way to store a lot of code that I have rotting in my laptop, which I think is relevant and useful. The models will follow the Bayesian tradition of the “spike-and-slab” prior for sparsity (Mitchell and Beauchamp 1988), so do not expect to see the so called “Bayesian Lasso” here because [it doesn’t work](https://andrewgelman.com/2017/11/02/king-must-die/). The idea is that, as time permits, I will be adding more models to the package.

## Implemented models

- A Bayesian sparse linear regression model (`sblm`)
    - Similar in spirit to the spike-and-slab, but with a slightly different approach. I describe this model in [this blog entry](https://miguelbiron.github.io).
- Non-parametric Sparse Factor Analysis (Knowles and Ghahramani 2011).
    - Check out [this blog entry](https://miguelbiron.github.io/2018/12/20/sbglm-sparse-bayes-generalized-linear-models/) where you can find more information on this model.

## TODO

- SBLM
    - Should we treat an intercept differently?
- NSFA
    - Use pre-allocated matrices instead of current naive solution of varying size dynamically.
- Add more models.
