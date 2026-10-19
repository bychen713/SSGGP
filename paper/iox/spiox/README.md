# Inside-out cross-covariance for multivariate GPs

This package implements Bayesian hierarchical models for fitting GPs with IOX cross-covariance. 


## Install

In R, `devtools::install_github("mkln/spiox")` to compile from source.

## Inside-out cross-covariance for spatial multivariate data


As the spatial features of multivariate data are increasingly central in researchers' applied problems, 
there is a growing demand for novel spatially-aware methods that are flexible, easily interpretable, 
and scalable to large data. We develop inside-out cross-covariance (IOX) models for 
multivariate spatial likelihood-based inference. IOX leads to valid cross-covariance matrix functions 
which we interpret as inducing spatial dependence on independent replicates of a correlated random vector. 
The resulting sample cross-covariance matrices are "inside-out" relative to the ubiquitous 
linear model of coregionalization (LMC). However, unlike LMCs, our methods offer direct marginal inference, 
easy prior elicitation of covariance parameters, the ability to model outcomes with unequal smoothness, 
and flexible dimension reduction. As a covariance model for a q-variate Gaussian process, IOX leads to 
scalable models for noisy vector data as well as flexible latent models. For large n cases, IOX complements 
Vecchia approximations and related process-based methods based on sparse graphical models. 
We demonstrate superior performance of IOX on synthetic datasets as well as on colorectal cancer proteomics data.


#### Links to the article

- [doi.org/10.1080/01621459.2026.2640644](https://doi.org/10.1080/01621459.2026.2640644)
- [arXiv:2412.12407](https://arxiv.org/abs/2412.12407) 


### Citation

```
@article{,
  author  = {Michele Peruzzi},
  title   = {Inside-out cross-covariance for spatial multivariate data},
  journal = {Journal of the American Statistical Association (Forthcoming).},
  year    = {2026},
  doi     = {10.1080/01621459.2026.2640644}
}
```
