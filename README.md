# Primary Objective

The repo is to document the spatial disease risk analysis on several datasets, including `PBCshp`, containing PBC count and indices of multiple deprivations from the package `SDALGCP` <sup>[1](#fnt1)</sup>, which is first published in Taylor et al. (2015)<sup>[2](#fnt2)</sup>; HIV datasets from DHS surveys in Sub-Saharan African countries (for access of the data, please register at [DHS](https://www.dhsprogram.com/)), and overlayed population density data from HRSL<sup>[3](#fnt3)</sup>.  

The analysis concerns primarily the implementation of the Kernel methods on aggregated LGCP model (at regional level). The objective is to develop a scalable and reproducible procedure in doing modelling and inferences on continuous spatial risk (our primary concern). In doing this, we make comparison between the Markov Random Field based areal models (ICAR, BYM, BYM2), and discretised LGCP model in Johnson et al. (2019)<sup>[1](#fnt1)</sup>. 

The implementation is done in Bayesian frameworks with Stan and INLA (possibly TMB). We consider two approaches for the GP based model. The first concerns the weight space view of GP, from which GLM model can be directly implemented with its basis functions, which, in our implementation, is approximated by a finite set of Random Fourier Features (RFF). The second approach borrows from the BYM2 model structure, where spatial random effect component employs GP (MVN for a sample) instead of ICAR. Notably, when constructing the kernel, grid integration is used for approximating kernel aggregation on regions (instead of sampling quadrature points from point process).

Other assumptions are made for the practical analysis due to the limitation of the data and storage capacity, 
* Continuity is assumed despite discrete representation on a grid. e.g. population density with certain granularity (on a grid) are assumed to be aggregated from continuous ground truth; 
* Things to add.

With the belief of continuous latent risk surface as ground truth, simulation study is to conducted for
1. Checking models behaves as expected, recovering preset parameters and interval coverage etc.
2. Exploring effect of assumption violation in a systematic way. 

Further, we consider Continuous-block cross-validation to aid model selection with varying smoothing lengthscales for the base kernel.

# Repo Structure

The functions for construction of the kernel can be found in `/ArealKernel/RegKernelFunc.R`, and wrapper functions for directly producing the kernel can be found in `/ArealKernel/WrapperFunc.R`.  The analysis using models with Stan and INLA can be currently found in `/ArealKernel/FullAnalysis.R` and `/ArealKernel/ModelsFunc.r`, however, for finalised analysis, please see `nutUK_pop.r` and `malawi_pop.r`, which are still to be completed.
There is also ongoing implementation with spatial random effect model (BYM3 as described before) in `/ArealKernel/SpatialRE`. Ongoing implementation of Continuous-block CV can be found in `/CBCV/CBCV.r`.

The exploratory analysis of the block-structure incured from the kernel can be found in `/ArealKernel/FindBlocks.R`, and the produced plots can be found in `/ArealKernel/BlockPlots`. 

## Concerns Raised

* Although direct implementation of BYM3 in INLA gives great result, the specification requires, however, the precision matrix of a PSD covariance matrix. 
    * By using basis approximation (RFF), the covariance matrix is rendered non-PSD (and not full rank), therefore no inverse existed. 
    * It may be better using a more flexible framework, i.e. Stan. However, the stochasticity nature of HMC may be less efficient for inference comparing to INLA.
* How to set prior in Stan?? Reconstruct the PC prior?
* Concern regarding CBCV: by selecting contiguous blocks, we want to minimise the dependency between the neighbouring regions. However, what about the neighbours of the left out neighbouring regions? On the other hand full Bayesian inference is expensive, any other way to tune hyperparameters?

## References
<a name="fnt1">1</a>: Johnson, O. O., Diggle, P., & Chicas, E. G. (2019). A Spatially Discrete Approximation to Log-Gaussian Cox Processes for Modelling Aggregated Disease Count Data. Retrieved from https://arxiv.org/pdf/1901.09551.pdf

<a name="fnt2">2</a>: Taylor, B., Davies, T., Rowlingson, B., & Diggle, P. (2015). Bayesian inference and data augmentation schemes for spatial, spatiotemporal and multivariate log-Gaussian Cox processes in R. Journal of Statistical Software, 63, 1-48.

<a name="fnt3">3</a>: Facebook Connectivity Lab and Center for International Earth Science Information Network - CIESIN - Columbia University. 2016. High Resolution Settlement Layer (HRSL). Source imagery for HRSL © 2016 DigitalGlobe. 
