# Primary Objective

This repo is to document the spatial disease risk analysis on several datasets, including `PBCshp`, containing PBC count and index of multiple deprivation, from the package `SDALGCP` <sup>[1](#fnt1)</sup>, and is first published in Taylor et al. (2015)<sup>[2](#fnt2)</sup>; HIV datasets from DHS surveys in Sub-Sahara countries (for access of the data, please register at [DHS](https://www.dhsprogram.com/)), and overlayed population density data is from HRSL<sup>[3](#fnt3)</sup>.  

The analysis concerns primarily the implementation of the Kernel methods on aggregated LGCP model (at regional level). The objective of doing the analysis is to develop a scalable and reproducible procudure in doing modelling and inferences on continuous spatial risk (our primary concern). In doing this, we make comparison between the Markov random field based areal models (ICAR, BYM, BYM2), and discretised LGCP model in Johnson et al. (2019)<sup>[1](#fnt1)</sup>. 

The implementation are done in Bayesian frameworks with Stan and INLA (possibly TMB). We consider two approaches for the GP based model. The first concerns the weight space view of GP, then GLM model can be directly implemented with its basis functions, which, in our implementation, is approximated by a finite set of Random Fourier Features (RFF). The second approach borrows from the BYM2 model structure, where spatial random effect component employs GP (MVN for a sample) instead of ICAR. We therefore refer this approach as BYM3. Notably, when constructing the kernel, grid integration is used for approximating kernel aggregation on regions (instead of sampling quadrature points from point process).

Other assumptions are made for the practical analysis due to the limitation of the data and storage capacity, 
* Common assumption of doing spatial modelling;
* Continuity is assumed despite discrete representation on a grid. e.g. population density with certain granularity (on a grid) are assumed to be aggregated from continuous ground truth; 
* Things to add.

With the belief of continuous latent risk surface as ground truth, simulation study is to conducted for
1. Checking models behaves as expected, recovering preset parameters and interval coverage etc.
2. Exploring effect of assumption violation in a systematic way. 

Further, we consider Continuous-block cross-validation to aid the model selection with varying smoothing length-scale for the base kernel.

# Repo Structure

The functions for construction of the kernel can be found in `/ArealKernel/RegKernelFunc.R`, and wrapper functions for directly producing the kernel can be found in `/ArealKernel/WrapperFunc.R`.  The analysis using models with Stan and INLA can be currently found in `/ArealKernel/FullAnalysis.R` and `/ArealKernel/ModelsFunc.r`, however, for finalised analysis, please see `nutUK_pop.r` and `malawi_pop.r`, which are still to be completed.
There is also ongoing implementation with spatial random effect model (BYM3 as described before) in `/ArealKernel/SpatialRE`. Ongoing implementation of Continuous-block CV can be found in `/CBCV`.

The exploratory analysis of the block-structure incured from the kernel can be found in `/ArealKernel/FindBlocks.R`, and the produced plots can be found in `/ArealKernel/BlockPlots`. 


## Task List
- [ ] CBCV: 
    - [x] First, extending LOO strategy to include also neighbours;
    - [x] Try leaving out different percentages of the regions;
    - [x] Stratified CBCV, exploting the block structure of kernel.
    - [ ] Implement with modelling procedure.
- [ ] Spatial random effect: 
    - [ ] Inlcude the GP as an MVN in modelling the random effect;
        - [x] INLA (fitting directly with the precision matrix of aggregated kernel works well, but is it robust?)
        - [ ] Stan (Catalina breaks rstan, so sad. Will try running on virtual machine with Mojave.)
- [ ] HIV data: Malawi 
    - [x] fine grained population density of the region: data souce -- Facebook HRSL.
    - [x] test with the new data set, see if has established reproducible procedure.
- [ ] Simulation study: 
    - [ ] continuous LGCP or inhomogenous Poission?
    - [ ] Implement ICAR and BYM2 in INLA and Stan (straightforward code copying and pasting).
    - [ ] Model comparison with CBCV and Marginal Loglik, DIC.


## References
<a name="fnt1">1</a>: Johnson, O. O., Diggle, P., & Chicas, E. G. (2019). A Spatially Discrete Approximation to Log-Gaussian Cox Processes for Modelling Aggregated Disease Count Data. Retrieved from https://arxiv.org/pdf/1901.09551.pdf

<a name="fnt2">2</a>: Taylor, B., Davies, T., Rowlingson, B., & Diggle, P. (2015). Bayesian inference and data augmentation schemes for spatial, spatiotemporal and multivariate log-Gaussian Cox processes in R. Journal of Statistical Software, 63, 1-48.

<a name="fnt3">3</a>: Facebook Connectivity Lab and Center for International Earth Science Information Network - CIESIN - Columbia University. 2016. High Resolution Settlement Layer (HRSL). Source imagery for HRSL © 2016 DigitalGlobe. 
