# Primary Objective

This repo is to document the spatial disease risk analysis on the dataset `PBCshp`, containing PBC count and index of multiple deprivation, from the package `SDALGCP` [^fnt1], and is first published in Taylor et al. (2015)[^fnt2].

The analysis concerns primarily the implementation of the Kernel methods on aggregated LGCP model (at regional level). The objective of doing the analysis is to develop a scalable and reproducible procudure in doing modelling and inferences on continuous spatial risk (our primary concern). In doing this, we make comparison between the Markov random field based areal models (ICAR, BYM2), and discretised LGCP model in Johnson et al. (2019). 

The implementation are done in Bayesian frameworks with Stan and INLA (possibly TMB). The following numerical approximation are used:
* Basis approach (Random Fourier Features (RFF)) for modelling GP;
* Grid integration for kernel aggregation on regions;

Other assumptions are made for the practical analysis due to the limitation of the data and storage capacity, 
* Common assumption of doing spatial modelling;
* Continuity is assumed despite discrete representation on a grid. e.g. population density with certain granularity (on a grid) are assumed to be aggregated from continuous ground truth; 

Further, we consider Continuous-block cross-validation to aid the model selection with varying smoothing length-scale for the base Kernel.

# Repo Structure

The functions for main construction steps can be found in `/ArealKernel/RFFfunc.R`. The analysis using models with Stan and INLA can be found in `/ArealKernel/FullAnalysis.R`. The exploratory analysis of the block-structure incured from the kernel can be found in `/ArealKernel/FindBlocks.R`. There is also ongoing implementation with spatial random effect model in `/ArealKernel/SpatialRE`, borrwoing the same construction idea for BYM model. Ongoing implementation of Continuous-block CV can be found in `/ArealKernel/CBCV`.



[^fnt1]: Johnson, O. O., Diggle, P., & Chicas, E. G. (2019). A Spatially Discrete Approximation to Log-Gaussian Cox Processes for Modelling Aggregated Disease Count Data. Retrieved from https://arxiv.org/pdf/1901.09551.pdf

[^fnt2]: Taylor, B., Davies, T., Rowlingson, B., & Diggle, P. (2015). Bayesian inference and data augmentation schemes for spatial, spatiotemporal and multivariate log-Gaussian Cox processes in R. Journal of Statistical Software, 63, 1-48.
