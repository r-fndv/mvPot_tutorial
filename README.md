## mvPot: A Case Study

## Estimating the Spatial Dependence of Extreme Temperature over the Red Sea

Classical spatial models for extremes rely on block maxima.
This approach has been limited in applications by computational considerations and by its limited relevance when modelling single extreme events.
Functional peaks-over-threshold analysis does not suffer from such limitations.
Dombry and Ribatet (2015) extend the notion of exceedance to functions and introduce the r-Pareto process whose mathematical expression is simpler, compared to its max-stable counterpart, making it more suited for spatio-temporal applications on large environmental datasets.

In this workshop, the theoretical background for functional peaks-over-threshold analysis is presented through a case study where we estimate the spatial dependence of extreme temperature anomalies over the Red Sea.
Special attention is given to a process derived from log-Gaussian random functions, which relies on dependence models from the literature on spatio-temporal statistics.
Parametric inference is done using the efficient implementation of the censored likelihood function (Wadsworth and Tawn, 2014)  available in the `R` package `mvPot`.
We also present an alternative procedure based on the gradient score (de Fondeville and Davison, 2018) tractable in much higher dimensions.

# Required Software

* Mandatory `R` packages: `sp`, `ggplot2`,`evd`, `mvPot`.
