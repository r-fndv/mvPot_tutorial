# mvPot: A Case Study

## Estimating the Spatial Dependence of Extreme Temperature over the Red Sea

Classical spatial models for extremes rely on block maxima.
This approach has been limited in applications by computational considerations and by its limited relevance when modelling single extreme events.
Functional peaks-over-threshold analysis does not suffer from such limitations.
Dombry and Ribatet (2015) extend the notion of exceedance to functions and introduce the *r*-Pareto process whose mathematical expression is simpler, compared to its max-stable counterpart, making it more suited for spatio-temporal applications on large environmental datasets.

In this workshop, the theoretical background for functional peaks-over-threshold analysis is presented through a case study where we estimate the spatial dependence of extreme temperature anomalies over the Red Sea.
Special attention is given to a process derived from log-Gaussian random functions, which relies on dependence models from the literature on spatio-temporal statistics.
Parametric inference is done using the efficient implementation of the censored likelihood function (Wadsworth and Tawn, 2014)  available in the `R` package `mvPot`.
We also present an alternative procedure based on the gradient score (de Fondeville and Davison, 2018) tractable in much higher dimensions.

The Red Sea data set (Donlon et al., 2012) was kindly provided by GHRSST, Met Office and CMEMS.

## Required Software

* Mandatory `R` packages: `sp`, `ggplot2`,`evd`, `mvPot`.
* Suggested `R` packages: `ggmap`. To install the last version use:
```devtools::install_github("dkahle/ggmap")```
To download background maps and use  `ggmap`, you will need a google API key. To obtain one, you can find instructions here:
```
\href{https://developers.google.com/maps/documentation/geocoding/get-api-key}{https://developers.google.com/maps/documentation/geocoding/get-api-key} 
```
Once you get an API key, copy/paste it in the file `Code/google_api_key.txt`.
* For easier presentation, the workshop is given with Jupyter Notebook. To get Jupyter, you can simply install Anaconda; instructions for linux can be found here:
```
\href{http://docs.anaconda.com/anaconda/install/linux/}{http://docs.anaconda.com/anaconda/install/linux/} 
```
You can also find instructions for MaxOS and Windows on the same website. By default, Jupyter Notebook is set up for `Python`. To add the `R` kernel, follow the instructions here:
```
\href{https://docs.anaconda.com/anaconda/navigator/tutorials/r-lang/}{https://docs.anaconda.com/anaconda/navigator/tutorials/r-lang/} 
```
Once this is done, you can start Jupyter by just typing `jupyter notebook` in a console (preferentially from the folder containing the notebook). This will open a page on your web browser with Jupyter interface.

## References

C.J. Donlon, M. Martin,J.D. Stark, J. Roberts-Jones, E. Fiedler, W. Wimmer. The operational sea surface temperature and sea ice analysis (OSTIA) system. *Remote Sensing Environ.*, 116 (2012), pp. 140-158 http://dx.doi.org/10.1016/j.rse.2010.10.017".

de Fondeville, R. and Davison, A. C. (2018). High-dimensional Peaks-over- threshold Inference. *Biometrika*, 105(3):575–592.

Dombry, C. and Ribatet, M. (2015). Functional Regular Variations, Pareto Processes and Peaks Over Thresholds. *Statistics and Its Interface*, 8(1):9–17.

Wadsworth, J. L. and Tawn, J. A. (2014). Efficient Inference for Spatial Extreme Value Processes Associated to Log-Gaussian Random Functions. *Biometrika*, 101(1):1–15.
