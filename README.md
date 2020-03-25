# FAIRforecast
*Automatic interpretable sales forecasting*

FAIR is a forecasting tool developed to support decision making in a retail environment. It provides multi-step-ahead sales forecasts at the category-store level, which are based on an interpretable and transparent model. These aspects make it an objective tool with which different promotional strategies (scenarios) can be compared and insights can be generated. Additionally, the FAIRforecast package provides plotting functions that generate figures showing the relative strength of the interactions between categories and important promotional variables. For more information, see the [vignette](https://rgurlek.github.io/FAIRforecast/) and [paper](http://home.ku.edu.tr/~oali/Automatic%20Interpretable%20Retail%20Forecasting%20with%20Promotional%20Scenarios.pdf).

# Installation
``` r
# Simple installation
devtools::install_github("https://github.com/rgurlek/FAIRforecast")
# Install and build the vignette
devtools::install_github("https://github.com/rgurlek/FAIRforecast", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```
