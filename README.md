
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iteval-sims

<!-- badges: start -->
<!-- badges: end -->

The goal of the **iteval-sims** [GitHub](https://github.com/) repository
is to enable replication of the simulation study in Hoogland et
al. [DOI 10.1002/sim.10186](https://doi.org/10.1002/sim.10186). The R scripts are not meant to be used for
other purposes and are not annotated. R package
[**iteval**](https://github.com/jeroenhoogland/iteval) is available for
application of the methods in practice and contains example
illustrations and help for all functions.

This README file describes the contents of the **iteval-sims**
repository and the required steps for replication of the simulations,
tables, and figures. Note that this repository does not have the
structure of an R package and does not have to be installed.

First, `population.R` simulates the population data and stores the
required representations in `population.RData`. Second, `sim.RData`
performs the simulation study. Note that this script depends on R
package **Hmisc**, **rms**, **MatchIt** and **iteval** and sources
`support_functions.R` (which only contains helper functions specific to
the simulation).

Just running `sim.RData` for a number of simulation provides some
replications to check against the pre-computed results in the
`replicate` folder. Note that the simulations take a considerable amount
of time, and that all results are available in pre-computed form. The
`replicate.R` script provides exact replicates of the Figures and Tables
in Hoogland et al. \[yet to insert DOI\]. Results have been checked on R
version 4.2.0 and using the following package versions.

Supplementary simulations on model misspecification have also been added
in the `sup` folder, where the `sim*.R` scripts can run the simulations,
`replicate_suppl.R` again provides replication of key figures, and the
`supple_sim_*` folders contain pre-computed simulation results.

``` r
library(Hmisc)
#> 
#> Attaching package: 'Hmisc'
#> The following objects are masked from 'package:base':
#> 
#>     format.pval, units
library(rms)
library(MatchIt)
library(iteval)
library(ggplot2)
sessionInfo()
#> R version 4.2.0 (2022-04-22)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur/Monterey 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_3.5.0 iteval_0.1.0  MatchIt_4.5.5 rms_6.7-1     Hmisc_5.1-1  
#> 
#> loaded via a namespace (and not attached):
#>  [1] zoo_1.8-12         tidyselect_1.2.0   xfun_0.42          splines_4.2.0     
#>  [5] lattice_0.20-45    colorspace_2.1-0   vctrs_0.6.5        generics_0.1.3    
#>  [9] htmltools_0.5.7    yaml_2.3.8         base64enc_0.1-3    utf8_1.2.4        
#> [13] survival_3.4-0     rlang_1.1.3        pillar_1.9.0       withr_3.0.0       
#> [17] foreign_0.8-84     glue_1.7.0         multcomp_1.4-25    lifecycle_1.0.4   
#> [21] stringr_1.5.1      MatrixModels_0.5-3 munsell_0.5.0      gtable_0.3.4      
#> [25] mvtnorm_1.2-3      htmlwidgets_1.6.4  codetools_0.2-18   evaluate_0.23     
#> [29] knitr_1.45         fastmap_1.1.1      SparseM_1.81       quantreg_5.97     
#> [33] fansi_1.0.6        htmlTable_2.4.2    Rcpp_1.0.12        TH.data_1.1-2     
#> [37] scales_1.3.0       backports_1.4.1    checkmate_2.3.1    gridExtra_2.3     
#> [41] digest_0.6.34      stringi_1.8.3      polspline_1.1.23   dplyr_1.1.4       
#> [45] grid_4.2.0         cli_3.6.2          tools_4.2.0        sandwich_3.1-0    
#> [49] magrittr_2.0.3     tibble_3.2.1       Formula_1.2-5      cluster_2.1.4     
#> [53] pkgconfig_2.0.3    MASS_7.3-58.1      Matrix_1.6-5       data.table_1.15.2 
#> [57] rmarkdown_2.26     rstudioapi_0.15.0  R6_2.5.1           rpart_4.1.19      
#> [61] nnet_7.3-18        nlme_3.1-163       compiler_4.2.0
```
