
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iteval-sims

<!-- badges: start -->
<!-- badges: end -->

The goal of the **iteval-sims** [GitHub](https://github.com/) repository
is to enable replication of the simulation study in Hoogland et al. 2022
\[yet to insert DOI\]. The R scripts are not meant to be used for other
purposes and are not annotated. R package
[**iteval**](https://github.com/jeroenhoogland/iteval) is available for
application of the methods in practice and contains example
illustrations and help for all functions.

This README file describes the contents of the **iteval-sims**
repository and the required steps for replication of the simulations,
tables, and figures. Note that this repository does not have the
structure of an R package and does not have to be installed.

First, `population.R` simulates the population data and stored the
required representations in `population.RData`. Second, `sim.RData`
performs a single run of the simulation study. Note that this script
depends on R package **Hmisc**, **rms**, and **iteval** and sources
`support_functions.R` (which only contains helper functions specific to
the simulation).

Just running `sim.RData` for a number of simulation provides replication
with random noise. For exact replication, the seeds used were 1:500
(Mersenne-Twister type). To avoid the computation time, the results are
stored in the `replicate` folder. The `replicate.R` script provides
exact replicates of the Figures and Tables in Hoogland et al. 2022 \[yet
to insert DOI\]. Results have been checked on R version 3.6.3 using
`Hmisc_4.7-0`, `rms_6.3-0`, `ggplot2_3.3.5`, and `iteval0_1.0` and R
version 4.2.0 using `Hmisc_4.7-1`, `rms_6.3-0`, `ggplot2_3.3.6`, and
`iteval0_1.0`
