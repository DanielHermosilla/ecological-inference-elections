# *Methods for "A Fast Ecological Inference Algorithm for the RxC case".*

<!-- badges: start -->

[![R-CMD-check-Ubuntu](https://img.shields.io/github/actions/workflow/status/DanielHermosilla/ecological-inference-elections/R-CMD-check-ubuntu.yaml?branch=main&label=Ubuntu)](https://github.com/DanielHermosilla/ecological-inference-elections/actions)
[![R-CMD-check-Windows](https://img.shields.io/github/actions/workflow/status/DanielHermosilla/ecological-inference-elections/R-CMD-check-windows.yaml?branch=main&label=Windows)](https://github.com/DanielHermosilla/ecological-inference-elections/actions)
[![R-CMD-check-macOS](https://img.shields.io/github/actions/workflow/status/DanielHermosilla/ecological-inference-elections/R-CMD-check-macos.yaml?branch=main&label=macOS)](https://github.com/DanielHermosilla/ecological-inference-elections/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/fastei)](https://CRAN.R-project.org/package=fastei)

<!-- badges: end -->

The following library includes a method (`run_em`) to solve the R×C Ecological Inference problem for the non-parametric case by using the EM algorithm with different approximation methods for the E-Step. The standard deviation of the estimated probabilities can be computed using bootstrapping (bootstrap).

It also provides a function that generates synthetic election data (`simulate_election`) and a function that imports real election data (`chilean_election_2021`) from the Chilean first-round presidential election of 2021.

The setting in which the documentation presents the Ecological Inference problem is an election context where for a set of ballot-boxes we observe (i) the votes obtained by each candidate and (ii) the number of voters of each demographic group (for example, these can be defined by age ranges or sex). See [Thraves, C., Ubilla, P., Hermosilla, D. (2024): "A Fast Ecological Inference Algorithm for the R×C Case".](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834).

The methods to compute the conditional probabilities of the E-Step included in this package are the following:

- **Markov Chain Monte Carlo** (`mcmc`): Performs MCMC to sample vote outcomes for each ballot-box consistent with the observed data. This sample is used to estimate the conditional probability of the E-Step.

- **Multivariate Normal PDF** (`mvn_pdf`): Uses the PDF of a Multivariate Normal to approxi- mate the conditional probability.

- **Multivariate Normal CDF** (`mvn_cdf`): Uses the CDF of a Multivariate Normal to approxi- mate the conditional probability.

- **Multinomial** (`mult`): A single Multinomial is used to approximate the sum of Multinomial distributions.

- **Exact** (`exact`): Solves the E-Step exactly using the Total Probability Law, which requires enumerating an exponential number of terms.

On average, the Multinomial method is the most efficient and precise. Its precision matches the Exact method.

The documentation uses the following notation:

- `b`: number of ballot boxes.
- `g`: number of demographic groups.
- `c`: number of candidates.
- `a`: number of aggregated macro-groups.

To learn more about `fastei`, please consult the available vignettes:

```{r}
browseVignettes("fastei")
```

Aditionally, it is possible to browse the full documentation on the [library website](https://danielhermosilla.github.io/ecological-inference-elections/).
### Installation

It can either be installed from source by running `devtools::install_github("DanielHermosilla/ecological-inference-elections)` or by the CRAN repository `install.packages("fastei")`. Fortran support is required, usually, shipped by R itself. 
