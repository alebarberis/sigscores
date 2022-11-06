
<!-- README.md is generated from README.Rmd. Please edit the README.Rmd file -->
<!-- # ```{r, include = FALSE} -->
<!-- # knitr::opts_chunk$set( -->
<!-- #   collapse = TRUE, -->
<!-- #   comment = "#>", -->
<!-- #   # fig.path = "man/figures/README-", -->
<!-- #   fig.path = "articles/figures/", -->
<!-- #   # fig.path = "vignettes/figures/", -->
<!-- #   # fig.path = "man/figures/", -->
<!-- #   out.width = "100%" -->
<!-- # ) -->
<!-- # ``` -->

# sigscores

**sigscores** is an R package providing ready-to-use functions to
compute summary scores.

The idea behind **sigscores** is to provide in a single R package
different summary scores allowing for an easy computation and comparison
of the metrics.

<!-- # ```{r echo=FALSE, fig.align='center', fig.cap='', out.width='75%'} -->
<!-- # knitr::include_graphics(path = 'articles/figures/renoir_workflow.png', error = FALSE) -->
<!-- # ``` -->

Different strategies are currently supported, such as the common mean,
median, and mode scores.

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install latest development version from GitHub (requires
[devtools](https://github.com/hadley/devtools) package):

``` r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github(
  repo = "alebarberis/sigscores", 
  dependencies = TRUE, 
  build_vignettes = FALSE
)
```

## Getting started

If you are just getting started with **sigscores** we recommend starting
with [Getting Started](articles/sigscores.html) section of the site.
<!-- `vignette("sigscores")` -->
