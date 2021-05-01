
<!-- README.md is generated from README.Rmd. Please edit that file -->

# skydivide

<!-- badges: start -->
<!-- badges: end -->

Skydivide is a package for doing divide-and-conquer MCMC in
phylodynamics. User first split the genomic data into disjoint subsets,
run MCMC using each subsets in BEAST separately. Then skydivide combines
all the outputs of BEAST and reconstructs the effective population size.

## Installation

<!-- You can install the released version of skydivide from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("skydivde") -->
<!-- ``` -->

You can install the development version of skydivide from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("liusf15/skydivide")
```

The example can be found in the [vignette](vignettes/vignette.Rmd)

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- library(skydivide) -->
<!-- ## basic example code -->
<!-- ``` -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
