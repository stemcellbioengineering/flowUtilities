README
================

# `flowUtilities`

[Report a
bug](https://github.com/stemcellbioengineering/flowUtilities/issues)

`flowUtilities` is an R package that simplifies flow cytometry data
analysis by providing convenient wrapper functions and utilities for
common tasks in the flow cytometry analysis workflow. The package builds
on popular libraries within the
[Bioconductor](https://bioconductor.org/) ecosystem to streamline data
import, gating, visualization, and statistical analysis.

## Installation

You can install `flowUtilities` from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("stemcellbioengineering/flowUtilities")
```

## Dependencies

`flowUtilities` requires the following Bioconductor packages:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("flowCore", "flowWorkspace", "flowStats", "ggcyto"))
```

These should be installed before installing `flowUtilities`.

## Documentation

For a complete workflow example, see the package vignette:

``` r
browseVignettes("flowUtilities")
```

An `.html` format on Github can be viewed by appending the vignette url
to <http://htmlpreview.github.io/>?

## Contributions

Contributors are welcome! If you have an idea for a new feature you can
[issue a pull
request](https://github.com/stemcellbioengineering/flowUtilities/pulls).

## License

This project is licensed under the GNU General Public License v3.0 - see
the LICENSE file for details.
