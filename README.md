
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MultiRNAflow

<!-- useful internetlink
https://shields.io/badges/static-badge
https://simpleicons.org/?q=R
-->
<!-- badges: start -->

![R logo](https://img.shields.io/badge/R-code-276DC3?style=flat&logo=R)
![GitHub R package
version](https://img.shields.io/github/r-package/v/loubator/MultiRNAflow)
[![GPLv3
License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
![Static Badge](https://img.shields.io/badge/test_coverage-90%25-green)
[<img src="https://www.bioconductor.org/images/logo/jpg/bioconductor_logo_rgb.jpg" width="200" align="right"/>](https://bioconductor.org/)
<!-- badges: end -->

## Introduction

The MultiRNAflow package is aimed at biologists and bioinformaticians
who wish to automatically analyze RNAseq datasets with multiple times
and/or multiple biological conditions. The MultiRNAflow suite gathers in
a unified framework methodological tools found in various existing
packages allowing to perform:

1.  Exploratory (unsupervised) analysis of the data.
2.  Statistical (supervised) analysis of dynamic transcriptional
    expression (DE genes), based on DESeq2 R package.
3.  Functional and GO analysis of subsets of genes automatically
    selected by the package, such as specific genes or genes with a
    given DE temporal pattern.

The package automates a commonly used workflow of analysis for studying
complex biological phenomena.

Our package takes as input a transcriptional dataset with an
experimental design that includes multiple groups of individuals and/or
multiple times, t0 (reference time) and ti (1 $\leq$ i $\leq$ n). The
dataset is a table of raw counts where lines correspond to genes and
columns correspond to samples. Each sample show the raw counts of an
individual sequencing, corresponding to a biological condition, an
individual sampling in this biological condition and a time.

## Installation

Before installing the necessary packages, you must install (or update)
the R software in a version superior or equal to 4.2.1 “Funny-Looking
Kid” (released on 2022/06/23) from .

Then, in order to use the MultiRNAflow package, the following R packages
must be installed:

-   From [CRAN](https://cran.r-project.org) repository:
    [reshape2](https://cran.r-project.org/web/packages/reshape2) (&gt;=
    1.4.4), [ggplot2](https://ggplot2.tidyverse.org) (&gt;= 3.4.0),
    [ggalluvial](https://cran.r-project.org/web/packages/ggalluvial)
    (&gt;= 0.12.3),
    [ggrepel](https://cran.r-project.org/web/packages/ggrepel) (&gt;=
    0.9.2),
    [FactoMineR](https://cran.r-project.org/web/packages/FactoMineR)
    (&gt;= 2.6),
    [factoextra](https://cran.r-project.org/web/packages/factoextra)
    (&gt;= 1.0.7),
    [plot3D](https://cran.r-project.org/web/packages/plot3D) (&gt;=
    1.4), [plot3Drgl](https://cran.r-project.org/web/packages/plot3Drgl)
    (&gt;= 1.0.3),
    [ggplotify](https://cran.r-project.org/web/packages/ggplotify)
    (&gt;= 0.1.2),
    [UpSetR](https://cran.r-project.org/web/packages/UpSetR) (&gt;=
    1.4.0),
    [gprofiler2](https://cran.r-project.org/web/packages/gprofiler2)
    (&gt;= 0.2.1)
-   From [CRAN](https://cran.r-project.org) repository and usually
    already included by default in R: graphics (&gt;= 4.2.2), grDevices
    (&gt;= 4.2.2), grid (&gt;= 4.2.2), utils (&gt;= 4.2.2), stats (&gt;=
    4.2.2).
-   From [Bioconductor](https://bioconductor.org) repository:
    [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
    (&gt;= 1.28.0),
    [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
    (&gt;= 1.38.1),
    [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
    (&gt;= 2.14.0),
    [Mfuzz](https://www.bioconductor.org/packages/release/bioc/html/Mfuzz.html)
    (&gt;= 2.58.0).

Before installing a package, for instance the package FactoMineR, the
user must check if the package is already installed with the command .
If the package has not been previously installed, the user must use the
command (packages from CRAN). For beginners in programming, we recommend
to follow the steps below for importing CRAN and Bioconductor packages.

For the packages which must be download from
[CRAN](https://cran.r-project.org),

``` r
Cran.pck <- c("reshape2", "ggplot2", "ggrepel", "ggalluvial",
              "FactoMineR", "factoextra",
              "plot3D", "plot3Drgl", "ggplotify", "UpSetR", "gprofiler2")
```

the user can copy and paste the following lines of code for each package
in order to download the missing packages.

``` r
Select.package.CRAN <- "FactoMineR"
if (!require(package=Select.package.CRAN,
             quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)) {
    install.packages(pkgs=Select.package.CRAN, dependencies=TRUE)
}
```

If the package is already installed (for instance here
[FactoMineR](https://cran.r-project.org/web/packages/FactoMineR)), the
previous lines of code will return nothing.

For the packages which must be download from
[Bioconductor](https://bioconductor.org),

``` r
Bioconductor.pck <- c("SummarizedExperiment", "S4Vectors", "DESeq2",
                      "Mfuzz", "ComplexHeatmap")
```

the user must first copy and paste the following lines of code in order
to install [BiocManager](https://www.bioconductor.org/install)

``` r
if (!require(package="BiocManager",
             quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)) {
    install.packages("BiocManager")
}# if(!require(package="BiocManager", quietly=TRUE, character.only=TRUE))
```

then copy and paste the following lines of code in order to install the
version 3.18 of [Bioconductor](https://bioconductor.org) (it works with
R version 4.2.0)

``` r
BiocManager::install(version="3.18")
```

and then copy and paste the following lines of code for each package in
order to download the missing packages.

``` r
Select.package.Bioc <- "DESeq2"
if (!require(package=Select.package.Bioc,
             quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)) {
    BiocManager::install(pkgs=Select.package.Bioc)
}
```

If the package is already installed (for instance here
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)),
the previous lines of code will return nothing.

Once all packages have been installed, the user may load our package.

``` r
## BiocManager::install(pkgs="MultiRNAflow") ## to install
library(MultiRNAflow)
```
