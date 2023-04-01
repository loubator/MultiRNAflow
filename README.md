
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MultiRNAflow

<!-- badges: start -->
<!-- badges: end -->

Our R package MultiRNAflow, built from the R package DESeq2 , provides
an easy to use unified framework allowing to automatically make both
unsupervised and supervised (DE) analysis for datasets with an arbitrary
number of biological conditions and time points. Specifically, our
package realizes:

In particular, our code makes a deep downstream analysis of DE
information, e.g. identifying temporal patterns across biological
conditions and DE genes which are specific to a biological condition for
each time.

Our package takes as input a transcriptional dataset with an
experimental design that includes multiple groups of individuals and/or
multiple times, t0 (reference time) and ti (1 $\leq$ i $\leq$ n). The
dataset is a table of raw counts where lines correspond to genes and
columns correspond to samples. Each sample show the raw counts of an
individual sequencing, corresponding to a biological condition, an
individual sampling in this biological condition and a time.

## Installation

Before installing the necessary packages, you must install (or update)
the R version 4.2.1 “Funny-Looking Kid” (released on 2022/06/23) from

Then, in order to use the MultiRNAflow package, the following R packages
must be installed:

Before installing a package, for instance the package FactoMineR, the
user must check if the package is already installed with the command .
If the package has not been previously installed, the user must use the
command (for CRAN). For beginners in programming, we recommend to follow
the steps below for importing CRAN and Bioconductor packages.

For the packages which must be download from CRAN,

``` r
Cran.pck<-c("scales", "reshape2", "plyr",
            "ggplot2", "ggrepel", "ggalluvial",
            "FactoMineR", "factoextra",
            "plot3D", "plot3Drgl",
            "UpSetR", "gprofiler2")
```

the user can copy and paste the following lines of code for each package
in order to download the missing packages.

``` rinstallcranpck
Select.package.CRAN<-"FactoMineR"
#
if(!require(package=Select.package.CRAN,
            quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)){
  install.packages(pkgs=Select.package.CRAN, dependencies=TRUE)
}# if(!require(package=Cran.pck[i], quietly=TRUE, character.only=TRUE))
```

If the package is already installed (for instance here “FactoMineR”),
the previous lines of code will return nothing.

For the packages which must be download from Bioconductor,

``` rbioconductorpck
Bioconductor.pck<-c("Mfuzz", "SummarizedExperiment", "DESeq2", "ComplexHeatmap")
```

the user must first copy and paste the following lines of code in order
to install “BiocManager”

``` r
if(!require(package="BiocManager",
            quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)){
  install.packages("BiocManager")
}# if(!require(package="BiocManager", quietly=TRUE, character.only=TRUE))
```

then copy and paste the following lines of code in order to install the
version 3.16 of bioconductor (it works with R version 4.2.0)

``` rbiocconductorversion
BiocManager::install(version = "3.16")
```

and then copy and paste the following lines of code for each package in
order to download the missing packages.

``` r
Select.package.CRAN<-"DESeq2"
if(!require(package=Select.package.CRAN,
            quietly=TRUE, character.only=TRUE, warn.conflicts=FALSE)){
  BiocManager::install(pkgs=Select.package.CRAN)
}# if(!require(package=Select.package.CRAN, quietly=TRUE, character.only=TRUE))
```

If the package is already installed (for instance here “DESeq2”), the
previous lines of code will return nothing.

Once all packages have been installed, you may install and load our
package.

``` r
# BiocManager::install(pkgs="MultiRNAflow") # to install
library(MultiRNAflow)
```
