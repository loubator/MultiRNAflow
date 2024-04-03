# MultiRNAflow 0.99.0

* Added a `NEWS.md` file to track changes to the package.

# Changes in versions: 0.99.1 ... 0.99.6

Bug fixes

* Necessary correction by the automated single package builder of bioconductor.org.

# Changes in versions: 0.99.7

Suggestions and remarks from a bioconductor team member

* Correction of `DESCRIPTION` file
* Vignette (`Running_analysis_with_MultiRNAflow.Rmd`)
  - code-style to highlight (function, variable, package names)
  - More details about the vignette
  - `eval=FALSE` replaced by `eval=TRUE`
* More unit tests

# Changes in versions: 0.99.8

Suggestions and remarks from a bioconductor team member

* Correction of `DESCRIPTION` file
* add a new function to improve inputs in others function
* modification of R functions depending on previous inputs
* Vignette (`Running_analysis_with_MultiRNAflow.Rmd`)
  - code-style to highlight (function, variable, package names)
  - keeping includegraphics only for the introduction part
* More unit tests (now 45 unit tests in 25 files)
* add the R file `MultiRNAflow-package.R` for man folder

# Changes in versions: 0.99.9

Suggestions and remarks from a bioconductor team member

* modification of outputs in order to replace list ouputs by SE class object
* Vignette (`Running_analysis_with_MultiRNAflow.Rmd`)
  - modification because of ouputs modification 
* More unit tests (now 197 unit tests in 31 files), coverage = 60/100

# Changes in versions: 0.99.10

* Dependency packages:
  - we remove unnecessary dependencies: `ggsci`, `methods`, `plyr`,
  `RcolorBrew`, `rlang`, `scales`
  - we add the R package `ggplotify` from CRAN
* The package now contains 321 tests and the coverage is now 90%
* R functions of the packages
  - `DEanalysisTimeAndGroup()`: Correction when there are two time points
  by adding `data.frame(X)` as input of the function `apply()`.
  - `DATAplotExpression1Gene()`: When data depends only on biological
  conditions, violin plots are plotted only if the number of individual
  per condition is >50.
  - For all the main functions, all results (plots too now) are now saved
  in the `SummaryExperiment` class object and can be accessed with
  the R function metadata (`S4Vectors` package).
  - `PCAgraphics()`: the 2D PCA graph with temporal links is now realized
  with `ggplot2` and all 3D PCA graph is saved as `ggplot2` object
  with the R package `ggplotify`
  - `MFUZZclustersNumber()`: the plot is now realized with `ggplot2`
  - `PCAgraphics()`, `PCAanalysis()` and `HCPCanalysis()`:
  the input `D3.mouvment` is now called `motion3D`
  - Many functions have been split into several parts to make the code easier
  to read
* Vignette (`Running_analysis_with_MultiRNAflow.Rmd`)
  - modification because of outputs and inputs modification
* Vignette (`MultiRNAflow_vignette-knitr.Rnw`)
  - modification because of outputs and inputs modification
  - update of the vignette in order to answer comment of reviewers
  of Bioinformatics
