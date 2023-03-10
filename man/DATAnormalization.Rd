% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DATAnormalization.R
\name{DATAnormalization}
\alias{DATAnormalization}
\title{Normalization of raw counts.}
\usage{
DATAnormalization(
  RawCounts,
  Column.gene,
  Group.position,
  Time.position,
  Individual.position,
  Normalization = "vst",
  Blind.rlog.vst = FALSE,
  Plot.Boxplot = TRUE,
  Colored.By.Factors = FALSE,
  Color.Group = NULL,
  Plot.genes = FALSE,
  path.result = NULL,
  Name.folder.norm = NULL
)
}
\arguments{
\item{RawCounts}{Data.frame with \eqn{N_g} rows and (\eqn{N_{s+k}}) columns,
where \eqn{N_g} is the number of genes,
\eqn{N_s} is the number of samples and
\eqn{k=1} if a column is used to specify gene names,
or \eqn{k=0} otherwise.
If \eqn{k=1}, the position of the column containing gene names is given
by \code{Column.gene}.
The data.frame contains non negative integers giving gene expressions of
each gene in each sample.
Column names of the data.frame must describe each sample's information
(individual, biological condition and time) and have the structure described
in the section \code{Details}.}

\item{Column.gene}{Integer indicating the column where gene names are given.
Set \code{Column.gene=NULL} if there is no such column.}

\item{Group.position}{Integer indicating the position of group information
in the string of characters in each sample names (see \code{Details}).
Set \code{Group.position=NULL} if there is only one or
no biological information in the string of character in each sample name.}

\item{Time.position}{Integer indicating the position of time measurement
information in the string of characters in each sample names
(see \code{Details}).
Set \code{Time.position=NULL} if there is only one or no time measurement
information in the string of character in each sample name.}

\item{Individual.position}{Integer indicating the position of the name of
the individual (e.g patient, replicate, mouse, yeasts culture ...)
in the string of characters in each sample names (see \code{Details}).
The names of different individuals must be all different.
Furthermore, if individual names are just numbers, they will be transform in
a vector of class "character" by \code{\link[=CharacterNumbers]{CharacterNumbers()}} and
a "r" will be added to each individual name ("r" for replicate).}

\item{Normalization}{"rle", "vst", "rlog".
Each corresponds to a method of normalization proposed by \code{DESeq2}
(see \code{\link[BiocGenerics:dge]{BiocGenerics::estimateSizeFactors()}} for "rle",
\code{\link[DESeq2:rlog]{DESeq2::rlog()}} for "rlog" and \code{\link[DESeq2:vst]{DESeq2::vst()}} for "vst").}

\item{Blind.rlog.vst}{TRUE or FALSE. See input 'blind' in \code{\link[DESeq2:rlog]{DESeq2::rlog()}}.
It is recommended to set \code{Blind.rlog.vst=FALSE} for downstream analysis.}

\item{Plot.Boxplot}{\code{TRUE} or \code{FALSE}. TRUE by default.
If \code{Plot.Boxplot=TRUE}, the function \code{\link[=DATAplotBoxplotSamples]{DATAplotBoxplotSamples()}} will be
called and boxplots will be plotted. Otherwise, no boxplots will be plotted.}

\item{Colored.By.Factors}{\code{TRUE} or \code{FALSE}.
If \code{TRUE}, boxplots will be colored with different colors for different
time measurements (if data were collected at different time points).
Otherwise, boxplots will be colored with different colors for different
biological conditions.}

\item{Color.Group}{\code{NULL} or a data.frame with \eqn{N_{bc}} rows and
two columns where \eqn{N_{bc}} is the number of biological conditions.
If \code{Color.Group} is a data.frame, the first column must contain
the name of each biological condition and the second column must contain
the colors associated to each biological condition.
If \code{Color.Group=NULL}, the function will automatically attribute a color
for each biological condition.
If samples belong to different time points only,
\code{Color.Group} will not be used.}

\item{Plot.genes}{\code{TRUE} or \code{FALSE}.
If \code{TRUE}, points representing gene expressions
(normalized or raw counts) will be plotted for each sample.
Otherwise, only boxplots will be plotted.}

\item{path.result}{Character or \code{NULL}. Path to save all results.
If \code{path.result} contains a sub folder entitled
"1_Normalization_\code{Name.folder.norm}"
all results will be saved in the sub folder
"1_Normalization_\code{Name.folder.norm}".
Otherwise, a sub folder entitled "1_Normalization_\code{Name.folder.norm}"
will be created in \code{path.result} and all results will be saved in
"1_Normalization_\code{Name.folder.norm}".
If NULL, the results will not be saved in a folder. NULL as default.}

\item{Name.folder.norm}{Character or \code{NULL}.
If \code{Name.folder.norm} is a character,
the folder name which will contain the results will be
"1_Normalization_\code{Name.folder.norm}".
Otherwise, the folder name will be "1_Normalization".}
}
\value{
The function returns a normalized count data.frame and
plots a boxplot (if \code{Plot.Boxplot=TRUE}).
}
\description{
From raw counts, this function realizes one of
the three methods of normalization of the package \code{DESeq2}:
\itemize{
\item Relative Log Expression (rle) transformation
(see \code{\link[BiocGenerics:dge]{BiocGenerics::estimateSizeFactors()}})
\item Regularized Log (rlog) transformation (see \code{\link[DESeq2:rlog]{DESeq2::rlog()}})
\item Variance Stabilizing Transformation (vst) transformation
(see \code{\link[DESeq2:vst]{DESeq2::vst()}})
}
}
\details{
The column names of \code{ExprData} must be a vector of strings
of characters containing
\itemize{
\item a string of characters (if \eqn{k=1}) which is the label of
the column containing gene names.
\item \eqn{N_s} sample names which must be strings of characters containing
at least: the name of the individual (e.g patient, mouse, yeasts culture),
its biological condition (if there is at least two) and
the time where data have been collected if there is at least two;
(must be either 't0', 'T0' or '0' for time 0,
't1', 'T1' or '1' for time 1, ...).
}

All these sample information must be separated by underscores in
the sample name. For instance 'CLL_P_t0_r1',
corresponds to the patient 'r1' belonging to the biological condition 'P'
and where data were collected at time 't0'. I this example, 'CLL' describe
the type of cells (here chronic lymphocytic leukemia) and
is not used in our analysis.

In the string of characters 'CLL_P_t0_r1',
'r1' is localized after the third underscore, so \code{Individual.position=4},
'P' is localized after the first underscore, so \code{Group.position=2} and
't0' is localized after the second underscore, so \code{Time.position=3}.
}
\examples{
data(RawCounts_Antoszewski2022_MOUSEsub500)
res.Norm<-DATAnormalization(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
                            Column.gene=1,
                            Group.position=1,
                            Time.position=NULL,
                            Individual.position=2,
                            Normalization="rle",
                            Blind.rlog.vst=FALSE,
                            Plot.Boxplot=TRUE,
                            Colored.By.Factors=FALSE,
                            Color.Group=NULL,
                            Plot.genes=FALSE,
                            path.result=NULL,
                            Name.folder.norm=NULL)
print(res.Norm)
}
\seealso{
The \code{\link[=DATAnormalization]{DATAnormalization()}} function calls the R functions
\code{\link[BiocGenerics:dge]{BiocGenerics::estimateSizeFactors()}}, \code{\link[DESeq2:rlog]{DESeq2::rlog()}} and \code{\link[DESeq2:vst]{DESeq2::vst()}}
in order to realized the normalization.
}
