% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DATAplotBoxplotSamples.R
\name{DATAplotBoxplotSamples}
\alias{DATAplotBoxplotSamples}
\title{Visualization of the distribution of all gene expressions using
a boxplot for each sample.}
\usage{
DATAplotBoxplotSamples(
  ExprData,
  Column.gene,
  Group.position,
  Time.position,
  Individual.position,
  Log2.transformation = TRUE,
  Colored.By.Factors = FALSE,
  Color.Group = NULL,
  Plot.genes = FALSE,
  y.label = NULL
)
}
\arguments{
\item{ExprData}{Data.frame with \eqn{N_g} rows and (\eqn{N_{s+k}}) columns,
where \eqn{N_g} is the number of genes,
\eqn{N_s} is the number of samples and
\eqn{k=1} if a column is used to specify gene names, or \eqn{k=0} otherwise.
If \eqn{k=1}, the position of the column containing gene names is given
by \code{Column.gene}.
The data.frame contains numeric values giving gene expressions of each gene
in each sample.
Gene expressions can be raw counts or normalized raw counts.
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

\item{Log2.transformation}{\code{TRUE} or \code{FALSE}.
If \code{TRUE}, each numeric value \eqn{x} in \code{ExprData} will become
\eqn{log_2(x+1)} (see \code{Details}).}

\item{Colored.By.Factors}{\code{TRUE} or \code{FALSE}.
If \code{TRUE}, boxplots will be colored with different colors for different
time measurements (if data were collected at different time points).
Otherwise, boxplots will be colored with different colors for
different biological conditions.}

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
If \code{TRUE}, points representing gene expression
(normalized or raw counts) will be added for each sample.}

\item{y.label}{\code{NULL} or a character. \code{NULL} as default.
If \code{y.label} is a character, it will be the y label of the graph.
If \code{y.label=NULL}, the label will be either "log2(Gene expression +1)"
(if \code{Log2.transformation=TRUE}) either "Gene expression"
(if \code{Log2.transformation=FALSE}).}
}
\value{
The function returns a graph which plots the distribution of all gene
expressions using a boxplot for each sample (see \code{\link[ggplot2:geom_boxplot]{ggplot2::geom_boxplot()}}).
}
\description{
From a gene expression dataset
(raw counts or normalized raw counts),
the function plots the distribution of all gene expressions using
a boxplot for each sample.
}
\details{
The column names of \code{ExprData} must be a vector of strings
of characters containing
\itemize{
\item a string of characters (if \eqn{k=1}) which is the label of the column
containing gene names.
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
and where data were collected at time 't0'.
I this example, 'CLL' describe the type of cells
(here chronic lymphocytic leukemia) and is not used in our analysis.

In the string of characters 'CLL_P_t0_r1',
'r1' is localized after the third underscore, so \code{Individual.position=4},
'P' is localized after the first underscore, so \code{Group.position=2} and
't0' is localized after the second underscore, so \code{Time.position=3}.

The boxplot allows to visualize six summary statistics
(see \code{\link[ggplot2:geom_boxplot]{ggplot2::geom_boxplot()}}):
\itemize{
\item the median
\item two hinges: first and third quartiles denoted Q1 and Q3.
\item two whiskers: \eqn{W1:=Q1-1.5*IQR} and \eqn{W3:=Q3+1.5*IQR}
with \eqn{IQR=Q3-Q1}, the interquartile range.
\item outliers: data beyond the end of the whiskers are called "outlying" points
and are plotted in black.
}

For better visualization of the six summary statistics described above,
raw counts must be transformed using the function \eqn{log_2(x+1)}.
This transformation is automatically performed by other functions of
the package, such as \code{\link[=DATAnormalization]{DATAnormalization()}}.
\code{Log2.transformation} will be set as TRUE in \code{\link[=DATAnormalization]{DATAnormalization()}}
if \code{Normalization ="rle"}, otherwise \code{Log2.transformation=FALSE}.
}
\examples{
data(RawCounts_Antoszewski2022_MOUSEsub500)
DATAplotBoxplotSamples(ExprData=RawCounts_Antoszewski2022_MOUSEsub500,
                       Column.gene=1,
                       Group.position=1,
                       Time.position=NULL,
                       Individual.position=2,
                       Log2.transformation=TRUE,
                       Colored.By.Factors=TRUE,
                       Color.Group=NULL,
                       Plot.genes=FALSE,
                       y.label=NULL)
}
\seealso{
The \code{\link[=DATAplotBoxplotSamples]{DATAplotBoxplotSamples()}} function
\itemize{
\item is used by the following function of our package: \code{\link[=DATAnormalization]{DATAnormalization()}}.
\item calls the R functions \link[ggplot2:geom_boxplot]{ggplot2::geom_boxplot} and \link[ggplot2:geom_jitter]{ggplot2::geom_jitter}
in order to print the boxplot.
}
}
