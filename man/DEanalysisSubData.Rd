% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEanalysisSubData.R
\name{DEanalysisSubData}
\alias{DEanalysisSubData}
\title{Sub data of a data.frame}
\usage{
DEanalysisSubData(
  SEresDE,
  ColumnsCriteria = 1,
  Set.Operation = "union",
  Save.SubData = FALSE
)
}
\arguments{
\item{SEresDE}{A SummarizedExperiment class object. Output from
\code{\link[=DEanalysisGlobal]{DEanalysisGlobal()}}
(see \code{Examples}).}

\item{ColumnsCriteria}{A vector of integers where each integer indicates
a column of  \code{SummarizedExperiment::rowData(SEresDE)}.
These columns should either contain only binary values, or may contain other
numerical value, in which case extracted outputs from \code{SEresDE}
will be those with >0 values (see \code{Details}).}

\item{Set.Operation}{A character.
The user must choose between "union" (default), "intersect", "setdiff"
(see \code{Details}).}

\item{Save.SubData}{\code{TRUE} or \code{FALSE} or a Character.
\code{FALSE} as default.
If \code{TRUE}, two csv files (see \code{Value}) will be saved in the folder
"2_SupervisedAnalysis_\code{Name.folder.DE}"
(see \code{\link[=DEanalysisGlobal]{DEanalysisGlobal()}}).}
}
\value{
The function returns a SummarizeExperiment class object containing
\itemize{
\item sub data.frames of the different dataset included in \code{SEresDE}
containing only the rows specified by
\code{ColumnsCriteria} and \code{Set.Operation}.
\item the DE results saved in \code{SEresDE} of genes selected by
\code{ColumnsCriteria} and \code{Set.Operation}.
\item The genes specified by \code{ColumnsCriteria} and \code{Set.Operation}.
}
}
\description{
From the results from our function
\code{\link[=DEanalysisGlobal]{DEanalysisGlobal()}},
the function extracts from the SummarizedExperiment class outputs of
the subset of genes selected with the inputs \code{Set.Operation} and
\code{ColumnsCriteria}, and saves them in a SummarizeExperiment object.
}
\details{
We have the following three cases:
\itemize{
\item If \code{Set.Operation="union"} then the rows extracted from
the different datasets included in \code{SEresDE}
are those such that the sum of the selected columns of
\code{SummarizedExperiment::rowData(SEresDE)}
by \code{ColumnsCriteria} is >0.
For example, the selected genes can be those DE at least at t1 or t2
(versus the reference time t0).
\item If \code{Set.Operation="intersect"} then the rows extracted from
the different datasets included in \code{SEresDE}
are those such that the product of the selected columns of
\code{SummarizedExperiment::rowData(SEresDE)}
by \code{ColumnsCriteria} is >0.
For example, the selected genes can be those DE at times t1 and t2
(versus the reference time t0).
\item If \code{Set.Operation="setdiff"} then the rows extracted from
the different datasets included in \code{SEresDE}
are those such that only one element of the selected columns of
\code{SummarizedExperiment::rowData(SEresDE)}
by \code{ColumnsCriteria} is >0.
For example, the selected genes can be those DE at times t1 only and
at times t2 only (versus the reference time t0).
}
}
\examples{
## Simulation raw counts
resSIMcount <- RawCountsSimulation(Nb.Group=2, Nb.Time=1, Nb.per.GT=4,
                                   Nb.Gene=5)
## Preprocessing step
resDATAprepSE <- DATAprepSE(RawCounts=resSIMcount$Sim.dat,
                            Column.gene=1,
                            Group.position=1,
                            Time.position=NULL,
                            Individual.position=2)

##------------------------------------------------------------------------##
## Transformation of resDATAprepSE into results of DEanalysisGlobal
resultsExamples <- data.frame(Gene=paste0("Gene", seq_len(5)),
                              DE1=c(0, 1, 0, 0, 1),
                              DE2=c(0, 1, 0, 1, 0))
listPATHnameEx <- list(Path.result=NULL, Folder.result=NULL)

SummarizedExperiment::rowData(resDATAprepSE) <- resultsExamples
S4Vectors::metadata(resDATAprepSE)$DESeq2obj$pathNAME <- listPATHnameEx
S4Vectors::metadata(resDATAprepSE)$DESeq2obj$SEidentification<-"SEresultsDE"

##------------------------------------------------------------------------##
## results of DEanalysisSubData
resDEsub <- DEanalysisSubData(SEresDE=resDATAprepSE,
                              ColumnsCriteria=c(2, 3),
                              Set.Operation="union",
                              Save.SubData=FALSE)
}
