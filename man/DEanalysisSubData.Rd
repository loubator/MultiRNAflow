% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEanalysisSubData.R
\name{DEanalysisSubData}
\alias{DEanalysisSubData}
\title{Sub data of a data.frame}
\usage{
DEanalysisSubData(
  Data,
  Res.DE.analysis,
  ColumnsCriteria = 1,
  Set.Operation = "union",
  Save.SubData = FALSE
)
}
\arguments{
\item{Data}{Data.frame with \eqn{N_g} rows and (\eqn{N_s+k}) columns,
where \eqn{N_g} is the number of genes,
\eqn{N_s} is the number of samples and
\eqn{k=1} if a column is used to specify gene names, or \eqn{k=0} otherwise.
If \eqn{k=1}, the position of the column containing gene names is given by
\code{Column.gene}.
The data.frame contains numeric values giving gene expressions of each gene
in each sample. Gene expressions can be raw counts or normalized raw counts.
Data and \code{Res.DE.analysis} must have the same number of rows.}

\item{Res.DE.analysis}{A list containing a data.frame or a data.frame or
a binary vector.
If it is a list, it must be the outputs from \code{\link[=DEanalysisGlobal]{DEanalysisGlobal()}}
(see \code{Examples}).
If it is a data.frame, it must contains at least one binary column
(filled with 0 and 1).}

\item{ColumnsCriteria}{A vector of integers where each integer indicates
a column of \code{Res.DE.analysis}.
These columns should either contain only binary values, or may contain other
numerical value, in which case extracted rows from \code{Data} will be
those with >0 values (see \code{Details}).}

\item{Set.Operation}{A character.
The user must choose between "union" (default), "intersect", "setdiff"
(see \code{Details}).}

\item{Save.SubData}{\code{TRUE} or \code{FALSE} or a Character.
\code{NULL} as default.
If \code{TRUE}, two csv files (see \code{Value}) will be saved in the folder
"2_SupervisedAnalysis_\code{Name.folder.DE}" (see \code{\link[=DEanalysisGlobal]{DEanalysisGlobal()}}).}
}
\value{
The function returns
\itemize{
\item A sub data.frame of \code{Data} containing only the rows specified by
\code{ColumnsCriteria} and \code{Set.Operation}.
\item A sub data.frame of \code{Data} containing only the rows specified by
\code{ColumnsCriteria} and \code{Set.Operation}.
\item The rows specified by \code{ColumnsCriteria} and \code{Set.Operation}.
}
}
\description{
From an initial data.frame, the function extracts a sub
data.frame containing predefined rows.
}
\details{
If \code{Res.DE.analysis} is a data.frame() or the output from
\code{\link[=DEanalysisGlobal]{DEanalysisGlobal()}}, then
\itemize{
\item If \code{Set.Operation="union"} then the rows extracted from \code{Data}
are those such that the sum of the selected columns by \code{ColumnsCriteria}
in \code{Res.DE.analysis} is >0.
For example, if \code{Res.DE.analysis} is the outputs from
\code{\link[=DEanalysisGlobal]{DEanalysisGlobal()}}, the rows extracted from \code{Data} will be those DE
at least at one time ti (except the reference time t0).
\item If \code{Set.Operation="intersect"} then the rows extracted from
\code{Data} are those such that the product of the selected columns
by \code{ColumnsCriteria} in \code{Res.DE.analysis} is >0.
For example, if \code{Res.DE.analysis} is the outputs from
\code{\link[=DEanalysisGlobal]{DEanalysisGlobal()}}, the rows extracted from \code{Data} will be those DE
at all time ti (except the reference time t0).
\item If \code{Set.Operation="setdiff"} then the rows extracted from \code{Data}
are those such that only one element of the selected columns by
\code{ColumnsCriteria} in \code{Res.DE.analysis} is >0.
For example, if \code{Res.DE.analysis} is the outputs from
\code{\link[=DEanalysisGlobal]{DEanalysisGlobal()}}, the rows extracted from \code{Data} will be those DE
at only one time ti (except the reference time t0).
}

If \code{Res.DE.analysis} is a vector, \code{ColumnsCriteria} and
\code{Set.Operation} are not used.
}
\examples{
Data.EX<-data.frame(matrix(sample(x=1:150, size=20), nrow=5))
colnames(Data.EX)<-paste("Gene",1:4,sep="")
#--------------------------------------------------------------------------#
# Exemple of output from the function DEanalysisGlobal()
res.DEanalysisGlobal.Ex<-list(data.frame(Gene=paste("Gene",1:5,sep="."),
                                         DE1=c(0,1,0,0,1),
                                         DE2=c(0,1,0,1,0)))
names(res.DEanalysisGlobal.Ex)<-c("DE.results")
#--------------------------------------------------------------------------#
res.SubDE<-DEanalysisSubData(Data=Data.EX,
                             Res.DE.analysis=res.DEanalysisGlobal.Ex,
                             ColumnsCriteria=c(2,3),
                             Set.Operation="union",
                             Save.SubData=FALSE)
}
