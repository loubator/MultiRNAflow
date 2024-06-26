% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HCPCanalysis.R
\name{HCPCanalysis}
\alias{HCPCanalysis}
\title{Hierarchical clustering analysis with HCPC (Main function)}
\usage{
HCPCanalysis(
  SEresNorm,
  DATAnorm = TRUE,
  gene.deletion = NULL,
  sample.deletion = NULL,
  Plot.HCPC = FALSE,
  Color.Group = NULL,
  Phi = 25,
  Theta = 140,
  epsilon = 0.2,
  Cex.point = 0.7,
  Cex.label = 0.7,
  motion3D = FALSE,
  path.result = NULL,
  Name.folder.hcpc = NULL
)
}
\arguments{
\item{SEresNorm}{Results of the function
\code{\link[=DATAnormalization]{DATAnormalization()}}.}

\item{DATAnorm}{\code{TRUE} or \code{FALSE}. \code{TRUE} as default.
\code{TRUE} means the function uses the normalized data.
\code{FALSE} means the function uses the raw counts data.}

\item{gene.deletion}{\code{NULL} or a vector of characters or a vector of
integers. \code{NULL} as default.
If \code{gene.deletion} is a vector of characters, all genes with names in
\code{gene.deletion} will be deleted from the data set as input
\code{RawCounts} of our function
\code{\link[=DATAprepSE]{DATAprepSE()}}.
If \code{gene.deletion} is a vector of integers,
all the corresponding row numbers will be deleted from the data set as input
\code{RawCounts} of our function
\code{\link[=DATAprepSE]{DATAprepSE()}}.
If \code{gene.deletion=NULL} all genes will be used in the construction
of the PCA.}

\item{sample.deletion}{\code{NULL} or a vector of characters or
a vector of integers. \code{NULL} as default.
If \code{sample.deletion} is a vector of characters, all samples with names
in \code{sample.deletion} will not be used in the construction of the PCA.
If \code{sample.deletion} is a vector of integers,
all the corresponding column numbers will not be used in the construction
of the PCA from the data set as input \code{RawCounts} of our function
\code{\link[=DATAprepSE]{DATAprepSE()}}.
If \code{sample.deletion=NULL} all samples will be used
in the construction of the PCA.}

\item{Plot.HCPC}{\code{TRUE} or \code{FALSE}. \code{FALSE} as default.
If \code{TRUE}, all graphs will be plotted.
Otherwise no graph will be plotted.}

\item{Color.Group}{\code{NULL} or a data.frame with \eqn{N_{bc}} rows and
two columns where \eqn{N_{bc}} is the number of biological conditions.
If \code{Color.Group} is a data.frame, the first column must contain
the name of each biological condition and the second column must contain
the colors associated to each biological condition.
If \code{Color.Group=NULL}, the function will automatically attribute
a color for each biological condition.
If samples belong to different time points only,
\code{Color.Group} will not be used.}

\item{Phi}{Angle defining the colatitude direction for the 3D PCA plot
(see \code{Details} in
\code{\link[graphics:persp]{graphics::persp()}}).}

\item{Theta}{Angle defining the azimuthal direction for the 3D PCA plot
(see \code{Details} in
\code{\link[graphics:persp]{graphics::persp()}}).}

\item{epsilon}{Non negative numeric value giving the length between points
and their labels in all PCA plots which are not automatically plotted
by \code{\link[FactoMineR:PCA]{FactoMineR::PCA()}}.}

\item{Cex.point}{Non negative numeric value giving the size of points
in all PCA plots which are not automatically plotted by
\code{\link[FactoMineR:PCA]{FactoMineR::PCA()}}.}

\item{Cex.label}{Non negative numeric value giving the size of the labels
associated to each point of the all PCA graphs which are not automatically
plotted by \code{\link[FactoMineR:PCA]{FactoMineR::PCA()}}.}

\item{motion3D}{\code{TRUE} or \code{FALSE}.
If \code{TRUE}, the 3D PCA plots will also be plotted in a rgl window
(see \code{\link[plot3Drgl:plotrgl]{plot3Drgl::plotrgl()}})
allowing to interactively rotate and zoom.}

\item{path.result}{Character or \code{NULL}. Path to save all results.
If \code{path.result} contains a sub folder entitled
"1_UnsupervisedAnalysis_\code{Name.folder.hcpc}" and a sub sub folder,
"1-3_HCPCanalysis_\code{Name.folder.hcpc}"
all results will be saved in the sub folder
"1_UnsupervisedAnalysis_\code{Name.folder.hcpc}/
1-3_HCPCanalysis_\code{Name.folder.hcpc}".
Otherwise, a sub folder entitled
"1_UnsupervisedAnalysis_\code{Name.folder.hcpc}" and/or a sub sub folder
"1-3_HCPCanalysis_\code{Name.folder.hcpc}"
will be created in \code{path.result} and all results will be saved in
"1_UnsupervisedAnalysis_\code{Name.folder.hcpc}/
1-3_HCPCanalysis_\code{Name.folder.hcpc}".
If \code{NULL}, the results will not be saved in a folder.
\code{NULL} as default.}

\item{Name.folder.hcpc}{Character or \code{NULL}.
If \code{Name.folder.hcpc} is a character, the folder and sub folder names
which will contain the PCA graphs will respectively be
"1_UnsupervisedAnalysis_\code{Name.folder.hcpc}"
and "1-3_HCPCanalysis_\code{Name.folder.hcpc}".
Otherwise, the folder and sub folder names will respectively be
"1_UnsupervisedAnalysis" and "1-3_HCPCanalysis".}
}
\value{
The function returns the same SummarizedExperiment class object
\code{SEresNorm} with the outputs from the function
\code{\link[FactoMineR:HCPC]{FactoMineR::HCPC()}},
(saved in the metadata \code{Results[[1]][[3]]} of \code{SEresNorm})
\itemize{
\item a dendrogram (also called hierarchical tree) using the function
\code{\link[factoextra:fviz_dend]{factoextra::fviz_dend()}}
\item one 2D PCA and two 3D PCA produced by the function
\code{\link[=PCAgraphics]{PCAgraphics()}}
where samples are colored with different colors for different clusters.
The two 3D PCA graphs are identical but one of them will be opened
in a rgl window
(see \code{\link[plot3Drgl:plotrgl]{plot3Drgl::plotrgl()}})
allowing to interactively rotate and zoom.
The interactive 3D graph will be plotted only if \code{motion3D=TRUE}.
\item A graph indicating for each sample, its cluster and
the time and/or biological condition associated to the sample.
\item the outputs of
\code{\link[FactoMineR:HCPC]{FactoMineR::HCPC()}}.
}
}
\description{
The functions performs a hierarchical clustering on results
from a factor analysis with the R function
\code{\link[FactoMineR:HCPC]{FactoMineR::HCPC()}}.
}
\details{
All results are built from the results of our function
\code{\link[=DATAnormalization]{DATAnormalization()}}.

The number of clusters is automatically selected by
\code{\link[FactoMineR:HCPC]{FactoMineR::HCPC()}}
and is described in the section \code{Details} of
\code{\link[FactoMineR:HCPC]{FactoMineR::HCPC()}}.
}
\examples{
## Simulation raw counts
resSIMcount <- RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
                                   Nb.Gene=10)
## Preprocessing step
resDATAprepSE <- DATAprepSE(RawCounts=resSIMcount$Sim.dat,
                            Column.gene=1,
                            Group.position=1,
                            Time.position=2,
                            Individual.position=3)
## Normalization
resNorm <- DATAnormalization(SEres=resDATAprepSE,
                             Normalization="rle",
                             Plot.Boxplot=FALSE,
                             Colored.By.Factors=FALSE)
##------------------------------------------------------------------------##
resHCPCanalysis <- HCPCanalysis(SEresNorm=resNorm,
                                DATAnorm=TRUE,
                                sample.deletion=NULL,
                                gene.deletion=NULL,
                                Plot.HCPC=TRUE,
                                Color.Group=NULL,
                                Phi=25, Theta=140,
                                Cex.point=1, Cex.label=0.6, epsilon=0.4,
                                motion3D=FALSE,
                                path.result=NULL,
                                Name.folder.hcpc=NULL)
}
\seealso{
The function calls the functions
\code{\link[=PCArealization]{PCArealization()}} and
\code{\link[FactoMineR:HCPC]{FactoMineR::HCPC()}}.
The function
\code{\link[FactoMineR:HCPC]{FactoMineR::HCPC()}}
will take as input the output of
\code{\link[=PCArealization]{PCArealization()}}.
}
