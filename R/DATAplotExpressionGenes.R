#' @title Plot expression of a subset of genes.
#'
#' @description The function allows to plot gene expression profiles
#' according to time and/or biological conditions.
#'
#' @details The column names of \code{ExprData} must be a vector of strings
#' of characters containing
#' * a string of characters (if \eqn{k=1}) which is the label of the column
#' containing gene names.
#' * \eqn{N_s} sample names which must be strings of characters containing
#' at least : the name of the individual (e.g patient, mouse, yeasts culture),
#' its biological condition (if there is at least two) and
#' the time where data have been collected if there is at least two;
#' (must be either 't0', 'T0' or '0' for time 0,
#' 't1', 'T1' or '1' for time 1, ...).
#'
#' All these sample information must be separated by underscores in
#' the sample name. For instance 'CLL_P_t0_r1',
#' corresponds to the patient 'r1' belonging to the biological condition 'P'
#' and where data were collected at time 't0'.
#' I this example, 'CLL' describe the type of cells
#' (here chronic lymphocytic leukemia) and is not used in our analysis.
#'
#' In the string of characters 'CLL_P_t0_r1',
#' 'r1' is localized after the third underscore, so \code{Individual.position=4},
#' 'P' is localized after the first underscore, so \code{Group.position=2} and
#' 't0' is localized after the second underscore, so \code{Time.position=3}.
#'
#' @param ExprData Data.frame with \eqn{N_g} rows and (\eqn{N_{s+k}}) columns,
#' where \eqn{N_g} is the number of genes,
#' \eqn{N_s} is the number of samples and
#' \eqn{k=1} if a column is used to specify gene names, or \eqn{k=0} otherwise.
#' If \eqn{k=1}, the position of the column containing gene names is given
#' by \code{Column.gene}.
#' The data.frame contains numeric values giving gene expressions of each gene
#' in each sample.
#' Gene expressions can be raw counts or normalized raw counts.
#' Column names of the data.frame must describe each sample's information
#' (individual, biological condition and time) and have the structure described
#' in the section \code{Details}.
#' @param Vector.row.gene Vector of integer indicating the rows of the genes
#' to be plotted.
#' @param Column.gene Integer indicating the column where gene names are given.
#' Set \code{Column.gene=NULL} if there is no such column.
#' @param Group.position Integer indicating the position of group information
#' in the string of characters in each sample names (see \code{Details}).
#' Set \code{Group.position=NULL} if there is only one or no biological
#' information in the string of character in each sample name.
#' @param Time.position Integer indicating the position of time measurement
#' information in the string of characters in each sample names
#' (see \code{Details}).
#' Set \code{Time.position=NULL} if there is only one or no time measurement
#' information in the string of character in each sample name.
#' @param Individual.position Integer indicating the position of the name of
#' the individual (e.g patient, replicate, mouse, yeasts culture ...)
#' in the string of characters in each sample names (see \code{Details}).
#' The names of different individuals must be all different.
#' Furthermore, if individual names are just numbers, they will be transform in
#' a vector of class "character" by [CharacterNumbers()] and
#' a "r" will be added to each individual name ("r" for replicate).
#' @param Color.Group NULL or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute
#' a color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param path.result Character or \code{NULL}. Path to save all results.
#' If \code{path.result} contains a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}" and a sub sub folder,
#' "1-5_ProfileExpression_\code{Name.folder.profile}"
#' all results will be saved in the sub folder
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}/1-5_ProfileExpression_\code{Name.folder.profile}".
#' Otherwise, a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}" and/or a sub sub folder
#' "1-5_ProfileExpression_\code{Name.folder.profile}"
#' will be created in \code{path.result} and all results will be saved in
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}/1-5_ProfileExpression_\code{Name.folder.profile}".
#' If NULL, the results will not be saved in a folder. NULL as default.
#' @param Name.folder.profile Character or \code{NULL}.
#' If \code{Name.folder.profile} is a character, the folder and sub folder names
#' which will contain the PCA graphs will respectively be
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}" and
#' "1-5_ProfileExpression_\code{Name.folder.profile}".
#' Otherwise, the folder and sub folder names will respectively be
#' "1_UnsupervisedAnalysis" and "1-5_ProfileExpression".
#'
#' @return The function plots for each gene selected with
#' the input \code{Vector.row.gene}
#' * In the case where samples belong to different time points only :
#' the evolution of the expression of each replicate across time and
#' the evolution of the mean and the standard deviation of the expression
#' across time.
#' * In the case where samples belong to different biological conditions only :
#' a violin plot (see [ggplot2::geom_violin()]),
#' and error bars (standard deviation) (see [ggplot2::geom_errorbar()])
#' for each biological condition.
#' * In the case where samples belong to different time points and different
#' biological conditions : the evolution of the expression of each replicate
#' across time and the evolution of the mean and the standard deviation
#' of the expression across time for each biological condition.
#'
#' @seealso The function calls the function [DATAplotExpression1Gene()]
#' for each selected genes with \code{Vector.row.gene}.
#'
#' @importFrom grDevices pdf dev.off
#'
#' @export
#'
#' @examples
#' res.sim.count=RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
#'                                   Nb.Gene=10)
#' #
#' Res.evo=DATAplotExpressionGenes(ExprData=res.sim.count$Sim.dat,
#'                                 Vector.row.gene=c(1,3),
#'                                 Column.gene=1,
#'                                 Group.position=1,
#'                                 Time.position=2,
#'                                 Individual.position=3,
#'                                 Color.Group=NULL,
#'                                 path.result=NULL,
#'                                 Name.folder.profile=NULL)
#' print(Res.evo)

DATAplotExpressionGenes<-function(ExprData,
                                  Vector.row.gene,
                                  Column.gene,
                                  Group.position,
                                  Time.position,
                                  Individual.position,
                                  Color.Group=NULL,
                                  path.result=NULL,
                                  Name.folder.profile=NULL){
  #---------------------------------------------------------------------------#
  # Folder creation if no existence
  #---------------------------------------------------------------------------#
  if(is.null(Name.folder.profile)==TRUE){
    Name.folder.profile<-""
    SubFolder.name<-"1_UnsupervisedAnalysis"
  }else{
    Name.folder.profile<-paste("_",Name.folder.profile,sep="")
    SubFolder.name<-paste("1_UnsupervisedAnalysis", Name.folder.profile,sep="")
  }# if(is.null(Name.folder.profile)==TRUE)
  #
  if(is.null(path.result)==FALSE){
    if(SubFolder.name%in%dir(path = path.result)==FALSE){
      print("Folder creation")
      dir.create(path=paste(path.result,"/",SubFolder.name,sep=""))
      path.result.f<-paste(path.result,"/",SubFolder.name,sep="")
    }else{
      path.result.f<-paste(path.result,"/",SubFolder.name,sep="")
    }
  }else{
    path.result.f<-NULL
  }# if(is.null(path.result)==FALSE)
  #
  if(is.null(path.result.f)==FALSE){
    nom.dossier.result<-paste("1-5_ProfileExpressionAnalysis",
                              Name.folder.profile, sep="")
    if(nom.dossier.result%in%dir(path = path.result.f)==FALSE){
      dir.create(path=paste(path.result.f,"/",nom.dossier.result,sep=""))
      path.result.new<-paste(path.result.f,"/",nom.dossier.result,sep="")
    }else{
      path.result.new<-paste(path.result.f,"/",nom.dossier.result,sep="")
    }# if(nom.dossier.result%in%dir(path = path.result.f)==FALSE)
  }else{
    path.result.new<-NULL
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  # Name all genes
  if(is.null(Column.gene)==TRUE){
    if(is.null(row.names(ExprData))==TRUE){
      Name.G<-paste("Gene.",Vector.row.gene,sep="")
    }else{
      Name.G<-row.names(ExprData)[Vector.row.gene]
    }
  }else{
    Name.G<-ExprData[Vector.row.gene,Column.gene]
  }# if(is.null(Column.gene)==TRUE)
  #---------------------------------------------------------------------------#
  # Data with only expression
  if(is.null(Column.gene)==TRUE){
    ind.col.expr<-seq_len(ncol(ExprData))#c(1:ncol(ExprData))
  }else{
    ind.col.expr<-seq_len(ncol(ExprData))[-Column.gene]
  }# if(is.null(Column.gene)==TRUE)
  #---------------------------------------------------------------------------#
  # Data expression of the selected genes
  ExprData.f<-data.frame(Gene=Name.G,
                         as.matrix(ExprData[Vector.row.gene,ind.col.expr]))
  row.names(ExprData.f)<-Name.G
  colnames(ExprData.f)[-1]<-colnames(ExprData)[ind.col.expr]
  #---------------------------------------------------------------------------#
  List.All.G<-vector(mode="list",length=length(Vector.row.gene))
  names(List.All.G)<-Name.G
  #---------------------------------------------------------------------------#
  # Save of all graph in a pdf file
  if(is.null(path.result)==FALSE){
    grDevices::pdf(paste(path.result.new, "/PlotsProfileGeneExpression",
                         Name.folder.profile, ".pdf",sep=""),
                   width = 11, height = 8,
                   onefile = TRUE)
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  cpt<-0
  for(g.sel in Vector.row.gene){
    cpt<-cpt+1
    PlotExpr1G<-DATAplotExpression1Gene(ExprData=ExprData.f,
                                        row.gene=cpt,
                                        Column.gene=1,
                                        Group.position=Group.position,
                                        Time.position=Time.position,
                                        Individual.position=Individual.position,
                                        Color.Group=Color.Group)
    List.All.G[[cpt]]<-PlotExpr1G
    print(PlotExpr1G)
  }# for(g.sel in Vector.row.gene)
  #---------------------------------------------------------------------------#
  if(is.null(path.result)==FALSE){
    grDevices::dev.off()
  }
  #---------------------------------------------------------------------------#
  return(list(List.plots=List.All.G,
              DataForPlot=ExprData.f))
}# DATAplotExpressionGenes
