#' @title Reshaped dataset for factorial analysis.
#'
#' @description The function generates a dataset reshaped from the original
#' dataset, to be used by the function [FactoMineR::PCA()],
#' which performs the Principal Component Analysis (PCA).
#' This function is called by the function [PCArealization()],
#' which also calls the function [FactoMineR::PCA()].
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
#' All these sample information must be separated by underscores
#' in the sample name. For instance 'CLL_P_t0_r1',
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
#' @param Individual.position Integer indicating the position of the name
#' of the individual (e.g patient, replicate, mouse, yeasts culture ...)
#' in the string of characters in each sample names (see \code{Details}).
#' The names of different individuals must be all different.
#' Furthermore, if individual names are just numbers, they will be transform in
#' a vector of class "character" by [CharacterNumbers()] and
#' a "r" will be added to each individual name ("r" for replicate).
#'
#' @seealso The function is called by the function [PCArealization()].
#'
#' @return A reshape of the originally dataset which corresponds to
#' a data.frame with (\eqn{N_g+k}) columns and \eqn{N_s} rows, where
#' \eqn{N_g} is the number of genes, \eqn{N_s} is the number of samples and
#' * \eqn{k=1} if samples belong to different biological condition or
#' time points.
#' In that case, the first column will contain the biological condition
#' or the time point associated to each sample.
#' * \eqn{k=2} if samples belong to different biological condition
#' and time points.
#' In that case, the first column will contain the biological condition
#' and the second column the time point associated to each sample.
#'
#' The other \eqn{N_g} columns form a sub data.frame which is a transpose of
#' the data.frame composed of the \eqn{N_s} numeric columns of \code{ExprData}.
#'
#' @export
#'
#' @examples
#' res.sim.count<-RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
#'                                    Nb.Gene=10)
#' #--------------------------------------------------------------------------#
#' res.dat.PCA<-PCApreprocessing(ExprData=res.sim.count$Sim.dat,
#'                               Column.gene=1,
#'                               Group.position=1,
#'                               Time.position=2,
#'                               Individual.position=3)

PCApreprocessing<-function(ExprData,
                           Column.gene,
                           Group.position,
                           Time.position,
                           Individual.position){
  #---------------------------------------------------------------------------#
  res.Factors<-ColnamesToFactors(ExprData=ExprData,
                                 Column.gene=Column.gene,
                                 Group.position=Group.position,
                                 Time.position=Time.position,
                                 Individual.position=Individual.position)
  #
  Vector.group<-res.Factors$Group.Info
  Vector.time.ini<-res.Factors$Time.Info
  Vector.patient<-res.Factors$Individual.info
  #
  if(is.null(Vector.time.ini)==FALSE){
    Tt.Del<-gsub("t","",gsub("T","",as.character(Vector.time.ini)))
    Vector.time<-paste("t",Tt.Del,sep="")
  }else{
    Vector.time<-Vector.time.ini
  }# if(is.null(Vector.time.ini)==FALSE)
  #---------------------------------------------------------------------------#
  null.index.vector<-which(c(is.null(Vector.group), is.null(Vector.time)))
  if(length(null.index.vector)==2){
    stop("You need a qualitative variable")
  }# if(length(null.index.vector)==2)
  #---------------------------------------------------------------------------#
  if(length(null.index.vector)==0){
    final.list<-list(Quali.Sup.Group=as.factor(Vector.group),
                     Quali.Sup.Time=as.factor(Vector.time))
    paste.quali.var<-do.call("paste", c(final.list, sep = "_"))
  }else{
    final.list<-list(Quali.Sup.Group=as.factor(Vector.group),
                     Quali.Sup.Time=as.factor(Vector.time))[-null.index.vector]
    paste.quali.var<-as.character(unlist(final.list))
  }# if(length(null.index.vector)==0)
  #---------------------------------------------------------------------------#
  if(is.null(Column.gene)==TRUE){
    data.f<-cbind.data.frame(final.list, as.data.frame(t(ExprData)))
  }else{
    data.f<-cbind.data.frame(final.list,
                             as.data.frame(t(ExprData[,-Column.gene])))
  }# if(is.null(Column.gene)==TRUE)
  row.names(data.f)<-paste(Vector.patient, "_", paste.quali.var,sep="")
  #---------------------------------------------------------------------------#
  Id.colname.gene<--seq_len(ncol(data.f)-nrow(ExprData))
  # c(-1:(-ncol(data.f)+nrow(ExprData)))
  Nb.unique.gene<-length(unique(ExprData[,Column.gene]))
  if(is.null(Column.gene)==FALSE & Nb.unique.gene==nrow(ExprData)){
    colnames(data.f)[Id.colname.gene]<-as.character(ExprData[,Column.gene])
  }else{
    colnames(data.f)[Id.colname.gene]<-paste("Gene.",
                                             seq_len(nrow(ExprData)), sep="")
  }# if(is.null(Column.gene)==FALSE & Nb.unique.gene==nrow(ExprData))
  #---------------------------------------------------------------------------#
  order.row<-seq_len(length(Vector.patient))
  #---------------------------------------------------------------------------#
  return(list(data.to.pca=data.f[order.row,],
              nb.quali.var=length(final.list),
              List.Factors=list(Vector.group=Vector.group[order.row],
                                Vector.time=Vector.time[order.row],
                                Vector.patient=Vector.patient[order.row])))
}# PCApreprocessing()
