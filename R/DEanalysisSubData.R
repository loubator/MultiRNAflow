#' @title Sub data of a data.frame
#'
#' @description From an initial data.frame, the function extracts a sub
#' data.frame containing predefined rows.
#'
#' @details
#' If \code{Res.DE.analysis} is a data.frame() or the output from
#' [DEanalysisGlobal()], then
#' * If \code{Set.Operation="union"} then the rows extracted from \code{Data}
#' are those such that the sum of the selected columns by \code{ColumnsCriteria}
#' in \code{Res.DE.analysis} is >0.
#' For example, if \code{Res.DE.analysis} is the outputs from
#' [DEanalysisGlobal()], the rows extracted from \code{Data} will be those DE
#' at least at one time ti (except the reference time t0).
#'
#' * If \code{Set.Operation="intersect"} then the rows extracted from
#' \code{Data} are those such that the product of the selected columns
#' by \code{ColumnsCriteria} in \code{Res.DE.analysis} is >0.
#' For example, if \code{Res.DE.analysis} is the outputs from
#' [DEanalysisGlobal()], the rows extracted from \code{Data} will be those DE
#' at all time ti (except the reference time t0).
#'
#' * If \code{Set.Operation="setdiff"} then the rows extracted from \code{Data}
#' are those such that only one element of the selected columns by
#' \code{ColumnsCriteria} in \code{Res.DE.analysis} is >0.
#' For example, if \code{Res.DE.analysis} is the outputs from
#' [DEanalysisGlobal()], the rows extracted from \code{Data} will be those DE
#' at only one time ti (except the reference time t0).
#'
#' If \code{Res.DE.analysis} is a vector, \code{ColumnsCriteria} and
#' \code{Set.Operation} are not used.
#'
#' @param Data Data.frame with \eqn{N_g} rows and (\eqn{N_s+k}) columns,
#' where \eqn{N_g} is the number of genes,
#' \eqn{N_s} is the number of samples and
#' \eqn{k=1} if a column is used to specify gene names, or \eqn{k=0} otherwise.
#' If \eqn{k=1}, the position of the column containing gene names is given by
#' \code{Column.gene}.
#' The data.frame contains numeric values giving gene expressions of each gene
#' in each sample. Gene expressions can be raw counts or normalized raw counts.
#' Data and \code{Res.DE.analysis} must have the same number of rows.
#' @param Res.DE.analysis A list containing a data.frame or a data.frame or
#' a binary vector.
#' If it is a list, it must be the outputs from [DEanalysisGlobal()]
#' (see \code{Examples}).
#' If it is a data.frame, it must contains at least one binary column
#' (filled with 0 and 1).
#' @param ColumnsCriteria A vector of integers where each integer indicates
#' a column of \code{Res.DE.analysis}.
#' These columns should either contain only binary values, or may contain other
#' numerical value, in which case extracted rows from \code{Data} will be
#' those with >0 values (see \code{Details}).
#' @param Set.Operation A character.
#' The user must choose between "union" (default), "intersect", "setdiff"
#' (see \code{Details}).
#' @param Save.SubData \code{TRUE} or \code{FALSE} or a Character.
#' \code{NULL} as default.
#' If \code{TRUE}, two csv files (see \code{Value}) will be saved in the folder
#' "2_SupervisedAnalysis_\code{Name.folder.DE}" (see [DEanalysisGlobal()]).
#'
#' @return The function returns
#' * A sub data.frame of \code{Data} containing only the rows specified by
#' \code{ColumnsCriteria} and \code{Set.Operation}.
#' * A sub data.frame of \code{Data} containing only the rows specified by
#' \code{ColumnsCriteria} and \code{Set.Operation}.
#' * The rows specified by \code{ColumnsCriteria} and \code{Set.Operation}.
#'
#' @export
#'
#' @examples
#' Data.EX<-data.frame(matrix(sample(x=1:150, size=20), nrow=5))
#' colnames(Data.EX)<-paste("Gene",1:4,sep="")
#' #--------------------------------------------------------------------------#
#' # Exemple of output from the function DEanalysisGlobal()
#' res.DEanalysisGlobal.Ex<-list(data.frame(Gene=paste("Gene",1:5,sep="."),
#'                                          DE1=c(0,1,0,0,1),
#'                                          DE2=c(0,1,0,1,0)))
#' names(res.DEanalysisGlobal.Ex)<-c("DE.results")
#' #--------------------------------------------------------------------------#
#' res.SubDE<-DEanalysisSubData(Data=Data.EX,
#'                              Res.DE.analysis=res.DEanalysisGlobal.Ex,
#'                              ColumnsCriteria=c(2,3),
#'                              Set.Operation="union",
#'                              Save.SubData=FALSE)

DEanalysisSubData<-function(Data,
                            Res.DE.analysis,
                            ColumnsCriteria=1,
                            Set.Operation="union",
                            Save.SubData=FALSE){
  #---------------------------------------------------------------------------#
  if(Set.Operation%in%c("union", "intersect", "setdiff")==FALSE){
    stop("Set.Operation mut be 'union', 'intersect' or 'setdiff'")
  }# if(Set.Operation%in%c("union", "intersect", "setdiff")==FALSE)
  #---------------------------------------------------------------------------#
  if(is.list(Res.DE.analysis)==TRUE){
    DatRowSel<-data.frame(Res.DE.analysis$DE.results)
  }else{
    DatRowSel<-data.frame(Res.DE.analysis)
  }# if(is.list(Res.DE.analysis)==TRUE)
  ncol.DatRowSel<-ncol(DatRowSel)
  nrow.DatRowSel<-nrow(DatRowSel)
  #---------------------------------------------------------------------------#
  if(is.numeric(ColumnsCriteria)==TRUE){
    if(sum(abs(floor(ColumnsCriteria)-ColumnsCriteria))>0){
      stop("'ColumnsCriteria' must be integers or characters")
    }# if(sum(abs(floor(ColumnsCriteria)-ColumnsCriteria))>0)
  }else{
    if(is.character(ColumnsCriteria)==TRUE){
      ColumnsCriteria.2<-rep(NA, times=length(ColumnsCriteria))
      for(i in seq_len(length(ColumnsCriteria))){# 1:length(ColumnsCriteria)
        if(length(grep(pattern=ColumnsCriteria[i], x=colnames(DatRowSel)))==0){
          Stop.WrongNames<-paste("The element ", i,
                                 " of 'ColumnsCriteria' is not correct",sep="")
          stop(Stop.WrongNames)
        }
        ColumnsCriteria.2[i]<-grep(pattern=ColumnsCriteria[i],
                                   x=colnames(DatRowSel))
      }# for(i in 1:length(ColumnsCriteria))
      #
      ColumnsCriteria<-sort(ColumnsCriteria.2)
    }else{
      stop("'ColumnsCriteria' must be integers or characters")
    }# if(is.character(ColumnsCriteria)==TRUE)
  }# if(is.numeric(ColumnsCriteria)==TRUE)
  #---------------------------------------------------------------------------#
  if(ncol.DatRowSel==1){
    ColumnsCriteria<-c(1)
  }# if(ncol.DatRowSel==1)
  #
  if(max(ColumnsCriteria)>ncol.DatRowSel | min(ColumnsCriteria)<1){
    Stop.WrongIntegers<-paste("Integers of 'ColumnsCriteria' must be",
                              "between 1 and", ncol.DatRowSel, sep=" ")
    stop(Stop.WrongIntegers)
  }# if(max(ColumnsCriteria)>ncol.DatRowSel | min(ColumnsCriteria)<1)
  #---------------------------------------------------------------------------#
  Nb.rows.Data<-nrow(Data)
  if(Nb.rows.Data!=nrow.DatRowSel){
    Stop.nrow<-paste("The number of rows of 'Data' and, ",
                     "the length or the number of rows of 'Res.DE.analysis', ",
                     "must be identical.", sep= "")
    stop(Stop.nrow)
  }# if(Nb.rows.Data!=nrow.DatRowSel)
  #---------------------------------------------------------------------------#
  if(length(ColumnsCriteria)==1){
    DEsel<-which(as.numeric(DatRowSel[,ColumnsCriteria])>0)
  }# if(length(ColumnsCriteria)==1)
  #
  if(Set.Operation=="union" & length(ColumnsCriteria)>1){
    Sum.colsel<-apply(data.frame(DatRowSel[,ColumnsCriteria]), 1, sum)
    DEsel<-which(as.numeric(Sum.colsel)>0)
  }
  #
  if(Set.Operation=="intersect" & length(ColumnsCriteria)>1){
    Prod.colsel<-apply(data.frame(DatRowSel[,ColumnsCriteria]), 1, prod)
    DEsel<-which(as.numeric(Prod.colsel)>0)
  }
  #
  if(Set.Operation=="setdiff" & length(ColumnsCriteria)>1){
    Nb0.colsel<-apply(X=data.frame(DatRowSel[,ColumnsCriteria]), MARGIN=1,
                      FUN=function(x) length(which(x==0)))
    DEsel<-which(as.numeric(Nb0.colsel)==(length(ColumnsCriteria)-1))
  }
  #
  L.DEsel<-length(DEsel)
  #---------------------------------------------------------------------------#
  if(L.DEsel==0){
    print(paste("No selection because the column selected is full of 0.",
                "The original 'Data' is returned."), sep="")
    DataSub<-Data
  }# if(L.DEsel==0)
  if(L.DEsel==Nb.rows.Data){
    print(paste("All rows are selected because there is no 0",
                "in the column selected.", "The original 'Data' is returned.",
                sep=" "))
    DataSub<-Data
  }# if(L.DEsel==Nb.rows.Data)
  if(L.DEsel>0 & L.DEsel<Nb.rows.Data){
    DataSub<-Data[DEsel,]
  }# if(L.DEsel>0 & L.DEsel<Nb.rows.Data)
  #---------------------------------------------------------------------------#
  # Folder path and creation
  if(isFALSE(Save.SubData)==FALSE){
    if(Save.SubData==TRUE){
      path.result<-Res.DE.analysis$Path.result
    }else{
      path.result<-Save.SubData
    }# if(Save.SubData==TRUE)
    utils::write.table(DataSub,
                       file=paste(path.result,"/Sub_Data.csv",sep=""),
                       sep=";", row.names=FALSE)
    utils::write.table(DatRowSel,
                       file=paste(path.result,"/Sub_DEresults.csv",sep=""),
                       sep=";", row.names=FALSE)
  }else{
    path.result<-NULL
  }# if(isFALSE(Save.SubData)==FALSE)
  #---------------------------------------------------------------------------#
  return(list(SubData=DataSub,
              SubDataCriteria=DatRowSel,
              RowsSelected=DEsel))
}# DEanalysisSubData()
