#' @title Automatic creation of a DESeq2 object for DE analysis.
#'
#' @description
#' This function creates automatically a DESeq2 object from raw counts data
#' using the R function [DESeq2::DESeqDataSetFromMatrix()].
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
#' @param RawCounts Data.frame with \eqn{N_g} rows and (\eqn{N_{s+k}}) columns,
#' where \eqn{N_g} is the number of genes,
#' \eqn{N_s} is the number of samples and
#' \eqn{k=1} if a column is used to specify gene names, or \eqn{k=0} otherwise.
#' If \eqn{k=1}, the position of the column containing gene names is given
#' by \code{Column.gene}.
#' The data.frame contains non negative integers giving gene expressions of
#' each gene in each sample.
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
#' @param Individual.position Integer indicating the position of the name of
#' the individual (e.g patient, replicate, mouse, yeasts culture ...)
#' in the string of characters in each sample names (see \code{Details}).
#' The names of different individuals must be all different.
#'
#' @return The function returns the \code{DESeq2} object for DE analysis and
#' the raw counts data (a matrix of non-negative integers).
#'
#' @seealso The [DEanalysisPreprocessing()] function
#' * is used by the following functions of our package : [DATAnormalization()],
#' [DEanalysisGlobal()].
#' * calls the R function [DESeq2::DESeqDataSetFromMatrix()] in order to create
#' the DESeq2 object.
#'
#' @importFrom stats as.formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#'
#' @export
#'
#' @examples
#' Grp.test=rep(c("P", "NP"), each=27)
#' Tps.test=rep(paste0("t", 0:8), times=6)
#' Pat.test=rep(paste0("pcl", 1:6), each=9)
#'
#' Name.sel.test=paste(Grp.test, Pat.test, Tps.test, sep="_")
#' Mat.test.1=data.frame(Gene.name=paste0("Name", 1:10),
#'                       matrix(sample(1:100,length(Name.sel.test)*10,
#'                                     replace=TRUE),
#'                              ncol=length(Name.sel.test), nrow=10))
#' colnames(Mat.test.1)=c("Gene.name", Name.sel.test)
#' ##-------------------------------------------------------------------------#
#' DESeq2.info.test=DEanalysisPreprocessing(RawCounts=Mat.test.1,
#'                                          Column.gene=1,
#'                                          Group.position=1,
#'                                          Time.position=3,
#'                                          Individual.position=2)
#' print(DESeq2.info.test)

DEanalysisPreprocessing<-function(RawCounts,
                                  Column.gene,
                                  Group.position,
                                  Time.position,
                                  Individual.position){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if(is.null(Time.position) & is.null(Group.position)){
        stop("'Time.position' and 'Group.position' can not be both NULL")
    }# if(is.null(Time.position)==TRUE & is.null(Group.position)==TRUE)

    ##------------------------------------------------------------------------#
    # Pre-processing
    res.Factors<-ColnamesToFactors(ExprData=RawCounts,
                                   Column.gene=Column.gene,
                                   Group.position=Group.position,
                                   Time.position=Time.position,
                                   Individual.position=Individual.position)
    Vect.group<-res.Factors$Group.Info
    Vect.time<-res.Factors$Time.Info

    ##------------------------------------------------------------------------#
    ## Biological conditions and time present
    if(!is.null(Vect.group) & !is.null(Vect.time)){
        Vect.time<-gsub("T","",gsub("t","",as.character(Vect.time)))
        colData.DESeq2<-data.frame(Group=as.factor(Vect.group),
                                   Time=as.factor(Vect.time))
        design.DESeq2<-stats::as.formula(~ Time + Group + Time:Group)
    }# if(is.null(Vect.group)==FALSE & is.null(Vect.time)==FALSE)

    ##------------------------------------------------------------------------#
    ## Biological condition present & Time absent
    if(!is.null(Vect.group) & is.null(Vect.time)){
        colData.DESeq2<-data.frame(Group=as.factor(Vect.group))
        design.DESeq2<-stats::as.formula(~ Group)
    }# if(is.null(Vect.group)==FALSE & is.null(Vect.time)==TRUE)

    ##------------------------------------------------------------------------#
    ## Biological conditions absent & Time present
    if(is.null(Vect.group) & !is.null(Vect.time)){
        Vect.time<-gsub("T", "", gsub("t", "", as.character(Vect.time)))
        colData.DESeq2<-data.frame(Time=as.factor(Vect.time))
        design.DESeq2<-stats::as.formula(~ Time)
    }# if(is.null(Vect.group)==TRUE & is.null(Vect.time)==FALSE)

    ##------------------------------------------------------------------------#
    ## Biological conditions and time absent
    if(is.null(Vect.group) & is.null(Vect.time)){
        colData.DESeq2<-NULL
        design.DESeq2<-stats::as.formula(~ 1)
    }# if(is.null(Vect.group)==TRUE & is.null(Vect.time)==TRUE)

    ##------------------------------------------------------------------------#
    ## Data with only expression
    if(is.null(Column.gene)){
        ind.col.expr<-seq_len(ncol(RawCounts))
        RowNamesRawCounts<-paste0("Gene", seq_len(nrow(RawCounts)))
    }else{
        ind.col.expr<-seq_len(ncol(RawCounts))[-Column.gene]
        RowNamesRawCounts<-RawCounts[,Column.gene]
    }

    mat.Data<-as.matrix(RawCounts[,ind.col.expr])
    row.names(mat.Data)<-RowNamesRawCounts
    colnames(mat.Data)<-res.Factors$Final.Name

    ##------------------------------------------------------------------------#
    ## Creation of Deseq2 object
    dds<-DESeq2::DESeqDataSetFromMatrix(countData=mat.Data,
                                        colData=colData.DESeq2,
                                        design=design.DESeq2)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Outputs
    return(list(DESeq2.obj=dds,
                Data.Expression=mat.Data,
                Factors.Info=data.frame(colData.DESeq2,
                                        Samples=res.Factors$Individual.info)))
}# DEanalysisPreprocessing()
