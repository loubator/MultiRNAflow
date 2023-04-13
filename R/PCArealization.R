#' @title PCA realization
#'
#' @description From a gene expression dataset, the functions performs
#' the Principal Component Analysis (PCA) through the R function
#' [FactoMineR::PCA()].
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
#' 'r1' is localized after the third underscore,
#' so \code{Individual.position=4},
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
#' Set \code{Group.position=NULL} if there is only one or
#' no biological information in the string of character in each sample name.
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
#' @param gene.deletion \code{NULL} or a vector of characters or a vector of
#' integers. \code{NULL} as default.
#' If \code{gene.deletion} is a vector of characters, all genes with names in
#' \code{gene.deletion} will be deleted from \code{ExprData}.
#' If \code{gene.deletion} is a vector of integers,
#' all the corresponding row numbers of \code{ExprData} will be deleted.
#' If \code{gene.deletion=NULL} all genes of \code{ExprData} will be used
#' in the construction of the PCA.
#' @param sample.deletion \code{NULL} or a vector of characters or
#' a vector of integers. \code{NULL} as default.
#' If \code{sample.deletion} is a vector of characters, all samples with names
#' in \code{sample.deletion} will not be used in the construction of the PCA.
#' If \code{sample.deletion} is a vector of integers,
#' all the corresponding column numbers of \code{ExprData} will not be used
#' in the construction of the PCA.
#' If \code{sample.deletion=NULL} all samples will be used
#' in the construction of the PCA.
#' @param Supp.del.sample \code{TRUE} or \code{FALSE}.
#' If \code{FALSE}, the samples selected with \code{sample.deletion} will
#' be deleted.
#' If \code{TRUE}, the samples selected with \code{sample.deletion} will
#' be plotted.
#' These individuals are called supplementary individuals in
#' [FactoMineR::PCA()].
#'
#' @seealso The [PCArealization()] function
#' * is used by the following functions of our package :
#' [PCAanalysis()] and [HCPCanalysis()].
#' * calls the R function [PCApreprocessing()] for reshaping the data and
#' uses its output for performing a Principal Component (PCA)
#' with [FactoMineR::PCA()].
#'
#' @return The function returns the output of the [FactoMineR::PCA()] function
#' (see [FactoMineR::PCA()]).
#'
#' @export
#'
#' @importFrom FactoMineR PCA
#'
#' @examples
#' Sim.Dat.pca<-RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
#'                                  Nb.Gene=10)
#' ##-------------------------------------------------------------------------#
#' res.pca.ex<-PCArealization(ExprData=Sim.Dat.pca$Sim.dat,
#'                            Column.gene=1,
#'                            Group.position=1,
#'                            Time.position=2,
#'                            Individual.position=3,
#'                            gene.deletion=c(3,5),
#'                            sample.deletion=c("G1_t0_Ind2","G1_t1_Ind3"),
#'                            Supp.del.sample=FALSE)
#' ##-------------------------------------------------------------------------#
#' res.pca.ex<-PCArealization(ExprData=Sim.Dat.pca$Sim.dat,
#'                            Column.gene=1,
#'                            Group.position=1,
#'                            Time.position=2,
#'                            Individual.position=3,
#'                            gene.deletion=c("Gene3","Gene5"),
#'                            sample.deletion=c(3,8),
#'                            Supp.del.sample=TRUE)

PCArealization<-function(ExprData,
                         Column.gene,
                         Group.position,
                         Time.position,
                         Individual.position,
                         gene.deletion,
                         sample.deletion,
                         Supp.del.sample=FALSE){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    resPCAprepro<-PCApreprocessing(ExprData=ExprData,
                                   Column.gene=Column.gene,
                                   Group.position=Group.position,
                                   Time.position=Time.position,
                                   Individual.position=Individual.position)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## In this "if" section, we want to know if some samples must be deleted or
    ## be plotted as "supplementray", e.g. not used in the built of the axes
    ## of the PCA

    if(!is.null(sample.deletion)){
        if(is.numeric(sample.deletion)){
            if(is.null(Column.gene)){
                Ind.del.f<-sample.deletion
            }else{
                ColSplDel<-colnames(ExprData)[sample.deletion]
                Ind.del.f<-which(colnames(ExprData)[-Column.gene]%in%ColSplDel)
            }## if(is.null(Column.gene))
        }else{
            if(is.null(Column.gene)){
                Ind.del.f<-which(colnames(ExprData)%in%sample.deletion)
            }else{
                Ind.del.f<-which(colnames(ExprData[,-Column.gene])%in%sample.deletion)
            }## if(is.null(Column.gene))
        }## if(is.numeric(sample.deletion))

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        if(isFALSE(Supp.del.sample)){
            Supp.del.sample.f<-NULL
            data.pca<-resPCAprepro$data.to.pca[-Ind.del.f,]

            ListFactors.F<-resPCAprepro$List.Factors

            for(l in seq_len(length(ListFactors.F))){
                if(!is.null(ListFactors.F[[l]])){
                    VFct<-ListFactors.F[[l]]
                    ListFactors.F[[l]]<-VFct[-Ind.del.f]
                }## if(!is.null(ListFactors.F[[l]]))
            }## for(l in 1:length(ListFactors.F))

        }else{
            Supp.del.sample.f<-Ind.del.f
            data.pca<-resPCAprepro$data.to.pca
            ListFactors.F<-resPCAprepro$List.Factors
        }## if(isFALSE(Supp.del.sample))

        ##--------------------------------------------------------------------#
    }else{
        Supp.del.sample.f<-NULL
        data.pca<-resPCAprepro$data.to.pca
        ListFactors.F<-resPCAprepro$List.Factors
    }## if(!is.null(sample.deletion))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## In this "if" section, we want to know if some genes must be deleted
    if(!is.null(gene.deletion)){
        if(is.numeric(gene.deletion)){
            GeneDel.f<-gene.deletion
        }else{
            GeneDel.f<-which(colnames(resPCAprepro$data.to.pca)%in%gene.deletion)
        }## if(is.numeric(gene.deletion))
        data.pca.f<-data.pca[,-GeneDel.f]
    }else{
        data.pca.f<-data.pca
    }## if(!is.null(gene.deletion))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    res.pca<-FactoMineR::PCA(X=data.pca.f,
                             graph=FALSE,
                             quali.sup=seq_len(resPCAprepro$nb.quali.var),
                             ind.sup=Supp.del.sample.f)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(list(res.pca=res.pca,
                List.Factors=ListFactors.F))
}## PCArealization()
