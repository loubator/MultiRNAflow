#' @title PCA realization
#'
#' @description From a gene expression dataset, the functions performs
#' the Principal Component Analysis (PCA) through the R function
#' [FactoMineR::PCA()].
#'
#' @details All results are built from the results of our function
#' [DATAnormalization()].
#'
#' @param SEresNorm Results of the function
#' [DATAnormalization()].
#' @param DATAnorm \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' \code{TRUE} means the function uses the normalized data.
#' \code{FALSE} means the function uses the raw counts data.
#' @param gene.deletion \code{NULL} or a vector of characters or a vector of
#' integers. \code{NULL} as default.
#' If \code{gene.deletion} is a vector of characters, all genes with names in
#' \code{gene.deletion} will be deleted from the data set as input
#' \code{RawCounts} of our function
#' [DATAprepSE()].
#' If \code{gene.deletion} is a vector of integers,
#' all the corresponding row numbers will be deleted from the data set as input
#' \code{RawCounts} of our function
#' [DATAprepSE()].
#' If \code{gene.deletion=NULL} all genes will be used in the construction
#' of the PCA.
#' @param sample.deletion \code{NULL} or a vector of characters or
#' a vector of integers. \code{NULL} as default.
#' If \code{sample.deletion} is a vector of characters, all samples with names
#' in \code{sample.deletion} will not be used in the construction of the PCA.
#' If \code{sample.deletion} is a vector of integers,
#' all the corresponding column numbers will not be used in the construction
#' of the PCA from the data set as input \code{RawCounts} of our function
#' [DATAprepSE()].
#' If \code{sample.deletion=NULL} all samples will be used
#' in the construction of the PCA.
#' @param Supp.del.sample \code{TRUE} or \code{FALSE}. \code{FALSE} by default.
#' If \code{FALSE}, the samples selected with \code{sample.deletion} will
#' be deleted.
#' If \code{TRUE}, the samples selected with \code{sample.deletion} will
#' be plotted.
#' These individuals are called supplementary individuals in
#' [FactoMineR::PCA()].
#'
#' @seealso The [PCArealization()] function
#' * is used by the following functions of our package :
#' [PCAanalysis()] and
#' [HCPCanalysis()].
#' * calls the R function
#' [PCApreprocessing()]
#' for reshaping the data and
#' uses its output for performing a Principal Component (PCA)
#' with
#' [FactoMineR::PCA()].
#'
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEresNorm} but with the output of the
#' [FactoMineR::PCA()]
#' function (see
#' [FactoMineR::PCA()])
#' saved in the metadata \code{Results[[1]][[2]]} of \code{SEresNorm}.
#'
#' @export
#'
#' @importFrom FactoMineR PCA
#' @importFrom S4Vectors metadata
#'
#' @examples
#' ## Simulation raw counts
#' resSIMcount <- RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
#'                                    Nb.Gene=10)
#' ## Preprocessing step
#' resDATAprepSE <- DATAprepSE(RawCounts=resSIMcount$Sim.dat,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=2,
#'                             Individual.position=3)
#' ## Normalization
#' resNorm <- DATAnormalization(SEres=resDATAprepSE,
#'                              Normalization="rle",
#'                              Plot.Boxplot=FALSE,
#'                              Colored.By.Factors=FALSE)
#' ##------------------------------------------------------------------------##
#' resPCAex <- PCArealization(SEresNorm=resNorm,
#'                            DATAnorm=TRUE,
#'                            gene.deletion=c(3, 5),
#'                            sample.deletion=c("G1_t0_Ind2", "G1_t1_Ind3"),
#'                            Supp.del.sample=FALSE)
#' ##------------------------------------------------------------------------##
#' resPCAex2 <- PCArealization(SEresNorm=resNorm,
#'                             DATAnorm=TRUE,
#'                             gene.deletion=c("Gene3", "Gene5"),
#'                             sample.deletion=c(3, 8),
#'                             Supp.del.sample=TRUE)

PCArealization <- function(SEresNorm,
                           DATAnorm=TRUE,
                           gene.deletion=NULL,
                           sample.deletion=NULL,
                           Supp.del.sample=FALSE) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check inputs
    resErr <- ErrPCArealization(SEresNorm=SEresNorm,
                                DATAnorm=DATAnorm,
                                gene.deletion=gene.deletion,
                                sample.deletion=sample.deletion,
                                Supp.del.sample=Supp.del.sample)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    SEresPCAprepro <- PCApreprocessing(SEresNorm=SEresNorm,
                                       DATAnorm=DATAnorm)

    PCAlist <- S4Vectors::metadata(SEresPCAprepro)$Results[[1]][[2]]
    nPCAlist <- length(PCAlist)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## In this "if" section, we want to know if some samples must be deleted or
    ## be plotted as "supplementray", e.g. not used in the built of the axes
    ## of the PCA

    colG <- S4Vectors::metadata(SEresPCAprepro)$colGene
    RAWcolN <- S4Vectors::metadata(SEresPCAprepro)$RAWcolnames

    ## Default PCA: sample.deletion == NULL
    data.pca <- PCAlist$data.to.pca
    NBquali <- PCAlist$nb.quali.var
    ListFactors.F <- PCAlist$List.Factors

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## In this "if" section, we want to know if some genes must be deleted
    if (!is.null(gene.deletion)) {
        if (is.numeric(gene.deletion)) {
            GeneDel.f <- gene.deletion
        } else {
            GeneDel.f <- which(colnames(data.pca)%in%gene.deletion)

            if (length(GeneDel.f) != length(gene.deletion)) {
                NOhereG <- paste(which(!gene.deletion%in%colnames(data.pca)),
                                 collapse=", ")
                Err_gene <- paste("The genes", NOhereG,
                                  "selected are not present in your dataset.")
                stop(Err_gene)
            }
        }## if(is.numeric(gene.deletion))
        data.pca.f <- data.pca[, -GeneDel.f]
    } else {
        data.pca.f <- data.pca
    }## if(!is.null(gene.deletion))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    IDquali <- seq_len(NBquali)
    Supp.del.sample.f <- NULL

    if (!is.null(sample.deletion)) {
        if (is.numeric(sample.deletion)) {
            if (is.null(colG)) {
                Ind.del.f <- sample.deletion
            } else {
                ColSplDel <- RAWcolN[sample.deletion]
                Ind.del.f <- which(RAWcolN[-colG]%in%ColSplDel)
            }## if(is.null(colG))
        } else {
            if (is.null(colG)) {
                Ind.del.f <- which(RAWcolN%in%sample.deletion)
            } else {
                Ind.del.f <- which(RAWcolN[-colG]%in%sample.deletion)
            }## if(is.null(colG))
        }## if(is.numeric(sample.deletion))

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        if (isFALSE(Supp.del.sample)) {
            data.pca.f <- data.pca.f[-Ind.del.f,]

            for (l in seq_len(length(ListFactors.F))) {
                if (!is.null(ListFactors.F[[l]])) {
                    VFct <- ListFactors.F[[l]]
                    ListFactors.F[[l]] <- VFct[-Ind.del.f]
                }## if(!is.null(ListFactors.F[[l]]))
            }## for(l in 1:length(ListFactors.F))

        } else {
            Supp.del.sample.f <- Ind.del.f
        }## if(isFALSE(Supp.del.sample))
    }## if(!is.null(sample.deletion))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    resPCA <- FactoMineR::PCA(X=data.pca.f,
                              graph=FALSE,
                              quali.sup=IDquali,
                              ind.sup=Supp.del.sample.f)

    SEresPCA <- SEresPCAprepro
    PCAlist <- append(PCAlist, list(PCAresults=resPCA))
    PCAlist$List.Factors <- ListFactors.F

    S4Vectors::metadata(SEresPCA)$Results[[1]][[2]] <- PCAlist
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    return(SEobj=SEresPCA)
}## PCArealization()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrPCArealization <- function(SEresNorm,
                              DATAnorm=TRUE,
                              gene.deletion=NULL,
                              sample.deletion=NULL,
                              Supp.del.sample=FALSE) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check SEresNorm (cf DATAplotExpressionGenes())
    Err1_PCAres <- ErrSEresNorm(SEresNorm=SEresNorm,
                                DATAnorm=DATAnorm,
                                path.result=NULL)

    ##-----------------------------------------------------------------------##
    RAWcolN <- S4Vectors::metadata(SEresNorm)$RAWcolnames

    Err2_PCAres <- ErrDeletion(gene.deletion=gene.deletion,
                               sample.deletion=sample.deletion,
                               Supp.del.sample=Supp.del.sample,
                               RAWcolnames=RAWcolN)

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrPCArealization()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrDeletion <- function(gene.deletion=NULL,
                        sample.deletion=NULL,
                        Supp.del.sample=FALSE,
                        RAWcolnames="sample") {
    ##-----------------------------------------------------------------------##
    if (!isTRUE(Supp.del.sample) & !isFALSE(Supp.del.sample)) {
        stop("'Supp.del.sample' must be TRUE or FALSE.")
    }## if (!isTRUE(Supp.del.sample) & !isFALSE(Supp.del.sample))

    ##-----------------------------------------------------------------------##
    ## Check deletion
    if (!is.null(gene.deletion)) {
        ErrG_delChNum <- ErrDeletionChNum(gene.deletion, "gene.deletion")
    }## if (!is.null(gene.deletion))

    ## Check deletion
    if (!is.null(sample.deletion)) {
        ErrS_delChNum <- ErrDeletionChNum(sample.deletion, "sample.deletion")

        if (is.character(sample.deletion)) {
            noSmpl <- which(!sample.deletion%in%RAWcolnames)
            if (length(noSmpl) > 0) {
                no_delspl <- paste(sample.deletion[noSmpl], collapse="', '")
                Err_delspl <- paste0("The elements '", no_delspl, "' are not ",
                                     "present among the column names of the ",
                                     "raw count data.")
                stop(Err_delspl)
            }## if (length(noSmpl) > 0)
        }## if (!is.character(sample.deletion))
    }## if (!is.null(sample.deletion))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrDeletion()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrDeletionChNum <- function(TOdelete, NAMEdelete) {
    ##-----------------------------------------------------------------------##
    if (!is.character(TOdelete)) {
        Err_delint <- paste0("'", NAMEdelete, "' must be either NULL ",
                             "either character or non negative integers.")
        if (!is.numeric(TOdelete)) {
            stop(Err_delint)
        } else {
            if (sum(abs(floor(TOdelete) - TOdelete)) !=0) {
                stop(Err_delint)
            }## if (floor(TOdelete) != TOdelete)
            if (min(TOdelete) <= 0) {
                stop(Err_delint)
            }## if (min(TOdelete) <= 0)
        }## if (!is.numeric(TOdelete))
    }## if (!is.character(TOdelete))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrDeletionChNum()

