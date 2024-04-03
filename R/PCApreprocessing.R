#' @title Reshaped dataset for factorial analysis.
#'
#' @description The function generates a SummarizedExperiment class object
#' containing the dataset reshaped from the original dataset,
#' to be used by the function
#' [FactoMineR::PCA()],
#' which performs the Principal Component Analysis (PCA).
#' This function is called by the function
#' [PCArealization()],
#' which also calls the function
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
#'
#'
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEresNorm} with the different elements below
#' * information for the functions
#' \code{PCArealization()} and \code{PCAgraphics()}
#' * a reshape of the originally dataset for the PCA analysis
#' (realized by the function \code{PCArealization()})
#'
#' saved in the metadata \code{Results[[1]][[2]]} of \code{SEresNorm}.
#'
#' The reshaped dataset which corresponds to a data.frame with
#' (\eqn{N_g+k}) columns and \eqn{N_s} rows, where
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
#' @seealso The function is called by our function
#' [PCArealization()]
#' and uses our function
#' [DATAnormalization()].
#'
#' @importFrom SummarizedExperiment colData rownames assays
#' @importFrom S4Vectors metadata
#'
#' @export
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
#' resPCAdata <- PCApreprocessing(SEresNorm=resNorm,
#'                                DATAnorm=TRUE)

PCApreprocessing <- function(SEresNorm,
                             DATAnorm=TRUE) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check inputs
    resErr <- ErrPCApreprocessing(SEresNorm=SEresNorm, DATAnorm=DATAnorm)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## preprocssing
    if (isTRUE(DATAnorm)) {
        aSE <- 2
    } else {
        aSE <- 1
    }## if (isTRUE(DATAnorm))

    ExprData <- SummarizedExperiment::assays(SEresNorm)[[aSE]]
    ExprData <- data.frame(ExprData)

    NameG <- as.character(SummarizedExperiment::rownames(SEresNorm))
    cSEdat <- SummarizedExperiment::colData(SEresNorm)

    if (c("Group")%in%colnames(cSEdat)) {
        Vector.group <- cSEdat$Group
    } else {
        Vector.group <- NULL
    }## if (c("Group")%in%colnames(cSEdat))

    if (c("Time")%in%colnames(cSEdat)) {
        Vector.time.ini <- cSEdat$Time
    } else {
        Vector.time.ini <- NULL
    }## if (c("Time")%in%colnames(cSEdat))

    Vector.patient <- cSEdat$ID

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (!is.null(Vector.time.ini)) {
        Tt.Del <- gsub("t", "", gsub("T", "", as.character(Vector.time.ini)))
        Vector.time <- paste0("t", Tt.Del)
    }else{
        Vector.time <- Vector.time.ini
    }## if(!is.null(Vector.time.ini))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    null.index.vector <- which(c(is.null(Vector.group), is.null(Vector.time)))
    ## if (length(null.index.vector) == 2) {
    ##     stop("You need a qualitative variable")
    ## }## if(length(null.index.vector)==2)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    final.list <- list(Quali.Sup.Group=as.factor(Vector.group),
                       Quali.Sup.Time=as.factor(Vector.time))

    if (length(null.index.vector) == 0) {
        paste.quali.var <- do.call("paste", c(final.list, sep="_"))
    } else {
        final.list <- final.list[-null.index.vector]
        paste.quali.var <- as.character(unlist(final.list))
    }## if(length(null.index.vector) == 0)

    RownamesPCA <- paste0(Vector.patient, "_", paste.quali.var)

    data.f <- cbind.data.frame(final.list,
                               as.data.frame(t(ExprData)))
    row.names(data.f) <- RownamesPCA

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    Id.colname.gene <- - seq_len(ncol(data.f) - nrow(ExprData))
    colnames(data.f)[Id.colname.gene] <- NameG

    ## order.row <- seq_len(length(Vector.patient))
    listFCTRS <- list(Vector.group=Vector.group,
                      Vector.time=Vector.time,
                      Vector.patient=Vector.patient)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## SE object
    NBcolINFO <- length(unlist(S4Vectors::metadata(SEresNorm)$colDataINFO))

    SEprepPCA <- SEresNorm
    SummarizedExperiment::colData(SEprepPCA)$PCA.name <- RownamesPCA
    S4Vectors::metadata(SEprepPCA)$colDataINFO$colINFOnamePCA <- NBcolINFO + 1

    PCAlist <- list(data.to.pca=data.f,
                    nb.quali.var=length(final.list),
                    List.Factors=listFCTRS)
    S4Vectors::metadata(SEprepPCA)$Results[[1]][[2]] <- PCAlist

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    return(SEobj=SEprepPCA)
}## PCApreprocessing()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrPCApreprocessing <- function(SEresNorm,
                                DATAnorm=TRUE) {
    ##-----------------------------------------------------------------------##
    ## Check SEresNorm (cf DATAplotExpressionGenes())
    res_ErrSEresNorm <- ErrSEresNorm(SEresNorm=SEresNorm, DATAnorm=DATAnorm,
                                     path.result=NULL)
    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrPCApreprocessing()
