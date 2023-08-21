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
#' @return The function returns a SummarizedExperiment class object containing
#' information and a reshape of the originally dataset for the PCA analysis.
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
#' ##-------------------------------------------------------------------------#
#' resPCAdata <- PCApreprocessing(SEresNorm=resNorm,
#'                                DATAnorm=TRUE)

PCApreprocessing <- function(SEresNorm,
                             DATAnorm=TRUE) {
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check 1
    ## DATAprepSE
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    if (!is(SEresNorm, "SummarizedExperiment")) {
        stop(Err_SE)
    } else {
        codeDEres <- S4Vectors::metadata(SEresNorm)$SEidentification

        if (is.null(codeDEres)) {
            stop(Err_SE)
        }## if (is.null(codeDEres))

        if (codeDEres != "SEresNormalization") {
            stop(Err_SE)
        }## if (codeDEres != "SEresNormalization")
    }## if (!is(SEresNorm, "SummarizedExperiment"))

    if (!isTRUE(DATAnorm) & !isFALSE(DATAnorm)) {
        stop("'DATAnorm' must be TRUE or FALSE.")
    }## if (!isTRUE(DATAnorm) & !isFALSE(DATAnorm))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## preprocssing
    if (DATAnorm == TRUE) {
        aSE <- 2
    } else {
        aSE <- 1
    }## if (DATAnorm == TRUE)

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

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if (!is.null(Vector.time.ini)) {
        Tt.Del <- gsub("t", "", gsub("T", "", as.character(Vector.time.ini)))
        Vector.time <- paste0("t", Tt.Del)
    }else{
        Vector.time <- Vector.time.ini
    }## if(!is.null(Vector.time.ini))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    null.index.vector<-which(c(is.null(Vector.group), is.null(Vector.time)))
    ## if (length(null.index.vector) == 2) {
    ##     stop("You need a qualitative variable")
    ## }## if(length(null.index.vector)==2)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    final.list <- list(Quali.Sup.Group=as.factor(Vector.group),
                       Quali.Sup.Time=as.factor(Vector.time))

    if (length(null.index.vector) == 0) {
        paste.quali.var <- do.call("paste", c(final.list, sep="_"))
    } else {
        final.list <- final.list[-null.index.vector]
        paste.quali.var <- as.character(unlist(final.list))
    }## if(length(null.index.vector)==0)

    RownamesPCA <- paste0(Vector.patient, "_", paste.quali.var)

    data.f <- cbind.data.frame(final.list,
                               as.data.frame(t(ExprData)))

    row.names(data.f) <- RownamesPCA

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    Id.colname.gene <- - seq_len(ncol(data.f) - nrow(ExprData))
    Nb.unique.gene <- length(unique(NameG))

    if (Nb.unique.gene == nrow(ExprData)) {
        colnames(data.f)[Id.colname.gene] <- NameG
    } else {
        colnames(data.f)[Id.colname.gene] <- paste0("Gene.",
                                                    seq_len(nrow(ExprData)))
    }## if(Nb.unique.gene == nrow(ExprData))

    ## order.row <- seq_len(length(Vector.patient))
    listFCTRS <- list(Vector.group=Vector.group,
                      Vector.time=Vector.time,
                      Vector.patient=Vector.patient)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## SE object
    NBcolINFO <- length(unlist(S4Vectors::metadata(SEresNorm)$colDataINFO))

    SEprepPCA <- SEresNorm
    SummarizedExperiment::colData(SEprepPCA)$PCA.name <- RownamesPCA
    S4Vectors::metadata(SEprepPCA)$PCA <- list(data.to.pca=data.f,
                                               nb.quali.var=length(final.list),
                                               List.Factors=listFCTRS)
    S4Vectors::metadata(SEprepPCA)$colDataINFO$colINFOnamePCA <- NBcolINFO + 1

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(SEobj=SEprepPCA)
}## PCApreprocessing()
