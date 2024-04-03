#' @title Sub data of a data.frame
#'
#' @description From the results from our function
#' [DEanalysisGlobal()],
#' the function extracts from the SummarizedExperiment class outputs of
#' the subset of genes selected with the inputs \code{Set.Operation} and
#' \code{ColumnsCriteria}, and saves them in a SummarizeExperiment object.
#'
#' @details We have the following three cases:
#' * If \code{Set.Operation="union"} then the rows extracted from
#' the different datasets included in \code{SEresDE}
#' are those such that the sum of the selected columns of
#' \code{SummarizedExperiment::rowData(SEresDE)}
#' by \code{ColumnsCriteria} is >0.
#' For example, the selected genes can be those DE at least at t1 or t2
#' (versus the reference time t0).
#' * If \code{Set.Operation="intersect"} then the rows extracted from
#' the different datasets included in \code{SEresDE}
#' are those such that the product of the selected columns of
#' \code{SummarizedExperiment::rowData(SEresDE)}
#' by \code{ColumnsCriteria} is >0.
#' For example, the selected genes can be those DE at times t1 and t2
#' (versus the reference time t0).
#' * If \code{Set.Operation="setdiff"} then the rows extracted from
#' the different datasets included in \code{SEresDE}
#' are those such that only one element of the selected columns of
#' \code{SummarizedExperiment::rowData(SEresDE)}
#' by \code{ColumnsCriteria} is >0.
#' For example, the selected genes can be those DE at times t1 only and
#' at times t2 only (versus the reference time t0).
#'
#' @param SEresDE A SummarizedExperiment class object. Output from
#' [DEanalysisGlobal()]
#' (see \code{Examples}).
#' @param ColumnsCriteria A vector of integers where each integer indicates
#' a column of  \code{SummarizedExperiment::rowData(SEresDE)}.
#' These columns should either contain only binary values, or may contain other
#' numerical value, in which case extracted outputs from \code{SEresDE}
#' will be those with >0 values (see \code{Details}).
#' @param Set.Operation A character.
#' The user must choose between "union" (default), "intersect", "setdiff"
#' (see \code{Details}).
#' @param Save.SubData \code{TRUE} or \code{FALSE} or a Character.
#' \code{FALSE} as default.
#' If \code{TRUE}, two csv files (see \code{Value}) will be saved in the folder
#' "2_SupervisedAnalysis_\code{Name.folder.DE}"
#' (see [DEanalysisGlobal()]).
#'
#' @return The function returns a SummarizeExperiment class object containing
#' * sub data.frames of the different dataset included in \code{SEresDE}
#' containing only the rows specified by
#' \code{ColumnsCriteria} and \code{Set.Operation}.
#' * the DE results saved in \code{SEresDE} of genes selected by
#' \code{ColumnsCriteria} and \code{Set.Operation}.
#' * The genes specified by \code{ColumnsCriteria} and \code{Set.Operation}.
#'
#' @importFrom SummarizedExperiment rowData colData rownames assays
#' @importFrom S4Vectors metadata
#' @importFrom utils write.table

#'
#' @export
#'
#' @examples
#' ## Simulation raw counts
#' resSIMcount <- RawCountsSimulation(Nb.Group=2, Nb.Time=1, Nb.per.GT=4,
#'                                    Nb.Gene=5)
#' ## Preprocessing step
#' resDATAprepSE <- DATAprepSE(RawCounts=resSIMcount$Sim.dat,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=NULL,
#'                             Individual.position=2)
#'
#' ##------------------------------------------------------------------------##
#' ## Transformation of resDATAprepSE into results of DEanalysisGlobal
#' resultsExamples <- data.frame(Gene=paste0("Gene", seq_len(5)),
#'                               DE1=c(0, 1, 0, 0, 1),
#'                               DE2=c(0, 1, 0, 1, 0))
#' listPATHnameEx <- list(Path.result=NULL, Folder.result=NULL)
#'
#' SummarizedExperiment::rowData(resDATAprepSE) <- resultsExamples
#' S4Vectors::metadata(resDATAprepSE)$DESeq2obj$pathNAME <- listPATHnameEx
#' S4Vectors::metadata(resDATAprepSE)$DESeq2obj$SEidentification<-"SEresultsDE"
#'
#' ##------------------------------------------------------------------------##
#' ## results of DEanalysisSubData
#' resDEsub <- DEanalysisSubData(SEresDE=resDATAprepSE,
#'                               ColumnsCriteria=c(2, 3),
#'                               Set.Operation="union",
#'                               Save.SubData=FALSE)

DEanalysisSubData <- function(SEresDE,
                              ColumnsCriteria=1,
                              Set.Operation="union",
                              Save.SubData=FALSE){
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check 1
    ## DATAprepSE
    Err_SE <- paste0("'SEresDE' mut be the results of the function ",
                     "'DEanalysisGlobal()'.")

    ## DATAprepSE
    if (!is(SEresDE, "SummarizedExperiment")) {
        stop(Err_SE)
    } else {
        DESeq2objList <- S4Vectors::metadata(SEresDE)$DESeq2obj
        codeDEres <- DESeq2objList$SEidentification

        if (is.null(codeDEres)) {
            stop(Err_SE)
        }## if (is.null(codeDEres))

        if (codeDEres != "SEresultsDE") {
            stop(Err_SE)
        }## if (codeDEres != SEresultsDE))
    }## if (!is(SEresDE, "SummarizedExperiment"))

    ##-----------------------------------------------------------------------##
    ## Check
    if (!Set.Operation%in%c("union", "intersect", "setdiff")) {
        stop("Set.Operation mut be 'union', 'intersect' or 'setdiff'")
    }## if(Set.Operation%in%c("union", "intersect", "setdiff")==FALSE)

    ##-----------------------------------------------------------------------##
    ## S4Vectors::metadata(SEresDE)$DESeq2obj$pathNAME
    resDESeq2obj <- S4Vectors::metadata(SEresDE)$DESeq2obj
    resPATH <- resDESeq2obj$pathNAME
    DatRowSel <- SummarizedExperiment::rowData(SEresDE)
    nrow.DatRowSel <- Nb.rows.Data <- nrow(DatRowSel)
    ncol.DatRowSel <- ncol(DatRowSel)

    ##-----------------------------------------------------------------------##
    if (is.numeric(ColumnsCriteria)) {

        if (sum(abs(floor(ColumnsCriteria)-ColumnsCriteria))>0) {
            stop("'ColumnsCriteria' must be integers or characters")
        }## if(sum(abs(floor(ColumnsCriteria)-ColumnsCriteria))>0)

    } else {

        if (is.character(ColumnsCriteria)) {
            ColumnsCriteria.2 <- rep(NA, times=length(ColumnsCriteria))

            for (i in seq_len(length(ColumnsCriteria))) {
                lenCOLcriteria <- length(grep(pattern=ColumnsCriteria[i],
                                              x=colnames(DatRowSel),
                                              fixed=TRUE))
                if (lenCOLcriteria == 0) {
                    Stop.WrongNames <- paste("The element", i, "of",
                                             "'ColumnsCriteria'",
                                             "is not correct.")
                    stop(Stop.WrongNames)
                }## if(length(grep(pattern=ColumnsCriteria[i],
                ### x=colnames(DatRowSel)))==0)

                ColumnsCriteria.2[i] <- grep(pattern=ColumnsCriteria[i],
                                             x=colnames(DatRowSel),
                                             fixed=TRUE)
            }## for(i in 1:length(ColumnsCriteria))

            ColumnsCriteria <- sort(ColumnsCriteria.2)
        } else {
            stop("'ColumnsCriteria' must be integers or characters")
        }## if(is.character(ColumnsCriteria)==TRUE)
    }## if(is.numeric(ColumnsCriteria)==TRUE)

    ##-----------------------------------------------------------------------##
    ## if (ncol.DatRowSel == 1) {
    ##     ColumnsCriteria <- c(1)
    ## }## if(ncol.DatRowSel==1)

    if (max(ColumnsCriteria)>ncol.DatRowSel | min(ColumnsCriteria)<1) {
        Stop.WrongIntegers <- paste("Integers of 'ColumnsCriteria' must be",
                                    "between 1 and", ncol.DatRowSel)
        stop(Stop.WrongIntegers)
    }## if(max(ColumnsCriteria)>ncol.DatRowSel | min(ColumnsCriteria)<1)

    ##-----------------------------------------------------------------------##
    if (length(ColumnsCriteria) == 1) {
        DEsel <- which(as.numeric(DatRowSel[, ColumnsCriteria])>0)
    }## if(length(ColumnsCriteria)==1)

    if (Set.Operation == "union" & length(ColumnsCriteria)>1) {
        Sum.colsel <- apply(X=data.frame(DatRowSel[, ColumnsCriteria]),
                            MARGIN=1,
                            FUN=sum)
        DEsel <- which(as.numeric(Sum.colsel)>0)
    }## if(Set.Operation=="union" & length(ColumnsCriteria)>1)

    if (Set.Operation == "intersect" & length(ColumnsCriteria)>1) {
        Prod.colsel <- apply(X=data.frame(DatRowSel[, ColumnsCriteria]),
                             MARGIN=1,
                             FUN=prod)
        DEsel <- which(as.numeric(Prod.colsel)>0)
    }## if(Set.Operation=="intersect" & length(ColumnsCriteria)>1)

    if (Set.Operation == "setdiff" & length(ColumnsCriteria)>1) {
        Nb0.colsel <- apply(X=data.frame(DatRowSel[,ColumnsCriteria]),
                            MARGIN=1,
                            FUN=function(x) length(which(x==0)))
        DEsel <- which(as.numeric(Nb0.colsel) == (length(ColumnsCriteria)-1))
    }## if (Set.Operation == "setdiff" & length(ColumnsCriteria)>1)

    L.DEsel <- length(DEsel)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    listDEselected <- list(DEselectedGenes=DEsel,
                           ColumnsCriteria=ColumnsCriteria,
                           Set.Operation=Set.Operation)

    if (L.DEsel == 0) {
        print(paste("No selection because the column selected is full of 0.",
                    "The original 'SEresDE' is returned."))
        SEresSEsub <- SEresDE
        S4Vectors::metadata(SEresSEsub)$DEselection <- listDEselected
    }## if(L.DEsel == 0)

    if (L.DEsel == Nb.rows.Data) {
        print(paste("All rows are selected because there is no 0",
                    "in the column selected.",
                    "The original 'SEresDE' is returned."))
        SEresSEsub <- SEresDE
        S4Vectors::metadata(SEresSEsub)$DEselection <- listDEselected
    }## if(L.DEsel == Nb.rows.Data)

    if (L.DEsel>0 & L.DEsel<Nb.rows.Data) {
        colFCTRS <- S4Vectors::metadata(SEresDE)$colDataINFO$colINFOfactors
        colFCTRS <- as.numeric(colFCTRS)
        colSEsub <- SummarizedExperiment::colData(SEresDE)[, colFCTRS]
        geneNames1 <- SummarizedExperiment::rownames(SEresDE)

        ## COUNTsub <- SummarizedExperiment::assays(resDESeq2obj$DESeq2results)
        ## COUNTsub <- COUNTsub$counts
        COUNTsub <- SummarizedExperiment::assays(SEresDE)$counts
        COUNTsub <- data.frame(Gene=geneNames1[DEsel], COUNTsub[DEsel,])

        RLEsub <- SummarizedExperiment::assays(SEresDE)$rle
        colNAMESini <- S4Vectors::metadata(SEresDE)$RAWcolnames
        colGene1 <- S4Vectors::metadata(SEresDE)$colGene

        SEresSEsub <- DATAprepSE(RawCounts=COUNTsub,
                                 Column.gene=1,
                                 Group.position=NULL,
                                 Time.position=NULL,
                                 Individual.position=NULL,
                                 colData=data.frame(colSEsub))

        SummarizedExperiment::assays(SEresSEsub)$rle <- RLEsub[DEsel,]
        S4Vectors::metadata(SEresSEsub)$RAWcolnames <- colNAMESini
        S4Vectors::metadata(SEresSEsub)$colGene <- colGene1
        S4Vectors::metadata(SEresSEsub)$DEselection <- listDEselected
        S4Vectors::metadata(SEresSEsub)$SEidentification <- "SEresNormalization"

        SummarizedExperiment::rowData(SEresSEsub) <- DatRowSel[DEsel,]

    }## if(L.DEsel>0 & L.DEsel<Nb.rows.Data)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Folder path and creation
    if (!isFALSE(Save.SubData)) {
        if (isTRUE(Save.SubData)) {
            path.result <- resPATH$Path.result
        } else {
            path.result <- Save.SubData
        }## if(Save.SubData==TRUE)

        if (!is.null(path.result)) {
            utils::write.table(SummarizedExperiment::assays(SEresSEsub)$rle,
                               file=file.path(path.result, "Sub_DataRLE.csv"),
                               sep=";", row.names=FALSE)

            utils::write.table(SummarizedExperiment::assays(SEresSEsub)$counts,
                               file=file.path(path.result,
                                              "Sub_RawCountsData.csv"),
                               sep=";", row.names=FALSE)

            utils::write.table(DatRowSel,
                               file=file.path(path.result, "Sub_DEresults.csv"),
                               sep=";", row.names=FALSE)
        }## if (!is.null(path.result))
    } else {
        path.result <- NULL
    }## if(isFALSE(Save.SubData)==FALSE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    ## list(SubData, SubDataCriteria=DatRowSel, RowsSelected=DEsel)
    return(SEresSEsub)
}## DEanalysisSubData()
