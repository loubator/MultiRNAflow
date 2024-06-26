#' @title Data preparation for exploratory and statistical analysis
#' (Main Function)
#'
#' @description
#' This function creates automatically  a SummarizedExperiment (SE) object
#' from raw counts data to store
#' * information for exploratory (unsupervised) analysis using the R function
#' [SummarizedExperiment::SummarizedExperiment()]
#' * a DESeq2 object from raw counts data in order to store all information
#' for statistical (supervised) analysis using the R function
#' [DESeq2::DESeqDataSetFromMatrix()].
#'
#' @details The column names of \code{RawCounts} must be a vector of strings
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
#' If the user does not have all these sample information separated by
#' underscores in the sample name, the user can build the data.frame
#' \code{colData} describing the samples.
#'
#' @param RawCounts Data.frame with \eqn{N_g} rows and (\eqn{N_{s}+k}) columns,
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
#' Furthermore, if individual names are just numbers, they will be transform
#' in a vector of class "character" by
#' [CharacterNumbers()] and
#' a "r" will be added to each individual name ("r" for replicate).
#' @param colData \code{NULL} or data.frame with \eqn{N_s} rows and two or
#' three columns describing the samples. \code{NULL} as default.
#' Optional input (see \code{Details}).
#' If \code{Group.position}, \code{Time.position} and
#' \code{Individual.position} are filled, set \code{colData=NULL}.
#' * If samples belong to different times point and different biological
#' condition
#'   * the first column must contain the biological condition for each sample.
#'   The column name must be "Group".
#'   * the second column must contain the time measurement for each sample.
#'   The column name must be "Time".
#'   * The third column must contain the individual name for each sample.
#'   The column name must be "ID".
#' * If samples belong to different times point or different biological
#' condition
#'   * the first column must contain, either the biological condition
#'   for each sample, or the time measurement for each sample.
#'   The column name must be either "Group", or "Time".
#'   * The second column must contain the individual name for each sample.
#'   The column name must be "ID".
#' @param VARfilter Positive numeric value, 0 as default.
#' All rows of \code{RawCounts} which the variance of counts is strictly under
#' the threshold \code{VARfilter} are deleted
#' @param SUMfilter Positive numeric value, 0 as default.
#' All rows of \code{RawCounts} which the sum of counts is strictly under
#' the threshold \code{SUMfilter} are deleted.
#' @param RNAlength \code{NULL} or "hsapiens" or data.frame with two columns.
#' \code{NULL} as default.
#' * if \code{RNAlength} is a data.frame
#'   * the first column must contain gene names
#'   (similar to those of \code{RawCounts})
#'   * the second columns must contain the median of the transcript length
#'   of each gene of the first column
#'   and all rows of \code{RawCounts} whose genes are not included in
#'   the first column of \code{RNAlength} will be deleted.
#' * if \code{RNAlength=NULL}, no rows will be deleted.
#'
#' If \code{RNAlength} is either "hsapiens" or a data.frame,
#' \code{Column.gene} can not be \code{NULL}.
#'
#' @return The function returns a SummarizedExperiment object containing
#' all information for exploratory (unsupervised) analysis and
#' DE statistical analysis.
#'
#' @seealso The [DATAprepSE()] function
#' * is used by the following functions of our package :
#' [DATAnormalization()],
#' [DEanalysisGlobal()].
#' * calls the R function
#' [DESeq2::DESeqDataSetFromMatrix()]
#' in order to create the DESeq2 object and
#' [SummarizedExperiment::SummarizedExperiment()]
#' in order to create the SummarizedExperiment object
#'
#' @importFrom stats as.formula reformulate
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors metadata
#'
#' @export
#'
#' @examples
#' BgCdEx <- rep(c("P", "NP"), each=27)
#' TimeEx <- rep(paste0("t", seq_len(9) - 1), times=6)
#' IndvEx <- rep(paste0("pcl", seq_len(6)), each=9)
#'
#' SampleNAMEex <- paste(BgCdEx, IndvEx, TimeEx, sep="_")
#' RawCountEx <- data.frame(Gene.name=paste0("Name", seq_len(10)),
#'                          matrix(sample(seq_len(100),
#'                                        length(SampleNAMEex)*10,
#'                                        replace=TRUE),
#'                                 ncol=length(SampleNAMEex), nrow=10))
#' colnames(RawCountEx) <- c("Gene.name", SampleNAMEex)
#' ##------------------------------------------------------------------------##
#' resDATAprepSE <- DATAprepSE(RawCounts=RawCountEx,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=3,
#'                             Individual.position=2)
#' ##
#' ## colDataEx <- data.frame(Group=BgCdEx, Time=TimeEx, ID=IndvEx)

DATAprepSE <- function(RawCounts,
                       Column.gene,
                       Group.position,
                       Time.position,
                       Individual.position,
                       colData=NULL,
                       VARfilter=0,
                       SUMfilter=0,
                       RNAlength=NULL) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check inputs
    resErr <- ErrDATAprepSE(RawCounts=RawCounts, Column.gene=Column.gene,
                            Group.position=Group.position,
                            Time.position=Time.position,
                            Individual.position=Individual.position,
                            colData=colData, RNAlength=RNAlength,
                            VARfilter=VARfilter, SUMfilter=SUMfilter)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## RNA length filter
    resLength <- subRNAfilter(RawCounts=RawCounts, Column.gene=Column.gene,
                              RNAlength=RNAlength)
    RawCounts2 <- resLength$RawCounts
    RNAfilter2 <- resLength$RNAfilter

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## RNA count filter
    resCount <- subCOUNTfilter(RawCounts=RawCounts2, Column.gene=Column.gene,
                               VARfilter=VARfilter, SUMfilter=SUMfilter)
    RawCounts2 <- resCount$RawCounts
    RNAfilter3 <- append(resCount$COUNTfilter, RNAfilter2)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Pre-processing
    if (!is.null(colData)) {
        dataColToFct <- data.frame(RawCounts2)[c(1, 2),]
        colData1v <- as.character(apply(colData, 1, paste, collapse="_"))

        if (nrow(colData) == ncol(RawCounts2)) {
            colnames(dataColToFct) <- colData1v
        } else {

            if (nrow(colData) != ncol(RawCounts2)-1) {
                Err_nrowcolData <- paste("The number of rows of 'colData'",
                                         "must be equal to the number of",
                                         "samples (numbers of column (Nc) of",
                                         "'RawCounts' if 'Column.gene==NULL',",
                                         "Nc - 1 otherwise).")
                stop(Err_nrowcolData)
            }## if (nrow(colData) != ncol(RawCounts2)-1)

            ColNameColToFct <- c(colData1v, NA)
            SeqGeneFct <- seq(Column.gene, nrow(colData))
            ColNameColToFct[SeqGeneFct + 1] <- ColNameColToFct[SeqGeneFct]
            ColNameColToFct[Column.gene] <- "Gene"

            colnames(dataColToFct) <- ColNameColToFct
        }## if(nrow(colData) == ncol(RawCounts2))

        if (c("Group")%in%colnames(colData)) {
            Gposition <- which(colnames(colData)%in%c("Group"))
        } else {
            Gposition <- NULL
        }## if (c("Group")%in%colnames(colData))

        if (c("Time")%in%colnames(colData)) {
            Tposition <- which(colnames(colData)%in%c("Time"))
        } else {
            Tposition <- NULL
        }## if (c("Time")%in%colnames(colData))

        Iposition <- ncol(colData)

    } else {
        dataColToFct <- RawCounts2
        Gposition <- Group.position
        Tposition <- Time.position
        Iposition <- Individual.position
    }## if (!is.null(colData))

    res.Factors <- ColnamesToFactors(ExprData=dataColToFct,
                                     Column.gene=Column.gene,
                                     Group.position=Gposition,
                                     Time.position=Tposition,
                                     Individual.position=Iposition)
    Vect.group <- res.Factors$Group.Info
    Vect.time <- res.Factors$Time.Info
    Vect.id <- res.Factors$Individual.info

    ##-----------------------------------------------------------------------##
    ## Biological conditions and time present
    if (!is.null(Vect.group) & !is.null(Vect.time)) {
        Vect.time <- gsub("T" , "", gsub("t", "", as.character(Vect.time)))
        colData.DESeq2 <- data.frame(Group=as.factor(Vect.group),
                                     Time=as.factor(Vect.time))
        design.DESeq2 <- stats::as.formula(~ Time + Group + Time:Group)
    }## if(!is.null(Vect.group) & !is.null(Vect.time))

    ##-----------------------------------------------------------------------##
    ## Biological condition present & Time absent
    if (!is.null(Vect.group) & is.null(Vect.time)) {
        colData.DESeq2 <- data.frame(Group=as.factor(Vect.group))
        design.DESeq2 <- stats::as.formula(~ Group)
    }## if(!is.null(Vect.group) & is.null(Vect.time))

    ##-----------------------------------------------------------------------##
    ## Biological conditions absent & Time present
    if (is.null(Vect.group) & !is.null(Vect.time)) {
        Vect.time <- gsub("T", "", gsub("t", "", as.character(Vect.time)))
        colData.DESeq2 <- data.frame(Time=as.factor(Vect.time))
        design.DESeq2 <- stats::as.formula(~ Time)
    }## if(is.null(Vect.group) & !is.null(Vect.time))

    ##-----------------------------------------------------------------------##
    ## Data with only expression
    if (is.null(Column.gene)) {
        ind.col.expr <- seq_len(ncol(RawCounts2))
        RowNamesRawCounts <- paste0("Gene", seq_len(nrow(RawCounts2)))
    } else {
        ind.col.expr <- seq_len(ncol(RawCounts2))[-Column.gene]
        RowNamesRawCounts <- RawCounts2[, Column.gene]
    }## if (is.null(Column.gene))

    mat.Data <- as.matrix(RawCounts2[, ind.col.expr])
    row.names(mat.Data) <- RowNamesRawCounts
    colnames(mat.Data) <- res.Factors$Final.Name

    ##-----------------------------------------------------------------------##
    ## Preparation of SummarizedExperiment object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=mat.Data,
                                          colData=colData.DESeq2,
                                          design=design.DESeq2)

    colDataSE <- cbind(data.frame(colData.DESeq2, ID=factor(Vect.id)),
                       res.Factors$Data.code.names)


    SEformula <- design.DESeq2
    ## SEformula <- stats::reformulate(as.character(design.DESeq2)[-1],
    ##                                 response="counts")
    colINFOnb <- ncol(colData.DESeq2) + 1
    nameINFO <- seq(from=colINFOnb+1, to=colINFOnb+2, by=1)
    listColDataINFO <- list(colINFOfactors=seq_len(colINFOnb),
                            colINFOname=nameINFO)

    listResults_ini <- list(UnsupervisedAnalysis=list(Normalization=NULL,
                                                      PCA=NULL,
                                                      HCPC=NULL,
                                                      Mfuzz=NULL,
                                                      GenesExpression=NULL),
                            SupervisedAnalysis=list(Normalization=NULL,
                                                    DEanalysis=NULL,
                                                    VolcanoMAplots=NULL,
                                                    Heatmaps=NULL,
                                                    Rgprofiler2=NULL))

    ##-----------------------------------------------------------------------##
    ## Creation of SummarizedExperiment object ## design.DESeq2
    SEobj <- SEobjFUN(mat.Data, colDataSE)

    S4Vectors::metadata(SEobj)$RAWcolnames <- colnames(RawCounts)
    S4Vectors::metadata(SEobj)$colGene <- Column.gene
    S4Vectors::metadata(SEobj)$colDataINFO <- listColDataINFO
    S4Vectors::metadata(SEobj)$RNAfiltering <- RNAfilter3
    S4Vectors::metadata(SEobj)$DESeq2obj <- list(formula=SEformula,
                                                 DESeq2preproceesing=dds)
    S4Vectors::metadata(SEobj)$Results <- listResults_ini
    S4Vectors::metadata(SEobj)$SEidentification <- c("SEstep")

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Outputs ## Data.Expression=mat.Data,
    ## Factors.Info=data.frame(colData.DESeq2, Samples=Vect.id)
    ## Data.code.names=res.Factors$Data.code.names,
    return(SEobj=SEobj)
}## DATAprepSE()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

SEobjFUN <- function(x, y, z=NULL) {
    rSE <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=x),
                                                      colData=y,
                                                      rowData=z)
    return(resSE=rSE)
}## SEobjFUN()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

subCOUNTfilter <- function(RawCounts,
                           Column.gene,
                           VARfilter=0,
                           SUMfilter=0) {
    ##-----------------------------------------------------------------------##
    ## Step 1
    matCount <- RawCounts

    if (!is.null(Column.gene)) {
        matCount <- matCount[,-Column.gene]
        Genes <- as.character(matCount[,Column.gene])
    } else {
        Genes <- seq_len(nrow(matCount))
    }## if (is.null(Column.gene))

    COUNTfilter <- list(SUMdeleted=NULL, VARdeleted=NULL)

    ##-----------------------------------------------------------------------##
    ## Counts filter
    idVar <- idSum <- c()

    if (SUMfilter > 0) {
        idSum <- which(apply(matCount, 1 , sum) < SUMfilter)
    }## if (SUMfilter>0)

    if (VARfilter > 0) {
        idVar <- which(apply(matCount, 1 , var) < VARfilter)
    }## if (VARfilter>0)

    ##-----------------------------------------------------------------------##
    ## Filtering
    if (length(idSum) + length(idVar) > 0) {
        idFilter <- sort(unique(c(idSum, idVar)))
        RawCounts2 <- RawCounts[-idFilter,]

        if (length(idSum) > 0) {
            COUNTfilter[[1]] <- Genes[-idSum]
        }## if (length(idSum) > 0)
        if (length(idVar) > 0) {
            COUNTfilter[[2]] <- Genes[-idVar]
        }## if (length(idVar) > 0)

    } else {
        RawCounts2 <- RawCounts
    }## if (length(idSum) + length(idVar) > 0)

    ##-----------------------------------------------------------------------##
    return(list(RawCounts=RawCounts2,
                COUNTfilter=COUNTfilter))
}## subCOUNTfilter()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

subRNAfilter <- function(RawCounts, Column.gene, RNAlength=NULL) {
    ##-----------------------------------------------------------------------##
    ## Sub RawCounts
    if (!is.null(RNAlength)) {

        hspns <- "hsapiens"
        hspns01 <- identical(RNAlength, hspns)

        if (isTRUE(hspns01)) {
            RNAlength <- get(utils::data("Transcript_HomoSapiens_Database",
                                         package="MultiRNAflow"))
        }## if (isTRUE(hspns01))

        gene_rwctINtrlg <- which(RawCounts[,Column.gene]%in%RNAlength[,1])
        Nsimilar <- length(gene_rwctINtrlg)

        if (Nsimilar == 0) {
            Err_RNA5 <- paste0("No corresponding genes.")
            stop(Err_RNA5)
        }## if (Nsimilar == 0)

        if (Nsimilar == nrow(RawCounts) | is.null(RNAlength)) {
            RawCounts2 <- RawCounts
            delGenes <- "No deleted Genes"
        } else {
            RawCounts2 <- RawCounts[gene_rwctINtrlg,]
            delGenes <- as.character(RawCounts[-gene_rwctINtrlg, Column.gene])

            Ndel <- nrow(RawCounts) - Nsimilar
            print(paste0(Ndel, " genes deleted."))
        }## if (Nsimilar == nrow(RawCounts) | is.null(RNAlength))

        gene_trlgINrwct <- which(RNAlength[,1]%in%RawCounts2[,Column.gene])
        RNAlength2 <- RNAlength[gene_trlgINrwct,]

        RawCounts2 <- RawCounts2[order(RNAlength2[,Column.gene]),]
        RNAlength2 <- RNAlength2[order(RNAlength2[,1]),]

        RNAfilter <- list(LENGTHdeleted=delGenes, RNAlength=RNAlength2)
    } else {
        RawCounts2 <- RawCounts
        RNAlength2 <- NULL
        RNAfilter <- list(LENGTHdeleted=NULL, RNAlength=NULL)
    }## if (!is.null(RNAlength))

    ##-----------------------------------------------------------------------##
    return(list(RawCounts=RawCounts2,
                RNAfilter=RNAfilter))
}## subRNAfilter()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrcolData <- function(colData) {

    if (!is.data.frame(colData)) {
        stop("'colData' must be a matrix of class data.frame.")
    }## if (!is.data.frame(colData))
    colData <- data.frame(colData)

    if (ncol(colData) == 1 | ncol(colData) > 3) {
        stop("'colData' must have two or three coulumns")
    }## if (!is.null(colData))

    if (ncol(colData) == 3) {
        columnTest <- identical(colnames(colData),
                                c("Group" ,"Time", "ID"))
        ##
        if (!isTRUE(columnTest)) {
            Err_colname <- paste0("The column names of 'colData' must be ",
                                  "'Group' ,'Time', 'ID'")
            stop(Err_colname)
        }## if (!isTRUE(columnTest))
    }## if (ncol(colData) == 3)

    if (ncol(colData) == 2) {
        columnTestG <- identical(colnames(colData), c("Group" ,"ID"))
        columnTestT <- identical(colnames(colData), c("Time", "ID"))
        ##
        if (!isTRUE(columnTestG) & !isTRUE(columnTestT)) {
            Err_colname <- paste0("The column names of 'colData' must be ",
                                  "either 'Group', 'ID'",
                                  "either 'Time', 'ID'")
            stop(Err_colname)
        }## if (!isTRUE(columnTestG) & !isTRUE(columnTestT))
    }## if (ncol(colData) == 2)

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrcolData()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrFilter <- function(VARfilter=0,
                      SUMfilter=0) {
    ##-----------------------------------------------------------------------##
    if (!is.numeric(SUMfilter)) {
        Err_sum <- paste0("'SUMfilter' must be a positive numeric value")
        stop(Err_sum)
    }## if (!is.numeric(SUMfilter))

    if (SUMfilter < 0) {
        Err_sum <- paste0("'SUMfilter' must be a positive numeric value")
        stop(Err_sum)
    }## if (SUMfilter < 0)

    if (!is.numeric(VARfilter)) {
        Err_var <- paste0("'VARfilter' must be a positive numeric value")
        stop(Err_var)
    }## if (!is.numeric(VARfilter))

    if (VARfilter < 0) {
        Err_var <- paste0("'VARfilter' must be a positive numeric value")
        stop(Err_var)
    }## if (VARfilter < 0)

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrFilter()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrRNAlength <- function(Column.gene, RNAlength=NULL) {
    ##-----------------------------------------------------------------------##
    hspns01 <- identical(RNAlength, "hsapiens")

    if (!is.null(RNAlength) & !hspns01 & !is.data.frame(RNAlength)) {
        Err_RNA <- paste0("'RNAlength' must be either NULL, 'hsapiens' ",
                          "or a two columns data.frame.")
        stop(Err_RNA)
    }## if (!is.null(RNAlength) & !hspns01 & !is.data.frame(RNAlength))

    if (!is.null(RNAlength) & is.null(Column.gene)) {
        Err_RNA2 <- paste0("If 'RNAlength' is not NULL, ",
                           "'Column.gene' can not be NULL too.")
        stop(Err_RNA2)
    }## if (!is.null(RNAlength) & is.null(Column.gene))

    if (is.data.frame(RNAlength)) {
        if (ncol(RNAlength) != 2) {
            Err_RNA3 <- paste0("'RNAlength' must be a two columns data.frame.")
            stop(Err_RNA3)
        }## if (ncol(RNAlength) != 2)

        if (!is.character(RNAlength[,1]) | !is.numeric(RNAlength[,2])) {
            Err_RNA4 <- paste0("The first column of 'RNAlength' must be ",
                               "character and the second numeric.")
            stop(Err_RNA4)
        }## if (is.character(RNAlength[,1]) & is.numeric(RNAlength[,2]))
    }## if (is.data.frame(RNAlength))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrRNAlength()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrRawCounts <- function(RawCounts,
                         Column.gene){
    ##-----------------------------------------------------------------------##
    if (!is.data.frame(RawCounts)) {
        stop("'RawCounts' must be a matrix of class data.frame.")
    }## if (!is.data.frame(RawCounts))

    ##-----------------------------------------------------------------------##
    ## check duplicated row.names colnames
    if (!is.null(Column.gene)) {
        Genes <- as.character(RawCounts[, Column.gene])
        dupliGene <- which(duplicated(Genes))
        if (length(dupliGene) > 0) {
            Err_dupliGene <- paste0("There are ", length(dupliGene),
                                    " duplicated genes: ",
                                    unique(Genes[dupliGene]))
            stop(Err_dupliGene)
        }## if (length(dupliGene) > 0)
    }## if (!is.null(Column.gene))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrRawCounts()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrDATAprepSE <- function(RawCounts,
                          Column.gene,
                          Group.position,
                          Time.position,
                          Individual.position,
                          colData=NULL,
                          VARfilter=0,
                          SUMfilter=0,
                          RNAlength=NULL) {
    ##-----------------------------------------------------------------------##
    ## Check 1
    if (is.null(Time.position) & is.null(Group.position) & is.null(colData)) {
        Err_NULL <- paste0("'Time.position', 'Group.position' and 'colData' ",
                           "can not be all 'NULL'")
        stop(Err_NULL)
    }## if(is.null(Time.position) & is.null(Group.position) & is.null(colData))

    ##-----------------------------------------------------------------------##
    ## Check 2
    if (is.null(colData)) {
        ## Every sample must have an individual name
        res_ErrPosition <- ErrPosition(Column.gene=Column.gene,
                                       Group.position=Group.position,
                                       Time.position=Time.position,
                                       Individual.position=Individual.position)
    } else {
        res_colData <- ErrcolData(colData=colData)
    }## if (!is.null(colData))

    ##-----------------------------------------------------------------------##
    res_ErrRawCounts <- ErrRawCounts(RawCounts=RawCounts,
                                     Column.gene=Column.gene)

    ##-----------------------------------------------------------------------##
    ## Check RNAlength, SUMfilter, VARfilter
    res_ErrRNAlength <- ErrRNAlength(Column.gene=Column.gene,
                                     RNAlength=RNAlength)

    res_ErrFilter <- ErrFilter(VARfilter=VARfilter, SUMfilter=SUMfilter)

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrDATAprepSE()

