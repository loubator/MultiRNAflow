#' @title Normalization of raw counts (Main Function).
#'
#' @description From raw counts, this function realizes one of
#' the three methods of normalization of the package \code{DESeq2}:
#' * Relative Log Expression (rle) transformation
#' (see [BiocGenerics::estimateSizeFactors()])
#' * Regularized Log (rlog) transformation
#' (see [DESeq2::rlog()])
#' * Variance Stabilizing Transformation (vst) transformation
#' (see [DESeq2::vst()])
#'
#' @details All results are built from the results of the function
#' [DATAprepSE()].
#'
#' @param SEres Results of the function
#' [DATAprepSE()].
#' @param Normalization "rle", "vst", "rlog" and "rpkm".
#' "rle", "vst" and "rlog" correspond to a method of normalization proposed
#' by \code{DESeq2} (see
#' [BiocGenerics::estimateSizeFactors()] for "rle",
#' [DESeq2::rlog()] for "rlog" and
#' [DESeq2::vst()] for "vst").
#' "rpkm" corresponds to a RPKM normalization after a "rle" normalization.
#' @param Blind.rlog.vst \code{TRUE} or \code{FALSE}. \code{FALSE} by default.
#' See input 'blind' in
#' [DESeq2::rlog()].
#' It is recommended to set \code{Blind.rlog.vst=FALSE}
#' for downstream analysis.
#' @param Plot.Boxplot \code{TRUE} or \code{FALSE}. \code{TRUE} by default.
#' If \code{Plot.Boxplot=TRUE}, the function
#' [DATAplotBoxplotSamples()] will be
#' called and boxplots will be plotted. Otherwise, no boxplots will be plotted.
#' @param Colored.By.Factors \code{TRUE} or \code{FALSE}.
#' \code{FALSE} by default.
#' If \code{TRUE}, boxplots will be colored with different colors for different
#' time measurements (if data were collected at different time points).
#' Otherwise, boxplots will be colored with different colors for different
#' biological conditions.
#' @param Color.Group \code{NULL} or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' \code{NULL} by default.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute a
#' color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param Plot.genes \code{TRUE} or \code{FALSE}. \code{FALSE} by default.
#' If \code{TRUE}, points representing gene expressions
#' (normalized or raw counts) will be plotted for each sample.
#' Otherwise, only boxplots will be plotted.
#' @param path.result Character or \code{NULL}. \code{NULL} by default.
#' Path to save all results.
#' If \code{path.result} contains a sub folder entitled
#' "1_Normalization_\code{Name.folder.norm}"
#' all results will be saved in the sub folder
#' "1_Normalization_\code{Name.folder.norm}".
#' Otherwise, a sub folder entitled "1_Normalization_\code{Name.folder.norm}"
#' will be created in \code{path.result} and all results will be saved in
#' "1_Normalization_\code{Name.folder.norm}".
#' If \code{NULL}, the results will not be saved in a folder.
#' \code{NULL} as default.
#' @param Name.folder.norm Character or \code{NULL}. \code{NULL} by default.
#' If \code{Name.folder.norm} is a character,
#' the folder name which will contain the results will be
#' "1_Normalization_\code{Name.folder.norm}".
#' Otherwise, the folder name will be "1_Normalization".
#'
#' @return The function returns a SummarizedExperiment object (\code{SEresNorm})
#' identical as \code{SEres} but
#' * with the normalized count data saved in \code{assays(SEresNORM)[[2]]}
#' * with the boxplot of normalized count
#' saved in the metadata \code{Results[[1]][[1]]} of \code{SEresNorm}.
#'
#' The boxplot is plotted if \code{Plot.Boxplot=TRUE}.
#'
#' @seealso The [DATAnormalization()]
#' function calls our R function
#' [DATAprepSE()],
#' and the R functions
#' [BiocGenerics::estimateSizeFactors()],
#' [DESeq2::rlog()] and
#' [DESeq2::vst()]
#' in order to realized the normalization.
#'
#' @importFrom DESeq2 vst varianceStabilizingTransformation rlog
#' estimateSizeFactors counts
#' @importFrom SummarizedExperiment assay assays colData colnames rownames
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' data(RawCounts_Antoszewski2022_MOUSEsub500)
#' ##------------------------------------------------------------------------##
#' resDATAprepSE <- DATAprepSE(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=NULL,
#'                             Individual.position=2)
#' ##------------------------------------------------------------------------##
#' resNorm <- DATAnormalization(SEres=resDATAprepSE,
#'                              Normalization="rle",
#'                              Plot.Boxplot=TRUE,
#'                              Colored.By.Factors=TRUE)

DATAnormalization <- function(SEres,
                              Normalization="vst",
                              Blind.rlog.vst=FALSE,
                              Plot.Boxplot=TRUE,
                              Colored.By.Factors=FALSE,
                              Color.Group=NULL,
                              Plot.genes=FALSE,
                              path.result=NULL,
                              Name.folder.norm=NULL) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check inputs
    resErr1 <- ErrDATAnormalization1(SEres=SEres,
                                     Normalization=Normalization,
                                     Blind.rlog.vst=Blind.rlog.vst)

    resErr2 <- ErrDATAnormalization2(Plot.Boxplot=Plot.Boxplot,
                                     Colored.By.Factors=Colored.By.Factors,
                                     Color.Group=Color.Group,
                                     Plot.genes=Plot.genes,
                                     path.result=path.result,
                                     Name.folder.norm=Name.folder.norm)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Folder creation if no existence
    if (is.null(Name.folder.norm)) {
        Name.folder.norm <- ""
        SubFolder.name <- "1_UnsupervisedAnalysis"
    } else {
        Name.folder.norm <- paste0("_", Name.folder.norm)
        SubFolder.name <- paste0("1_UnsupervisedAnalysis", Name.folder.norm)
    }## if(is.null(Name.folder.norm))

    nameRESfolder <- paste0("1-1_Normalization", Name.folder.norm)

    if (!is.null(path.result)) {
        if (!SubFolder.name%in%dir(path=path.result)) {
            print("Folder creation")
            dir.create(path=file.path(path.result, SubFolder.name))
        }## if(!SubFolder.name%in%dir(path=path.result))
        ##
        path.result.f <- file.path(path.result, SubFolder.name)
    } else {
        path.result.f <- NULL
    }## if(!is.null(path.result))

    if (!is.null(path.result.f)) {
        if (!nameRESfolder%in%dir(path=path.result.f)) {
            dir.create(path=file.path(path.result.f, nameRESfolder))
        }## if(!nameRESfolder%in%dir(path=path.result.f))
        ##
        path.result.new <- file.path(path.result.f, nameRESfolder)
    } else {
        path.result.new<-NULL
    }## if(!is.null(path.result.f))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Preprocessing, SummarizedExperiment::rowData, dimnames(DESeq2.obj)
    ## colnames(RawCounts), row.names(RawCounts) ## FactorBoxplt
    RawCounts <- SummarizedExperiment::assays(SEres)[[1]]
    SEsampleName <- SummarizedExperiment::colnames(SEres)
    Name.G <- SummarizedExperiment::rownames(SEres)
    DESeq2.obj <- S4Vectors::metadata(SEres)$DESeq2obj$DESeq2preproceesing

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    # Normalization vst
    if (Normalization == "vst") {
        if (nrow(RawCounts) > 1350) {
            log2.data.prepare <- DESeq2::vst(DESeq2.obj, blind=Blind.rlog.vst)
        } else {
            log2.data.prepare <- vst2(DESeq2.obj, Blind.rlog.vst)
        }## if (nrow(RawCounts) > 1350)

        log2.data <- data.frame(SummarizedExperiment::assay(log2.data.prepare))
        log2.data <- round(log2.data, digits=3)
        ##
        Norm.dat <- data.frame(Gene=Name.G, log2.data)
        row.names(Norm.dat) <- Name.G
        ##
        Log2Trf <- FALSE
        YlabelNorm <- "vst normalized counts"
    }## if(Normalization == "vst")

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Normalization rlog
    if (Normalization == "rlog") {
        log2.data.prepare <- DESeq2::rlog(DESeq2.obj, blind=Blind.rlog.vst)
        log2.data <- data.frame(SummarizedExperiment::assay(log2.data.prepare))
        log2.data <- round(log2.data, digits=3)
        ##
        Norm.dat <- data.frame(Gene=Name.G, log2.data)
        row.names(Norm.dat) <- Name.G
        ##
        Log2Trf <- FALSE
        YlabelNorm <- "rlog normalized counts"
    }## if(Normalization == "rlog")

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Normalization rle
    if (Normalization == "rle") {
        dds.SF <- DESeq2::estimateSizeFactors(DESeq2.obj)
        rle.data <- round(DESeq2::counts(dds.SF, normalized=TRUE), digits=3)
        Norm.dat <- data.frame(Gene=Name.G, rle.data)
        row.names(Norm.dat) <- Name.G
        ##
        Log2Trf <- TRUE
        YlabelNorm <- "log2 (rle normalized counts + 1)"
    }## if(Normalization == "rle")

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Normalization rpkm
    if (Normalization == "rpkm") {
        RNAfilter <- S4Vectors::metadata(SEres)$RNAfiltering$RNAlength

        if (is.null(RNAfilter)) {
            Err_RNA <- paste0("The input 'RNAlength' of the function ",
                              "'DATAprepSE' can not be NULL.")
            stop(Err_RNA)
        }## if (is.null(RNAfilter))

        RNAnumBases <- RNAfilter[,2]
        rpkmdata <- myRPKM(ExprData=RawCounts, numBases=RNAnumBases)
        ##
        Norm.dat <- data.frame(Gene=Name.G, data.frame(rpkmdata))
        row.names(Norm.dat) <- Name.G
        ##
        Log2Trf <- TRUE
        YlabelNorm <- "log2 (rpkm normalized counts + 1)"
    }## if(Normalization == "rpkm")

    colnames(Norm.dat)[-1] <- SEsampleName

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## New SummarizeExperiment class object
    SEresNORM <- SEres
    SummarizedExperiment::assays(SEresNORM)[[2]] <- as.matrix(Norm.dat[, -1])
    names(SummarizedExperiment::assays(SEresNORM))[2] <- Normalization
    S4Vectors::metadata(SEresNORM)$SEidentification <- c("SEresNormalization")

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Boxplot
    resBXPLT <- DATAplotBoxplotSamples(SEres=SEresNORM,
                                       Log2.transformation=Log2Trf,
                                       Colored.By.Factors=Colored.By.Factors,
                                       Color.Group=Color.Group,
                                       Plot.genes=Plot.genes,
                                       y.label=YlabelNorm)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## SE object
    resLISTnorm <- list(normBoxplot=resBXPLT,
                        normMethod=Normalization)

    S4Vectors::metadata(SEresNORM)$Results[[1]][[1]] <- resLISTnorm

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Save of the boxplot graph
    csvNorm_filename <- paste0("NormalizedData_", Normalization, ".csv")
    pdfBox_filename <- paste0("BoxplotNormalization_", Normalization, ".pdf")

    if (!is.null(path.result)) {
        ## Save of the normalized data
        csvNorm <- file.path(path.result.new, csvNorm_filename)
        pdfBox <- file.path(path.result.new, pdfBox_filename)

        utils::write.table(Norm.dat, file=csvNorm, sep=";", row.names=FALSE)

        grDevices::pdf(file=pdfBox, width=11, height=8)
        print(resBXPLT)
        grDevices::dev.off()

    }## if(!is.null(path.result))

    if (isTRUE(Plot.Boxplot)) {
        print(resBXPLT)
    }## if(isTRUE(Plot.Boxplot))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output ## NormalizedBoxplot=resBXPLT ## NormalizedData=Norm.dat
    return(SEobj=SEresNORM)
}## DATAnormalization()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

vst2 <- function(x, y) {
    resVST2 <- DESeq2::varianceStabilizingTransformation(x, blind=y)
    return(resVST2=resVST2)
}## vst2()

myRPKM <- function(ExprData, numBases) {
    totalReads_inv <- 1/as.numeric(apply(data.frame(ExprData), 2, sum))
    RNAlength_inv <- 1/numBases

    Ngenes <- nrow(ExprData)
    Nsamples <- length(totalReads_inv)

    RPKMdata <- ExprData*(10^9)

    for (s in seq_len(Nsamples)) {
        RPKMdata[, s] <- RPKMdata[, s]*totalReads_inv[s]
    }## for (g in seq_len(Nsamples))

    for (g in seq_len(Ngenes)) {
        RPKMdata[g,] <- RPKMdata[g,]*RNAlength_inv[g]
    }## for (g in seq_len(Ngenes))

    return(RPKMdata)
}## myRPKM()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrDATAnormalization1 <- function(SEres,
                                  Normalization="vst",
                                  Blind.rlog.vst=FALSE) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## DATAprepSE
    if (!is(SEres, "SummarizedExperiment")) {
        stop("'SEres' mut be the results of the function 'DATAprepSE()'")
    } else {
        codeDEres <- S4Vectors::metadata(SEres)$SEidentification

        if (is.null(codeDEres)) {
            stop("'SEres' mut be the results of the function 'DATAprepSE()'")
        }## if (is.null(codeDEres))

        if (codeDEres != "SEstep") {
            stop("'SEres' mut be the results of the function 'DATAprepSE()'")
        }## if (codeDEres != "SEstep")
    }## if (!is(SEres, "SummarizedExperiment"))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Different normalization authorized
    if (!Normalization%in%c("rlog", "vst", "rle", "rpkm")) {
        stop("'Normalization' mut be 'vst', 'rlog', 'rle' or 'rpkm'")
    }## if(!Normalization%in%c("rlog","vst","rle"))

    ## Blind.rlog.vst
    if (!isTRUE(Blind.rlog.vst) & !isFALSE(Blind.rlog.vst)) {
        stop("'Blind.rlog.vst' must be TRUE or FALSE.")
    }## if (!isTRUE(Blind.rlog.vst) & !isFALSE(Blind.rlog.vst))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrDATAnormalization1()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrDATAnormalization2 <- function(Plot.Boxplot=TRUE,
                                  Colored.By.Factors=FALSE,
                                  Color.Group=NULL,
                                  Plot.genes=FALSE,
                                  path.result=NULL,
                                  Name.folder.norm=NULL) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check 1: TRUE FALSE
    if (!isTRUE(Plot.Boxplot) & !isFALSE(Plot.Boxplot)) {
        stop("'Plot.Boxplot' must be TRUE or FALSE.")
    }## if (!isTRUE(Plot.Boxplot) & !isFALSE(Plot.Boxplot))

    if (!isTRUE(Colored.By.Factors) & !isFALSE(Colored.By.Factors)) {
        stop("'Colored.By.Factors' must be TRUE or FALSE.")
    }## if (!isTRUE(Colored.By.Factors) & !isFALSE(Colored.By.Factors))

    if (!isTRUE(Plot.genes) & !isFALSE(Plot.genes)) {
        stop("'Plot.genes' must be TRUE or FALSE.")
    }## if (!isTRUE(Plot.genes) & !isFALSE(Plot.genes))

    if (!is.null(Color.Group)) {
        if (!is.data.frame(Color.Group)) {
            stop("'Color.Group' must be NULL or a data.frame.")
        }## if (!is.data.frame(Color.Group))
    }## if (!is.null(Color.Group))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check 2: path.result, Name.folder.norm
    if (!is.null(path.result)) {
        if (!is.character(path.result)) {
            stop("'path.result' must be NULL or a character.")
        }## if (!is.character(path.result))
    }## if (!is.null(path.result))

    if (!is.null(Name.folder.norm)) {
        if (!is.character(Name.folder.norm)) {
            stop("'Name.folder.norm' must be NULL or a character.")
        }## if (!is.character(Name.folder.norm))
    }## if (!is.null(Name.folder.norm))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrDATAnormalization2()

