#' @title Realization of the DE analysis (Main Function).
#'
#' @description The function realizes the DE analysis in three cases:
#' either samples belonging to different time measurements,
#' or samples belonging to different biological conditions,
#' or samples belonging to different time measurements
#' and different biological conditions.
#'
#' @details All results are built from the results of either our R function
#' [DATAprepSE()],
#' or our R function
#' [DATAnormalization()].
#'
#' @param SEres Results of either our R function
#' [DATAprepSE()],
#' or our R function
#' [DATAnormalization()].
#' @param pval.min Numeric value between 0 and 1. A gene will be considered as
#' differentially expressed (DE) between two biological conditions if
#' its Benjamini-Hochberg adjusted p-value
#' (see [stats::p.adjust()])
#' is below the threshold \code{pval.min}. Default value is 0.05.
#' @param pval.vect.t \code{NULL} or vector of dimension \eqn{T-1} filled with
#' numeric values between 0 and 1, with \eqn{T} the number of
#' time measurements.
#' A gene will be considered as differentially expressed (DE) between
#' the time ti and the reference time t0 if its Benjamini-Hochberg adjusted
#' p-value (see [stats::p.adjust()])
#' is below the i-th threshold of \code{pval.vect.t}.
#' If \code{NULL}, \code{pval.vect.t} will be vector of dimension \eqn{T-1}
#' filled with \code{pval.min}.
#' @param log.FC.min Non negative numeric value.
#' If the \eqn{log_2} fold change between biological conditions or times
#' has an absolute value below the threshold \code{log.FC.min},
#' then the gene is not selected even if is considered as DE.
#' Default value is 1.
#' If \code{log.FC.min=0}, all DE genes will be kept.
#' @param LRT.supp.info \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the algorithm realizes another statistical test in order to
#' detect if, among all biological conditions and/or times,
#' at least one has a different behavior than the others
#' (see the input \code{test} in [DESeq2::DESeq()]).
#' @param Plot.DE.graph \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, all graphs will be plotted.
#' Otherwise no graph will be plotted.
#' @param path.result Character or \code{NULL}. Path to save all results.
#' If \code{path.result} contains a sub folder entitled
#' "DEanalysis_\code{Name.folder.DE}" all results will be saved in
#' the sub folder "DEanalysis_\code{Name.folder.DE}".
#' Otherwise, a sub folder entitled "DEanalysis_\code{Name.folder.DE}"
#' will be created in \code{path.result} and all results will be saved in
#' "DEanalysis_\code{Name.folder.DE}".
#' If \code{NULL}, the results will not be saved in a folder.
#' \code{NULL} as default.
#' @param Name.folder.DE Character or \code{NULL}.
#' If \code{Name.folder.DE} is a character, the folder names which will
#' contain all results will be "DEanalysis_\code{Name.folder.DE}".
#' Otherwise, the folder name will be "DEanalysis".
#'
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEres} with the rle normalized count data (cf [DATAnormalization()])
#' automatically realized by
#' [DESeq2::DESeq()]
#' and saved in \code{assays(SEresNORM)$rle}, and with the following results
#' saved in the metadata \code{Results[[2]][[2]]} of \code{SEres},
#' depending on the experimental design.
#'
#' * If samples belong to different biological conditions only
#' (see [DEanalysisGroup()]),
#' the function returns
#'   * a data.frame (output \code{rowData(SEres)}) which contains
#'     * pvalues, log2 fold change and DE genes between each pairs of
#'     biological conditions.
#'     * a binary column (1 and 0) where 1 means the gene is DE between
#'     at least one pair of biological conditions.
#'     * \eqn{N_{bc}} binary columns,
#'     where \eqn{N_{bc}} is the number of biological conditions,
#'     which gives the specific genes for each biological condition.
#'     A '1' in one of these columns means the gene is specific to
#'     the biological condition associated to the given column. 0 otherwise.
#'     A gene is called specific to a given biological condition BC1,
#'     if the gene is DE between BC1 and any other biological conditions,
#'     but not DE between any pair of other biological conditions.
#'     * \eqn{N_{bc}} columns filled with -1, 0 and 1, one per biological
#'     condition.
#'     A '1' in one of these columns means the gene is up-regulated
#'     (or over-expressed) for the biological condition associated to the
#'     given column.
#'     A gene is called up-regulated for a given biological condition BC1 if
#'     the gene is specific to the biological condition BC1 and expressions
#'     in BC1 are higher than in the other biological conditions.
#'     A '-1' in one of these columns means the gene is down-regulated
#'     (or under-expressed) for the biological condition associated to the
#'     given column.
#'     A gene is called down-regulated for a given biological condition BC1 if
#'     the gene is specific to the biological condition BC1 and expressions
#'     in BC1 are lower than in the other biological conditions.
#'     A '0' in one of these columns means the gene is not specific to the
#'     biological condition associated to the given column.
#'   * an UpSet plot (Venn diagram displayed as a barplot) which gives the
#'   number of genes for each possible intersection
#'   (see [DEplotVennBarplotGroup()]).
#'   We consider that a set of pairs of biological conditions forms an
#'   intersection if there is at least one gene which is DE for each of these
#'   pairs of biological conditions, but not for the others.
#'   * a barplot which gives the number of genes categorized as "Upregulated"
#'   and "DownRugulated", per biological condition
#'   (see [DEplotBarplot()]).
#'   * a barplot which gives the number of genes categorized as "Upregulated",
#'   "DownRugulated" and "Other", per biological condition
#'   (see [DEplotBarplot()]).
#'   A gene is categorized as 'Other', for a given biological condition,
#'   if the gene is not specific to the given biological condition.
#'   So this barplot, only plotted when there are strictly more than two
#'   biological conditions, is similar to the previous barplot but with
#'   the category "Other".
#'   * a list (output \code{List.Glossary}) containing the glossary of
#'   the column names of \code{DE.results}.
#'   * a list (output \code{Summary.Inputs}) containing a summary of sample
#'   information and inputs of
#'   [DEanalysisGlobal()].
#'
#' * If data belong to different time points only
#' (see [DEanalysisTime()]),
#' the function returns
#'   * a data.frame (output \code{rowData(SEres)}) which contains
#'     * gene names
#'     * pvalues, log2 fold change and DE genes between each time ti versus
#'     the reference time t0.
#'     * a binary column (1 and 0) where 1 means the gene is DE at at least
#'     between one time ti versus the reference time t0.
#'     * a column where each element is succession of 0 and 1.
#'     The positions of '1' indicate the set of times ti such that the gene
#'     is DE between ti and the reference time t0.
#'   * an alluvial graph of differentially expressed (DE) genes
#'   (see [DEplotAlluvial()])
#'   * a graph showing the number of DE genes as a function of time for
#'   each temporal group
#'   (see [DEplotAlluvial()]).
#'   By temporal group, we mean the sets of genes which are first DE at
#'   the same time.
#'   * a barplot which gives the number of DE genes per time
#'   (see [DEplotBarplotTime()])
#'   * an UpSet plot which gives the number of genes per temporal pattern
#'   (see [DEplotVennBarplotTime()]).
#'   By temporal pattern, we mean the set of times ti such that the gene is
#'   DE between ti and the reference time t0.
#'   * a similar UpSet plot where each bar is split in different colors
#'   corresponding to all possible numbers of DE times where genes are over
#'   expressed in a given temporal pattern.
#'   * a list (output \code{List.Glossary}) containing the glossary of
#'   the column names of \code{DE.results}.
#'   * a list (output \code{Summary.Inputs}) containing a summary of sample
#'   information and inputs of
#'   [DEanalysisGlobal()].
#'
#' * If data belong to different time points and different biological
#' conditions
#' (see [DEanalysisTimeAndGroup()]),
#' the function returns
#'   * a data.frame (output \code{rowData(SEres)}) which contains
#'     * gene names
#'     * Results from the temporal statistical analysis
#'       * pvalues, log2 fold change and DE genes between each pairs of
#'       biological conditions for each fixed time.
#'       * \eqn{N_{bc}} binary columns (0 and 1), one per biological condition
#'       (with \eqn{N_{bc}} the number of biological conditions).
#'       A 1 in one of these two columns means the gene is DE at least between
#'       one time ti versus the reference time t0, for the biological condition
#'       associated to the given column.
#'       * \eqn{N_{bc}} columns, one per biological condition, where each
#'       element is succession of 0 and 1. The positions of 1 in one of these
#'       two columns, indicate the set of times ti such that the gene is DE
#'       between ti and the reference time t0, for the biological condition
#'       associated to the given column.
#'     * Results from the statistical analysis by biological condition
#'       * pvalues, log2 fold change and DE genes between each time ti
#'       and the reference time t0 for each biological condition.
#'       * \eqn{T} binary columns (0 and 1), one per time
#'       (with \eqn{T} the number of time measurements).
#'       A 1 in one of these columns, means the gene is DE between at least
#'       one pair of biological conditions, for the fixed time associated
#'       to the given column.
#'       * \eqn{T \times N_{bc}} binary columns, which give the genes specific
#'       for each biological condition at each time ti.
#'       A 1 in one of these columns means the gene is specific to the
#'       biological condition at a fixed time associated to the given column.
#'       0 otherwise. A gene is called specific to a given biological condition
#'       BC1 at a time ti, if the gene is DE between BC1 and any other
#'       biological conditions at time ti, but not DE between any pair of
#'       other biological conditions at time ti.
#'       * \eqn{T \times N_{bc}} columns filled with -1, 0 and 1.
#'       A 1 in one of these columns means the gene is up-regulated
#'       (or over-expressed) for the biological condition at a fixed time
#'       associated to the given column. A gene is called up-regulated for a
#'       given biological condition BC1 at time ti if the gene is specific to
#'       the biological condition BC1 at time ti and expressions in BC1 at time
#'       ti are higher than in the other biological conditions at time ti.
#'       A -1 in one of these columns means the gene is down-regulated
#'       (or under-expressed) for the biological condition at a fixed time
#'       associated to the given column. A gene is called down-regulated for a
#'       given biological condition at a time ti BC1 if the gene is specific to
#'       the biological condition BC1 at time ti and expressions in BC1 at time
#'       ti are lower than in the other biological conditions at time ti.
#'       A 0 in one of these columns means the gene is not specific to the
#'       biological condition at a fixed time associated to the given column.
#'       * \eqn{N_{bc}} binary columns (0 and 1). A 1 in one of these columns
#'       means the gene is specific at at least one time ti, for the biological
#'       condition associated to the given column. 0 otherwise.
#'     * Results from the combination of temporal and biological statistical
#'     analysis
#'       * \eqn{T \times N_{bc}} binary columns, which give the signatures
#'       genes for each biological condition at each time ti.
#'       A 1 in one of these columns means the gene is signature gene to the
#'       biological condition at a fixed time associated to the given column.
#'       0 otherwise. A gene is called signature of a biological condition
#'       BC1 at a given time ti, if the gene is specific to the biological
#'       condition BC1 at time ti and DE between ti versus the reference time
#'       t0 for the biological condition BC1.
#'       * \eqn{N_{bc}} binary columns (0 and 1). A 1 in one of these columns
#'       means the gene is signature at at least one time ti, for the
#'       biological condition associated to the given column. 0 otherwise.
#'   * the following plots from the temporal statistical analysis
#'     * a barplot which gives the number of DE genes between ti and the
#'     reference time t0, for each time ti (except the reference time t0) and
#'     biological condition
#'     (see [DEplotBarplotFacetGrid()]).
#'     * \eqn{N_{bc}} alluvial graphs of DE genes
#'     (see [DEplotAlluvial()]),
#'     one per biological condition.
#'     * \eqn{N_{bc}} graphs showing the number of DE genes as a function of
#'     time for each temporal group, one per biological condition.
#'     By temporal group, we mean the sets of genes which are first DE at the
#'     same time.
#'     * \eqn{2\times N_{bc}} UpSet plot showing the number of DE genes
#'     belonging to each DE temporal pattern, for each biological condition.
#'     By temporal pattern, we mean the set of times ti such that the gene is
#'     DE between ti and the reference time t0
#'     (see [DEplotVennBarplotTime()]).
#'     * an alluvial graph for DE genes which are DE at least one time for
#'     each group.
#'   * the following plots from the statistical analysis by biological
#'   condition
#'     * a barplot which gives the number of specific DE genes for each
#'     biological condition and time
#'     (see [DEplotBarplotFacetGrid()]).
#'     * \eqn{N_{bc}(N_{bc}-1)/2} UpSet plot which give the number of genes
#'     for each possible intersection (set of pairs of biological conditions),
#'     one per time
#'     (see [DEplotVennBarplotGroup()]).
#'     * an alluvial graph of genes which are specific at least one time
#'     (see [DEplotAlluvial()]).
#'   * the following plots from the combination of temporal and biological
#'   statistical analysis
#'     * a barplot which gives the number of signature genes for each
#'     biological condition and time
#'     (see [DEplotBarplotFacetGrid()]).
#'     * a barplot showing the number of genes which are DE at at least one
#'     time, specific at at least one time and signature at at least one time,
#'     for each biological condition.
#'     * an alluvial graph of genes which are signature at least one time
#'     (see [DEplotAlluvial()]).
#'
#'
#' @importFrom DESeq2 DESeq counts
#' @importFrom SummarizedExperiment assays colnames rownames colData rowData
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' data(RawCounts_Antoszewski2022_MOUSEsub500)
#' ## No time points. We take only two groups for the speed of the example
#' RawCounts_T1Wt <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200),
#'                                                         seq_len(7)]
#' ##------------------------------------------------------------------------##
#' ## Preprocessing
#' resDATAprepSE <- DATAprepSE(RawCounts=RawCounts_T1Wt,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=NULL,
#'                             Individual.position=2)
#' ##------------------------------------------------------------------------##
#' ## DE analysis
#' resDE <- DEanalysisGlobal(SEres=resDATAprepSE,
#'                           pval.min=0.05,
#'                           pval.vect.t=NULL,
#'                           log.FC.min=1,
#'                           LRT.supp.info=FALSE,
#'                           Plot.DE.graph=TRUE,
#'                           path.result=NULL,
#'                           Name.folder.DE=NULL)

DEanalysisGlobal <- function(SEres,
                             pval.min=0.05,
                             pval.vect.t=NULL,
                             log.FC.min=1,
                             LRT.supp.info=FALSE,
                             Plot.DE.graph=TRUE,
                             path.result=NULL,
                             Name.folder.DE=NULL){
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check 1
    ## DATAprepSE
    Err_SE <- paste0("'SEres' mut be the results of either the function ",
                     "'DATAprepSE()' or 'DATAnormalization()'.")

    ## DATAprepSE
    if (!is(SEres, "SummarizedExperiment")) {
        stop(Err_SE)
    } else {
        codeDEres <- S4Vectors::metadata(SEres)$SEidentification

        if (is.null(codeDEres)) {
            stop(Err_SE)
        }## if (is.null(codeDEres))

        if (!codeDEres%in%c("SEstep", "SEresNormalization")) {
            stop(Err_SE)
        }## if (!codeDEres%in%c("SEstep", "SEresNormalization"))
    }## if (!is(SEres, "SummarizedExperiment"))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Creation folder if no existence
    if (is.null(Name.folder.DE)) {
        Name.folder.DE.ini <- NULL
        Name.folder.DE <- ""
        SubFolder.name <- "2_SupervisedAnalysis"
    } else {
        Name.folder.DE.ini <- Name.folder.DE
        Name.folder.DE <- paste0("_", Name.folder.DE.ini)
        SubFolder.name <- paste0("2_SupervisedAnalysis", Name.folder.DE)
    }## if (is.null(Name.folder.DE))

    name.folder.result1 <- paste0("2-1_RLEnormalizedDATA", Name.folder.DE)

    if (!is.null(path.result)) {
        if (!SubFolder.name%in%dir(path=path.result)){
            print("Folder creation")
            dir.create(path=file.path(path.result, SubFolder.name))
        }## if(SubFolder.name%in%dir(path=path.result) == FALSE)
        path.result.f <- file.path(path.result, SubFolder.name)
    } else {
        path.result.f <- NULL
    }## if(is.null(path.result)==FALSE)

    ## Folder for RLE normalized count data
    if (!is.null(path.result)) {
        if (!name.folder.result1%in%dir(path=path.result.f)) {
            dir.create(path=file.path(path.result.f, name.folder.result1))
        }## if (!name.folder.result1%in%dir(path = path.result.f))
        path.result.new1 <- file.path(path.result.f, name.folder.result1)
    } else {
        path.result.new1 <- NULL
    }## if(is.null(path.result)==FALSE)

    ## Folder for DE results csv
    if (!is.null(path.result)) {
        name.folder.result3 <- paste0("2-3_CSV_file_DEanalysis", Name.folder.DE)
        if (!name.folder.result3%in%dir(path=path.result.f)) {
            dir.create(path=file.path(path.result.f, name.folder.result3))
        }## if(name.folder.result2%in%dir(path = path.result.f)==FALSE)

        path.result.new3 <- file.path(path.result.f, name.folder.result3)
    } else {
        path.result.new3 <- NULL
    }## if(is.null(path.result)==FALSE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Preprocessing
    print("Preprocessing")
    DESeq2.obj <- S4Vectors::metadata(SEres)$DESeq2obj$DESeq2preproceesing

    RawCounts <- SummarizedExperiment::assays(SEres)[[1]]
    SEsampleName <- SummarizedExperiment::colnames(SEres)
    Name.G <- as.character(SummarizedExperiment::rownames(SEres))
    cDat <- data.frame(SummarizedExperiment::colData(DESeq2.obj))
    FactorBoxplt <- data.frame(SummarizedExperiment::colData(SEres))

    colFCTRS <- S4Vectors::metadata(SEres)$colDataINFO
    colFCTRS <- as.numeric(colFCTRS$colINFOfactors)

    FactorInfo.f <- FactorBoxplt[, colFCTRS]
    names(FactorInfo.f)[ncol(FactorInfo.f)] <- "Samples"

    ##-----------------------------------------------------------------------##
    if (c("Group")%in%colnames(cDat)) {
        Vector.group <- as.character(cDat$Group)
        LvlsGROUP <- levels(as.factor(Vector.group))
    } else {
        Vector.group <- NULL
    }## if (c("Group")%in%colnames(cDat))

    ##-----------------------------------------------------------------------##
    if (c("Time")%in%colnames(cDat)) {
        Vector.time <- as.character(cDat$Time)
        Levels.time <- levels(as.factor(Vector.time))
        Nb.time <- length(Levels.time)

        Tnumeric <- gsub("t", "", gsub("T", "", FactorInfo.f$Time, fixed=TRUE),
                         fixed=TRUE)
        TlevNumeric <- gsub("t", "", gsub("T", "", Levels.time, fixed=TRUE),
                            fixed=TRUE)
        FactorInfo.f$Time <- paste("t", Tnumeric, sep="")
        TlevNumeric <- paste("t", Tnumeric, sep="")
    } else {
        Vector.time <- NULL
    }## if (c("Time")%in%colnames(cDat))

    ##-----------------------------------------------------------------------##
    if (!is.null(Vector.time)) {
        if (!is.null(pval.vect.t)) {
            if (length(pval.vect.t) > Nb.time - 1) {
                pval.vect.t <- pval.vect.t[seq_len(Nb.time-1)]
            }## if(length(pval.vect.t)>Nb.time-1)

            if (length(pval.vect.t) < Nb.time - 1) {
                repPvalmin <- Nb.time - length(pval.vect.t) - 1
                pval.vect.t <- c(pval.vect.t, rep(pval.min, times=repPvalmin))
            }## if(length(pval.vect.t)<Nb.time-1)
        } else {
            pval.vect.t <- rep(pval.min, times=Nb.time-1)
        }## if(is.null(pval.vect.t)==FALSE)
    }## if(is.null(Vector.time)==FALSE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Names folder for DE graph
    if (is.null(Vector.time) & !is.null(Vector.group)) {
        name.folder.result2 <- paste0("2-2_Group_DEanalysis", Name.folder.DE)
    }## if (is.null(Vector.time) & !is.null(Vector.group))

    if (!is.null(Vector.time) & is.null(Vector.group)) {
        name.folder.result2 <- paste0("2-2_Temporal_DEanalysis", Name.folder.DE)
    }## if (!is.null(Vector.time) & is.null(Vector.group))

    if (!is.null(Vector.time) & !is.null(Vector.group)) {
        name.folder.result2 <- paste0("2-2_DEanalysis" , Name.folder.DE)
    }## if (!is.null(Vector.time) & !is.null(Vector.group))

    ## Folder for DE graphs
    if(!is.null(path.result)) {
        if (!name.folder.result2%in%dir(path=path.result.f)) {
            dir.create(path=file.path(path.result.f, name.folder.result2))
        }## if (!name.folder.result2%in%dir(path=path.result.f))
        path.result.new2 <- file.path(path.result.f, name.folder.result2)
    } else {
        path.result.new2 <- NULL
    }## if(!is.null(path.result)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Differential expression
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    print("Differential expression step with DESeq2::DESeq()")

    if (isTRUE(LRT.supp.info)) {
        dds.norm.diff <- DESeq2::DESeq(DESeq2.obj, betaPrior=FALSE,
                                       test="LRT", reduced=~1)
    } else {
        dds.norm.diff <- DESeq2::DESeq(DESeq2.obj, betaPrior=FALSE,
                                       test="Wald")
    }## if (isTRUE(LRT.supp.info))

    ##-----------------------------------------------------------------------##
    ScaledData <- round(DESeq2::counts(dds.norm.diff, normalized=TRUE),
                        digits=3)
    RLEdata <- data.frame(Gene=Name.G, ScaledData)

    ##-----------------------------------------------------------------------##
    ## SE object
    SEresDE <- SEres
    S4Vectors::metadata(SEresDE)$SEidentification <- c("SEresNormalization")
    SummarizedExperiment::assays(SEresDE)$rle <- ScaledData
    S4Vectors::metadata(SEresDE)$DESeq2obj$DESeq2results <- dds.norm.diff

    ##-----------------------------------------------------------------------##
    yRLEboxplot <- "log2 (rle normalized counts + 1)"
    rleBXPLT <- DATAplotBoxplotSamples(SEres=SEresDE,
                                       Log2.transformation=TRUE,
                                       Colored.By.Factors=FALSE,
                                       Color.Group=NULL,
                                       Plot.genes=FALSE,
                                       y.label=yRLEboxplot)

    rleLISTnorm <- list(normBoxplot=rleBXPLT, normMethod="rle")
    S4Vectors::metadata(SEresDE)$Results[[2]][[1]] <- rleLISTnorm

    if (!is.null(path.result)) {
        RLEname <- paste0("NormalizedData_rle",  Name.folder.DE, ".csv")
        utils::write.table(RLEdata, file=file.path(path.result.new1, RLEname),
                           sep=";", row.names = FALSE)

        BXPLTname <- paste0("BoxplotNormalization_rle", ".pdf")
        grDevices::pdf(file=file.path(path.result.new1, BXPLTname),
                       width=11, height=8)## width = 8, height = 11
        print(rleBXPLT)
        grDevices::dev.off()
    }## if (!is.null(path.result))

    if (isTRUE(Plot.DE.graph)) {
        print(rleBXPLT)
    }## if (isTRUE(Plot.DE.graph))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Results from Differential expression step with DESeq2::DESeq()
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Case 1 analysis DE : Biological conditions only
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##

    if (is.null(Vector.time) & !is.null(Vector.group)) {
        print("Case 1 analysis : Biological conditions only")

        resDEbioncond <- DEanalysisGroup(DESeq.result=dds.norm.diff,
                                         LRT.supp.info=LRT.supp.info,
                                         log.FC.min=log.FC.min,
                                         pval.min=pval.min,
                                         Plot.DE.graph=Plot.DE.graph,
                                         path.result=path.result.new2,
                                         SubFile.name=SubFolder.name)

        resDEgroup <- S4Vectors::metadata(resDEbioncond)$DEresultsGroup
        resDEsummary <- data.frame(resDEgroup$DEsummary)

        ##-------------------------------------------------------------------##
        ## SE final
        SumInfo <- list(ExprCond=c("Group"),
                        FactorsInfo=FactorInfo.f,
                        GroupLevels=LvlsGROUP,
                        logFCmin=log.FC.min,
                        pvalGroup=pval.min)

        resGlossary <- Glossary(path.result.new3, Case=1)

        listDEresults <- append(resDEgroup[-c(1, 2, 3)],
                                list(Glossary=resGlossary),
                                after=0)
    }## if (is.null(Vector.time) & !is.null(Vector.group))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Case 2 analysis DE : Time only
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##

    if (!is.null(Vector.time) & is.null(Vector.group)) {
        print("Case 2 analysis : Time only")

        resDEtimeSE <- DEanalysisTime(DESeq.result=dds.norm.diff,
                                      LRT.supp.info=LRT.supp.info,
                                      log.FC.min=log.FC.min,
                                      pval.min=pval.min,
                                      pval.vect.t=pval.vect.t,
                                      Plot.DE.graph=Plot.DE.graph,
                                      path.result=path.result.new2,
                                      SubFile.name=SubFolder.name)

        resDEtime <- S4Vectors::metadata(resDEtimeSE)$DEresultsTime
        resDEsummary <- data.frame(resDEtime$DEsummary)

        ##-------------------------------------------------------------------##
        ## SE final
        SumInfo <- list(ExprCond=c("Time"),
                        FactorsInfo=FactorInfo.f,
                        TimeLevels=levels(factor(TlevNumeric)),
                        logFCmin=log.FC.min,
                        pvalsTime=pval.vect.t)

        resGlossary <- Glossary(path.result.new3, Case=2)

        listDEresults <- append(resDEtime[-c(1,2)], list(Glossary=resGlossary),
                                after=0)
    }## if(!is.null(Vector.time) & is.null(Vector.group))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Case 3 analysis DE : Time and Biological conditions
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (!is.null(Vector.time) & !is.null(Vector.group)) {
        print("Case 3 analysis : Biological conditions and Times.")

        resDEtg <- DEanalysisTimeAndGroup(DESeq.result=dds.norm.diff,
                                          LRT.supp.info=LRT.supp.info,
                                          log.FC.min=log.FC.min,
                                          pval.min=pval.min,
                                          pval.vect.t=pval.vect.t,
                                          Plot.DE.graph=Plot.DE.graph,
                                          path.result=path.result.new2,
                                          SubFile.name=SubFolder.name)

        DEresTG <- S4Vectors::metadata(resDEtg)$DEresultsTimeGroup
        resDEsummary <- data.frame(DEresTG$DEsummary)

        ##-------------------------------------------------------------------##
        ## SE final
        SumInfo <- list(ExprCond=c("Time","Group"),
                        FactorsInfo=FactorInfo.f,
                        TimeLevels=levels(factor(TlevNumeric)),
                        GroupLevels=LvlsGROUP,
                        logFCmin=log.FC.min,
                        pvalsTime=pval.vect.t,
                        pvalGroup=pval.min)

        resGlossary <- Glossary(path.result.new3, Case=3)

        listDEresults <- append(DEresTG[-c(1, 2)], list(Glossary=resGlossary),
                                after=0)
    }## if (!is.null(Vector.time) & !is.null(Vector.group))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## SE common information
    SummarizedExperiment::rowData(SEresDE) <- resDEsummary
    S4Vectors::metadata(SEresDE)$Results[[2]][[2]] <- listDEresults

    listPATHname <- list(Path.result=path.result.f,
                         Folder.result=Name.folder.DE.ini)
    S4Vectors::metadata(SEresDE)$DESeq2obj$Summary.Inputs <- SumInfo
    S4Vectors::metadata(SEresDE)$DESeq2obj$pathNAME <- listPATHname
    S4Vectors::metadata(SEresDE)$DESeq2obj$SEidentification <- "SEresultsDE"
    ## S4Vectors::metadata(SEresDE)$DESeq2obj$List.Glossary <- resGlossary
    ## S4Vectors::metadata(SEresDE)$DESeq2obj$RLEdata <- RLEdata

    ## Tables which contain all results
    if (!is.null(path.result)) {
        DEsum_file <- paste0("ALLresults_DEanalysis", Name.folder.DE, ".csv")
        utils::write.table(resDEsummary,
                           file=file.path(path.result.new3, DEsum_file),
                           sep=";", row.names=FALSE)
    }## if (!is.null(path.result))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    return(SEresDE)
}## DEanalysisGlobal()


##---------------------------------------------------------------------------##
##===========================================================================##
##===========================================================================##
##---------------------------------------------------------------------------##

Glossary <- function(path.result, Case){
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## List Glossary general information
    ListPrepro <- vector(mode="list", length=8)
    paste100rep <- paste(rep("=", times=100), collapse="")
    paste15rep <- paste(rep("-", times=15), collapse="")
    paste15rep <- paste0("=", paste15rep, "=", collapse="")

    ListPrepro[[2]] <- ListPrepro[[1]] <- paste100rep
    ListPrepro[[5]] <- ListPrepro[[4]] <- paste100rep
    ListPrepro[[8]] <- paste100rep

    ListPrepro[[3]] <- paste(paste15rep, "GLOSSARY of column names in the",
                             "cvs file 'ALLresults_DEanalysis'",
                             paste15rep)
    ListPrepro[[6]] <- c("\n")
    ListPrepro[[7]] <- paste("The cvs file 'ALLresults_DEanalysis' gathers",
                             "all the results of MultiRNAflow analysis.",
                             "The goal of this glossary is to clarify",
                             "the meaning of the column names in the file.",
                             "The first column gives the names of all genes.",
                             "In the list below, we specify the meaning of",
                             "each entry in the corresponding column.\n")

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Glossary when samples depends on biological condition only
    if (Case == 1) {
        ListGlossary <- vector(mode="list", length=6)
        names(ListGlossary) <- c("Log2FoldChange.Group2.versus.Group1",
                                 "Pvalue.adjusted.Group2.versus.Group1",
                                 "DE.Group2.versus.Group1",
                                 "DE.1pair.of.Group.minimum",
                                 "Specific.genes_Group1",
                                 "OverUnder.regulated.genes_Group1")

        ListGlossary[[1]] <- paste("Log2 fold change between",
                                   "the biological condition Group2 and",
                                   "the biological condition Group1.")
        ListGlossary[[2]] <- paste("Adjusted p-value between",
                                   "the biological condition Group2 and",
                                   "the biological condition Group1.")
        ListGlossary[[3]] <- paste("Binary number (0 or 1) where 1 means",
                                   "that the gene is DE and 0 means",
                                   "that it is not DE between",
                                   "the biological condition Group2 and",
                                   "the biological condition Group1.",
                                   "The value 1 is given if both",
                                   "'abs(Log2FoldChange.Group2.versus.Group1)",
                                   ">log.FC.min' and",
                                   "'Pvalue.adjusted.Group2.versus.Group1",
                                   "<pval.min', where",
                                   "'log.FC.min' and 'pval.min' are inputs",
                                   "of the function DEanalysisGlobal()")
        ListGlossary[[4]] <- paste("Binary number (0 or 1) where 1 means that",
                                   "the gene is DE between at least one pair",
                                   "of biological conditions.")
        ListGlossary[[5]] <- paste("Binary number (0 or 1)",
                                   "where 1 means that the gene is specific",
                                   "for the biological condition Group1.",
                                   "This means that the gene is DE between",
                                   "Group1 and any other biological conditions",
                                   "but not DE between any pairs of other",
                                   "biological conditions.")
        ListGlossary[[6]] <- paste("Number (-1, 0 or 1).",
                                   "1 means the gene is specific",
                                   "('Specific.genes_Group1=1') and the",
                                   "gene is up-regulated (or over expressed)",
                                   "in Group1 versus",
                                   "the other biological conditions.",
                                   "-1' means the gene is specific",
                                   "('Specific.genes_Group1=1') and",
                                   "the gene is down-regulated",
                                   "(or under expressed) in Group1",
                                   "versus the other biological conditions.",
                                   "0 otherwise.")
    }## if(Case==1)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Glossary when samples depends on time only
    if (Case == 2) {
        ListGlossary <- vector(mode="list", length=5)
        names(ListGlossary) <- c("Log2FoldChange.ti.versus.t0",
                                 "Pvalue.adjusted.ti.versus.t0",
                                 "DE.ti.versus.t0",
                                 "DE.Temporal.Pattern",
                                 "DE.1time.minimum")

        ListGlossary[[1]] <- paste("Log2 fold change between the time ti and",
                                   "the reference time t0.")
        ListGlossary[[2]] <- paste("Adjusted pvalue between the time ti and",
                                   "the reference time t0.")
        ListGlossary[[3]] <- paste0("Binary number (0 or 1) ",
                                    "where 1 means that the gene is DE ",
                                    "and 0 means it is not DE between ",
                                    "the time ti and the reference time t0. ",
                                    "The value 1 is given if both ",
                                    "'abs(Log2FoldChange.ti.versus.t0)>",
                                    "log.FC.min' and 'Pvalue.adjusted.",
                                    "ti.versus.t0<pval.vect.t[i]'.")
        ListGlossary[[4]] <- paste("Vector of 0 and 1 corresponding to times",
                                   "t1 to tn. The values 1 correspond to the",
                                   "times ti such that the gene is DE between",
                                   "ti and the reference time t0.")
        ListGlossary[[5]] <- paste("Binary number (0 or 1) where 1 means that",
                                   "the gene is DE at least between one time",
                                   " ti versus the reference time t0.")
    }## if(Case==2)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Glossary when samples depends on time and biological condition
    if (Case == 3) {
        ListGlossary <- vector(mode="list", length=14)
        names(ListGlossary) <- c("Log2FoldChange.ti.versus.t0_Group1",
                                 "Pvalue.adjusted.ti.versus.t0_Group1",
                                 "DE.ti.versus.t0_Group1",
                                 "DE.Temporal.Pattern_Group1",
                                 "DE.1time.minimum_Group1",
                                 "Log2FoldChange.Group2.versus.Group1_Time.ti",
                                 "Pvalue.adjusted.Group2.versus.Group1_Time.ti",
                                 "DE.Group2.versus.Group1_Time.ti",
                                 "DE.1pair.of.Group.minimum_Time.ti",
                                 "Specific.genes_Group1_Time.ti",
                                 "OverUnder.regulated.genes_Group1_Time.ti",
                                 "Specific.genes_Group1_1t.minimum",
                                 "Signature.genes_Group.Group1_Time.ti",
                                 "Signature.genes_Group.Group1_1time.minimum")

        ## Time for each group
        ListGlossary[[1]] <- paste("Log2 fold change between the time ti and",
                                   "the reference time t0,",
                                   "for the biological condition Group1.")
        ListGlossary[[2]] <- paste("Adjusted pvalue between the time ti and ",
                                   "the reference time t0, ",
                                   "for the biological condition Group1.",
                                   sep="")
        ListGlossary[[3]] <- paste("Binary number (0 or 1) where 1 means that ",
                                   "the gene is DE and 0 means it is not DE, ",
                                   "for the biological condition Group1. ",
                                   "The value 1 is given if both ",
                                   "'abs(Log2FoldChange.ti.versus.t0_Group1)>",
                                   "log.FC.min' and ",
                                   "'Pvalue.adjusted.ti.versus.t0_Group1<",
                                   "pval.vect.t[i]', for the group Group1.",
                                   sep="")
        ListGlossary[[4]] <- paste("Vector of 0 and 1 corresponding to times ",
                                   "t1 to tn. The values 1 correspond to the ",
                                   "times ti such that the gene is DE between ",
                                   "ti and the reference time t0, ",
                                   "for the biological condition Group1.",
                                   sep="")
        ListGlossary[[5]] <- paste("Binary number (0 or 1) ",
                                   "where 1 means that the ",
                                   "gene is DE at least between one time ti ",
                                   "versus the reference time t0, ",
                                   "for the group Group1.", sep="")

        ## Group for each time
        ListGlossary[[6]] <- paste("Log2 fold change between ",
                                   "the biological condition Group2 and ",
                                   "the biological condition Group1, ",
                                   "at time ti.", sep="")
        ListGlossary[[7]] <- paste("Adjusted p-value between ",
                                   "the biological condition Group2 and ",
                                   "the biological condition Group1, ",
                                   "at time ti.", sep="")
        ListGlossary[[8]] <- paste("Binary number (0 or 1) where ",
                                   "1 means that the gene is DE and 0 means ",
                                   "that it is not DE between ",
                                   "the biological condition Group2 and the ",
                                   "biological condition Group1, at time ti. ",
                                   "The value 1 is given if both 'abs(Pvalue.",
                                   "adjusted.Group2.versus.Group1_Time.ti)",
                                   ">log.FC.min' and 'Pvalue.adjusted.",
                                   "Group2.versus.Group1_Time.ti < pval.min', ",
                                   "where 'log.FC.min' and 'pval.min' are ",
                                   "inputs of the function DEanalysisGlobal()",
                                   sep="")
        ListGlossary[[9]] <- paste("Binary number (0 or 1) where 1 means that ",
                                   "the gene is DE between at least one pair ",
                                   "of biological conditions, at time ti.",
                                   sep="")
        ListGlossary[[10]] <- paste("Binary number (0 or 1)",
                                    "where 1 means that the gene",
                                    "is specific for the biological condition",
                                    "Group1, at time ti. This means that the",
                                    "gene is DE between Group1 and any other",
                                    "biological conditions but not DE between",
                                    "any pairs of other biological conditions,",
                                    "at time ti.")
        ListGlossary[[11]] <- paste("Number (-1, 0 or 1). ",
                                    "1 means the gene is specific ",
                                    "('Specific.genes_Group1_Time.ti=1') ",
                                    "at time ti and the gene is up-regulated ",
                                    "(or over expressed) in Group1 versus ",
                                    "the other biological conditions, ",
                                    "for the time ti. ",
                                    "-1' means the gene is specific ",
                                    "('Specific.genes_Group1_Time.ti=1') ",
                                    "at time ti and the gene is down-regulated",
                                    " (or under expressed) in Group1 versus ",
                                    "the other biological conditions, ",
                                    "for the time ti. ",
                                    "0 otherwise.", sep="")
        ListGlossary[[12]] <- paste("Binary number (0 or 1). ",
                                    "1 means that the gene is specific ",
                                    "for the biological condition Group1, ",
                                    "at at least one time ti. ",
                                    "0 otherwise.", sep="")

        ## Time and Group
        ListGlossary[[13]] <- paste("Binary number (0 or 1). ",
                                    "1 means that the gene ",
                                    "is a signature gene for the biological ",
                                    "condition Group1 at time ti. ",
                                    "This means that the gene is specific ",
                                    "for the biological condition Group1 ",
                                    "at time ti and the gene is ",
                                    "DE between the time ti and the reference ",
                                    "time t0 for the group Group1.", sep="")
        ListGlossary[[14]] <- paste("Binary number (0 or 1). ",
                                    "1 means that the gene ",
                                    "is a signature gene for the biological ",
                                    "condition Group1 at at least one time ti.",
                                    " 0 otherwise.", sep="")
    }## if(Case==3)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Specific case 3 and ending list
    if (Case == 3) {
        ListCase3 <- vector(mode="list", length=7)
        ListCase3[[1]] <- paste("==-----= ", "Temporal statistical analysis",
                                "\n", collapse="")
        ListCase3[[2]] <- c("\n")
        ListCase3[[3]] <- paste100rep
        ListCase3[[4]] <- paste("==-----= ",
                                "Statistical analysis by biological condition",
                                "\n", collapse="")
        ListCase3[[5]] <- c("\n")
        ListCase3[[6]] <- paste100rep
        ListCase3[[7]] <- paste("==-----= Combination of",
                                "temporal and condition analysis",
                                "\n", collapse="")
    }## if (Case == 3)


    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (!is.null(path.result)) {
        fileRdme <- paste(path.result,"/", "Glossary.txt", sep="")
        cat(ListPrepro[[1]], sep="", file=fileRdme, append=FALSE)

        for (l in seq(2, 8, 1)) {
            cat(ListPrepro[[l]], sep="", file=fileRdme, append=TRUE)
        }## for (l in seq(2, 8, 1))

        for(i in seq_len(length(ListGlossary))){
            if (i%in%c(1, 6, 13) & Case == 3) {
                if (i == 1) {
                    cat(ListCase3[[1]], sep="", file=fileRdme, append=TRUE)
                }## if(i==1)

                if (i == 6) {
                    cat(ListCase3[[2]], sep="", file=fileRdme, append=TRUE)
                    cat(ListCase3[[3]], sep="", file=fileRdme, append=TRUE)
                    cat(ListCase3[[4]], sep="", file=fileRdme, append=TRUE)
                }## if(i == 6)

                if (i == 13) {
                    cat(ListCase3[[5]], sep="", file=fileRdme, append=TRUE)
                    cat(ListCase3[[6]], sep="", file=fileRdme, append=TRUE)
                    cat(ListCase3[[7]], sep="", file=fileRdme, append=TRUE)
                }## if(i==13)

                ## cat("\n", sep="", file=fileRdme, append=TRUE)
                ## cat(rep("=",times=100), "\n", sep="", file=fileRdme,
                ##     append=TRUE)
            }## if(i%in%c(6,12) & Case==3)

            lineEND <- paste(" ** ", names(ListGlossary)[i], " : ",
                             ListGlossary[[i]], "\n", collapse="")

            cat(c("\n"), sep="", file=fileRdme, append=TRUE)
            cat(lineEND, sep="", file=fileRdme, append=TRUE)
        }## for(i in 1:length(ListGlossary))
    }## if(is.null(path.result)==FALSE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    return(ListGlossary)
}## Glossary()
