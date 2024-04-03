#' @title Volcano and MA graphs
#'
#' @description The function returns Volcano plots and MA plots from
#' the results of our function
#' [DEanalysisGlobal()].
#'
#' @details
#' * If data belong to different time points only, the function returns
#' \eqn{T-1} volcano and MA plots
#' (with \eqn{T} the number of time measurements), corresponding to the
#' \eqn{log_2} fold change between each time ti and the reference time t0,
#' for all \eqn{i>0}.
#' * If data belong to different biological conditions only,
#' the function returns \eqn{(N_{bc}*(N_{bc}-1))/2} volcano and MA plots
#' (with \eqn{N_{bc}} the number of biological conditions),
#' corresponding to the \eqn{log_2} fold change between each pair of
#' biological condition.
#' * If data belong to different biological conditions and time points,
#' the function returns
#'   * \eqn{(T-1)*N_{bc}} volcano and MA plots,
#'   corresponding to the \eqn{log_2} fold change between each time ti and
#'   the reference time t0, for all biological condition.
#'   * \eqn{((T-1)*N_{bc}*(N_{bc}-1))/2} volcano and MA plots,
#'   corresponding to the \eqn{log_2} fold change between
#'   each pair of biological conditions, for all fixed time point.
#'
#' @param SEresDE A SummarizedExperiment class object. Output from
#' [DEanalysisGlobal()]
#' (see \code{Examples}).
#' @param NbGene.plotted Non negative integer. The algorithm computes the sum
#' of all the absolute \eqn{log_2} fold change present in the element
#' \code{DE.results} of \code{Res.DE.analysis} for each gene.
#' Only the highest \code{NbGene.plotted} genes are plotted in the volcano and
#' MA plots. By default, \code{NbGene.plotted}=2.
#' @param SizeLabel Numeric. Give the size of the names of plotted genes.
#' By default, \code{SizeLabel}=3.
#' @param Display.plots \code{TRUE} or \code{FALSE}. \code{FALSE} as default.
#' If \code{TRUE}, all graphs will be plotted.
#' Otherwise no graph will be plotted.
#' @param Save.plots \code{TRUE} or \code{FALSE} or a Character.
#' \code{FALSE} as default. Path to save the Volcano and MA plots.
#' If \code{NULL}, the Volcano and MA plots will not be saved in a sub folder
#' in \code{path.result}.
#'
#' If \code{path.result} contains a sub folder entitled "VolcanoPlots",
#' all the Volcano plots will be saved in the sub folder "VolcanoPlots".
#' Otherwise, a sub folder entitled "VolcanoPlots" will be created
#' in \code{path.result} and all the Volcano plots will be saved in
#' the sub folder created.
#'
#' If \code{path.result} contains a sub folder entitled "MAplots",
#' all the MA plots will be saved in the sub folder "MAplots".
#' Otherwise, a sub folder entitled "MAplots" will be created in
#' \code{path.result} and all the MA plots will be saved in
#' the sub folder created.
#'
#' @importFrom SummarizedExperiment rowData rownames assays
#' @importFrom S4Vectors metadata
#' @importFrom ggplot2 ggplot geom_point aes scale_shape_manual
#' scale_color_manual scale_color_manual geom_hline ggtitle xlab theme_minimal
#' @importFrom ggrepel geom_label_repel
#'
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEresDE} with Volcano plots and MA plots from the results of
#' our function [DEanalysisGlobal()],
#' all saved in the metadata \code{Results[[2]][[3]]} of \code{SEresDE}.
#'
#' @seealso The function calls the output of
#' [DEanalysisGlobal()].
#'
#' @export
#'
#' @examples
#' ## data importation
#' data(RawCounts_Antoszewski2022_MOUSEsub500)
#' ## No time points. We take only two groups for the speed of the example
#' dataT1wt <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200), seq_len(7)]
#'
#' ## Preprocessing with Results of DEanalysisGlobal()
#' resDATAprepSE <- DATAprepSE(RawCounts=dataT1wt,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=NULL,
#'                             Individual.position=2)
#' ##------------------------------------------------------------------------##
#' ## DE analysis
#' resDET1wt <- DEanalysisGlobal(SEres=resDATAprepSE,
#'                               pval.min=0.05,
#'                               pval.vect.t=NULL,
#'                               log.FC.min=1,
#'                               LRT.supp.info=FALSE,
#'                               Plot.DE.graph=FALSE,
#'                               path.result=NULL,
#'                               Name.folder.DE=NULL)
#'
#' ##------------------------------------------------------------------------##
#' ## Volcano MA
#' resVolcanoMA <- DEplotVolcanoMA(SEresDE=resDET1wt,
#'                                 NbGene.plotted=5,
#'                                 Display.plots=TRUE,
#'                                 Save.plots=FALSE)

DEplotVolcanoMA <- function(SEresDE,
                            NbGene.plotted=2,
                            SizeLabel=3,
                            Display.plots=FALSE,
                            Save.plots=FALSE) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check
    resErr <- ErrVolcanoMA(SEresDE=SEresDE,
                           NbGene.plotted=NbGene.plotted,
                           SizeLabel=SizeLabel,
                           Display.plots=Display.plots,
                           Save.plots=Save.plots)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Preprocessing
    DEresults <- SummarizedExperiment::rowData(SEresDE)
    Name.G <- as.character(SummarizedExperiment::rownames(SEresDE))
    dataRLE <- data.frame(SummarizedExperiment::assays(SEresDE)$rle)
    resPATH <- S4Vectors::metadata(SEresDE)$DESeq2obj$pathNAME
    resInputs <- S4Vectors::metadata(SEresDE)$DESeq2obj$Summary.Inputs

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## To avoid "no visible binding for global variable" with devtools::check()
    log2FoldChange <- padj <- Gene <- ColpFC <- Shape.pFC <- NULL

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Folder path and creation
    if (!isFALSE(Save.plots)) {
        if (isTRUE(Save.plots)) {
            path.result <- resPATH$Path.result
        } else {
            path.result <- Save.plots
        }## if(isTRUE(Save.plots))

        if (!is.null(path.result)) {
            if (!is.null(resPATH$Folder.result)) {
                SufixDE <- paste0("_", resPATH$Folder.result)
            } else {
                SufixDE <- NULL
            }## if(!is.null(resPATH$Folder.result))

            SuppPlotFolder <- paste0("2-4_Supplementary_Plots", SufixDE)

            if (!SuppPlotFolder%in%dir(path=path.result)) {
                print("Folder creation")
                dir.create(path=file.path(path.result, SuppPlotFolder))
            }## if(!SuppPlotFolder%in%dir(path = path.result))

            path.result.f <- file.path(path.result, SuppPlotFolder)

        } else {
            path.result.f <- NULL
        }## if(!is.null(path.result))

        if (!is.null(path.result.f)) {
            if (!"Plots_Volcano"%in%dir(path=path.result.f)) {
                dir.create(path=file.path(path.result.f, "Plots_Volcano"))
            }## if(!"Plots_Volcano"%in%dir(path=path.result.f))
            path.result.vol <- file.path(path.result.f, "Plots_Volcano")

            if (!"Plots_RatioIntensity(MA)"%in%dir(path=path.result.f)) {
                dir.create(path=file.path(path.result.f,
                                          "Plots_RatioIntensity(MA)"))
            }## if(!"Plots_RatioIntensity(MA)"%in%dir(path=path.result.f))
            path.result.ma <- file.path(path.result.f,
                                        "Plots_RatioIntensity(MA)")
        } else {
            path.result.vol <- NULL
            path.result.ma <- NULL
            path.result <- NULL
        }## if(!is.null(path.result.f))
    } else {
        path.result <- NULL
    }## if(!isFALSE(Save.plots))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Distinction case when samples belong to different group and/or time
    if (length(resInputs$ExprCond) == 2) {
        ##-------------------------------------------------------------------##
        NbGroup <- length(resInputs$GroupLevels)
        NbTime <- length(resInputs$TimeLevels)
        pevo <- c(rep(resInputs$pvalsTime, times=NbGroup),
                  rep(resInputs$pvalGroup, times=NbTime*NbGroup*(NbGroup-1)/2))

        Title.g <- rep(NA,times=NbGroup*(NbGroup-1)/2)
        cptBC <- 0

        ##-------------------------------------------------------------------##
        for (i in seq_len(NbGroup-1)) {
            for (k in seq(from=(i+1), to=NbGroup, by=1)) {
                cptBC <- cptBC + 1
                Title.g[cptBC] <- paste(resInputs$GroupLevels[k], "vs",
                                        resInputs$GroupLevels[i])
            }## for(k in (i+1):NbGroup)
        }## for(i in seq_len(NbGroup-1))

        ##-------------------------------------------------------------------##
        fixedTime <- rep(paste0("t", seq_len(NbTime)-1),
                         each=NbGroup*(NbGroup-1)/2)
        fixedBC <- rep(resInputs$GroupLevels, each=NbTime-1)

        Title.g2 <- paste0("(", fixedTime, ") ", Title.g)

        Title.vMA <- c(paste0("(", fixedBC, ") t",
                              rep(seq_len(NbTime-1), times=NbGroup), " vs t0"),
                       Title.g2)
        TitleMAv_R <- gsub("(", "",
                           gsub(")", "", gsub(" ", "_", Title.vMA, fixed=TRUE),
                                fixed=TRUE),
                           fixed=TRUE)

    } else {
        ##-------------------------------------------------------------------##
        if ("Time"%in%resInputs$ExprCond) {
            NbTime <- length(resInputs$TimeLevels)
            pevo <- resInputs$pvalsTime

            Title.vMA <- paste0("t", rep(seq_len(NbTime-1), times=1), " vs t0")
        } else {
            NbGroup <- length(resInputs$GroupLevels)
            pevo <- rep(resInputs$pvalGroup, times=NbGroup*(NbGroup-1))

            Title.vMA <- rep(NA, times=NbGroup*(NbGroup-1)/2)
            cptBC <- 0

            for (i in seq_len(NbGroup-1)) {
                for (k in seq(from=(i+1), to=NbGroup, by=1)) {
                    cptBC <- cptBC + 1
                    Title.vMA[cptBC] <- paste(resInputs$GroupLevels[k], "vs",
                                              resInputs$GroupLevels[i])
                }## for (k in seq(from=(i+1), to=NbGroup, by=1))
            }## for (i in seq_len(NbGroup-1))
        }## if("Time"%in%resInputs$ExprCond)

        TitleMAv_R <- Title.vMA
    }## if(length(resInputs$ExprCond)==2)

    AbsLog2Fcmin <- resInputs$logFCmin

    ##-----------------------------------------------------------------------##
    log2Mean <- log2(apply(dataRLE, 1 , mean) + 1)
    IdLog2FC <- grep(pattern="Log2FoldChange",x=colnames(DEresults),fixed=TRUE)

    ##-----------------------------------------------------------------------##
    Nlog2FC <- length(IdLog2FC)

    listVolcanos <- vector(mode="list", length=Nlog2FC)
    listMAplots <- vector(mode="list", length=Nlog2FC)
    bestGenes <- vector(mode="list", length=Nlog2FC)

    TitleMAv_Rplot <- gsub(" ", "_", TitleMAv_R, fixed=TRUE)
    names(listMAplots) <- names(listVolcanos) <- TitleMAv_Rplot
    names(bestGenes) <- TitleMAv_Rplot

    ##------------------------------------------------------------------------#
    ## Volcano and MA plots for each case
    for (i in seq_len(Nlog2FC)) {
        DatVolMA <- data.frame(Gene=Name.G, log2Mean=log2Mean,
                               log2FoldChange=NA, padj=NA,
                               DEpFC=rep("Not significant",
                                         times=length(log2Mean)),
                               ColpFC="#999999", Shape.pFC=18)
        row.names(DatVolMA) <- NULL

        ##-------------------------------------------------------------------##
        DatVolMA$log2FoldChange <- DEresults[, IdLog2FC[i]]
        DatVolMA$padj <- DEresults[, IdLog2FC[i]+1]

        ##-------------------------------------------------------------------##
        if (min(DatVolMA$padj) == 0 & length(unique(DatVolMA$padj))>1) {
            sortVolMApadj <- unique(sort(DatVolMA$padj))[2]*0.01
            DatVolMA$padj[which(DatVolMA$padj == 0)] <- sortVolMApadj
        }## if(min(DatVolMA$padj)==0 & length(unique(DatVolMA$padj))>1)

        DatVolMA$DEpFC[which(DatVolMA$padj<pevo[i])] <- "Pvalue adjusted"
        SelpvFC <- which(DatVolMA$padj < pevo[i] &
                             abs(DatVolMA$log2FoldChange) > AbsLog2Fcmin)
        DatVolMA$DEpFC[SelpvFC] <- "Pvalue adjusted & Log2FC"

        ##-------------------------------------------------------------------##
        Gene.ns <- which(DatVolMA$DEpFC == "Not significant")
        Gene.p <- which(DatVolMA$DEpFC == "Pvalue adjusted")
        Gene.pFC <- which(DatVolMA$DEpFC == "Pvalue adjusted & Log2FC")

        if (length(Gene.p)>0) {
            DatVolMA$ColpFC[Gene.p] <- "#56B4E9"
            DatVolMA$Shape.pFC[Gene.p] <- 17
        }## if(length(Gene.p)>0)

        if (length(Gene.pFC)>0) {
            DatVolMA$ColpFC[Gene.pFC] <- "#E69F00"
            DatVolMA$Shape.pFC[Gene.pFC] <- 16
        }## if(length(Gene.pFC)>0)

        ##-------------------------------------------------------------------##
        DatVolMA$Shape.pFC <- as.factor(DatVolMA$Shape.pFC)
        DatVolMA$ColpFC <- as.factor(DatVolMA$ColpFC)

        ##-------------------------------------------------------------------##
        Eucl.pFC.Volc<-sqrt(DatVolMA$log2FoldChange^2 + log10(DatVolMA$padj)^2)

        IdPlotted.Volc <- c(intersect(order(Eucl.pFC.Volc,decreasing=TRUE),
                                      Gene.pFC),
                            intersect(order(Eucl.pFC.Volc, decreasing=TRUE),
                                      Gene.p),
                            intersect(order(Eucl.pFC.Volc, decreasing=TRUE),
                                      Gene.ns))[seq_len(NbGene.plotted)]

        bestGenes[[i]] <- DatVolMA$Gene[IdPlotted.Volc]

        ##-------------------------------------------------------------------##
        options(warn = - 1)

        gVolca <- ggplot2::ggplot(data=DatVolMA,
                                  ggplot2::aes(x=log2FoldChange,
                                               y=-log10(padj),
                                               label=Gene)) +
            ggrepel::geom_label_repel(data=DatVolMA[IdPlotted.Volc,],
                                      col=DatVolMA$ColpFC[IdPlotted.Volc],
                                      box.padding=0.5,
                                      max.overlaps=NbGene.plotted+2,
                                      min.segment.length=0, size=SizeLabel) +
            ggplot2::geom_point(ggplot2::aes(color=ColpFC,shape=Shape.pFC)) +
            ggplot2::scale_shape_manual(values=c(18, 17, 16),
                                        name="Criteria",
                                        breaks=c(18, 17, 16),
                                        labels=c("Not significant",
                                                 "Pvalue adjusted",
                                                 "Pvalue adjusted & Log2FC")) +
            ggplot2::scale_color_manual(values=c("#999999", "#56B4E9",
                                                 "#E69F00"),
                                        name="Criteria",
                                        breaks=c("#999999", "#56B4E9",
                                                 "#E69F00"),
                                        labels=c("Not significant",
                                                 "Pvalue adjusted",
                                                 "Pvalue adjusted & Log2FC")) +
            ggplot2::geom_vline(xintercept=c(-AbsLog2Fcmin, AbsLog2Fcmin),
                                linetype="dashed",
                                col="black") +
            ggplot2::geom_hline(yintercept=-log10(pevo[i]),
                                linetype="dashed",
                                col="black") +
            ggplot2::ggtitle(paste0("Volcano plot : ", Title.vMA[i])) +
            ggplot2::xlab("log2 (Fold Change)") +
            ggplot2::ylab("-log10 (Pvalue adjusted)") +
            ggplot2::theme_minimal()

        listVolcanos[[i]] <- gVolca

        options(warn = 0)

        ##-------------------------------------------------------------------##
        options(warn = - 1)

        gMA <- ggplot2::ggplot(data=DatVolMA,
                               ggplot2::aes(x=log2Mean,
                                            y=log2FoldChange,
                                            label=Gene)) +
            ggrepel::geom_label_repel(data=DatVolMA[IdPlotted.Volc,],
                                      col=DatVolMA$ColpFC[IdPlotted.Volc],
                                      box.padding=0.5,
                                      max.overlaps=NbGene.plotted+2,
                                      min.segment.length=0,
                                      size=SizeLabel)+
            ggplot2::geom_point(ggplot2::aes(color=ColpFC,
                                             shape=Shape.pFC)) +
            ggplot2::scale_shape_manual(values=c(18, 17, 16),
                                        name="Criteria", ## "Attribute"
                                        breaks=c(18, 17, 16),
                                        labels=c("Not significant",
                                                 "Pvalue adjusted",
                                                 "Pvalue adjusted & Log2FC")) +
            ggplot2::scale_color_manual(values=c("#999999", "#56B4E9",
                                                 "#E69F00"),
                                        name="Criteria",
                                        breaks=c("#999999", "#56B4E9",
                                                 "#E69F00"),
                                        labels=c("Not significant",
                                                 "Pvalue adjusted",
                                                 "Pvalue adjusted & Log2FC"))+
            ggplot2::geom_hline(yintercept=-AbsLog2Fcmin,
                                linetype="dashed", col="black") +
            ggplot2::geom_hline(yintercept=0, col="black") +
            ggplot2::geom_hline(yintercept=AbsLog2Fcmin,
                                linetype="dashed",
                                col="black") +
            ggplot2::ggtitle(paste0("MA plot : ", Title.vMA[i])) +
            ggplot2::xlab("log2 (Mean rle +1)") +
            ggplot2::ylab("log2 (Fold Change)") +
            ggplot2::theme_minimal()

        listMAplots[[i]] <- gMA

        options(warn=0)

        ##-------------------------------------------------------------------##
        options(warn=-1)

        if (!is.null(path.result)) {
            if (isTRUE(Save.plots) & !is.null(resPATH$Folder.result)) {
                SuffixSaveplots <- paste0("_", resPATH$Folder.result)
            } else {
                SuffixSaveplots <- ""
            }## if(isTRUE(Save.plots)& !is.null(resPATH$Folder.result))

            TitleVol <- paste0("Volcano_",
                               gsub(" ", "_", x=Title.vMA[i], fixed=TRUE),
                               SuffixSaveplots, ".pdf")
            TitleMA <- paste0("RatioIntensityPlot(MAplot)_",
                              gsub(" ", "_", x=Title.vMA[i], fixed=TRUE),
                              SuffixSaveplots, ".pdf")

            grDevices::pdf(file=file.path(path.result.vol, TitleVol),
                           width=11, height=8)
            print(gVolca)
            grDevices::dev.off()

            grDevices::pdf(file=file.path(path.result.ma, TitleMA),
                           width=11, height=8)
            print(gMA)
            grDevices::dev.off()

        } else {
            if (isTRUE(Display.plots)) {
                print(gVolca)
                print(gMA)
            }## if(isTRUE(Display.plots))
        }## if(is.null(path.result)==FALSE)

        options(warn = 0)
        ##-------------------------------------------------------------------##
    }## for(i in seq_len(length(IdLog2FC)))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (length(resInputs$ExprCond) == 2) {
        lVOL <- split(listVolcanos, c(fixedBC, fixedTime))
        lMA <- split(listMAplots, c(fixedBC, fixedTime))
        lGenes <- split(bestGenes, c(fixedBC, fixedTime))
    } else {
        lVOL <- listVolcanos
        lMA <- listMAplots
        lGenes <- bestGenes
    }## if (length(resInputs$ExprCond) == 2)

    listMAvolcano <- list(lGenes, Volcano=lVOL, MAplots=lMA)
    names(listMAvolcano)[1] <- paste0("highest", NbGene.plotted, "DEgenes")

    SEresDE_MA <- SEresDE
    S4Vectors::metadata(SEresDE_MA)$Results[[2]][[3]] <- listMAvolcano

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    return(SEresDE=SEresDE_MA)
}## DEplotVolcanoMA()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrSEresDE <- function(SEresDE) {
    ##-----------------------------------------------------------------------##
    Err_SE <- paste0("'SEresDE' mut be the results of the function ",
                     "'DEanalysisGlobal()'.")

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
    return(Message="No error")
}## ErrSEresDE()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrSaveDisplayPlot <- function(Display.plots=FALSE,
                               Save.plots=FALSE) {
    ##-----------------------------------------------------------------------##
    ## Display.plots
    if (!isTRUE(Display.plots) & !isFALSE(Display.plots)) {
        stop("'Display.plots' must be TRUE or FALSE.")
    }## if (!isTRUE(Display.plots) & !isFALSE(Display.plots))

    ## Save.plots
    if (!isTRUE(Save.plots) & !isFALSE(Save.plots)) {
        stop("'Save.plots' must be TRUE or FALSE.")
    }## if (!isTRUE(Save.plots) & !isFALSE(Save.plots))
    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrSaveDisplayPlot()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrVolcanoMA <- function(SEresDE,
                         NbGene.plotted=2,
                         SizeLabel=3,
                         Display.plots=FALSE,
                         Save.plots=FALSE) {
    ##-----------------------------------------------------------------------##
    ## SEresDE
    Err1 <- ErrSEresDE(SEresDE=SEresDE)
    Err2 <- ErrNNI(NNI=NbGene.plotted, NNIname="NbGene.plotted")
    Err3 <- ErrSaveDisplayPlot(Display.plots=Display.plots,
                               Save.plots=Save.plots)

    ##-----------------------------------------------------------------------##
    if (!is.numeric(SizeLabel)) {
        stop("'SizeLabel' must be strictly positive numeric value")
    }## if (!is.numeric(SizeLabel))

    if (SizeLabel <= 0) {
        stop("'SizeLabel' must be strictly positive numeric value")
    }## if (SizeLabel <= 0)

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrVolcanoMA()

