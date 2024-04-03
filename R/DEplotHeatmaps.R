#' @title Heatmaps of DE genes
#'
#' @description The function returns two heatmaps:
#' one heatmap of gene expressions between samples and selected genes and
#' a correlation heatmap between samples from the output of
#' [DEanalysisGlobal()].
#'
#' @details We have the following three cases:
#' * If \code{Set.Operation="union"} then the rows extracted from
#' the different datasets (raw counts, normalized data and
#' \code{SummarizedExperiment::rowData(SEresDE)})
#' included in the SummarizedExperiment class object \code{SEresDE}
#' are those such that the sum of the selected columns of
#' \code{SummarizedExperiment::rowData(SEresDE)}
#' given in \code{ColumnsCriteria} is >0.
#' This means that the selected genes are those having at least one ’1’
#' in one of the selected columns.
#' * If \code{Set.Operation="intersect"} then the rows extracted from
#' the different datasets (raw counts, normalized data and
#' \code{SummarizedExperiment::rowData(SEresDE)})
#' included in the SummarizedExperiment class object \code{SEresDE}
#' are those such that the product of the selected columns of
#' \code{SummarizedExperiment::rowData(SEresDE)}
#' given in \code{ColumnsCriteria} is >0.
#' This means that the selected genes are those having ’1’
#' in all of the selected columns.
#' * If \code{Set.Operation="setdiff"} then the rows extracted from
#' the different datasets (raw counts, normalized data and
#' \code{SummarizedExperiment::rowData(SEresDE)})
#' included in the SummarizedExperiment class object \code{SEresDE}
#' are those such that only one element of the selected columns of
#' \code{SummarizedExperiment::rowData(SEresDE)}
#' given in \code{ColumnsCriteria} is >0.
#' This means that the selected genes are those having ’1’
#' in only one of the selected columns.
#'
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
#' @param NbGene.analysis An integer or \code{NULL}.
#' If it is an integer, the heatmaps will be plotted with the
#' \code{NbGene.analysis} genes which have the highest sum of absolute
#' log2 fold change, among the DE genes selected using \code{ColumnsCriteria}
#' and \code{Set.Operation}.
#' If \code{NULL}, all the DE selected genes will be used for both heatmaps.
#' @param Color.Group \code{NULL} or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute
#' a color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param SizeLabelRows Numeric >0.
#' Size of the labels for the genes in the heatmaps.
#' @param SizeLabelCols Numeric >0.
#' Size of the labels for the samples in the heatmaps.
#' @param Display.plots \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, all graphs will be plotted.
#' Otherwise no graph will be plotted.
#' @param Save.plots TRUE or FALSE or a Character.
#' If \code{Save.plots=FALSE}, the different files will not be saved.
#' If \code{Save.plots=TRUE} and the \code{path.result} of
#' [DEanalysisGlobal()]
#' is not NULL, all files will be saved in
#' "2_SupervisedAnalysis_\code{Name.folder.DE}/
#' 2-4_Supplementary_Plots_\code{Name.folder.DE}/Plots_Heatmaps".
#' If \code{Save.plots} is a character, it must be a path and
#' all files will be saved in the sub-folder "Plots_Heatmaps".
#'
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEresDE} with two heatmaps
#' saved in the metadata \code{Results[[2]][[4]]} of \code{SEresDE}
#' * a correlation heatmap between samples (correlation heatmap)
#' * a heatmap across samples and genes called Zscore heatmap,
#' for a subset of genes that can be selected by the user.
#'
#' The two heatmaps are plotted if \code{Display.plots=TRUE}.
#' The second heatmap is build from the normalized
#' count data after being both centered and scaled (Zscore).
#'
#'
#' @seealso The function calls the function
#' [ComplexHeatmap::Heatmap()]
#' in order to plot the Heatmaps.
#'
#' @importFrom stats cor
#' @importFrom FactoMineR HCPC
#' @importFrom ComplexHeatmap HeatmapAnnotation rowAnnotation Heatmap
#' @importFrom grid gpar
#' @importFrom grDevices pdf dev.off
#'
#' @export
#'
#' @examples
#' ## data importation
#' data("RawCounts_Antoszewski2022_MOUSEsub500")
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
#' resHeatmap <- DEplotHeatmaps(SEresDE=resDET1wt,
#'                              ColumnsCriteria=3, ## Specific genes N1haT1ko
#'                              Set.Operation="union",
#'                              NbGene.analysis=20,
#'                              Color.Group=NULL,
#'                              SizeLabelRows=5,
#'                              SizeLabelCols=5,
#'                              Display.plots=TRUE,
#'                              Save.plots=FALSE)

DEplotHeatmaps <- function(SEresDE,
                           ColumnsCriteria=2,
                           Set.Operation="union",
                           NbGene.analysis=20,
                           Color.Group=NULL,
                           SizeLabelRows=5,
                           SizeLabelCols=5,
                           Display.plots=TRUE,
                           Save.plots=FALSE){
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check
    resErr <- ErrHeatmaps(SEresDE=SEresDE,
                          ColumnsCriteria=ColumnsCriteria,
                          Set.Operation=Set.Operation,
                          NbGene.analysis=NbGene.analysis,
                          Color.Group=Color.Group,
                          SizeLabelRows=SizeLabelRows,
                          SizeLabelCols=SizeLabelCols,
                          Display.plots=Display.plots,
                          Save.plots=Save.plots)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Preprocessing
    DEresults <- SummarizedExperiment::rowData(SEresDE)
    Name.G <- as.character(SummarizedExperiment::rownames(SEresDE))
    RleDat <- data.frame(SummarizedExperiment::assays(SEresDE)$rle)
    resPATH <- S4Vectors::metadata(SEresDE)$DESeq2obj$pathNAME
    resInputs <- S4Vectors::metadata(SEresDE)$DESeq2obj$Summary.Inputs

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Sub data preprocessing for both heatmaps
    if (!is.null(ColumnsCriteria)) {
        ##-------------------------------------------------------------------##
        ResSubDE <- DEanalysisSubData(SEresDE=SEresDE,
                                      ColumnsCriteria=ColumnsCriteria,
                                      Set.Operation=Set.Operation)

        rowSel <- S4Vectors::metadata(ResSubDE)$DEselection$DEselectedGenes

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## Selection of genes with the highest sum of absolute log2 fold change
        IdLog2FC <- grep(pattern="Log2FC", x=colnames(DEresults), fixed=TRUE)
        SubLog2FC <- data.frame(DEresults[rowSel, IdLog2FC])
        SumLog2FC <- apply(SubLog2FC, 1, function(x) sum(abs(x)))
        OrderLog2FC.ini <- order(SumLog2FC, decreasing=TRUE)

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ##as.matrix(ResSubDE$SubData[,-1])
        RleSubCount <- SummarizedExperiment::assays(ResSubDE)$rle
        NbSelectedRows <- length(rowSel)

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        if (!is.null(NbGene.analysis)) {
            NbGene.analysis.f <- min(NbGene.analysis, NbSelectedRows)
            OrderLog2FC <- OrderLog2FC.ini[seq_len(NbGene.analysis.f)]
        } else {
            OrderLog2FC <- OrderLog2FC.ini
        }## if(!is.null(NbGene.analysis))

        RleScaled_orderFC <- RleScaled <- t(scale(t(RleSubCount)))

        if (nrow(RleScaled) < 3) {
            rowStop <- paste("Only", nrow(RleScaled), "genes",
                             "have been selected.", "The function will work",
                             "with at least three genes.")
            stop(rowStop)
        }## if (nrow(RleScaled) < 3)

        CorSample <- stats::cor(RleSubCount)

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        RleLog2corInd <- FactoMineR::HCPC(data.frame(t(RleScaled)),
                                          method="ward", graph=FALSE)
        RleLog2corVar <- FactoMineR::HCPC(data.frame(RleScaled),
                                          method="ward", graph=FALSE, min=2)
    } else {
        RleSubCount <- as.matrix(RleDat)
        SumLog2FC <- apply(data.frame(RleSubCount), 1, var)
        OrderLog2FC.ini <- order(SumLog2FC, decreasing=TRUE)
        NbSelectedRows <- length(SumLog2FC)

        ##-------------------------------------------------------------------##
        if (!is.null(NbGene.analysis)) {
            NbGene.analysis.f <- min(NbGene.analysis, NbSelectedRows)
            OrderLog2FC <- OrderLog2FC.ini[seq_len(NbGene.analysis.f)]
        } else {
            OrderLog2FC <- OrderLog2FC.ini
        }## if(!is.null(NbGene.analysis))

        RleScaled <- t(scale(t(RleSubCount)))
        CorSample <- stats::cor(RleSubCount)
        RleScaled_orderFC <- RleScaled[OrderLog2FC,]

        ##-------------------------------------------------------------------##
        RleLog2corInd<-FactoMineR::HCPC(data.frame(t(RleScaled_orderFC)),
                                        method="ward", graph=FALSE)
        RleLog2corVar<-FactoMineR::HCPC(data.frame(RleScaled_orderFC),
                                        method="ward", graph=FALSE)
    }## if(is.null(ColumnsCriteria))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (length(resInputs$ExprCond) == 2) {
        ##-------------------------------------------------------------------##
        vectTIMEini <- resInputs$FactorsInfo$Time
        Vect.Time <- paste0("t", gsub("T", "", gsub("t", "", vectTIMEini)))
        Tlevels <- levels(factor(Vect.Time))
        NbTime <- length(Tlevels)

        ##-------------------------------------------------------------------##
        MypaletteT <- myPaletteT(Nt=NbTime)
        Color.Time <- data.frame(Name=Tlevels, Col=MypaletteT)

        ## Color.Time <- NULL
        ## if (is.null(Color.Time)) {
        ##     ColorTime0 <- scales::hue_pal()(length(Tlevels)-1)
        ##     Color.Time <- data.frame(Name=Tlevels,
        ##                              Col=c("#737373", ColorTime0))
        ## } else {
        ##     Id.LevelColT <- order(Color.Time[, 1])
        ##     Color.Time <- data.frame(Name=Tlevels,
        ##                              Col=Color.Time[Id.LevelColT, 2])
        ## }## if(is.null(Color.Time))
        ##-------------------------------------------------------------------##
        fa1.t <- eval(parse(text=paste0("c(",
                                        paste0("'", Tlevels, "'='",
                                               Color.Time[, 2], "'",
                                               collapse=","),
                                        ")")))

        ##-------------------------------------------------------------------##
        Vect.Group <- resInputs$FactorsInfo$Group
        Glevels <- resInputs$GroupLevels
        NbGroup <- length(Glevels)

        ##-------------------------------------------------------------------##
        if (is.null(Color.Group)) {
            MypaletteG <- myPaletteBC(Nbc=NbGroup)
            Color.Group <- data.frame(Name=Glevels, Col=MypaletteG)
        } else {
            Id.LevelCol.G <- order(Color.Group[, 1])
            Color.Group <- data.frame(Name=Glevels,
                                      Col=Color.Group[Id.LevelCol.G, 2])
        }## if(is.null(Color.Group))

        ##-------------------------------------------------------------------##
        fa2.g <- eval(parse(text=paste0("c(", paste0("'", Glevels, "'='",
                                                     Color.Group[, 2], "'",
                                                     collapse=","),
                                        ")")))

        listGTtitle <- list(Condition=list(title="Biological condition"))

        ##-------------------------------------------------------------------##
        AnnotCplxHeat<-ComplexHeatmap::HeatmapAnnotation(
            Time=Vect.Time,
            Condition=Vect.Group,
            col=list(Time=fa1.t, Condition=fa2.g),
            annotation_legend_param=listGTtitle)
        AnnotCplxHeat2 <- ComplexHeatmap::rowAnnotation(
            Time=Vect.Time, Condition=Vect.Group,
            col=list(Time=fa1.t, Condition=fa2.g),
            annotation_legend_param=listGTtitle, show_legend=FALSE)
    } else {
        ##-------------------------------------------------------------------##
        if ("Time"%in%resInputs$ExprCond) {
            vectTIMEini <- resInputs$FactorsInfo$Time
            Vect.Time <- paste0("t", gsub("T", "", gsub("t", "", vectTIMEini)))
            Tlevels <- levels(factor(Vect.Time))
            NbTime <- length(Tlevels)

            ##---------------------------------------------------------------##
            MypaletteT <- myPaletteT(Nt=NbTime)
            Color.Time <- data.frame(Name=Tlevels, Col=MypaletteT)

            ## Color.Time <- NULL
            ## if (is.null(Color.Time)) {
            ##     ColorTime0 <- scales::hue_pal()(length(Tlevels)-1)
            ##     Color.Time <- data.frame(Name=Tlevels,
            ##                              Col=c("#737373", ColorTime0))
            ## } else {
            ##     Id.LevelColT <- order(Color.Time[,1])
            ##     Color.Time <- data.frame(Name=Tlevels,
            ##                              Col=Color.Time[Id.LevelColT,2])
            ## }## if(is.null(Color.Time))

            ##---------------------------------------------------------------##
            fa1.t <- eval(parse(text=paste0("c(", paste0("'", Tlevels, "'='",
                                                         Color.Time[,2], "'",
                                                         collapse=","),
                                            ")")))
            listTfa1 <- list(Time=fa1.t)

            AnnotCplxHeat <- ComplexHeatmap::HeatmapAnnotation(Time=Vect.Time,
                                                               col=listTfa1)
            AnnotCplxHeat2 <- ComplexHeatmap::rowAnnotation(Time=Vect.Time,
                                                            col=listTfa1,
                                                            show_legend=FALSE)
        } else {
            Vect.Group <- resInputs$FactorsInfo$Group
            Glevels <- resInputs$GroupLevels
            NbGroup <- length(Glevels)

            ##---------------------------------------------------------------##
            if (is.null(Color.Group)) {
                MypaletteG <- myPaletteBC(Nbc=NbGroup)
                Color.Group <- data.frame(Name=Glevels, Col=MypaletteG)
            } else {
                Id.LevelCol.G <- order(Color.Group[,1])
                Color.Group <- data.frame(Name=Glevels,
                                          Col=Color.Group[Id.LevelCol.G, 2])
            }## if(is.null(Color.Group))

            ##---------------------------------------------------------------##
            fa2.g <- eval(parse(text=paste0("c(", paste0("'", Glevels, "'='",
                                                         Color.Group[, 2],
                                                         "'", collapse=","),
                                            ")")))

            listGfa2 <- list(Condition=fa2.g)
            listGtitle <- list(Condition=list(title="Biological condition"))

            AnnotCplxHeat <- ComplexHeatmap::HeatmapAnnotation(
                Condition=Vect.Group,
                col=listGfa2,
                annotation_legend_param=listGtitle)
            AnnotCplxHeat2 <- ComplexHeatmap::rowAnnotation(
                Condition=Vect.Group,
                show_legend=FALSE,
                col=listGfa2,
                annotation_legend_param=listGtitle)
        }## if("Time"%in%resInputs$ExprCond)
    }## if(length(resInputs$ExprCond)==2)

    ##-----------------------------------------------------------------------##
    gp_fontsize_Cols <- grid::gpar(fontsize=SizeLabelCols)
    gp_fontsize_Rows <- grid::gpar(fontsize=SizeLabelRows)

    ## Heatmaps of the selected genes
    H.SG <- ComplexHeatmap::Heatmap(RleScaled_orderFC, name="Z-score",
                                    top_annotation=AnnotCplxHeat,
                                    row_split=RleLog2corVar$call$t$nb.clust,
                                    column_split=RleLog2corInd$call$t$nb.clust,
                                    column_title="Samples", row_title="Genes",
                                    column_names_gp=gp_fontsize_Cols,
                                    row_names_gp=gp_fontsize_Rows)

    ## Correlation heatmaps
    H.cor <- ComplexHeatmap::Heatmap(CorSample, name="Correlation",
                                     top_annotation=AnnotCplxHeat,
                                     left_annotation=AnnotCplxHeat2,
                                     row_split=RleLog2corInd$call$t$nb.clust,
                                     column_split=RleLog2corInd$call$t$nb.clust,
                                     column_title="Samples",
                                     row_title="Samples",
                                     column_names_gp=gp_fontsize_Cols,
                                     row_names_gp=gp_fontsize_Cols)

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
            }## if(!SuppPlotFolder%in%dir(path=path.result))

            path.result.f <- file.path(path.result, SuppPlotFolder)

        } else {
            path.result.f <- NULL
        }## if(!is.null(path.result))

        if (!is.null(path.result.f)) {
            if (!"Plots_Heatmaps"%in%dir(path=path.result.f)) {
                dir.create(path=file.path(path.result.f, "Plots_Heatmaps"))
            }## if(!"Plots_Heatmaps"%in%dir(path = path.result.f))
            path.result.Heat <- file.path(path.result.f, "Plots_Heatmaps")
        } else {
            path.result.f <- path.result <- NULL
        }## if (!is.null(path.result.f))
    } else {
        path.result <- NULL
    }## if (!isFALSE(Save.plots))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (!is.null(path.result)) {
        grDevices::pdf(file=file.path(path.result.Heat,
                                      paste0("Heatmap_Zscore_samples_vs_genes",
                                             "_most_expressed_genes", ".pdf")),
                       width=11, height=8)
        print(H.SG)
        grDevices::dev.off()
        #
        grDevices::pdf(file=file.path(path.result.Heat,
                                      paste0("Heatmap_Correlation_samples",
                                             ".pdf")),
                       width=11, height=8)
        print(H.cor)
        grDevices::dev.off()
    }## if(!is.null(path.result))

    if (isTRUE(Display.plots)) {
        print(H.SG)
        print(H.cor)
    }## if(isTRUE(Display.plots))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    resHM <- list(Zscores=RleScaled_orderFC, CorrelationMatrix=CorSample,
                  Heatmap_Zscore=H.SG, Heatmap_Correlation=H.cor)

    SEresDE_HM <- SEresDE
    S4Vectors::metadata(SEresDE_HM)$Results[[2]][[4]] <- resHM

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    return(SEresDE_HM)
}## DEplotHeatmaps()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrNNIvector <- function(NNIs, NNIname) {
    ##-----------------------------------------------------------------------##
    vNNImessage <- paste0("'", NNIname, "'", " must be a non-negative integer,",
                          " a vector of non-negative integers or 'NULL'")

    if (!is.null(NNIs)) {
        if (!is.numeric(NNIs)) {
            stop(vNNImessage)
        }## if (!is.numeric(NNIs))

        sumNNI <- sum(NNIs)
        minNNI <- min(NNIs)

        if (floor(sumNNI) != sumNNI | minNNI < 1) {
            stop(vNNImessage)
        }## if (floor(sumNNI) != sumNNI | minNNI < 1)
    }## if (!is.null(NNIs))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrNonNegativeInteger()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##


ErrColumnsCriteria <- function(ColumnsCriteria=2,
                               Set.Operation="union") {
    ##-----------------------------------------------------------------------##
    ## Check ColumnsCriteria
    Err1 <- ErrNNIvector(NNIs=ColumnsCriteria, NNIname="ColumnsCriteria")

    ## Different Set.Operation
    if (!Set.Operation%in%c("union", "intersect", "setdiff")) {
        stop("'Set.Operation' mut be 'union', 'intersect' or 'setdiff'")
    }## if (!Set.Operation%in%c("union", "intersect", "setdiff"))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrColumnsCriteria

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrHeatmaps <- function(SEresDE,
                        ColumnsCriteria=2,
                        Set.Operation="union",
                        NbGene.analysis=20,
                        Color.Group=NULL,
                        SizeLabelRows=5,
                        SizeLabelCols=5,
                        Display.plots=TRUE,
                        Save.plots=FALSE) {
    ##-----------------------------------------------------------------------##
    ## SEresDE
    Err1 <- ErrSEresDE(SEresDE=SEresDE)
    Err2 <- ErrColumnsCriteria(ColumnsCriteria=ColumnsCriteria,
                               Set.Operation=Set.Operation)
    Err3 <- ErrNNI(NNI=NbGene.analysis, NNIname="NbGene.analysis")
    Err4 <- ErrSaveDisplayPlot(Display.plots=Display.plots,
                               Save.plots=Save.plots)

    ##-----------------------------------------------------------------------##
    ## Color.Group
    if (!is.null(Color.Group)) {
        if (!is.data.frame(Color.Group)) {
            stop("'Color.Group' must be NULL or a data.frame.")
        }## if (!is.data.frame(Color.Group))
    }## if (!is.null(Color.Group))

    ## 'SizeLabelRows' and 'SizeLabelCols'
    Err_size <- paste("'SizeLabelRows' and 'SizeLabelCols' must be",
                      "strictly positive numeric values")

    if (!is.numeric(SizeLabelRows) | !is.numeric(SizeLabelCols)) {
        stop(Err_size)
    }## if (!is.numeric(SizeLabelRows) | !is.numeric(SizeLabelCols))

    if (min(c(SizeLabelRows, SizeLabelCols)) < 0) {
        stop(Err_size)
    }## if (min(c(SizeLabelRows, SizeLabelCols))<0)

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrVolcanoMA()
