#' @title GSEA analysis with gprofiler2
#'
#' @description The function realizes, from the outputs of
#' [DEanalysisGlobal()],
#' an enrichment analysis (GSEA) of a subset of genes with
#' the R package \code{gprofiler2}.
#'
#' @details If \code{ColumnsLog2ordered} is a vector of integers,
#' the rows of \code{Res.DE.analysis$DE.results} (corresponding to genes)
#' will be decreasingly ordered according to the sum of absolute \eqn{log_2}
#' fold change (the selected columns must contain \eqn{log_2} fold change
#' values) before the enrichment analysis.
#' The enrichment analysis will take into account the genes order as
#' the first genes will be considered to have the highest biological importance
#' and the last genes the lowest.
#' See the input \code{ordered_query} of
#' [gprofiler2::gost()]
#' and the vignette of \code{gprofiler2} for more details.
#'
#' We have the following three cases:
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
#' @param SEresDE A SummarizedExperiment class object. Output from
#' [DEanalysisGlobal()]
#' (see \code{Examples}).
#' @param ColumnsCriteria A vector of integers where each integer indicates
#' a column of  \code{SummarizedExperiment::rowData(SEresDE)}.
#' These columns should either contain only binary values, or may contain other
#' numerical value, in which case extracted outputs from \code{SEresDE}
#' will be those with >0 values (see \code{Details}).
#' @param Internet.Connection \code{TRUE} or \code{FALSE}.
#' \code{FALSE} by default. If the user is sure to have an internet connection,
#' the user must set \code{Internet.Connection=TRUE}, otherwise, the algorithm
#' will not run.
#' @param Set.Operation A character. The user must choose between "union"
#' (default), "intersect", "setdiff" (see \code{Details}).
#' @param ColumnsLog2ordered \code{NULL} or a vector of integers.
#' If \code{ColumnsLog2ordered} is a vector of integers, it corresponds to
#' the columns number of \code{Res.DE.analysis$DE.results}, the output of
#' [DEanalysisGlobal()],
#' which must contains \eqn{log_2} fold change values (see \code{Details}).
#' @param Organism A character indicating the organism
#' where data were taken from.
#' See vignette of the R package \code{gprofiler2} for supported organisms.
#' See [gprofiler2::gost()].
#' @param Background \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the statistical enrichment analysis to find
#' over-representation of functions from Gene Ontology (GO) and
#' biological pathways (e.g. KEGG) will be done by comparing the functions and
#' biological pathways among the selected DE genes with those associated with
#' all genes in \code{Res.DE.analysis$DE.results}.
#' If \code{FALSE}, the statistical enrichment analysis will be done
#' by comparing the functions and biological pathways among the selected
#' DE genes with all functions and biological pathways included in the database
#' of \code{gprofiler2} (link in \code{See Also}).
#' See also [gprofiler2::gost()].
#' @param MaxNumberGO An integer.
#' The user can select the \code{MaxNumberGO} most important Gene Ontology
#' (GO) names to be plotted in a lollipop graph.
#' By default, \code{MaxNumberGO=20}.
#' @param Display.plots \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, all graphs will be plotted.
#' Otherwise no graph will be plotted.
#' @param Save.plots \code{TRUE} or \code{FALSE} or a Character.
#' If \code{Save.plots=TRUE} and the output \code{path.result} of
#' [DEanalysisGlobal()] is not \code{NULL}, all files will be saved in
#' "2_SupervisedAnalysis_\code{Name.folder.DE}/
#' 2-5_Enrichment_analysis_\code{Name.folder.DE}/
#' 2-5-1_gprofiler2_results_\code{Name.folder.DE}", with \code{Name.folder.DE}
#' an input of
#' [DEanalysisGlobal()].
#' If \code{Save.plots} is a character, it must be a path and all files
#' will be saved in the sub-folder "gprofiler2_results_\code{Name.folder.DE}".
#' Otherwise, the different files will not be saved.
#'
#'
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEresDE} with
#' * a data.frame which contains the outputs of
#' [gprofiler2::gost()]
#' * a Manhattan plot showing all GO names according to their pvalue
#' * a lollipop graph showing the \code{MaxNumberGO} most important GO.
#'
#' saved in the metadata \code{Results[[2]][[5]]} of \code{SEresDE}.
#'
#' The Manhattan plot and the lollipop graph are plotted if
#' \code{Display.plots=TRUE}.
#'
#' @importFrom SummarizedExperiment rowData assays rownames
#' @importFrom S4Vectors metadata
#' @importFrom gprofiler2 gost gostplot publish_gosttable
#' @importFrom grDevices png dev.off
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_hline geom_bar
#' scale_fill_manual scale_color_manual scale_shape_manual geom_segment
#' scale_x_continuous scale_y_continuous xlab ylab ggtitle theme_classic theme
#' element_blank element_text unit rel coord_flip labs theme_bw geom_area
#' @importFrom ggrepel geom_label_repel
#'
#' @seealso The function uses the R package \code{gprofiler2}
#' \url{https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html}.
#'
#' The R package \code{gprofiler2} provides an R interface to the web toolset
#' g:Profiler \url{https://biit.cs.ut.ee/gprofiler/gost}.
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
#'
#' ##------------------------------------------------------------------------##
#' ## Internet is needed in order to run the following lines of code because
#' ## gprofileR2 needs an internet connection
#' ## DE analysis
#' # resDET1wt <- DEanalysisGlobal(SEres=resDATAprepSE,
#' #                               pval.min=0.05,
#' #                               pval.vect.t=NULL,
#' #                               log.FC.min=1,
#' #                               LRT.supp.info=FALSE,
#' #                               Plot.DE.graph=FALSE,
#' #                               path.result=NULL,
#' #                               Name.folder.DE=NULL)
#' #########
#' # resGs <- GSEAQuickAnalysis(Internet.Connection=TRUE,
#' #                            SEresDE=resDET1wt,
#' #                            ColumnsCriteria=3,
#' #                            ColumnsLog2ordered=NULL,
#' #                            Set.Operation="union",
#' #                            Organism="mmusculus",
#' #                            MaxNumberGO=20,
#' #                            Background=FALSE,
#' #                            Display.plots=TRUE,
#' #                            Save.plots=FALSE)

GSEAQuickAnalysis <- function(Internet.Connection=FALSE,
                              SEresDE,
                              ColumnsCriteria=1,
                              ColumnsLog2ordered=NULL,
                              Set.Operation="union",
                              Organism="hsapiens",
                              MaxNumberGO=20,
                              Background=FALSE,
                              Display.plots=TRUE,
                              Save.plots=FALSE){
    ##-----------------------------------------------------------------------##
    ## To avoid "no visible binding for global variable" with devtools::check()
    term_name <- term_id <- p_value <- significant <- xmh <- NULL

    LollipopPlot <- GSEAtablepng <- GSEAtable <- ManhPlot <- NULL
    path.result.GSEA <- path.result <- NULL

    GOmat <- gManhattan <- glolipop <- NULL

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check
    resErr <- ErrGprofiler2(Internet.Connection=Internet.Connection,
                            SEresDE=SEresDE,
                            ColumnsCriteria=ColumnsCriteria,
                            ColumnsLog2ordered=ColumnsLog2ordered,
                            Set.Operation=Set.Operation,
                            MaxNumberGO=MaxNumberGO,
                            Background=Background,
                            Display.plots=Display.plots,
                            Save.plots=Save.plots)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Preprocessing
    scaled.data <- round(SummarizedExperiment::assays(SEresDE)$rle, digits=3)
    GeneNames <- as.character(SummarizedExperiment::rownames(SEresDE))
    resPATH <- S4Vectors::metadata(SEresDE)$DESeq2obj$pathNAME
    resInputs <- S4Vectors::metadata(SEresDE)$DESeq2obj$Summary.Inputs

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Folder path and creation
    if (!isFALSE(Save.plots)) {
        if (isTRUE(Save.plots)) {
            path.result <- resPATH$Path.result
        } else {
            path.result <- Save.plots
        }## if (isTRUE(Save.plots))

        if (!is.null(path.result)) {
            if (!is.null(resPATH$Folder.result)) {
                SufixDE <- paste0("_", resPATH$Folder.result)
            } else {
                SufixDE <- NULL
            }## if(is.null(resPATH$Folder.result)==FALSE)

            if (isTRUE(Save.plots)) {
                SuppPlotFolder <- paste0("2-5_Enrichment_analysis", SufixDE)

                if (!SuppPlotFolder%in%dir(path=path.result)) {
                    print("Folder creation")
                    dir.create(path=file.path(path.result, SuppPlotFolder))
                }## if(SuppPlotFolder%in%dir(path = path.result)==FALSE)

                path.result.f <- file.path(path.result, SuppPlotFolder)

            } else {
                path.result.f <- path.result
            }## if (isTRUE(Save.plots))

        } else {
            path.result.f <- NULL
        }## if(is.null(path.result)==FALSE)

        if (!is.null(path.result)) {
            if (isTRUE(Save.plots)) {
                nom.dossier.result1 <- paste0("2-5-1_gprofiler2_results",
                                              SufixDE)
            } else {
                nom.dossier.result1 <- paste0("gprofiler2_results",
                                              SufixDE)
            }## if(Save.plots==TRUE)

            if (!nom.dossier.result1%in%dir(path=path.result.f)) {
                dir.create(path=file.path(path.result.f,
                                          nom.dossier.result1))
            }## if(nom.dossier.result1%in%dir(path = path.result.f)==FALSE)

            path.result.GSEA <- file.path(path.result.f,
                                          nom.dossier.result1)

            LollipopPlot <- file.path(path.result.GSEA, "LollipopChart.pdf")
            ManhPlot <- file.path(path.result.GSEA, "ManhattanPlot.pdf")
            GSEAtable <- file.path(path.result.GSEA, "GSEAtable.csv")
            GSEAtablepng <- file.path(path.result.GSEA, "GSEAtable.pdf")
        }## else{path.result.GSEA<-NULL} ## if(is.null(path.result)==FALSE)

    }## else{} ## if(isFALSE(Save.plots)==FALSE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Selection of genes according to Columns criteria and Set.operation
    ResSubDE <- DEanalysisSubData(SEresDE=SEresDE,
                                  ColumnsCriteria=ColumnsCriteria,
                                  Set.Operation=Set.Operation)

    ##-----------------------------------------------------------------------##
    ## Gene order considered as important or not
    if (!is.null(ColumnsLog2ordered)) {
        Log2FCchosen <- data.frame(SummarizedExperiment::rowData(ResSubDE))
        Log2FCchosen <- Log2FCchosen[, ColumnsLog2ordered]

        orderLog2FC <- order(abs(Log2FCchosen))
        GeneSelected.i <- rownames(ResSubDE)
        GeneSelected <- GeneSelected.i[orderLog2FC]
        QueryOrder <- TRUE
    } else {
        GeneSelected <- rownames(ResSubDE)
        QueryOrder <- FALSE
    }## if(is.null(ColumnsLog2ordered)==FALSE)

    ##-----------------------------------------------------------------------##
    ## Inputs of gprofiler2::gost()
    if (isFALSE(Background)) {
        Backgene <- NULL
        DomainScoped <- "annotated"
    } else {
        Backgene <- GeneNames
        DomainScoped <- "custom"
    }## if(Background==FALSE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (isFALSE(Internet.Connection)) {
        MessageInternet <- paste("Once the user is sure to have",
                                 "an internet connection, the user must set",
                                 "'Internet.Connection=TRUE' in order to",
                                 "realize the enrichment analysis", sep=" ")

        MessageUser <- paste("Enrichment analysis were not realized because",
                             "'Internet.Connection=FALSE'.",
                             "The user must set 'Internet.Connection=TRUE'",
                             "in order to realize the enrichment analysis",
                             sep=" ")

        listGSEA <- MessageUser

        message(MessageInternet)
    } else {
        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## gprofiler2 analysis
        gostres <- gprofiler2::gost(query=GeneSelected, organism=Organism,
                                    ordered_query=QueryOrder,
                                    multi_query=FALSE, significant=FALSE,
                                    user_threshold=0.05, exclude_iea=FALSE,
                                    measure_underrepresentation=FALSE,
                                    evcodes=TRUE, correction_method="g_SCS",
                                    domain_scope=DomainScoped,
                                    custom_bg=Backgene,
                                    sources=c("GO:BP","GO:MF","GO:CC","KEGG"),
                                    numeric_ns="", as_short_link=FALSE)
        ## term_size - number of genes that are annotated to the term (GO)
        ## query_sizes
        #### - number of genes that were included in the query in the order of
        ####   input queries
        ## intersection_sizes
        #### - the number of genes in the input query that are annotated to
        ####   the corresponding term in the order of input queries
        ## precision
        #### - the proportion of genes in the input list that are annotated to
        ####   the function (defined as intersection_size/query_size)
        ## recall
        #### - the proportion of functionally annotated genes that the query
        ####   recovers (defined as intersection_size/term_size)
        ## p <- gprofiler2::gostplot(RESgost, capped=TRUE, interactive=FALSE)

        ##-------------------------------------------------------------------##
        ## Creation of GOmat, a data.frame containing result
        orderdPva <- order(gostres$result$p_value)
        NewgoRes <- gostres$result[orderdPva,]
        row.names(NewgoRes) <- NULL
        gostres$result <- NewgoRes
        GOmat <- NewgoRes[, c(9:11, 3:2, 4:8)]

        GeneID <- paste0("(", NewgoRes$intersection, ")")
        lapply_GOparents <- lapply(NewgoRes$parents, function(x) myCollapse(x))
        GOmat$GOparents <- unlist(lapply_GOparents)
        GOmat$Gene_id <- gsub(",", ")_(", GeneID, fixed=TRUE)

        gproList <- RESgprofiler2(RESgost=GOmat, MaxNumberGO=MaxNumberGO)

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        if (!is.null(path.result)) {
            grDevices::pdf(file=ManhPlot, width=11, height=8)
            print(gManhattan)
            grDevices::dev.off()

            grDevices::pdf(file=LollipopPlot, width=11, height=8)
            print(glolipop)
            grDevices::dev.off()

            utils::write.table(GOmat, file=GSEAtable, sep=";", row.names=FALSE)

        }## if (!is.null(path.result))

        if (isTRUE(Display.plots)) {
            print(gManhattan)
            print(glolipop)
        }## if (isTRUE(Display.plots))

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## SE object ## gproList
        listGSEA <- list(GSEAresults=gproList$GOmat,
                         selectedGenes=GeneSelected,
                         lollipopChart=gproList$glolipop,
                         manhattanPlot=gproList$gManhattan)
    }## if(Internet.Connection==FALSE)

    SEresGSEA <- SEresDE
    S4Vectors::metadata(SEresGSEA)$Results[[2]][[5]] <- listGSEA

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    return(SEresGSEA)
}## GSEAQuickAnalysis()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

RESgprofiler2 <- function(RESgost, MaxNumberGO) {
    ##-----------------------------------------------------------------------##
    ## To avoid "no visible binding for global variable" with devtools::check()
    p_value <- significant <- term_name <- xmh <- NULL
    ##-----------------------------------------------------------------------##

    ## preprocessing
    GOmat <- RESgost

    GOsource <- c("GO:BP", "GO:CC", "GO:MF", "KEGG")
    GOcolors <- c("#CB2027", "#059748", "#EE7733", "#7B3A96")

    rel6 <- ggplot2::rel(0.6)
    rel7 <- ggplot2::rel(0.7)
    rel8 <- ggplot2::rel(0.8)

    ##-----------------------------------------------------------------------##
    ## Preprocessing data.frame for graphs
    GOmat$term_id <- factor(GOmat$term_id, levels=rev(GOmat$term_id))
    GOmat$significant <- factor(GOmat$significant, levels=c("TRUE", "FALSE"))
    GOmat$source <- factor(GOmat$source, levels=GOsource)

    if (max(GOmat$p_value) == 1) {
        GOmat <- GOmat[-which(GOmat$p_value == 1),]
    }## if(max(GOmat$p_value)==1)

    ##-----------------------------------------------------------------------##
    ## data.frame for lolipop graph
    DatLolipop <- GOmat[seq_len(min(c(MaxNumberGO, nrow(GOmat)))),
                        c(1, 2, 4, 5, 3)]
    ## we can use either base::substr or stringr::str_sub()
    DatLolipop$term_name <- substr(DatLolipop$term_name, 1, 50)

    index_dupliname <- which(duplicated(DatLolipop$term_name))
    N_dupliname <- length(index_dupliname)
    if (N_dupliname > 0) {
        strlast <- sample(c(47, 48, 49), size=N_dupliname, replace=TRUE)
        str_sample <- substr(DatLolipop$term_name[index_dupliname], 1, strlast)
        DatLolipop$term_name[index_dupliname] <- str_sample
    }## if (length(index_dupliname) > 0)

    DatLolipop$term_name <- factor(DatLolipop$term_name,
                                   levels=rev(DatLolipop$term_name))
    ## length(which(duplicated(DatLolipop$term_name)))

    ## Graph option
    IntegerMaxLogPval <- ceiling(max(-log10(DatLolipop$p_value)))
    SizeGOlabel <- min(1, (0.5*4*15)/nrow(DatLolipop))
    ggaty_size <- ggplot2::rel(SizeGOlabel)

    ##-----------------------------------------------------------------------##
    lvlsSignificantL <- levels(DatLolipop$significant)
    uniqueSignificantL <- unique(DatLolipop$significant)
    LevelSignificant <- which(lvlsSignificantL%in%uniqueSignificantL)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Lolipop graph ## term_name <-> term_id,
    glolipop <- ggplot2::ggplot(DatLolipop,
                                ggplot2::aes(x=term_name, y=-log10(p_value),
                                             ymin=0,
                                             ymax=IntegerMaxLogPval)) +
        ggplot2::geom_point(ggplot2::aes(colour=significant),
                            size=1.9, shape=19) +
        ggplot2::coord_flip() +
        ggplot2::geom_segment(ggplot2::aes(x=term_name, xend=term_name,
                                           y=-log10(p_value), yend=0,
                                           colour=significant),
                              linewidth=1) +
        ggplot2::geom_bar(stat="identity", color="black", size=0.5,
                          ggplot2::aes(fill=source,
                                       y=-0.10*(IntegerMaxLogPval*0.55))) +
        ggplot2::geom_hline(yintercept=-log10(0.05), linetype="dashed",
                            size=0.6) +
        ggplot2::scale_color_manual(values=c("#E69F00",
                                             "#56B4E9")[LevelSignificant]) +
        ggplot2::scale_fill_manual(values=GOcolors, name="Source") +
        ggplot2::xlab("") + ggplot2::ylab("-log10 (pvalue)") +
        ggplot2::ggtitle("gprofiler2 enrichments analysis") +
        ggplot2::labs(color = "Significant\n  (<0.05)") +
        ggplot2::theme_bw() + ## theme_linedraw()+
        ggplot2::theme(legend.position="right", legend.box="vertical",
                       axis.text.y=ggplot2::element_text(face="bold",
                                                         size=ggaty_size),
                       legend.title=ggplot2::element_text(face="bold",
                                                          size=rel8),
                       legend.text=ggplot2::element_text(size=rel6),
                       legend.key.size=ggplot2::unit(0.4, "cm"))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## data.frame for Manhattan graph
    if (length(which(GOmat$p_value == 1)) > 0) {
        DatManhattan <- GOmat[-which(GOmat$p_value == 1), c(1, 2, 4, 5, 3)]
    } else {
        DatManhattan <- GOmat[, c(1, 2, 4, 5, 3)]
    }## if(length(which(GOmat$p_value == 1)) > 0)

    DatManhattan <- DatManhattan[order(DatManhattan$source,
                                       DatManhattan$p_value),]
    DatManhattan$xmh <- seq_len(nrow(DatManhattan))

    DatManhattan$Color <- factor(DatManhattan$source, levels=GOsource)
    levels(DatManhattan$Color) <- GOcolors

    DatManhattan$Point <- DatManhattan$significant
    levels(DatManhattan$Point) <- c(8, 16)

    ##-----------------------------------------------------------------------##
    ## Option for the Manhattan graph
    NbPerGOw <- as.numeric(table(DatManhattan$source))

    if (length(which(NbPerGOw == 0)) > 0) {
        NbPerGOw <- NbPerGOw[-which(NbPerGOw == 0)]
    }## if(length(which(NbPerGOw == 0))>0)

    xGO <- cumsum(NbPerGOw) + 1 - ((NbPerGOw+1)/2)
    SelMh <- seq_len(max(ceiling(nrow(DatManhattan)*0.05), 10))
    ## 1:max(ceiling(nrow(DatManhattan)*0.05),10)
    ## min(c(MaxNumberGO,nrow(DatManhattan)))

    ##-----------------------------------------------------------------------##
    lvlsSignificantM <- levels(DatManhattan$significant)
    uniqueSignificantM <- unique(DatManhattan$significant)
    LvlsSignificantManh <- which(lvlsSignificantM%in%uniqueSignificantM)
    DatManhattan_repel <- DatManhattan[order(DatManhattan$p_value),][SelMh,]

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Manhattan plot ## term_name <-> term_id,
    gManhattan <- ggplot2::ggplot(DatManhattan,
                                  ggplot2::aes(x=xmh, y=-log10(p_value),
                                               label=term_name, ymin=0,
                                               ymax=IntegerMaxLogPval)) +
        ggplot2::geom_point(ggplot2::aes(color=source, shape=significant),
                            size=1.7) + ##alpha=0.8,
        ggplot2::geom_line(ggplot2::aes(color=as.character(source)), size=0.4)+
        ggplot2::geom_hline(yintercept=-log10(0.05), linetype="dashed",
                            size=0.6) +
        ggplot2::geom_area(ggplot2::aes(fill=as.character(source),
                                        group=as.character(source)),
                           alpha=0.2, position='identity') +
        ggplot2::geom_bar(stat="identity", size=1.5,
                          ggplot2::aes(fill=source,color=source,
                                       y=-0.06*(IntegerMaxLogPval*0.55))) +
        ggplot2::scale_fill_manual(values=GOcolors, name="Source")+
        ggplot2::scale_color_manual(values=GOcolors, name="Source") +
        ggplot2::scale_shape_manual(values=c(8, 16)[LvlsSignificantManh],
                                    name="Significant\n   (<0.05)") +
        ggrepel::geom_label_repel(data=DatManhattan_repel,
                                  col=DatManhattan_repel$Color, size=2,
                                  max.overlaps=nrow(DatManhattan),
                                  box.padding=0.2, min.segment.length=0)+
        ggplot2::scale_x_continuous(label=unique(DatManhattan$source),
                                    breaks=xGO,
                                    guide=ggplot2::guide_axis(angle=90)) +
        ggplot2::scale_y_continuous(limits=c(-0.06*(IntegerMaxLogPval*0.55),
                                             IntegerMaxLogPval),
                                    expand=c(0, 0)) +
        ggplot2::xlab("") + ggplot2::ylab("-log10 (pvalue)") +
        ggplot2::ggtitle("Manhattan plot") +
        ggplot2::theme_classic()+ ## theme_bw() +
        ggplot2::theme(legend.position="right", legend.box="vertical",
                       panel.border=ggplot2::element_blank(),
                       panel.grid.major.x=ggplot2::element_blank(),
                       panel.grid.minor.x=ggplot2::element_blank(),
                       legend.title=ggplot2::element_text(face="bold",
                                                          size=rel8),
                       legend.text=ggplot2::element_text(size=rel7),
                       legend.key.size=ggplot2::unit(0.5, "cm"))

    ##-----------------------------------------------------------------------##
    return(list(GSEAresults=GOmat,
                glolipop=glolipop,
                gManhattan=gManhattan))
}## RESgprofiler2()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

myCollapse <- function(x) {
    resCollapse <- paste0("(", paste(x, collapse=")_("), ")")
    return(resCollapse)
}## myCollapse()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrGprofiler2 <- function(Internet.Connection=FALSE,
                          SEresDE,
                          ColumnsCriteria=1, ColumnsLog2ordered=NULL,
                          Set.Operation="union",
                          MaxNumberGO=20, Background=FALSE,
                          Display.plots=TRUE, Save.plots=FALSE) {
    ##-----------------------------------------------------------------------##
    ## SEresDE
    Err1 <- ErrSEresDE(SEresDE=SEresDE)
    Err2 <- ErrColumnsCriteria(ColumnsCriteria=ColumnsCriteria,
                               Set.Operation=Set.Operation)
    Err3 <- ErrNNIvector(NNIs=ColumnsLog2ordered, NNIname="ColumnsLog2ordered")
    Err4 <- ErrNNI(NNI=MaxNumberGO, NNIname="MaxNumberGO")
    Err5 <- ErrSaveDisplayPlot(Display.plots=Display.plots,
                               Save.plots=Save.plots)

    ##-----------------------------------------------------------------------##
    if (!isTRUE(Internet.Connection) & !isFALSE(Internet.Connection)) {
        stop("'Internet.Connection' must be TRUE or FALSE.")
    }## if (!isTRUE(Internet.Connection) & !isFALSE(Internet.Connection))

    if (!isTRUE(Background) & !isFALSE(Background)) {
        stop("'Background' must be TRUE or FALSE.")
    }## if (!isTRUE(Background) & !isFALSE(Background))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrGprofiler2()

