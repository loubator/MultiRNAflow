#' @title Heatmaps of DE genes
#'
#' @description The function returns two heatmaps:
#' one heatmap of gene expressions between samples and selected genes and
#' a correlation heatmap between samples from the output of
#' [DEanalysisGlobal()].
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
#' @return The function returns two heatmaps:
#' one heatmap of gene expressions between samples and selected genes;
#' and a correlation heatmap between samples.
#'
#' @seealso The function calls the function
#' [ComplexHeatmap::Heatmap()]
#' in order to plot the Heatmaps.
#'
#' @importFrom stats cor
#' @importFrom FactoMineR HCPC
#' @importFrom ComplexHeatmap HeatmapAnnotation rowAnnotation Heatmap
#' @importFrom grid gpar
#' @importFrom scales hue_pal
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices pdf dev.off
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
#' ##-------------------------------------------------------------------------#
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
#' ##-------------------------------------------------------------------------#
#' resHeatmap <- DEplotHeatmaps(SEresDE=resDET1wt,
#'                              ColumnsCriteria=3, ## Specific genes N1haT1ko
#'                              Set.Operation="union",
#'                              NbGene.analysis=20,
#'                              Color.Group=NULL,
#'                              SizeLabelRows=5,
#'                              SizeLabelCols=5,
#'                              Display.plots=TRUE,
#'                              Save.plots=FALSE)

DEplotHeatmaps<-function(SEresDE,
                         ColumnsCriteria=2,
                         Set.Operation="union",
                         NbGene.analysis=20,
                         Color.Group=NULL,
                         SizeLabelRows=5,
                         SizeLabelCols=5,
                         Display.plots=TRUE,
                         Save.plots=FALSE){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check 1
    ## DATAprepSE
    Err_SE <- paste0("'SEresDE' mut be the results of the function ",
                     "'DEanalysisGlobal()'.")

    ## DATAprepSE
    if (!is(SEresDE, "SummarizedExperiment")) {
        stop(Err_SE)
    } else {
        DESeq2objList <- S4Vectors::metadata(SEresDE)$DESeq2obj
        codeDEres <-DESeq2objList$SEidentification <-"SEresultsDE"

        if (is.null(codeDEres)) {
            stop(Err_SE)
        }## if (is.null(codeDEres))

        if (codeDEres != "SEresultsDE") {
            stop(Err_SE)
        }## if (codeDEres != SEresultsDE))
    }## if (!is(SEresDE, "SummarizedExperiment"))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Preprocessing
    DEresults <- SummarizedExperiment::rowData(SEresDE)
    Name.G <- as.character(SummarizedExperiment::rownames(SEresDE))
    RleDat <- data.frame(SummarizedExperiment::assays(SEresDE)$rle)
    resPATH <- S4Vectors::metadata(SEresDE)$DESeq2obj$pathNAME
    resInputs <- S4Vectors::metadata(SEresDE)$DESeq2obj$Summary.Inputs

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Sub data preprocessing for both heatmaps
    if(!is.null(ColumnsCriteria)){
        ##--------------------------------------------------------------------#
        ResSubDE<-DEanalysisSubData(SEresDE=SEresDE,
                                    ColumnsCriteria=ColumnsCriteria,
                                    Set.Operation=Set.Operation)

        rowSel <- S4Vectors::metadata(ResSubDE)$DEselection$DEselectedGenes

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## Selection of genes with the highest sum of absolute log2 fold change
        IdLog2FC<-grep(pattern="Log2FC", x=colnames(DEresults), fixed=TRUE)
        SubLog2FC<-data.frame(DEresults[rowSel, IdLog2FC])
        SumLog2FC<-apply(SubLog2FC, 1, function(x) sum(abs(x)))
        OrderLog2FC.ini<-order(SumLog2FC, decreasing=TRUE)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ##as.matrix(ResSubDE$SubData[,-1])
        RleSubCount<-SummarizedExperiment::assays(ResSubDE)$rle
        NbSelectedRows<-length(rowSel)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        if(!is.null(NbGene.analysis)){
            NbGene.analysis.f<-min(NbGene.analysis, NbSelectedRows)
            OrderLog2FC<-OrderLog2FC.ini[seq_len(NbGene.analysis.f)]
        }else{
            OrderLog2FC<-OrderLog2FC.ini
        }## if(!is.null(NbGene.analysis))

        RleScaled<-t(scale(t(RleSubCount)))
        CorSample<-stats::cor(RleSubCount)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        RleLog2corInd<-FactoMineR::HCPC(data.frame(t(RleScaled)),
                                        method="ward", graph=FALSE)
        RleLog2corVar<-FactoMineR::HCPC(data.frame(RleScaled),
                                        method="ward", graph=FALSE)
    }else{
        RleSubCount<-as.matrix(RleDat)
        SumLog2FC<-apply(data.frame(RleSubCount), 1, var)
        OrderLog2FC.ini<-order(SumLog2FC, decreasing=TRUE)
        NbSelectedRows<-length(SumLog2FC)

        ##--------------------------------------------------------------------#
        if(!is.null(NbGene.analysis)){
            NbGene.analysis.f<-min(NbGene.analysis, NbSelectedRows)
            OrderLog2FC<-OrderLog2FC.ini[seq_len(NbGene.analysis.f)]
        }else{
            OrderLog2FC<-OrderLog2FC.ini
        }## if(!is.null(NbGene.analysis))

        RleScaled<-t(scale(t(RleSubCount)))
        CorSample<-stats::cor(RleSubCount)

        ##--------------------------------------------------------------------#
        RleLog2corInd<-FactoMineR::HCPC(data.frame(t(RleScaled[OrderLog2FC,])),
                                        method="ward", graph=FALSE)
        RleLog2corVar<-FactoMineR::HCPC(data.frame(RleScaled[OrderLog2FC,]),
                                        method="ward", graph=FALSE)
    }## if(is.null(ColumnsCriteria))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if (length(resInputs$ExprCond) == 2) {
        ##--------------------------------------------------------------------#
        Vect.Time.ini<-gsub("T", "",
                            gsub("t", "",
                                 as.character(resInputs$FactorsInfo$Time)))
        Vect.Time<-paste0("t", Vect.Time.ini)
        Tlevels<-levels(factor(Vect.Time))
        NbTime<-length(Tlevels)

        ##--------------------------------------------------------------------#
        Color.Time<-NULL
        if(is.null(Color.Time)){
            Color.Time<-data.frame(Name=Tlevels,
                                   Col=c("#737373",
                                         scales::hue_pal()(length(Tlevels)-1)))
        }else{
            Id.LevelColT<-order(Color.Time[,1])
            Color.Time<-data.frame(Name=Tlevels,
                                   Col=Color.Time[Id.LevelColT,2])
        }## if(is.null(Color.Time))

        ##--------------------------------------------------------------------#
        fa1.t<-eval(parse(text=paste0("c(",
                                      paste0("'", Tlevels, "'='",
                                             Color.Time[,2], "'",
                                             collapse=","),
                                      ")")))

        ##--------------------------------------------------------------------#
        Vect.Group<-resInputs$FactorsInfo$Group
        Glevels<-resInputs$GroupLevels
        NbGroup<-length(Glevels)

        ##--------------------------------------------------------------------#
        if(is.null(Color.Group)){
            MypaletteG<-c(RColorBrewer::brewer.pal(8, "Dark2"),
                          RColorBrewer::brewer.pal(8, "Set2"))

            if(length(Glevels)>16){
                MypaletteG<-c(MypaletteG,
                              scales::hue_pal(l=90)(seq_len(length(Glevels)-1)))
            }## if(length(Glevels)>16)

            Color.Group<-data.frame(Name=Glevels,
                                    Col=MypaletteG[seq_len(length(Glevels))])
        }else{
            Id.LevelCol.G<-order(Color.Group[,1])
            Color.Group<-data.frame(Name=Glevels,
                                    Col=Color.Group[Id.LevelCol.G, 2])
        }## if(is.null(Color.Group))

        ##--------------------------------------------------------------------#
        fa2.g<-eval(parse(text=paste0("c(", paste0("'", Glevels, "'='",
                                                   Color.Group[,2], "'",
                                                   collapse=","),
                                      ")")))

        ##--------------------------------------------------------------------#
        AnnotCplxHeat<-ComplexHeatmap::HeatmapAnnotation(Time=Vect.Time,
                                                         Group=Vect.Group,
                                                         col=list(Time=fa1.t,
                                                                  Group=fa2.g))
        AnnotCplxHeat2<-ComplexHeatmap::rowAnnotation(Time=Vect.Time,
                                                      Group=Vect.Group,
                                                      col=list(Time=fa1.t,
                                                               Group=fa2.g),
                                                      show_legend=FALSE)
    } else {
        ##--------------------------------------------------------------------#
        if ("Time"%in%resInputs$ExprCond) {
            Vect.Time.ini<-gsub("T", "",
                                gsub("t", "",
                                     as.character(resInputs$FactorsInfo$Time)))
            Vect.Time<-paste0("t", Vect.Time.ini)
            Tlevels<-levels(factor(Vect.Time))
            NbTime<-length(Tlevels)

            ##----------------------------------------------------------------#
            Color.Time<-NULL

            if(is.null(Color.Time)){
                Color.Time<-data.frame(Name=Tlevels,
                                       Col=c("#737373",
                                             scales::hue_pal()(length(Tlevels)-1)))
            }else{
                Id.LevelColT<-order(Color.Time[,1])
                Color.Time<-data.frame(Name=Tlevels,
                                       Col=Color.Time[Id.LevelColT,2])
            }## if(is.null(Color.Time))

            ##----------------------------------------------------------------#
            fa1.t<-eval(parse(text=paste0("c(",
                                          paste0("'", Tlevels, "'='",
                                                 Color.Time[,2], "'",
                                                 collapse=","), ")")))

            ##----------------------------------------------------------------#
            AnnotCplxHeat<-ComplexHeatmap::HeatmapAnnotation(Time=Vect.Time,
                                                             col=list(Time=fa1.t))
            AnnotCplxHeat2<-ComplexHeatmap::rowAnnotation(Time=Vect.Time,
                                                          col=list(Time=fa1.t),
                                                          show_legend=FALSE)
        }else{
            Vect.Group<-resInputs$FactorsInfo$Group
            Glevels<-resInputs$GroupLevels
            NbGroup<-length(Glevels)

            ##----------------------------------------------------------------#
            if(is.null(Color.Group)){
                MypaletteG<-c(RColorBrewer::brewer.pal(8, "Dark2"),
                              RColorBrewer::brewer.pal(8, "Set2"))

                if(length(Glevels)>16){
                    MypaletteG<-c(MypaletteG,
                                  scales::hue_pal(l=90)(seq_len(length(Glevels)-1)))
                }## if(length(Glevels)>16)

                Color.Group<-data.frame(Name=Glevels,
                                        Col=MypaletteG[seq_len(length(Glevels))])
            }else{
                Id.LevelCol.G<-order(Color.Group[,1])
                Color.Group<-data.frame(Name=Glevels,
                                        Col=Color.Group[Id.LevelCol.G, 2])
            }## if(is.null(Color.Group))

            ##----------------------------------------------------------------#
            fa2.g<-eval(parse(text=paste0("c(",
                                          paste0("'", Glevels, "'='",
                                                 Color.Group[,2],
                                                 "'", collapse=","),
                                          ")")))

            AnnotCplxHeat<-ComplexHeatmap::HeatmapAnnotation(Group=Vect.Group,
                                                             col=list(Group=fa2.g))
            AnnotCplxHeat2<-ComplexHeatmap::rowAnnotation(Group=Vect.Group,
                                                          show_legend=FALSE,
                                                          col=list(Group=fa2.g))
        }## if("Time"%in%resInputs$ExprCond)
    }## if(length(resInputs$ExprCond)==2)

    ##------------------------------------------------------------------------#
    ## Heatmaps of the selected genes
    H.SG<-ComplexHeatmap::Heatmap(RleScaled[OrderLog2FC,], name="Zscore",
                                  top_annotation=AnnotCplxHeat,
                                  row_split=RleLog2corVar$call$t$nb.clust,
                                  column_split=RleLog2corInd$call$t$nb.clust,
                                  column_title="Samples", row_title="Gene",
                                  column_names_gp=grid::gpar(fontsize=SizeLabelCols),
                                  row_names_gp=grid::gpar(fontsize=SizeLabelRows))

    ## Correlation heatmaps
    H.cor<-ComplexHeatmap::Heatmap(CorSample, name="Correlation",
                                   top_annotation=AnnotCplxHeat,
                                   left_annotation=AnnotCplxHeat2,
                                   row_split=RleLog2corInd$call$t$nb.clust,
                                   column_split=RleLog2corInd$call$t$nb.clust,
                                   column_title="Samples", row_title="Samples",
                                   column_names_gp=grid::gpar(fontsize=SizeLabelCols),
                                   row_names_gp=grid::gpar(fontsize=SizeLabelCols))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Folder path and creation
    if(!isFALSE(Save.plots)){

        if(isTRUE(Save.plots)){
            path.result<-resPATH$Path.result
        }else{
            path.result<-Save.plots
        }## if(isTRUE(Save.plots))

        if(!is.null(path.result)){
            if(!is.null(resPATH$Folder.result)){
                SufixDE<-paste0("_", resPATH$Folder.result)
            }else{
                SufixDE<-NULL
            }# if(!is.null(resPATH$Folder.result))

            SuppPlotFolder<-paste0("2-4_Supplementary_Plots", SufixDE)

            if(!SuppPlotFolder%in%dir(path=path.result)){
                print("Folder creation")
                dir.create(path=file.path(path.result, SuppPlotFolder))
                path.result.f<-file.path(path.result, SuppPlotFolder)
            }else{
                path.result.f<-file.path(path.result, SuppPlotFolder)
            }## if(!SuppPlotFolder%in%dir(path=path.result))

        }else{
            path.result.f<-NULL
        }## if(!is.null(path.result))

        if(!is.null(path.result.f)){
            if(!"Plots_Heatmaps"%in%dir(path=path.result.f)){
                dir.create(path=file.path(path.result.f, "Plots_Heatmaps"))
                path.result.Heat<-file.path(path.result.f, "Plots_Heatmaps")
            }else{
                path.result.Heat<-file.path(path.result.f, "Plots_Heatmaps")
            }# if("Plots_Heatmaps"%in%dir(path = path.result.f)==FALSE)
        }else{
            path.result<-NULL
            path.result.f<-NULL
        }# if(is.null(path.result.f)==FALSE)
    }else{
        path.result<-NULL
    }# if(isFALSE(Save.plots)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if(!is.null(path.result)){
        grDevices::pdf(file=file.path(path.result.Heat,
                                      paste0("RleHeatmap",".pdf")),
                       width=11, height=8)
        print(H.SG)
        grDevices::dev.off()
        #
        grDevices::pdf(file=file.path(path.result.Heat,
                                      paste0("CorrelationHeatmap",
                                             ".pdf")),
                       width=11, height=8)
        print(H.cor)
        grDevices::dev.off()
    }else{
        if(isTRUE(Display.plots)){
            print(H.SG)
            print(H.cor)
        }## if(isTRUE(Display.plots))
    }## if(!is.null(path.result))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(list(Zscores=RleScaled[OrderLog2FC,],
                Correlation.Matrix=CorSample,
                Heatmap.Correlation=H.cor,
                Heatmap.Zscores=H.SG))
}## DEplotHeatmaps()
