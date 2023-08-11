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
#' * If \code{Set.Operation="union"} then the rows (so genes) extracted from
#' \code{Res.DE.analysis$DE.results} are those such that the sum of
#' the selected columns in \code{Res.DE.analysis} is >0.
#' For example, the rows extracted from \code{Res.DE.analysis$DE.results}
#' will be those DE at least at one time.
#'
#' * If \code{Set.Operation="intersect"} then the rows extracted from
#' \code{Res.DE.analysis$DE.results} are those such that the product of
#' the selected columns in \code{Res.DE.analysis$DE.results} is >0.
#' For example, the rows extracted from \code{Res.DE.analysis$DE.results}
#' will be those DE at all time ti (except the reference time t0).
#'
#' * If \code{Set.Operation="setdiff"} then the rows extracted from data
#' are those such that only one element of the selected columns in
#' \code{Res.DE.analysis$DE.results} is >0.
#' For example, the rows extracted from \code{Res.DE.analysis$DE.results}
#' will be those DE at only one time ti (except the reference time t0).
#'
#' @param Internect.Connection \code{TRUE} or \code{FALSE}.
#' \code{FALSE} as default.
#' If \code{TRUE}, the function realizes an enrichment analysis.
#' If \code{FALSE}, the function returns a message indicating that
#' the user must have an internet connection.
#' @param Res.DE.analysis A list corresponding to the output of
#' [DEanalysisGlobal()].
#' @param ColumnsCriteria A vector of integers where each integer indicates
#' a column of \code{Res.DE.analysis$DE.results} to be selected for
#' GSEA analysis. These columns should either contain only binary values,
#' or may contain other numerical value, in which case extracted rows
#' from \code{Res.DE.analysis$DE.results} will be those with >0 values
#' (see \code{Details}).
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
#' @return The function returns
#' * a data.frame which contains the outputs of
#' [gprofiler2::gost()]
#' * a Manhattan plot showing all GO names according to their pvalue.
#' * A lollipop graph showing the \code{MaxNumberGO} most important GO.
#'
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
#' data(Results_DEanalysis_sub500)
#' ## Results of DEanalysisGlobal() with the dataset of Antoszewski
#' res.all<-Results_DEanalysis_sub500$DE_Antoszewski2022_MOUSEsub500
#'
#' ##-------------------------------------------------------------------------#
#' ## Internet is needed in order to run the following lines of code
#' ## resGpA<-GSEAQuickAnalysis(Internect.Connection=FALSE,
#' ##                           Res.DE.analysis=res.all,
#' ##                           ColumnsCriteria=2,
#' ##                           ColumnsLog2ordered=NULL,
#' ##                           Set.Operation="union",
#' ##                           Organism="mmusculus",
#' ##                           MaxNumberGO=20,
#' ##                           Background=FALSE,
#' ##                           Display.plots=TRUE,
#' ##                           Save.plots=FALSE)
#' ##-------------------------------------------------------------------------#
#'
#' ##-------------------------------------------------------------------------#
#' ## The results res.all of DEanalysisGlobal with the dataset Antoszewski2022
#' ## data(RawCounts_Antoszewski2022_MOUSEsub500)
#' ## res.all<-DEanalysisGlobal(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
#' ##                           Column.gene=1, Group.position=1,
#' ##                           Time.position=NULL, Individual.position=2,
#' ##                           pval.min=0.05, log.FC.min=1,LRT.supp.info=FALSE,
#' ##                           path.result=NULL, Name.folder.DE=NULL)

GSEAQuickAnalysis<-function(Internect.Connection=FALSE,
                            Res.DE.analysis,
                            ColumnsCriteria=1,
                            ColumnsLog2ordered=NULL,
                            Set.Operation="union",
                            Organism="hsapiens",
                            MaxNumberGO=20,
                            Background=FALSE,
                            Display.plots=TRUE,
                            Save.plots=FALSE){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if(!is.list(Res.DE.analysis) & !is(Res.DE.analysis, 'DESeqDataSet')){
        stop("Res.DE.analysis must be a list or a 'DESeqDataSet' object")
    }## if(!is.list(Res.DE.analysis) & !is(classDeseq2, 'DESeqDataSet'))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if(isFALSE(Internect.Connection)){
        GOmat<-gManhattan<-glolipop<-NULL

        MessageInternet<-paste("Once the user is sure to have",
                               "an internet connection, the user must set",
                               "'Internect.Connection=TRUE' in order to realize",
                               "the enrichment analysis", sep=" ")
        message(MessageInternet)
    }else{
        ##--------------------------------------------------------------------#
        ## To avoid "no visible binding for global variable" with
        ## devtools::check()
        term_id<-p_value<-significant<-xmh<-NULL

        ##--------------------------------------------------------------------#
        ## Folder path and creation
        if(!isFALSE(Save.plots)){
            if(Save.plots==TRUE){
                path.result<-Res.DE.analysis$Path.result
            }else{
                path.result<-Save.plots
            }## if(Save.plots==TRUE)

            if(!is.null(path.result)){
                if(!is.null(Res.DE.analysis$Folder.result)){
                    SufixDE<-paste0("_", Res.DE.analysis$Folder.result)
                }else{
                    SufixDE<-NULL
                }## if(is.null(Res.DE.analysis$Folder.result)==FALSE)

                if(Save.plots==TRUE){
                    SuppPlotFolder<-paste0("2-5_Enrichment_analysis", SufixDE)

                    if(!SuppPlotFolder%in%dir(path=path.result)){
                        print("Folder creation")
                        dir.create(path=file.path(path.result, SuppPlotFolder))
                        path.result.f<-file.path(path.result, SuppPlotFolder)
                    }else{
                        path.result.f<-file.path(path.result, SuppPlotFolder)
                    }## if(SuppPlotFolder%in%dir(path = path.result)==FALSE)

                }else{
                    path.result.f<-path.result
                }## if(Save.plots==TRUE)

            }else{
                path.result.f<-NULL
            }## if(is.null(path.result)==FALSE)

            if(!is.null(path.result)){
                if(Save.plots==TRUE){
                    nom.dossier.result1<-paste0("2-5-1_gprofiler2_results",
                                                SufixDE)
                }else{
                    nom.dossier.result1<-paste0("gprofiler2_results",
                                                SufixDE)
                }

                if(!nom.dossier.result1%in%dir(path=path.result.f)){
                    dir.create(path=file.path(path.result.f,
                                              nom.dossier.result1))
                    path.result.GSEA<-file.path(path.result.f,
                                                nom.dossier.result1)
                }else{
                    path.result.GSEA<-file.path(path.result.f,
                                                nom.dossier.result1)
                }## if(nom.dossier.result1%in%dir(path = path.result.f)==FALSE)

                LollipopPlot<-file.path(path.result.GSEA, "LollipopChart.pdf")
                ManhPlot<-file.path(path.result.GSEA, "ManhattanPlot.pdf")
                GSEAtable<-file.path(path.result.GSEA, "GSEAtable.csv")
                GSEAtablepng<-file.path(path.result.GSEA, "GSEAtable.pdf")
            }else{
                ## To avoid "no visible binding for global variable"
                ## with devtools::check()
                LollipopPlot<-GSEAtablepng<-GSEAtable<-ManhPlot<-NULL
                path.result.GSEA<-path.result<-NULL
                ## path.result.GSEA<-NULL
            }## if(is.null(path.result)==FALSE)

        }else{
            ## To avoid "no visible binding for global variable" with
            ## devtools::check()
            LollipopPlot<-GSEAtablepng<-GSEAtable<-ManhPlot<-NULL
            path.result.GSEA<-path.result<-NULL
        }## if(isFALSE(Save.plots)==FALSE)

        ##--------------------------------------------------------------------#
        ## Selection of genes according to Columns criteria and Set.operation
        ResSubDE<-DEanalysisSubData(Data=Res.DE.analysis$RLEdata,
                                    Res.DE.analysis=Res.DE.analysis,
                                    ColumnsCriteria=ColumnsCriteria,
                                    Set.Operation=Set.Operation)

        GeneNames<-Res.DE.analysis$RLEdata$Gene

        ##--------------------------------------------------------------------#
        ## Gene order considered as important or not
        if(!is.null(ColumnsLog2ordered)){
            Log2FCchosen<-Res.DE.analysis$DE.results[ResSubDE$RowsSelected,
                                                     ColumnsLog2ordered]

            orderLog2FC<-order(abs(Log2FCchosen))
            GeneSelected.i<-GeneNames[ResSubDE$RowsSelected]
            GeneSelected<-GeneSelected.i[orderLog2FC]
            QueryOrder<-TRUE
        }else{
            GeneSelected<-GeneNames[ResSubDE$RowsSelected]
            QueryOrder<-FALSE
        }## if(is.null(ColumnsLog2ordered)==FALSE)

        ##--------------------------------------------------------------------#
        ## Inputs of gprofiler2::gost()
        if(isFALSE(Background)){
            Backgene<-NULL
            DomainScoped<-"annotated"
        }else{
            Backgene<-GeneNames
            DomainScoped<-"custom"
        }# if(Background==FALSE)

        ##--------------------------------------------------------------------#
        ## gprofiler2 analysis
        gostres<-gprofiler2::gost(query=GeneSelected,
                                  organism=Organism,
                                  ordered_query=QueryOrder,
                                  multi_query=FALSE,
                                  significant=FALSE,
                                  exclude_iea=FALSE,
                                  measure_underrepresentation=FALSE,
                                  evcodes=TRUE,
                                  user_threshold=0.05,
                                  correction_method="g_SCS",
                                  domain_scope=DomainScoped,
                                  custom_bg=Backgene,
                                  numeric_ns="",
                                  sources=c("GO:BP", "GO:MF", "GO:CC", "KEGG"),
                                  as_short_link=FALSE)
        # term_size - number of genes that are annotated to the term (GO)
        # query_sizes
        ### - number of genes that were included in the query in the order of
        ###   input queries
        # intersection_sizes
        ### - the number of genes in the input query that are annotated to
        ###   the corresponding term in the order of input queries
        # precision
        ### - the proportion of genes in the input list that are annotated to
        ###   the function (defined as intersection_size/query_size)
        # recall
        ### - the proportion of functionally annotated genes that the query
        ###   recovers (defined as intersection_size/term_size)

        ##--------------------------------------------------------------------#
        ## Creation of GOmat, a data.frame containing result
        orderdPva<-order(gostres$result$p_value)
        NewgoRes<-gostres$result[orderdPva,]
        row.names(NewgoRes)<-NULL
        gostres$result<-NewgoRes
        GOmat<-NewgoRes[,c(9:11, 3:2, 4:8)]

        ##--------------------------------------------------------------------#
        ## Preprocessing data.frame for graphs
        GeneID<-paste0("(", NewgoRes$intersection, ")")
        GOmat$Gene_id<-gsub(",",")_(",GeneID, fixed=TRUE)
        GOmat$GOparents<-unlist(lapply(NewgoRes$parents,
                                       function(x) paste0("(",
                                                          paste(x,
                                                                collapse=")_("),
                                                          ")")))

        GOmat$term_id<-factor(GOmat$term_id, levels=rev(GOmat$term_id))
        GOmat$significant<-factor(GOmat$significant,
                                  levels=c("TRUE", "FALSE"))
        GOmat$source<-factor(GOmat$source,
                             levels=c("GO:BP", "GO:CC", "GO:MF" ,"KEGG"))

        if(max(GOmat$p_value)==1){
            GOmat<-GOmat[-which(GOmat$p_value==1),]
        }## if(max(GOmat$p_value)==1)

        ##--------------------------------------------------------------------#
        if(!is.null(path.result) | Display.plots==TRUE){
            ## data.frame for lolipop graph
            DatLolipop<-GOmat[seq_len(min(c(MaxNumberGO, nrow(GOmat)))),
                              c(1, 2, 4, 5)]

            ## Graph option
            IntegerMaxLogPval<-ceiling(max(-log10(DatLolipop$p_value)))
            SizeGOlabel<-(0.5*4*15)/nrow(DatLolipop)

            ##----------------------------------------------------------------#
            LevelSignificant<-which(levels(DatLolipop$significant)%in%unique(DatLolipop$significant)==TRUE)

            ##----------------------------------------------------------------#
            ## Lolipop graph
            glolipop<-ggplot2::ggplot(DatLolipop,
                                      ggplot2::aes(x=term_id,
                                                   y=-log10(p_value),
                                                   ymin=0,
                                                   ymax=IntegerMaxLogPval))+
                ggplot2::geom_point(ggplot2::aes(colour=significant),
                                    size=1.9,
                                    shape=19)+
                ggplot2::coord_flip() +
                ggplot2::geom_segment(ggplot2::aes(x=term_id, xend=term_id,
                                                   y=-log10(p_value), yend=0,
                                                   colour=significant),
                                      size=1) +
                ggplot2::geom_bar(stat="identity", color="black", size=0.5,
                                  ggplot2::aes(fill=source,
                                               y=-0.10*(IntegerMaxLogPval*0.55))) +
                ggplot2::geom_hline(yintercept=-log10(0.05),
                                    linetype="dashed",
                                    size=0.6)+
                ggplot2::scale_color_manual(values=c("#E69F00",
                                                     "#56B4E9")[LevelSignificant])+
                ggplot2::scale_fill_manual(values=c("#CB2027","#059748",
                                                    "#EE7733","#7B3A96"),
                                           name="Source")+
                ggplot2::xlab("") + ggplot2::ylab("-log10(pvalue)")+
                ggplot2::ggtitle("gprofiler2 enrichments analysis")+
                ggplot2::labs(color = "Significant\n  (<0.05)")+
                ggplot2::theme_bw()+#theme_linedraw()+
                ggplot2::theme(legend.position="right",#"bottom",
                               legend.box="vertical",
                               axis.text.y=ggplot2::element_text(face="bold",
                                                                 size=ggplot2::rel(SizeGOlabel)),
                               legend.title=ggplot2::element_text(face="bold",
                                                                  size=ggplot2::rel(0.8)),
                               legend.text=ggplot2::element_text(size=ggplot2::rel(0.6)),
                               legend.key.size = ggplot2::unit(0.4, "cm"))

            ##----------------------------------------------------------------#
            ## data.frame for Manhattan graph
            if(length(which(GOmat$p_value==1))>0){
                DatManhattan<-GOmat[-which(GOmat$p_value==1), c(1, 2, 4, 5)]
            }else{
                DatManhattan<-GOmat[,c(1,2,4,5)]
            }## if(length(which(GOmat$p_value==1))>0)

            DatManhattan<-DatManhattan[order(DatManhattan$source,
                                             DatManhattan$p_value),]
            DatManhattan$xmh<-seq_len(nrow(DatManhattan))
            DatManhattan$Color<-factor(DatManhattan$source,
                                       levels=c("GO:BP", "GO:CC",
                                                "GO:MF", "KEGG"))

            levels(DatManhattan$Color)<-c("#CB2027", "#059748",
                                          "#EE7733", "#7B3A96")
            DatManhattan$Point<-DatManhattan$significant
            levels(DatManhattan$Point)<-c(8, 16)

            ##----------------------------------------------------------------#
            ## Option for the Manhattan graph
            NbPerGOw<-as.numeric(table(DatManhattan$source))
            if(length(which(NbPerGOw==0))>0){
                NbPerGOw<-NbPerGOw[-which(NbPerGOw == 0)]
            }## if(length(which(NbPerGOw==0))>0)

            xGO<-cumsum(NbPerGOw)+1-((NbPerGOw+1)/2)
            SelMh<-seq_len(max(ceiling(nrow(DatManhattan)*0.05),10))
            ## 1:max(ceiling(nrow(DatManhattan)*0.05),10)
            ## min(c(MaxNumberGO,nrow(DatManhattan)))
            LevelSignificantManh<-which(levels(DatManhattan$significant)%in%unique(DatManhattan$significant)==TRUE)

            ##----------------------------------------------------------------#
            ## Manhattan plot
            gManhattan<-ggplot2::ggplot(DatManhattan,
                                        ggplot2::aes(x=xmh,
                                                     y=-log10(p_value),
                                                     label=term_id,
                                                     ymin=0,
                                                     ymax=IntegerMaxLogPval))+
                ggplot2::geom_point(ggplot2::aes(color=source,
                                                 shape=significant),
                                    size=1.7) + #alpha=0.8,
                ggplot2::geom_line(ggplot2::aes(color=as.character(source)),
                                   size=0.4)+
                ggplot2::geom_hline(yintercept=-log10(0.05),
                                    linetype="dashed",
                                    size=0.6)+
                ggplot2::geom_area(ggplot2::aes(fill=as.character(source),
                                                group=as.character(source)),
                                   alpha=0.2, position='identity')+
                ggplot2::geom_bar(stat="identity", size=1.5,
                                  ggplot2::aes(fill=source,color=source,
                                               y=-0.06*(IntegerMaxLogPval*0.55))) +
                ggplot2::scale_fill_manual(values=c("#CB2027", "#059748",
                                                    "#EE7733", "#7B3A96"),
                                           name="Source")+
                ggplot2::scale_color_manual(values=c("#CB2027", "#059748",
                                                     "#EE7733", "#7B3A96"),
                                            name="Source")+
                ggplot2::scale_shape_manual(values=c(8,
                                                     16)[LevelSignificantManh],
                                            name="Significant\n   (<0.05)")+
                ggrepel::geom_label_repel(data=DatManhattan[order(DatManhattan$p_value),][SelMh,],
                                          col=DatManhattan$Color[order(DatManhattan$p_value)][SelMh],
                                          box.padding=0.2,
                                          max.overlaps=nrow(DatManhattan),
                                          min.segment.length=0, size=2)+
                ggplot2::scale_x_continuous(label=unique(DatManhattan$source),
                                            breaks=xGO,
                                            guide=ggplot2::guide_axis(angle=90))+
                ggplot2::scale_y_continuous(limits=c(-0.06*(IntegerMaxLogPval*0.55),
                                                     IntegerMaxLogPval),
                                            expand=c(0,0))+
                ggplot2::xlab("")+ ggplot2::ylab("-log10(pvalue)")+
                ggplot2::ggtitle("Manhattan plot")+
                ggplot2::theme_classic()+#theme_bw() +
                ggplot2::theme(legend.position="right",# legend.position="bottom",
                               legend.box="vertical",# legend.box = "vertical",
                               panel.border = ggplot2::element_blank(),
                               panel.grid.major.x = ggplot2::element_blank(),
                               panel.grid.minor.x = ggplot2::element_blank(),
                               legend.title=ggplot2::element_text(face="bold",
                                                                  size=ggplot2::rel(0.8)),
                               legend.text=ggplot2::element_text(size=ggplot2::rel(0.7)),
                               legend.key.size = ggplot2::unit(0.5, "cm"))
            ##
            ## p<-gprofiler2::gostplot(gostres, capped = TRUE, interactive = FALSE)
        }else{
            gManhattan<-glolipop<-NULL
        }

        if(!is.null(path.result)){
            ## grDevices::png(file = ManhPlot, width = 800, height = 700)
            grDevices::pdf(file=ManhPlot, width=11, height=8)
            print(gManhattan)
            grDevices::dev.off()

            ## grDevices::png(file = LollipopPlot, width = 800, height = 700)
            grDevices::pdf(file=LollipopPlot, width=11, height=8)
            print(glolipop)
            grDevices::dev.off()

            utils::write.table(GOmat, file=GSEAtable, sep=";",row.names=FALSE)

            if(Display.plots==TRUE){
                print(gManhattan)
                print(glolipop)
            }

        }else{

            if(Display.plots==TRUE){
                print(gManhattan)
                print(glolipop)
            }

        }## if(is.null(path.result)==FALSE)
    }## if(Internect.Connection==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(list(GSEAresults=GOmat,
                LollipopChart=glolipop,
                ManhattanPlot=gManhattan))
}## GSEAQuickAnalysis()
