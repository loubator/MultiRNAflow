#' @title Heatmaps of DE genes
#'
#' @description The function returns two heatmaps:
#' one heatmap of gene expressions between samples and selected genes and
#' a correlation heatmap between samples
#' from the output of [DEanalysisGlobal()].
#'
#' @details
#' * If \code{Set.Operation="union"} then the rows extracted from \code{Data}
#' are those such that the sum of the selected columns in \code{Res.DE.analysis}
#' is >0.
#' For example, if \code{Res.DE.analysis} is the outputs from
#' [DEanalysisGlobal()], the rows extracted from \code{Data} will be those DE
#' at least at one time.
#'
#' * If \code{Set.Operation="intersect"} then the rows extracted from data are
#' those such that the product of the selected columns in \code{Res.DE.analysis}
#' is >0.
#' For example, if \code{Res.DE.analysis} is the outputs from
#' [DEanalysisGlobal()], the rows extracted from \code{Data} will be those DE
#' at all time ti (except the reference time t0).
#'
#' * If \code{Set.Operation="setdiff"} then the rows extracted from data are
#' those such that only one element of the selected columns in
#' \code{Res.DE.analysis} is >0.
#' For example, if \code{Res.DE.analysis} is the outputs from
#' [DEanalysisGlobal()], the rows extracted from \code{Data} will be those DE
#' at only one time ti (except the reference time t0).
#'
#' The time in the first column of \code{Color.Time} must be either
#' 't0', 'T0' or '0' for time 0,
#' 't1', 'T1' or '1' for time 1, ...
#'
#' @param Res.DE.analysis A list. Output from [DEanalysisGlobal()].
#' @param ColumnsCriteria A vector of integers where each integer indicates
#' column of the output \code{DE.results} from \code{Res.DE.analysis}.
#' Columns should contained binary values, but if there are other
#' integers or other values, rows selected are those which values are >0
#' (see \code{Details}).
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
#' @param Save.plots TRUE or FALSE or a Character.
#' If \code{Save.plots=FALSE}, the different files will not be saved.
#' If \code{Save.plots=TRUE} and the \code{path.result} of [DEanalysisGlobal()]
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
#' @seealso The function calls the function [ComplexHeatmap::Heatmap()]
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
#' data(Results_DEanalysis_sub500)
#' # Results of DEanalysisGlobal() with the dataset of Antoszewski
#' res.all<-Results_DEanalysis_sub500$DE_Antoszewski2022_MOUSEsub500
#' #
#' DEplotHeatmaps(Res.DE.analysis=res.all,
#'                ColumnsCriteria=2,
#'                Set.Operation="union",
#'                NbGene.analysis=20,
#'                Color.Group=NULL,
#'                SizeLabelRows=5,
#'                SizeLabelCols=5,
#'                Save.plots=FALSE)
#' #
#' #--------------------------------------------------------------------------#
#' ## The results res.all of DEanalysisGlobal with the dataset Antoszewski2022
#' # data(RawCounts_Antoszewski2022_MOUSEsub500)
#' # res.all<-DEanalysisGlobal(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
#' #                           Column.gene=1, Group.position=1,
#' #                           Time.position=NULL, Individual.position=2,
#' #                           pval.min=0.05, log.FC.min=1,LRT.supp.info=FALSE,
#' #                           path.result=NULL, Name.folder.DE=NULL)

DEplotHeatmaps<-function(Res.DE.analysis,
                         ColumnsCriteria=2,
                         Set.Operation="union",
                         NbGene.analysis=20,
                         Color.Group=NULL,
                         SizeLabelRows=5,
                         SizeLabelCols=5,
                         Save.plots=FALSE){
  #---------------------------------------------------------------------------#
  # RLE count data
  RleDat<-Res.DE.analysis$List.Datas$RLEdata
  # Sub data preprocessing for both heatmaps
  if(is.null(ColumnsCriteria)==FALSE){
    ResSubDE<-DEanalysisSubData(Data=RleDat,
                                Res.DE.analysis=Res.DE.analysis,
                                ColumnsCriteria=ColumnsCriteria,
                                Set.Operation=Set.Operation)
    #-------------------------------------------------------------------------#
    # Selection of genes with the highest sum of absolute log2 fold change
    IdLog2FC<-grep(pattern="Log2FC", x=colnames(Res.DE.analysis$DE.results),
                   fixed=TRUE)
    SubLog2FC<-data.frame(Res.DE.analysis$DE.results[ResSubDE$RowsSelected,
                                                     IdLog2FC])
    SumLog2FC<-apply(SubLog2FC, 1, function(x) sum(abs(x)))
    OrderLog2FC.ini<-order(SumLog2FC, decreasing=TRUE)
    #-------------------------------------------------------------------------#
    RleSubCount<-as.matrix(ResSubDE$SubData[,-1])
    NbSelectedRows<-length(ResSubDE$RowsSelected)
    #-------------------------------------------------------------------------#
    if(is.null(NbGene.analysis)==FALSE){
      NbGene.analysis.f<-min(NbGene.analysis, NbSelectedRows)
      OrderLog2FC<-OrderLog2FC.ini[seq_len(NbGene.analysis.f)]
    }else{# seq_len(NbGene.analysis.f) == 1:NbGene.analysis.f
      OrderLog2FC<-OrderLog2FC.ini
    }# if(is.null(NbGene.analysis)==FALSE)
    #
    RleScaled<-t(scale(t(RleSubCount)))
    CorSample<-stats::cor(RleSubCount)
    #-------------------------------------------------------------------------#
    RleLog2corInd<-FactoMineR::HCPC(data.frame(t(RleScaled)),
                                    method="ward", graph=FALSE)
    RleLog2corVar<-FactoMineR::HCPC(data.frame(RleScaled),
                                    method="ward", graph=FALSE)
  }else{
    RleSubCount<-as.matrix(RleDat[,-1])
    SumLog2FC<-apply(data.frame(RleSubCount), 1, var)
    OrderLog2FC.ini<-order(SumLog2FC, decreasing=TRUE)
    NbSelectedRows<-length(SumLog2FC)
    #-------------------------------------------------------------------------#
    if(is.null(NbGene.analysis)==FALSE){
      NbGene.analysis.f<-min(NbGene.analysis, NbSelectedRows)
      OrderLog2FC<-OrderLog2FC.ini[seq_len(NbGene.analysis.f)]
    }else{ # [seq_len(NbGene.analysis.f) == 1:NbGene.analysis.f
      OrderLog2FC<-OrderLog2FC.ini
    }# if(is.null(NbGene.analysis)==FALSE)
    #
    RleScaled<-t(scale(t(RleSubCount)))
    CorSample<-stats::cor(RleSubCount)
    #-------------------------------------------------------------------------#
    RleLog2corInd<-FactoMineR::HCPC(data.frame(t(RleScaled[OrderLog2FC,])),
                                    method="ward", graph=FALSE)
    RleLog2corVar<-FactoMineR::HCPC(data.frame(RleScaled[OrderLog2FC,]),
                                    method="ward", graph=FALSE)
  }# if(is.null(ColumnsCriteria)==FALSE)
  #---------------------------------------------------------------------------#
  if(length(Res.DE.analysis$Summary.Inputs$ExprCond)==2){
    Vect.Time.ini<-gsub("T","",
                        gsub("t","",
                             as.character(Res.DE.analysis$Summary.Inputs$FactorsInfo$Time)))
    Vect.Time<-paste("t", Vect.Time.ini, sep="")
    Tlevels<-levels(factor(Vect.Time))
    NbTime<-length(Tlevels)
    #
    Color.Time<-NULL
    if(is.null(Color.Time)==TRUE){
      Color.Time<-data.frame(Name=Tlevels,
                             Col=c("#737373",
                                   scales::hue_pal()(length(Tlevels)-1)))
    }else{
      Id.LevelColT<-order(Color.Time[,1])
      Color.Time<-data.frame(Name=Tlevels, Col=Color.Time[Id.LevelColT,2])
    }# if(is.null(Color.Time)==TRUE)
    #
    fa1.t<-eval(parse(text=paste("c(",
                                 paste("'", Tlevels, "'='", Color.Time[,2], "'",
                                       sep="", collapse=","),
                                 ")", sep="")))
    #
    Vect.Group<-Res.DE.analysis$Summary.Inputs$FactorsInfo$Group
    Glevels<-Res.DE.analysis$Summary.Inputs$GroupLevels
    NbGroup<-length(Glevels)
    #
    if(is.null(Color.Group)==TRUE){
      MypaletteG<-c(RColorBrewer::brewer.pal(8,"Dark2"),
                    RColorBrewer::brewer.pal(8,"Set2"))
      #
      if(length(Glevels)>16){
        MypaletteG<-c(MypaletteG,
                      hue_pal(l=90)(seq_len(length(Glevels)-1)))
      }# if(length(Glevels)>16)
      #
      Color.Group<-data.frame(Name=Glevels,
                              Col=MypaletteG[seq_len(length(Glevels))])
    }else{
      Id.LevelCol.G<-order(Color.Group[,1])
      Color.Group<-data.frame(Name=Glevels, Col=Color.Group[Id.LevelCol.G,2])
    }# if(is.null(Color.Group)==TRUE)
    #
    fa2.g<-eval(parse(text=paste("c(", paste("'", Glevels, "'='",
                                             Color.Group[,2], "'",
                                             sep="",collapse=","),
                                 ")", sep="")))
    #
    AnnotCplxHeat<-ComplexHeatmap::HeatmapAnnotation(Time=Vect.Time,
                                                     Group=Vect.Group,
                                                     col=list(Time=fa1.t,
                                                              Group=fa2.g))
    AnnotCplxHeat2<-ComplexHeatmap::rowAnnotation(Time=Vect.Time,
                                                  Group=Vect.Group,
                                                  col=list(Time=fa1.t,
                                                           Group=fa2.g),
                                                  show_legend=FALSE)
  }else{
    if("Time"%in%Res.DE.analysis$Summary.Inputs$ExprCond){
      Vect.Time.ini<-gsub("T","",
                          gsub("t","",
                               as.character(Res.DE.analysis$Summary.Inputs$FactorsInfo$Time)))
      Vect.Time<-paste("t",Vect.Time.ini,sep="")
      Tlevels<-levels(factor(Vect.Time))
      NbTime<-length(Tlevels)
      #
      Color.Time<-NULL
      if(is.null(Color.Time)==TRUE){
        Color.Time<-data.frame(Name=Tlevels,
                               Col=c("#737373",
                                     scales::hue_pal()(length(Tlevels)-1)))
      }else{
        Id.LevelColT<-order(Color.Time[,1])
        Color.Time<-data.frame(Name=Tlevels, Col=Color.Time[Id.LevelColT,2])
      }# if(is.null(Color.Time)==TRUE)
      #
      fa1.t<-eval(parse(text=paste("c(",
                                   paste("'",Tlevels,"'='", Color.Time[,2], "'",
                                         sep="", collapse=","), ")", sep="")))
      #
      AnnotCplxHeat<-ComplexHeatmap::HeatmapAnnotation(Time=Vect.Time,
                                                       col=list(Time=fa1.t))
      AnnotCplxHeat2<-ComplexHeatmap::rowAnnotation(Time=Vect.Time,
                                                    col=list(Time=fa1.t),
                                                    show_legend=FALSE)
    }else{
      Vect.Group<-Res.DE.analysis$Summary.Inputs$FactorsInfo$Group
      Glevels<-Res.DE.analysis$Summary.Inputs$GroupLevels
      NbGroup<-length(Glevels)
      #
      if(is.null(Color.Group)==TRUE){
        MypaletteG<-c(RColorBrewer::brewer.pal(8,"Dark2"),
                      RColorBrewer::brewer.pal(8,"Set2"))
        #
        if(length(Glevels)>16){
          MypaletteG<-c(MypaletteG,
                        hue_pal(l=90)(seq_len(length(Glevels)-1)))
        }# if(length(Glevels)>16)
        #
        Color.Group<-data.frame(Name=Glevels,
                                Col=MypaletteG[seq_len(length(Glevels))])
      }else{
        Id.LevelCol.G<-order(Color.Group[,1])
        Color.Group<-data.frame(Name=Glevels, Col=Color.Group[Id.LevelCol.G,2])
      }# if(is.null(Color.Group)==TRUE)
      #
      fa2.g<-eval(parse(text=paste("c(",
                                   paste("'", Glevels, "'='", Color.Group[,2],
                                         "'", sep="",collapse=","),
                                   ")", sep="")))
      #
      AnnotCplxHeat<-ComplexHeatmap::HeatmapAnnotation(Group = Vect.Group,
                                                       col=list(Group = fa2.g))
      AnnotCplxHeat2<-ComplexHeatmap::rowAnnotation(Group = Vect.Group,
                                                    show_legend = FALSE,
                                                    col = list(Group = fa2.g))
    }# if("Time"%in%Res.DE.analysis$Summary.Inputs$ExprCond)
  }# if(length(Res.DE.analysis$Summary.Inputs$ExprCond)==2)
  #---------------------------------------------------------------------------#
  # Heatmaps of the selected genes
  H.SG<-ComplexHeatmap::Heatmap(RleScaled[OrderLog2FC,], name="Zscore",
                                top_annotation=AnnotCplxHeat,
                                row_split=RleLog2corVar$call$t$nb.clust,
                                column_split=RleLog2corInd$call$t$nb.clust,
                                column_title="Samples", row_title="Gene",
                                column_names_gp=grid::gpar(fontsize=SizeLabelCols),
                                row_names_gp=grid::gpar(fontsize=SizeLabelRows))
  # heatmap_legend_param = list(title_position = "leftcenter-rot"),
  # Correlation heatmaps
  H.cor<-ComplexHeatmap::Heatmap(CorSample, name="Correlation",
                                 top_annotation=AnnotCplxHeat,
                                 left_annotation=AnnotCplxHeat2,
                                 row_split=RleLog2corInd$call$t$nb.clust,
                                 column_split=RleLog2corInd$call$t$nb.clust,
                                 column_title="Samples", row_title="Samples",
                                 column_names_gp=grid::gpar(fontsize=SizeLabelCols),
                                 row_names_gp=grid::gpar(fontsize=SizeLabelCols))
  #---------------------------------------------------------------------------#
  # Folder path and creation
  if(isFALSE(Save.plots)==FALSE){
    if(Save.plots==TRUE){
      path.result<-Res.DE.analysis$Path.result
    }else{
      path.result<-Save.plots
    }# if(Save.plots==TRUE)
    #
    if(is.null(path.result)==FALSE){
      if(is.null(Res.DE.analysis$Folder.result)==FALSE){
        SufixDE<-paste("_",Res.DE.analysis$Folder.result,sep="")
      }else{
        SufixDE<-NULL
      }
      SuppPlotFolder<-paste("2-4_Supplementary_Plots",SufixDE,sep="")
      if(SuppPlotFolder%in%dir(path = path.result)==FALSE){
        print("Folder creation")
        dir.create(path=paste(path.result,"/",SuppPlotFolder,sep=""))
        path.result.f<-paste(path.result,"/",SuppPlotFolder,sep="")
      }else{
        path.result.f<-paste(path.result,"/",SuppPlotFolder,sep="")
      }# if(SuppPlotFolder%in%dir(path = path.result)==FALSE)
    }else{
      path.result.f<-NULL
    }# if(is.null(path.result)==FALSE)
    #
    if(is.null(path.result.f)==FALSE){
      if("Plots_Heatmaps"%in%dir(path = path.result.f)==FALSE){
        dir.create(path=paste(path.result.f,"/","Plots_Heatmaps",sep=""))
        path.result.Heat<-paste(path.result.f,"/","Plots_Heatmaps",sep="")
      }else{
        path.result.Heat<-paste(path.result.f,"/","Plots_Heatmaps",sep="")
      }# if("Plots_Heatmaps"%in%dir(path = path.result.f)==FALSE)
    }else{
      path.result<-NULL
      path.result.f<-NULL
    }# if(is.null(path.result.f)==FALSE)
  }else{
    path.result<-NULL
  }# if(isFALSE(Save.plots)==FALSE)
  #---------------------------------------------------------------------------#
  if(is.null(path.result)==FALSE){
    grDevices::pdf(file=paste(path.result.Heat,"/","RleHeatmap",".pdf",sep=""),
                   width = 11, height = 8)#width = 8, height = 11
    print(H.SG)
    grDevices::dev.off()
    #
    grDevices::pdf(file=paste(path.result.Heat,"/","CorrelationHeatmap",".pdf",
                              sep=""),
                   width = 11, height = 8)#width = 8, height = 11
    print(H.cor)
    grDevices::dev.off()
  }else{
    print(H.SG)
    print(H.cor)
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  return(list(Zscores=RleScaled[OrderLog2FC,],
              Heatmap.Correlation=H.cor,
              Heatmap.Zscores=H.SG))
}# DEplotHeatmaps()
