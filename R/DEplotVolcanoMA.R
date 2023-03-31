#' @title Volcano and MA graphs
#'
#' @description The function returns Volcano plots and MA plots from
#' the results of our function [DEanalysisGlobal()].
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
#' @param Res.DE.analysis A list. Output from [DEanalysisGlobal()]
#' (see \code{Examples}).
#' @param NbGene.plotted Non negative integer. The algorithm computes the sum of
#' all the absolute \eqn{log_2} fold change present in the element
#' \code{DE.results} of \code{Res.DE.analysis} for each gene.
#' Only the highest \code{NbGene.plotted} genes are plotted in the volcano and
#' MA plots. By default, \code{NbGene.plotted}=2.
#' @param SizeLabel Numeric. Give the size of the names of plotted genes.
#' By default, \code{SizeLabel}=3.
#' @param Display.plots \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
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
#' @importFrom ggplot2 ggplot geom_point aes scale_shape_manual
#' scale_color_manual scale_color_manual geom_hline ggtitle xlab theme_minimal
#' @importFrom ggrepel geom_label_repel
#'
#' @return The function returns Volcano plots and MA plots
#' from the results of our function [DEanalysisGlobal()].
#'
#' @seealso The function calls the output of [DEanalysisGlobal()].
#'
#' @export
#'
#' @examples
#' data(Results_DEanalysis_sub500)
#' # Results of DEanalysisGlobal() with the dataset of Antoszewski
#' res.all<-Results_DEanalysis_sub500$DE_Antoszewski2022_MOUSEsub500
#' #
#' resVolcanoMA<-DEplotVolcanoMA(res.all,
#'                               NbGene.plotted=5,
#'                               Display.plots=TRUE,
#'                               Save.plots=FALSE)
#' #
#' #--------------------------------------------------------------------------#
#' ## The results res.all of DEanalysisGlobal with the dataset Antoszewski2022
#' # data(RawCounts_Antoszewski2022_MOUSEsub500)
#' # res.all<-DEanalysisGlobal(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
#' #                           Column.gene=1, Group.position=1,
#' #                           Time.position=NULL, Individual.position=2,
#' #                           pval.min=0.05, log.FC.min=1,LRT.supp.info=FALSE,
#' #                           path.result=NULL, Name.folder.DE=NULL)

DEplotVolcanoMA<-function(Res.DE.analysis,
                          NbGene.plotted=2,
                          SizeLabel=3,
                          Display.plots=TRUE,
                          Save.plots=FALSE){
  #---------------------------------------------------------------------------#
  # To avoid "no visible binding for global variable" with devtools::check()
  log2FoldChange<-padj<-Gene<-ColpFC<-Shape.pFC<-NULL
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
      }# if(is.null(Res.DE.analysis$Folder.result)==FALSE)
      SuppPlotFolder<-paste("2-4_Supplementary_Plots", SufixDE,sep="")
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
      if("Plots_Volcano"%in%dir(path = path.result.f)==FALSE){
        dir.create(path=paste(path.result.f,"/","Plots_Volcano",sep=""))
        path.result.vol<-paste(path.result.f,"/","Plots_Volcano",sep="")
      }else{
        path.result.vol<-paste(path.result.f,"/","Plots_Volcano",sep="")
      }# if("Plots_Volcano"%in%dir(path = path.result.f)==FALSE)
      #
      if("Plots_RatioIntensity(MA)"%in%dir(path = path.result.f)==FALSE){
        dir.create(path=paste(path.result.f,"/","Plots_RatioIntensity(MA)",
                              sep=""))
        path.result.ma<-paste(path.result.f,"/","Plots_RatioIntensity(MA)",
                              sep="")
      }else{
        path.result.ma<-paste(path.result.f,"/","Plots_RatioIntensity(MA)",
                              sep="")
      }# if("Plots_RatioIntensity(MA)"%in%dir(path = path.result.f)==FALSE)
    }else{
      path.result.vol<-NULL
      path.result.ma<-NULL
      path.result<-NULL
    }# if(is.null(path.result)==FALSE)
  }else{
    path.result<-NULL
  }# if(isFALSE(Save.plots)==FALSE)
  #---------------------------------------------------------------------------#
  # Distinction case when samples belong to different group and/or time
  if(length(Res.DE.analysis$Summary.Inputs$ExprCond)==2){
    NbGroup<-length(Res.DE.analysis$Summary.Inputs$GroupLevels)
    NbTime<-length(Res.DE.analysis$Summary.Inputs$TimeLevels)
    pevo<-c(rep(Res.DE.analysis$Summary.Inputs$pvalsTime, times=NbGroup),
            rep(Res.DE.analysis$Summary.Inputs$pvalGroup,
                times=NbTime*NbGroup*(NbGroup-1)/2))
    #
    Title.g<-rep(NA,times=NbGroup*(NbGroup-1)/2)
    cptBC<-0
    for(i in seq_len(NbGroup-1)){# seq_len(32)
      for(k in seq(from=(i+1), to=NbGroup, by=1)){# (i+1):NbGroup
        cptBC<-cptBC+1
        Title.g[cptBC]<-paste(Res.DE.analysis$Summary.Inputs$GroupLevels[k],
                              "vs",
                              Res.DE.analysis$Summary.Inputs$GroupLevels[i])
      }# for(k in (i+1):NbGroup)
    }# for(i in 1:(NbGroup-1))
    #
    Title.g2<-paste("(", rep(paste("t",0:(NbTime-1),sep=""),
                             each=NbGroup*(NbGroup-1)/2),
                    ") ", Title.g,sep="")
    #
    Title.vMA<-c(paste("(",rep(Res.DE.analysis$Summary.Inputs$GroupLevels,
                               each=NbTime-1),
                       ") t", rep(seq_len(NbTime-1), times=NbGroup),
                       " vs t0", sep=""),
                 Title.g2)
    #
  }else{
    if("Time"%in%Res.DE.analysis$Summary.Inputs$ExprCond){
      NbTime<-length(Res.DE.analysis$Summary.Inputs$TimeLevels)
      pevo<-Res.DE.analysis$Summary.Inputs$pvalsTime
      #
      Title.vMA<-paste("t",rep(seq_len(NbTime-1), times=1)," vs t0",sep="")
    }else{
      NbGroup<-length(Res.DE.analysis$Summary.Inputs$GroupLevels)
      pevo<-rep(Res.DE.analysis$Summary.Inputs$pvalGroup,
                times=NbGroup*(NbGroup-1))
      #
      Title.vMA<-rep(NA, times=NbGroup*(NbGroup-1)/2)
      cptBC<-0
      for(i in seq_len(NbGroup-1)){# 1:(NbGroup-1)
        for(k in seq(from=(i+1), to=NbGroup, by=1)){
          cptBC<-cptBC+1
          Title.vMA[cptBC]<-paste(Res.DE.analysis$Summary.Inputs$GroupLevels[k],
                                  "vs",
                                  Res.DE.analysis$Summary.Inputs$GroupLevels[i])
        }# for(k in (i+1):NbGroup)
      }# for(i in 1:(NbGroup-1))
    }# if("Time"%in%Res.DE.analysis$Summary.Inputs$ExprCond)
  }# if(length(Res.DE.analysis$Summary.Inputs$ExprCond)==2)
  #
  AbsLog2Fcmin<-Res.DE.analysis$Summary.Inputs$logFCmin
  #---------------------------------------------------------------------------#
  log2Mean<-log2(apply(Res.DE.analysis$List.Datas$RLEdata[,-1], 1 , mean)+1)
  IdLog2FC<-grep(pattern="Log2FoldChange",
                 x=colnames(Res.DE.analysis$DE.results),
                 fixed = TRUE)
  #---------------------------------------------------------------------------#
  Nlog2FC<-length(IdLog2FC)
  #
  List.plot.Volcano<-vector(mode="list", length=Nlog2FC)
  List.plot.MA<-vector(mode="list", length=Nlog2FC)
  names(List.plot.MA)<-names(List.plot.Volcano)<-gsub(" ","",Title.vMA)
  #---------------------------------------------------------------------------#
  # Volcano and MA plots for each case
  for(i in seq_len(Nlog2FC)){# 1:length(IdLog2FC)
    DatVolMA<-data.frame(Gene=Res.DE.analysis$List.Datas$RLEdata$Gene,
                         log2Mean=log2Mean,
                         log2FoldChange=NA,
                         padj=NA,
                         DEpFC=rep("Not significant", times=length(log2Mean)),
                         ColpFC="#999999",
                         Shape.pFC=18)
    row.names(DatVolMA)<-NULL
    #-------------------------------------------------------------------------#
    DatVolMA$log2FoldChange<-Res.DE.analysis$DE.results[,IdLog2FC[i]]
    DatVolMA$padj<-Res.DE.analysis$DE.results[,IdLog2FC[i]+1]
    #-------------------------------------------------------------------------#
    if(min(DatVolMA$padj)==0 & length(unique(DatVolMA$padj))>1){
      DatVolMA$padj[which(DatVolMA$padj==0)]<-unique(sort(DatVolMA$padj))[2]*0.01
    }# if(min(DatVolMA$padj)==0 & length(unique(DatVolMA$padj))>1)
    #
    DatVolMA$DEpFC[which(DatVolMA$padj<pevo[i])]<-"Pvalue-adj"
    SelpvFC<-which(DatVolMA$padj<pevo[i] & abs(DatVolMA$log2FoldChange)>AbsLog2Fcmin)
    DatVolMA$DEpFC[SelpvFC]<-"Pvalue-adj & Log2FC"
    #-------------------------------------------------------------------------#
    Gene.ns<-which(DatVolMA$DEpFC=="Not significant")
    Gene.p<-which(DatVolMA$DEpFC=="Pvalue-adj")
    Gene.pFC<-which(DatVolMA$DEpFC=="Pvalue-adj & Log2FC")
    #
    if(length(Gene.p)>0){
      DatVolMA$ColpFC[Gene.p]<-"#56B4E9"
      DatVolMA$Shape.pFC[Gene.p]<-17
    }# if(length(Gene.p)>0)
    #
    if(length(Gene.pFC)>0){
      DatVolMA$ColpFC[Gene.pFC]<-"#E69F00"
      DatVolMA$Shape.pFC[Gene.pFC]<-16
    }# if(length(Gene.pFC)>0)
    #-------------------------------------------------------------------------#
    DatVolMA$Shape.pFC<-as.factor(DatVolMA$Shape.pFC)
    DatVolMA$ColpFC<-as.factor(DatVolMA$ColpFC)
    #-------------------------------------------------------------------------#
    Eucl.pFC.Volc<-sqrt(DatVolMA$log2FoldChange^2 + log10(DatVolMA$padj)^2)
    #
    IdPlotted.Volc<-c(intersect(order(Eucl.pFC.Volc,decreasing=TRUE),Gene.pFC),
                      intersect(order(Eucl.pFC.Volc, decreasing=TRUE),Gene.p),
                      intersect(order(Eucl.pFC.Volc, decreasing=TRUE),
                                Gene.ns))[seq_len(NbGene.plotted)]
    #-------------------------------------------------------------------------#
    options(warn = - 1)
    #
    gVolca<-ggplot2::ggplot(data=DatVolMA,
                            ggplot2::aes(x=log2FoldChange, y=-log10(padj),
                                         label=Gene)) +
      ggrepel::geom_label_repel(data=DatVolMA[IdPlotted.Volc,],
                                col=DatVolMA$ColpFC[IdPlotted.Volc],
                                box.padding=0.5, max.overlaps=NbGene.plotted+2,
                                min.segment.length=0, size=SizeLabel)+
      ggplot2::geom_point(ggplot2::aes(color=ColpFC,shape=Shape.pFC))+
      ggplot2::scale_shape_manual(values=c(18, 17, 16),
                                  name="Criteria",
                                  breaks=c(18, 17, 16),
                                  labels=c("Not significant", "Pvalue-adj",
                                           "Pvalue-adj & Log2FC"))+
      ggplot2::scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00"),
                                  name="Criteria",
                                  breaks=c("#999999", "#56B4E9", "#E69F00"),
                                  labels=c("Not significant", "Pvalue-adj",
                                           "Pvalue-adj & Log2FC"))+
      ggplot2::geom_vline(xintercept=c(-AbsLog2Fcmin, AbsLog2Fcmin),
                          linetype="dashed",col="black") +
      ggplot2::geom_hline(yintercept=-log10(pevo[i]),
                          linetype="dashed",col="black")+
      ggplot2::ggtitle(paste("Volcano plot : ",Title.vMA[i],sep=""))+
      ggplot2::theme_minimal()
    #
    List.plot.Volcano[[i]]<-gVolca
    #
    options(warn = 0)
    #-------------------------------------------------------------------------#
    options(warn = - 1)
    #
    gMA<-ggplot2::ggplot(data=DatVolMA,
                         ggplot2::aes(x=log2Mean,y=log2FoldChange,label=Gene))+
      ggrepel::geom_label_repel(data=DatVolMA[IdPlotted.Volc,],
                                col=DatVolMA$ColpFC[IdPlotted.Volc],
                                box.padding=0.5, max.overlaps=NbGene.plotted+2,
                                min.segment.length=0, size = SizeLabel)+
      ggplot2::geom_point(ggplot2::aes(color=ColpFC,shape=Shape.pFC))+
      ggplot2::scale_shape_manual(values=c(18, 17, 16),
                                  name="Criteria",
                                  breaks=c(18, 17, 16),
                                  labels=c("Not significant", "Pvalue-adj",
                                           "Pvalue-adj & Log2FC"))+
      ggplot2::scale_color_manual(values=c("#999999", "#56B4E9", "#E69F00"),
                                  name="Criteria",
                                  breaks=c("#999999", "#56B4E9", "#E69F00"),
                                  labels=c("Not significant", "Pvalue-adj",
                                           "Pvalue-adj & Log2FC"))+
      ggplot2::geom_hline(yintercept=-AbsLog2Fcmin,
                          linetype="dashed", col="black")+
      ggplot2::geom_hline(yintercept=0, col="black")+
      ggplot2::geom_hline(yintercept=AbsLog2Fcmin,
                          linetype="dashed",col="black")+
      ggplot2::ggtitle(paste("MA plot : ",Title.vMA[i],sep=""))+
      ggplot2::xlab("log2(Mean rle +1)")+
      ggplot2::theme_minimal()
    #
    List.plot.MA[[i]]<-gMA
    #
    options(warn = 0)
    #-------------------------------------------------------------------------#
    options(warn = - 1)
    #
    if(is.null(path.result)==FALSE){
      if(Save.plots==TRUE & is.null(Res.DE.analysis$Folder.result)==FALSE){
        SuffixSaveplots<-paste("_",Res.DE.analysis$Folder.result,sep="")
      }else{
        SuffixSaveplots<-""
      }# if(Save.plots==TRUE & is.null(Res.DE.analysis$Folder.result)==FALSE)
      grDevices::pdf(file=paste(path.result.vol,"/","Volcano_",
                                gsub(" ", "_", x=Title.vMA[i],fixed = TRUE),
                                SuffixSaveplots, ".pdf",sep=""),
                     width=11, height=8)
      print(gVolca)
      grDevices::dev.off()
      #
      grDevices::pdf(file=paste(path.result.ma,"/",
                                "RatioIntensityPlot(MAplot)_",
                                gsub(" ", "_", x=Title.vMA[i],fixed = TRUE),
                                SuffixSaveplots, ".pdf",sep=""),
                     width=11, height=8)
      print(gMA)
      grDevices::dev.off()
    }else{
      if(Display.plots==TRUE){
        print(gVolca)
        print(gMA)
      }# if(Display.plots==TRUE)
    }# if(is.null(path.result)==FALSE)
    #
    options(warn = 0)
    #-------------------------------------------------------------------------#
  }# for(i in 1:length(IdLog2FC))
  #---------------------------------------------------------------------------#
  return(list(List.plot.Volcano=List.plot.Volcano,
              List.plot.MA=List.plot.MA))
}# DEplotVolcanoMA()
