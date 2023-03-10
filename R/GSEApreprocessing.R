#' @title GSEA preprocessing for official software and online tools.
#'
#' @description The function returns, from the output of [DEanalysisGlobal()],
#' specific files designed to be used as input for several online online tools
#' and software given in the section \code{Value}.
#'
#' @details
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
#' * If \code{Set.Operation="setdiff"} then the rows extracted from data are
#' those such that only one element of the selected columns in
#' \code{Res.DE.analysis$DE.results} is >0.
#' For example, the rows extracted from \code{Res.DE.analysis$DE.results}
#' will be those DE at only one time ti (except the reference time t0).
#'
#'
#' @param Res.DE.analysis A list corresponding to the output from
#' [DEanalysisGlobal()] (see \code{Examples}).
#' @param ColumnsCriteria A vector of integers where each integer indicates a
#' column of \code{Res.DE.analysis$DE.results}.
#' These columns should either contain only binary values, or may contain other
#' numerical value, in which case extracted rows from
#' \code{Res.DE.analysis$DE.results} will be those with >0 values
#' (see \code{Details}).
#' @param Set.Operation A character. The user must choose between "union"
#' (default), "intersect", "setdiff" (see \code{Details}).
#' @param Save.files \code{TRUE} or \code{FALSE} or a Character.
#' If \code{Save.files=TRUE} and the \code{path.result} of [DEanalysisGlobal()]
#' is not NULL, all files will be saved in
#' "2_SupervisedAnalysis_\code{Name.folder.DE}/2-5_Enrichment_analysis_\code{Name.folder.DE}/
#' 2-5-2_EnrichmentGO_software_preprocessing".
#' If \code{Save.files} is a character, it must be a path and all files
#' will be saved in the sub-folder "EnrichmentGO_software_preprocessing".
#' Otherwise, the different files will not be saved.
#'
#' @return The function returns
#' * A vector of character containing gene names specified by
#' \code{ColumnsCriteria} and \code{Set.Operation}.
#' * A vector of character containing all gene names
#' * And, in case where \code{Save.files=TRUE} and the \code{path.result}
#' of [DEanalysisGlobal()] is not NULL, specific files designed to be used
#' as input for the following online tools and software :
#'   * GSEA : \url{https://www.gsea-msigdb.org/gsea/index.jsp}
#'   * DAVID : \url{https://david.ncifcrf.gov/tools.jsp}
#'   * WebGestalt : \url{http://www.webgestalt.org}
#'   * gProfiler : \url{https://biit.cs.ut.ee/gprofiler/gost}
#'   * Panther : \url{http://www.pantherdb.org}
#'   * ShinyGO : \url{http://bioinformatics.sdstate.edu/go/}
#'   * Enrichr : \url{https://maayanlab.cloud/Enrichr/}
#'   * GOrilla : \url{http://cbl-gorilla.cs.technion.ac.il}.
#'
#' @export
#'
#' @importFrom stats aggregate
#' @importFrom utils write.table
#'
#' @examples
#' data(Results_DEanalysis_sub500)
#' # Results of DEanalysisGlobal() with the dataset of Antoszewski
#' res.all<-Results_DEanalysis_sub500$DE_Antoszewski2022_MOUSEsub500
#' #
#' resGp<-GSEApreprocessing(Res.DE.analysis=res.all,
#'                          ColumnsCriteria=2,
#'                          Set.Operation="union",
#'                          Save.files=FALSE)
#' #
#' #--------------------------------------------------------------------------#
#' ## The results res.all of DEanalysisGlobal with the dataset Antoszewski2022
#' # data(RawCounts_Antoszewski2022_MOUSEsub500)
#' # res.all<-DEanalysisGlobal(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
#' #                           Column.gene=1, Group.position=1,
#' #                           Time.position=NULL, Individual.position=2,
#' #                           pval.min=0.05, log.FC.min=1,LRT.supp.info=FALSE,
#' #                           path.result=NULL, Name.folder.DE=NULL)

GSEApreprocessing<-function(Res.DE.analysis,
                            ColumnsCriteria,
                            Set.Operation,
                            Save.files=FALSE){
  #---------------------------------------------------------------------------#
  # RLE normalized data
  scaled.data<-round(Res.DE.analysis$List.Datas$RLEdata[,-1], digits=3)
  #---------------------------------------------------------------------------#
  # Selection of DE gene
  ReDEsel<-DEanalysisSubData(Data=data.frame(scaled.data),
                             Res.DE.analysis=Res.DE.analysis,
                             ColumnsCriteria=ColumnsCriteria,
                             Set.Operation=Set.Operation)
  Norm.dat<-ReDEsel$SubData
  NameAllG<-Res.DE.analysis$List.Datas$RLEdata$Gene#RawCounts[,column.gene]
  Name.G<-NameAllG[ReDEsel$RowsSelected]
  row.names(Norm.dat)<-Name.G
  #---------------------------------------------------------------------------#
  if(isFALSE(Save.files)==FALSE & is.null(Res.DE.analysis$Path.result)==FALSE){
    if(Save.files==TRUE){
      path.result<-Res.DE.analysis$Path.result
    }else{
      path.result<-Save.files
    }# if(Save.files==TRUE)
    #
    if(is.null(Res.DE.analysis$Folder.result)==FALSE){
      SufixDE<-paste("_",Res.DE.analysis$Folder.result,sep="")
    }else{
      SufixDE<-NULL
    }# if(is.null(Res.DE.analysis$Folder.result)==FALSE)
    #
    if(Save.files==TRUE){
      SuppPlotFolder<-paste("2-5_Enrichment_analysis", SufixDE,sep="")
      if(SuppPlotFolder%in%dir(path=path.result)==FALSE){
        print("Folder creation")
        dir.create(path=paste(path.result,"/",SuppPlotFolder,sep=""))
        path.result.f<-paste(path.result,"/",SuppPlotFolder,sep="")
      }else{
        path.result.f<-paste(path.result,"/",SuppPlotFolder,sep="")
      }# if(SuppPlotFolder%in%dir(path = path.result)==FALSE)
      #
      if("2-5-2_EnrichmentGO_software_preprocessing"%in%dir(path=path.result.f)==FALSE){
        # print("Folder creation")
        dir.create(path=paste(path.result.f,"/",
                              "2-5-2_EnrichmentGO_software_preprocessing",
                              sep=""))
      }
      path.resultf<-paste(path.result.f,"/",
                          "2-5-2_EnrichmentGO_software_preprocessing",sep="")
    }else{
      SuppPlotFolder<-paste("EnrichmentGO_software_preprocessing",SufixDE,
                            sep="")
      if(SuppPlotFolder%in%dir(path = path.result)==FALSE){# GOFolder%in%
        # print("Folder creation")
        dir.create(path=paste(path.result,"/",SuppPlotFolder,sep=""))
      }# if(SuppPlotFolder%in%dir(path=path.result)==FALSE)
      path.resultf<-paste(path.result,"/",SuppPlotFolder,sep="")
    }# if(Save.files==TRUE)
    #-------------------------------------------------------------------------#
    if("GSEA"%in%dir(path = path.resultf)==FALSE){
      dir.create(path=paste(path.resultf,"/","GSEA",sep=""))
    }
    #
    if("GSEApreranked"%in%dir(path=paste(path.resultf,"/","GSEA",sep=""))==FALSE){
      dir.create(path=paste(path.resultf,"/","/GSEA/GSEApreranked",sep=""))
    }
    #
    if("Webgestalt"%in%dir(path = path.resultf)==FALSE){
      dir.create(path=paste(path.resultf,"/","Webgestalt",sep=""))
    }
    path.resultw<-paste(path.resultf,"/","Webgestalt",sep="")
    #
    if("ORA_NTA_analysis"%in%dir(path = path.resultw)==FALSE){
      dir.create(path=paste(path.resultw,"/","ORA_NTA_analysis",sep=""))
    }
    #
    if("GSEA_analysis"%in%dir(path = path.resultw)==FALSE){
      dir.create(path=paste(path.resultw,"/","GSEA_analysis",sep=""))
    }
    #
    OtherEnrichGO<-"DAVID_gProfiler_Panther_ShinyGO_Enrichr_GOrilla"
    if("DAVID_gProfiler_Panther_ShinyGO_Enrichr_GOrilla"%in%dir(path = path.resultf)==FALSE){
      dir.create(path=paste(path.resultf,"/",OtherEnrichGO,sep=""))
    }
    #-------------------------------------------------------------------------#
    if(ncol(Res.DE.analysis$Summary.Inputs$FactorsInfo)==3){
      TnumFct<-gsub("t","",gsub("T","",
                                Res.DE.analysis$Summary.Inputs$FactorsInfo[,2],
                                fixed=TRUE))
      FactorBoxplt<-as.factor(paste(Res.DE.analysis$Summary.Inputs$FactorsInfo[,1],
                                    paste("t",TnumFct,sep=""),
                                    sep="."))
    }else{
      FactorBoxplt<-as.factor(Res.DE.analysis$Summary.Inputs$FactorsInfo[,1])
    }# if(ncol(Res.DE.analysis$Summary.Inputs$FactorsInfo)==3)
    #-------------------------------------------------------------------------#
    IndexOrder<-order(as.character(FactorBoxplt))
    Norm.dat<-Norm.dat[,IndexOrder]
    #
    # GCT file for GSEA
    options(warn=-1)
    filename.Norm <- paste(path.resultf,"/GSEA/ScaledData", ".gct", sep="")
    cat("#1.2", "\n", sep="\t", file=filename.Norm)
    cat(nrow(Norm.dat), ncol(Norm.dat), "\n", sep="\t",
        file=filename.Norm, append=TRUE)
    utils::write.table(cbind(Name.G, Name.G,Norm.dat),
                       file=filename.Norm, row.names=FALSE,
                       quote=FALSE, sep="\t", append=TRUE)
    options(warn=0)# suppressWarnings(utils::write.table())
    #-------------------------------------------------------------------------#
    # CLS file for GSEA
    LevelsFct<-levels(FactorBoxplt)
    VectFct<-as.vector(FactorBoxplt)[IndexOrder]
    level<-LevelsFct[match(LevelsFct, unique(VectFct))]
    filenameFct<-paste(path.resultf,"/GSEA/Factors", ".cls", sep="")
    cat(c(length(VectFct), length(level), 1), "\n", sep=" ", file=filenameFct)
    cat(c("#", level), "\n", sep=" ", file=filenameFct, append=TRUE)
    cat(VectFct, sep=" ", file=filenameFct, append=TRUE)
    #-------------------------------------------------------------------------#
    # txt files for OtherEnrichGO
    filepathOtherEnrichGO<-paste(path.resultf,"/",OtherEnrichGO,sep="")
    utils::write.table(data.frame(Name.G),
                       file=paste(filepathOtherEnrichGO,"/TargetDEGene",".txt",
                                  sep=""),
                       row.names=FALSE, col.names=FALSE,quote = FALSE)
    utils::write.table(data.frame(NameAllG),
                       file=paste(filepathOtherEnrichGO,"/BackgroungGene",
                                  ".txt", sep=""),
                       row.names=FALSE, col.names=FALSE,quote = FALSE)
    #-------------------------------------------------------------------------#
    # txt file for webgestalt (ORA)
    filenameORA<-paste(path.resultw,
                       "/ORA_NTA_analysis/WebgestaltORANTATargetGene",".txt",
                       sep="")
    utils::write.table(data.frame(Name.G), file=filenameORA,
                       row.names=FALSE, col.names=FALSE,quote = FALSE)
    utils::write.table(data.frame(NameAllG),
                       file=paste(path.resultw,
                                  "/ORA_NTA_analysis/BackgroungGene",".txt",
                                  sep=""),
                       row.names=FALSE, col.names=FALSE,quote = FALSE)
    #-------------------------------------------------------------------------#
    # rnk file for webgestalt (GSEA)
    # http://www.webgestalt.org
    MNallg<-stats::aggregate(t(Norm.dat), by=list(VectFct), FUN=mean)
    SDallg<-stats::aggregate(t(Norm.dat), by=list(VectFct), FUN=sd)
    NperFct<-as.numeric(table(VectFct))
    for(fc in seq_len(length(LevelsFct)-1)){
      for(fnc in seq(from=(fc+1), to=length(LevelsFct), by=1)){
        Dmn<-as.numeric(MNallg[fnc,-1]-MNallg[fc,-1])
        Dsd<-as.numeric(sqrt((SDallg[fnc,-1]^2)/NperFct[fnc]+(SDallg[fc,-1]^2)/NperFct[fc]))
        filenameWGSEA<-paste(path.resultw,
                             "/GSEA_analysis/WebgestaltGSEATargetGene_",
                             paste("Tscore",LevelsFct[fnc],"versus",
                                   LevelsFct[fc],sep="_"),
                             ".rnk", sep="")
        utils::write.table(data.frame(Gene=Name.G,
                                      Score=round(Dmn/Dsd,digits=4)),
                           sep="\t",
                           file=filenameWGSEA, row.names=FALSE,
                           col.names=FALSE, quote = FALSE)
        #
        Prerankedpath<-paste(path.resultf,
                             "/GSEA/GSEApreranked/WebgestaltGSEATargetGene_",
                             paste("Tscore",LevelsFct[fnc],"versus",
                                   LevelsFct[fc],sep="_"),
                             ".rnk", sep="")
        #
        utils::write.table(data.frame(Gene=Name.G,
                                      Score=round(Dmn/Dsd,digits=4)),
                           sep="\t",
                           file=Prerankedpath, row.names=FALSE,
                           col.names=FALSE, quote = FALSE)
      }# end for(fnc in (fc+1):length(LevelsFct))
    }# end for(fc in 1:(length(LevelsFct)-1))
    #-------------------------------------------------------------------------#
    #-------------------------------------------------------------------------#
    # Creation file Readme
    fileRdme<-paste(path.resultf,"/","ReferenceEnrichmentGO.txt",sep="")
    ReferenceEnrichmentGOtxt(path.fileRdme=fileRdme)
  }# if(is.FALSE(Save.files)==FALSE)
  #---------------------------------------------------------------------------#
  return(list(TargetDEGene=Name.G,
              BackgroungGene=NameAllG))
}# GSEApreprocessing()


ReferenceEnrichmentGOtxt<-function(path.fileRdme){
  #
  fileRdme<-path.fileRdme
  #---------------------------------------------------------------------------#
  cat(rep("=",times=75), "\n", sep="", file=fileRdme)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c(rep("=",times=5),"   GSEA: http://www.gsea-msigdb.org/gsea/index.jsp"),
      "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(c("* Link of the article associated with GSEA:"),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  cat(c("  - Title: Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles."),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - Number of citations : 31586 (from 25/10/2005 to 22/04/2022)",
        "mean rate = 1914 citations per year"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - https://pubmed.ncbi.nlm.nih.gov/16199517/"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 1 from the youtube channel 'IUSM-CCBB'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","B4B: Module 4 - Functional Analysis - GSEA (theory)"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=ibyyUpgG4Wo", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsGSEA1<-"No real youtube description."
  cat(c("  - Youtube description:",AbsGSEA1),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 2 from the youtube channel 'IUSM-CCBB'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","B4B: Module 4 - Functional Analysis - GSEA (hands on)"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=JIMZipAPKkk", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsGSEA2<-"No real youtube description."
  cat(c("  - Youtube description:",AbsGSEA2),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 3 from the youtube channel 'Genomics Gurus'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","HOW TO PERFORM GSEA - A tutorial on gene set enrichment analysis for RNA-seq"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=KY6SS4vRchY", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsGSEA3<-"In this tutorial, we explain what gene set enrichment analysis (GSEA) is and what it offers you. We show you how to run the analysis on your computer and take you through how to interpret the outputs. The tutorial also covers leading edge analysis and analysis of gene networks with Cytoscape."
  cat(c("  - Youtube description:",AbsGSEA3),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 4 from the youtube channel",
        "'Online Faculty Mentoring Network to Develop Video Tutorials for Computational Genomics'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","GSEAtheory"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=bT00oJh2x_4", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsGSEA4<-"Brief introduction to how gene set enrichment analysis works."
  cat(c("  - Youtube description:",AbsGSEA4),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c(rep("=",times=5),"   DAVID: https://david.ncifcrf.gov/"), "\n",
      sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(c("* Link of the article associated with DAVID:"),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  cat(c("  - Title: Systematic and integrative analysis of large gene lists using DAVID bioinformatics resources."),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - Number of citations : 23137 (from 18/12/2008 to 22/04/2022)",
        "mean rate = 1629 citations per year"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - https://www.nature.com/articles/nprot.2008.211"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 1 from the youtube channel 'Sanbomics'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","Simple gene ontology and pathway enrichment from a gene list"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=XLRA0A5qsoE", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsDavid1<-"Here I show you how to use the DAVID functional annotation resource to easily do gene ontology and KEGG pathway enrichment from a gene list. This method takes genes in almost any format and is incredibly easy to use. This requires no programatic skills and barely any thinking"
  cat(c("  - Youtube description:",AbsDavid1),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 2 from the youtube channel 'Bioinformatics'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","How to do Gene Ontology (GO) enrichment analysis with DAVID"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=9zFqtviLUzQ", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsDavid2<-"Basic usage of DAVID: finding  GO Direct pathways of Biological Process (BP), Cellular Component (CC), and Molecular Function (MF). How to select a set of single list of genes on spread sheet, and paste the list on DAVID."
  cat(c("  - Youtube description:",AbsDavid2),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 3 from the youtube channel 'Bioinformatics'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("TITLE:","using DAVID to convert gene IDs to official gene symbols."),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=NwPMVlqzNUg", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsDavid3<-"No real youtube description"
  cat(c("  - Youtube description:",AbsDavid3),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c(rep("=",times=5),"   WebGestalt: http://www.webgestalt.org"), "\n",
      sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(c("* Link of the article associated with WebGestalt:"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - Title: WebGestalt 2019: gene set analysis toolkit with revamped UIs and APIs."),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - Number of citations : 1079 (from 02/07/2019 to 22/04/2022)",
        "mean rate = 600 citations per year"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - https://academic.oup.com/nar/article/47/W1/W199/5494758"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 1 from the youtube channel 'IIT Bombay July 2018'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","Lecture 49 : WebGestalt - I"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=zDNWc9YiKpo", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsWebGSEA1<-"No real youtube description."
  cat(c("  - Youtube description:",AbsWebGSEA1),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 2 from the youtube channel 'IIT Bombay July 2018'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","Lecture 50 : WebGestalt - II"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=OCmO_H6MoLo", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsWebGSEA2<-"No real youtube description."
  cat(c("  - Youtube description:",AbsWebGSEA2),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c(rep("=",times=5),"   Panther: http://pantherdb.org/"), "\n",
      sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(c("* Link of the article associated with Panther:"),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  cat(c("  - Title: PANTHER version 16: a revised family classification, tree-based classification tool, enhancer regions and extensive API."),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - Number of citations : 323 (from 08/01/2021 to 22/04/2022)",
        "mean rate = 323 citations per year"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - https://academic.oup.com/nar/article/49/D1/D394/6027812?login=false"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 1 from the youtube channel 'Bioinformatics World'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","PANTHER database | PANTHER tutorial | Bioinformatics tutorial |"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=Xqp-sNb8PBY", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  Abspanther1<-"No real youtube description."
  cat(c("  - Youtube description:",Abspanther1),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c(rep("=",times=5),"   g:profiler: https://biit.cs.ut.ee/gprofiler/gost"),
      "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(c("* Link of the article associated with g:profiler:"),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  cat(c("  - Title: g:Profiler: a web server for functional enrichment analysis and conversions of gene lists"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - Number of citations : 834 (from 02/07/2019 to 22/04/2022)",
        "mean rate = 298 citations per year"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - https://pubmed.ncbi.nlm.nih.gov/31066453/"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 1 from the youtube channel 'Saniya Khullar'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","Overview of gProfiler Tools for Bioinformatics/Biological Analysis"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=QH8dD0tcNpg", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  Absgprofile1<-"Please note that gProfiler is a tool for functional enrichment analysis, gene identifier conversion, mapping homologous genes across related organisms, and mapping SNPs to functional effects and target genes. There are 4 key tools available on gProfiler, such as:  GOSt, Convert, Orth, and SNPense.  Please note that Saniya has dedicated videos on each of these 4 tools."
  cat(c("  - Youtube description:",Absgprofile1),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 2 from the youtube channel 'Saniya Khullar'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","g:Profiler for Gene-Set Enrichment Analysis"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=Ymu7EruyIeE&t=0s", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  Absgprofile2<-"Please note that in this video, Saniya goes through how you can use g:Profiler (https://biit.cs.ut.ee/gprofiler/gost) to perform gene-set enrichment analysis (GSEA), along with other functions, and also download gene data sources as .gmt files.  That is, for a given list of genes (gene names and/or gene Entrez IDs or Ensembl IDs), you can use g:Profiler to uncover biological pathways, networks, functions, associations, tissues, etc. that are enriched (have a statistically significant overlap with our gene list) for our list of genes.  In short, through Metascape (1 of many tools out there), we can query hundreds of different biological gene lists (across a variety of databases and atlases like Gene Ontology, biological pathways (e.g. KEGG, Reactome, WikiPathways), regulatory motifs in DNA (e.g. TRANSFAC, miRTarBase), protein databases (like Human Protein Atlas or CORUM), human phenotype ontology (e.g. HP) to determine what significant biological results are there for our gene list.  Saniya also explains how results differ for different gene sets."
  cat(c("  - Youtube description:",Absgprofile2),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c(rep("=",times=5),"   Enrichr: https://maayanlab.cloud/Enrichr/"), "\n",
      sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(c("* Link of the article associated with Enrichr:"),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  cat(c("  - Title: Enrichr: a comprehensive gene set enrichment analysis web server"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - Number of citations : 4598 (from 03/05/2016 to 22/04/2022)",
        "mean rate = 766 citations per year"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4987924/"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 1 from the youtube channel '
de.STAIR'"),
"\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","Enrichment analysis with Enrichr"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=qTfOXAObNwo", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsEnrich1<-" In this short video we show the use of Enrichr (Chen et al. 2013, Kuleshov et al. 2016) to carry out the enrichment analysis of a list of differtially expressed genes, which we obtained from the RNA-Seq analysis of two NGS Human breast cancer data sets."
  cat(c("  - Youtube description:",AbsEnrich1),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c(rep("=",times=5),"   ShinyGO: http://bioinformatics.sdstate.edu/go/"),
      "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(c("* Link of the article associated with ShinyGO:"),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  cat(c("  - Title: ShinyGO: a graphical gene-set enrichment tool for animals and plants"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - Number of citations : 424 (from 15/04/2009 to 22/04/2022)",
        "mean rate = 212 citations per year"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - https://academic.oup.com/bioinformatics/article/36/8/2628/5688742"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - https://www.biorxiv.org/content/biorxiv/suppl/2018/05/04/315150.DC1/315150-1.pdf"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("-",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c("* VIDEO 1 from the youtube channel 'Dr. Asif's Mol. Biology'"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - TITLE:","Gene ontology : GO enrichment analysis | shiny GO | Web tool"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat("  - https://www.youtube.com/watch?v=Xqp-sNb8PBY", "\n", sep=" ",
      file=fileRdme, append=TRUE)
  AbsShinyGO1<-"In this video, I have explained how can we use an online tool for generating gene ontology enrichment graphs? Go shiny is a web tool and can make the GO enrichment graphs, PPI and KEGG pathways in seconds."
  cat(c("  - Youtube description:",AbsShinyGO1),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(c(rep("=",times=5),"   GOrilla: http://cbl-gorilla.cs.technion.ac.il"),
      "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(c("* Link of the article associated with GOrilla:"),"\n", sep=" ",
      file=fileRdme, append=TRUE)
  cat(c("  - Title: GOrilla: A Tool For Discovery And Visualization of Enriched GO Terms in Ranked Gene Lists"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - Number of citations : 3044 (from 03/02/2009 to 22/04/2022)",
        "mean rate = 154 citations per year"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  cat(c("  - https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-48"),
      "\n", sep=" ", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  cat(rep("=",times=75), "\n", sep="", file=fileRdme, append=TRUE)
  #---------------------------------------------------------------------------#
}# ReferenceEnrichmentGOtxt
