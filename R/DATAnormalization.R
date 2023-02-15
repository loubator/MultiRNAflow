#' @title Normalization of raw counts.
#'
#' @description From raw counts, this function realizes one of
#' the three methods of normalization of the package \code{DESeq2}:
#' * Relative Log Expression (rle) transformation
#' (see [BiocGenerics::estimateSizeFactors()])
#' * Regularized Log (rlog) transformation (see [DESeq2::rlog()])
#' * Variance Stabilizing Transformation (vst) transformation
#' (see [DESeq2::vst()])
#'
#' @details The column names of \code{ExprData} must be a vector of strings
#' of characters containing
#' * a string of characters (if \eqn{k=1}) which is the label of
#' the column containing gene names.
#' * \eqn{N_s} sample names which must be strings of characters containing
#' at least: the name of the individual (e.g patient, mouse, yeasts culture),
#' its biological condition (if there is at least two) and
#' the time where data have been collected if there is at least two;
#' (must be either 't0', 'T0' or '0' for time 0,
#' 't1', 'T1' or '1' for time 1, ...).
#'
#' All these sample information must be separated by underscores in
#' the sample name. For instance 'CLL_P_t0_r1',
#' corresponds to the patient 'r1' belonging to the biological condition 'P'
#' and where data were collected at time 't0'. I this example, 'CLL' describe
#' the type of cells (here chronic lymphocytic leukemia) and
#' is not used in our analysis.
#'
#' In the string of characters 'CLL_P_t0_r1',
#' 'r1' is localized after the third underscore, so \code{Individual.position=4},
#' 'P' is localized after the first underscore, so \code{Group.position=2} and
#' 't0' is localized after the second underscore, so \code{Time.position=3}.
#'
#' @param RawCounts Data.frame with \eqn{N_g} rows and (\eqn{N_{s+k}}) columns,
#' where \eqn{N_g} is the number of genes,
#' \eqn{N_s} is the number of samples and
#' \eqn{k=1} if a column is used to specify gene names,
#' or \eqn{k=0} otherwise.
#' If \eqn{k=1}, the position of the column containing gene names is given
#' by \code{Column.gene}.
#' The data.frame contains non negative integers giving gene expressions of
#' each gene in each sample.
#' Column names of the data.frame must describe each sample's information
#' (individual, biological condition and time) and have the structure described
#' in the section \code{Details}.
#' @param Column.gene Integer indicating the column where gene names are given.
#' Set \code{Column.gene=NULL} if there is no such column.
#' @param Group.position Integer indicating the position of group information
#' in the string of characters in each sample names (see \code{Details}).
#' Set \code{Group.position=NULL} if there is only one or
#' no biological information in the string of character in each sample name.
#' @param Time.position Integer indicating the position of time measurement
#' information in the string of characters in each sample names
#' (see \code{Details}).
#' Set \code{Time.position=NULL} if there is only one or no time measurement
#' information in the string of character in each sample name.
#' @param Individual.position Integer indicating the position of the name of
#' the individual (e.g patient, replicate, mouse, yeasts culture ...)
#' in the string of characters in each sample names (see \code{Details}).
#' The names of different individuals must be all different.
#' Furthermore, if individual names are just numbers, they will be transform in
#' a vector of class "character" by [CharacterNumbers()] and
#' a "r" will be added to each individual name ("r" for replicate).
#' @param Normalization "rle", "vst", "rlog".
#' Each corresponds to a method of normalization proposed by \code{DESeq2}
#' (see [BiocGenerics::estimateSizeFactors()] for "rle",
#' [DESeq2::rlog()] for "rlog" and [DESeq2::vst()] for "vst").
#' @param Blind.rlog.vst TRUE or FALSE. See input 'blind' in [DESeq2::rlog()].
#' It is recommended to set \code{Blind.rlog.vst=FALSE} for downstream analysis.
#' @param Plot.Boxplot \code{TRUE} or \code{FALSE}. TRUE by default.
#' If \code{Plot.Boxplot=TRUE}, the function [DATAplotBoxplotSamples()] will be
#' called and boxplots will be plotted. Otherwise, no boxplots will be plotted.
#' @param Colored.By.Factors \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, boxplots will be colored with different colors for different
#' time measurements (if data were collected at different time points).
#' Otherwise, boxplots will be colored with different colors for different
#' biological conditions.
#' @param Color.Group \code{NULL} or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute a color
#' for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param Plot.genes \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, points representing gene expressions
#' (normalized or raw counts) will be plotted for each sample.
#' Otherwise, only boxplots will be plotted.
#' @param path.result Character or \code{NULL}. Path to save all results.
#' If \code{path.result} contains a sub folder entitled
#' "1_Normalization_\code{Name.folder.norm}"
#' all results will be saved in the sub folder
#' "1_Normalization_\code{Name.folder.norm}".
#' Otherwise, a sub folder entitled "1_Normalization_\code{Name.folder.norm}"
#' will be created in \code{path.result} and all results will be saved in
#' "1_Normalization_\code{Name.folder.norm}".
#' If NULL, the results will not be saved in a folder. NULL as default.
#' @param Name.folder.norm Character or \code{NULL}.
#' If \code{Name.folder.norm} is a character,
#' the folder name which will contain the results will be
#' "1_Normalization_\code{Name.folder.norm}".
#' Otherwise, the folder name will be "1_Normalization".
#'
#' @return The function returns a normalized count data.frame and
#' plots a boxplot (if \code{Plot.Boxplot=TRUE}).
#'
#' @seealso The [DATAnormalization()] function calls the R functions
#' [BiocGenerics::estimateSizeFactors()], [DESeq2::rlog()] and [DESeq2::vst()]
#' in order to realized the normalization.
#'
#' @importFrom DESeq2 vst varianceStabilizingTransformation rlog
#' estimateSizeFactors counts
#' @importFrom SummarizedExperiment assay
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' data(RawCounts_Antoszewski2022_MOUSEsub500)
#' res.Norm<-DATAnormalization(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=NULL,
#'                             Individual.position=2,
#'                             Normalization="rle",
#'                             Blind.rlog.vst=FALSE,
#'                             Plot.Boxplot=TRUE,
#'                             Colored.By.Factors=FALSE,
#'                             Color.Group=NULL,
#'                             Plot.genes=FALSE,
#'                             path.result=NULL,
#'                             Name.folder.norm=NULL)
#' print(res.Norm)

DATAnormalization<-function(RawCounts,
                            Column.gene,
                            Group.position,
                            Time.position,
                            Individual.position,
                            Normalization="vst",
                            Blind.rlog.vst=FALSE,
                            Plot.Boxplot=TRUE,
                            Colored.By.Factors=FALSE,
                            Color.Group=NULL,
                            Plot.genes=FALSE,
                            path.result=NULL,
                            Name.folder.norm=NULL){
  #---------------------------------------------------------------------------#
  # Sample must belong to different biological condition or time points
  if(is.null(Time.position)==TRUE & is.null(Group.position)==TRUE){
    stop("'Time.position' and 'Group.position' can not be both NULL")
  }# (is.null(Time.position)==TRUE & is.null(Group.position)==TRUE)
  #---------------------------------------------------------------------------#
  # Different normalization authorized
  if(Normalization%in%c("rlog","vst","rle")==FALSE){
    stop("Normalization mut be 'vst', 'rlog' or 'rle'")
  }# (Normalization%in%c("rlog","vst","rle")==FALSE)
  #---------------------------------------------------------------------------#
  # Folder creation if no existence
  #---------------------------------------------------------------------------#
  if(is.null(Name.folder.norm)==TRUE){
    Name.folder.norm<-""
    SubFolder.name<-"1_UnsupervisedAnalysis"
  }else{
    Name.folder.norm<-paste("_",Name.folder.norm,sep="")
    SubFolder.name<-paste("1_UnsupervisedAnalysis",Name.folder.norm,sep="")
  }# if(is.null(Name.folder.norm)==TRUE)
  #
  if(is.null(path.result)==FALSE){
    if(SubFolder.name%in%dir(path=path.result)==FALSE){
      print("Folder creation")
      dir.create(path=paste(path.result,"/",SubFolder.name,sep=""))
      path.result.f<-paste(path.result,"/",SubFolder.name,sep="")
    }else{
      path.result.f<-paste(path.result,"/",SubFolder.name,sep="")
    }# if(SubFolder.name%in%dir(path = path.result)==FALSE)
  }else{
    path.result.f<-NULL
  }# if(is.null(path.result)==FALSE)
  #
  if(is.null(path.result.f)==FALSE){
    nom.dossier.result<-paste("1-1_Normalization",Name.folder.norm,sep="")
    if(nom.dossier.result%in%dir(path=path.result.f)==FALSE){
      dir.create(path=paste(path.result.f,"/",nom.dossier.result,sep=""))
      path.result.new<-paste(path.result.f,"/",nom.dossier.result,sep="")
    }else{
      path.result.new<-paste(path.result.f,"/",nom.dossier.result,sep="")
    }# if(nom.dossier.result%in%dir(path = path.result.f)==FALSE)
  }else{
    path.result.new<-NULL
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Name of all genes
  if(is.null(Column.gene)==TRUE){
    if(is.null(row.names(RawCounts))==TRUE){
      Name.G<-as.character(seq_len(nrow(RawCounts)))# 1:nrow(RawCounts)
    }else{
      Name.G<-row.names(RawCounts)
    }# is.null(row.names(RawCounts))==TRUE
  }else{
    Name.G<-RawCounts[,Column.gene]
  }# if(is.null(Column.gene)==TRUE)
  #---------------------------------------------------------------------------#
  # Columns with only expression
  if(is.null(Column.gene)==TRUE){
    ind.col.expr<-seq_len(ncol(RawCounts))#c(1:ncol(RawCounts))
  }else{
    ind.col.expr<-seq_len(ncol(RawCounts))[-Column.gene]
  }# if(is.null(Column.gene)==TRUE)
  #---------------------------------------------------------------------------#
  # Preprocessing
  resPreProcessing<-DEanalysisPreprocessing(RawCounts=RawCounts,
                                            Column.gene=Column.gene,
                                            Group.position=Group.position,
                                            Time.position=Time.position,
                                            Individual.position=Individual.position)
  DESeq2.obj<-resPreProcessing$DESeq2.obj
  FactorBoxplt<-resPreProcessing$DESeq2.obj@colData
  #---------------------------------------------------------------------------#
  # Normalization vst
  if(Normalization=="vst"){
    if(nrow(RawCounts)>1000){
      log2.data.prepare<-DESeq2::vst(DESeq2.obj,blind=Blind.rlog.vst)
    }else{
      log2.data.prepare<-DESeq2::varianceStabilizingTransformation(DESeq2.obj,
                                                                   blind=Blind.rlog.vst)
    }# if(nrow(RawCounts)>1000)
    log2.data<-round(SummarizedExperiment::assay(log2.data.prepare),digits=3)
    #
    Norm.dat<-data.frame(Gene=Name.G,log2.data)
    row.names(Norm.dat)<-Name.G
    #
    Log2Trf<-FALSE
    YlabelNorm<-"vst normalized counts"
  }# if(Normalization=="vst")
  #---------------------------------------------------------------------------#
  # Normalization rlog
  if(Normalization=="rlog"){
    log2.data.prepare <- DESeq2::rlog(DESeq2.obj,blind=Blind.rlog.vst)
    log2.data<-round(SummarizedExperiment::assay(log2.data.prepare),digits=3)
    #
    Norm.dat<-data.frame(Gene=Name.G,log2.data)
    row.names(Norm.dat)<-Name.G
    #
    Log2Trf<-FALSE
    YlabelNorm<-"rlog normalized counts"
  }# if(Normalization=="rlog")
  #---------------------------------------------------------------------------#
  # Normalization rle
  if(Normalization=="rle"){
    dds.SF<-DESeq2::estimateSizeFactors(DESeq2.obj)
    rle.data<-round(DESeq2::counts(dds.SF, normalized=TRUE),digits=3)
    Norm.dat<-data.frame(Gene=Name.G,rle.data)
    row.names(Norm.dat)<-Name.G
    #
    Log2Trf<-TRUE
    YlabelNorm<-"log2(rle normalized counts+1)"
  }# if(Normalization=="rle")
  colnames(Norm.dat)[-1]<-colnames(RawCounts)[ind.col.expr]
  #---------------------------------------------------------------------------#
  # Boxplot
  res.bxplt<-DATAplotBoxplotSamples(ExprData=Norm.dat,
                                    Column.gene=Column.gene,
                                    Group.position=Group.position,
                                    Time.position=Time.position,
                                    Individual.position=Individual.position,
                                    Log2.transformation=Log2Trf,
                                    Colored.By.Factors=Colored.By.Factors,
                                    Color.Group=Color.Group,
                                    Plot.genes=Plot.genes,
                                    y.label=YlabelNorm)
  #---------------------------------------------------------------------------#
  # Save of the boxplot graph
  if(is.null(path.result)==FALSE){
    grDevices::pdf(file = paste(path.result.new,"/",Normalization,
                                "_NormalizedBoxplotSamples.pdf",sep=""),
                   width = 11, height = 8)#width = 8, height = 11
    print(res.bxplt)
    grDevices::dev.off()
    #
    if(Plot.Boxplot==TRUE){
      print(res.bxplt)
    }# if(Plot.Boxplot==TRUE)
    # Save of the normalized data
    utils::write.table(Norm.dat,
                       file=paste(path.result.new,"/",
                                  Normalization,"_NormalizedData",
                                  ".csv",sep=""),
                       sep=";",row.names = FALSE)
  }else{
    if(Plot.Boxplot==TRUE){
      print(res.bxplt)
    }# if(Plot.Boxplot==TRUE)
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  return(list(NormalizedData=Norm.dat,
              NormalizedBoxplot=res.bxplt))
}# DATAnormalization()
