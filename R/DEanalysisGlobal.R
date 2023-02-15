#' @title Realization of the DE analysis (Main Function).
#'
#' @description The function realizes the DE analysis in three cases:
#' either samples belonging to different time measurements,
#' or samples belonging to different biological conditions,
#' or samples belonging to different time measurements
#' and different biological conditions.
#'
#' @details The column names of \code{ExprData} must be a vector of strings
#' of characters containing
#' * a string of characters (if \eqn{k=1}) which is the label of the column
#' containing gene names.
#' * \eqn{N_s} sample names which must be strings of characters
#' containing at least :
#' the name of the individual (e.g patient, mouse, yeasts culture),
#' its biological condition (if there is at least two) and
#' the time where data have been collected if there is at least two;
#' (must be either 't0', 'T0' or '0' for time 0,
#' 't1', 'T1' or '1' for time 1, ...).
#'
#' All these sample information must be separated by underscores
#' in the sample name. For instance 'CLL_P_t0_r1',
#' corresponds to the patient 'r1' belonging to the biological condition 'P'
#' and where data were collected at time 't0'.
#' I this example, 'CLL' describe the type of cells
#' (here chronic lymphocytic leukemia) and is not used in our analysis.
#'
#' In the string of characters 'CLL_P_t0_r1',
#' 'r1' is localized after the third underscore, so \code{Individual.position=4},
#' 'P' is localized after the first underscore, so \code{Group.position=2} and
#' 't0' is localized after the second underscore, so \code{Time.position=3}.
#'
#' @param RawCounts Data.frame with \eqn{N_g} rows and (\eqn{N_{s+k}}) columns,
#' where \eqn{N_g} is the number of genes,
#' \eqn{N_s} is the number of samples and
#' \eqn{k=1} if a column is used to specify gene names, or \eqn{k=0} otherwise.
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
#' @param pval.min Numeric value between 0 and 1. A gene will be considered as
#' differentially expressed (DE) between two biological conditions if
#' its Benjamini-Hochberg adjusted p-value (see [stats::p.adjust()]) is below
#' the threshold \code{pval.min}. Default value is 0.05.
#' @param pval.vect.t \code{NULL} or vector of dimension \eqn{T-1} filled with
#' numeric values between 0 and 1, with \eqn{T} the number of time measurements.
#' A gene will be considered as differentially expressed (DE) between
#' the time ti and the reference time t0 if its Benjamini-Hochberg adjusted
#' p-value (see [stats::p.adjust()]) is below the i-th threshold
#' of \code{pval.vect.t}.
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
#' If \code{Name.folder.DE} is a character, the folder names which will contain
#' all results will be "DEanalysis_\code{Name.folder.DE}".
#' Otherwise, the folder name will be "DEanalysis".
#'
#' @return The function [DEanalysisGlobal()] returns first the raw counts and
#' the rle normalized data automatically realized by [DESeq2::DESeq()]
#' (output \code{List.Datas}). Then
#' * If samples belong to different biological conditions, the function returns
#'   * a data.frame (output \code{DE.results}) which contains
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
#'     condition. A '1' in one of these columns means the gene is up-regulated
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
#'   and "DownRugulated", per biological condition (see [DEplotBarplot()]).
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
#'   information and inputs of [DEanalysisGlobal()].
#'
#' * If data belong to different time points only, the function returns
#'   * a data.frame (output \code{Results}) which contains
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
#'   each temporal group (see [DEplotAlluvial()]).
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
#'   information and inputs of [DEanalysisGlobal()].
#'
#'
#'
#'
#'
#'

#'
#' * If data belong to different time points and different biological conditions,
#' the function returns
#'   * a data.frame (output \code{Results}) which contains
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
#'     biological condition (see [DEplotBarplotFacetGrid()]).
#'     * \eqn{N_{bc}} alluvial graphs of DE genes (see [DEplotAlluvial()]),
#'     one per biological condition.
#'     * \eqn{N_{bc}} graphs showing the number of DE genes as a function of
#'     time for each temporal group, one per biological condition.
#'     By temporal group, we mean the sets of genes which are first DE at the
#'     same time.
#'     * \eqn{2\times N_{bc}} UpSet plot showing the number of DE genes
#'     belonging to each DE temporal pattern, for each biological condition.
#'     By temporal pattern, we mean the set of times ti such that the gene is
#'     DE between ti and the reference time t0 (see [DEplotVennBarplotTime()]).
#'     * an alluvial graph for DE genes which are DE at least one time for
#'     each group.
#'   * the following plots from the statistical analysis by biological
#'   condition
#'     * a barplot which gives the number of specific DE genes for each
#'     biological condition and time (see [DEplotBarplotFacetGrid()]).
#'     * \eqn{N_{bc}(N_{bc}-1)/2} UpSet plot which give the number of genes
#'     for each possible intersection (set of pairs of biological conditions),
#'     one per time (see [DEplotVennBarplotGroup()]).
#'     * an alluvial graph of genes which are specific at least one time
#'     (see [DEplotAlluvial()]).
#'   * the following plots from the combination of temporal and biological
#'   statistical analysis
#'     * a barplot which gives the number of signature genes for each
#'     biological condition and time (see [DEplotBarplotFacetGrid()]).
#'     * a barplot showing the number of genes which are DE at at least one
#'     time, specific at at least one time and signature at at least one time,
#'     for each biological condition.
#'     * an alluvial graph of genes which are signature at least one time
#'     (see [DEplotAlluvial()]).
#'
#'
#' @importFrom DESeq2 DESeq counts varianceStabilizingTransformation
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' data(RawCounts_Antoszewski2022_MOUSEsub500)
#' # No time points. We take only two groups for the speed of the example
#' RawCounts_T1Wt<-RawCounts_Antoszewski2022_MOUSEsub500[,1:7]
#' #
#' res.all=DEanalysisGlobal(RawCounts=RawCounts_T1Wt,
#'                          Column.gene=1,
#'                          Group.position=1,
#'                          Time.position=NULL,
#'                          Individual.position=2,
#'                          pval.min=0.05,
#'                          pval.vect.t=NULL,
#'                          log.FC.min=1,
#'                          LRT.supp.info=FALSE,
#'                          path.result=NULL,
#'                          Name.folder.DE=NULL)

DEanalysisGlobal<-function(RawCounts,
                           Column.gene,
                           Group.position,
                           Time.position,
                           Individual.position,
                           pval.min=0.05,
                           pval.vect.t=NULL,
                           log.FC.min=1,
                           LRT.supp.info=FALSE,
                           path.result=NULL,
                           Name.folder.DE=NULL){
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Creation folder if no existence
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  if(is.null(Name.folder.DE)==TRUE){
    Name.folder.DE.ini<-NULL
    Name.folder.DE<-""
    SubFolder.name<-"2_SupervisedAnalysis"
  }else{
    Name.folder.DE.ini<-Name.folder.DE
    Name.folder.DE<-paste("_",Name.folder.DE.ini,sep="")
    SubFolder.name<-paste("2_SupervisedAnalysis", Name.folder.DE,sep="")
  }# if(is.null(Name.folder.DE)==TRUE)
  #
  if(is.null(path.result)==FALSE){
    if(SubFolder.name%in%dir(path=path.result)==FALSE){
      print("Folder creation")
      dir.create(path=paste(path.result,"/",SubFolder.name,sep=""))
      path.result.f<-paste(path.result,"/",SubFolder.name,sep="")
    }else{
      path.result.f<-paste(path.result,"/",SubFolder.name,sep="")
    }# if(SubFolder.name%in%dir(path=path.result)==FALSE)
  }else{
    path.result.f<-NULL
  }# if(is.null(path.result)==FALSE)
  # Folder for RLE normalized count data
  if(is.null(path.result)==FALSE){
    name.folder.result1<-paste("2-1_RLEnormalizedDATA",Name.folder.DE,sep="")
    if(name.folder.result1%in%dir(path = path.result.f)==FALSE){
      dir.create(path=paste(path.result.f,"/",name.folder.result1,sep=""))
      path.result.new1<-paste(path.result.f,"/",name.folder.result1,sep="")
    }else{
      path.result.new1<-paste(path.result.f,"/",name.folder.result1,sep="")
    }# if(name.folder.result1%in%dir(path = path.result.f)==FALSE)
  }else{
    path.result.new1<-NULL
  }# if(is.null(path.result)==FALSE)
  # Folder for DE results csv
  if(is.null(path.result)==FALSE){
    name.folder.result3<-paste("2-3_CSV_file_DEanalysis",Name.folder.DE,sep="")
    if(name.folder.result3%in%dir(path = path.result.f)==FALSE){
      dir.create(path=paste(path.result.f,"/",name.folder.result3,sep=""))
      path.result.new3<-paste(path.result.f,"/",name.folder.result3,sep="")
    }else{
      path.result.new3<-paste(path.result.f,"/",name.folder.result3,sep="")
    }# if(name.folder.result2%in%dir(path = path.result.f)==FALSE)
  }else{
    path.result.new3<-NULL
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Pre-processing
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  print("Preprocessing")
  resPreProcessing<-DEanalysisPreprocessing(RawCounts=RawCounts,
                                            Column.gene=Column.gene,
                                            Group.position=Group.position,
                                            Time.position=Time.position,
                                            Individual.position=Individual.position)
  DESeq2.obj<-resPreProcessing$DESeq2.obj
  #---------------------------------------------------------------------------#
  # columns with only expression
  if(is.null(Column.gene)==TRUE){
    ind.col.expr<-seq_len(ncol(RawCounts))#c(1:ncol(RawCounts))
    Data.f<-cbind(Gene=paste("Gene",seq_len(nrow(RawCounts)),sep=""),RawCounts)
  }else{
    ind.col.expr<-seq_len(ncol(RawCounts))[-Column.gene]
    Data.f<-cbind(Gene=RawCounts[,Column.gene], RawCounts[,-Column.gene])
  }# if(is.null(Column.gene)==TRUE)
  row.names(Data.f)<-Data.f[,1]
  #---------------------------------------------------------------------------#
  All.Datas<-vector("list", length = 2)
  names(All.Datas)<-c("RawCounts", "RLEdata")
  All.Datas[[1]]<-Data.f
  #---------------------------------------------------------------------------#
  # if(is.null(path.result)==FALSE){
  #   utils::write.table(Data.f,
  #                      file=paste(path.result,"/",SubFolder.name,"/",
  #                                 "DATA_RawCounts",Name.folder.DE,".csv",sep=""),
  #                      sep=";",row.names = FALSE)
  # }
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Differential expression
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  print("Differential expression step with DESeq2::DESeq()")
  if(LRT.supp.info==TRUE){
    dds.norm.diff<-DESeq2::DESeq(DESeq2.obj,
                                 betaPrior=FALSE, test="LRT", reduced=~1)
  }else{
    dds.norm.diff<-DESeq2::DESeq(DESeq2.obj, betaPrior=FALSE, test="Wald")
  }# if(LRT.supp.info==TRUE)
  #---------------------------------------------------------------------------#
  ScaledData<-round(DESeq2::counts(dds.norm.diff, normalized=TRUE), digits=3)
  ScaledData.f<-data.frame(Gene=Data.f[,1], ScaledData)
  All.Datas[[2]]<-ScaledData.f
  colnames(ScaledData.f)<-c("Gene", colnames(RawCounts)[ind.col.expr])
  #---------------------------------------------------------------------------#
  res.bxplt<-DATAplotBoxplotSamples(ExprData=ScaledData.f,
                                    Column.gene=1,
                                    Group.position=Group.position,
                                    Time.position=Time.position,
                                    Individual.position=Individual.position,
                                    Log2.transformation=TRUE,
                                    Colored.By.Factors=FALSE,
                                    Color.Group=NULL,
                                    Plot.genes=FALSE,
                                    y.label="log2(rle normalized counts+1)")
  #
  if(is.null(path.result)==FALSE){
    utils::write.table(ScaledData.f,
                       file=paste(path.result.new1,"/",
                                  "DATA_RLE",Name.folder.DE,".csv",sep=""),
                       sep=";", row.names = FALSE)
    #
    grDevices::pdf(file=paste(path.result.new1,"/BoxplotSamples_RLE",
                              Name.folder.DE,".pdf",sep=""),
                   width=11, height=8)#width = 8, height = 11
    print(res.bxplt)
    grDevices::dev.off()
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  if(is.null(Time.position)==FALSE){
    Levels.time<-levels(as.factor(dds.norm.diff@colData$Time))
    Nb.time<-length(Levels.time)
    #
    if(is.null(pval.vect.t)==FALSE){
      if(length(pval.vect.t)>Nb.time-1){
        pval.vect.t<-pval.vect.t[seq_len(Nb.time-1)]
      }# if(length(pval.vect.t)>Nb.time-1)
      if(length(pval.vect.t)<Nb.time-1){
        pval.vect.t<-c(pval.vect.t,
                       rep(pval.min,times=Nb.time-length(pval.vect.t)-1))
      }# if(length(pval.vect.t)<Nb.time-1)
    }else{
      pval.vect.t<-rep(pval.min, times=Nb.time-1)
    }# if(is.null(pval.vect.t)==FALSE)
  }# if(is.null(Time.position)==FALSE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Results from Differential expression step with DESeq2::DESeq()
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  if(is.null(Time.position)==TRUE & is.null(Group.position)==TRUE){
    stop("'Time.position' and 'Group.position' can not be both NULL")
  }# if(is.null(Time.position)==TRUE & is.null(Group.position)==TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Case 1 analysis DE : Biological conditions only
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  if(is.null(Time.position)==TRUE & is.null(Group.position)==FALSE){
    print("Case 2 analysis : Biological conditions only")
    # Folder for DE graphs
    if(is.null(path.result)==FALSE){
      name.folder.result2<-paste("2-2_Group_DEanalysis", Name.folder.DE,sep="")
      if(name.folder.result2%in%dir(path = path.result.f)==FALSE){
        dir.create(path=paste(path.result.f,"/",name.folder.result2,sep=""))
        path.result.new2<-paste(path.result.f,"/",name.folder.result2,sep="")
      }else{
        path.result.new2<-paste(path.result.f,"/",name.folder.result2,sep="")
      }# if(name.folder.result2%in%dir(path = path.result.f)==FALSE)
    }else{
      path.result.new2<-NULL
    }# if(is.null(path.result)==FALSE)
    #
    Res.DE.BC<-DEanalysisGroup(DESeq.result=dds.norm.diff,
                               LRT.supp.info=LRT.supp.info,
                               log.FC.min=log.FC.min,
                               pval.min=pval.min,
                               path.result=path.result.new2,
                               SubFile.name=SubFolder.name)
    #-------------------------------------------------------------------------#
    # Table which contains all results.
    #-------------------------------------------------------------------------#
    if(is.null(path.result)==FALSE){
      utils::write.table(data.frame(Res.DE.BC$Results),
                         file=paste(path.result.new3,"/",
                                    "ALLresults_DEanalysis", Name.folder.DE,
                                    ".csv", sep=""),
                         sep=";", row.names=FALSE)
    }# if(is.null(path.result)==FALSE)
    #
    SumInfo<-list(ExprCond=c("Group"),
                  FactorsInfo=resPreProcessing$Factors.Info,
                  GroupLevels=levels(resPreProcessing$DESeq2.obj@colData$Group),
                  logFCmin=log.FC.min,
                  pvalGroup=pval.min)
    #
    resGlossary<-Glossary(path.result.new3, Case=1)
    #
    return(list(List.Datas=All.Datas,
                Summary.Inputs=SumInfo,
                DE.results=Res.DE.BC$Results,
                Path.result=path.result.f,
                Folder.result=Name.folder.DE.ini,
                List.Glossary=resGlossary,
                List.Plots.DE.Analysis=Res.DE.BC$List.Plots.DE.Group,
                DESeq.dds=dds.norm.diff))#SubFolder.name
  }# if(is.null(Time.position)==TRUE & is.null(Group.position)==FALSE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Case 2 analysis DE : Time only
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  if(is.null(Time.position)==FALSE & is.null(Group.position)==TRUE){
    print("Case 1 analysis : Time only")
    # Folder for DE graphs
    if(is.null(path.result)==FALSE){
      name.folder.result2<-paste("2-2_Temporal_DEanalysis", Name.folder.DE,
                                 sep="")
      if(name.folder.result2%in%dir(path = path.result.f)==FALSE){
        dir.create(path=paste(path.result.f,"/",name.folder.result2,sep=""))
        path.result.new2<-paste(path.result.f,"/",name.folder.result2,sep="")
      }else{
        path.result.new2<-paste(path.result.f,"/",name.folder.result2,sep="")
      }# if(name.folder.result2%in%dir(path = path.result.f)==FALSE)
    }else{
      path.result.new2<-NULL
    }# if(is.null(path.result)==FALSE)
    #
    Res.DE.Time<-DEanalysisTime(DESeq.result=dds.norm.diff,
                                LRT.supp.info=LRT.supp.info,
                                log.FC.min=log.FC.min,
                                pval.min=pval.min,
                                pval.vect.t=pval.vect.t,
                                path.result=path.result.new2,
                                SubFile.name=SubFolder.name)
    #-------------------------------------------------------------------------#
    # Table which contains all results
    #-------------------------------------------------------------------------#
    if(is.null(path.result)==FALSE){
      utils::write.table(data.frame(Res.DE.Time$Results),
                         file=paste(path.result.new3,"/",
                                    "ALLresults_DEanalysis", Name.folder.DE,
                                    ".csv", sep=""),
                         sep=";", row.names=FALSE)
    }# if(is.null(path.result)==FALSE)
    #
    FactorInfo.f<-resPreProcessing$Factors.Info
    Tnumeric<-gsub("t","",gsub("T","",FactorInfo.f$Time,fixed=TRUE),fixed=TRUE)
    TlevNumeric<-gsub("t","",
                      gsub("T","",
                           levels(resPreProcessing$DESeq2.obj@colData$Time),
                           fixed=TRUE), fixed=TRUE)
    FactorInfo.f$Time<-paste("t",Tnumeric,sep="")
    TlevNumeric<-paste("t",Tnumeric,sep="")
    #
    SumInfo<-list(ExprCond=c("Time"),
                  FactorsInfo=FactorInfo.f,
                  TimeLevels=levels(factor(TlevNumeric)),
                  logFCmin=log.FC.min,
                  pvalsTime=pval.vect.t)
    #
    resGlossary<-Glossary(path.result.new3, Case=2)
    #
    return(list(List.Datas=All.Datas,
                Summary.Inputs=SumInfo,
                DE.results=Res.DE.Time$Results,
                Path.result=path.result.f,
                Folder.result=Name.folder.DE.ini,
                List.Glossary=resGlossary,
                List.Plots.DE.Analysis=Res.DE.Time$List.Plots.DE.Time,
                DESeq.dds=dds.norm.diff))#SubFolder.name
  }# if(is.null(Time.position)==FALSE & is.null(Group.position)==TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Case 3 analysis DE : Time and Biological conditions
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  if(is.null(Time.position)==FALSE | is.null(Group.position)==FALSE){
    print("Case 3 analysis : Biological conditions and Times.")
    # Folder for DE graphs
    if(is.null(path.result)==FALSE){
      name.folder.result2<-paste("2-2_DEanalysis",Name.folder.DE,sep="")
      if(name.folder.result2%in%dir(path = path.result.f)==FALSE){
        dir.create(path=paste(path.result.f,"/",name.folder.result2,sep=""))
        path.result.new2<-paste(path.result.f,"/",name.folder.result2,sep="")
      }else{
        path.result.new2<-paste(path.result.f,"/",name.folder.result2,sep="")
      }# if(name.folder.result2%in%dir(path = path.result.f)==FALSE)
    }else{
      path.result.new2<-NULL
    }# if(is.null(path.result)==FALSE)
    #
    Res.DE.T.G<-DEanalysisTimeAndGroup(DESeq.result=dds.norm.diff,
                                       LRT.supp.info=LRT.supp.info,
                                       log.FC.min=log.FC.min,
                                       pval.min=pval.min,
                                       pval.vect.t=pval.vect.t,
                                       path.result=path.result.new2,
                                       SubFile.name=SubFolder.name)
    #-------------------------------------------------------------------------#
    # Tables which contain all results
    #-------------------------------------------------------------------------#
    if(is.null(path.result)==FALSE){
      utils::write.table(data.frame(Res.DE.T.G$Results),
                         file=paste(path.result.new3,"/",
                                    "ALLresults_DEanalysis", Name.folder.DE,
                                    ".csv", sep=""),
                         sep=";", row.names=FALSE)
    }# if(is.null(path.result)==FALSE)
    #
    FactorInfo.f<-resPreProcessing$Factors.Info
    Tnumeric<-gsub("t","",gsub("T","",FactorInfo.f$Time,fixed=TRUE),fixed=TRUE)
    TlevNumeric<-gsub("t","",
                      gsub("T","",
                           levels(resPreProcessing$DESeq2.obj@colData$Time),
                                  fixed=TRUE), fixed=TRUE)
    FactorInfo.f$Time<-paste("t",Tnumeric,sep="")
    TlevNumeric<-paste("t",Tnumeric,sep="")
    #
    SumInfo<-list(ExprCond=c("Time","Group"),
                  FactorsInfo=FactorInfo.f,
                  TimeLevels=levels(factor(TlevNumeric)),
                  GroupLevels=levels(resPreProcessing$DESeq2.obj@colData$Group),
                  logFCmin=log.FC.min,
                  pvalsTime=pval.vect.t,
                  pvalGroup=pval.min)
    #
    resGlossary<-Glossary(path.result.new3, Case=3)
    #-------------------------------------------------------------------------#
    return(list(List.Datas=All.Datas,
                Summary.Inputs=SumInfo,
                DE.results=data.frame(Res.DE.T.G$Results),
                Path.result=path.result.f,
                Folder.result=Name.folder.DE.ini,
                List.Glossary=resGlossary,
                List.Plots.DE.Analysis=Res.DE.T.G$List.Plots.DE.Time.Group,
                DESeq.dds=dds.norm.diff))
  }# if(is.null(Time.position)==FALSE | is.null(Group.position)==FALSE)
}# DEanalysisGlobal()


#-----------------------------------------------------------------------------#
#=============================================================================#
#-----------------------------------------------------------------------------#

Glossary<-function(path.result, Case){
  #---------------------------------------------------------------------------#
  if(is.null(path.result)==FALSE){
    # path for glossary
    fileRdme<-paste(path.result,"/", "Glossary.txt",sep="")
    # Title of the text file
    cat(rep("=",times=100), "\n", sep="", file=fileRdme, append=FALSE)
    cat(rep("=",times=100), "\n", sep="", file=fileRdme, append=TRUE)
    cat(c("=", rep("-",times=15), "= GLOSSARY ",
          "of column names in the cvs file 'ALLresults_DEanalysis' =",
          rep("-",times=15),"=\n"),
        sep="", file=fileRdme, append=TRUE)
    cat(rep("=",times=100), "\n", sep="", file=fileRdme, append=TRUE)
    cat(rep("=",times=100), "\n", sep="", file=fileRdme, append=TRUE)
    #
    cat("\n", sep="", file=fileRdme, append=TRUE)
    cat(c("The cvs file 'ALLresults_DEanalysis' gathers all the results of ",
          "MultiRNAflow analysis. The goal of this glossary is to clarify ",
          "the meaning of the column names in the file. ",
          "The first column gives the names of all genes. ",
          "In the list below, we specify the meaning of each entry ",
          "in the corresponding column.\n"),
        sep="", file=fileRdme, append=TRUE)
    cat(rep("=",times=100), "\n", sep="", file=fileRdme, append=TRUE)
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Glossary when samples depends on biological condition only
  if(Case==1){
    ListGlossary<-vector(mode="list", length=6)
    names(ListGlossary)<-c("Log2FoldChange.Group2.versus.Group1",
                           "Pvalue.adjusted.Group2.versus.Group1",
                           "DE.Group2.versus.Group1",
                           "DE.1pair.of.Group.minimum",
                           "Specific.genes_Group1",
                           "OverUnder.regulated.genes_Group1")
    #
    ListGlossary[[1]]<-paste("Log2 fold change between ",
                             "the biological condition Group2 and ",
                             "the biological condition Group1.", sep="")
    ListGlossary[[2]]<-paste("Adjusted p-value between ",
                             "the biological condition Group2 and ",
                             "the biological condition Group1.", sep="")
    ListGlossary[[3]]<-paste("Binary number (0 or 1) where 1 means that ",
                             "the gene is DE and 0 means that it is not DE ",
                             "between the biological condition Group2 and ",
                             "the biological condition Group1. ",
                             "The value 1 is given if both ",
                             "'abs(Log2FoldChange.Group2.versus.Group1)",
                             ">log.FC.min' and ",
                             "'Pvalue.adjusted.Group2.versus.Group1<pval.min',",
                             " where 'log.FC.min' and 'pval.min' are inputs ",
                             "of the function DEanalysisGlobal()", sep="")
    ListGlossary[[4]]<-paste("Binary number (0 or 1) where 1 means that ",
                             "the gene is DE between at least one pair of ",
                             "biological conditions.", sep="")
    ListGlossary[[5]]<-paste("Binary number (0 or 1) where 1 means that the ",
                             "gene is specific for the biological condition ",
                             "Group1. This means that the gene is DE between ",
                             "Group1 and any other biological conditions ",
                             "but not DE between any pairs of other ",
                             "biological conditions.", sep="")
    ListGlossary[[6]]<-paste("Number (-1, 0 or 1). 1 means the gene is ",
                             "specific ('Specific.genes_Group1=1') and the ",
                             "gene is up-regulated (or over expressed) in ",
                             "Group1 versus the other biological conditions. ",
                             "-1' means the gene is specific ",
                             "('Specific.genes_Group1=1') and the gene is ",
                             "down-regulated (or under expressed) in Group1 ",
                             "versus the other biological conditions. ",
                             "0 otherwise.", sep="")
  }# if(Case==1)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Glossary when samples depends on time only
  if(Case==2){
    ListGlossary<-vector(mode="list", length=5)
    names(ListGlossary)<-c("Log2FoldChange.ti.versus.t0",
                           "Pvalue.adjusted.ti.versus.t0",
                           "DE.ti.versus.t0",
                           "DE.Temporal.Pattern",
                           "DE.1time.minimum")
    #
    ListGlossary[[1]]<-paste("Log2 fold change between the time ti and ",
                             "the reference time t0.", sep="")
    ListGlossary[[2]]<-paste("Adjusted pvalue between the time ti and ",
                             "the reference time t0.", sep="")
    ListGlossary[[3]]<-paste("Binary number (0 or 1) where 1 means that the ",
                             "gene is DE and 0 means it is not DE between ",
                             "the time ti and the reference time t0. ",
                             "The value 1 is given if both ",
                             "'abs(Log2FoldChange.ti.versus.t0)>log.FC.min' ",
                             "and ",
                             "'Pvalue.adjusted.ti.versus.t0< pval.vect.t[i]'.",
                             sep="")
    ListGlossary[[4]]<-paste("Vector of 0 and 1 corresponding to times ",
                             "t1 to tn. The values 1 correspond to the times ",
                             "ti such that the gene is DE between ",
                             "ti and the reference time t0.", sep="")
    ListGlossary[[5]]<-paste("Binary number (0 or 1) where 1 means that the ",
                             "gene is DE at least between one time ti ",
                             "versus the reference time t0.", sep="")
  }#if(Case==2)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # Glossary when samples depends on time and biological condition
  if(Case==3){
    ListGlossary<-vector(mode="list", length=14)
    names(ListGlossary)<-c("Log2FoldChange.ti.versus.t0_Group1",
                           "Pvalue.adjusted.ti.versus.t0_Group1",
                           "DE.ti.versus.t0_Group1",
                           "DE.Temporal.Pattern_Group1",
                           "DE.1time.minimum_Group1",#
                           "Log2FoldChange.Group2.versus.Group1_Time.ti",
                           "Pvalue.adjusted.Group2.versus.Group1_Time.ti",
                           "DE.Group2.versus.Group1_Time.ti",
                           "DE.1pair.of.Group.minimum_Time.ti",
                           "Specific.genes_Group1_Time.ti",
                           "OverUnder.regulated.genes_Group1_Time.ti",
                           "Specific.genes_Group1_1t.minimum",
                           "Signature.genes_Group.Group1_Time.ti",
                           "Signature.genes_Group.Group1_1time.minimum")
    #
    ListGlossary[[1]]<-paste("Log2 fold change between the time ti and the ",
                             "reference time t0, for the biological condition",
                             " Group1.", sep="")
    ListGlossary[[2]]<-paste("Adjusted pvalue between the time ti and the ",
                             "reference time t0, for the biological condition",
                             " Group1.", sep="")
    ListGlossary[[3]]<-paste("Binary number (0 or 1) where 1 means that the ",
                             "gene is DE and 0 means it is not DE, for the ",
                             "biological condition Group1. ",
                             "The value 1 is given if both ",
                             "'abs(Log2FoldChange.ti.versus.t0_Group1)>",
                             "log.FC.min' and ",
                             "'Pvalue.adjusted.ti.versus.t0_Group1<",
                             "pval.vect.t[i]', for the group Group1.", sep="")
    ListGlossary[[4]]<-paste("Vector of 0 and 1 corresponding to times ",
                             "t1 to tn. The values 1 correspond to the times ",
                             "ti such that the gene is DE between ",
                             "ti and the reference time t0, ",
                             "for the biological condition Group1.", sep="")
    ListGlossary[[5]]<-paste("Binary number (0 or 1) where 1 means that the ",
                             "gene is DE at least between one time ti ",
                             "versus the reference time t0, ",
                             "for the group Group1.", sep="")
    #
    ListGlossary[[6]]<-paste("Log2 fold change between ",
                             "the biological condition Group2 and ",
                             "the biological condition Group1, ",
                             "at time ti.",sep="")
    ListGlossary[[7]]<-paste("Adjusted p-value between ",
                             "the biological condition Group2 and ",
                             "the biological condition Group1, ",
                             "at time ti.", sep="")
    ListGlossary[[8]]<-paste("Binary number (0 or 1) where 1 means that ",
                             "the gene is DE and 0 means that it is not DE ",
                             "between the biological condition Group2 and ",
                             "the biological condition Group1, at time ti. ",
                             "The value 1 is given if both 'abs(",
                             "Pvalue.adjusted.Group2.versus.Group1_Time.ti)",
                             ">log.FC.min' and ",
                             "'Pvalue.adjusted.Group2.versus.Group1_Time.ti",
                             "<pval.min', where 'log.FC.min' and 'pval.min' ",
                             "are inputs of the function DEanalysisGlobal()",
                             sep="")
    ListGlossary[[9]]<-paste("Binary number (0 or 1) where 1 means that ",
                             "the gene is DE between at least one pair of ",
                             "biological conditions, at time ti.",
                             sep="")
    ListGlossary[[10]]<-paste("Binary number (0 or 1) where 1 means that the ",
                              "gene is specific for the biological condition ",
                              "Group1, at time ti. This means that the ",
                              "gene is DE between Group1 and any other ",
                              "biological conditions but not DE between ",
                              "any pairs of other biological conditions, ",
                              "at time ti.",
                              sep="")
    ListGlossary[[11]]<-paste("Number (-1, 0 or 1). 1 means the gene is ",
                              "specific ('Specific.genes_Group1_Time.ti=1') ",
                              "at time ti and the gene is up-regulated ",
                              "(or over expressed) in Group1 versus the ",
                              "other biological conditions, for the time ti. ",
                              "-1' means the gene is specific ",
                              "('Specific.genes_Group1_Time.ti=1') ",
                              "at time ti and the gene is down-regulated ",
                              "(or under expressed) in Group1 versus the ",
                              "other biological conditions, for the time ti. ",
                              "0 otherwise.", sep="")
    ListGlossary[[12]]<-paste("Binary number (0 or 1). 1 means that the gene ",
                              "is specific for the biological condition ",
                              "Group1, at at least one time ti.",
                              "0 otherwise.", sep="")
    #
    ListGlossary[[13]]<-paste("Binary number (0 or 1). 1 means that the gene ",
                              " is a signature gene for the biological ",
                              "condition Group1 at time ti. This means that ",
                              "the gene is specific for the biological ",
                              "condition Group1 at time ti and the gene is ",
                              "DE between the time ti and the reference time ",
                              "t0 for the group Group1.", sep="")
    ListGlossary[[14]]<-paste("Binary number (0 or 1). 1 means that the gene ",
                              "is a signature gene for the biological ",
                              "condition Group1 at at least one time ti. ",
                              "0 otherwise.", sep="")
  }# if(Case==3)
  #---------------------------------------------------------------------------#
  if(is.null(path.result)==FALSE){
    for(i in seq_len(length(ListGlossary))){# 1:length(ListGlossary)
      if(i%in%c(1,6,13) & Case==3){
        if(i==1){
          cat(c("==-----= ", "Temporal statistical analysis", "\n"),
              sep="", file=fileRdme, append=TRUE)
        }
        if(i==6){
          cat("\n", sep="", file=fileRdme, append=TRUE)
          cat(rep("=",times=100), "\n", sep="", file=fileRdme, append=TRUE)
          cat(c("==-----= ", "Statistical analysis by biological condition",
                "\n"), sep="", file=fileRdme, append=TRUE)
        }
        if(i==13){
          cat("\n", sep="", file=fileRdme, append=TRUE)
          cat(rep("=",times=100), "\n", sep="", file=fileRdme, append=TRUE)
          cat(c("==-----= ", "Combination of temporal and condition analysis",
                "\n"), sep="", file=fileRdme, append=TRUE)
        }
        # cat("\n", sep="", file=fileRdme, append=TRUE)
        # cat(rep("=",times=100), "\n", sep="", file=fileRdme, append=TRUE)
      }# if(i%in%c(6,12) & Case==3)
      cat("\n", sep="", file=fileRdme, append=TRUE)
      cat(c(" ** ",names(ListGlossary)[i], " : ", ListGlossary[[i]], "\n"),
          sep="", file=fileRdme, append=TRUE)
    }# for(i in 1:length(ListGlossary))
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  return(ListGlossary)
}# Glossary()
