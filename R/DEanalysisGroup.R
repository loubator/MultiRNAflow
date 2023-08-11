#' @title DE Analysis when samples belong to different biological conditions.
#'
#' @description The function realizes from the
#' [DESeq2::DESeq()]
#' output the analysis of DE genes between all pairs of biological conditions.
#'
#' @param DESeq.result Output from the function
#' [DESeq2::DESeq()].
#' @param pval.min Numeric value between 0 and 1. A gene will be considered as
#' differentially expressed (DE) between two biological conditions if
#' its Benjamini-Hochberg adjusted p-value
#' (see [stats::p.adjust()])
#' is below the threshold \code{pval.min}. Default value is 0.05.
#' @param log.FC.min Non negative numeric value.
#' If the log2 fold change between biological conditions or times has
#' an absolute value below the threshold \code{log.FC.min}, then the gene is
#' not selected even if is considered as DE. Default value is 1.
#' If \code{log.FC.min=0}, all DE genes will be kept.
#' @param LRT.supp.info \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the algorithm realizes another statistical test in order
#' to detect if, among all biological conditions and/or times, at least one
#' has a different behavior than the others (see the input 'test' in
#' [DESeq2::DESeq()]).
#' @param Plot.DE.graph \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, all graphs will be plotted.
#' Otherwise no graph will be plotted.
#' @param path.result Character or \code{NULL}.
#' If \code{path.result} is a character, it must be a path to a folder,
#' all graphs will be saved in \code{path.result}.
#' If \code{NULL}, the results will not be saved in a folder.
#' \code{NULL} as default.
#' @param SubFile.name Character or \code{NULL}.
#' If \code{SubFile.name} is a character, each saved file names will contain
#' the strings of characters "_\code{SubFile.name}".
#' If \code{NULL}, no suffix will be added.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom grDevices dev.off pdf
#'
#' @return The function returns
#' * a data.frame (output \code{Results}) which contains
#'   * gene names
#'   * pvalues, log2 fold change and DE genes between each pairs of
#'   biological conditions.
#'   * a binary column (1 and 0) where 1 means the gene is DE between at
#'   least one pair of biological conditions.
#'   * \eqn{N_{bc}} binary columns, where \eqn{N_{bc}} is the number of
#'   biological conditions, which gives the specific genes for each
#'   biological condition.
#'   A '1' in one of these columns means the gene is specific to the
#'   biological condition associated to the given column. 0 otherwise.
#'   A gene is called specific to a given biological condition BC1,
#'   if the gene is DE between BC1 and any other biological conditions,
#'   but not DE between any pair of other biological conditions.
#'   * \eqn{N_{bc}} columns filled with -1, 0 and 1, one per biological
#'   condition. A '1' in one of these columns means the gene is up-regulated
#'   (or over-expressed) for the biological condition associated to the
#'   given column. A gene is called up-regulated for a given biological
#'   condition BC1 if the gene is specific to the biological condition BC1
#'   and expressions in BC1 are higher than in the other
#'   biological conditions.
#'   A '-1' in one of these columns means the gene is down-regulated
#'   (or under-expressed) for the biological condition associated to the
#'   given column.
#'   A gene is called down-regulated for a given biological condition BC1 if
#'   the gene is specific to the biological condition BC1 and expressions in
#'   BC1 are lower than in the other biological conditions.
#'   A '0' in one of these columns means the gene is not specific to the
#'   biological condition associated to the given column.
#' * A contingency matrix (output \code{Summary.DEanalysis}) which gives
#' for each biological condition the number of genes categorized as
#' "Upregulated", "DownRugulated" and "Other".
#' A gene is categorized as 'Other', for a given biological condition,
#' if the gene is not specific to the given biological condition.
#' The category 'Other' does not exist when there are only two biological
#' conditions.
#' * an UpSet plot (Venn diagram displayed as a barplot) which gives the
#' number of genes for each possible intersection
#' (see [DEplotVennBarplotGroup()]).
#' We consider that a set of pairs of biological conditions forms an
#' intersection if there is at least one gene which is DE for each of
#' these pairs of biological conditions, but not for the others.
#' * a barplot which gives the number of genes categorized as "Upregulated"
#' and "DownRugulated", per biological condition
#' (see [DEplotBarplot()]).
#' * a barplot which gives the number of genes categorized as "Upregulated",
#' "DownRugulated" and "Other", per biological condition
#' (see [DEplotBarplot()]).
#' So this barplot, only plotted when there are strictly more than
#' two biological conditions, is similar to the previous barplot but with
#' the category "Other".
#'
#' @seealso The outputs of the function are used by the main function
#' [DEanalysisGlobal()].
#'
#' @export
#'
#' @examples
#' ## Data
#' data(RawCounts_Antoszewski2022_MOUSEsub500)
#' ## No time points. We take only two groups for the speed of the example
#' RawCounts_T1Wt<-RawCounts_Antoszewski2022_MOUSEsub500[,1:7]
#'
#' ## Preprocessing step
#' resDATAprepSEmus2<- DATAprepSE(RawCounts=RawCounts_T1Wt,
#'                                Column.gene=1,
#'                                Group.position=1,
#'                                Time.position=NULL,
#'                                Individual.position=2)
#'
#' ##------------------------------------------------------------------------#
#' dds.DE.G<-DESeq2::DESeq(resDATAprepSEmus2$DESeq2.obj,
#'                         quiet=TRUE, betaPrior=FALSE)
#' ##
#' res.sum.group<-DEanalysisGroup(DESeq.result=dds.DE.G,
#'                                pval.min=0.01,
#'                                log.FC.min=1,
#'                                LRT.supp.info=FALSE,
#'                                Plot.DE.graph=TRUE,
#'                                path.result=NULL,
#'                                SubFile.name=NULL)

DEanalysisGroup<-function(DESeq.result,
                          pval.min=0.05,
                          log.FC.min=1,
                          LRT.supp.info=TRUE,
                          Plot.DE.graph=TRUE,
                          path.result=NULL,
                          SubFile.name=NULL){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if(!is(DESeq.result, 'DESeqDataSet')){
        stop("Res.DE.analysis must be a 'DESeqDataSet' object")
    }## if(!is(classDeseq2, 'DESeqDataSet'))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 1) Summary DESeq2 results
    Fct.group <- data.frame(SummarizedExperiment::colData(DESeq.result))[,1]
    Fct.group<-as.factor(as.character(Fct.group))
    nb.group<-length(levels(Fct.group))
    nb.pair.of.group<-(nb.group*(nb.group-1))/2

    Sum.DE.analysis.G<-DEresultGroup(DESeq.result=DESeq.result,
                                     LRT.supp.info=LRT.supp.info,
                                     log.FC.min=log.FC.min,
                                     pval.min=pval.min)

    Cont.per.group<-Sum.DE.analysis.G$Contingence.per.group

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 2) Barplot and Upset plot

    if(nb.pair.of.group>1){
        List.plots.DE.group<-vector(mode="list", length=3)
        names(List.plots.DE.group)<-c("VennBarplot",
                                      "NumberDEgenes_SpecificAndNoSpecific",
                                      "NumberDEgenes_SpecificGenes")

        G.Upset<-DEplotVennBarplotGroup(Mat.DE.pair.group=Sum.DE.analysis.G$DE.per.pair.G)
        Spe.NoSpe.Barplot<-DEplotBarplot(Cont.per.group, dodge=FALSE)
        Spe.Barplot<-DEplotBarplot(Cont.per.group[-3,], dodge=FALSE)

        List.plots.DE.group[[1]]<-G.Upset$Upset.global
        List.plots.DE.group[[2]]<-Spe.NoSpe.Barplot
        List.plots.DE.group[[3]]<-Spe.Barplot
    }else{
        List.plots.DE.group<-vector(mode="list", length=1)
        names(List.plots.DE.group)<-c("NumberDEgenes_UpDownRegulated")
        #
        Spe.NoSpe.Barplot<-DEplotBarplot(Cont.per.group, dodge=FALSE)
        #
        List.plots.DE.group[[1]]<-Spe.NoSpe.Barplot
    }# if(nb.pair.of.group>1)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 6) Save
    if(!is.null(SubFile.name)){
        SubFile.name<-paste0("_", SubFile.name)
    }# if(is.null(SubFile.name)==TRUE)

    if(nb.pair.of.group==1){
        if(!is.null(path.result)){
            OvUnd2Gfile<-paste0("Plot_NumberDEgenes_UpDownRegulated",
                                SubFile.name, "_OverUnder_DE_2groups", ".pdf")

            grDevices::pdf(file=file.path(path.result, OvUnd2Gfile),
                           width=11, height=8)
            print(Spe.NoSpe.Barplot)
            grDevices::dev.off()
        }else{
            if(Plot.DE.graph==TRUE){
                print(Spe.NoSpe.Barplot)
            }
        }# if(is.null(path.result)==FALSE)
    }# if(nb.pair.of.group==1)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if(nb.pair.of.group>1){
        if(!is.null(path.result)){
            VennFile<-paste0("Plot_VennBarplot",  SubFile.name, ".pdf")

            grDevices::pdf(file=file.path(path.result, VennFile),
                           width=11, height=8)
            print(G.Upset$Upset.global)
            grDevices::dev.off()

            ##----------------------------------------------------------------#
            SpeFile<-paste0("Plot_NumberDEgenes_",
                            "SpecificAndNoSpecific_perBiologicalCondition",
                            SubFile.name, ".pdf")

            grDevices::pdf(file=file.path(path.result, SpeFile),
                           width=11, height=8)
            print(Spe.NoSpe.Barplot)
            grDevices::dev.off()

            ##----------------------------------------------------------------#
            UDregulates<-paste0("Plot_NumberSpecificGenes_",
                                "UpDownRegulated_perBiologicalCondition",
                                SubFile.name, ".pdf")

            grDevices::pdf(file=file.path(path.result, UDregulates),
                           width=11, height=8)
            print(Spe.Barplot)
            grDevices::dev.off()
        }else{
            if(Plot.DE.graph==TRUE){
                print(G.Upset$Upset.global)
                # print(G.Upset$Upset.threshold)
                print(Spe.NoSpe.Barplot)
                print(Spe.Barplot)
            }
        }# if(is.null(path.result)==FALSE)
    }# if(nb.pair.of.group>1)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 6) End
    ## List.Plots.DE.Group=List.plots.DE.group
    return(list(Results=Sum.DE.analysis.G$Results,
                Summary.DEanalysis=Sum.DE.analysis.G$Contingence.per.group))
}## DEanalysisGroup()
