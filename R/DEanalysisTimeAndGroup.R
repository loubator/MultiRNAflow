#' @title DE analysis when samples belong to different biological condition
#' and time points.
#'
#' @description The function realizes from the [DESeq2::DESeq()] output
#' the analysis of :
#' * DE genes between all pairs of biological conditions for each fixed time.
#' * DE genes between all times ti and the reference time t0for each
#' biological condition.
#'
#' @param DESeq.result Output from the function [DESeq2::DESeq()].
#' @param pval.min Numeric value between 0 and 1. A gene will be considered as
#' differentially expressed (DE) between two biological conditions
#' if its Benjamini-Hochberg adjusted p-value (see [stats::p.adjust()])
#' is below the threshold \code{pval.min}. Default value is 0.05
#' @param pval.vect.t \code{NULL} or vector of dimension \eqn{T-1} filled with
#' numeric values between 0 and 1, with \eqn{T} the number of time measurements.
#' A gene will be considered as differentially expressed (DE) between the
#' time ti and the reference time t0 if its Benjamini-Hochberg adjusted p-value
#' (see [stats::p.adjust()]) is below the i-th threshold of \code{pval.vect.t}.
#' If NULL, \code{pval.vect.t} will be vector of dimension \eqn{T-1} filled
#' with \code{pval.min}.
#' @param log.FC.min Non negative numeric value.
#' If the log2 fold change between biological conditions or times has an
#' absolute value below the threshold \code{log.FC.min}, then the gene is
#' not selected even if is considered as DE. Default value is 1.
#' If \code{log.FC.min=0}, all DE genes will be kept.
#' @param LRT.supp.info \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the algorithm realizes another statistical test in order
#' to detect if, among all biological conditions and/or times, at least one has
#' a different behavior than the others
#' (see the input \code{test} in [DESeq2::DESeq()]).
#' @param Plot.DE.graph \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, all graphs will be plotted.
#' Otherwise no graph will be plotted.
#' @param path.result Character or \code{NULL}.
#' If \code{path.result} is a character, it must be a path to a folder,
#' all graphs will be saved in different sub-folders in \code{path.result}.
#' If \code{NULL}, the results will not be saved. \code{NULL} as default.
#' @param SubFile.name Character or \code{NULL}.
#' If \code{SubFile.name} is a character, each saved file names and
#' created folders names will contain the strings of characters
#' "_\code{SubFile.name}". If \code{NULL}, no suffix will be added.
#'
#' @return The function returns
#' * a data.frame (output \code{Results}) which contains
#'   * gene names
#'   * Results from the temporal statistical analysis
#'     * pvalues, log2 fold change and DE genes between each pairs of
#'     biological conditions for each fixed time.
#'     * \eqn{N_{bc}} binary columns (0 and 1), one per biological condition
#'     (with \eqn{N_{bc}} the number of biological conditions).
#'     A 1 in one of these two columns means the gene is DE at least between
#'     one time ti versus the reference time t0, for the biological condition
#'     associated to the given column.
#'     * \eqn{N_{bc}} columns, one per biological condition, where each element
#'     is succession of 0 and 1. The positions of 1 in one of these two columns,
#'     indicate the set of times ti such that the gene is DE between ti and
#'     the reference time t0, for the biological condition associated to
#'     the given column.
#'   * Results from the statistical analysis by biological condition
#'     * pvalues, log2 fold change and DE genes between each time ti
#'     and the reference time t0 for each biological condition.
#'     * \eqn{T} binary columns (0 and 1), one per time
#'     (with \eqn{T} the number of time measurements).
#'     A 1 in one of these columns, means the gene is DE between at least
#'     one pair of biological conditions, for the fixed time associated
#'     to the given column.
#'     * \eqn{T \times N_{bc}} binary columns, which give the genes specific
#'     for each biological condition at each time ti.
#'     A 1 in one of these columns means the gene is specific to the biological
#'     condition at a fixed time associated to the given column. 0 otherwise.
#'     A gene is called specific to a given biological condition BC1 at a
#'     time ti, if the gene is DE between BC1 and any other biological
#'     conditions at time ti, but not DE between any pair of other biological
#'     conditions at time ti.
#'     * \eqn{T \times N_{bc}} columns filled with -1, 0 and 1.
#'     A 1 in one of these columns means the gene is up-regulated
#'     (or over-expressed) for the biological condition at a fixed time
#'     associated to the given column. A gene is called up-regulated for a
#'     given biological condition BC1 at time ti if the gene is specific to
#'     the biological condition BC1 at time ti and expressions in BC1 at time
#'     ti are higher than in the other biological conditions at time ti.
#'     A -1 in one of these columns means the gene is down-regulated
#'     (or under-expressed) for the biological condition at a fixed time
#'     associated to the given column. A gene is called down-regulated for a
#'     given biological condition at a time ti BC1 if the gene is specific to
#'     the biological condition BC1 at time ti and expressions in BC1 at time
#'     ti are lower than in the other biological conditions at time ti.
#'     A 0 in one of these columns means the gene is not specific to the
#'     biological condition at a fixed time associated to the given column.
#'     * \eqn{N_{bc}} binary columns (0 and 1). A 1 in one of these columns
#'     means the gene is specific at at least one time ti, for the biological
#'     condition associated to the given column. 0 otherwise.
#'   * Results from the combination of temporal and biological statistical
#'   analysis
#'     * \eqn{T \times N_{bc}} binary columns, which give the signatures genes
#'     for each biological condition at each time ti. A 1 in one of these
#'     columns means the gene is signature gene to the biological condition at
#'     a fixed time associated to the given column. 0 otherwise.
#'     A gene is called signature of a biological condition BC1 at a given time
#'     ti, if the gene is specific to the biological condition BC1 at time ti
#'     and DE between ti versus the reference time t0 for the biological
#'     condition BC1.
#'     * \eqn{N_{bc}} binary columns (0 and 1). A 1 in one of these columns
#'     means the gene is signature at at least one time ti, for the biological
#'     condition associated to the given column. 0 otherwise.
#' * the following plots from the temporal statistical analysis
#'   * a barplot which gives the number of DE genes between ti and the
#'   reference time t0, for each time ti (except the reference time t0) and
#'   biological condition (see [DEplotBarplotFacetGrid()]).
#'   * \eqn{N_{bc}} alluvial graphs of DE genes (see [DEplotAlluvial()]),
#'   one per biological condition.
#'   * \eqn{N_{bc}} graphs showing the number of DE genes as a function of time
#'   for each temporal group, one per biological condition. By temporal group,
#'   we mean the sets of genes which are first DE at the same time.
#'   * \eqn{2\times N_{bc}} UpSet plot showing the number of DE genes
#'   belonging to each DE temporal pattern, for each biological condition.
#'   By temporal pattern, we mean the set of times ti such that the gene is
#'   DE between ti and the reference time t0 (see [DEplotVennBarplotTime()]).
#'   * an alluvial graph for DE genes which are DE at least one time for
#'   each group.
#' * the following plots from the statistical analysis by biological condition
#'   * a barplot which gives the number of specific DE genes for each
#'   biological condition and time (see [DEplotBarplotFacetGrid()]).
#'   * \eqn{N_{bc}(N_{bc}-1)/2} UpSet plot which give the number of genes
#'   for each possible intersection (set of pairs of biological conditions),
#'   one per time (see [DEplotVennBarplotGroup()]).
#'   * an alluvial graph of genes which are specific at least one time
#'   (see [DEplotAlluvial()]).
#' * the following plots from the combination of temporal and biological
#' statistical analysis
#'   * a barplot which gives the number of signature genes for each biological
#'   condition and time (see [DEplotBarplotFacetGrid()]).
#'   * a barplot showing the number of genes which are DE at at least one time,
#'   specific at at least one time and signature at at least one time,
#'   for each biological condition.
#'   * an alluvial graph of genes which are signature at least one time
#'   (see [DEplotAlluvial()]).
#'
#' @seealso The outputs of the function will be used by the main
#' function [DEanalysisGlobal()].
#'
#' @importFrom DESeq2 results resultsNames
#' @importFrom grDevices dev.off png
#'
#' @export
#'
#' @examples
#' data(RawCounts_Schleiss2021_CLLsub500)
#' ## We take only the first three times (both group) for the speed of
#' ## the example
#' Index3t<-c(2:4,11:13,20:22, 29:31,38:40,47:49)
#' RawCounts_P_NP_3t<-RawCounts_Schleiss2021_CLLsub500[,c(1,Index3t)]
#' DESeq2.info<-DEanalysisPreprocessing(RawCounts=RawCounts_P_NP_3t,
#'                                      Column.gene=1,
#'                                      Group.position=2,
#'                                      Time.position=4,
#'                                      Individual.position=3)
#' ##
#' dds.DE<-DESeq2::DESeq(DESeq2.info$DESeq2.obj)
#' ##
#' res.G.T<-DEanalysisTimeAndGroup(DESeq.result=dds.DE,
#'                                 LRT.supp.info=FALSE,
#'                                 pval.min=0.05,
#'                                 pval.vect.t=NULL,
#'                                 log.FC.min=0.1,
#'                                 Plot.DE.graph=TRUE,
#'                                 path.result=NULL,
#'                                 SubFile.name=NULL)

DEanalysisTimeAndGroup<-function(DESeq.result,
                                 LRT.supp.info=TRUE,
                                 log.FC.min,pval.min,
                                 pval.vect.t,
                                 Plot.DE.graph=TRUE,
                                 path.result,
                                 SubFile.name){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if(class(DESeq.result)[1]!="DESeqDataSet"){
        stop("Res.DE.analysis must a 'DESeqDataSet' object")
    }## if(class(DESeq.result)[1]!="DESeqDataSet")

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Folder creation
    if(!is.null(SubFile.name)){
        SubFile.name<-paste0("_", SubFile.name)
    }# if(is.null(SubFile.name)==TRUE)

    if(!is.null(path.result)){
        Folder_TIMEperBC<-"2-2-1_Temporal_DEanalysis_perGroup"

        if(!Folder_TIMEperBC%in%dir(path=path.result)){
            dir.create(path=file.path(path.result, Folder_TIMEperBC))
            PathResult_TIMEperBC<-file.path(path.result, Folder_TIMEperBC)
        }else{
            PathResult_TIMEperBC<-file.path(path.result, Folder_TIMEperBC)
        }
    }# if(is.null(path.result)==FALSE)

    if(!is.null(path.result)){
        Folder_BCperTime<-paste0("2-2-2_Group_DEanalysis_perTime")#Other.t

        if(!Folder_BCperTime%in%dir(path=path.result)){
            dir.create(path=file.path(path.result, Folder_BCperTime))
            PathResult_BCperTime<-file.path(path.result, Folder_BCperTime)
        }else{
            PathResult_BCperTime<-file.path(path.result, Folder_BCperTime)
        }
    }# if(is.null(path.result)==FALSE)

    if(!is.null(path.result)){
        Folder_TIMEandBC<-"2-2-3_Time_and_Group_DEanalysis"

        if(!Folder_TIMEandBC%in%dir(path=path.result)){
            dir.create(path=file.path(path.result, Folder_TIMEandBC))
            PathResult_TIMEandBC<-file.path(path.result, Folder_TIMEandBC)
        }else{
            PathResult_TIMEandBC<-file.path(path.result, Folder_TIMEandBC)
        }
    }# if(is.null(path.result)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 1) Parameters

    # Gene names and number of genes
    # res.dds<-DESeq.result@rowRanges@partitioning@NAMES
    res.dds<-dimnames(DESeq.result)[[1]]
    # DESeq2::results(DESeq.result)

    if(is.null(res.dds)){
        Row.names.f<-paste0("Gene", seq_len(length(res.dds)))
    }else{
        Row.names.f<-res.dds
    }# if(is.null(res.dds)==TRUE)

    Nb.gene<-length(Row.names.f)
    # nrow(DESeq2::results(DESeq.result))
    # row.names(SummarizedExperiment::assay(DESeq.result))
    # Vector.group<-as.factor(DESeq.result@colData$Group)
    Vector.group<-as.factor(SummarizedExperiment::colData(DESeq.result)$Group)
    Nb.group<-length(levels(Vector.group))
    nb.pair.of.group<-(Nb.group*(Nb.group-1))/2
    Levels.group<-levels(Vector.group)
    # Levels.group<-levels(as.factor(DESeq.result@colData$Group))


    vTimeDESeq<-as.factor(SummarizedExperiment::colData(DESeq.result)$Time)
    vTimeDESeq<-as.factor(gsub("T", "", gsub("t", "", vTimeDESeq)))
    Levels.time<-levels(as.factor(vTimeDESeq))
    Nb.time<-length(Levels.time)
    # levels(DESeq.result@colData$Time)<-levels(as.factor(gsub("T","",
    #                                                          gsub("t","",
    #                                                               DESeq.result@colData$Time))))
    # Levels.time<-levels(as.factor(DESeq.result@colData$Time))
    # Nb.time<-length(Levels.time)


    ref.level.group<-Levels.group[1]
    ref.level.time<-Levels.time[1]
    Id.ref.T<-which(vTimeDESeq==ref.level.time)
    Vect.T.no.ref<-vTimeDESeq[-Id.ref.T]
    Other.t<-Levels.time[-1]
    # DESeq.result@colData$Time

    timeline.basis.num<-as.numeric(gsub(x=ref.level.time, pattern="t",
                                        replacement=""))
    Other.t.num<-as.numeric(gsub(x=Other.t, pattern="t", replacement=""))

    ##------------------------------------------------------------------------#
    ## pvalues
    if(LRT.supp.info==TRUE){
        res.LRT<-DESeq2::results(DESeq.result, test="LRT")
        padj.LRT<-res.LRT$padj
        if(length(which(is.na(padj.LRT)))>0){
            padj.LRT[which(is.na(padj.LRT))]<-1
        }# if(length(which(is.na(padj.LRT)))>0)
    }#if(LRT.supp.info==TRUE)

    if(!is.null(pval.vect.t)){
        if(length(pval.vect.t)>Nb.time-1){
            pval.vect.t<-pval.vect.t[seq_len(Nb.time-1)]# 1:(Nb.time-1)
        }# if(length(pval.vect.t)>Nb.time-1)
        if(length(pval.vect.t)<Nb.time-1){
            pval.vect.t<-c(pval.vect.t,
                           rep(pval.min,
                               times=Nb.time-length(pval.vect.t)-1))
        }# if(length(pval.vect.t)<Nb.time-1)
    }else{
        pval.vect.t<-rep(pval.min, times=Nb.time-1)
    }# if(is.null(pval.vect.t)==FALSE)

    ##------------------------------------------------------------------------#
    # Save plots in R variable
    # m=0 or 1. a(1-m)+e(m)=1+m(e-a).
    # m=min(length(Other.t)-1,1) or min(Nb.group-2,1)
    Nb.Graph.TperG<-2 + min(length(Other.t)-1,1)*((4*Nb.group)+2-2)
    Nb.Graph.GperT<-1 + min(Nb.group-2,1)*(Nb.time+3-1)
    Nb.Graph.GandT<-2 + min(Nb.group-2,1)*(3-2)

    Nb.Graph.all<-Nb.Graph.TperG+Nb.Graph.GperT+Nb.Graph.GandT
    List.plots.DE.time.group<-vector(mode="list", length=Nb.Graph.all)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 2) DE time analysis for each biological condition

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    #### 2.1) Data.frame and list for all time for each biological condition

    List.M.sum.DE.t0<-vector(mode="list", length=Nb.group)
    List.M.time.DE<-vector(mode="list", length=Nb.group)
    List.M.log2FC<-vector(mode="list", length=Nb.group)
    #
    Cont.G.T<-as.data.frame(matrix(0, nrow=Nb.group, ncol=length(Other.t)))
    colnames(Cont.G.T)<-paste0("t_", seq_len(length(Other.t)))
    #
    M.DE.1tmin<-as.data.frame(matrix(0, nrow=Nb.gene, ncol=Nb.group))
    #
    row.names(M.DE.1tmin)<-Row.names.f
    colnames(M.DE.1tmin)<-rep(paste0("DE.1tmin_", Levels.group))

    ##------------------------------------------------------------------------#
    #### 2.2) Data.frame and list filling
    print("DE time analysis for each biological condition.")
    #
    for(g in seq_len(Nb.group)){# 1:Nb.group
        groupe.g<-Levels.group[g]

        M.sum.DE.t0.g<-as.data.frame(matrix(data=0,
                                            nrow=Nb.gene,
                                            ncol=length(Other.t)*3))
        colnames(M.sum.DE.t0.g)<-paste0(rep(c("Log2FoldChange_",
                                              "Pvalue.adjusted_",
                                              "DE_"),
                                            times=length(Other.t)),
                                        rep(paste0("t", Other.t, ".versus.t",
                                                   ref.level.time), each=3),
                                        "_", groupe.g)

        Row.name.assay<-Row.names.f
        row.names(M.sum.DE.t0.g)<-Row.name.assay

        M.time.DE.g<-matrix(data=0, nrow=Nb.gene, ncol=length(Other.t))
        row.names(M.time.DE.g)<-Row.name.assay
        colnames(M.time.DE.g)<-paste0("t", Other.t.num, "_", groupe.g)

        M.log2FC.g<-matrix(data=0, nrow=Nb.gene, ncol=length(Other.t))
        row.names(M.log2FC.g)<-Row.name.assay
        colnames(M.log2FC.g)<-paste0("t", Other.t.num, "_", groupe.g)

        gene.DE.1tmin.g<-c()

        ##------------------------------------------------------------------------#
        #### 2.3) Data.frame filling and DE genes
        for(i in seq_len(length(Other.t))){# 1:length(Other.t)
            Index.contrast<-rep(0,
                                times=length(DESeq2::resultsNames(DESeq.result)))
            Ind.tivst0<-grep(pattern=paste0("Time_", i),
                             x=DESeq2::resultsNames(DESeq.result))

            if(groupe.g!=ref.level.group){
                Ind.tivst0.g<-grep(pattern=paste0("Time", i,
                                                  ".Group", groupe.g),
                                   x=DESeq2::resultsNames(DESeq.result))
                Ind.tivst0.g.f<-c(Ind.tivst0, Ind.tivst0.g)
            }else{
                Ind.tivst0.g.f<-Ind.tivst0
            }# if(groupe.g!=ref.level.group)

            Index.contrast[Ind.tivst0.g.f]<-1
            # print(Index.contrast)

            ##----------------------------------------------------------------#
            res.dds.diff.t<-DESeq2::results(DESeq.result,
                                            contrast=Index.contrast,
                                            test="Wald")

            Log2.FC<-res.dds.diff.t$log2FoldChange
            Log2.FC[which(is.na(Log2.FC))]<-0
            Pvalue.adj<-res.dds.diff.t$padj
            Pvalue.adj[which(is.na(Pvalue.adj))]<-1
            M.sum.DE.t0.g[,(i-1)*3+1]<-round(Log2.FC,digits=3)
            M.sum.DE.t0.g[,(i-1)*3+2]<-Pvalue.adj#round(Pvalue.adj,digits=3)
            M.log2FC.g[,i]<-Log2.FC

            ind.logFC.sel<-which(abs(Log2.FC)>log.FC.min)
            ind.padj.sel<-which(Pvalue.adj<pval.vect.t[i])#which(Pvalue.adj<pval.min)

            if(LRT.supp.info==TRUE){
                ind.padj.LRT.sel<-which(padj.LRT<pval.min)
                gene.DE.i<-sort(intersect(intersect(ind.logFC.sel,ind.padj.sel),
                                          ind.padj.LRT.sel))
            }else{
                gene.DE.i<-sort(intersect(ind.logFC.sel, ind.padj.sel))
            }# if(LRT.supp.info==TRUE)

            gene.DE.1tmin.g<-c(gene.DE.1tmin.g, gene.DE.i)
            M.time.DE.g[gene.DE.i,i]<-1
            M.sum.DE.t0.g[gene.DE.i,(i-1)*3+3]<-1
            Cont.G.T[g,i]<-length(gene.DE.i)
            row.names(Cont.G.T)[g]<-groupe.g
        }# for(i in 1:length(Other.t))

        ##--------------------------------------------------------------------#
        #### 2.4) Results
        M.DE.1tmin[sort(unique(gene.DE.1tmin.g)), g]<-1
        #
        Pattern.tps.g<-apply(M.time.DE.g,
                             MARGIN=1,
                             FUN=function(x) paste(as.character(x),
                                                   collapse=""))
        #
        M.sum.DE.t0.g.f<-data.frame(DE.1t.min.g=M.DE.1tmin[,g],
                                    Pattern.DE.g=paste0(".",
                                                       Pattern.tps.g, "."),
                                    M.sum.DE.t0.g)
        colnames(M.sum.DE.t0.g.f)<-c(paste0(c("DE.1time.minimum_",
                                             "Pattern.DE_"), groupe.g),
                                     colnames(M.sum.DE.t0.g))

        List.M.sum.DE.t0[[g]]<-M.sum.DE.t0.g.f
        List.M.time.DE[[g]]<-M.time.DE.g
        List.M.log2FC[[g]]<-M.log2FC.g

        ##--------------------------------------------------------------------#
        ## 2.5) Graphs
        # Alluvial graph et nb de gene per group per time
        Col.M.time.DE.g<-colnames(M.time.DE.g)
        colnames(M.time.DE.g)<-paste0("t", seq_len(ncol(M.time.DE.g)))

        if(length(Other.t)>1){
            title.alluvial.TG<-paste0("Alluvial graph of genes which are DE at",
                                     "at least one time,\nfor the group ",
                                     groupe.g)
            title.evolution.TG<-paste0("Time evolution of the number of genes which ",
                                      "are DE at at least one time within each ",
                                      "temporal group,\nfor the group ",
                                      groupe.g)

            res.allu.g<-DEplotAlluvial(table.DE.time=M.time.DE.g,
                                       Temporal.Group=TRUE,
                                       title.alluvial=title.alluvial.TG,
                                       title.evolution=title.evolution.TG)
            colnames(M.time.DE.g)<-Col.M.time.DE.g

            ##----------------------------------------------------------------#
            # Upset Venn graph with or not over vs t0
            res.VennBarplot<-DEplotVennBarplotTime(table.DE.time=M.time.DE.g,
                                                   Log2.FC.matrix=M.log2FC.g)

            ##----------------------------------------------------------------#
            List.plots.DE.time.group[[2+(g-1)*4+1]]<-res.allu.g$g.alluvial
            List.plots.DE.time.group[[2+(g-1)*4+2]]<-res.allu.g$g.alluvial.freq
            List.plots.DE.time.group[[2+(g-1)*4+3]]<-res.VennBarplot$Upset.graph
            List.plots.DE.time.group[[2+(g-1)*4+4]]<-res.VennBarplot$Upset.graph.with.nb.over

            names(List.plots.DE.time.group)[2+(g-1)*4+1]<-paste0("Alluvial.graph.Group_",
                                                                 groupe.g)
            names(List.plots.DE.time.group)[2+(g-1)*4+2]<-paste0("NumberDEgenes.acrossTime",
                                                                 ".perTemporalGroup",
                                                                 ".Group_",
                                                                 groupe.g)
            names(List.plots.DE.time.group)[2+(g-1)*4+3]<-paste0("VennBarplot.Group_",
                                                                 groupe.g)
            names(List.plots.DE.time.group)[2+(g-1)*4+4]<-paste0("VennBarplot.",
                                                                 "withNumberUpRegulated",
                                                                 ".Group_",
                                                                 groupe.g)
        }# if(length(Other.t)>1)

        ##--------------------------------------------------------------------#
        ### 2.6) Save
        if(!is.null(path.result)){
            if(length(Other.t)>1){
                # Save alluvial graph

                AlluGfile<-paste0("Group_", groupe.g, "_AlluvialGraph",
                                  SubFile.name, ".pdf")

                grDevices::pdf(file=file.path(PathResult_TIMEperBC, AlluGfile),
                               width=11, height=8)
                print(res.allu.g$g.alluvial)
                grDevices::dev.off()

                # Time evolution per temporal group
                NbDEgfile<-paste0("Group_", groupe.g,
                                  "_NumberDEgenes_acrossTime_perTemporalGroup",
                                  SubFile.name, ".pdf")

                grDevices::pdf(file=file.path(PathResult_TIMEperBC, NbDEgfile),
                               width=11, height=8)
                print(res.allu.g$g.alluvial.freq)
                grDevices::dev.off()

                # Upset graph
                Venngfile<-paste0("Group_", groupe.g, "_VennBarplot",
                                  SubFile.name, ".pdf")

                grDevices::pdf(file=file.path(PathResult_TIMEperBC, Venngfile),
                               width=11, height=8)
                print(res.VennBarplot$Upset.graph)
                grDevices::dev.off()

                VenngUPfile<-paste0("Group_", groupe.g,
                                    "_VennBarplot_withNumberUpRegulated",
                                    SubFile.name, ".pdf")

                grDevices::pdf(file=file.path(PathResult_TIMEperBC,
                                              VenngUPfile),
                               width=11, height=8)
                print(res.VennBarplot$Upset.graph.with.nb.over)
                grDevices::dev.off()
            }# if(length(Other.t)>1)
        }else{
            if(length(Other.t)>1 & Plot.DE.graph==TRUE){
                print(res.allu.g$g.alluvial)
                print(res.allu.g$g.alluvial.freq)
                print(res.VennBarplot$Upset.graph)
                print(res.VennBarplot$Upset.graph.with.nb.over)
            }# if(length(Other.t)>1)
        }# if(is.null(path.result)==FALSE)
    }# for(g in 1:Nb.group)

    ##------------------------------------------------------------------------#
    ### 2.7) Alluvial of DE genes at least one time for each group

    title.alluvial.DEg<-paste("Alluvial graph of genes which are DE at",
                              "at least one time for each group", sep=" ")

    Col.M.DE.1tmin<-colnames(M.DE.1tmin)
    colnames(M.DE.1tmin)<-gsub("DE.1tmin_","",Col.M.DE.1tmin)

    res.allu.g1tmin<-DEplotAlluvial(table.DE.time=M.DE.1tmin,
                                    Temporal.Group=FALSE,
                                    title.alluvial=title.alluvial.DEg)
    colnames(M.DE.1tmin)<-Col.M.DE.1tmin

    ##------------------------------------------------------------------------#
    ### 2.8) Summary of DE time analysis per group data
    Mat.facet<-data.frame(Attribute=rep(c("UpRegulated", "DownRegulated"),
                                        times=(Nb.time-1)*Nb.group),
                          Group=rep(Levels.group, each=(Nb.time-1)*2),
                          Time=rep(rep(paste0("t", seq_len(Nb.time-1)),
                                       each=2),
                                   times=Nb.group),# 1:(Nb.time-1)
                          value=rep(0, times=2*(Nb.time-1)*Nb.group))

    for(g in seq_len(length(List.M.time.DE))){# 1:length(List.M.time.DE)
        prodlogpval<-List.M.time.DE[[g]]*sign(List.M.log2FC[[g]])
        nb.up<-apply(prodlogpval,2,function(x) length(which(x==1)))
        nb.down<-apply(prodlogpval,2,function(x) length(which(x==-1)))
        id.Matfacet2<-seq_len(2*(Nb.time-1)) + (g-1)*(2*(Nb.time-1))
        #
        UpDown<-rep(nb.up,each=2)
        UpDown[2*seq_len(Nb.time-1)]<-nb.down# (1:(Nb.time-1))
        #
        Mat.facet[id.Matfacet2,4]<-UpDown
    }# for(g in seq_len(length(List.M.time.DE)))

    Color.Leg.T<-data.frame(Levels=c("UpRegulated", "DownRegulated"),
                            Color=c("#EA8626", "#8BEAF0"))
    #c("#E41A1C", "steelblue")

    res.facet.BC<-DEplotBarplotFacetGrid(Data=Mat.facet,
                                         Abs.col=3,
                                         Legend.col=1,
                                         Facet.col=2,
                                         Value.col=4,
                                         Color.Legend=Color.Leg.T)

    ##------------------------------------------------------------------------#
    List.plots.DE.time.group[[1]]<-res.facet.BC
    List.plots.DE.time.group[[2]]<-res.allu.g1tmin
    names(List.plots.DE.time.group)[1]<-"AlluvialGraph_DE.1tmin_perGroup"
    names(List.plots.DE.time.group)[2]<-paste0("NumberDEgenes_UpDownRegulated",
                                               "_perTimeperGroup")

    ##------------------------------------------------------------------------#
    ### 2.9) Summary of DE time analysis per group Graph
    if(!is.null(path.result)){
        Allu1tmin<-paste0("TimeResultsPerGroup_",
                          "AlluvialGraph_DE1Tmin_perGroup",
                          SubFile.name, ".pdf")

        grDevices::pdf(file=file.path(PathResult_TIMEperBC, Allu1tmin),
                       width=11, height=8)
        print(res.allu.g1tmin)
        grDevices::dev.off()

        NbUDperGT<-paste0("TimeResultsPerGroup_",
                          "NumberDEgenes_UpDownRegulated_perTimeperGroup",
                          SubFile.name, ".pdf")

        grDevices::pdf(file=file.path(PathResult_TIMEperBC, NbUDperGT),
                       width=11, height=8)
        print(res.facet.BC)
        grDevices::dev.off()
    }else{
        if(Plot.DE.graph==TRUE){
            print(res.allu.g1tmin)
            print(res.facet.BC)
        }
    }# if(is.null(path.result)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ### 3) DE group analysis for each time
    print("DE group analysis for each time measurement.")

    res.sum.G.T<-DEresultGroupPerTime(DESeq.result=DESeq.result,
                                      LRT.supp.info=LRT.supp.info,
                                      log.FC.min=log.FC.min,
                                      pval.min=pval.min)

    ##------------------------------------------------------------------------#
    ### 3.1) Specific and no specific DE genes per time and
    ###      Specific genes at least one time for each group

    ## Data summary DE group analysis per time
    Cont.Group.per.Time<-res.sum.G.T$Sum.cont.G.per.T
    Cont.Group.per.Time$Time<-paste0("t",
                                     gsub("t", "",
                                          x=Cont.Group.per.Time$Time,
                                          fixed=TRUE))
    ## Color
    if(Nb.group>2){
        Color.Leg.Spe<-data.frame(Levels=c("Other",
                                           "UpRegulated",
                                           "DownRegulated"),
                                  Color=c("#999999", "#E41A1C", "steelblue"))
    }else{
        Color.Leg.Spe<-data.frame(Levels=c("UpRegulated", "DownRegulated"),
                                  Color=c("#E41A1C", "steelblue"))
    }# if(Nb.group>2)

    # Graph
    res.dodge<-DEplotBarplotFacetGrid(Data=Cont.Group.per.Time,
                                      Abs.col=2,
                                      Legend.col=1,
                                      Facet.col=3,
                                      Value.col=4,
                                      Color.Legend=Color.Leg.Spe)

    if(Nb.group>2){
        Ind.cont.other<-which(Cont.Group.per.Time$Attribute=="Other")
        Sum.cont.G.per.T.spe.sign<-Cont.Group.per.Time[-Ind.cont.other,]
        res.dodge.spe.sign<-DEplotBarplotFacetGrid(Data=Sum.cont.G.per.T.spe.sign,
                                                   Abs.col=2,
                                                   Legend.col=1,
                                                   Facet.col=3,
                                                   Value.col=4,
                                                   Color.Legend=Color.Leg.Spe[-1,])

        ##--------------------------------------------------------------------#
        # Specific genes at least one time for each group
        SpeDat1tmin<-data.frame(res.sum.G.T$Spe.G.1t.min)
        colnames(SpeDat1tmin)<-Levels.group

        TitleSpe.g<-paste("Alluvial graph of specific genes at at least one time",
                          "for each group", sep=" ")

        res.allu.Spe1tmin<-DEplotAlluvial(SpeDat1tmin, FALSE, TitleSpe.g)

        ##--------------------------------------------------------------------#
        List.plots.DE.time.group[[Nb.Graph.TperG+1]]<-res.dodge
        List.plots.DE.time.group[[Nb.Graph.TperG+2]]<-res.dodge.spe.sign
        List.plots.DE.time.group[[Nb.Graph.TperG+3]]<-res.allu.Spe1tmin
        names(List.plots.DE.time.group)[Nb.Graph.TperG+1]<-paste0("NumberDEgenes_",
                                                                  "SpecificAndNoSpecific_",
                                                                  "perBiologicalCondition")
        names(List.plots.DE.time.group)[Nb.Graph.TperG+2]<-paste0("NumberSpecificGenes_",
                                                                  "UpDownRegulated_",
                                                                  "perBiologicalCondition")
        names(List.plots.DE.time.group)[Nb.Graph.TperG+3]<-paste0("AlluvialGraph_",
                                                                  "SpecificGenes1tmin_",
                                                                  "perBiologicalCondition")
    }else{
        List.plots.DE.time.group[[Nb.Graph.TperG+1]]<-res.dodge
        names(List.plots.DE.time.group)[Nb.Graph.TperG+1]<-paste0("NumberSpecificGenes_",
                                                                  "UpDownRegulated_",
                                                                  "perBiologicalCondition")
    }# if(Nb.group>2)

    ##------------------------------------------------------------------------#
    ### 3.2) Save Venn barplot group analysis for each time and previous graph

    # List.DE.per.pair.G.per.T, List.DE.per.pair.G.per.T,
    # List.DE.per.pair.G.per.T
    if(Nb.group>2){
        for(t in seq_len(Nb.time)){# 1:Nb.time
            G.Upset.t<-DEplotVennBarplotGroup(res.sum.G.T$List.DE.per.pair.G.per.T[[t]])

            Timet<-paste0("t",
                          gsub(x=gsub(x=Levels.time[t], pattern="T",
                                      replacement=""),
                               pattern="t", replacement=""))

            List.plots.DE.time.group[[Nb.Graph.TperG+3+t]]<-G.Upset.t$Upset.global
            names(List.plots.DE.time.group)[Nb.Graph.TperG+3+t]<-paste0("VennBarplot_",
                                                                        "BiologicalConditions",
                                                                        "_atTime.",
                                                                        Timet)

            # Other.t.num<-as.numeric(gsub(x=Other.t, pattern="t", replacement=""))
            if(!is.null(path.result)){
                VennTfile<-paste0("VennBarplot_BiologicalConditions_atTime.",
                                  Timet, SubFile.name, ".pdf")

                grDevices::pdf(file=file.path(PathResult_BCperTime, VennTfile),
                               width=11, height=8)
                print(G.Upset.t$Upset.global)
                grDevices::dev.off()
            }else{
                if(Plot.DE.graph==TRUE){
                    print(G.Upset.t$Upset.global)
                }
            }# if(is.null(path.result)==FALSE)
        }# for(t in 1:Nb.time)
    }# if(Nb.group>2)

    if(!is.null(path.result)){
        if(Nb.group>2){
            SpeCondfile<-paste0("GroupResultsPerTime_",
                                "NumberDEgenes_SpecificAndNoSpecific_",
                                "perBiologicalCondition",
                                SubFile.name, ".pdf")

            grDevices::pdf(file=file.path(PathResult_BCperTime, SpeCondfile),
                           width=11, height=8)
            print(res.dodge)
            grDevices::dev.off()

            SpeCondfile2<-paste0("GroupResultsPerTime_",
                                 "NumberSpecificGenes_UpDownRegulated_",
                                 "perBiologicalCondition",
                                 SubFile.name, ".pdf")

            grDevices::pdf(file=file.path(PathResult_BCperTime, SpeCondfile2),
                           width=11, height=8)
            print(res.dodge.spe.sign)
            grDevices::dev.off()


            AlluSpeG<-paste0("GroupResultsPerTime_",
                             "_AlluvialGraph_SpecificGenes1tmin_perBiologicalCondition",
                             SubFile.name, ".pdf")

            grDevices::pdf(file=file.path(PathResult_BCperTime, AlluSpeG),
                           width=11, height=8)
            print(res.allu.Spe1tmin)
            grDevices::dev.off()
        }else{
            SpeCondfile3<-paste0("GroupResultsPerTime_",
                                 "NumberSpecificGenes_UpDownRegulated_",
                                 "perBiologicalCondition",
                                 SubFile.name, ".pdf")

            grDevices::pdf(file=file.path(PathResult_BCperTime, SpeCondfile3),
                           width=11, height=8)
            print(res.dodge)
            grDevices::dev.off()
        }# if(Nb.group>2)
    }else{
        if(Plot.DE.graph==TRUE){
            print(res.dodge)
            if(Nb.group>2){
                print(res.dodge.spe.sign)
                print(res.allu.Spe1tmin)
            }# if(Nb.group>2)
        }
    }# if(is.null(path.result)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    #### 4) DE group and time analysis
    print("Combined time and group results.")

    ##------------------------------------------------------------------------#
    ##### 4.0) Preprocessing

    M.sum.1tmin.all.g<-do.call(cbind, List.M.sum.DE.t0)
    DE1timeAllg<-do.call(cbind, List.M.time.DE)

    DEspeAllti<-res.sum.G.T$OverUnder.per.G.per.T
    DEspeAllti.not0<-DEspeAllti[,-grep("t0", colnames(DEspeAllti), fixed=TRUE)]

    DEtxt0<-do.call(cbind, List.M.time.DE)

    Intersect.G.all.t<-M.DE.1tmin*res.sum.G.T$Spe.G.1t.min

    if(length(which(Intersect.G.all.t>1))>1){
        Intersect.G.all.t[which(Intersect.G.all.t>1, arr.ind=TRUE)]<-1
    }# if(length(which(Intersect.G.all.t>1))>1)
    colnames(Intersect.G.all.t)<-paste0("Intersect_", Levels.group, "_T")

    ##------------------------------------------------------------------------#
    #### 4.1) Signature Data
    Signature.Mat.Sign<-data.frame(matrix(0,
                                          nrow=Nb.gene,
                                          ncol=Nb.group*Nb.time))

    for(g in seq_len(Nb.group)){# 1:Nb.group
        Group.g<-Levels.group[g]
        # for(t in seq_len(Nb.time-1)){# 1:(Nb.time-1)
        #   Time.t<-paste("t",Levels.time[t+1],sep="")
        #   t.sel<-grep(paste(Time.t,"_",sep=""),
        #               paste(colnames(DE1timeAllg),"_",sep=""),
        #               fixed=TRUE)
        # }# for(t in 1:(Nb.time-1))
        g.sel.Spe<-grep(paste0(".", Group.g, "_"),
                        colnames(DEspeAllti.not0), fixed=TRUE)
        g.sel.Time<-grep(paste0("_", Group.g),
                         colnames(DE1timeAllg), fixed=TRUE)

        Signature.g.t<-DE1timeAllg[,g.sel.Time]*DEspeAllti.not0[,g.sel.Spe]
        Sum.Signature.g.t<-apply(abs(Signature.g.t), 1, sum)

        if(max(Sum.Signature.g.t)>1){
            Sum.Signature.g.t[Sum.Signature.g.t>1]<-1
        }# if(max(Sum.Signature.g.t)>1)

        ## 1:(Nb.time-1)
        IdColSignature<-Nb.group+(Nb.time-1)*(g-1)+seq_len(Nb.time-1)
        Signature.Mat.Sign[,IdColSignature]<-Signature.g.t
        colnames(Signature.Mat.Sign)[IdColSignature]<-paste0("Signature.genes_",
                                                             "Group.",
                                                             Group.g, "_Time.t",
                                                             seq_len(Nb.time-1))
        Signature.Mat.Sign[,g]<-Sum.Signature.g.t
        colnames(Signature.Mat.Sign)[g]<-paste0("Signature.genes_Group.",
                                                Group.g,
                                                "_1time.minimum")
    }# for(g in 1:Nb.group)
    ## do.call(cbind,
    ## List.M.time.DE)#abs(res.sum.G.T$OverUnder.per.G.per.T[,-1])

    ##------------------------------------------------------------------------#
    # 4.2) Data containing the summary DE group and time analysis
    Sum.final.g.t<-cbind(Gene=res.sum.G.T$M.sum.DE.pair.G.per.T[,1],
                         abs(Signature.Mat.Sign), #Intersect.G.all.t,
                         M.sum.1tmin.all.g,
                         res.sum.G.T$M.sum.DE.pair.G.per.T[,-1])

    ##------------------------------------------------------------------------#
    # 4.3) Alluvial graph of signature genes at least one time for each group
    if(Nb.group>2){
        alluvialSignature<-matrix(0, ncol=Nb.group,
                                  nrow=nrow(abs(Signature.Mat.Sign)))
        Signature.All.TG<-abs(Signature.Mat.Sign)[,-seq_len(Nb.group)]

        for(g in seq_len(Nb.group)){
            Index.g.signature<-(g-1)*(Nb.time-1)+seq_len(Nb.time-1)
            alluvialSignature[,g]<-apply(Signature.All.TG[,Index.g.signature],1,sum)
        }# for(g in 1:Nb.group)

        alluvialSignature[which(alluvialSignature>1, arr.ind=TRUE)]<-1
        colnames(alluvialSignature)<-Levels.group

        title.allu.signature<-paste0("Alluvial graph of signature genes at ",
                                     "at least one time for each group")
        res.allu.signature<-DEplotAlluvial(alluvialSignature,
                                           FALSE,
                                           title.allu.signature)
    }# if(Nb.group>2)

    ##------------------------------------------------------------------------#
    ### 4.3) Number of DE genes per time for each group with
    ###      signature information

    Mat.facet.S<-data.frame(Attribute=rep(c("UpRegulated", "DownRegulated"),
                                          times=(Nb.time-1)*Nb.group),
                            Group=rep(rep(Levels.group, each=2),
                                      times=(Nb.time-1)),
                            Time=rep(paste0("t", seq_len(Nb.time-1)),
                                     each=2*Nb.group),
                            value=rep(0, times=2*(Nb.time-1)*Nb.group))
    #
    for(t in seq_len(Nb.time-1)){
        t.sel.Sign<-grep(paste0("_Time.t",Levels.time[t+1]),
                         colnames(Signature.Mat.Sign),
                         fixed=TRUE)
        t.sel.txt0<-grep(paste0("t", Levels.time[t+1], "_"),
                         colnames(DEtxt0), fixed=TRUE)

        prodlogpval.S<-Signature.Mat.Sign[,t.sel.Sign]*DEtxt0[,t.sel.txt0]
        nb.up.S<-apply(prodlogpval.S, 2, function(x) length(which(x==1)))
        nb.down.S<-apply(prodlogpval.S, 2, function(x) length(which(x==-1)))
        id.Matfacet2.S<-seq_len(2*Nb.group) + (t-1)*(2*Nb.group)

        UpDown.S<-as.numeric(rep(nb.up.S, each=2))
        UpDown.S[2*seq_len(Nb.group)]<-as.numeric(nb.down.S)

        Mat.facet.S[id.Matfacet2.S,4]<-as.numeric(UpDown.S)
    }# for(g in 1:(Nb.time-1))

    # Data for graph
    Mat.facet.all<-rbind(Mat.facet, Mat.facet.S)
    Mat.facet.all$Attribute<-as.factor(c(Mat.facet$Attribute,
                                         paste0("Signature",
                                               Mat.facet.S$Attribute)))

    Dfacet<-Mat.facet.S[order(Mat.facet.S$Group),]
    Mat.facet.all$value[seq_len(nrow(Mat.facet))]<-abs(Mat.facet$value-Dfacet$value)

    Color.Leg.S<-data.frame(Levels=c("SignatureUpRegulated", "UpRegulated",
                                     "DownRegulated",
                                     "SignatureDownRegulated"),
                            Color=c("#E41A1C", "#EA8626",
                                    "#8BEAF0","steelblue"))

    res.facet.BC.signature<-DEplotBarplotFacetGrid(Data=Mat.facet.all,
                                                   Abs.col=3,
                                                   Legend.col=1,
                                                   Facet.col=2,
                                                   Value.col=4,
                                                   Color.Legend=Color.Leg.S)

    ##------------------------------------------------------------------------#
    #### 4.4) Summary DE, specific and signature at least one time
    Mat.Sum<-data.frame(Attribute=rep(c("DE.1time.minimum",
                                        "Specific.1time.minimum",
                                        "Signature.1time.minimum"),
                                      times=Nb.group),
                        Group=rep(Levels.group, each=3),
                        value=rep(0, times=3*Nb.group))

    Mat.Sum$Attribute<-factor(Mat.Sum$Attribute,
                              levels=c("DE.1time.minimum",
                                       "Specific.1time.minimum",
                                       "Signature.1time.minimum"))

    Nb.DE1tminperG<-apply(M.DE.1tmin, 2, sum)
    Nb.Spe.1tmin<-apply(res.sum.G.T$Spe.G.1t.min, 2, sum)
    Nb.Signature.1tmin<-apply(abs(Signature.Mat.Sign[,seq_len(Nb.group)]),
                              2,
                              sum)

    #apply(Intersect.G.all.t, 2, sum)
    for(g in seq_len(Nb.group)){# 1:Nb.group
        Mat.Sum$value[(g-1)*3+c(1, 2, 3)]<-c(Nb.DE1tminperG[g],
                                             Nb.Spe.1tmin[g],
                                             Nb.Signature.1tmin[g])
    }# for(g in seq_len(Nb.group))

    res.Sum.signature<-DEplotBarplotFacetGrid(Data=Mat.Sum,
                                              Abs.col=1,
                                              Legend.col=1,
                                              Facet.col=2,
                                              Value.col=3,
                                              Color.Legend=NULL)

    ##------------------------------------------------------------------------#
    Nb.GandT<-Nb.Graph.TperG+Nb.Graph.GperT
    List.plots.DE.time.group[[Nb.GandT+1]]<-res.facet.BC.signature
    List.plots.DE.time.group[[Nb.GandT+2]]<-res.Sum.signature
    names(List.plots.DE.time.group)[Nb.GandT+1]<-paste0("Number_DEgenes_SignatureGenes_",
                                                       "UpDownRegulated_perTimeperGroup")
    names(List.plots.DE.time.group)[Nb.GandT+2]<-paste0("Number_DEgenes1TimeMinimum_",
                                                       "Specific1TimeMinimum_",
                                                       "Signature1TimeMinimum_",
                                                       "perBiologicalCondition")
    if(Nb.group>2){
        List.plots.DE.time.group[[Nb.GandT+3]]<-res.Sum.signature
        names(List.plots.DE.time.group)[Nb.GandT+3]<-paste0("Alluvial_SignatureGenes_",
                                                           "1TimeMinimum_perGroup")
    }# if(Nb.group>2)

    ##------------------------------------------------------------------------#
    ##### 4.5) Save graph
    if(!is.null(path.result)){
        SignaFile<-paste0("Number_DEgenes_SignatureGenes_",
                          "UpDownRegulated_perTimeperGroup", SubFile.name,
                          ".pdf")

        grDevices::pdf(file=file.path(PathResult_TIMEandBC, SignaFile),
                       width=11, height=8)
        print(res.facet.BC.signature)
        grDevices::dev.off()

        Sum1tminFile<-paste0("Number_DEgenes1TimeMinimum_",
                             "Specific1TimeMinimum_",
                             "Signature1TimeMinimum_",
                             "perBiologicalCondition", SubFile.name, ".pdf")

        grDevices::pdf(file=file.path(PathResult_TIMEandBC, Sum1tminFile),
                       width=11, height=8)
        print(res.Sum.signature)
        grDevices::dev.off()

        if(Nb.group>2){
            AlluSignaFile<-paste0("Alluvial_SignatureGenes_",
                                  "1TimeMinimum_perGroup",
                                  SubFile.name, ".pdf")

            grDevices::pdf(file=file.path(PathResult_TIMEandBC, AlluSignaFile),
                           width=11, height=8)
            print(res.allu.signature)
            grDevices::dev.off()
        }# if(Nb.group>2)
        #
    }else{
        if(Plot.DE.graph==TRUE){
            print(res.facet.BC.signature)
            print(res.Sum.signature)
            if(Nb.group>2){
                print(res.allu.signature)
            }# if(Nb.group>2)
        }
    }# if(is.null(path.result)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ##### 5) Output
    # List.Plots.DE.Time.Group=List.plots.DE.time.group
    return(list(Results=Sum.final.g.t,
                Summary.DEanalysis=Mat.facet.all))
}## DEanalysisTimeAndGroup()
