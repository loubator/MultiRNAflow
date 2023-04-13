#' @title Intermediate analysis when samples belong to different
#' biological conditions and different time points.
#'
#' @description This function realizes the intermediate steps of the analysis
#' of the function [DEanalysisTimeAndGroup()].
#'
#' @param DESeq.result Output from the function [DESeq2::DESeq()].
#' @param pval.min Numeric value between 0 and 1. A gene will be considered as
#' differentially expressed (DE) between two biological conditions if its
#' Benjamini-Hochberg adjusted p-value (see [stats::p.adjust()]) is below
#' the threshold \code{pval.min}. Default value is 0.05.
#' @param log.FC.min Non negative numeric value.
#' If the log2 fold change between biological conditions or times has
#' an absolute value below the threshold \code{log.FC.min},
#' then the gene is not selected even if is considered as DE.
#' Default value is 1.
#' If \code{log.FC.min=0}, all DE genes will be kept.
#' @param LRT.supp.info \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the algorithm realizes another statistical test in order
#' to detect if, among all biological conditions and/or times, at least one has
#' a different behavior than the others
#' (see the input \code{test} in [DESeq2::DESeq()]).
#'
#' @importFrom DESeq2 results resultsNames
#'
#' @return The function returns
#' * a data.frame (output 'Global') which contains
#'   * pvalues, log2 fold change and DE genes between each pairs
#'     of biological conditions for a fixed time ti
#'     (except the reference time t0).
#'   * DE specific genes per biological condition for a fixed time ti
#'   (except the reference time t0).
#' * inputs for the functions : [DEplotBarplot()], [DEplotBarplotTime()],
#' [DEplotVennBarplotGroup()], [DEplotVennBarplotTime()],
#' [DEplotBarplotFacetGrid()], [DEplotAlluvial()].
#'
#' @seealso The outputs of the function are used by the main function
#' [DEanalysisTimeAndGroup()].
#'
#' @export
#'
#' @examples
#' data(RawCounts_Schleiss2021_CLLsub500)
#' ## We take only the first three time (both group)
#' ## for the speed of the example
#' Index3t<-c(2:4,11:13,20:22, 29:31,38:40,47:49)
#' RawCounts_P_NP_3t<-RawCounts_Schleiss2021_CLLsub500[,c(1,Index3t)]
#' DESeq2.info<-DEanalysisPreprocessing(RawCounts=RawCounts_P_NP_3t,
#'                                      Column.gene=1,
#'                                      Group.position=2,
#'                                      Time.position=4,
#'                                      Individual.position=3)
#' ##
#' dds.DE <- DESeq2::DESeq(DESeq2.info$DESeq2.obj)
#' ##
#' res.G.T.2<-DEresultGroupPerTime(DESeq.result=dds.DE,
#'                                 LRT.supp.info=FALSE,
#'                                 log.FC.min=1,
#'                                 pval.min=0.05)

DEresultGroupPerTime<-function(DESeq.result,
                               LRT.supp.info=TRUE,
                               pval.min=0.05,
                               log.FC.min=1){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check
    if(class(DESeq.result)[1]!="DESeqDataSet"){
        stop("Res.DE.analysis must a 'DESeqDataSet' object")
    }## if(class(DESeq.result)[1]!="DESeqDataSet")

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 1) Parameters
    Nb.gene<-nrow(DESeq2::results(DESeq.result)) # length(DESeq.result)

    # Vector.group<-as.factor(DESeq.result@colData$Group)
    Vector.group<-as.factor(SummarizedExperiment::colData(DESeq.result)$Group)
    Nb.group<-length(levels(Vector.group))
    nb.pair.of.group<-(Nb.group*(Nb.group-1))/2
    Levels.group<-levels(as.factor(Vector.group))

    # Levels.time<-levels(as.factor(DESeq.result@colData$Time))
    Levels.time<-levels(as.factor(SummarizedExperiment::colData(DESeq.result)$Time))
    Nb.time<-length(Levels.time)
    # time.order=sort(Levels.time)
    # timeline.basis=time.order[1]

    ref.level.group<-Levels.group[1]
    ref.level.time<-Levels.time[1]
    Other.t<-Levels.time[-1]
    timeline.basis.num<-as.numeric(gsub(x=ref.level.time, pattern="t",
                                        replacement=""))
    Other.t.num<-as.numeric(gsub(x=Other.t, pattern="t", replacement=""))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 2) Initialisation data.frame/list all group per time
    # print("DE between biological conditions for each time")
    List.M.sum.DE.pair.g<-vector(mode="list", length=Nb.time)#length(Other.t)
    List.Over.Under.spe.g<-vector(mode="list", length=Nb.time)#length(Other.t)
    List.M.sum.DE.Cont<-vector(mode = "list", length=Nb.time)#length(Other.t)
    names(List.M.sum.DE.Cont)<-Levels.time#Other.t
    List.M.sum.DE.Cont.melt<-vector(mode="list", length=Nb.time)#length(Other.t)
    names(List.M.sum.DE.Cont.melt)<-Levels.time#Other.t
    List.M.DE.pair.g.per.t<-vector(mode="list", length=Nb.time)#length(Other.t)
    List.spe.all.T<-vector(mode="list", length=Nb.time)#length(Other.t)


    if(LRT.supp.info==TRUE){# ind.padj.LRT.sel=which(padj.LRT<pval.min)
        res.LRT<-DESeq2::results(DESeq.result, test="LRT")
        padj.LRT<-res.LRT$padj
        if(length(which(is.na(padj.LRT)))>0){
            padj.LRT[which(is.na(padj.LRT))]<-1
        }# if(length(which(is.na(padj.LRT)))>0)
    }# if(LRT.supp.info==TRUE)


    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 3) Filling data.frames for each time

    for(t in seq_len(Nb.time)){# 1:length(Other.t) # as.numeric(Levels.time)
        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 3.1)  Differential expression between biological conditions

        ## print(paste("DE between biological conditions.",
        ## " Time t",t-1,sep=""))
        Vect.fight.group<-rep(NA, nb.pair.of.group)

        DE.per.2BC<-data.frame(matrix(0, ncol=3*nb.pair.of.group, nrow=Nb.gene))
        Bin.mat.DE<-data.frame(matrix(0, ncol=nb.pair.of.group, nrow=Nb.gene))

        cpt<-0
        gene.DE<-c()

        ##--------------------------------------------------------------------#
        for(i in seq_len(Nb.group-1)){# 1:(Nb.group-1)
            for(k in seq(from=(i+1), to=Nb.group, by=1)){# (i+1):Nb.group
                cpt<-cpt+1
                Index.contrast<-rep(0,times=length(DESeq2::resultsNames(DESeq.result)))
                if(Levels.group[i]==ref.level.group){
                    Ind.givsgk<-grep(pattern=paste0("_", Levels.group[k],
                                                    "_", "vs_",
                                                    Levels.group[i]),
                                     x=DESeq2::resultsNames(DESeq.result))
                    #
                    if(t==0){
                        Index.contrast[Ind.givsgk]<-1
                    }else{
                        Ind.givsgk.t<-grep(pattern=paste0("Time", t-1,
                                                          ".Group",
                                                          Levels.group[k]),
                                           x=DESeq2::resultsNames(DESeq.result))
                        ##t ->t-1
                        Ind.givsgk.t.f<-c(Ind.givsgk,Ind.givsgk.t)
                        Index.contrast[Ind.givsgk.t.f]<-1
                    }

                    fight.group<-paste0(".", levels(Vector.group)[k], "..",
                                        levels(Vector.group)[i], ".")
                }else{
                    Ind.givsgk.p<-grep(pattern=paste0("_", Levels.group[i],
                                                     "_", "vs_",
                                                     ref.level.group),
                                       x=DESeq2::resultsNames(DESeq.result))

                    Ind.givsgk.n<-grep(pattern=paste0("_", Levels.group[k],
                                                     "_", "vs_",
                                                     ref.level.group),
                                       x=DESeq2::resultsNames(DESeq.result))

                    if(t==0){
                        Index.contrast[c(Ind.givsgk.p)]<-1
                        Index.contrast[c(Ind.givsgk.n)]<--1
                    }else{
                        Ind.givsgk.t.p<-grep(pattern=paste0("Time", t-1,
                                                            ".Group",
                                                            Levels.group[i]),
                                             x=DESeq2::resultsNames(DESeq.result))
                        ##t ->t-1
                        Ind.givsgk.t.n<-grep(pattern=paste0("Time", t-1,
                                                            ".Group",
                                                            Levels.group[k]),
                                             x=DESeq2::resultsNames(DESeq.result))
                        ##t ->t-1
                        Index.contrast[c(Ind.givsgk.p,Ind.givsgk.t.p)]<-1
                        Index.contrast[c(Ind.givsgk.n,Ind.givsgk.t.n)]<--1
                    }# if(t==0)

                    fight.group<-paste0(".", levels(Vector.group)[i], "..",
                                        levels(Vector.group)[k], ".")
                }# if(Levels.group[i]==ref.level.group)
                #
                res.dds.group<-DESeq2::results(DESeq.result,
                                               contrast=Index.contrast,
                                               test="Wald")

                # Condition Pvalue<0.05 et |log2Fold-Change|>1
                Padj.i.VS.k<-res.dds.group$padj
                if(length(which(is.na(Padj.i.VS.k)))>0){
                    Padj.i.VS.k[which(is.na(Padj.i.VS.k))]<-1
                }# if(length(which(is.na(Padj.i.VS.k)))>0)

                Log2.FC.i.VS.k<-res.dds.group$log2FoldChange
                if(length(which(is.na(Log2.FC.i.VS.k)))>0){
                    Log2.FC.i.VS.k[which(is.na(Log2.FC.i.VS.k))]<-0
                }# if(length(which(is.na(Log2.FC.i.VS.k)))>0)

                ind.logFC.sel<-which(abs(Log2.FC.i.VS.k)>log.FC.min)
                ind.padj.sel<-which(Padj.i.VS.k<pval.min)

                if(LRT.supp.info==TRUE){
                    ind.padj.LRT.sel<-which(padj.LRT<pval.min)
                    criteria<-sort(intersect(intersect(ind.logFC.sel,ind.padj.sel),
                                             ind.padj.LRT.sel))
                }else{
                    criteria<-sort(intersect(ind.logFC.sel,ind.padj.sel))
                }# if(LRT.supp.info==TRUE)
                #
                gene.DE<-c(gene.DE,criteria)
                pvalue.log2FoldChange<-rep(0,nrow(res.dds.group))
                pvalue.log2FoldChange[criteria]<-1
                #
                Bin.mat.DE[,cpt]<-pvalue.log2FoldChange
                colnames(Bin.mat.DE)[cpt]<-paste0(fight.group, "_t",
                                                  Levels.time[t])#Other.t[t]

                Vect.fight.group[cpt]<-fight.group

                DE.per.2BC[,3*(cpt-1)+c(1,2,3)]<-data.frame(Log2FC=round(Log2.FC.i.VS.k,
                                                                         digits=3),
                                                            Pvalue=Padj.i.VS.k,
                                                            Condition=pvalue.log2FoldChange)
                # round(Padj.i.VS.k,digits=4),
                colnames(DE.per.2BC)[3*(cpt-1)+c(1,2,3)]<-paste0(c("Log2FoldChange.",
                                                                  "Pvalue.adjusted.",
                                                                  "DE."),
                                                                 fight.group)

                data.pair.of.groups<-data.frame(Log2FC=round(Log2.FC.i.VS.k,digits=3),
                                                Pvalue=Padj.i.VS.k,
                                                Condition=pvalue.log2FoldChange)
                # round(Padj.i.VS.k,digits=4),
            }# end for group i
        }# end for group k>i

        ##--------------------------------------------------------------------#
        row.names(Bin.mat.DE)<-row.names(res.dds.group)
        row.names(DE.per.2BC)<-row.names(res.dds.group)

        List.M.DE.pair.g.per.t[[t]]<-Bin.mat.DE
        gene.DE<-sort(unique(gene.DE))
        DE.1tmin.g<-rep(0,times=Nb.gene)
        DE.1tmin.g[DE.1tmin.g]<-1

        if(is.null(row.names(res.dds.group))==TRUE){
            Row.name.res<-paste0("Gene", seq_len(nrow(res.dds.group)))
        }else{# 1:nrow(res.dds.group)
            Row.name.res<-row.names(res.dds.group)
        }## if(is.null(row.names(res.dds.group))==TRUE)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 3.2) Specific genes for each biological condition for each time
        if(nb.pair.of.group>1){
            # print("Specific genes per biological condition")
            Gene.spe.per.group<-matrix(0, nrow=Nb.gene,
                                       ncol=length(levels(Vector.group)))
            colnames(Gene.spe.per.group)<-paste0("Specific.genes_",
                                                levels(Vector.group))
            row.names(Gene.spe.per.group)<-row.names(res.dds.group)

            Nb.DE.per.group<-rep(NA, times=Nb.group)
            names(Nb.DE.per.group)<-levels(Vector.group)

            for(g in seq_len(Nb.group)){# 1:Nb.group
                group.sel<-paste0(".", levels(Vector.group)[g], ".")
                index.group.selec<-grep(pattern=group.sel,
                                        x=Vect.fight.group,
                                        fixed=TRUE)

                Id.Column.spe<-c(seq_len(nb.pair.of.group)*3)[index.group.selec]
                Id.Column.nospe<-c(seq_len(nb.pair.of.group)*3)[-index.group.selec]
                # (1:nb.pair.of.group)

                sum.row.spe<-apply(X=data.frame(DE.per.2BC[,Id.Column.spe]),
                                   MARGIN=1, FUN=sum)
                sum.row.no.spe<-apply(X=data.frame(DE.per.2BC[,Id.Column.nospe]),
                                      MARGIN=1, FUN=sum)

                Nb.DE.per.group[g]<-length(which(sum.row.spe>0))

                Spe.vect<-rep(0,times=nrow(res.dds.group))

                for(i in seq_len(length(Spe.vect))){# 1:length(Spe.vect)
                    if(sum.row.no.spe[i]==0 & sum.row.spe[i]==length(index.group.selec)){
                        Spe.vect[i]<-1
                    }
                }## for(i in 1:length(Spe.vect))
                Gene.spe.per.group[,g]<-Spe.vect
            }## for(g in 1:Nb.group)
        }else{
            Gene.spe.per.group<-cbind(DE.per.2BC[,3], DE.per.2BC[,3])
            colnames(Gene.spe.per.group)<-paste0("Specific.genes_",
                                                 levels(Vector.group))
            row.names(Gene.spe.per.group)<-row.names(res.dds.group)
            Nb.DE.per.group<-rep(length(which(DE.per.2BC[,3]>0)), times=2)
        }## if(Nb.group>1)

        List.spe.all.T[[t]]<-Gene.spe.per.group

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        # 3.3) Over and Under expressed genes in all groups
        OverUnder.expr.per.g<-matrix(0, nrow=nrow(res.dds.group),
                                     ncol=Nb.group)
        colnames(OverUnder.expr.per.g)<-paste0("OverUnder.regulated.genes_Group.",
                                              Levels.group, "_Time.t",
                                              Levels.time[t])
        row.names(OverUnder.expr.per.g)<-Row.name.res #row.names(res.dds.group)

        if(nb.pair.of.group>1){
            for(g in seq_len(Nb.group)){# 1:Nb.group
                group.sel<-paste0(".", levels(Vector.group)[g], ".")
                index.group.sel<-grep(pattern=group.sel,
                                      x=Vect.fight.group,
                                      fixed=TRUE)
                #
                Id.Column.pval<-c(seq_len(nb.pair.of.group)*3)[index.group.sel]
                Id.Column.log2fc<-c(seq_len(nb.pair.of.group)*3-2)[index.group.sel]
                # (1:nb.pair.of.group)
                sign.matrix<-apply(X=DE.per.2BC[,Id.Column.pval]*DE.per.2BC[,Id.Column.log2fc],
                                   MARGIN=2, FUN=function(x) sign(x))

                ##------------------------------------------------------------#
                Position.log0<-matrix(unlist(strsplit(Vect.fight.group[index.group.sel],
                                                      split=".." ,
                                                      fixed=TRUE)),
                                      nrow=2)
                Position.log1<-gsub(".", "", Position.log0, fixed=TRUE)
                Position.log<-matrix(paste0(".", Position.log1, "."),
                                     ncol=length(index.group.sel), byrow=FALSE)
                vec.Position.log<-apply(Position.log, MARGIN=2,
                                        FUN=function(x) which(x==group.sel))
                vec.Position.log[which(vec.Position.log==2)]<--1

                ##------------------------------------------------------------#
                sign.matrix2<-sign.matrix*matrix(rep(vec.Position.log,
                                                     times=nrow(sign.matrix)),
                                                 nrow=nrow(sign.matrix),
                                                 byrow=TRUE)
                sum.sign.matrix<-as.numeric(apply(X=sign.matrix2,
                                                  MARGIN=1,
                                                  FUN=sum))
                sum.sign.matrix.f<-sum.sign.matrix*Gene.spe.per.group[,g]

                ##------------------------------------------------------------#
                for(i in seq_len(nrow(OverUnder.expr.per.g))){
                    sign.gene.i<-sum.sign.matrix.f[i]
                    if(sign.gene.i==length(index.group.sel) & !is.na(sign.gene.i)){
                        OverUnder.expr.per.g[i,g]<-1
                    }
                    if(sign.gene.i==-length(index.group.sel) & !is.na(sign.gene.i)){
                        OverUnder.expr.per.g[i,g]<- -1
                    }
                }# end for(i in 1:nrow(OverUnder.expr.per.g))
            }# end for(g in 1:Nb.group)
        }else{
            OverUnder.expr.per.g[,1]<--sign(pvalue.log2FoldChange*Log2.FC.i.VS.k)
            OverUnder.expr.per.g[,2]<-sign(pvalue.log2FoldChange*Log2.FC.i.VS.k)
        }# if(nb.pair.of.group>1)
        #-------------------------------------------------------------------------#
        List.Over.Under.spe.g[[t]]<-OverUnder.expr.per.g
        contin.spe.g.ini<-rbind(apply(OverUnder.expr.per.g, 2,
                                      function(x) length(which(x==1))),
                                apply(OverUnder.expr.per.g, 2,
                                      function(x) length(which(x==-1))))
        delta.spe.sign.spe<-apply(Gene.spe.per.group,2,
                                  sum) - apply(contin.spe.g.ini,2,sum)
        #
        contin.spe.g<-rbind(rbind(contin.spe.g.ini, delta.spe.sign.spe),
                            Nb.DE.per.group - apply(rbind(contin.spe.g.ini,
                                                          delta.spe.sign.spe)
                                                    , 2, sum))
        colnames(contin.spe.g)<-levels(Vector.group)
        row.names(contin.spe.g)<-c("UpRegulated", "DownRegulated",
                                   "No.specific", "Other")
        #
        if(nb.pair.of.group>1){
            contin.spe.g.f<-data.frame(Attribute=row.names(contin.spe.g),
                                       contin.spe.g)[-3,]
        }else{
            contin.spe.g.f<-data.frame(Attribute=row.names(contin.spe.g),
                                       contin.spe.g)[-c(3,4),]
        }# if(nb.pair.of.group>1)
        #
        List.M.sum.DE.Cont[[t]]<-contin.spe.g.f
        melt.cont<-reshape2::melt(contin.spe.g.f, id.vars="Attribute")
        List.M.sum.DE.Cont.melt[[t]]<-cbind(melt.cont,
                                            Time=rep(Levels.time[t],#Other.t[t]
                                                     times=nrow(melt.cont)))[,c(1,2,4,3)]
        colnames(List.M.sum.DE.Cont.melt[[t]])<-c("Attribute", "Group","Time",
                                                  "value")

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 3.4) Table with all results for one time
        # print("Summary all steps")
        Resultats.DEseq2.groups<-data.frame(DE.1pair.of.Group.minimum=DE.1tmin.g,
                                            cbind(Gene.spe.per.group,
                                                  OverUnder.expr.per.g,
                                                  DE.per.2BC))
        #
        ColNew<-c(seq_len(ncol(Gene.spe.per.group)+1),
                  seq(from=ncol(Gene.spe.per.group)+ncol(OverUnder.expr.per.g)+2,
                      to=ncol(Resultats.DEseq2.groups), by=1))
        #
        colnames(Resultats.DEseq2.groups)[ColNew]<-paste0(colnames(Resultats.DEseq2.groups)[ColNew],
                                                         "_Time.t",
                                                         Levels.time[t])
        row.names(Resultats.DEseq2.groups)<-Row.name.res

        List.M.sum.DE.pair.g[[t]]<-Resultats.DEseq2.groups
    }# for(t in 1:Nb.time)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    # 4) Table with all results for all times
    M.DE.pair.g.per.t<-cbind(Gene=row.names(Resultats.DEseq2.groups),
                             do.call(cbind, List.M.DE.pair.g.per.t))

    Sum.cont.t.g<-do.call(rbind, List.M.sum.DE.Cont.melt)
    colnames(Sum.cont.t.g)<-c("Attribute", "Group", "Time", "value")

    Spe.1t.min<-Reduce('+', List.spe.all.T)
    colnames(Spe.1t.min)<-paste0(colnames(Spe.1t.min), "_1time.minimum")
    if(length(which(Spe.1t.min>1))>0){
        Spe.1t.min[which(Spe.1t.min>1, arr.ind=TRUE)]<-1
    }

    M.sum.DE.pair.g.per.t<-cbind(Gene=row.names(Resultats.DEseq2.groups),
                                 Spe.1t.min,
                                 do.call(cbind, List.M.sum.DE.pair.g))

    M.Over.Under.spe.g<-data.frame(Gene=row.names(Resultats.DEseq2.groups),
                                   do.call(cbind, List.Over.Under.spe.g))
    #
    colnames(M.sum.DE.pair.g.per.t)<-gsub("Log2FoldChange..",
                                          "Log2FoldChange_",
                                          colnames(M.sum.DE.pair.g.per.t),
                                          fixed=TRUE)
    colnames(M.sum.DE.pair.g.per.t)<-gsub("Pvalue.adjusted..",
                                          "Pvalue.adjusted_",
                                          colnames(M.sum.DE.pair.g.per.t),
                                          fixed=TRUE)
    colnames(M.sum.DE.pair.g.per.t)<-gsub("DE..", "DE_",
                                          colnames(M.sum.DE.pair.g.per.t),
                                          fixed=TRUE)
    colnames(M.sum.DE.pair.g.per.t)<-gsub("..", ".versus.",
                                          colnames(M.sum.DE.pair.g.per.t),
                                          fixed=TRUE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    # 5) Output
    return(list(M.sum.DE.pair.G.per.T=M.sum.DE.pair.g.per.t,
                Spe.G.1t.min=Spe.1t.min,
                Sum.cont.G.per.T=Sum.cont.t.g,
                List.DE.per.pair.G.per.T=List.M.DE.pair.g.per.t,
                OverUnder.per.G.per.T=M.Over.Under.spe.g))
}# DEresultGroupPerTime()
