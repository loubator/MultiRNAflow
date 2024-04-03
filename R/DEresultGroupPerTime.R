#' @title Intermediate analysis when samples belong to different
#' biological conditions and different time points.
#'
#' @description This function realizes the intermediate steps of the analysis
#' of the function [DEanalysisTimeAndGroup()].
#'
#' @param DESeq.result Output from the function
#' [DESeq2::DESeq()].
#' @param pval.min Numeric value between 0 and 1. A gene will be considered as
#' differentially expressed (DE) between two biological conditions if its
#' Benjamini-Hochberg adjusted p-value
#' (see [stats::p.adjust()])
#' is below the threshold \code{pval.min}. Default value is 0.05.
#' @param log.FC.min Non negative numeric value.
#' If the log2 fold change between biological conditions or times has
#' an absolute value below the threshold \code{log.FC.min},
#' then the gene is not selected even if is considered as DE.
#' Default value is 1.
#' If \code{log.FC.min=0}, all DE genes will be kept.
#' @param LRT.supp.info \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the algorithm realizes another statistical test in order
#' to detect if, among all biological conditions and/or times, at least one has
#' a different behavior than the others (see the input \code{test} in
#' [DESeq2::DESeq()]).
#'
#' @importFrom DESeq2 results resultsNames
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom reshape2 melt
#'
#' @return The function returns the same DESeqDataSet class object
#' \code{DESeq.result} with the following results,
#' saved in the metadata \code{DEresultsTimeGroup} of \code{DESeq.result}:
#' * a data.frame (output \code{DEsummary} of \code{DEresultsTimeGroup})
#' which contains
#'   * pvalues, log2 fold change and DE genes between each pairs
#'     of biological conditions for a fixed time ti
#'     (except the reference time t0).
#'   * DE specific genes per biological condition for a fixed time ti
#'   (except the reference time t0).
#' * inputs for the functions :
#' [DEplotBarplot()],
#' [DEplotBarplotTime()],
#' [DEplotVennBarplotGroup()],
#' [DEplotVennBarplotTime()],
#' [DEplotBarplotFacetGrid()],
#' [DEplotAlluvial()].
#'
#' @seealso The output of the function are used by the main function
#' [DEanalysisTimeAndGroup()].
#'
#' @export
#'
#' @examples
#' data("RawCounts_Schleiss2021_CLLsub500")
#' ## We take only the first three times (both group) for the speed of
#' ## the example
#' Index3t<-c(2:4,11:13,20:22, 29:31,38:40,47:49)
#' RawCounts_3t<-RawCounts_Schleiss2021_CLLsub500[seq_len(200), c(1,Index3t)]
#'
#' ## Preprocessing step
#' resDATAprepSEleuk <- DATAprepSE(RawCounts=RawCounts_3t,
#'                                 Column.gene=1,
#'                                 Group.position=2,
#'                                 Time.position=4,
#'                                 Individual.position=3)
#'
#' DESeq2preprocess <- S4Vectors::metadata(resDATAprepSEleuk)$DESeq2obj
#' DESeq2obj <- DESeq2preprocess$DESeq2preproceesing
#'
#' ##------------------------------------------------------------------------##
#' dds.DE <- DESeq2::DESeq(DESeq2obj)
#' ##
#' res.G.T.2 <- DEresultGroupPerTime(DESeq.result=dds.DE,
#'                                   LRT.supp.info=FALSE,
#'                                   log.FC.min=1,
#'                                   pval.min=0.05)

DEresultGroupPerTime <- function(DESeq.result,
                                 LRT.supp.info=TRUE,
                                 pval.min=0.05,
                                 log.FC.min=1){
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check
    if(!is(DESeq.result, 'DESeqDataSet')){
        stop("Res.DE.analysis must be a 'DESeqDataSet' object")
    }## if(!is(classDeseq2, 'DESeqDataSet'))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 1) Parameters
    Nb.gene<-nrow(DESeq2::results(DESeq.result)) ## length(DESeq.result)
    namesCoefDESeq2 <- DESeq2::resultsNames(DESeq.result)


    Vector.group <- SummarizedExperiment::colData(DESeq.result)$Group
    Vector.group <- as.factor(Vector.group)

    Levels.group <- levels(Vector.group)
    ref.level.group <- Levels.group[1]
    Nb.group <- length(Levels.group)
    NpairGroups <- (Nb.group*(Nb.group-1))/2


    Vector.time <- SummarizedExperiment::colData(DESeq.result)$Time
    Vector.time <- as.factor(Vector.time)

    Levels.time <- levels(Vector.time)
    ref.level.time <- Levels.time[1]
    Other.t <- Levels.time[-1]
    Nb.time <- length(Levels.time)
    ## time.order=sort(Levels.time) ## timeline.basis=time.order[1]

    timeline.basis.num <- as.numeric(gsub(x=ref.level.time, pattern="t",
                                          replacement=""))
    Other.t.num <- as.numeric(gsub(x=Other.t, pattern="t", replacement=""))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 2) Initialisation data.frame/list all group per time
    ## print("DE between biological conditions for each time")
    ## length=Nb.time ## before: length=length(Other.t)
    List.M.sum.DE.pair.g <- vector(mode="list", length=Nb.time)
    List.Over.Under.spe.g <- vector(mode="list", length=Nb.time)
    List.M.sum.DE.Cont <- vector(mode = "list", length=Nb.time)
    names(List.M.sum.DE.Cont) <- Levels.time ##Other.t
    List.M.sum.DE.Cont.melt <- vector(mode="list", length=Nb.time)
    names(List.M.sum.DE.Cont.melt)<-Levels.time ##Other.t
    listDEpairGperT <- vector(mode="list", length=Nb.time)
    List.spe.all.T <- vector(mode="list", length=Nb.time)

    if (isTRUE(LRT.supp.info)) {## ind.padj.LRT.sel=which(padj.LRT<pval.min)
        res.LRT <- DESeq2::results(DESeq.result, test="LRT")
        padj.LRT <- res.LRT$padj
        if (length(which(is.na(padj.LRT))) > 0) {
            padj.LRT[which(is.na(padj.LRT))] <- 1
        }# if(length(which(is.na(padj.LRT)))>0)
    }# if(isTRUE(LRT.supp.info))


    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 3) Filling data.frames for each time

    for(t in seq_len(Nb.time)){
        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 3.1)  Differential expression between biological conditions

        ## print(paste("DE between biological conditions.",
        ## " Time t",t-1,sep=""))
        Vect.fight.group <- rep(NA, NpairGroups)

        DEper2BCt <- data.frame(matrix(0, ncol=3*NpairGroups, nrow=Nb.gene))
        Bin.mat.DE <- data.frame(matrix(0, ncol=NpairGroups, nrow=Nb.gene))

        cpt <- 0
        gene.DE <- c()

        ##-------------------------------------------------------------------##
        for (i in seq_len(Nb.group-1)) {# 1:(Nb.group-1)
            for (k in seq(from=(i+1), to=Nb.group, by=1)) {
                cpt <- cpt+1
                Index.contrast <- rep(0, times=length(namesCoefDESeq2))

                if (Levels.group[i] == ref.level.group) {
                    Ind.givsgk <- grep(pattern=paste0("_", Levels.group[k],
                                                      "_", "vs_",
                                                      Levels.group[i]),
                                       x=namesCoefDESeq2)

                    if (t == 0) {
                        Index.contrast[Ind.givsgk] <- 1
                    } else {
                        Ind.givsgk.t <- grep(pattern=paste0("Time", t-1,
                                                            ".Group",
                                                            Levels.group[k]),
                                             x=namesCoefDESeq2)
                        ##t ->t-1
                        Ind.givsgk.t.f <- c(Ind.givsgk,Ind.givsgk.t)
                        Index.contrast[Ind.givsgk.t.f] <- 1
                    }## if (t == 0)

                    fight.group <- paste0(".", levels(Vector.group)[k],
                                          "..", levels(Vector.group)[i], ".")
                } else {
                    Ind.givsgk.p <- grep(pattern=paste0("_", Levels.group[i],
                                                        "_", "vs_",
                                                        ref.level.group),
                                         x=namesCoefDESeq2)

                    Ind.givsgk.n <- grep(pattern=paste0("_", Levels.group[k],
                                                        "_", "vs_",
                                                        ref.level.group),
                                         x=namesCoefDESeq2)

                    if (t == 0) {
                        Index.contrast[c(Ind.givsgk.p)] <- 1
                        Index.contrast[c(Ind.givsgk.n)] <- -1
                    } else {
                        Ind.givsgk.t.p <- grep(pattern=paste0("Time", t-1,
                                                              ".Group",
                                                              Levels.group[i]),
                                               x=namesCoefDESeq2)
                        ##t ->t-1
                        Ind.givsgk.t.n <- grep(pattern=paste0("Time", t-1,
                                                              ".Group",
                                                              Levels.group[k]),
                                               x=namesCoefDESeq2)
                        ##t ->t-1
                        Index.contrast[c(Ind.givsgk.p,Ind.givsgk.t.p)] <- 1
                        Index.contrast[c(Ind.givsgk.n,Ind.givsgk.t.n)] <- -1
                    }## if(t==0)

                    fight.group <- paste0(".", levels(Vector.group)[i],
                                          "..", levels(Vector.group)[k], ".")
                }## if(Levels.group[i]==ref.level.group)

                res.dds.group <- DESeq2::results(DESeq.result,
                                                 contrast=Index.contrast,
                                                 test="Wald")

                ##-----------------------------------------------------------##
                ## Condition Pvalue<0.05 et |log2Fold-Change|>1
                Padj.i.VS.k <- res.dds.group$padj
                if (length(which(is.na(Padj.i.VS.k))) > 0) {
                    Padj.i.VS.k[which(is.na(Padj.i.VS.k))] <- 1
                }## if(length(which(is.na(Padj.i.VS.k)))>0)

                Log2FC.i.VS.k <- res.dds.group$log2FoldChange
                if (length(which(is.na(Log2FC.i.VS.k))) > 0) {
                    Log2FC.i.VS.k[which(is.na(Log2FC.i.VS.k))]<-0
                }## if(length(which(is.na(Log2FC.i.VS.k)))>0)

                ind.logFC.sel <- which(abs(Log2FC.i.VS.k)>log.FC.min)
                ind.padj.sel <- which(Padj.i.VS.k<pval.min)

                ##-----------------------------------------------------------##
                criteria <- sort(intersect(ind.logFC.sel, ind.padj.sel))

                if (isTRUE(LRT.supp.info)) {
                    ind.padj.LRT.sel <- which(padj.LRT<pval.min)
                    criteria <- sort(intersect(criteria, ind.padj.LRT.sel))
                }## if(LRT.supp.info==TRUE)

                ##-----------------------------------------------------------##
                gene.DE <- c(gene.DE, criteria)
                pval_log2FC <- rep(0, nrow(res.dds.group))
                pval_log2FC[criteria] <- 1

                Bin.mat.DE[,cpt] <- pval_log2FC
                colnames(Bin.mat.DE)[cpt] <- paste0(fight.group, "_t",
                                                    Levels.time[t])#Other.t[t]

                Vect.fight.group[cpt] <- fight.group

                IDcol2BCt <- 3*(cpt - 1) + c(1, 2, 3)
                DEper2BCt[,IDcol2BCt] <- data.frame(Log2FC=round(Log2FC.i.VS.k,
                                                                 digits=3),
                                                    Pvalue=Padj.i.VS.k,
                                                    Condition=pval_log2FC)
                ## round(Padj.i.VS.k,digits=4),
                colnames(DEper2BCt)[IDcol2BCt] <- paste0(c("Log2FoldChange.",
                                                           "Pvalue.adjusted.",
                                                           "DE."),
                                                         fight.group)

                data.pair.of.groups <- data.frame(Log2FC=round(Log2FC.i.VS.k,
                                                               digits=3),
                                                  Pvalue=Padj.i.VS.k,
                                                  Condition=pval_log2FC)
                ## round(Padj.i.VS.k,digits=4),
            }## end for group i
        }## end for group k>i

        ##-------------------------------------------------------------------##
        row.names(Bin.mat.DE) <- row.names(res.dds.group)
        row.names(DEper2BCt) <- row.names(res.dds.group)

        listDEpairGperT[[t]] <- Bin.mat.DE
        gene.DE <- sort(unique(gene.DE))
        DE.1tmin.g <- rep(0,times=Nb.gene)
        DE.1tmin.g[DE.1tmin.g] <- 1

        if (is.null(row.names(res.dds.group)) == TRUE) {
            Row.name.res <- paste0("Gene", seq_len(nrow(res.dds.group)))
        } else {
            Row.name.res <- row.names(res.dds.group)
        }## if(is.null(row.names(res.dds.group))==TRUE)

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 3.2) Specific genes for each biological condition for each time
        if (NpairGroups>1) {
            ## print("Specific genes per biological condition")
            SPEgenesPerGroup <- matrix(0, nrow=Nb.gene,
                                         ncol=length(levels(Vector.group)))
            colnames(SPEgenesPerGroup) <- paste0("Specific.genes_",
                                                   levels(Vector.group))
            row.names(SPEgenesPerGroup) <- row.names(res.dds.group)

            Nb.DE.per.group <- rep(NA, times=Nb.group)
            names(Nb.DE.per.group) <- levels(Vector.group)

            for (g in seq_len(Nb.group)) {
                group.sel <- paste0(".", levels(Vector.group)[g], ".")
                IDgroupSelect <- grep(pattern=group.sel,
                                      x=Vect.fight.group,
                                      fixed=TRUE)

                IDcolSPE <- c(seq_len(NpairGroups)*3)[IDgroupSelect]
                IDcolNOspe <- c(seq_len(NpairGroups)*3)[-IDgroupSelect]

                SUMrowSPE <- apply(X=data.frame(DEper2BCt[,IDcolSPE]),
                                   MARGIN=1, FUN=sum)
                SUMrowNOspe <- apply(X=data.frame(DEper2BCt[,IDcolNOspe]),
                                     MARGIN=1, FUN=sum)

                Nb.DE.per.group[g] <- length(which(SUMrowSPE>0))

                Spe.vect <- rep(0,times=nrow(res.dds.group))

                for (i in seq_len(length(Spe.vect))) {
                    if(SUMrowNOspe[i]==0 & SUMrowSPE[i]==length(IDgroupSelect)){
                        Spe.vect[i] <- 1
                    }
                }## for(i in 1:length(Spe.vect))

                SPEgenesPerGroup[,g] <- Spe.vect

            }## for(g in 1:Nb.group)
        } else {
            SPEgenesPerGroup <- cbind(DEper2BCt[,3], DEper2BCt[,3])
            colnames(SPEgenesPerGroup) <- paste0("Specific.genes_",
                                                 levels(Vector.group))
            row.names(SPEgenesPerGroup) <- row.names(res.dds.group)
            Nb.DE.per.group <- rep(length(which(DEper2BCt[,3]>0)), times=2)
        }## if(Nb.group>1)

        List.spe.all.T[[t]] <- SPEgenesPerGroup

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        # 3.3) Over and Under expressed genes in all groups
        OverUnder.expr.per.g <- matrix(0, nrow=nrow(res.dds.group),
                                       ncol=Nb.group)
        colnames(OverUnder.expr.per.g) <- paste0("OverUnder.regulated.genes",
                                                 "_Group.", Levels.group,
                                                 "_Time.t", Levels.time[t])
        row.names(OverUnder.expr.per.g) <- Row.name.res
        ## row.names(res.dds.group)

        if (NpairGroups>1) {
            for (g in seq_len(Nb.group)) {
                group.sel <- paste0(".", levels(Vector.group)[g], ".")
                IDgroupSelect <- grep(pattern=group.sel,
                                      x=Vect.fight.group,
                                      fixed=TRUE)

                IDcolPVAL <- c(seq_len(NpairGroups)*3)[IDgroupSelect]
                IDcolLog2FC <- c(seq_len(NpairGroups)*3-2)[IDgroupSelect]

                DEper2BCt_DE <- DEper2BCt[,IDcolPVAL]*DEper2BCt[,IDcolLog2FC]
                sign.matrix <- apply(X=DEper2BCt_DE,
                                     MARGIN=2,
                                     FUN=function(x) sign(x))

                ##-----------------------------------------------------------##
                Position.log0_str <- strsplit(Vect.fight.group[IDgroupSelect],
                                              split=".." ,
                                              fixed=TRUE)
                Position.log0 <- matrix(unlist(Position.log0_str), nrow=2)
                Position.log1 <- gsub(".", "", Position.log0, fixed=TRUE)

                Position.log <- matrix(paste0(".", Position.log1, "."),
                                       ncol=length(IDgroupSelect), byrow=FALSE)
                vec.Position.log <- apply(Position.log, MARGIN=2,
                                          FUN=function(x) which(x==group.sel))
                vec.Position.log[which(vec.Position.log == 2)] <- -1

                ##-----------------------------------------------------------##
                sign.matrix2 <- sign.matrix*matrix(rep(vec.Position.log,
                                                       times=nrow(sign.matrix)),
                                                   nrow=nrow(sign.matrix),
                                                   byrow=TRUE)
                sum.sign.matrix <- as.numeric(apply(X=sign.matrix2,
                                                    MARGIN=1,
                                                    FUN=sum))
                sum.sign.matrix.f <- sum.sign.matrix*SPEgenesPerGroup[,g]

                ##-----------------------------------------------------------##
                for (i in seq_len(nrow(OverUnder.expr.per.g))) {
                    sign.gene.i <- sum.sign.matrix.f[i]
                    NidGroupSel <- length(IDgroupSelect)
                    if(sign.gene.i == NidGroupSel & !is.na(sign.gene.i)){
                        OverUnder.expr.per.g[i, g] <- 1
                    }
                    if(sign.gene.i == -NidGroupSel & !is.na(sign.gene.i)){
                        OverUnder.expr.per.g[i, g] <- -1
                    }
                }## for (i in seq_len(nrow(OverUnder.expr.per.g)))
            }## for(g in seq_len(Nb.group))
        }else{
            OverUnder.expr.per.g[,1] <- -sign(pval_log2FC*Log2FC.i.VS.k)
            OverUnder.expr.per.g[,2] <- sign(pval_log2FC*Log2FC.i.VS.k)
        }## if(NpairGroups>1)
        ##-------------------------------------------------------------------##
        List.Over.Under.spe.g[[t]] <- OverUnder.expr.per.g
        contin.spe.g.ini <- rbind(apply(OverUnder.expr.per.g, 2,
                                        function(x) length(which(x == 1))),
                                  apply(OverUnder.expr.per.g, 2,
                                        function(x) length(which(x == -1))))
        delta.spe.sign.spe <- apply(SPEgenesPerGroup,2,
                                    sum) - apply(contin.spe.g.ini, 2, sum)

        contin.spe.g <- rbind(rbind(contin.spe.g.ini, delta.spe.sign.spe),
                              Nb.DE.per.group - apply(rbind(contin.spe.g.ini,
                                                            delta.spe.sign.spe)
                                                      , 2, sum))
        colnames(contin.spe.g) <- levels(Vector.group)
        row.names(contin.spe.g) <- c("UpRegulated", "DownRegulated",
                                     "No.specific", "Other")

        if (NpairGroups>1) {
            contin.spe.g.f <- data.frame(Attribute=row.names(contin.spe.g),
                                         contin.spe.g)[-3,]
        } else {
            contin.spe.g.f <- data.frame(Attribute=row.names(contin.spe.g),
                                         contin.spe.g)[-c(3, 4),]
        }## if(NpairGroups>1)

        List.M.sum.DE.Cont[[t]] <- contin.spe.g.f
        melt.cont <- reshape2::melt(contin.spe.g.f, id.vars="Attribute")

        meltContTime <- rep(Levels.time[t], times=nrow(melt.cont)) ##Other.t[t]
        List.M.sum.DE.Cont.melt[[t]] <- cbind(melt.cont,
                                              Time=meltContTime)[, c(1,2,4,3)]
        colnames(List.M.sum.DE.Cont.melt[[t]]) <- c("Attribute", "Group",
                                                    "Time", "value")

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 3.4) Table with all results for one time
        ## print("Summary all steps")
        ResultsDEseq2groups <- data.frame(DE.1pair.of.Group.minimum=DE.1tmin.g,
                                          cbind(SPEgenesPerGroup,
                                                OverUnder.expr.per.g,
                                                DEper2BCt))

        SUMcolSpOvUn <- ncol(SPEgenesPerGroup) + ncol(OverUnder.expr.per.g) + 2
        ColNew <- c(seq_len(ncol(SPEgenesPerGroup)+1),
                    seq(from=SUMcolSpOvUn, to=ncol(ResultsDEseq2groups), by=1))
        ColnamesSelect <- colnames(ResultsDEseq2groups)[ColNew]

        colnames(ResultsDEseq2groups)[ColNew] <- paste0(ColnamesSelect,
                                                        "_Time.t",
                                                        Levels.time[t])
        row.names(ResultsDEseq2groups) <- Row.name.res

        List.M.sum.DE.pair.g[[t]] <- ResultsDEseq2groups
    }## for(t in 1:Nb.time)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 4) Table with all results for all times
    M.DE.pair.g.per.t <- cbind(Gene=row.names(ResultsDEseq2groups),
                               do.call(cbind, listDEpairGperT))

    Sum.cont.t.g <- do.call(rbind, List.M.sum.DE.Cont.melt)
    colnames(Sum.cont.t.g) <- c("Attribute", "Group", "Time", "value")

    Spe.1t.min <- Reduce('+', List.spe.all.T)
    colnames(Spe.1t.min) <- paste0(colnames(Spe.1t.min), "_1time.minimum")
    if (length(which(Spe.1t.min>1)) > 0) {
        Spe.1t.min[which(Spe.1t.min>1, arr.ind=TRUE)] <- 1
    }## if (length(which(Spe.1t.min>1)) > 0)

    summaryDEpairGperT <- cbind(Gene=row.names(ResultsDEseq2groups),
                                Spe.1t.min,
                                do.call(cbind, List.M.sum.DE.pair.g))

    M.Over.Under.spe.g <- data.frame(Gene=row.names(ResultsDEseq2groups),
                                     do.call(cbind, List.Over.Under.spe.g))

    colnames(summaryDEpairGperT) <- gsub("Log2FoldChange..",
                                         "Log2FoldChange_",
                                         colnames(summaryDEpairGperT),
                                         fixed=TRUE)
    colnames(summaryDEpairGperT) <- gsub("Pvalue.adjusted..",
                                         "Pvalue.adjusted_",
                                         colnames(summaryDEpairGperT),
                                         fixed=TRUE)
    colnames(summaryDEpairGperT) <- gsub("DE..", "DE_",
                                         colnames(summaryDEpairGperT),
                                         fixed=TRUE)
    colnames(summaryDEpairGperT) <- gsub("..", ".versus.",
                                         colnames(summaryDEpairGperT),
                                         fixed=TRUE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## SE object
    listDEresGT <- list(summaryDEpairGperT=summaryDEpairGperT,
                        Spe.G.1t.min=Spe.1t.min,
                        Sum.cont.G.per.T=Sum.cont.t.g,
                        listDEpairGperT=listDEpairGperT,
                        OverUnder.per.G.per.T=M.Over.Under.spe.g)

    DESeqclass <- DESeq.result
    S4Vectors::metadata(DESeqclass)$DEresultsTimeGroup <- listDEresGT

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 5) Output
    return(DESeqclass)
}# DEresultGroupPerTime()
