#' @title DE analysis when samples belong to different time points only.
#'
#' @description The function realizes from the
#' [DESeq2::DESeq()]
#' output the analysis of DE genes between each time versus
#' the reference time t0.
#'
#' @param DESeq.result Output from the function
#' [DESeq2::DESeq()].
#' @param pval.min Numeric value between 0 and 1. A gene will be considered as
#' differentially expressed (DE) between two biological conditions if
#' its Benjamini-Hochberg adjusted p-value
#' (see [stats::p.adjust()])
#' is below the threshold pval.min. Default value is 0.05
#' @param pval.vect.t \code{NULL} or vector of dimension \eqn{T-1} filled with
#' numeric values between 0 and 1, with \eqn{T} the number of time
#' measurements.
#' A gene will be considered as differentially expressed (DE) between the time
#' ti and the reference time t0 if its Benjamini-Hochberg adjusted p-value
#' (see [stats::p.adjust()])
#' is below the i-th threshold of \code{pval.vect.t}.
#' If NULL, \code{pval.vect.t} will be vector of dimension \eqn{T-1} filled
#' with \code{pval.min}.
#' @param log.FC.min Non negative numeric value.
#' If the \eqn{log_2} fold change between biological conditions or times has
#' an absolute value below the threshold \code{log.FC.min}, then the gene is
#' not selected even if is considered as DE. Default value is 1.
#' If \code{log.FC.min=0}, all DE genes will be kept.
#' @param LRT.supp.info \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the algorithm realizes another statistical test in order to
#' detect if, among all biological conditions and/or times, at least one has
#' a different behavior than the others (see the input 'test' in
#' [DESeq2::DESeq()]).
#' @param Plot.DE.graph \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, all graphs will be plotted.
#' Otherwise no graph will be plotted.
#' @param path.result Character or \code{NULL}.
#' If \code{path.result} is a character, it must be a path to a folder,
#' all graphs will be saved in \code{path.result}.
#' If \code{NULL}, the results will not be saved in a folder. NULL as default.
#' @param SubFile.name Character or \code{NULL}
#' If \code{SubFile.name} is a character, each saved file names will contain
#' the strings of characters "_\code{SubFile.name}".
#' If \code{NULL}, no suffix will be added.
#'
#' @importFrom DESeq2 results
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom ggplot2 ggplot aes xlab ylab geom_bar scale_fill_manual
#' theme_minimal
#' @importFrom grDevices pdf dev.off
#'
#' @return The function returns the same DESeqDataSet class object
#' \code{DESeq.result} with the following results,
#' saved in the metadata \code{DEresultsTime} of \code{DESeq.result}:
#' * a data.frame (output \code{DEsummary} of \code{DEresultsTime})
#' which contains
#'   * gene names
#'   * pvalues, log2 fold change and DE genes between each time ti versus the
#'   reference time t0.
#'   * a binary column (1 and 0) where 1 means the gene is DE at least at
#'   between one time ti versus the reference time t0.
#'   * a column where each element is succession of 0 and 1.
#'   The positions of '1' indicate the set of times ti such that the gene is
#'   DE between ti and the reference time t0.
#' * an alluvial graph of differentially expressed (DE) genes
#' (see [DEplotAlluvial()])
#' * a graph showing the number of DE genes as a function of time for each
#' temporal group
#' (see [DEplotAlluvial()]).
#' By temporal group, we mean the sets of genes which are first DE at
#' the same time.
#' * a barplot which gives the number of DE genes per time
#' (see [DEplotBarplotTime()])
#' * an UpSet plot (Venn diagram displayed as a barplot) which gives the number
#' of genes per temporal pattern
#' (see [DEplotVennBarplotTime()]).
#' By temporal pattern, we mean the set of times ti such that the gene is
#' DE between ti and the reference time t0.
#' * a similar UpSet plot where each bar is split in different colors
#' corresponding to all possible numbers of DE times where genes are over
#' expressed in a given temporal pattern.
#'
#' @seealso The outputs of the function will be used by the main function
#' [DEanalysisGlobal()].
#'
#' @export
#'
#' @examples
#' data(RawCounts_Leong2014_FISSIONsub500wt)
#' ## We take only the first three time for the speed of the example
#' RawCounts_Fission_3t<-RawCounts_Leong2014_FISSIONsub500wt[seq_len(200),
#'                                                           seq_len(10)]
#'
#' ## Preprocessing step
#' resDATAprepSEfission <- DATAprepSE(RawCounts=RawCounts_Fission_3t,
#'                                    Column.gene=1,
#'                                    Group.position=NULL,
#'                                    Time.position=2,
#'                                    Individual.position=3)
#'
#' DESeq2preprocess <- S4Vectors::metadata(resDATAprepSEfission)$DESeq2obj
#' DESeq2obj <- DESeq2preprocess$DESeq2preproceesing
#'
#' ##------------------------------------------------------------------------##
#' dds.DE.T <- DESeq2::DESeq(DESeq2obj, quiet=TRUE, betaPrior=FALSE)
#' ##
#' res.T <- DEanalysisTime(DESeq.result=dds.DE.T,
#'                         pval.min=0.05,
#'                         pval.vect.t=c(0.01,0.05,0.05),
#'                         log.FC.min=1,
#'                         LRT.supp.info=FALSE,
#'                         Plot.DE.graph=TRUE,
#'                         path.result=NULL,
#'                         SubFile.name=NULL)

DEanalysisTime <- function(DESeq.result,
                           pval.min=0.05,
                           pval.vect.t=NULL,
                           log.FC.min=1,
                           LRT.supp.info=FALSE,
                           Plot.DE.graph=TRUE,
                           path.result=NULL,
                           SubFile.name=NULL){
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (!is(DESeq.result, 'DESeqDataSet')) {
        stop("Res.DE.analysis must be a 'DESeqDataSet' object")
    }## if(!is(classDeseq2, 'DESeqDataSet'))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 1) Parameters
    ## Gene names and number of genes

    Gene.Names <- dimnames(DESeq.result)[[1]] ## DESeq2::results(DESeq.result)

    if (is.null(Gene.Names)) {
        Row.name.res <- paste0("Gene", seq_len(length(Gene.Names)))
    } else {
        Row.name.res <- Gene.Names
    }## if(is.null(row.names(res.dds.group))==TRUE)

    Nb.gene <- length(Row.name.res)

    ## Time
    Fct.time <- data.frame(SummarizedExperiment::colData(DESeq.result))[,1]
    Fct.time <- as.factor(as.character(Fct.time))
    Levels.time <- levels(Fct.time)
    ref.level.time <- Levels.time[1]
    Other.t <- Levels.time[-1]
    Nb.time <- length(Levels.time)
    Other.t.num <- as.numeric(gsub("T", "", gsub("t", "", Other.t)))
    ref.level.time.num <- as.numeric(gsub("T", "",
                                          gsub("t", "", ref.level.time)))

    ## pval.vect.t adjustment if necessary
    if (!is.null(pval.vect.t)) {
        if (length(pval.vect.t)>Nb.time-1) {
            pval.vect.t <- pval.vect.t[seq_len(Nb.time-1)]
        }## if(length(pval.vect.t)>Nb.time-1)

        if (length(pval.vect.t)<Nb.time-1) {
            pval.vect.t <- c(pval.vect.t,
                           rep(pval.min, times=Nb.time-length(pval.vect.t)-1))
        }## if(length(pval.vect.t)<Nb.time-1)
    } else {
        pval.vect.t <- rep(pval.min, times=Nb.time-1)
    }## if(is.null(pval.vect.t)==FALSE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 2) Data.frames initialisation
    ## Table with all results
    Table.Results <- as.data.frame(matrix(data=0, nrow=Nb.gene,
                                          ncol=(Nb.time-1)*3))
    colnames(Table.Results) <- paste0(rep(c("Log2FoldChange.",
                                            "Pvalue.adjusted.", "DE."),
                                          times=Nb.time-1),
                                      rep(paste0("t", Other.t.num, ".versus.t",
                                                 ref.level.time.num),
                                          each=3))
    row.names(Table.Results) <- Row.name.res

    ## Table which indicates the times where each gene is DE
    table.DE.Times <- matrix(data=0, nrow=Nb.gene, ncol=Nb.time-1)
    colnames(table.DE.Times) <- paste0("t", Other.t.num)
    row.names(table.DE.Times) <- Row.name.res

    ## Table which contains all log2 fold change
    table.log2FC <- matrix(data=0, nrow=Nb.gene, ncol=Nb.time-1)
    colnames(table.log2FC) <- paste0("t", Other.t.num)
    row.names(table.log2FC) <- Row.name.res

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 3) Filling of data.frames and DE genes
    if (isTRUE(LRT.supp.info)) {
        res.LRT <- DESeq2::results(DESeq.result, test="LRT")
        ## the argument 'reduced' is useless because there is only factor time

        padj.LRT <- res.LRT$padj

        if (length(which(is.na(padj.LRT)))>0) {
            padj.LRT[which(is.na(padj.LRT))] <- 1
        }## if(length(which(is.na(padj.LRT)))>0)
    }## if(isTRUE(LRT.supp.info))

    DE.gene.all.ti <- c()

    for (i in seq_len(length(Other.t))) {
        res.dds.norm.diff.i <- DESeq2::results(DESeq.result, test="Wald",
                                               contrast=c("Time", Other.t[i],
                                                          ref.level.time))

        Pvalue.adj <- res.dds.norm.diff.i$padj

        if (length(which(is.na(Pvalue.adj)))>0) {
            Pvalue.adj[which(is.na(Pvalue.adj))] <- 1
        }## if(length(which(is.na(Pvalue.adj)))>0)

        Log2.FC <- res.dds.norm.diff.i$log2FoldChange

        if (length(which(is.na(Log2.FC)))>0) {
            Log2.FC[which(is.na(Log2.FC))] <- 0
        }## if(length(which(is.na(Log2.FC)))>0)

        ind.logFC.sel <- which(abs(Log2.FC)>log.FC.min)
        ind.padj.sel <- which(Pvalue.adj<pval.vect.t[i])

        if (isTRUE(LRT.supp.info)) {
            ind.padj.LRT.sel <- which(padj.LRT<pval.min)
            DE.gene.ti <- sort(intersect(intersect(ind.logFC.sel,
                                                   ind.padj.sel),
                                         ind.padj.LRT.sel))
        } else {
            DE.gene.ti <- sort(intersect(ind.logFC.sel, ind.padj.sel))
        }## if (LRT.supp.info == TRUE)

        Table.Results[,(i-1)*3+1] <- round(Log2.FC,digits=3)
        Table.Results[,(i-1)*3+2] <- Pvalue.adj ## round(Pvalue.adj,digits=3)
        table.log2FC[, i] <- Log2.FC
        DE.gene.all.ti <- c(DE.gene.all.ti, DE.gene.ti)
        table.DE.Times[DE.gene.ti, i] <- 1
        Table.Results[DE.gene.ti, (i-1)*3+3] <- 1
    }## end for(i in 1:length(Other.t))

    DE.at.least.1t <- rep(0, times=nrow(Table.Results))
    DE.at.least.1t[sort(unique(DE.gene.all.ti))] <- 1

    Pattern.tps <- apply(table.DE.Times,
                         MARGIN=1,
                         FUN=function(x) paste(as.character(x), collapse=""))

    Table.Results.f <- data.frame(Gene=Row.name.res,
                                  DE.1time.minimum=DE.at.least.1t,
                                  DE.Temporal.Pattern=paste0(".", Pattern.tps,
                                                             "."),
                                  Table.Results)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 4) Graphs
    Nb.Graph <- (1 - min(Nb.time, 1)) + 5*min(Nb.time, 1)
    listDEplotsTIME <- vector(mode="list", length=Nb.Graph)

    ##-----------------------------------------------------------------------##
    if (Nb.time == 2) {
        Pos.2t <- length(which(table.DE.Times*table.log2FC>0))
        Neg.2t <- length(which(table.DE.Times*table.log2FC<0))

        ## To avoid " no visible binding for global variable"
        ## in devtools::check()
        Number.DE <- Attribute <- NULL
        data.sign.2t <- data.frame(Attribute=c("UpRegulated", "DownRegulated"),
                                   Number.DE=c(Pos.2t, Neg.2t))

        LVLSattribute <- rev(levels(factor(data.sign.2t$Attribute)))
        data.sign.2t$Attribute <- factor(data.sign.2t$Attribute,
                                         levels=LVLSattribute)

        barplot2t <- ggplot2::ggplot(data=data.sign.2t,
                                     ggplot2::aes(x="t1 vs t0", y=Number.DE,
                                                  fill=Attribute)) +
            ggplot2::xlab("DE") + ggplot2::ylab("Number of DE genes")+
            ggplot2::geom_bar(stat="identity") +
            ggplot2::scale_fill_manual(values=c("#E41A1C", "steelblue")) +
            ggplot2::theme_minimal()
        ## ggplot2::scale_fill_brewer(palette="Paired")

        listDEplotsTIME[[1]] <- barplot2t
        names(listDEplotsTIME)[1] <- "Number.Up.Down.Regulated"
    }## if(Nb.time==2)

    ##-----------------------------------------------------------------------##
    if (Nb.time>2) {
        title.alluvial.T <- paste("Alluvial graph of genes which are DE",
                                  "at least at one time",
                                  sep=" ")
        title.evolution.T <- paste("Time evolution of the number of genes",
                                   "which are DE at least at one time within",
                                   "each temporal group",
                                   sep=" ")

        res.allu.g <- DEplotAlluvial(table.DE.time=table.DE.Times,
                                     Temporal.Group=TRUE,
                                     title.alluvial=title.alluvial.T,
                                     title.evolution=title.evolution.T)

        listDEplotsTIME[[1]] <- res.allu.g$g.alluvial
        names(listDEplotsTIME)[1] <- "Alluvial.graph"

        listDEplotsTIME[[2]] <- res.allu.g$g.alluvial.freq
        names(listDEplotsTIME)[2] <- paste0("NumberDEgenes_acrossTime",
                                            "_perTemporalGroup")

        ##-------------------------------------------------------------------##
        ## Number DE gene per time with or without sign
        res.nb.DEPerTime <- DEplotBarplotTime(table.DE.time=table.DE.Times,
                                              Log2.FC.matrix=table.log2FC)

        listDEplotsTIME[[3]] <- res.nb.DEPerTime$g.nb.DEPerTime.sign
        names(listDEplotsTIME)[3] <- "NumberDEgenes_UpDownRegulated_perTime"

        ##-------------------------------------------------------------------##
        ## Upset Venn graph with or not over vs t0
        res.VennBarplot <- DEplotVennBarplotTime(table.DE.time=table.DE.Times,
                                                 Log2.FC.matrix=table.log2FC)

        listDEplotsTIME[[4]] <- res.VennBarplot$Upset.graph
        names(listDEplotsTIME)[4] <- "VennBarplot"

        listDEplotsTIME[[5]] <- res.VennBarplot$Upset.graph.with.nb.over
        names(listDEplotsTIME)[5] <- "VennBarplot_withNumberUpRegulated"
    }## if(Nb.time>2)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 5) Save
    if (!is.null(SubFile.name)) {
        SubFile.name <- paste0("_", SubFile.name)
    }## if(is.null(SubFile.name)==TRUE)

    ##-----------------------------------------------------------------------##
    if (Nb.time == 2) {
        NbUDfile <- paste0("Plot_NumberDEgenes_UpDownRegulated",
                           SubFile.name, ".pdf")

        if (!is.null(path.result)) {
            grDevices::pdf(file=file.path(path.result, NbUDfile),
                           width=11, height=8)
            print(barplot2t)
            grDevices::dev.off()
        }## if(is.null(path.result)==FALSE)

        if (isTRUE(Plot.DE.graph)) {
            print(barplot2t)
        }## if (isTRUE(Plot.DE.graph))

    }## if(Nb.time==2)

    ##-----------------------------------------------------------------------##
    if (Nb.time>2) {
        Allufile <- paste0("Plot_AlluvialGraph", SubFile.name, ".pdf")
        NbaccrosTfile <- paste0("Plot_NumberDEgenes_acrossTime",
                                "_perTemporalGroup", SubFile.name, ".pdf")
        NbUDtfile <- paste0("Plot_NumberDEgenes_UpDownRegulated_perTime",
                            SubFile.name, ".pdf")
        Venn1file <- paste0("Plot_VennBarplot", SubFile.name, ".pdf")
        Venn2file <- paste0("Plot_VennBarplot_withNumberUpRegulated",
                            SubFile.name, ".pdf")

        if (!is.null(path.result)) {
            ##---------------------------------------------------------------##
            ## Save alluvial graph
            grDevices::pdf(file=file.path(path.result, Allufile),
                           width=11, height=8)
            print(res.allu.g$g.alluvial)
            grDevices::dev.off()

            ##---------------------------------------------------------------##
            ## Save graph nb de gene per group per time
            grDevices::pdf(file=file.path(path.result, NbaccrosTfile),
                           width=11, height=8)
            print(res.allu.g$g.alluvial.freq)
            grDevices::dev.off()

            ##---------------------------------------------------------------##
            ## Save graph NB DE per time with or not sign
            grDevices::pdf(file=file.path(path.result, NbUDtfile),
                           width=11, height=8)
            print(res.nb.DEPerTime$g.nb.DEPerTime.sign)
            grDevices::dev.off()

            ##---------------------------------------------------------------##
            grDevices::pdf(file=file.path(path.result, Venn1file),
                           width=11, height=8)
            print(res.VennBarplot$Upset.graph)
            grDevices::dev.off()

            ##---------------------------------------------------------------##
            grDevices::pdf(file=file.path(path.result, Venn2file),
                           width=11, height=8)
            print(res.VennBarplot$Upset.graph.with.nb.over)
            grDevices::dev.off()
            ##---------------------------------------------------------------##
        }## if(is.null(path.result)==FALSE)

        if (isTRUE(Plot.DE.graph)) {
            print(res.allu.g$g.alluvial)
            print(res.allu.g$g.alluvial.freq)

            print(res.nb.DEPerTime$g.nb.DEPerTime.sign)

            print(res.VennBarplot$Upset.graph)
            print(res.VennBarplot$Upset.graph.with.nb.over)
        }## if (isTRUE(Plot.DE.graph))

    }## if(Nb.time>2)

    ##-----------------------------------------------------------------------##
    if (Nb.time>2) {
        DE.Temporal.Pattern.f <- res.VennBarplot$DE.pattern.t.01.sum
    } else {
        DE.Temporal.Pattern.f <- data.frame(matrix(NA, nrow=1, ncol=4))
        colnames(DE.Temporal.Pattern.f) <- c("DE.Temporal.Pattern",
                                             "0.ti.over.t0", "1.ti.over.t0",
                                             "Total")
        ##
        DE.Temporal.Pattern.f[1,] <- c(1, data.sign.2t$Number.DE[c(2, 1)],
                                       sum(data.sign.2t$Number.DE[c(2, 1)]))
    }## if(Nb.time>2)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## SE final
    listDEresultTime <- append(list(DEsummary=Table.Results.f,
                                    DE_TemporalPattern=DE.Temporal.Pattern.f),
                               listDEplotsTIME)

    DESeqclass <- DESeq.result
    S4Vectors::metadata(DESeqclass)$DEresultsTime <- listDEresultTime

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 6) End et outputs
    return(DESeqclass=DESeqclass)
}## DEanalysisTime()
