#' @title Intermediate analysis when samples belong to different
#' biological conditions
#'
#' @description This function realizes the intermediary steps of the analysis
#' of the function
#' [DEanalysisGroup()].
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
#' an absolute value below the threshold \code{log.FC.min},
#' then the gene is not selected even if is considered as DE.
#' Default value is 1. If \code{log.FC.min=0}, all DE genes will be kept.
#' @param LRT.supp.info \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the algorithm realizes another statistical test in order to
#' detect if, among all biological conditions and/or times, at least one has
#' a different behavior than the others (see the input \code{test} in
#' [DESeq2::DESeq()]).
#'
#'
#' @return The function returns the same DESeqDataSet class object
#' \code{DESeq.result} with the following results,
#' saved in the metadata \code{DEresultsGroup} of \code{DESeq.result}:
#' * a data.frame (output \code{DEsummary} of \code{DEresultsGroup})
#' which contains
#'   * gene names
#'   * pvalues, log2 fold change and DE genes between each pairs of
#'   biological conditions.
#'   * a binary column (1 and 0) where 1 means the gene is DE between at least
#'   one pair of biological conditions.
#'   * \eqn{N_{bc}} binary columns, where \eqn{N_{bc}} is the number of
#'   biological conditions, which gives the specific genes for each
#'   biological condition.
#'   A '1' in one of these columns means the gene is specific to the
#'   biological condition associated to the given column. 0 otherwise.
#'   A gene is called specific to a given biological condition BC1,
#'   if the gene is DE between BC1 and any other biological conditions,
#'   but not DE between any pair of other biological conditions.
#'   * \eqn{N_{bc}} columns filled with -1, 0 and 1, one per
#'   biological condition.
#'   A '1' in one of these columns means the gene is up-regulated
#'   (or over-expressed) for the biological condition associated
#'   to the given column.
#'   A gene is called up-regulated for a given biological condition BC1 if
#'   the gene is specific to the biological condition BC1 and expressions in
#'   BC1 are higher than in the other biological conditions.
#'   A '-1' in one of these columns means the gene is down-regulated
#'   (or under-expressed) for the biological condition associated to the given
#'   column.
#'   A gene is called regulated for a given biological condition BC1 if
#'   the gene is specific to the biological condition BC1 and expressions in
#'   BC1 are lower than in the other biological conditions.
#'   A '0' in one of these columns means the gene is not specific to
#'   the biological condition associated to the given column.
#' * a data.frame (output \code{DE.per.pair.G} of \code{DEresultsGroup})
#' with \eqn{N_g} rows and \eqn{((N_{bc}-1)\times N_{bc})/2} columns
#' with \eqn{N_g} the number of genes
#' and \eqn{N_{bc}} the number of biological conditions.
#' The number of 1 in the n-th row gives the number of pairs of
#' biological conditions where the gene \eqn{n} is DE.
#' The output \code{DE.per.pair.G} will be the input of the function
#' [DEplotVennBarplotGroup()].
#' * a contingency matrix (output \code{Contingence.per.group}
#' of \code{DEresultsGroup}) which gives for each biological condition
#' the number of genes categorized as
#' "Upregulated", "DownRugulated" and "Other".
#' A gene is categorized as 'Other', for a given biological condition BC1,
#' if the gene is not specific to the biological condition BC1.
#' The category 'Other' does not exist when there are only two
#' biological conditions.
#'
#' The output \code{Contingence.per.group} will be the input of the function
#' [DEplotBarplot()].
#'
#' @importFrom DESeq2 results
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#'
#' @seealso The output of the function are used by the main function
#' [DEanalysisGroup()].
#'
#' @export
#'
#' @examples
#' ## Data
#' data("RawCounts_Antoszewski2022_MOUSEsub500")
#' ## No time points. We take only two groups for the speed of the example
#' RawCounts_T1Wt <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200),
#'                                                         seq_len(7)]
#'
#' ## Preprocessing step
#' resDATAprepSEmus1 <- DATAprepSE(RawCounts=RawCounts_T1Wt,
#'                                 Column.gene=1,
#'                                 Group.position=1,
#'                                 Time.position=NULL,
#'                                 Individual.position=2)
#'
#' DESeq2preprocess <- S4Vectors::metadata(resDATAprepSEmus1)$DESeq2obj
#' DESeq2obj <- DESeq2preprocess$DESeq2preproceesing
#'
#' ##------------------------------------------------------------------------##
#' dds.DE.G <- DESeq2::DESeq(DESeq2obj)
#'
#' res.sum.G <- DEresultGroup(DESeq.result=dds.DE.G,
#'                            LRT.supp.info=FALSE,
#'                            log.FC.min=1,
#'                            pval.min=0.05)

DEresultGroup <- function(DESeq.result,
                          LRT.supp.info=TRUE,
                          pval.min=0.05,
                          log.FC.min=1) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (!is(DESeq.result, 'DESeqDataSet')) {
        stop("Res.DE.analysis must be a 'DESeqDataSet' object")
    }## if(!is(classDeseq2, 'DESeqDataSet'))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 0)  Parameters
    ## Gene names and number of genes
    geneNames <- dimnames(DESeq.result)[[1]] ## DESeq2::results(DESeq.result)

    if (is.null(geneNames)) {
        Row.name.res <- paste0("Gene", seq_len(length(geneNames)))
    } else {
        Row.name.res <- geneNames
    }## if(is.null(row.names(resDDSgroup))==TRUE)

    Nb.gene <- length(Row.name.res)

    ## Biological conditions
    Vector.group <- data.frame(SummarizedExperiment::colData(DESeq.result))[,1]
    Vector.group <- as.factor(as.character(Vector.group))

    nb.group <- length(levels(Vector.group))
    NpairGroups <- (nb.group*(nb.group-1))/2

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 1)  Differential expression between each pair of biological condition
    Vect.fight.group <- rep(NA, NpairGroups)

    DEper2BC <- data.frame(matrix(0, ncol=3*NpairGroups, nrow=Nb.gene))
    Bin.mat.DE <- data.frame(matrix(0, ncol=NpairGroups, nrow=Nb.gene))

    ##-----------------------------------------------------------------------##
    if (LRT.supp.info==TRUE) {
        res.LRT <- DESeq2::results(DESeq.result, test="LRT")
        padj.LRT <- res.LRT$padj

        if (length(which(is.na(padj.LRT)))>0) {
            padj.LRT[which(is.na(padj.LRT))] <- 1
        }## if(length(which(is.na(padj.LRT)))>0)
    }## if(LRT.supp.info==TRUE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    NameColDat <- names(SummarizedExperiment::colData(DESeq.result))[1]

    cpt <- 0
    gene.DE <- c()

    for (i in seq_len(nb.group-1)) {
        for (k in seq(from=(i+1), to=nb.group, by=1)) {
            cpt <- cpt+1
            resDDSgroup <- DESeq2::results(DESeq.result, test="Wald",
                                           contrast=c(NameColDat,
                                                      levels(Vector.group)[k],
                                                      levels(Vector.group)[i]))

            fight.group <- paste0(".", levels(Vector.group)[k], "..",
                                  levels(Vector.group)[i], ".")

            ##---------------------------------------------------------------##
            Padj.i.VS.k <- resDDSgroup$padj
            if (length(which(is.na(Padj.i.VS.k)))>0) {
                Padj.i.VS.k[which(is.na(Padj.i.VS.k))] <- 1
            }## if (length(which(is.na(Padj.i.VS.k)))>0)

            Log2.FC.i.VS.k <- resDDSgroup$log2FoldChange
            if (length(which(is.na(Log2.FC.i.VS.k)))>0) {
                Log2.FC.i.VS.k[which(is.na(Log2.FC.i.VS.k))] <- 0
            }## if (length(which(is.na(Log2.FC.i.VS.k)))>0)

            ##---------------------------------------------------------------##
            criteria <- sort(intersect(which(abs(Log2.FC.i.VS.k)>log.FC.min),
                                       which(Padj.i.VS.k<pval.min)))
            if (isTRUE(LRT.supp.info)) {
                criteria <- sort(intersect(criteria, which(padj.LRT<pval.min)))
            }## if(LRT.supp.info==TRUE)

            ##---------------------------------------------------------------##
            gene.DE <- c(gene.DE,criteria)
            pvalue.log2FoldChange <- rep(0,nrow(resDDSgroup))
            pvalue.log2FoldChange[criteria] <- 1

            Bin.mat.DE[cpt] <- pvalue.log2FoldChange
            colnames(Bin.mat.DE)[cpt] <- fight.group

            Vect.fight.group[cpt] <- fight.group

            IDcol2BC <- 3*(cpt - 1) + c(1, 2, 3)
            DEper2BC[,IDcol2BC] <- data.frame(Log2FC=round(Log2.FC.i.VS.k,
                                                           digits=3),
                                              Pvalue=Padj.i.VS.k,
                                              Condition=pvalue.log2FoldChange)
            colnames(DEper2BC)[IDcol2BC] <- paste0(c("Log2FoldChange.",
                                                     "Pvalue.adjusted.",
                                                     "DE."),
                                                   fight.group)
            ## round(Padj.i.VS.k,digits=4),
        }# end for group i
    }# end for group k>i

    ##-----------------------------------------------------------------------##
    row.names(Bin.mat.DE) <- row.names(resDDSgroup)
    row.names(DEper2BC) <- row.names(resDDSgroup)
    gene.DE <- sort(unique(gene.DE))
    DE.1min.pair.g <- rep(0, times=Nb.gene)
    DE.1min.pair.g[gene.DE] <- 1

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 2) Specific genes for each biological condition

    if (NpairGroups>1) {
        ## print("Specific genes per biological condition")
        SPEgenesPerGroup <- matrix(0, nrow=Nb.gene,
                                   ncol=length(levels(Vector.group)))
        colnames(SPEgenesPerGroup) <- paste0("Specific.genes_",
                                             levels(Vector.group))
        row.names(SPEgenesPerGroup) <- row.names(resDDSgroup)

        Nb.DE.per.group <- rep(NA, times=nb.group)
        names(Nb.DE.per.group) <- levels(Vector.group)

        for(g in seq_len(nb.group)){
            group.sel <- paste0(".", levels(Vector.group)[g], ".")
            IDgroupSelect <- grep(pattern=group.sel,
                                  x=Vect.fight.group,
                                  fixed=TRUE)

            IDcolSPE <- c(seq_len(NpairGroups)*3)[IDgroupSelect]
            IDcolNOspe <- c(seq_len(NpairGroups)*3)[-IDgroupSelect]

            SUMrowSPE <- apply(X=data.frame(DEper2BC[,IDcolSPE]),
                               MARGIN=1, FUN=sum)
            SUMrowNOspe <- apply(X=data.frame(DEper2BC[,IDcolNOspe]),
                                 MARGIN=1, FUN=sum)

            Nb.DE.per.group[g] <- length(which(SUMrowSPE>0))

            Spe.vect <- rep(0, times=nrow(resDDSgroup))
            for (i in seq_len(length(Spe.vect))) {
                if(SUMrowNOspe[i] == 0 & SUMrowSPE[i] == length(IDgroupSelect)){
                    Spe.vect[i] <- 1
                }## if()
            }## for (i in seq_len(length(Spe.vect)))
            SPEgenesPerGroup[,g] <- Spe.vect
        }## for(g in seq_len(nb.group))

    } else {
        SPEgenesPerGroup <- cbind(DEper2BC[,3], DEper2BC[,3])
        colnames(SPEgenesPerGroup) <- paste0("Specific.genes_",
                                             levels(Vector.group))
        row.names(SPEgenesPerGroup) <- row.names(resDDSgroup)
        Nb.DE.per.group <- rep(length(which(DEper2BC[,3]>0)), times=2)
    }## if(nb.group>1)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 3) Over and under expressed genes per biological condition
    OverUnder.expr.per.g <- matrix(0, nrow=Nb.gene, ncol=nb.group)
    colnames(OverUnder.expr.per.g) <- paste0("OverUnder.regulated.genes_",
                                             levels(Vector.group))
    row.names(OverUnder.expr.per.g) <- row.names(resDDSgroup)

    if (NpairGroups>1) {
        for (g in seq_len(nb.group)) {
            group.sel <- paste0(".", levels(Vector.group)[g], ".")
            IDgroupSelect <- grep(pattern=group.sel,
                                  x=Vect.fight.group,
                                  fixed=TRUE)

            IDcolPVAL <- c(seq_len(NpairGroups)*3)[IDgroupSelect]
            IDcolLog2FC <- c(seq_len(NpairGroups)*3-2)[IDgroupSelect]

            sign.matrix <- apply(X=DEper2BC[,IDcolPVAL]*DEper2BC[,IDcolLog2FC],
                                 MARGIN=2, FUN=function(x) sign(x))

            ##---------------------------------------------------------------##
            Position.log0_strsplit <- strsplit(Vect.fight.group[IDgroupSelect],
                                               split=".." , fixed=TRUE)
            Position.log0 <- matrix(unlist(Position.log0_strsplit), nrow=2)
            Position.log1 <- gsub(".", "", Position.log0, fixed=TRUE)

            Position.log <- matrix(paste0(".", Position.log1, "."),
                                   ncol=length(IDgroupSelect), byrow=FALSE)
            vec.Position.log <- apply(Position.log, MARGIN=2,
                                      FUN=function(x) which(x == group.sel))
            vec.Position.log[which(vec.Position.log == 2)] <- -1

            ##---------------------------------------------------------------##
            sign.matrix2 <- sign.matrix*matrix(rep(vec.Position.log,
                                                   times=nrow(sign.matrix)),
                                               nrow=nrow(sign.matrix),
                                               byrow=TRUE)
            sum.sign.matrix <- as.numeric(apply(X=sign.matrix2, MARGIN=1,
                                                FUN=sum))
            sum.sign.matrix.f <- sum.sign.matrix*SPEgenesPerGroup[,g]

            for (i in seq_len(nrow(OverUnder.expr.per.g))) {
                sign.gene.i <- sum.sign.matrix.f[i]
                if(sign.gene.i == length(IDgroupSelect) & !is.na(sign.gene.i)){
                    OverUnder.expr.per.g[i, g] <- 1
                }
                if(sign.gene.i == -length(IDgroupSelect) & !is.na(sign.gene.i)){
                    OverUnder.expr.per.g[i, g] <- -1
                }
            }# end for(i in 1:nrow(OverUnder.expr.per.g))
        }# end for (g in 1:nb.group)
    } else {
        OverUnder.expr.per.g[,1] <- -sign(pvalue.log2FoldChange*Log2.FC.i.VS.k)
        OverUnder.expr.per.g[,2] <- sign(pvalue.log2FoldChange*Log2.FC.i.VS.k)
    }# if(NpairGroups>1)

    ##-----------------------------------------------------------------------##
    contin.spe.g.ini <- rbind(apply(OverUnder.expr.per.g, 2,
                                    function(x) length(which(x==1))),
                              apply(OverUnder.expr.per.g, 2,
                                    function(x) length(which(x==-1))))
    delta.spe.sign.spe <- apply(SPEgenesPerGroup, 2,
                                sum) - apply(contin.spe.g.ini, 2, sum)

    contin.spe.g <- rbind(rbind(contin.spe.g.ini, delta.spe.sign.spe),
                          Nb.DE.per.group - apply(rbind(contin.spe.g.ini,
                                                        delta.spe.sign.spe),
                                                  2, sum))
    colnames(contin.spe.g) <- levels(Vector.group)
    row.names(contin.spe.g) <- c("UpRegulated", "DownRegulated",
                                 "No.specific", "Other")

    if (NpairGroups>1) {
        contin.spe.g.f <- contin.spe.g[-3,]
    } else {
        contin.spe.g.f <- contin.spe.g[-c(3, 4),]
    }## if(NpairGroups>1)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 4) Table with all results
    ## print("Summary all steps")
    ResultsDEseq2groups <- data.frame(Gene=Row.name.res,
                                      DE.1pair.of.Group.minimum=DE.1min.pair.g,
                                      cbind(SPEgenesPerGroup,
                                            OverUnder.expr.per.g,
                                            DEper2BC))
    row.names(ResultsDEseq2groups) <- row.names(resDDSgroup)

    colnames(ResultsDEseq2groups) <- gsub("Log2FoldChange..",
                                          "Log2FoldChange_",
                                          colnames(ResultsDEseq2groups),
                                          fixed=TRUE)
    colnames(ResultsDEseq2groups) <- gsub("Pvalue.adjusted..",
                                          "Pvalue.adjusted_",
                                          colnames(ResultsDEseq2groups),
                                          fixed=TRUE)
    colnames(ResultsDEseq2groups) <- gsub("DE..", "DE_",
                                          colnames(ResultsDEseq2groups),
                                          fixed=TRUE)
    colnames(ResultsDEseq2groups) <- gsub("..", ".versus.",
                                          colnames(ResultsDEseq2groups),
                                          fixed=TRUE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## SE object
    listDEresultGroup <- list(DEsummary=ResultsDEseq2groups,
                              DE.per.pair.G=Bin.mat.DE,
                              Contingence.per.group=contin.spe.g.f)

    DESeqclass <- DESeq.result
    S4Vectors::metadata(DESeqclass)$DEresultsGroup <- listDEresultGroup

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    return(DESeqclass=DESeqclass)
}## DEresultGroup()
