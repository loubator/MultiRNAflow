#' @title Intermediate analysis when samples belong to different
#' biological conditions
#'
#' @description This function realizes the intermediary steps of the analysis
#' of the function [DEanalysisGroup()].
#'
#' @param DESeq.result Output from the function [DESeq2::DESeq()].
#' @param pval.min Numeric value between 0 and 1. A gene will be considered as
#' differentially expressed (DE) between two biological conditions if
#' its Benjamini-Hochberg adjusted p-value (see [stats::p.adjust()])
#' is below the threshold \code{pval.min}. Default value is 0.05.
#' @param log.FC.min Non negative numeric value.
#' If the log2 fold change between biological conditions or times has
#' an absolute value below the threshold \code{log.FC.min},
#' then the gene is not selected even if is considered as DE.
#' Default value is 1. If \code{log.FC.min=0}, all DE genes will be kept.
#' @param LRT.supp.info \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the algorithm realizes another statistical test in order to
#' detect if, among all biological conditions and/or times, at least one has
#' a different behavior than the others (see the input \code{test}
#' in [DESeq2::DESeq()]).
#'
#' @return The function returns
#' * a data.frame (output \code{Results}) which contains
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
#'   the gene is specific to the biological condition BC1 and expressions in BC1
#'   are higher than in the other biological conditions.
#'   A '-1' in one of these columns means the gene is down-regulated
#'   (or under-expressed) for the biological condition associated to the given
#'   column.
#'   A gene is called regulated for a given biological condition BC1 if
#'   the gene is specific to the biological condition BC1 and expressions in BC1
#'   are lower than in the other biological conditions.
#'   A '0' in one of these columns means the gene is not specific to
#'   the biological condition associated to the given column.
#' * a data.frame (output \code{DE.per.pair.G}) with \eqn{N_g} rows
#' and \eqn{((N_{bc}-1)\times N_{bc})/2} columns with
#' \eqn{N_g} the number of genes
#' and \eqn{N_{bc}} the number of biological conditions.
#' The number of 1 in the n-th row gives the number of pairs of
#' biological conditions where the gene \eqn{n} is DE.
#' The output \code{DE.per.pair.G} will be the input of the function
#' [DEplotVennBarplotGroup()].
#' * a contingency matrix (output \code{Contingence.per.group}) which gives
#' for each biological condition the number of genes categorized as
#' "Upregulated", "DownRugulated" and "Other".
#' A gene is categorized as 'Other', for a given biological condition BC1,
#' if the gene is not specific to the biological condition BC1.
#' The category 'Other' does not exist when there are only two
#' biological conditions.
#' The output \code{Contingence.per.group} will be the input of
#' the function [DEplotBarplot()].
#'
#' @importFrom DESeq2 results
#'
#' @seealso The outputs of the function are used by the main function
#' [DEanalysisGroup()].
#'
#' @export
#'
#' @examples
#' data(RawCounts_Antoszewski2022_MOUSEsub500)
#' # No time points. We take only two groups for the speed of the example
#' RawCounts_T1Wt<-RawCounts_Antoszewski2022_MOUSEsub500[,1:7]
#' DESeq2.info<-DEanalysisPreprocessing(RawCounts=RawCounts_T1Wt,
#'                                      Column.gene=1,
#'                                      Group.position=1,
#'                                      Time.position=NULL,
#'                                      Individual.position=2)
#' #
#' dds.DE.G<-DESeq2::DESeq(DESeq2.info$DESeq2.obj)# ,test = "LRT",reduced=~1
#' res.sum.G<-DEresultGroup(DESeq.result=dds.DE.G,
#'                          LRT.supp.info=FALSE,
#'                          log.FC.min=1,
#'                          pval.min=0.05)

DEresultGroup<-function(DESeq.result,
                        LRT.supp.info=TRUE,
                        pval.min=0.05,
                        log.FC.min=1){
  #---------------------------------------------------------------------------#
  # 0)  Parameters
  #---------------------------------------------------------------------------#
  # Gene names and number of genes
  # DESeq2::results(DESeq.result)
  # Gene.Names<-row.names(SummarizedExperiment::assay(DESeq.result))
  res.dds<-DESeq.result@rowRanges@partitioning@NAMES
  #
  if(is.null(res.dds)==TRUE){# is.null(row.names(res.dds))==TRUE)
    Row.name.res<-paste("Gene", seq_len(length(res.dds)), sep="")
  }else{# # 1:length(res.dds), nrow(res.dds)
    Row.name.res<-res.dds#row.names(res.dds)
  }# if(is.null(row.names(res.dds.group))==TRUE)
  #
  Nb.gene<-length(Row.name.res)#nrow(DESeq2::results(DESeq.result))
  # Biological conditions
  Vector.group<-as.factor(DESeq.result@colData@listData[[1]])
  nb.group<-length(levels(Vector.group))
  nb.pair.of.group<-(nb.group*(nb.group-1))/2
  #---------------------------------------------------------------------------#
  # 1)  Differential expression between each pair of biological condition
  #---------------------------------------------------------------------------#
  Vect.fight.group<-rep(NA, nb.pair.of.group)
  #
  DE.per.2BC<-data.frame(matrix(0, ncol=3*nb.pair.of.group, nrow=Nb.gene))
  Bin.mat.DE<-data.frame(matrix(0, ncol=nb.pair.of.group, nrow=Nb.gene))
  #---------------------------------------------------------------------------#
  if(LRT.supp.info==TRUE){
    res.LRT<-DESeq2::results(DESeq.result, test="LRT")
    padj.LRT<-res.LRT$padj
    if(length(which(is.na(padj.LRT)))>0){
      padj.LRT[which(is.na(padj.LRT))]<-1
    }# if(length(which(is.na(padj.LRT)))>0)
  }# if(LRT.supp.info==TRUE)
  #---------------------------------------------------------------------------#
  cpt<-0
  gene.DE<-c()
  for(i in seq_len(nb.group-1)){
    for(k in seq(from=(i+1), to=nb.group, by=1)){
      cpt<-cpt+1
      res.dds.group<-DESeq2::results(DESeq.result, test="Wald",
                                     contrast=c(names(DESeq.result@colData@listData)[1],
                                                levels(Vector.group)[k],
                                                levels(Vector.group)[i]))
      #
      fight.group<-paste(".", levels(Vector.group)[k], "..",
                         levels(Vector.group)[i], ".", sep="")
      #-----------------------------------------------------------------------#
      Padj.i.VS.k<-res.dds.group$padj
      if(length(which(is.na(Padj.i.VS.k)))>0){
        Padj.i.VS.k[which(is.na(Padj.i.VS.k))]<-1
      }
      Log2.FC.i.VS.k<-res.dds.group$log2FoldChange
      if(length(which(is.na(Log2.FC.i.VS.k)))>0){
        Log2.FC.i.VS.k[which(is.na(Log2.FC.i.VS.k))]<-0
      }
      #-----------------------------------------------------------------------#
      if(LRT.supp.info==TRUE){
        criteria<-sort(intersect(intersect(which(abs(Log2.FC.i.VS.k)>log.FC.min),
                                           which(Padj.i.VS.k<pval.min)),
                                 which(padj.LRT<pval.min)))
      }else{
        criteria<-sort(intersect(which(abs(Log2.FC.i.VS.k)>log.FC.min),
                                 which(Padj.i.VS.k<pval.min)))
      }# if(LRT.supp.info==TRUE)
      #-----------------------------------------------------------------------#
      gene.DE<-c(gene.DE,criteria)
      pvalue.log2FoldChange<-rep(0,nrow(res.dds.group))
      pvalue.log2FoldChange[criteria]<-1
      #
      Bin.mat.DE[cpt]<-pvalue.log2FoldChange
      colnames(Bin.mat.DE)[cpt]<-fight.group
      #
      Vect.fight.group[cpt]<-fight.group
      # round(Padj.i.VS.k,digits=4),
      DE.per.2BC[,3*(cpt-1)+c(1,2,3)]<-data.frame(Log2FC=round(Log2.FC.i.VS.k,
                                                               digits=3),
                                                  Pvalue=Padj.i.VS.k,
                                                  Condition=pvalue.log2FoldChange)
      colnames(DE.per.2BC)[3*(cpt-1)+c(1,2,3)]<-paste(c("Log2FoldChange.",
                                                        "Pvalue.adjusted.",
                                                        "DE."), fight.group,
                                                      sep="")
    }# end for group i
  }# end for group k>i
  #---------------------------------------------------------------------------#
  row.names(Bin.mat.DE)<-row.names(res.dds.group)
  row.names(DE.per.2BC)<-row.names(res.dds.group)
  gene.DE<-sort(unique(gene.DE))
  DE.1min.pair.g<-rep(0, times=Nb.gene)
  DE.1min.pair.g[gene.DE]<-1
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  # 2) Specific genes for each biological condition
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  if(nb.pair.of.group>1){
    # print("Specific genes per biological condition")
    Gene.spe.per.group<-matrix(0, nrow=Nb.gene,
                               ncol=length(levels(Vector.group)))
    colnames(Gene.spe.per.group)<-paste("Specific.genes_",levels(Vector.group),
                                        sep="")
    row.names(Gene.spe.per.group)<-row.names(res.dds.group)
    #
    Nb.DE.per.group<-rep(NA, times=nb.group)
    names(Nb.DE.per.group)<-levels(Vector.group)
    #
    for(g in seq_len(nb.group)){# 1:nb.group
      group.sel<-paste(".", levels(Vector.group)[g], ".", sep="")
      index.group.selec<-grep(pattern=group.sel, x=Vect.fight.group,
                              fixed=TRUE)
      #
      Id.Column.spe<-c(seq_len(nb.pair.of.group)*3)[index.group.selec]
      Id.Column.nospe<-c(seq_len(nb.pair.of.group)*3)[-index.group.selec]
      # 1:nb.pair.of.group
      sum.row.spe<-apply(X=DE.per.2BC[,Id.Column.spe], MARGIN=1, FUN=sum)
      sum.row.no.spe<-apply(X=DE.per.2BC[,Id.Column.nospe], MARGIN=1, FUN=sum)
      #
      Nb.DE.per.group[g]<-length(which(sum.row.spe>0))
      #
      Spe.vect<-rep(0, times=nrow(res.dds.group))
      for(i in seq_len(length(Spe.vect))){# 1:length(Spe.vect)
        if(sum.row.no.spe[i]==0 & sum.row.spe[i]==length(index.group.selec)){
          Spe.vect[i]<-1
        }
      }# for(i in 1:length(Spe.vect))
      Gene.spe.per.group[,g]<-Spe.vect
    }# for(g in 1:nb.group)
  }else{
    Gene.spe.per.group<-cbind(DE.per.2BC[,3], DE.per.2BC[,3])
    colnames(Gene.spe.per.group)<-paste("Specific.genes_",
                                        levels(Vector.group),
                                        sep="")
    row.names(Gene.spe.per.group)<-row.names(res.dds.group)
    Nb.DE.per.group<-rep(length(which(DE.per.2BC[,3]>0)), times=2)
  }# if(nb.group>1)
  #---------------------------------------------------------------------------#
  # 3) Over and under expressed genes per biological condition
  #---------------------------------------------------------------------------#
  OverUnder.expr.per.g<-matrix(0, nrow=Nb.gene, ncol=nb.group)
  colnames(OverUnder.expr.per.g)<-paste("OverUnder.regulated.genes_",
                                        levels(Vector.group), sep="")
  row.names(OverUnder.expr.per.g)<-row.names(res.dds.group)
  #
  if(nb.pair.of.group>1){
    for(g in seq_len(nb.group)){# 1:nb.group
      group.sel<-paste(".", levels(Vector.group)[g], ".", sep="")
      index.group.sel<-grep(pattern=group.sel, x=Vect.fight.group, fixed=TRUE)
      #
      Id.Column.pval<-c(seq_len(nb.pair.of.group)*3)[index.group.sel]
      Id.Column.log2fc<-c(seq_len(nb.pair.of.group)*3-2)[index.group.sel]
      # (1:nb.pair.of.group)
      sign.matrix<-apply(X=DE.per.2BC[,Id.Column.pval]*DE.per.2BC[,Id.Column.log2fc],
                         MARGIN=2, FUN=function(x) sign(x))
      ###
      Position.log0<-matrix(unlist(strsplit(Vect.fight.group[index.group.sel],
                                            split=".." , fixed=TRUE)), nrow=2)
      Position.log1<-gsub(".", "", Position.log0, fixed=TRUE)
      Position.log<-matrix(paste(".", Position.log1, "." ,sep=""),
                           ncol=length(index.group.sel), byrow=FALSE)
      vec.Position.log<-apply(Position.log, MARGIN=2,
                              FUN=function(x) which(x==group.sel))
      vec.Position.log[which(vec.Position.log==2)]<--1
      ###
      sign.matrix2<-sign.matrix*matrix(rep(vec.Position.log,
                                           times=nrow(sign.matrix)),
                                       nrow=nrow(sign.matrix), byrow=TRUE)
      sum.sign.matrix<-as.numeric(apply(X=sign.matrix2, MARGIN=1, FUN=sum))
      sum.sign.matrix.f<-sum.sign.matrix*Gene.spe.per.group[,g]
      #
      for(i in seq_len(nrow(OverUnder.expr.per.g))){
        sign.gene.i<-sum.sign.matrix.f[i]
        if(sign.gene.i==length(index.group.sel) & is.na(sign.gene.i)==FALSE){
          OverUnder.expr.per.g[i,g]<-1
        }
        if(sign.gene.i==-length(index.group.sel) & is.na(sign.gene.i)==FALSE){
          OverUnder.expr.per.g[i,g]<- -1
        }
      }# end for(i in 1:nrow(OverUnder.expr.per.g))
    }# end for (g in 1:nb.group)
  }else{
    OverUnder.expr.per.g[,1]<--sign(pvalue.log2FoldChange*Log2.FC.i.VS.k)
    OverUnder.expr.per.g[,2]<-sign(pvalue.log2FoldChange*Log2.FC.i.VS.k)
  }# if(nb.pair.of.group>1)
  #
  contin.spe.g.ini<-rbind(apply(OverUnder.expr.per.g, 2,
                                function(x) length(which(x==1))),
                          apply(OverUnder.expr.per.g, 2,
                                function(x) length(which(x==-1))))
  delta.spe.sign.spe<-apply(Gene.spe.per.group,2,sum) - apply(contin.spe.g.ini,2,sum)
  #
  contin.spe.g<-rbind(rbind(contin.spe.g.ini,delta.spe.sign.spe),
                      Nb.DE.per.group-apply(rbind(contin.spe.g.ini,delta.spe.sign.spe),
                                            2, sum))
  colnames(contin.spe.g)<-levels(Vector.group)
  row.names(contin.spe.g)<-c("UpRegulated","DownRegulated",
                             "No.specific","Other")
  #
  if(nb.pair.of.group>1){
    contin.spe.g.f<-contin.spe.g[-3,]
  }else{
    contin.spe.g.f<-contin.spe.g[-c(3,4),]
  }# if(nb.pair.of.group>1)
  #---------------------------------------------------------------------------#
  # 4) Table with all results
  #---------------------------------------------------------------------------#
  # print("Summary all steps")
  Resultats.DEseq2.groups<-data.frame(Gene=Row.name.res,
                                      DE.1pair.of.Group.minimum=DE.1min.pair.g,
                                      cbind(Gene.spe.per.group,
                                            OverUnder.expr.per.g,
                                            DE.per.2BC))
  row.names(Resultats.DEseq2.groups)<-row.names(res.dds.group)
  #
  colnames(Resultats.DEseq2.groups)<-gsub("Log2FoldChange..","Log2FoldChange_",
                                          colnames(Resultats.DEseq2.groups),
                                          fixed=TRUE)
  colnames(Resultats.DEseq2.groups)<-gsub("Pvalue.adjusted..",
                                          "Pvalue.adjusted_",
                                          colnames(Resultats.DEseq2.groups),
                                          fixed=TRUE)
  colnames(Resultats.DEseq2.groups)<-gsub("DE..","DE_",
                                          colnames(Resultats.DEseq2.groups),
                                          fixed=TRUE)
  colnames(Resultats.DEseq2.groups)<-gsub("..",".versus.",
                                          colnames(Resultats.DEseq2.groups),
                                          fixed=TRUE)
  #
  #---------------------------------------------------------------------------#
  # 5) END
  #---------------------------------------------------------------------------#
  return(list(Results=Resultats.DEseq2.groups,
              DE.per.pair.G=Bin.mat.DE,
              Contingence.per.group=contin.spe.g.f))
}# DEresultGroup()
