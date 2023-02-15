#' @title Barplot of DE genes per time
#'
#' @description The function takes as input two tables
#' * a binary table with \eqn{N_g} rows corresponding to genes and
#' \eqn{T-1} columns corresponding to times
#' (with \eqn{T} the number of time points).
#' A '1' in the n-th row and i-th column means that the n-th gene is
#' differentially expressed (DE) at time ti,
#' compared with the reference time t0.
#' * a numeric matrix with positive and negative values with
#' \eqn{N_g} rows corresponding to genes and \eqn{T-1} columns corresponding
#' to times.
#' The element in n-th row and i-th column corresponds to the log2 fold change
#' between the time ti and the reference time t0 for the n-th gene.
#' If the gene is DE and the sign is positive, then the gene n will be
#' considered as over-expressed (up-regulated) at the time ti.
#' If the gene is DE and the sign is negative, then the gene n will be
#' considered as under-expressed (down-regulated) at the time ti.
#'
#' The function plots two graphs: a barplot showing the number of DE genes
#' per time and a barplot showing the number of under- and over-expressed
#' genes per times.
#'
#' @param table.DE.time Binary matrix (table filled with 0 and 1) with
#' \eqn{N_g} rows and \eqn{T-1} columns with
#' \eqn{N_g} the number of genes and \eqn{T} the number of time points.
#' @param Log2.FC.matrix Numeric matrix with positive and negative with
#' \eqn{N_g} rows and \eqn{T-1} columns.
#'
#' @return The function plots two graphs:
#' a barplot showing the number of DE genes per time and
#' a barplot showing the number of under and over expressed genes per times.
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_brewer xlab ylab
#' theme_minimal labs
#'
#' @export
#'
#' @examples
#' set.seed(1994)
#' Dat1.FTP=matrix(sample(c(0,1),replace=TRUE,size=120,c(0.3,0.7)),ncol=3)
#' Dat2.FTP=matrix(round(rnorm(n=120, mean=0, sd=1),digits=2), ncol=3)
#' colnames(Dat1.FTP)=paste("t",1:3,sep="")
#' colnames(Dat2.FTP)=paste("t",1:3,sep="")
#' #--------------------------------------------------------------------------#
#' res.DE.all.t=DEplotBarplotTime(table.DE.time=Dat1.FTP,
#'                                Log2.FC.matrix=Dat2.FTP)
#' print(res.DE.all.t$g.nb.DEPerTime)
#' print(res.DE.all.t$g.nb.DEPerTime.sign)

DEplotBarplotTime<-function(table.DE.time,
                            Log2.FC.matrix){
  #---------------------------------------------------------------------------#
  # 1) Parameters
  #---------------------------------------------------------------------------#
  # DE gene
  if(sum(abs(table.DE.time))==0){stop("NO DE genes.")}
  gene.DE.number<-which(apply(table.DE.time,MARGIN=1,function(x) sum(x))!=0)
  table.DE.time.DE.gene<-table.DE.time[gene.DE.number,]
  if(is.null(Log2.FC.matrix)==FALSE){
    Log2.FC.matrix.DE<-Log2.FC.matrix[gene.DE.number,]
  }# if(is.null(Log2.FC.matrix)==FALSE)
  # To avoid "no visible binding for global variable" with devtools::check()
  Freq<-Time<-Characteristic<-NULL
  #---------------------------------------------------------------------------#
  # 2) Data pour Graph
  #---------------------------------------------------------------------------#
  # Data for computing the number of DE genes per time
  Sum.per.time<-apply(table.DE.time.DE.gene, MARGIN=2, FUN=sum)
  Nb.DE.genes.per.time<-data.frame(Time=colnames(table.DE.time.DE.gene),
                                   Freq=Sum.per.time)
  if(is.null(Log2.FC.matrix)==FALSE){
    freq.sign<-matrix(data=0, ncol=ncol(table.DE.time.DE.gene), nrow=2)
    product.m<-Log2.FC.matrix.DE*table.DE.time.DE.gene
    colnames(freq.sign)<-colnames(table.DE.time.DE.gene)
    #
    for(i in seq_len(ncol(table.DE.time.DE.gene))){
      Nb.negative<-length(which(sign(product.m)[,i]==-1))
      Nb.positive<-length(which(sign(product.m)[,i]==1))
      freq.sign[c(1,2),i]<-c(Nb.negative,Nb.positive)
    }# for(i in seq_len(ncol(table.DE.time.DE.gene)))
    #
    Sign.freq.per.time<-data.frame(Characteristic=c("DownRegulated",
                                                    "UpRegulated"),
                                   freq.sign)
    freq.sign.mat<-reshape2::melt(data=Sign.freq.per.time,
                                 id.vars=c("Characteristic"))
    colnames(freq.sign.mat)<-c("Characteristic", "Time", "Freq")
    #
    Characteristic.Levels<-rev(levels(factor(freq.sign.mat$Characteristic)))
    freq.sign.mat$Characteristic<-factor(freq.sign.mat$Characteristic,
                                         levels=Characteristic.Levels)
  }# if(is.null(Log2.FC.matrix)==FALSE)
  #---------------------------------------------------------------------------#
  # 3) Graphs
  #---------------------------------------------------------------------------#
  g.barplot<-ggplot2::ggplot(Nb.DE.genes.per.time,
                             ggplot2::aes(y=Freq, x=Time))+
    ggplot2::geom_bar(position="stack", stat="identity", fill="#999999",
                      color="black") +
    ggplot2::xlab("Time")+
    ggplot2::ylab("Number of DE gene")+
    ggplot2::theme_minimal()
  if(is.null(Log2.FC.matrix)==FALSE){
    ggplot2Title<-"Number of up- and down-regulated genes per time versus t0"
    #
    g.barplot.sign<-ggplot2::ggplot(freq.sign.mat,
                                    ggplot2::aes(fill=Characteristic,
                                                 y=Freq, x=Time))+
      ggplot2::geom_bar(position="stack", stat="identity", color="black") +
      ggplot2::scale_fill_manual(values=c("#E41A1C", "steelblue"))+
      ggplot2::xlab("Time")+
      ggplot2::ylab("Number of genes")+
      ggplot2::ggtitle(ggplot2Title)+
      ggplot2::labs(fill='Attribute') +
      ggplot2::theme_minimal()
  }else{g.barplot.sign<-NULL}
  #---------------------------------------------------------------------------#
  # 4) End
  #---------------------------------------------------------------------------#
  return(list(g.nb.DEPerTime=g.barplot,
              g.nb.DEPerTime.sign=g.barplot.sign))
}# DEplotBarplotTime()
