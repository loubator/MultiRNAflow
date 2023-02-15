#' @title Venn barplot of DE genes across time.
#'
#' @description The function takes as input two matrix or data.frame
#' * a binary matrix or data.frame with \eqn{N_g} rows corresponding to genes
#' and \eqn{T-1} columns corresponding to times
#' (with \eqn{T} the number of time points).
#' A '1' in the n-th row and i-th column means that the n-th gene is
#' differentially expressed (DE) at time ti, compared with
#' the reference time t0.
#' * a numeric matrix or data.frame with \eqn{N_g} rows corresponding to genes
#' and \eqn{T-1} columns corresponding to times.
#' The element in n-th row and i-th column corresponds to the \eqn{log_2}
#' fold change between the time ti and the reference time t0 for the n-th gene.
#' If the gene is DE and the sign is positive, then the gene n will be
#' considered as over-expressed (up-regulated) at time ti.
#' If the gene is DE and the sign is negative, then the gene n will be
#' considered as under-expressed (down-regulated) at time ti.
#'
#' @param table.DE.time Binary matrix or data.frame (table filled with 0 and 1)
#' with \eqn{N_g} rows and \eqn{T-1} columns with \eqn{N_g} the number of genes
#' and \eqn{T} the number of time points.
#' @param Log2.FC.matrix Numeric matrix or data.frame with \eqn{N_g} rows and
#' \eqn{T-1} columns.
#'
#' @return The function plots
#' * the number of genes per time patterns in an UpSet plot (Venn diagram
#' displayed as a barplot) with the R function [UpSetR::upset()].
#' By temporal pattern, we mean the set of times ti such that the gene is
#' DE between ti and the reference time t0.
#' * a similar UpSet plot where each bar is split in different colors
#' corresponding to all possible numbers of DE times where genes are over
#' expressed in a given temporal pattern.
#'
#' @seealso The function
#' * calls the function [UpSetR::upset()] in order to plot the UpSet plot.
#' * is called by the functions [DEanalysisTime()] and
#' [DEanalysisTimeAndGroup()].
#'
#' @importFrom UpSetR elements upset
#' @importFrom scales hue_pal
#' @importFrom grDevices rainbow
#'
#' @export
#'
#' @examples
#' set.seed(1994)
#' Nb.Time<-4 # Number of time measurement
#' #--------------------------------------------------------------------------#
#' table.DE.time.ex=matrix(sample(c(0,1), replace=TRUE,
#'                                size=40*(Nb.Time-1),c(0.2,0.8)),
#'                         ncol=Nb.Time-1)
#' colnames(table.DE.time.ex)=paste("t",1:(Nb.Time-1),sep="")
#' #--------------------------------------------------------------------------#
#' Log2.FC.matrix.ex=matrix(round(rnorm(n=40*(Nb.Time-1), mean=0, sd=1),
#'                                digits=2),
#'                          ncol=(Nb.Time-1))
#' colnames(Log2.FC.matrix.ex)=paste("t",1:(Nb.Time-1),sep="")
#' #--------------------------------------------------------------------------#
#' res.test.VennBarplot=DEplotVennBarplotTime(table.DE.time=table.DE.time.ex,
#'                                            Log2.FC.matrix=Log2.FC.matrix.ex)
#' print(res.test.VennBarplot$Upset.graph)
#' print(res.test.VennBarplot$Upset.graph.with.nb.over)
#' res.test.VennBarplot$DE.pattern.t.01.sum

DEplotVennBarplotTime<-function(table.DE.time,
                                Log2.FC.matrix){
  #---------------------------------------------------------------------------#
  # 1) DE genes
  #---------------------------------------------------------------------------#
  gene.DE.number<-which(apply(table.DE.time, MARGIN=1, function(x) sum(x))!=0)
  table.DE.time.DE.gene<-table.DE.time[gene.DE.number,]
  if(is.null(Log2.FC.matrix)==FALSE){
    Log2.FC.matrix.DE<-Log2.FC.matrix[gene.DE.number,]
  }# if(is.null(Log2.FC.matrix)==FALSE)
  #---------------------------------------------------------------------------#
  # 2) Data for the function UpSetR::upset()
  #---------------------------------------------------------------------------#
  if(is.null(Log2.FC.matrix)==FALSE){
    Sign.DE.t0<-sign(Log2.FC.matrix[gene.DE.number,]*table.DE.time.DE.gene)
    Sign.pat.t<-apply(Sign.DE.t0, MARGIN=1, FUN=function(x) length(which(x==1)))
    Over.Under.pattern<-data.frame(Nb.over=Sign.pat.t, abs(Sign.DE.t0))
    #-------------------------------------------------------------------------#
    id.over<-sort(unique(Over.Under.pattern$Nb.over))
    TpsNoOver<-which(c(0:ncol(table.DE.time))%in%id.over==FALSE)
    #
    if(length(TpsNoOver)>0){
      Add.OUpat<-matrix(0, nrow=length(TpsNoOver), ncol=ncol(table.DE.time)+1)
      Add.OUpat[,1]<-TpsNoOver-1
      colnames(Add.OUpat)<-colnames(Over.Under.pattern)
      #
      OU.pat.f<-rbind(Over.Under.pattern, Add.OUpat)
    }else{
      OU.pat.f<-as.data.frame(Over.Under.pattern)
      row.names(OU.pat.f)<-NULL
    }# if(length(TpsNoOver)>0)
    #-------------------------------------------------------------------------#
    # Temporal pattern
    Vect.t.01<-apply(Over.Under.pattern[,-1], MARGIN=1,
                     FUN=function(x) paste(as.character(x), collapse=""))
    #
    Pat.t.01<-table(Vect.t.01, Over.Under.pattern$Nb.over)
    #
    Sum.pat.t.01<-apply(Pat.t.01, MARGIN=1, FUN=sum)
    #
    DE.temporal.pattern<-data.frame(Binary.time.DE=row.names(Pat.t.01),
                                    cbind(Pat.t.01, Sum.pat.t.01))
    colnames(DE.temporal.pattern)<-c("Temporal.Pattern",
                                     paste(colnames(Pat.t.01),".ti.over.t0",
                                           sep=""),
                                     "Total")
    #
    DE.temporal.pattern.f<-DE.temporal.pattern[order(DE.temporal.pattern$Total,
                                                     decreasing=TRUE),]
    row.names(DE.temporal.pattern.f)<-NULL
    #-------------------------------------------------------------------------#
    # Queries preprocessing
    n.over<-0
    mylist.upsetg<-list()
    vectorUniqueValue<-id.over# 0:ncol(table.DE.time.DE.gene)
    # Colors
    RainBowCol<-as.character(rev(grDevices::rainbow(100)))[-seq_len(32)]
    SeqCol<-floor(seq(1, length(RainBowCol),
                      length.out=ncol(table.DE.time.DE.gene)+1))#length(id.over)
    colors<-RainBowCol[SeqCol]
    #-------------------------------------------------------------------------#
    # Queries for UpSetR::upset
    while(length(vectorUniqueValue)>0){
      n.over<-n.over+1
      mylist.upsetg[[n.over]]<-list(query=UpSetR::elements,
                                    color=colors[n.over], active=TRUE,
                                    params=list("Nb.over",
                                                as.character(vectorUniqueValue)),
                                    query.name=paste(as.character(id.over)[n.over],
                                                     ".ti.over.t0",sep=""))
      vectorUniqueValue<-vectorUniqueValue[-1]
    }# while(length(vectorUniqueValue)>0)
  }else{
    DE.temporal.pattern.f
  }# if(is.null(Log2.FC.matrix)==FALSE)
  #---------------------------------------------------------------------------#
  # 3) Graphs
  #---------------------------------------------------------------------------#
  Ti.Sum.sup0<-which(apply(as.data.frame(table.DE.time),2,sum)>0)
  #
  g.upset<-UpSetR::upset(as.data.frame(table.DE.time), order.by="freq",
                         sets=colnames(table.DE.time)[Ti.Sum.sup0],
                         keep.order=TRUE,
                         number.angles=20, mb.ratio=c(0.7, 0.3),
                         sets.bar.color = "#56B4E9")
  # nsets = ncol(table.DE.time) # print(g.upset)
  if(is.null(Log2.FC.matrix)==FALSE){
    g.upset.over<-UpSetR::upset(OU.pat.f, order.by="freq",
                                sets=colnames(table.DE.time)[Ti.Sum.sup0],
                                keep.order=TRUE,
                                queries=mylist.upsetg, query.legend="top",
                                number.angles=20, mb.ratio=c(0.7, 0.3))
    # nsets = ncol(table.DE.time),
  }# if(is.null(Log2.FC.matrix)==FALSE)
  #---------------------------------------------------------------------------#
  # 4) END
  #---------------------------------------------------------------------------#
  return(list(DE.pattern.t.01.sum=DE.temporal.pattern.f,
              Upset.graph=g.upset,
              Upset.graph.with.nb.over=g.upset.over))
}# DEplotVennBarplotTime()
