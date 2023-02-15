#' @title Venn barplot of DE genes across pairs of biological conditions.
#'
#' @description The function takes as input a binary matrix or data.frame with
#' \eqn{N_g} rows and \eqn{((N_{bc}-1)\times N_{bc})/2} columns with
#' \eqn{N_g} the number of genes and \eqn{N_{bc}}
#' the number of biological conditions.
#' The number of 1 in the n-th row gives the number of pairs of
#' biological conditions where the gene \eqn{n} is DE.
#' We consider that a set of pairs of biological conditions forms
#' an intersection if there is at least one gene which is DE for each of
#' these pairs of biological conditions, but not for the others.
#'
#' The function calls the [UpSetR::upset()] function in order to plot
#' the number of genes for each possible intersection in an UpSet plot
#' (Venn diagram displayed as a barplot).
#'
#' @param Mat.DE.pair.group Binary matrix or data.frame with \eqn{N_g} rows and
#' \eqn{((N_{bc}-1)*N_{bc})/2} columns with
#' \eqn{N_{bc}} the number of biological conditions.
#'
#' @return The function plots the number of genes for each possible
#' intersection in a UpSet plot.
#'
#' @seealso The function
#' * calls the function [UpSetR::upset()] in order to plot the UpSet plot.
#' * is called by the functions [DEanalysisGroup()] and
#' [DEanalysisTimeAndGroup()].
#'
#' @importFrom UpSetR upset
#' @importFrom plyr count
#'
#' @export
#'
#' @examples
#' set.seed(1994)
#' #--------------------------------------------------------------------------#
#' # Binary matrix
#' Bin.Table.G=matrix(c(sample(c(0,1),replace=TRUE,size=240,c(0.75,0.35)),
#'                      sample(c(0,1),replace=TRUE,size=240,c(0.3,0.7)),
#'                      rep(0,18)),
#'                    ncol=6, byrow = TRUE)
#' colnames(Bin.Table.G)=c(".A..B.",".A..C.",".A..D.",
#'                         ".B..C.",".B..D.",".C..D.")
#' #--------------------------------------------------------------------------#
#' # Results
#' res.t.upset=DEplotVennBarplotGroup(Mat.DE.pair.group=Bin.Table.G)
#' print(res.t.upset$Upset.global)
#' print(res.t.upset$Upset.threshold)

DEplotVennBarplotGroup<-function(Mat.DE.pair.group){
  #---------------------------------------------------------------------------#
  # 1) Data for graphs
  #---------------------------------------------------------------------------#
  Nb.pair.group<-ncol(Mat.DE.pair.group)
  Mat.DE.pair.group<-as.data.frame(Mat.DE.pair.group)
  #
  Bin.mat.DE.f<-Mat.DE.pair.group[which(apply(Mat.DE.pair.group,
                                              MARGIN=1,
                                              FUN=sum)>0),]
  row.names(Bin.mat.DE.f)<-NULL
  #
  Freq.pat<-plyr::count(Bin.mat.DE.f, vars=colnames(Bin.mat.DE.f))
  Freq.unique.pat<-nrow(plyr::count(Mat.DE.pair.group,
                                    vars=colnames(Mat.DE.pair.group)))
  s.upset<-nrow(Bin.mat.DE.f)/(Freq.unique.pat-1)
  # 0,...,0 excluded (Freq.unique.pat*2)
  Freq.pat.del<-Freq.pat[which(Freq.pat$freq<s.upset),-(Nb.pair.group+1)]
  Pattern.del<-as.character(apply(Freq.pat.del,
                                  MARGIN=1,
                                  FUN=function(x) paste(x,collapse="")))
  All.patern<-as.character(apply(Bin.mat.DE.f,
                                 MARGIN=1,
                                 FUN=function(x) paste(x,collapse="")))
  #---------------------------------------------------------------------------#
  # 2) graphs
  #---------------------------------------------------------------------------#
  if(Nb.pair.group>1){
    g.upset1<-UpSetR::upset(Bin.mat.DE.f,
                            nsets=Nb.pair.group, number.angles=20,
                            order.by="freq", keep.order=TRUE,
                            mb.ratio=c(0.7, 0.3),
                            sets.bar.color="#56B4E9")
    # print(g.upset1)
    #-------------------------------------------------------------------------#
    if(length(Pattern.del)>0){
      g.upset2<-UpSetR::upset(Bin.mat.DE.f[-which(All.patern%in%Pattern.del),],
                              nsets=Nb.pair.group, number.angles=20,
                              order.by="freq", keep.order=TRUE,
                              mb.ratio=c(0.7, 0.3),
                              sets.bar.color="#56B4E9")
      # print(g.upset2)
    }else{
      g.upset2<-NULL
    }# if(length(Pattern.del)>0)
  }else{
    g.upset1<-NULL
    g.upset2<-NULL
  }# if(Nb.pair.group>1)
  # output
  return(list(Upset.global=g.upset1,
              Upset.threshold=g.upset2))
}# DEplotVennBarplotGroup
