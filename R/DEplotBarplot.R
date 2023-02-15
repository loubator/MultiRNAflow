#' @title Barplot of DE genes from a contingency table.
#'
#' @description From a contingency table between two variables,
#' the function plots a barplot of the frequency distribution of one variable
#' against the other (see \code{Details}).
#'
#' @details
#' A contingency table (or cross-tabulation) is a table that displays the
#' frequency distribution of two variables (each containing several levels),
#' i.e. the number of observation recorded per pair of levels.
#' The function plots a single barplot from \code{ContingencyTable}.
#'
#' This function is called by [DEanalysisGroup()] and
#' [DEanalysisTimeAndGroup()].
#' These two functions produce several contingency tables,
#' giving information about specific and particular
#' DE genes, as described below.
#'
#' First, we look for all genes that are DE between at least two biological
#' conditions.
#' A gene will be called specific to a given biological condition BC1,
#' if the gene is DE between BC1 and any other biological conditions,
#' but not DE between any pair of other biological conditions.
#' Then each DE gene will be categorized as follow:
#' * If a gene is not specific, the gene will be categorized as 'Other'.
#' The category 'Other' does not exist when there are only two biological
#' conditions.
#' * If a gene is specific to a given biological condition BC1 and expressions
#' in BC1 are higher than in the other biological conditions,
#' the gene will be categorized as 'Upregulated'.
#' * If a gene is specific to a given biological condition BC1 and expressions
#' in BC1 are lower than in the other biological conditions,
#' the gene will be categorized as 'Downregulated'.
#'
#' The functions [DEanalysisGroup()] and [DEanalysisTimeAndGroup()] produce two
#' contingency table that allow to plot both
#' * the number of genes categorized as 'Other', 'Upregulated'
#' and 'Downregulated'
#' (only when there are strictly more than two biological conditions).
#' * the number of genes categorized 'Upregulated' and 'Downregulated'.
#'
#' Second, we look for all genes that are DE between at least one time point
#' (except t0) and t0 for each biological condition.
#' A gene will be categorized as 'particular' to a given biological condition
#' BC1 for a given time point ti (except t0), if the gene is DE between ti and
#' t0 for the biological condition BC1, but not DE between ti and t0 for
#' the other biological conditions.
#' A gene will be categorized as 'common' to all biological conditions,
#' if the gene is DE between ti and t0 for all biological conditions.
#' Otherwise, a gene will categorized as 'Other'.
#'
#' The function [DEanalysisTimeAndGroup()] produces a contingency table that
#' allow to plot the number of 'specific', 'common' and 'other' genes for
#' each ti (except t0).
#'
#' @param ContingencyTable A numeric data.frame,
#' corresponding to a contingency table, of dimension N1*N2, with N1 and N2,
#' respectively the number of levels in the first and second variable
#' (see examples and details).
#' @param dodge \code{TRUE} or \code{FALSE}.
#' \code{FALSE} means multiple bars in the barplot
#' (one per level of the first variable)
#' one for each fixed level of the other variable.
#' \code{TRUE} means multiple bars will be dodged side-to-side
#' (see [ggplot2::geom_bar()]).
#'
#' @return A barplot using [ggplot2] (see details).
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_bar position_dodge xlab ylab guide_axis
#'
#' @seealso The [DEplotBarplot()] function
#' * is used by the following functions of our package: [DEanalysisGroup()]
#' and [DEanalysisTimeAndGroup()].
#' * calls the R package [ggplot2] in order to plot the barplot.
#'
#' @export
#'
#' @examples
#' # Data simulation
#' CrossTabulation<-matrix(c(75,30,10,5, 5,35,5,20, 220,235,285,275),
#'                         ncol=4, byrow=TRUE)
#' colnames(CrossTabulation)=c("A","B","C","D")
#' row.names(CrossTabulation)=c("Spe.Pos", "Spe.Neg", "Other")
#' #--------------------------------------------------------------------------#
#' res.dodgeTRUE=DEplotBarplot(ContingencyTable=CrossTabulation,dodge=FALSE)
#' res.dodgeTRUE
#' #
#' res.dodgeFALSE=DEplotBarplot(ContingencyTable=CrossTabulation,dodge=TRUE)
#' res.dodgeFALSE

DEplotBarplot<-function(ContingencyTable,
                        dodge=TRUE){
  #---------------------------------------------------------------------------#
  # Data preprocessing for graph
  data.plot.i<-data.frame(Attribute=row.names(ContingencyTable),
                          ContingencyTable)
  data.plot.f<-reshape2::melt(data.plot.i, id.vars=c("Attribute"))
  #---------------------------------------------------------------------------#
  # Graph preprocessing
  ColLegend<-data.plot.f$Attribute
  lev1<-levels(factor(ColLegend))[1]
  levLast<-levels(factor(ColLegend))[-1]
  data.plot.f$Attribute<-factor(data.plot.f$Attribute, levels=c(levLast,lev1))
  # To avoid "no visible binding for global variable" with devtools::check()
  variable<-value<-Attribute<-NULL
  #---------------------------------------------------------------------------#
  # Graph
  if(dodge==TRUE){
    q.cont.g.t<-ggplot2::ggplot(data=data.plot.f,
                                ggplot2::aes(x=variable, y=value,
                                             fill=Attribute)) +
      ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge(),
                        color="black")
    #
  }else{
    q.cont.g.t<-ggplot2::ggplot(data=data.plot.f,
                                ggplot2::aes(x=variable, y=value,
                                             fill=Attribute)) +
      ggplot2::geom_bar(stat="identity", color="black")
  }# if(dodge==TRUE)
  #
  q.cont.g.t<-q.cont.g.t+
    ggplot2::ylab("Number of Genes") + ggplot2::xlab("")+
    ggplot2::scale_x_discrete(guide=ggplot2::guide_axis(angle = 45))
  #
  if(length(unique(data.plot.f$Attribute))==3){
    q.cont.g.t<- q.cont.g.t +
      ggplot2::scale_fill_manual(values=c("#999999", "#E41A1C", "steelblue"))
  }# if(length(unique(data.plot.f$Attribute))==3)
  #
  if(length(unique(data.plot.f$Attribute))==2){
    q.cont.g.t<- q.cont.g.t +
      ggplot2::scale_fill_manual(values=c("#E41A1C", "steelblue"))
  }# if(length(unique(data.plot.f$Attribute))==2)
  #---------------------------------------------------------------------------#
  return(graph.cont.g.t=q.cont.g.t)
}# DEplotBarplot()
