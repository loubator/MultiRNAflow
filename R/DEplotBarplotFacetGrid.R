#' @title Faceted barplot of specific DE genes
#'
#' @description The function creates a faceted barplot from a data.frame
#' containing two or three qualitative variables and one quantitative variable.
#'
#' @param Data Data.frame containing three or four columns.
#' One must contain quantitative variable and the other qualitative variables.
#' @param Abs.col Integer indicating the column of \code{Data} which will be
#' used for the x-axis.
#' The selected column must be one of the qualitative variables and must be
#' identical to \code{Legend.col} if there are only two qualitative variables.
#' Otherwise, \code{Abs.col} and \code{Legend.col} must be different.
#' @param Legend.col Integer indicating the column of \code{Data} which is used
#' for the color of the barplots.
#' The selected column must be one of the qualitative variables and must be
#' identical to \code{Abs.col} if there are only two qualitative variables.
#' Otherwise, \code{Abs.col} and \code{Legend.col} must be different.
#' @param Facet.col Integer indicating the column of \code{Data} which is used
#' for separating barplots in different panels, one per level of the
#' qualitative variable.
#' The selected column must be one of the qualitative variables.
#' @param Value.col Integer indicating the column of \code{Data} which contains
#' numeric values.
#' @param Color.Legend Data.frame or \code{NULL}.
#' If \code{Color.Legend} is a data.frame, the data.frame must have two columns
#' and \eqn{N_{bc}} rows where \eqn{N_{bc}} is the number of biological
#' conditions. The first column must contain the name of the \eqn{N_{bc}}
#' different biological conditions and the second column must the color
#' associated to each biological condition.
#' If \code{Color.Legend=NULL}, the function will automatically attribute
#' a color for each biological condition.
#'
#' @return The function will plot a facet grid barplot.
#' The function is called by our function [DEanalysisTimeAndGroup()]
#' in order to plot the number of specific (up- or down-regulated) DE genes
#' per biological condition for each time points.
#'
#' @seealso The function
#' * is called by the function [DEanalysisTimeAndGroup()]
#' * calls the R functions [ggplot2::facet_grid()] and [ggplot2::geom_bar()].
#'
#' @importFrom ggplot2 ggplot aes_string facet_grid geom_bar element_blank
#' theme element_rect element_text
#' @importFrom stats as.formula
#'
#' @export
#'
#' @examples
#' Group.ex=c('G1', 'G2','G3')
#' Time.ex=c('t1', 't2','t3', 't4')
#' Spe.sign.ex=c("Pos","Neg")
#' Nb.Spe=sample(3:60, length(Group.ex)*length(Time.ex), replace=FALSE)
#' Nb.Spe.sign=sample(3:60, length(Group.ex)*length(Time.ex)*2, replace=FALSE)
#' #--------------------------------------------------------------------------#
#' Melt.Dat.1=data.frame(Group=rep(Group.ex,times=length(Time.ex)),
#'                       Time=rep(Time.ex,each=length(Group.ex)),
#'                       Nb.Spe.DE=Nb.Spe)
#' #
#' DEplotBarplotFacetGrid(Data=Melt.Dat.1,Abs.col=2,Legend.col=2,
#'                        Facet.col=1,Value.col=3,
#'                        Color.Legend=NULL)
#' DEplotBarplotFacetGrid(Data=Melt.Dat.1,Abs.col=1,Legend.col=1,
#'                        Facet.col=2,Value.col=3,
#'                        Color.Legend=NULL)
#' #--------------------------------------------------------------------------#
#' Melt.Dat.2=data.frame(Group=rep(Group.ex,times=length(Time.ex)*2),
#'                       Time=rep(Time.ex,each=length(Group.ex)*2),
#'                       Spe.sign=rep(Spe.sign.ex,
#'                                    times=length(Time.ex)*length(Group.ex)*2),
#'                       Nb.Spe.DE=Nb.Spe.sign)
#'                       #
#' DEplotBarplotFacetGrid(Data=Melt.Dat.2,
#'                        Abs.col=1,
#'                        Legend.col=3,
#'                        Facet.col=2,
#'                        Value.col=4,
#'                        Color.Legend=NULL)

DEplotBarplotFacetGrid<-function(Data,
                                 Abs.col,
                                 Legend.col,
                                 Facet.col,
                                 Value.col,
                                 Color.Legend=NULL){
  #---------------------------------------------------------------------------#
  # Data preprocessing for graph if 'Abs.col!=Legend.col'
  if(Abs.col!=Legend.col){
    if(is.null(Color.Legend)==FALSE){
      Data[,Legend.col]<-factor(Data[,Legend.col], levels=Color.Legend[,1])
    }else{
      Data[,Legend.col]<-factor(Data[,Legend.col])
    }# if(is.null(Color.Legend)==FALSE)
  }# if(Abs.col!=Legend.col)
  #---------------------------------------------------------------------------#
  # Graph
  q.dodged<-ggplot2::ggplot(Data, fill=Data[,Legend.col],
                            ggplot2::aes_string(colnames(Data)[Abs.col],
                                                colnames(Data)[Value.col]))+
    ggplot2::facet_grid(stats::as.formula(paste(". ~",
                                                colnames(Data)[Facet.col])))+
    ggplot2::xlab("") + ggplot2::ylab("Number of genes") +
    ggplot2::theme(panel.background=ggplot2::element_rect(colour="dark grey"))+
    ggplot2::theme(strip.text.x=ggplot2::element_text(size=22, face="bold"),
                   strip.background=ggplot2::element_rect(colour="dark grey"))
  # ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
  #                axis.text.x = ggplot2::element_text(size=5)) +
  if(Abs.col!=Legend.col){
    q.dodged<-q.dodged+
      ggplot2::geom_bar(ggplot2::aes_string(fill=colnames(Data)[Legend.col]),
                        stat="identity", color="black")+ #position = "dodge",
      ggplot2::scale_x_discrete(guide=ggplot2::guide_axis(angle=45))
    #
    if(is.null(Color.Legend)==FALSE){
      q.dodged<-q.dodged+
        ggplot2::scale_fill_manual(values=as.character(Color.Legend[,2]))
    }# if(is.null(Color.Legend)==FALSE)
  }else{
    q.dodged<-q.dodged+
      ggplot2::geom_bar(fill="#E69F00", color="black", stat="identity")+
      ggplot2::scale_x_discrete(guide=ggplot2::guide_axis(angle=90))
  }# if(Abs.col!=Legend.col)
  #---------------------------------------------------------------------------#
  # Output
  return(Barplot.dodged.G.T=q.dodged)
}# DEplotBarplotFacetGrid()
