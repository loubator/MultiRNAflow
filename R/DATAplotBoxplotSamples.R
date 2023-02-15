#' @title Visualization of the distribution of all gene expressions using
#' a boxplot for each sample.
#'
#' @description From a gene expression dataset
#' (raw counts or normalized raw counts),
#' the function plots the distribution of all gene expressions using
#' a boxplot for each sample.
#'
#' @details The column names of \code{ExprData} must be a vector of strings
#' of characters containing
#' * a string of characters (if \eqn{k=1}) which is the label of the column
#' containing gene names.
#' * \eqn{N_s} sample names which must be strings of characters containing
#' at least: the name of the individual (e.g patient, mouse, yeasts culture),
#' its biological condition (if there is at least two) and
#' the time where data have been collected if there is at least two;
#' (must be either 't0', 'T0' or '0' for time 0,
#' 't1', 'T1' or '1' for time 1, ...).
#'
#' All these sample information must be separated by underscores in
#' the sample name. For instance 'CLL_P_t0_r1',
#' corresponds to the patient 'r1' belonging to the biological condition 'P'
#' and where data were collected at time 't0'.
#' I this example, 'CLL' describe the type of cells
#' (here chronic lymphocytic leukemia) and is not used in our analysis.
#'
#' In the string of characters 'CLL_P_t0_r1',
#' 'r1' is localized after the third underscore, so \code{Individual.position=4},
#' 'P' is localized after the first underscore, so \code{Group.position=2} and
#' 't0' is localized after the second underscore, so \code{Time.position=3}.
#'
#' The boxplot allows to visualize six summary statistics
#' (see [ggplot2::geom_boxplot()]):
#' * the median
#' * two hinges: first and third quartiles denoted Q1 and Q3.
#' * two whiskers: \eqn{W1:=Q1-1.5*IQR} and \eqn{W3:=Q3+1.5*IQR}
#' with \eqn{IQR=Q3-Q1}, the interquartile range.
#' * outliers: data beyond the end of the whiskers are called "outlying" points
#' and are plotted in black.
#'
#' For better visualization of the six summary statistics described above,
#' raw counts must be transformed using the function \eqn{log_2(x+1)}.
#' This transformation is automatically performed by other functions of
#' the package, such as [DATAnormalization()].
#' \code{Log2.transformation} will be set as TRUE in [DATAnormalization()]
#' if \code{Normalization ="rle"}, otherwise \code{Log2.transformation=FALSE}.
#'
#' @param ExprData Data.frame with \eqn{N_g} rows and (\eqn{N_{s+k}}) columns,
#' where \eqn{N_g} is the number of genes,
#' \eqn{N_s} is the number of samples and
#' \eqn{k=1} if a column is used to specify gene names, or \eqn{k=0} otherwise.
#' If \eqn{k=1}, the position of the column containing gene names is given
#' by \code{Column.gene}.
#' The data.frame contains numeric values giving gene expressions of each gene
#' in each sample.
#' Gene expressions can be raw counts or normalized raw counts.
#' Column names of the data.frame must describe each sample's information
#' (individual, biological condition and time) and have the structure described
#' in the section \code{Details}.
#' @param Column.gene Integer indicating the column where gene names are given.
#' Set \code{Column.gene=NULL} if there is no such column.
#' @param Group.position Integer indicating the position of group information
#' in the string of characters in each sample names (see \code{Details}).
#' Set \code{Group.position=NULL} if there is only one or
#' no biological information in the string of character in each sample name.
#' @param Time.position Integer indicating the position of time measurement
#' information in the string of characters in each sample names
#' (see \code{Details}).
#' Set \code{Time.position=NULL} if there is only one or no time measurement
#' information in the string of character in each sample name.
#' @param Individual.position Integer indicating the position of the name of
#' the individual (e.g patient, replicate, mouse, yeasts culture ...)
#' in the string of characters in each sample names (see \code{Details}).
#' The names of different individuals must be all different.
#' Furthermore, if individual names are just numbers, they will be transform in
#' a vector of class "character" by [CharacterNumbers()] and
#' a "r" will be added to each individual name ("r" for replicate).
#' @param Log2.transformation \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, each numeric value \eqn{x} in \code{ExprData} will become
#' \eqn{log_2(x+1)} (see \code{Details}).
#' @param Colored.By.Factors \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, boxplots will be colored with different colors for different
#' time measurements (if data were collected at different time points).
#' Otherwise, boxplots will be colored with different colors for
#' different biological conditions.
#' @param Color.Group \code{NULL} or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute a color
#' for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param Plot.genes \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, points representing gene expression
#' (normalized or raw counts) will be added for each sample.
#' @param y.label \code{NULL} or a character. \code{NULL} as default.
#' If \code{y.label} is a character, it will be the y label of the graph.
#' If \code{y.label=NULL}, the label will be either "log2(Gene expression +1)"
#' (if \code{Log2.transformation=TRUE}) either "Gene expression"
#' (if \code{Log2.transformation=FALSE}).
#'
#' @return The function returns a graph which plots the distribution of all gene
#' expressions using a boxplot for each sample (see [ggplot2::geom_boxplot()]).
#'
#' @seealso The [DATAplotBoxplotSamples()] function
#' * is used by the following function of our package: [DATAnormalization()].
#' * calls the R functions [ggplot2::geom_boxplot] and [ggplot2::geom_jitter]
#' in order to print the boxplot.
#'
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales hue_pal
#' @importFrom ggplot2 ggplot aes geom_boxplot theme labs guides
#' scale_x_discrete guide_axis element_text guide_legend scale_fill_manual
#' geom_jitter position_jitter
#'
#' @export
#'
#' @examples
#' data(RawCounts_Antoszewski2022_MOUSEsub500)
#' DATAplotBoxplotSamples(ExprData=RawCounts_Antoszewski2022_MOUSEsub500,
#'                        Column.gene=1,
#'                        Group.position=1,
#'                        Time.position=NULL,
#'                        Individual.position=2,
#'                        Log2.transformation=TRUE,
#'                        Colored.By.Factors=TRUE,
#'                        Color.Group=NULL,
#'                        Plot.genes=FALSE,
#'                        y.label=NULL)

DATAplotBoxplotSamples<-function(ExprData,
                                 Column.gene,
                                 Group.position,
                                 Time.position,
                                 Individual.position,
                                 Log2.transformation=TRUE,
                                 Colored.By.Factors=FALSE,
                                 Color.Group=NULL,
                                 Plot.genes=FALSE,
                                 y.label=NULL){
  #---------------------------------------------------------------------------#
  # To avoid "no visible binding for global variable" with devtools::check()
  BioCond<-Time<-NULL
  #---------------------------------------------------------------------------#
  # Preprocessing
  res.colnames<-ColnamesToFactors(ExprData=ExprData,
                                  Column.gene=Column.gene,
                                  Group.position=Group.position,
                                  Time.position=Time.position,
                                  Individual.position=Individual.position)
  N.spl<-length(res.colnames$Final.Name)
  #---------------------------------------------------------------------------#
  # columns with only expression
  if(is.null(Column.gene)==TRUE){
    ind.col.expr<-seq_len(ncol(ExprData))
  }else{
    ind.col.expr<-seq_len(ncol(ExprData))[-Column.gene]
  }# if(is.null(Column.gene)==TRUE)
  #---------------------------------------------------------------------------#
  # Data reshaped
  ExprData.f<-cbind(as.character(seq_len(nrow(ExprData))),
                    ExprData[,ind.col.expr])
  colnames(ExprData.f)<-c("Gene",res.colnames$Final.Name)
  #---------------------------------------------------------------------------#
  # Data used for boxplot
  Norm.dat.melt<-reshape2::melt(ExprData.f, id=c("Gene"),
                                value.name="Expr",
                                variable.name="Samples")
  colnames(Norm.dat.melt)<-c("Gene","Samples","Expression")
  Norm.dat.melt$Samples<-as.factor(Norm.dat.melt$Samples)
  #
  if(Log2.transformation==TRUE){
    Norm.dat.melt$Expression<-log2(Norm.dat.melt$Expression+1)
    ylab.epr<-"log2(Gene expression +1)"
  }else{
    ylab.epr<-"Gene expression"
  }# if(Log2.transformation==TRUE)
  #
  if(is.null(y.label)==FALSE){
    ylab.epr<-y.label
  }# if(is.null(y.label)==FALSE)
  #---------------------------------------------------------------------------#
  # Graph
  Samples<-Expression<-NULL
  if(Colored.By.Factors==FALSE){
    if(Plot.genes==FALSE){
      res.bxplt<-ggplot2::ggplot(Norm.dat.melt,
                                 ggplot2::aes(x=Samples, y=Expression)) +
        ggplot2::geom_boxplot(alpha=0.8, show.legend=FALSE, outlier.alpha=0.2,
                              color="black", fill="#E69F00")+
        ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))
    }else{
      Width.plt<-min(10/N.spl,0.3)
      res.bxplt<-ggplot2::ggplot(Norm.dat.melt,
                                 ggplot2::aes(x=Samples, y=Expression)) +
        ggplot2::geom_jitter(position=ggplot2::position_jitter(width=Width.plt,
                                                               height=0.001),
                             alpha=0.7, fill="#56B4E9", color="#56B4E9") +
        ggplot2::geom_boxplot(alpha=0.8, show.legend=FALSE, outlier.alpha=0.2,
                              color="black", fill="#E69F00")+
        ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))
    }# if(Plot.genes==FALSE)
    res.bxplt<-res.bxplt+ggplot2::labs(x="Samples", y=ylab.epr)
  }else{
    if(is.null(res.colnames$Group.Info)==FALSE){
      Expr.colG<-ExprData.f
      FactorBoxplG<-as.character(res.colnames$Group.Info)
      Num.FactorBC<-CharacterNumbers(Vect.number=seq_len(length(FactorBoxplG)))
      colnames(Expr.colG)<-c("Gene", paste(Num.FactorBC, FactorBoxplG, sep=""))
      #
      Max.digit.G<-floor(log10(abs(max(seq_len(length(FactorBoxplG))))))+ 1
      #
      Melt1F.G<-reshape2::melt(Expr.colG, id=c("Gene"),
                               value.name = "Expr",
                               variable.name = "BioCond")
      substrg.BC<-substring(Melt1F.G$BioCond,first=Max.digit.G+1,last=1000000L)
      Dat.bxplt.G<-data.frame(Norm.dat.melt[,-1], substrg.BC)
      colnames(Dat.bxplt.G)<-c("Samples","Expression","BioCond")
      #-----------------------------------------------------------------------#
      Glevels<-levels(factor(res.colnames$Group.Info))
      #
      if(is.null(Color.Group)==TRUE){
        MypaletteG<-c(RColorBrewer::brewer.pal(8,"Dark2"),
                      RColorBrewer::brewer.pal(8,"Set2"))
        if(length(Glevels)>16){
          MypaletteG<-c(MypaletteG,hue_pal(l=90)(seq_len(length(Glevels)-1)))
        }# if(length(Glevels)>16)
        Color.Group<-data.frame(Name=Glevels,
                                Col=MypaletteG[seq_len(length(Glevels))])
      }else{
        Id.LevelCol.G<-order(Color.Group[,1])
        Color.Group<-data.frame(Name=Glevels,
                                Col=Color.Group[Id.LevelCol.G,2])
      }# if(is.null(Color.Group)==TRUE)
      VcolG<-factor(res.colnames$Group.Info)
      levels(VcolG)<-Color.Group$Col
      #
      LegendTitle<-"Biological \nconditions"
      #-----------------------------------------------------------------------#
      if(Plot.genes==FALSE){
        res.bxplt<-ggplot2::ggplot(Dat.bxplt.G,
                                   ggplot2::aes(x=factor(Samples),
                                                y=Expression, fill=BioCond)) +
          ggplot2::geom_boxplot(alpha = 1, outlier.alpha = 0.2) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(size=9,
                                                              color="black",
                                                              face="bold")) +
          ggplot2::labs(x="Samples", y=ylab.epr)+
          ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))+
          ggplot2::guides(fill=ggplot2::guide_legend(LegendTitle))+
          ggplot2::scale_fill_manual(values=levels(VcolG))
      }else{
        JitterCol<-factor(Dat.bxplt.G[,3])
        levels(JitterCol)<-levels(VcolG)
        #
        res.bxplt<-ggplot2::ggplot(Dat.bxplt.G,
                                   ggplot2::aes(x=factor(Samples),
                                                y=Expression,
                                                fill=BioCond)) +
          ggplot2::geom_jitter(position=ggplot2::position_jitter(width=0.3,
                                                                 height=0.2),
                               colour=JitterCol,
                               alpha=0.9) +
          ggplot2::geom_boxplot(alpha = 1, outlier.alpha = 0.2) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(size=9,
                                                              color="black",
                                                              face="bold")) +
          ggplot2::labs(x="Samples", y=ylab.epr)+
          ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))+
          ggplot2::guides(fill=ggplot2::guide_legend(LegendTitle))+
          ggplot2::scale_fill_manual(values=levels(VcolG))
      }# if(Plot.genes==FALSE)
    }# if(is.null(res.colnames$Group.Info)==FALSE)
    #-------------------------------------------------------------------------#
    if(is.null(res.colnames$Group.Info)==TRUE & is.null(res.colnames$Time.Info)==FALSE){
      Expr.colT<-ExprData.f
      FactorBoxplT<-as.character(res.colnames$Time.Info)
      Num.FactorT<-CharacterNumbers(Vect.number=seq_len(length(FactorBoxplT)))
      colnames(Expr.colT)<-c("Gene", paste(Num.FactorT, FactorBoxplT, sep=""))
      Max.digit.G<-floor(log10(abs(max(seq_len(length(FactorBoxplT))))))+ 1
      #
      Melt1F.T<-reshape2::melt(Expr.colT, id=c("Gene"),
                               value.name="Expr", variable.name="Time")
      substrg.T<-substring(Melt1F.T$Time, first=Max.digit.G+1, last=1000000L)
      substrg.T2<-gsub("T","", gsub("t","", as.character(substrg.T)))
      substrg.Tf<-paste("t", substrg.T2, sep="")
      Dat.bxplt.T<-data.frame(Norm.dat.melt[,-1], substrg.Tf)
      colnames(Dat.bxplt.T)<-c("Samples","Expression","Time")
      #-----------------------------------------------------------------------#
      Tlevels<-levels(factor(substrg.Tf))
      NbTime<-length(Tlevels)
      #
      Color.Time<-NULL
      if(is.null(Color.Time)==TRUE){
        Color.Time<-data.frame(Name=Tlevels,
                               Col=c("#737373",# "#252525"
                                     scales::hue_pal()(length(Tlevels)-1)))
      }else{
        Id.LevelColT<-order(Color.Time[,1])
        Color.Time<-data.frame(Name=Tlevels,
                               Col=Color.Time[Id.LevelColT,2])
      }# if(is.null(Color.Time)==TRUE)
      VcolT<- factor(paste("t",gsub("T","",
                                    gsub("t","",
                                         as.character(res.colnames$Time.Info))),
                           sep=""))
      levels(VcolT)<-Color.Time$Col
      #-----------------------------------------------------------------------#
      if(Plot.genes==FALSE){
        res.bxplt<-ggplot2::ggplot(Dat.bxplt.T,
                                   ggplot2::aes(x=factor(Samples),
                                                y=Expression,fill=Time)) +
          ggplot2::geom_boxplot(alpha = 1, outlier.alpha = 0.2) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(size=9,
                                                              color="black",
                                                              face="bold")) +
          ggplot2::labs(x="Samples", y=ylab.epr)+
          ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))+
          ggplot2::guides(color=ggplot2::guide_legend("Time"))+
          ggplot2::scale_fill_manual(values=levels(VcolT))
      }else{
        JitterCol<-factor(Dat.bxplt.T[,3])
        levels(JitterCol)<-levels(VcolT)
        #
        res.bxplt<-ggplot2::ggplot(Dat.bxplt.T,
                                   ggplot2::aes(x=factor(Samples),
                                                y=Expression,fill=Time)) +
          ggplot2::geom_jitter(position=ggplot2::position_jitter(width=0.3,
                                                                 height=0.2),
                               colour=JitterCol,
                               alpha=0.9) +
          ggplot2::geom_boxplot(alpha = 1, outlier.alpha = 0.2) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(size=9,
                                                              color="black",
                                                              face="bold")) +
          ggplot2::labs(x="Samples", y=ylab.epr)+
          ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))+
          ggplot2::guides(color=ggplot2::guide_legend("Time"))+
          ggplot2::scale_fill_manual(values=levels(VcolT))
      }# if(Plot.genes==FALSE)
    }# if(is.null(res.colnames$Group.Info & res.colnames$Time.Info)==TRUE)
  }# if(Colored.By.Factors==FALSE)
  res.bxplt<-res.bxplt + ggplot2::ylim(min=min(Norm.dat.melt$Expression),
                                       max=max(Norm.dat.melt$Expression)+0.1)
  return(res.bxplt=res.bxplt)
}# DATAplotBoxplotSamples
