#' @title Visualization of the distribution of all gene expressions using
#' a boxplot for each sample.
#'
#' @description From the results of either our R function
#' [DATAprepSE()]
#' or our R function
#' [DATAnormalization()]
#' (raw counts or normalized raw counts),
#' the function plots the distribution of all gene expressions using
#' a boxplot for each sample.
#'
#' @details The boxplot allows to visualize six summary statistics
#' (see [ggplot2::geom_boxplot()]):
#' * the median
#' * two hinges: first and third quartiles denoted Q1 and Q3.
#' * two whiskers: \eqn{W1:=Q1-1.5*IQR} and \eqn{W3:=Q3+1.5*IQR}
#' with \eqn{IQR=Q3-Q1}, the interquartile range.
#' * outliers: data beyond the end of the whiskers are called "outlying"
#' points and are plotted in black.
#'
#' For better visualization of the six summary statistics described above,
#' raw counts must be transformed using the function \eqn{log_2(x+1)}.
#' This transformation is automatically performed by other functions of
#' the package, such as
#' [DATAnormalization()].
#' \code{Log2.transformation} will be set as TRUE in
#' [DATAnormalization()]
#' if \code{Normalization ="rle"}, otherwise \code{Log2.transformation=FALSE}.
#'
#' @param SEres Results of the function
#' [DATAprepSE()] or
#' [DATAnormalization()].
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
#' If \code{Color.Group=NULL}, the function will automatically attribute a
#' color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param Plot.genes \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, points representing gene expression
#' (normalized or raw counts) will be added for each sample.
#' @param y.label \code{NULL} or a character. \code{NULL} as default.
#' If \code{y.label} is a character, it will be the y label of the graph.
#' If \code{y.label=NULL}, the label will be either "log2(Gene expression +1)"
#' (if \code{Log2.transformation=TRUE}) or "Gene expression"
#' (if \code{Log2.transformation=FALSE}).
#'
#' @return The function returns a graph which plots the distribution of all
#' gene expressions using a boxplot for each sample
#' (see [ggplot2::geom_boxplot()]).
#'
#' @seealso The [DATAplotBoxplotSamples()] function
#' * is used by the following function of our package:
#' [DATAnormalization()].
#' * calls the R functions
#' [ggplot2::geom_boxplot] and
#' [ggplot2::geom_jitter]
#' in order to print the boxplot.
#'
#' @importFrom SummarizedExperiment assays colData
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
#' ##-------------------------------------------------------------------------#
#' resDATAprepSE <- DATAprepSE(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=NULL,
#'                             Individual.position=2)
#' ##-------------------------------------------------------------------------#
#' DATAplotBoxplotSamples(SEres=resDATAprepSE,
#'                        Log2.transformation=TRUE,
#'                        Colored.By.Factors=TRUE,
#'                        Color.Group=NULL,
#'                        Plot.genes=FALSE,
#'                        y.label=NULL)

DATAplotBoxplotSamples <- function(SEres,
                                   Log2.transformation=TRUE,
                                   Colored.By.Factors=FALSE,
                                   Color.Group=NULL,
                                   Plot.genes=FALSE,
                                   y.label=NULL) {
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check
    ## DATAprepSE
    Err_SE <- paste0("'SEres' mut be the results of either the function ",
                     "'DATAprepSE()' or 'DATAnormalization()'.")
    if (is.null(SEres$SEidentification)) {
        stop(Err_SE)
    } else {
        if (!SEres$SEidentification%in%c("SEstep", "SEresNormalization")) {
            stop(Err_SE)
        }
    }## if (SEres$SEidentification!="SEstep")

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## To avoid "no visible binding for global variable" with devtools::check()
    BioCond <- Time <- NULL

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Preprocessing

    Ndata <- length(SummarizedExperiment::assays(SEres$SEobj))
    ExprData <- SummarizedExperiment::assays(SEres$SEobj)[[Ndata]]

    cDat <- data.frame(SummarizedExperiment::colData(SEres$SEobj))

    ## which(colnames(cDat)%in%c("Group")), SEobj
    if (c("Group")%in%colnames(cDat)) {
        SEinfoGroup <- as.character(cDat$Group)
    } else {
        SEinfoGroup <- NULL
    }## if (c("Group")%in%colnames(cDat))

    if (c("Time")%in%colnames(cDat)) {
        SEinfoTime <- as.character(cDat$Time)
    } else {
        SEinfoTime <- NULL
    }## if (c("Time")%in%colnames(cDat))

    SEfinalname <- SEres$Data.code.names$Final.Name
    N.spl <- length(SEfinalname)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Data reshaped
    ExprData.f <- data.frame(Gene=as.character(seq_len(nrow(ExprData))),
                             ExprData)
    colnames(ExprData.f) <- c("Gene", SEfinalname)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Data used for boxplot
    Norm.dat.melt <- reshape2::melt(ExprData.f,
                                    id=c("Gene"),
                                    value.name="Expr",
                                    variable.name="Samples")
    colnames(Norm.dat.melt) <- c("Gene", "Samples", "Expression")
    Norm.dat.melt$Samples <- as.factor(Norm.dat.melt$Samples)

    if (isTRUE(Log2.transformation)) {
        Norm.dat.melt$Expression <- log2(Norm.dat.melt$Expression + 1)
        ylab.epr <- "log2(Gene expression +1)"
    } else {
        ylab.epr <- "Gene expression"
    }## if(isTRUE(Log2.transformation))

    if (!is.null(y.label)) {
        ylab.epr <- y.label
    }## if(!is.null(y.label))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Graph
    Samples <- Expression <- NULL

    if (isFALSE(Colored.By.Factors)) {
        if (isFALSE(Plot.genes)) {
            res.bxplt <- ggplot2::ggplot(Norm.dat.melt,
                                         ggplot2::aes(x=Samples,
                                                      y=Expression)) +
                ggplot2::geom_boxplot(alpha=0.8,
                                      show.legend=FALSE, outlier.alpha=0.2,
                                      color="black", fill="#E69F00") +
                ggplot2::scale_x_discrete(guide=ggplot2::guide_axis(angle=90))
        } else {
            Width.plt <- min(10/N.spl, 0.3)

            res.bxplt <- ggplot2::ggplot(Norm.dat.melt,
                                         ggplot2::aes(x=Samples,
                                                      y=Expression)) +
                ggplot2::geom_jitter(position=ggplot2::position_jitter(
                    width=Width.plt, height=0.001),
                    alpha=0.7, fill="#56B4E9", color="#56B4E9") +
                ggplot2::geom_boxplot(alpha=0.8,
                                      show.legend=FALSE, outlier.alpha=0.2,
                                      color="black", fill="#E69F00") +
                ggplot2::scale_x_discrete(guide=ggplot2::guide_axis(angle=90))
        }## if(isFALSE(Plot.genes))

        res.bxplt <- res.bxplt + ggplot2::labs(x="Samples", y=ylab.epr)

    } else {
        if (!is.null(SEinfoGroup)) {
            ##----------------------------------------------------------------#
            ## Setting
            FactorBoxplG <- as.character(SEinfoGroup)
            SEQlenFCT <- seq_len(length(FactorBoxplG))
            Num.FactorBC <- CharacterNumbers(Vect.number=SEQlenFCT)

            Expr.colG <- ExprData.f
            colnames(Expr.colG)<-c("Gene", paste0(Num.FactorBC, FactorBoxplG))

            Max.digit.G<-floor(log10(abs(max(SEQlenFCT))))
            Max.digit.G <- Max.digit.G + 1

            Melt1F.G <- reshape2::melt(Expr.colG, id=c("Gene"),
                                       value.name="Expr",
                                       variable.name="BioCond")
            substrg.BC <- substring(Melt1F.G$BioCond,
                                    first=Max.digit.G + 1,
                                    last=1000000L)
            Dat.bxplt.G <- data.frame(Norm.dat.melt[, -1],
                                      substrg.BC)
            colnames(Dat.bxplt.G) <- c("Samples", "Expression", "BioCond")

            ##----------------------------------------------------------------#
            Glevels <- levels(factor(SEinfoGroup))

            if (is.null(Color.Group)) {
                MypaletteG <- c(RColorBrewer::brewer.pal(8, "Dark2"),
                                RColorBrewer::brewer.pal(8, "Set2"))
                if (length(Glevels) > 16) {
                    MypaletteG <- c(MypaletteG,
                                    scales::hue_pal(l=90)(
                                        seq_len(length(Glevels) - 1)))
                }## if(length(Glevels)>16)

                Color.Group <- data.frame(Name=Glevels,
                                          Col=MypaletteG[seq_len(
                                              length(Glevels))])
            } else {
                Id.LevelCol.G <- order(Color.Group[, 1])
                Color.Group <- data.frame(Name=Glevels,
                                          Col=Color.Group[Id.LevelCol.G, 2])
            }## if(is.null(Color.Group))

            VcolG <- factor(SEinfoGroup)
            levels(VcolG) <- Color.Group$Col

            LegendTitle <- "Biological \nconditions"

            ##----------------------------------------------------------------#
            if (isFALSE(Plot.genes)) {
                res.bxplt <- ggplot2::ggplot(Dat.bxplt.G,
                                             ggplot2::aes(x=factor(Samples),
                                                          y=Expression,
                                                          fill=BioCond)) +
                    ggplot2::geom_boxplot(alpha=1, outlier.alpha=0.2) +
                    ggplot2::theme(strip.text.x=ggplot2::element_text(
                        size=9, color="black", face="bold")) +
                    ggplot2::labs(x="Samples", y=ylab.epr)+
                    ggplot2::scale_x_discrete(
                        guide=ggplot2::guide_axis(angle=90)) +
                    ggplot2::guides(fill=ggplot2::guide_legend(LegendTitle)) +
                    ggplot2::scale_fill_manual(values=levels(VcolG))
            } else {
                JitterCol <- factor(Dat.bxplt.G[, 3])
                levels(JitterCol) <- levels(VcolG)

                res.bxplt <- ggplot2::ggplot(Dat.bxplt.G,
                                             ggplot2::aes(x=factor(Samples),
                                                          y=Expression,
                                                          fill=BioCond)) +
                    ggplot2::geom_jitter(
                        position=ggplot2::position_jitter(width=0.3,
                                                          height=0.2),
                        colour=JitterCol,
                        alpha=0.9) +
                    ggplot2::geom_boxplot(alpha=1, outlier.alpha=0.2) +
                    ggplot2::theme(
                        strip.text.x=ggplot2::element_text(size=9,
                                                           color="black",
                                                           face="bold")) +
                    ggplot2::labs(x="Samples", y=ylab.epr)+
                    ggplot2::scale_x_discrete(
                        guide=ggplot2::guide_axis(angle=90))+
                    ggplot2::guides(fill=ggplot2::guide_legend(LegendTitle))+
                    ggplot2::scale_fill_manual(values=levels(VcolG))
            }# if(Plot.genes==FALSE)
        }# if(is.null(SEinfoGroup)==FALSE)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        if (is.null(SEinfoGroup) & !is.null(SEinfoTime)) {
            ##----------------------------------------------------------------#
            Expr.colT <- ExprData.f
            FactorBoxplT <- as.character(SEinfoTime)
            Num.FactorT <- CharacterNumbers(
                Vect.number=seq_len(length(FactorBoxplT)))
            colnames(Expr.colT) <- c("Gene", paste0(Num.FactorT, FactorBoxplT))
            Max.digit.G<-floor(log10(abs(max(seq_len(length(FactorBoxplT))))))
            Max.digit.G <- Max.digit.G + 1

            ##----------------------------------------------------------------#
            Melt1F.T <- reshape2::melt(Expr.colT, id=c("Gene"),
                                       value.name="Expr", variable.name="Time")
            substrg.T <- substring(Melt1F.T$Time, first=Max.digit.G + 1,
                                   last=1000000L)
            substrg.T2 <- gsub("T", "", gsub("t", "", as.character(substrg.T)))
            substrg.Tf <- paste0("t", substrg.T2)
            Dat.bxplt.T <- data.frame(Norm.dat.melt[, -1], substrg.Tf)
            colnames(Dat.bxplt.T)<-c("Samples", "Expression", "Time")

            ##----------------------------------------------------------------#
            Tlevels <- levels(factor(substrg.Tf))
            NbTime <- length(Tlevels)
            Color.Time <- NULL

            if (is.null(Color.Time)) {
                Color.Time <- data.frame(Name=Tlevels,
                                         Col=c("#737373",# "#252525"
                                               scales::hue_pal()(
                                                   length(Tlevels)-1)))
            } else {
                Id.LevelColT <- order(Color.Time[, 1])
                Color.Time <- data.frame(Name=Tlevels,
                                         Col=Color.Time[Id.LevelColT, 2])
            }## if(is.null(Color.Time))

            VcolT <- paste0("t",
                            gsub("T", "",
                                 gsub("t", "", as.character(SEinfoTime))))
            VcolT <- factor(VcolT)
            levels(VcolT) <- Color.Time$Col

            ##----------------------------------------------------------------#
            if (isFALSE(Plot.genes)) {
                res.bxplt <- ggplot2::ggplot(Dat.bxplt.T,
                                             ggplot2::aes(x=factor(Samples),
                                                          y=Expression,
                                                          fill=Time)) +
                    ggplot2::geom_boxplot(alpha=1,
                                          outlier.alpha=0.2) +
                    ggplot2::theme(
                        strip.text.x=ggplot2::element_text(size=9,
                                                           color="black",
                                                           face="bold")) +
                    ggplot2::labs(x="Samples", y=ylab.epr)+
                    ggplot2::scale_x_discrete(
                        guide=ggplot2::guide_axis(angle=90))+
                    ggplot2::guides(color=ggplot2::guide_legend("Time"))+
                    ggplot2::scale_fill_manual(values=levels(VcolT))
            } else {
                JitterCol <- factor(Dat.bxplt.T[,3])
                levels(JitterCol) <- levels(VcolT)

                res.bxplt <- ggplot2::ggplot(Dat.bxplt.T,
                                             ggplot2::aes(x=factor(Samples),
                                                          y=Expression,
                                                          fill=Time)) +
                    ggplot2::geom_jitter(
                        position=ggplot2::position_jitter(width=0.3,
                                                          height=0.2),
                        colour=JitterCol,
                        alpha=0.9) +
                    ggplot2::geom_boxplot(alpha=1, outlier.alpha=0.2) +
                    ggplot2::theme(
                        strip.text.x=ggplot2::element_text(size=9,
                                                           color="black",
                                                           face="bold")) +
                    ggplot2::labs(x="Samples", y=ylab.epr)+
                    ggplot2::scale_x_discrete(
                        guide=ggplot2::guide_axis(angle=90))+
                    ggplot2::guides(color=ggplot2::guide_legend("Time"))+
                    ggplot2::scale_fill_manual(values=levels(VcolT))
            }## if(isFALSE(Plot.genes))
        }##if(is.null(SEinfoGroup)&!is.null(SEinfoTime))
    }## if(isFALSE(Colored.By.Factors))

    ##------------------------------------------------------------------------#
    ## Final ggplot2
    res.bxplt <- res.bxplt +
        ggplot2::ylim(min=min(Norm.dat.melt$Expression),
                      max=max(Norm.dat.melt$Expression) + 0.1)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(res.bxplt=res.bxplt)
}## DATAplotBoxplotSamples()
