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
#' \code{TRUE} by default.
#' If \code{TRUE}, each numeric value \eqn{x} in \code{ExprData} will become
#' \eqn{log_2(x+1)} (see \code{Details}).
#' @param Colored.By.Factors \code{TRUE} or \code{FALSE}.
#' \code{FALSE} by default.
#' If \code{TRUE}, boxplots will be colored with different colors for different
#' time measurements (if data were collected at different time points).
#' Otherwise, boxplots will be colored with different colors for
#' different biological conditions.
#' @param Color.Group \code{NULL} or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' \code{NULL} by default.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute a
#' color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param Plot.genes \code{TRUE} or \code{FALSE}. \code{FALSE} by default.
#' If \code{TRUE}, points representing gene expression
#' (normalized or raw counts) will be added for each sample.
#' @param y.label \code{NULL} or a character. \code{NULL} by default.
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
#' @importFrom ggplot2 ggplot aes geom_boxplot theme labs guides
#' scale_x_discrete guide_axis element_text guide_legend scale_fill_manual
#' geom_jitter position_jitter
#'
#' @export
#'
#' @examples
#' data(RawCounts_Antoszewski2022_MOUSEsub500)
#' ##------------------------------------------------------------------------##
#' resDATAprepSE <- DATAprepSE(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=NULL,
#'                             Individual.position=2)
#' ##------------------------------------------------------------------------##
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
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check inputs
    resErr <- ErrDATAplotBoxplotSamples(SEres=SEres,
                                        Log2.transformation=Log2.transformation,
                                        Colored.By.Factors=Colored.By.Factors,
                                        Color.Group=Color.Group,
                                        Plot.genes=Plot.genes,
                                        y.label=y.label)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## To avoid "no visible binding for global variable" with devtools::check()
    BioCond <- Time <- NULL

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Preprocessing

    Ndata <- length(SummarizedExperiment::assays(SEres))
    ExprData <- SummarizedExperiment::assays(SEres)[[Ndata]]
    cDat <- data.frame(SummarizedExperiment::colData(SEres))

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

    SEfinalname <- as.character(cDat$Final.Name)
    N.spl <- length(SEfinalname)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Data reshaped
    ExprData.f <- data.frame(Gene=as.character(seq_len(nrow(ExprData))),
                             ExprData)
    colnames(ExprData.f) <- c("Gene", SEfinalname)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Data used for boxplot
    Norm.dat.melt <- reshape2::melt(ExprData.f,
                                    id=c("Gene"),
                                    value.name="Expr",
                                    variable.name="Samples")
    colnames(Norm.dat.melt) <- c("Gene", "Samples", "Expression")
    Norm.dat.melt$Samples <- as.factor(Norm.dat.melt$Samples)

    if (isTRUE(Log2.transformation)) {
        Norm.dat.melt$Expression <- log2(Norm.dat.melt$Expression + 1)
        ylab.epr <- "log2 (Gene expression +1)"
    } else {
        ylab.epr <- "Gene expression"
    }## if(isTRUE(Log2.transformation))

    if (!is.null(y.label)) {
        ylab.epr <- y.label
    }## if(!is.null(y.label))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Graph
    Samples <- Expression <- NULL

    if (isFALSE(Colored.By.Factors)) {
        res.bxplt <- ggplot2::ggplot(Norm.dat.melt,
                                     ggplot2::aes(x=Samples, y=Expression))

        if (isTRUE(Plot.genes)) {
            Width.plt <- min(10/N.spl, 0.3)
            gg_jitter <- ggplot2::position_jitter(width=Width.plt, height=0.001)

            res.bxplt <- res.bxplt +
                ggplot2::geom_jitter(position=gg_jitter,
                                     alpha=0.7, fill="#56B4E9", color="#56B4E9")
        }## if (isTRUE(Plot.genes))

        res.bxplt <- res.bxplt +
            ggplot2::geom_boxplot(alpha=0.8, color="black", fill="#E69F00",
                                  show.legend=FALSE, outlier.alpha=0.2)

    } else {
        if (!is.null(SEinfoGroup)) {
            ##---------------------------------------------------------------##
            ## Setting
            FactorBoxplG <- as.character(SEinfoGroup)
            SEQlenFCT <- seq_len(length(FactorBoxplG))
            Num.FactorBC <- CharacterNumbers(Vect.number=SEQlenFCT)

            Expr.colG <- ExprData.f
            colnames(Expr.colG)<-c("Gene", paste0(Num.FactorBC, FactorBoxplG))

            Max.digit.G <- floor(log10(abs(max(SEQlenFCT))))
            Max.digit.G <- Max.digit.G + 1

            Melt1F.G <- reshape2::melt(Expr.colG, id=c("Gene"),
                                       value.name="Expr",
                                       variable.name="BioCond")
            substrg.BC <- substring(Melt1F.G$BioCond,
                                    first=Max.digit.G + 1,
                                    last=1000000L)
            dataBoxplot <- data.frame(Norm.dat.melt[, -1],
                                      substrg.BC)
            colnames(dataBoxplot) <- c("Samples", "Expression", "BioCond")

            ##---------------------------------------------------------------##
            Glevels <- levels(factor(SEinfoGroup))

            if (is.null(Color.Group)) {
                MypaletteG <- myPaletteBC(Nbc=length(Glevels))
                Color.Group <- data.frame(Name=Glevels, Col=MypaletteG)

            } else {
                Id.LevelCol.G <- order(Color.Group[, 1])
                Color.Group <- data.frame(Name=Glevels,
                                          Col=Color.Group[Id.LevelCol.G, 2])
            }## if(is.null(Color.Group))

            VARcolor <- factor(SEinfoGroup)
            levels(VARcolor) <- Color.Group$Col

            ##---------------------------------------------------------------##
            LegendTitle <- "Biological \nconditions"

            res.bxplt <- ggplot2::ggplot(dataBoxplot,
                                         ggplot2::aes(x=factor(Samples),
                                                      y=Expression,
                                                      fill=BioCond))
        }# if (!is.null(SEinfoGroup))

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        if (is.null(SEinfoGroup) & !is.null(SEinfoTime)) {
            ##---------------------------------------------------------------##
            Expr.colT <- ExprData.f
            FactorBoxplT <- as.character(SEinfoTime)
            FactorBoxplT_seq <- seq_len(length(FactorBoxplT))

            Num.FactorT <- CharacterNumbers(Vect.number=FactorBoxplT_seq)
            colnames(Expr.colT) <- c("Gene", paste0(Num.FactorT, FactorBoxplT))

            Max.digit.G <- floor(log10(abs(max(FactorBoxplT_seq))))
            Max.digit.G <- Max.digit.G + 1

            ##---------------------------------------------------------------##
            Melt1F.T <- reshape2::melt(Expr.colT, id=c("Gene"),
                                       value.name="Expr", variable.name="Time")
            substrg.T <- substring(Melt1F.T$Time, first=Max.digit.G + 1,
                                   last=1000000L)
            substrg.T2 <- gsub("T", "", gsub("t", "", as.character(substrg.T)))
            substrg.Tf <- paste0("t", substrg.T2)
            dataBoxplot <- data.frame(Norm.dat.melt[, -1], substrg.Tf)
            colnames(dataBoxplot) <- c("Samples", "Expression", "Time")

            ##---------------------------------------------------------------##
            Tlevels <- levels(factor(substrg.Tf))
            NbTime <- length(Tlevels)

            MypaletteT <- myPaletteT(NbTime)
            Color.Time <- data.frame(Name=Tlevels, Col=MypaletteT)

            VARcolor <- gsub("T", "", gsub("t", "", as.character(SEinfoTime)))
            VARcolor <- paste0("t", VARcolor)
            VARcolor <- factor(VARcolor)
            levels(VARcolor) <- Color.Time$Col

            ##---------------------------------------------------------------##
            LegendTitle <- "Time"

            res.bxplt <- ggplot2::ggplot(dataBoxplot,
                                         ggplot2::aes(x=factor(Samples),
                                                      y=Expression,
                                                      fill=Time))
        }##if(is.null(SEinfoGroup)&!is.null(SEinfoTime))

        ##-------------------------------------------------------------------##
        gg_strip <- ggplot2::element_text(size=9, color="black", face="bold")
        gg_jitter <- ggplot2::position_jitter(width=0.3, height=0.2)

        res.bxplt <- res.bxplt +
            ggplot2::theme(strip.text.x=gg_strip) +
            ggplot2::scale_fill_manual(values=levels(VARcolor))

        if (isTRUE(Plot.genes)) {
            JitterCol <- factor(dataBoxplot[,3])
            levels(JitterCol) <- levels(VARcolor)

            res.bxplt <- res.bxplt +
                ggplot2::geom_jitter(position=gg_jitter, colour=JitterCol,
                                     alpha=0.9)
        }## if (isTRUE(Plot.genes))

        res.bxplt <- res.bxplt +
            ggplot2::geom_boxplot(alpha=1, outlier.alpha=0.2) +
            ggplot2::guides(color=ggplot2::guide_legend(LegendTitle),
                            fill=ggplot2::guide_legend(LegendTitle))

    }## if(isFALSE(Colored.By.Factors))

    ##-----------------------------------------------------------------------##
    ## Final ggplot2
    if (isFALSE(Plot.genes)) {
        minGENEplot <- min(Norm.dat.melt$Expression)
    } else {
        minGENEplot <- min(Norm.dat.melt$Expression) - 0.1
    }## if (isFALSE(Plot.genes))

    res.bxplt <- res.bxplt +
        ggplot2::scale_x_discrete(guide=ggplot2::guide_axis(angle=90)) +
        ggplot2::labs(x="Samples", y=ylab.epr) +
        ggplot2::ylim(min=minGENEplot, max=max(Norm.dat.melt$Expression) + 0.1)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    return(res.bxplt=res.bxplt)
}## DATAplotBoxplotSamples()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

myPaletteBC <- function(Nbc=4) {
    ## paletteBC <- c(RColorBrewer::brewer.pal(8, "Dark2"),
    ##                RColorBrewer::brewer.pal(8, "Set2"))

    paletteBC <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                   "#66A61E", "#E6AB02", "#A6761D", "#666666",
                   "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                   "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")

    if (Nbc > 16) {
        ## pal2 <- scales::hue_pal(l=90)(Nbc - 1)
        pal2 <- gg_color_hue(n=Nbc-1, l=90)
        paletteBC <- c(paletteBC, pal2)
    }## if (Nbc > 16)

    paletteBC <- paletteBC[seq_len(Nbc)]

    return(paletteBC)
}## myPaletteBC()

myPaletteT <- function(Nt=4) {
    ## paletteT <- c("#737373", scales::hue_pal()(Nt - 1)) ## "#252525"
    paletteT <- c("#737373", gg_color_hue(n=Nt-1, l=65)) ## "#252525"
    return(paletteT)
}## myPaletteT()

gg_color_hue <- function(n=10, l=65) {
    ## gg_color_hue(20) == scales::hue_pal()(20)
    hues <- seq(15, 375, length= n + 1)
    resHUES <- grDevices::hcl(h=hues, l=l, c=100)[seq_len(n)]
    ## grDevices::hcl(h=hues, l=65, c=100)[seq_len(n)] == scales::hue_pal()(n)
    ## grDevices::hcl(h=hues, l=90, c=100)[1:n] == scales::hue_pal(l=90)(n)
    return(resHUES)
}## gg_color_hue

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrDATAplotBoxplotSamples <- function(SEres,
                                      Log2.transformation=TRUE,
                                      Colored.By.Factors=FALSE,
                                      Color.Group=NULL,
                                      Plot.genes=FALSE,
                                      y.label=NULL) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check DATAprepSE
    Err_SE <- paste0("'SEres' mut be the results of either the function ",
                     "'DATAprepSE()' or 'DATAnormalization()'.")

    if (!is(SEres, "SummarizedExperiment")) {
        stop(Err_SE)
    } else {
        codeDEres <- S4Vectors::metadata(SEres)$SEidentification
        if (!codeDEres%in%c("SEstep", "SEresNormalization")) {
            stop(Err_SE)
        }## if (!codeDEres%in%c("SEstep", "SEresNormalization"))
    }## if (!is(SEres, "SummarizedExperiment"))

    if (!isTRUE(Log2.transformation) & !isFALSE(Log2.transformation)) {
        stop("'Log2.transformation' must be TRUE or FALSE.")
    }## if (!isTRUE(Log2.transformation) & !isFALSE(Log2.transformation))

    if (!isTRUE(Colored.By.Factors) & !isFALSE(Colored.By.Factors)) {
        stop("'Colored.By.Factors' must be TRUE or FALSE.")
    }## if (!isTRUE(Colored.By.Factors) & !isFALSE(Colored.By.Factors))

    if (!isTRUE(Plot.genes) & !isFALSE(Plot.genes)) {
        stop("'Plot.genes' must be TRUE or FALSE.")
    }## if (!isTRUE(Plot.genes) & !isFALSE(Plot.genes))

    if (!is.null(Color.Group)) {
        if (!is.data.frame(Color.Group)) {
            stop("'Color.Group' must be NULL or a data.frame.")
        }## if (!is.data.frame(Color.Group))
    }## if (!is.null(Color.Group))

    if (!is.null(y.label)) {
        if (!is.character(y.label)) {
            stop("'y.label' must be NULL or a character.")
        }## if (!is.character(y.label))
    }## if (!is.null(y.label))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrDATAnormalization()
