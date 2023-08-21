#' @title Plot expression of one gene.
#'
#' @description The function allows to plot the gene expression profile of
#' one gene only according to time and/or biological conditions.
#'
#' @details All results are built from either the results of our R function
#' [DATAprepSE()]
#' or the results of our R function
#' [DATAnormalization()].
#'
#' @param SEres Results of either our R function
#' [DATAprepSE()],
#' or our R function
#' [DATAnormalization()].
#' @param row.gene Non negative integer indicating the row of the gene
#' to be plotted.
#' @param Color.Group NULL or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute
#' a color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#'
#' @importFrom SummarizedExperiment assays colData rownames
#' @importFrom S4Vectors metadata
#' @importFrom reshape2 melt
#' @importFrom stats var sd
#' @importFrom scales hue_pal
#' @importFrom ggplot2 ggplot xlab ylab aes geom_errorbar position_dodge
#' geom_line geom_point position_nudge geom_jitter position_jitter
#' geom_violin geom_boxplot geom_dotplot
#'
#' @return The function plots for the gene selected with
#' the input \code{row.gene}
#' * In the case where samples belong to different time points only :
#' the evolution of the expression of each replicate across time and
#' the evolution of the mean and the standard deviation of
#' the expression across time.
#' * In the case where samples belong to different biological conditions only:
#' a violin plot
#' (see [ggplot2::geom_violin()]),
#' and error bars (standard deviation)
#' (see [ggplot2::geom_errorbar()])
#' for each biological condition.
#' * In the case where samples belong to different time points and
#' different biological conditions : the evolution of the expression of
#' each replicate across time and the evolution of the mean and
#' the standard deviation of the expression across time
#' for each biological condition.
#'
#' @seealso The [DATAplotExpression1Gene()]
#' function is used by the following function of our package:
#' [DATAplotExpressionGenes()].
#'
#' @export
#'
#' @examples
#' ## Simulation raw counts
#' resSIMcount <- RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
#'                                    Nb.Gene=10)
#' ## Preprocessing step
#' resDATAprepSE <- DATAprepSE(RawCounts=resSIMcount$Sim.dat,
#'                             Column.gene=1,
#'                             Group.position=1,
#'                             Time.position=2,
#'                             Individual.position=3)
#' ##-------------------------------------------------------------------------#
#' resEVO1gene <- DATAplotExpression1Gene(SEres=resDATAprepSE,
#'                                        row.gene=1,
#'                                        Color.Group=NULL)
#' print(resEVO1gene)

DATAplotExpression1Gene <- function(SEres,
                                    row.gene=1,
                                    Color.Group=NULL) {
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check 1
    ## DATAprepSE
    Err_SE <- paste0("'SEres' mut be the results of either the function ",
                     "'DATAprepSE()' or 'DATAnormalization()'.")

    ## DATAprepSE
    if (!is(SEres, "SummarizedExperiment")) {
        stop(Err_SE)
    } else {
        codeDEres <- S4Vectors::metadata(SEres)$SEidentification

        if (is.null(codeDEres)) {
            stop(Err_SE)
        }## if (is.null(codeDEres))

        if (!codeDEres%in%c("SEstep", "SEresNormalization")) {
            stop(Err_SE)
        }## if (!codeDEres%in%c("SEstep", "SEresNormalization"))
    }## if (!is(SEres, "SummarizedExperiment"))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check 2
    if (!is.numeric(row.gene) & !is.integer(row.gene)) {
        stop("'row.gene' must be a non negative integer.")
    } else {
        if(floor(row.gene) != row.gene){
            stop("'row.gene' must be a non negative integer.")
        }## if(floor(Individual.position) != Individual.position)

        if (row.gene <= 0) {
            stop("'row.gene' must be a non negative integer.")
        }## if (row.gene <= 0)
    }## if(is.null(Individual.position))

    if (!is.null(Color.Group)) {
        if (!is.data.frame(Color.Group)) {
            stop("'Color.Group' must be NULL or a data.frame.")
        }## if (!is.data.frame(Color.Group))
    }## if (!is.null(Color.Group))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## preprocessing
    Ndata <- length(SummarizedExperiment::assays(SEres))
    ExprData <- SummarizedExperiment::assays(SEres)[[Ndata]]
    cDat <- data.frame(SummarizedExperiment::colData(SEres))
    NameG <- as.character(SummarizedExperiment::rownames(SEres))

    if (c("Group")%in%colnames(cDat)) {
        Vector.group <- as.character(cDat$Group)
    } else {
        Vector.group <- NULL
    }## if (c("Group")%in%colnames(cDat))

    if (c("Time")%in%colnames(cDat)) {
        Vector.time <- as.character(cDat$Time)
    } else {
        Vector.time <- NULL
    }## if (c("Time")%in%colnames(cDat))

    Vector.patient <- cDat$ID

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    Nb.per.cond <- as.numeric(table(Vector.patient))
    Var.per.cond <- stats::var(as.numeric(table(Vector.patient)))

    if (Var.per.cond == 0 & Nb.per.cond[1] > 1 & !is.null(Vector.time)) {
        Lign.draw <- TRUE
    }else{
        Lign.draw <- FALSE
    }## if(Var.per.cond == 0 & Nb.per.cond[1] > 1 & !is.null(Vector.time))

    GENEtitle <- paste0("Expression of gene ", NameG[row.gene])

    Expr.1G <- as.matrix(ExprData[row.gene,])

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## When samples belong to several biological conditions and time points
    if (!is.null(Vector.time) & !is.null(Vector.group)) {
        ##--------------------------------------------------------------------#
        Vector.groupF <- as.factor(Vector.group)
        Vector.timeF <- as.factor(Vector.time)
        ## To avoid "no visible binding for global variable" with
        ## devtools::check()
        Time <- Group <- Mean <- Sd <- value <- Patient <- NULL

        ##--------------------------------------------------------------------#
        Dat.1G.to.melt <- data.frame(GeneExpr=as.numeric(Expr.1G),
                                     Time=Vector.timeF,
                                     Group=Vector.groupF,
                                     Patient=Vector.patient)
        ##
        Dat.1G.melted <- reshape2::melt(Dat.1G.to.melt,
                                        id=c("Time", "Group", "Patient"))
        ##
        by.Dat.1G <- data.frame(Time=rep(levels(Vector.timeF),
                                         times=length(levels(Vector.groupF))),
                                Group=rep(levels(Vector.groupF),
                                          each=length(levels(Vector.timeF))),
                                Mean=as.numeric(by(Dat.1G.melted[, 5],
                                                   Dat.1G.melted[, c(1, 2)],
                                                   mean)),
                                Sd=as.numeric(by(Dat.1G.melted[, 5],
                                                 Dat.1G.melted[, c(1, 2)],
                                                 stats::sd)))

        ##--------------------------------------------------------------------#
        Glevels <- levels(factor(Vector.groupF))
        NbGroup <- length(Glevels)

        if (is.null(Color.Group)) {
            MypaletteG <- c(RColorBrewer::brewer.pal(8, "Dark2"),
                            RColorBrewer::brewer.pal(8, "Set2"))

            if (length(Glevels) > 16) {
                MypaletteG <- c(MypaletteG,
                                scales::hue_pal(l=90)(
                                    seq_len(length(Glevels)-1)))
            }## if(length(Glevels)>16)

            Color.Group <- data.frame(Name=Glevels,
                                      Col=MypaletteG[seq_len(length(Glevels))])
        } else {
            Id.LevelCol.G <- order(Color.Group[, 1])
            Color.Group <- data.frame(Name=Glevels,
                                      Col=Color.Group[Id.LevelCol.G, 2])
        }## if(is.null(Color.Group))

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## Graph profile expression
        Expr.plot <- ggplot2::ggplot(data=by.Dat.1G, ymin=0,
                                     ggplot2::aes(x=factor(Time),
                                                  y=as.numeric(Mean),
                                                  group=Group,
                                                  color=Group,
                                                  linetype=Group)) +
            ggplot2::theme(legend.position="bottom") +
            ggplot2::xlab("Time") + ggplot2::ylab("Expression") +
            ggplot2::ggtitle(GENEtitle) +
            ggplot2::geom_errorbar(data=by.Dat.1G, width=0.2, linewidth=0.7,
                                   ggplot2::aes(x=factor(Time),
                                                ymin=Mean-Sd, ymax=Mean+Sd,
                                                group=Group, color=Group,
                                                linetype=Group),
                                   position=ggplot2::position_dodge(0.2)) +
            ggplot2::geom_line(data=by.Dat.1G, linewidth=1.1,
                               ggplot2::aes(x=factor(Time),
                                            y=as.numeric(Mean),
                                            group=Group,
                                            color=Group,linetype=Group),
                               position=ggplot2::position_dodge(0.2)) +
            ggplot2::geom_point(data=by.Dat.1G, size=2.2,
                                ggplot2::aes(x=factor(Time), y=Mean,
                                             group=Group,
                                             color=Group,
                                             shape=Group),
                                position=ggplot2::position_dodge(0.2)) +
            ggplot2::scale_color_manual(values=as.character(Color.Group$Col))

        if (isTRUE(Lign.draw)) {
            Expr.plot <- Expr.plot+
                ggplot2::geom_line(data=Dat.1G.melted, alpha=0.5,
                                   ggplot2::aes(x=factor(Time),
                                                y=as.numeric(value),
                                                group=Patient,
                                                color=Group,
                                                linetype=Group),
                                   position=ggplot2::position_dodge(0.01))+
                ggplot2::geom_point(data=Dat.1G.melted, size=1.1,
                                    ggplot2::aes(x=factor(Time), y=value,
                                                 group=Group, color=Group,
                                                 shape=Group),
                                    position=ggplot2::position_dodge(0.01))
        } else {
            Expr.plot <- Expr.plot+
                ggplot2::geom_point(data=Dat.1G.melted, size=1.1,
                                    ggplot2::aes(x=factor(Time),
                                                 y=value,
                                                 group=Group, color=Group,
                                                 shape=Group),
                                    position=ggplot2::position_dodge(0.01))
        }## if(isTRUE(Lign.draw))
    }## if(is.null(Vector.time)==FALSE & is.null(Vector.group)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## When samples belong to several time points but one biological condition
    if (!is.null(Vector.time) & is.null(Vector.group)) {
        ##--------------------------------------------------------------------#
        Vector.groupF <- as.factor(rep("G1",times=length(as.numeric(Expr.1G))))
        Vector.timeF <- as.factor(Vector.time)
        ## To avoid "no visible binding for global variable" with
        ## devtools::check()
        Time <- Mean <- Sd <- value <- Patient <- NULL

        ##--------------------------------------------------------------------#
        Dat.1G.to.melt <- data.frame(GeneExpr=as.numeric(Expr.1G),
                                     Time=Vector.timeF,
                                     Patient=Vector.patient)
        Dat.1G.melted <- reshape2::melt(Dat.1G.to.melt,
                                        id=c("Time","Patient"))
        by.Dat.1G <- data.frame(Time=levels(Vector.timeF),
                                Mean=as.numeric(by(Dat.1G.melted[,4],
                                                   Dat.1G.melted[,1],
                                                   mean)),
                                Sd=as.numeric(by(Dat.1G.melted[,4],
                                                 Dat.1G.melted[,1],
                                                 stats::sd)))

        Col.grougF <- Vector.groupF
        levels(Col.grougF) <- "#E76BF3"

        ##--------------------------------------------------------------------#
        Expr.plot <- ggplot2::ggplot(ymin=0) +
            ggplot2::theme(legend.position="bottom") +
            ggplot2::xlab("Time")+ ggplot2::ylab("Expression") +
            ggplot2::ggtitle(GENEtitle) +
            ggplot2::geom_errorbar(data=by.Dat.1G,
                                   linetype="dashed", width=.2, linewidth=0.7,
                                   ggplot2::aes(x=factor(Time),
                                                ymin=Mean-Sd,
                                                ymax=Mean+Sd),
                                   position=ggplot2::position_nudge(x=0.1,
                                                                    y=0)) +
            ggplot2::geom_line(data=by.Dat.1G, color="black", linewidth=1.1,
                               ggplot2::aes(x=factor(Time),
                                            y=as.numeric(Mean), group=1),
                               position=ggplot2::position_nudge(x=0.1, y=0)) +
            ggplot2::geom_point(data=by.Dat.1G, size=2.2, color="black",
                                ggplot2::aes(x=factor(Time), y=Mean),
                                position=ggplot2::position_nudge(x=0.1, y=0))

        if (isTRUE(Lign.draw)) {
            Expr.plot <- Expr.plot+
                ggplot2::geom_line(data=Dat.1G.melted,
                                   colour=Col.grougF, alpha=0.4,
                                   ggplot2::aes(x=factor(Time),
                                                y=value, group=Patient),
                                   position=ggplot2::position_nudge(x=0, y=0))+
                ggplot2::geom_point(data=Dat.1G.melted, colour=Col.grougF,
                                    size=1.1, alpha=0.4,
                                    ggplot2::aes(x=factor(Time),
                                                 y=value, group=Patient),
                                    position=ggplot2::position_nudge(x=0, y=0))
        } else {
            Expr.plot<-Expr.plot+
                ggplot2::geom_point(data=Dat.1G.melted, colour=Col.grougF,
                                    size=1.1, alpha=0.4,
                                    ggplot2::aes(x=factor(Time),
                                                 y=value, group=Patient),
                                    position=ggplot2::position_nudge(x=0, y=0))
        }## if(Lign.draw==TRUE)
    }## if(!is.null(Vector.time) & is.null(Vector.group))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## When samples belong to several biological conditions but one time point
    if (is.null(Vector.time) & !is.null(Vector.group)) {
        ##--------------------------------------------------------------------#
        Vector.groupF <- as.factor(Vector.group)
        Glevels <- levels(factor(Vector.groupF))
        NbGroup <- length(Glevels)

        ##--------------------------------------------------------------------#
        if (is.null(Color.Group)) {
            MypaletteG <- c(RColorBrewer::brewer.pal(8,"Dark2"),
                            RColorBrewer::brewer.pal(8,"Set2"))
            #
            if (length(Glevels) > 16) {
                MypaletteG <- c(MypaletteG,
                                scales::hue_pal(l=90)(
                                    seq_len(length(Glevels)-1)))
            }## if(length(Glevels)>16)

            Color.Group <- data.frame(Name=Glevels,
                                      Col=MypaletteG[seq_len(length(Glevels))])
        } else {
            Id.LevelCol.G <- order(Color.Group[, 1])
            Color.Group <- data.frame(Name=Glevels,
                                      Col=Color.Group[Id.LevelCol.G, 2])
        }## if(is.null(Color.Group))

        ##--------------------------------------------------------------------#
        ## To avoid "no visible binding for global variable" with
        ## devtools::check()
        Mean <- Sd <- value <- Patient <- NULL
        #
        Dat.1G.to.melt <- data.frame(GeneExpr=as.numeric(Expr.1G),
                                     Group=Vector.groupF,
                                     Patient=Vector.patient)
        Dat.1G.melted <- reshape2::melt(Dat.1G.to.melt,
                                        id=c("Group", "Patient"))
        by.Dat.1G <- data.frame(Group=levels(Vector.groupF),
                                Mean=as.numeric(by(Dat.1G.melted[, 4],
                                                   Dat.1G.melted[, 1],
                                                   mean)),
                                Sd=as.numeric(by(Dat.1G.melted[, 4],
                                                 Dat.1G.melted[, 1],
                                                 stats::sd)))
        #
        size.pt <- (max(Dat.1G.melted$value) - min(Dat.1G.melted$value))/50

        ##--------------------------------------------------------------------#
        Expr.plot <- ggplot2::ggplot(ymin=0) +
            ggplot2::xlab("Biological condition") +
            ggplot2::ylab("Expression") +
            ggplot2::ggtitle(GENEtitle) +
            ggplot2::geom_violin(data=Dat.1G.melted, trim=TRUE,
                                 ggplot2::aes(x=Group, y=value, fill=Group))+
            ggplot2::geom_boxplot(data=Dat.1G.melted,
                                  ggplot2::aes(x=Group, y=value),width=0.1)+
            ggplot2::geom_dotplot(data=Dat.1G.melted, colour="black",
                                  ggplot2::aes(x=factor(Group),
                                               y=value,
                                               fill=Group),
                                  binaxis='y', stackdir='center',
                                  binwidth=size.pt)+
            ggplot2::geom_errorbar(data=by.Dat.1G, width=0.2, linewidth=0.9,
                                   ggplot2::aes(x=factor(Group),
                                                ymin=Mean-Sd, ymax=Mean+Sd),
                                   position = ggplot2::position_nudge(x=0,
                                                                      y=0)) +
            ggplot2::geom_point(data=by.Dat.1G,
                                size=2.2, shape=15, colour="blue", fill="blue",
                                ggplot2::aes(x=factor(Group), y=Mean),
                                position = ggplot2::position_nudge(x=0, y=0))+
            ggplot2::scale_fill_manual(values=as.character(Color.Group$Col))
    }## if(is.null(Vector.time) & !is.null(Vector.group))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(Expr.plot=Expr.plot)
}## DATAplotExpression1Gene()
