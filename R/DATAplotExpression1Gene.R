#' @title Plot expression of one gene.
#'
#' @description The function allows to plot the gene expression profile of
#' one gene only according to time and/or biological conditions.
#'
#' @details The column names of \code{ExprData} must be a vector of strings
#' of characters containing
#' * a string of characters (if \eqn{k=1}) which is the label of the column
#' containing gene names.
#' * \eqn{N_s} sample names which must be strings of characters containing
#' at least : the name of the individual (e.g patient, mouse, yeasts culture),
#' its biological condition (if there is at least two) and
#' the time where data have been collected if there is at least two;
#' (must be either 't0', 'T0' or '0' for time 0,
#' 't1', 'T1' or '1' for time 1, ...).
#'
#' All these sample information must be separated by underscores in
#' the sample name. For instance 'CLL_P_t0_r1',
#' corresponds to the patient 'r1' belonging to the biological condition 'P'
#' and where data were collected at time 't0'. I this example, 'CLL' describe
#' the type of cells (here chronic lymphocytic leukemia)
#' and is not used in our analysis.
#'
#' In the string of characters 'CLL_P_t0_r1',
#' 'r1' is localized after the third underscore,
#' so \code{Individual.position=4},
#' 'P' is localized after the first underscore, so \code{Group.position=2} and
#' 't0' is localized after the second underscore, so \code{Time.position=3}.
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
#' @param row.gene Non negative integer indicating the row of the gene
#' to be plotted.
#' @param Column.gene Integer indicating the column where gene names are given.
#' Set \code{Column.gene=NULL} if there is no such column.
#' @param Group.position Integer indicating the position of group information
#' in the string of characters in each sample names (see \code{Details}).
#' Set \code{Group.position=NULL} if there is only one or no biological
#' information in the string of character in each sample name.
#' @param Time.position Integer indicating the position of time measurement
#' information in the string of characters in each sample names
#' (see \code{Details}).
#' Set \code{Time.position=NULL} if there is only one or no time measurement
#' information in the string of character in each sample name.
#' @param Individual.position Integer indicating the position of the name of
#' the individual (e.g patient, replicate, mouse, yeasts culture ...)
#' in the string of characters in each sample names (see \code{Details}).
#' The names of different individuals must be all different.
#' Furthermore, if individual names are just numbers, they will be transform
#' in a vector of class "character" by [CharacterNumbers()] and
#' a "r" will be added to each individual name ("r" for replicate).
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
#' @importFrom reshape2 melt
#' @importFrom stats var sd
#' @importFrom scales hue_pal
#' @importFrom ggplot2 ggplot xlab ylab aes geom_errorbar position_dodge
#' geom_line geom_point position_nudge geom_jitter position_jitter
#' geom_violin geom_boxplot geom_dotplot
#'
#' @return The function plots for each gene selected with
#' the input \code{Vector.row.gene}
#' * In the case where samples belong to different time points only :
#' the evolution of the expression of each replicate across time and
#' the evolution of the mean and the standard deviation of
#' the expression across time.
#' * In the case where samples belong to different biological conditions only :
#' a violin plot (see [ggplot2::geom_violin()]),
#' and error bars (standard deviation) (see [ggplot2::geom_errorbar()])
#' for each biological condition.
#' * In the case where samples belong to different time points and
#' different biological conditions : the evolution of the expression of
#' each replicate across time and the evolution of the mean and
#' the standard deviation of the expression across time
#' for each biological condition.
#'
#' @export
#'
#' @examples
#' ## Simulation raw counts
#' res.sim.count<-RawCountsSimulation(Nb.Group=2,Nb.Time=3,Nb.per.GT=4,
#'                                    Nb.Gene=10)
#' ##
#' Res.evo<-DATAplotExpression1Gene(ExprData=res.sim.count$Sim.dat,
#'                                  row.gene=1,
#'                                  Column.gene=1,
#'                                  Group.position=1,
#'                                  Time.position=2,
#'                                  Individual.position=3,
#'                                  Color.Group=NULL)
#' print(Res.evo)

DATAplotExpression1Gene<-function(ExprData,
                                  row.gene,
                                  Column.gene,
                                  Group.position,
                                  Time.position,
                                  Individual.position,
                                  Color.Group=NULL){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    res.colname<-ColnamesToFactors(ExprData=ExprData,
                                   Column.gene=Column.gene,
                                   Group.position=Group.position,
                                   Time.position=Time.position,
                                   Individual.position=Individual.position)
    Vector.group<-res.colname$Group.Info
    Vector.time<-res.colname$Time.Info
    Vector.patient<-res.colname$Individual.info

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    Nb.per.cond<-as.numeric(table(Vector.patient))
    Var.per.cond<-stats::var(as.numeric(table(Vector.patient)))

    if(Var.per.cond == 0 & Nb.per.cond[1] > 1 & !is.null(Time.position)){
        Lign.draw<-TRUE
    }else{
        Lign.draw<-FALSE
    }## if(Var.per.cond == 0 & Nb.per.cond[1] > 1 & !is.null(Time.position))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check
    if(is.null(Vector.time) & is.null(Vector.group)){
        stop("Time and/or biological condition must be selected")
    }## if(is.null(Vector.time) & is.null(Vector.group))

    ## Data with only expression
    if(is.null(Column.gene)){
        ind.col.expr<-seq_len(ncol(ExprData))
    }else{
        ind.col.expr<-seq_len(ncol(ExprData))[-Column.gene]
    }## if(is.null(Column.gene))
    Expr.1G<-as.matrix(ExprData[row.gene, ind.col.expr])

    # Gene names
    if(is.null(Column.gene)){
        if(is.null(row.names(ExprData))){
            Name.G<-paste0("Expression of gene ", row.gene)
        }else{
            Name.G<-paste0("Expression of gene ",
                           row.names(ExprData)[row.gene])
        }## if(is.null(row.names(ExprData)))
    }else{
        Name.G<-paste0("Expression of gene ", ExprData[row.gene, Column.gene])
    }## if(is.null(Column.gene))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## When samples belong to several biological conditions and time points
    if(!is.null(Vector.time) & !is.null(Vector.group)){
        ##--------------------------------------------------------------------#
        Vector.groupF<-as.factor(Vector.group)
        Vector.timeF<-as.factor(Vector.time)
        ## To avoid "no visible binding for global variable" with
        ## devtools::check()
        Time<-Group<-Mean<-Sd<-value<-Patient<-NULL

        ##--------------------------------------------------------------------#
        Dat.1G.to.melt<-data.frame(GeneExpr=as.numeric(Expr.1G),
                                   Time=Vector.timeF,
                                   Group=Vector.groupF,
                                   Patient=Vector.patient)
        #
        Dat.1G.melted<-reshape2::melt(Dat.1G.to.melt,
                                      id=c("Time", "Group", "Patient"))
        #
        by.Dat.1G<-data.frame(Time=rep(levels(Vector.timeF),
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
        Glevels<-levels(factor(Vector.groupF))
        NbGroup<-length(Glevels)

        if(is.null(Color.Group)){
            MypaletteG<-c(RColorBrewer::brewer.pal(8, "Dark2"),
                          RColorBrewer::brewer.pal(8, "Set2"))

            if(length(Glevels)>16){
                MypaletteG<-c(MypaletteG,
                              scales::hue_pal(l=90)(seq_len(length(Glevels)-1)))
            }## if(length(Glevels)>16)

            Color.Group<-data.frame(Name=Glevels,
                                    Col=MypaletteG[seq_len(length(Glevels))])
        }else{
            Id.LevelCol.G<-order(Color.Group[, 1])
            Color.Group<-data.frame(Name=Glevels,
                                    Col=Color.Group[Id.LevelCol.G, 2])
        }## if(is.null(Color.Group))

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## Graph profile expression
        Expr.plot<-ggplot2::ggplot(data=by.Dat.1G, ymin=0,
                                   ggplot2::aes(x=factor(Time),
                                                y=as.numeric(Mean),
                                                group=Group,
                                                color=Group,
                                                linetype=Group)) +
            ggplot2::theme(legend.position="bottom") +
            ggplot2::xlab("Time") + ggplot2::ylab("Expression") +
            ggplot2::ggtitle(Name.G) +
            ggplot2::geom_errorbar(data=by.Dat.1G, width=.2, size=0.7,
                                   ggplot2::aes(x=factor(Time),
                                                ymin=Mean-Sd, ymax=Mean+Sd,
                                                group=Group, color=Group,
                                                linetype=Group),
                                   position=ggplot2::position_dodge(0.2)) +
            ggplot2::geom_line(data=by.Dat.1G, size=1.1,
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

        if(isTRUE(Lign.draw)){
            Expr.plot<-Expr.plot+
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
        }else{
            Expr.plot<-Expr.plot+
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
    if(!is.null(Vector.time) & is.null(Vector.group)){
        ##--------------------------------------------------------------------#
        Vector.groupF<-as.factor(rep("G1",times=length(as.numeric(Expr.1G))))
        Vector.timeF<-as.factor(Vector.time)
        ## To avoid "no visible binding for global variable" with
        ## devtools::check()
        Time<-Mean<-Sd<-value<-Patient<-NULL

        ##--------------------------------------------------------------------#
        Dat.1G.to.melt<-data.frame(GeneExpr=as.numeric(Expr.1G),
                                   Time=Vector.timeF,
                                   Patient=Vector.patient)
        Dat.1G.melted <- reshape2::melt(Dat.1G.to.melt, id=c("Time","Patient"))
        by.Dat.1G<-data.frame(Time=levels(Vector.timeF),
                              Mean=as.numeric(by(Dat.1G.melted[,4],
                                                 Dat.1G.melted[,1],
                                                 mean)),
                              Sd=as.numeric(by(Dat.1G.melted[,4],
                                               Dat.1G.melted[,1],
                                               stats::sd)))

        Col.grougF<-Vector.groupF
        levels(Col.grougF)<-"#E76BF3"

        ##--------------------------------------------------------------------#
        Expr.plot<-ggplot2::ggplot(ymin=0) +
            ggplot2::theme(legend.position="bottom") +
            ggplot2::xlab("Time")+ ggplot2::ylab("Expression") +
            ggplot2::ggtitle(Name.G) +
            ggplot2::geom_errorbar(data=by.Dat.1G,
                                   linetype="dashed", width=.2, size=0.7,
                                   ggplot2::aes(x=factor(Time),
                                                ymin=Mean-Sd,
                                                ymax=Mean+Sd),
                                   position=ggplot2::position_nudge(x=0.1,
                                                                    y=0)) +
            ggplot2::geom_line(data=by.Dat.1G, color="black", size=1.1,
                               ggplot2::aes(x=factor(Time),
                                            y=as.numeric(Mean), group=1),
                               position=ggplot2::position_nudge(x=0.1, y=0)) +
            ggplot2::geom_point(data=by.Dat.1G, size=2.2, color="black",
                                ggplot2::aes(x=factor(Time), y=Mean),
                                position=ggplot2::position_nudge(x=0.1, y=0))

        if(isTRUE(Lign.draw)){
            Expr.plot<-Expr.plot+
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
        }else{
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
    if(is.null(Vector.time) & !is.null(Vector.group)){
        ##--------------------------------------------------------------------#
        Vector.groupF<-as.factor(Vector.group)
        Glevels<-levels(factor(Vector.groupF))
        NbGroup<-length(Glevels)

        ##--------------------------------------------------------------------#
        if(is.null(Color.Group)){
            MypaletteG<-c(RColorBrewer::brewer.pal(8,"Dark2"),
                          RColorBrewer::brewer.pal(8,"Set2"))
            #
            if(length(Glevels)>16){
                MypaletteG<-c(MypaletteG,
                              scales::hue_pal(l=90)(seq_len(length(Glevels)-1)))
            }## if(length(Glevels)>16)

            Color.Group<-data.frame(Name=Glevels,
                                    Col=MypaletteG[seq_len(length(Glevels))])
        }else{
            Id.LevelCol.G<-order(Color.Group[, 1])
            Color.Group<-data.frame(Name=Glevels,
                                    Col=Color.Group[Id.LevelCol.G, 2])
        }## if(is.null(Color.Group))

        ##--------------------------------------------------------------------#
        ## To avoid "no visible binding for global variable" with
        ## devtools::check()
        Mean<-Sd<-value<-Patient<-NULL
        #
        Dat.1G.to.melt<-data.frame(GeneExpr=as.numeric(Expr.1G),
                                   Group=Vector.groupF,
                                   Patient=Vector.patient)
        Dat.1G.melted<-reshape2::melt(Dat.1G.to.melt,
                                      id=c("Group", "Patient"))
        by.Dat.1G<-data.frame(Group=levels(Vector.groupF),
                              Mean=as.numeric(by(Dat.1G.melted[, 4],
                                                 Dat.1G.melted[, 1],
                                                 mean)),
                              Sd=as.numeric(by(Dat.1G.melted[, 4],
                                               Dat.1G.melted[, 1],
                                               stats::sd)))
        #
        size.pt<-(max(Dat.1G.melted$value)-min(Dat.1G.melted$value))/50

        ##--------------------------------------------------------------------#
        Expr.plot<-ggplot2::ggplot(ymin=0) +
            ggplot2::xlab("Biological condition") +
            ggplot2::ylab("Expression") +
            ggplot2::ggtitle(Name.G) +
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
            ggplot2::geom_errorbar(data=by.Dat.1G, width=.2, size=0.9,
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
