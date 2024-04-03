#' @title Alluvial graphs of differentially expressed (DE) genes
#'
#' @description
#' The function takes as input a binary table with \eqn{N_g} lines
#' corresponding to genes and
#' * if \code{Temporal.Group=TRUE} : \eqn{T-1} columns corresponding to times
#' (with \eqn{T} the number of time points).
#' A '1' in the n-th row and t-th column means that the n-th gene is
#' differentially expressed (DE) at time t, compared with
#' the reference time t0.
#' * if \code{Temporal.Group=FALSE} :
#' \eqn{G} columns corresponding to the number of group.
#' A '1' in the \eqn{n}-th row and \eqn{g}-th column means
#' that the n-th gene is
#'   * DE at least one time ti, compared with the reference time t0,
#'   for the group \eqn{g}.
#'   * specific at least one time ti, compared with the reference time t0,
#'   for the group \eqn{g} (see [DEanalysisTimeAndGroup()]
#'   for the notion 'specific').
#'   * a signature gene at least one time ti, compared with the reference time
#'   t0, for the group \eqn{g} (see [DEanalysisTimeAndGroup()]
#'   for the notion 'signature').
#'
#' The function plots
#' * if \code{Temporal.Group=TRUE}, two graphs: an alluvial graph and
#' a plot showing the time evolution of the number of DE genes within
#' each temporal group. By temporal group, we mean the sets of genes which
#' are first DE at the same time.
#' * if \code{Temporal.Group=FALSE} : an alluvial graph.
#'
#' @param table.DE.time Binary matrix (table filled with 0 and 1) with
#' \eqn{N_g} rows and \eqn{T-1} columns with \eqn{N_g} the number of genes and
#' \eqn{T-1} the number of time points.
#' @param Temporal.Group \code{TRUE} or \code{FALSE},
#' \code{FALSE} as default (see \code{Description}).
#' @param title.alluvial String of characters or \code{NULL},
#' \code{NULL} as default. The input \code{title.allluvial} corresponds
#' to the title of the alluvial graph.
#' If \code{title} is a string of characters,
#' \code{title} will be the title of the alluvial graph.
#' If \code{title=NULL}, the title of the alluvial graph will be
#' 'Alluvial graph'.
#' @param title.evolution String of characters or \code{NULL},
#' \code{NULL} as default. Only applied if \code{Temporal.Group=TRUE}.
#' The input \code{title.evolution} corresponds to the title of
#' the second graph (see \code{Description}).
#' If \code{title} is a string of characters, it will be to the title of
#' the second graph.
#' If \code{title=NULL}, the title of the second graph will be
#' 'Time evolution of the number of DE genes within each temporal group'.
#'
#' @details The names of the columns of the table will be the axis labels
#' in the plots.
#' If the table has no column names, the function will automatically create
#' column names (t1,t2,...).
#'
#' @return The function returns, as described in \code{description}
#' * if \code{Temporal.Group=TRUE}, two graphs: an alluvial graph and
#' a plot showing the time evolution of the number of DE genes within
#' each temporal group.
#' By temporal group, we mean the sets of genes which are first DE
#' at the same time.
#' * if \code{Temporal.Group=FALSE} : an alluvial graph.
#'
#' @seealso The [DEplotAlluvial()] function
#' * is used by the following functions of our package : [DEanalysisTime()]
#' and [DEanalysisTimeAndGroup()].
#' * calls the R package [ggplot2] in order to plot the two graphs.
#'
#' @importFrom ggalluvial stat_stratum stat_alluvium
#' @importFrom ggplot2 ggplot aes theme scale_x_continuous ggtitle labs guides
#' geom_line geom_point guide_legend geom_area scale_size_manual xlab ylab
#' scale_fill_manual
#' @importFrom stats aggregate
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#' set.seed(1994)
#'
#' NbTime.vst0 <- 4
#' BinTable <- matrix(sample(c(0,1),replace=TRUE,
#'                           size=NbTime.vst0*120,c(0.60,0.40)),
#'                    ncol=NbTime.vst0)
#' colnames(BinTable) <- paste0("t", 1:NbTime.vst0)
#'
#' ##------------------------------------------------------------------------##
#' res.alluvial <- DEplotAlluvial(table.DE.time=BinTable)
#' print(res.alluvial$g.alluvial)
#' print(res.alluvial$g.alluvial.freq)

DEplotAlluvial <- function(table.DE.time,
                           Temporal.Group=TRUE,
                           title.alluvial=NULL,
                           title.evolution=NULL){
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## To avoid "no visible binding for global variable" with
    ## devtools::check()
    variable <- Freq <- Gene <- value <- Category <- NULL
    Attribute <- Time_group <- NULL

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 1) Parameters : column names
    nTimes <- ncol(table.DE.time)

    if (is.null(colnames(table.DE.time))) {
        DATAcolnames <- paste0("t", seq_len(nTimes))
    } else {
        DATAcolnames <- colnames(table.DE.time)
    }## if(is.null(colnames(table.DE.time)))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 2) Data for graphics, common to Temporal.Group=TRUE or FALSE
    ## Only DE genes at least at one times versus t0 are kept
    DEgenes_number <- which(apply(table.DE.time, 1, function(x) sum(x)) != 0)
    nDEgenes <- length(DEgenes_number)

    if (nDEgenes == 1) {
        tableDEtime_DEonly <- t(as.matrix(table.DE.time[DEgenes_number,]))
    } else {
        tableDEtime_DEonly <- table.DE.time[DEgenes_number,]
    }## if (nDEgenes==1)

    colnames(tableDEtime_DEonly) <- DATAcolnames

    if (is.null(row.names(tableDEtime_DEonly))) {
        Gnames <- as.character(seq_len(nDEgenes))
    } else {
        Gnames <- row.names(tableDEtime_DEonly)
    }# if(is.null(row.names(tableDEtime_DEonly)))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 3) graphics parameters common to Temporal.Group=TRUE or FALSE: colors...

    colorguideList <- list(colour=c("black"), fill=c('black', 'grey'))
    colorLegend <- ggplot2::guide_legend(override.aes=colorguideList)

    alluvLabels <- DATAcolnames

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (isTRUE(Temporal.Group)) {
        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 2) Data for graphics, specific to Temporal.Group=TRUE

        ##-------------------------------------------------------------------##
        ## Position of DE genes in table.group.gene.DE.peak. Temporal group is
        ## calculated and corresponds to the first time a gene is DE.
        Time.first.peak <- apply(as.matrix(tableDEtime_DEonly),
                                 MARGIN=1,
                                 FUN=function(x) which(cumsum(x) == 1)[1])

        DE.genes.info <- data.frame(Index=DEgenes_number,
                                    Name=Gnames,
                                    T.group=as.numeric(Time.first.peak))

        ##-------------------------------------------------------------------##
        ## Data for alluvial graph per gene :
        ## Patterns (1 for DE vs t0), time group & Freq
        Freq1ini <- rep(1, times=nDEgenes)
        alluvData <- cbind(tableDEtime_DEonly,
                           data.frame(Time_group=DE.genes.info$T.group,
                                      Freq=Freq1ini))
        alluvData_agg <- alluvData[,seq_len(nTimes+1)]

        ## Freq=rep(1, nrow(alluvData)) = Freq1ini
        allu2DataTini <- stats::aggregate(x=list(Freq=Freq1ini),
                                          by=alluvData_agg,
                                          FUN=length)
        sortTgr <- sort(unique(allu2DataTini$Time_group))

        ##-------------------------------------------------------------------##
        allu2Data.t <- allu2DataTini

        if (length(unique(allu2DataTini$Time_group)) <= nTimes) {
            add.TG <- cbind(diag(rep(1,times=nTimes)),
                            seq_len(nTimes),
                            rep(0,times=nTimes))
            colnames(add.TG) <- colnames(allu2DataTini)
            allu2Data.t <- rbind(allu2Data.t, add.TG[-sortTgr,])
        }## if(length(unique(allu2DataTini$Time_group)) <= nTimes)

        allu2Data.t <- cbind(Gene=seq_len(nrow(allu2Data.t)), allu2Data.t)

        allu2Data.tf <- reshape2::melt(data=allu2Data.t,
                                       id.vars=c("Gene", "Time_group", "Freq"))

        allu2Data.tf$Time_group <- as.factor(allu2Data.tf$Time_group)
        allu2Data.tf$Gene <- as.factor(allu2Data.tf$Gene)
        allu2Data.tf$variable <- as.numeric(allu2Data.tf$variable)*50

        allu2Data.tf$value <- as.factor(allu2Data.tf$value)
        allu2valueKeep <- as.numeric(levels(allu2Data.tf$value)) + 1
        levels(allu2Data.tf$value) <- c("no", "yes")[allu2valueKeep]

        ##-------------------------------------------------------------------##
        ## Data for alluvial graph : number DE gene per temporal group
        ## at each time
        allu2Data_prep <- stats::aggregate(alluvData[,seq_len(nTimes)],
                                           by=list(Category=alluvData$Time_group),
                                           FUN=sum)

        if (length(unique(allu2Data.t$Time_group)) <= nTimes) {
            add.TG.2 <- cbind(seq_len(nTimes),
                              diag(rep(0,times=nTimes)))
            colnames(add.TG.2) <- colnames(allu2Data_prep)

            allu2Data_prep <- rbind(allu2Data_prep, add.TG.2[-sortTgr,])
        }## if (length(unique(allu2Data.t$Time_group)) <= nTimes)

        allu2Data <- reshape2::melt(data=allu2Data_prep, id.vars=c("Category"))
        Category_str <- as.character(allu2Data$Category)

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 3) graphics parameters: colors, labels ...

        allu2DataTini_apply <- allu2DataTini[,seq_len(nTimes)]
        Id.T.0DE <- which(apply(allu2DataTini_apply, 2, FUN=function(x) 0%in%x))
        Id.T.1DE <- which(apply(allu2DataTini_apply, 2, FUN=function(x) 1%in%x))
        colorbarFILLini <- rep(c('black', 'grey'), times=nTimes)

        if (length(Id.T.1DE) == 0 | length(Id.T.0DE) == 0) {
            if (length(Id.T.1DE) == 0) {
                colorbarFILLkeep <- 2*(as.numeric(Id.T.0DE) - 1) + 1
            } else {
                colorbarFILLkeep <- 2*as.numeric(Id.T.1DE)
            }## if (length(Id.T.1DE) == 0)
        } else {
            colorbarFILLkeep <- sort(c(2*(as.numeric(Id.T.0DE) - 1) + 1,
                                       2*as.numeric(Id.T.1DE)))
        }## if(length(Id.T.1DE)==0 | length(Id.T.0DE)==0)

        colorbarFILL <- colorbarFILLini[colorbarFILLkeep]

        qFreq_guides <- ggplot2::guide_legend(title="Temporal group",
                                              title.position="left")

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 3) Graphs
        ## Graph alluvial 1
        q.alluvial <- ggplot2::ggplot(allu2Data.tf,
                                      ggplot2::aes(x=variable, y=Freq,
                                                   fill=Time_group,
                                                   alluvium=Gene,
                                                   stratum=value)) +
            ggalluvial::stat_alluvium(geom="flow", lode.guidance="forward",
                                      width=5, reverse=FALSE) +
            ggalluvial::stat_stratum(width=5,
                                     mapping=ggplot2::aes(size=value),
                                     fill=colorbarFILL, reverse=FALSE,
                                     alpha=0.7) +
            ggplot2::labs(x="Time", y="Number of DE genes",
                          fill="Temporal group") +
            ggplot2::scale_x_continuous(breaks=seq_len(nTimes)*50,
                                        labels=alluvLabels) +
            ggplot2::scale_size_manual("DE tx versus t0", values=c(0.5 ,0.5),
                                       guide=colorLegend)

        if (is.null(title.alluvial)) {
            q.alluvial <- q.alluvial + ggplot2::ggtitle("Alluvial graph")
        } else {
            q.alluvial <- q.alluvial + ggplot2::ggtitle(title.alluvial)
        }## if(is.null(title.alluvial))

        ##-------------------------------------------------------------------##
        ## Graph alluvial 2
        qFreq_perTperBC <- ggplot2::ggplot(data=allu2Data,
                                           ggplot2::aes(x=variable,
                                                        y=value,
                                                        group=Category_str)) +
            ggplot2::guides(color=qFreq_guides, fill=qFreq_guides) +
            ggplot2::geom_line(ggplot2::aes(color=Category_str), linewidth=1) +
            ggplot2::geom_point(ggplot2::aes(color=Category_str), size=2) +
            ggplot2::geom_area(aes(fill=Category_str, group=Category_str),
                               alpha=0.2, position='identity') +
            ggplot2::theme(axis.text.x=ggplot2::element_text(angle=0,
                                                             hjust=0.5),
                           legend.position="bottom") +
            ggplot2::xlab("Time") +
            ggplot2::ylab("Number of DE genes")

        if (is.null(title.evolution)) {
            qFreq_perTperBC <- qFreq_perTperBC +
                ggplot2::ggtitle(paste0("Time evolution of the number of ",
                                        "DE genes within each temporal group."))
        } else {
            qFreq_perTperBC <- qFreq_perTperBC +
                ggplot2::ggtitle(title.evolution)
        }## if(is.null(title.evolution))

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 4) Output
        return(list(g.alluvial=q.alluvial,
                    g.alluvial.freq=qFreq_perTperBC))
    } else {

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 2) Data for graphs, specific to Temporal.Group=FALSE

        ##-------------------------------------------------------------------##
        ## Position of DE genes in table.group.gene.DE.peak
        ## Temporal group is calculated and corresponds to the first time
        ## a gene is DE.
        DE.nb.G <- apply(as.data.frame(tableDEtime_DEonly), MARGIN=1, FUN=sum)
        DE.nb.G[-which(DE.nb.G%in%c(1, nTimes))] <- nTimes + 10
        DE.nb.G <- as.factor(as.numeric(DE.nb.G))
        IdAttribute <- which(c(1, nTimes, nTimes + 10)%in%levels(DE.nb.G))
        levels(DE.nb.G) <- c("Specific", "Common", "Other")[IdAttribute]

        DE.genes.info <- data.frame(Index=DEgenes_number,
                                    Name=Gnames,
                                    Attribute=DE.nb.G)

        ##-------------------------------------------------------------------##
        ## Data for alluvial graph : Patterns (1 for DE vs t0),
        ## time group & Freq
        alluvData <- cbind(tableDEtime_DEonly,
                           data.frame(Attribute=DE.genes.info$Attribute,
                                      Freq=rep(1, times=nDEgenes)))

        allu2DataTini <- stats::aggregate(x=list(Freq=rep(1, nrow(alluvData))),
                                          by=alluvData[,seq_len(nTimes + 1)],
                                          FUN=length)

        allu2Data.t <- allu2DataTini
        allu2Data.t <- cbind(Gene=seq_len(nrow(allu2Data.t)), allu2Data.t)

        allu2Data.tf <- reshape2::melt(data=allu2Data.t,
                                       id.vars=c("Gene", "Attribute", "Freq"))
        allu2Data.tf[!duplicated(allu2Data.tf), ]

        allu2Data.tf$Attribute <- as.factor(allu2Data.tf$Attribute)
        allu2Data.tf$Gene <- as.factor(allu2Data.tf$Gene)
        allu2Data.tf$variable <- as.numeric(allu2Data.tf$variable)*50

        allu2Data.tf$value <- as.factor(allu2Data.tf$value)
        ## levels(allu2Data.tf$value) <- c("grey", "black")
        allu2valueKeep <- as.numeric(levels(allu2Data.tf$value)) + 1
        levels(allu2Data.tf$value)<-c("no", "yes")[allu2valueKeep]

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 3) graphics parameters: colors, labels ...
        allu2DataTini_apply <- allu2DataTini[,seq_len(nTimes)]
        Id.T.0DE <- which(apply(allu2DataTini_apply, 2, FUN=function(x) 0%in%x))
        Id.T.1DE <- which(apply(allu2DataTini_apply, 2, FUN=function(x) 1%in%x))
        colorbarFILLini <- rep(c('black','grey'), times=nTimes)

        if(length(Id.T.1DE)==0 | length(Id.T.0DE)==0){
            if (length(Id.T.1DE) == 0) {
                colorbarFILLkeep <- 2*(as.numeric(Id.T.0DE) - 1) + 1
            } else {
                colorbarFILLkeep <- 2*as.numeric(Id.T.1DE)
            }## if(length(Id.T.1DE)==0)
        } else {
            colorbarFILLkeep <- sort(c(2*(as.numeric(Id.T.0DE) - 1) + 1,
                                       2*as.numeric(Id.T.1DE)))
        }## if(length(Id.T.1DE)==0 | length(Id.T.0DE)==0)

        colorbarFILL <- colorbarFILLini[colorbarFILLkeep]

        ColorAllu <- allu2Data.tf$Attribute
        levels(ColorAllu)<-c("#0099B4FF", "#42B540FF", "#F39B7FB2")[IdAttribute]
        ## 4DBBD5B2, #00A087B2, #42B54099

        alphaAllu <- allu2Data.tf$Attribute
        levels(alphaAllu) <- c(0.7, 1, 0.6)[IdAttribute]
        alphaAllu <- as.numeric(as.character(alphaAllu))

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## 4) Graphs ## Graph alluvial 1 ## stat_alluvium(color="grey")
        g.alluvial <- ggplot2::ggplot(allu2Data.tf,
                                      ggplot2::aes(x=variable, y=Freq,
                                                   alluvium=Gene,
                                                   fill=Attribute,
                                                   stratum=value)) +
            ggalluvial::stat_alluvium(geom="flow", lode.guidance="forward",
                                      width=5, reverse=FALSE, alpha=alphaAllu) +
            ggplot2::scale_fill_manual(values=levels(ColorAllu)) +
            ggalluvial::stat_stratum(width=5, mapping=ggplot2::aes(size=value),
                                     fill=colorbarFILL,
                                     reverse=FALSE, alpha=0.7) +
            ggplot2::labs(x="", y="Number of DE genes") +
            ggplot2::scale_x_continuous(breaks=seq_len(nTimes)*50,
                                        labels=alluvLabels,
                                        guide=guide_axis(angle=45)) +
            ggplot2::scale_size_manual("DE tx versus t0", values=c(0.5, 0.5),
                                       guide=colorLegend)

        if (is.null(title.alluvial)) {
            g.alluvial <- g.alluvial + ggplot2::ggtitle("Alluvial graph")
        } else {
            g.alluvial <- g.alluvial + ggplot2::ggtitle(title.alluvial)
        }## if(is.null(title.alluvial))

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        ## Output
        return(g.alluvial)
    }## if(isTRUE(Temporal.Group))
}## DEplotAlluvial()
