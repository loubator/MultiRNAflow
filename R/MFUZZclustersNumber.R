#' @title Automatic choice of the number of clusters to use for
#' the Mfuzz analysis
#'
#' @description The function uses
#' [stats::kmeans()] or
#' [FactoMineR::HCPC()]
#' in order to compute the number of cluster for the
#' [Mfuzz::mfuzz()] analysis.
#'
#' @details All results are built from the results of our function
#' [DATAnormalization()].
#'
#' The \code{Mfuzz} package works with datasets where rows correspond to genes
#' and columns correspond to times.
#' If \code{RawCounts} (input of our function
#' [DATAprepSE()])
#' contains several replicates per time,
#' the algorithm computes the mean of replicates for each gene before using
#' [Mfuzz::mfuzz()].
#' When there are several biological conditions, the algorithm realizes
#' the [Mfuzz::mfuzz()]
#' analysis for each biological condition.
#'
#' The kmeans method or the hierarchical clustering method,
#' respectively included in
#' [stats::kmeans()] and
#' [FactoMineR::HCPC()],
#' is used in order to compute the optimal number of clusters.
#' If there are several biological conditions, the algorithm computes
#' one optimal number of clusters per biological condition.
#'
#' @param SEresNorm Results of the function
#' [DATAnormalization()].
#' @param DATAnorm \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' \code{TRUE} means the function uses the normalized data.
#' \code{FALSE} means the function uses the raw counts data.
#' @param Method "kmeans" or "hcpc". The method used for selecting the number
#' of cluster to be used for the temporal cluster analysis (see \code{Details}).
#' \code{Method="kmeans"} is advised for large number of genes.
#' @param Max.clust Integer strictly superior to 1 indicating
#' the maximum number of clusters. The default is \code{Max.clust=10}.
#' @param Min.std Numeric positive value. All genes where their
#' standard deviations are smaller than the threshold Min.std will be excluded.
#' @param Plot.Cluster \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, the output graph will be plotted.
#' Otherwise the graph will be plotted.
#' @param path.result Character or \code{NULL}.
#' Path to save the plot described in the section \code{Value}.
#' If \code{NULL}, the graph will not be saved in a folder.
#' \code{NULL} as default.
#'
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEresNorm} with the different elements below,
#' saved in the metadata \code{Results[[1]][[4]]} of \code{SEresNorm},
#' * the optimal number of clusters for each biological condition
#' (between 2 and \code{Max.clust}).
#' * a data.frame with (\eqn{N_{bc}+1}) columns and \code{Max.clust} rows
#' with \eqn{N_{bc}} the number of biological conditions.
#'   * If \code{Method="kmeans"}, the ith rows and the jth column correspond
#'     to the within-cluster intertia (see \code{tot.withinss} from
#'     [stats::kmeans()])
#'     dividing by the sum of the variance of each row of \code{ExprData}
#'     of the (j-1)th biological condition computed by
#'     [stats::kmeans()]
#'     with i clusters.
#'     When there is only one cluster, the within-cluster intertia
#'     corresponds to the sum of the variance of each row of
#'     \code{ExprData} (see \code{Details}).
#'     The first column contains integers between 1 and \code{Max.clust}
#'     which corresponds to the number of clusters selected for the
#'     [stats::kmeans()]
#'     analysis.
#'   * If \code{Method="hcpc"}, the jth column correspond to the clustering
#'   heights (see the output \code{height} from
#'   [FactoMineR::HCPC()])
#'   dividing by the maximum value of \code{height}.
#'   The first column contains integers between 1 and \code{Max.clust}
#'   which corresponds to the number of clusters selected for the
#'   [stats::kmeans()]
#'   analysis.
#' * a plot which gives
#'   * If \code{Method="kmeans"}, the evolution of the weighted
#'   within-cluster intertia per number of clusters
#'   (from 1 to \code{Max.clust}) for each biological condition.
#'   The optimal number of cluster for each biological condition
#'   will be colored in blue.
#'   * If \code{Method="hcpc"}, the evolution of the scaled height per
#'   number of clusters (from 1 to \code{Max.clust})
#'   for each biological condition.
#'   The optimal number of cluster for each biological condition will be
#'   colored in blue.
#'
#' @seealso The function is called by
#' [MFUZZanalysis()].
#'
#' @importFrom stats kmeans
#' @importFrom ggplot2 ggplot aes geom_line geom_point ylim scale_color_manual
#' ylab guides guide_legend theme
#' @importFrom FactoMineR HCPC
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors metadata
#' @importFrom graphics lines legend
#' @importFrom grDevices pdf dev.off
#'
#' @export
#'
#' @examples
#' ## Data simulation
#' set.seed(33)
#' DATAclustSIM <- matrix(rnorm(12*10*3, sd=0.2,
#'                              mean=rep(c(rep(c(1, 6, 9, 4, 3, 1,
#'                                               6.5, 0.7, 10), times=2),
#'                                         rep(c(2, 3.6, 3.7, 5, 7.9, 8,
#'                                               7.5, 3.5, 3.4), times=2)),
#'                                       each=10)),
#'                        nrow=30, ncol=12)
#' DATAclustSIM <- floor(DATAclustSIM*100)
#' ##
#' colnames(DATAclustSIM) <- c("G1_t0_r1", "G1_t1_r1", "G1_t2_r1",
#'                             "G1_t0_r2", "G1_t1_r2", "G1_t2_r2",
#'                             "G2_t0_r3", "G2_t1_r3", "G2_t2_r3",
#'                             "G2_t0_r4", "G2_t1_r4", "G2_t2_r4")
#' ##------------------------------------------------------------------------##
#' ## Plot the temporal expression of each individual
#' graphics::matplot(t(rbind(DATAclustSIM[, 1:3], DATAclustSIM[, 4:6],
#'                           DATAclustSIM[, 7:9], DATAclustSIM[, 10:12])),
#'                   col=rep(c("black", "red"), each=6*10),
#'                   xlab="Time", ylab="Gene expression", type=c("b"), pch=19)
#'
#' ##------------------------------------------------------------------------##
#' ## Preprocessing step
#' DATAclustSIM <- data.frame(DATAclustSIM)
#'
#' resDATAprepSE <- DATAprepSE(RawCounts=DATAclustSIM,
#'                             Column.gene=NULL,
#'                             Group.position=1,
#'                             Time.position=2,
#'                             Individual.position=3)
#' ## Normalization
#' resNorm <- DATAnormalization(SEres=resDATAprepSE,
#'                              Normalization="rle",
#'                              Plot.Boxplot=FALSE,
#'                              Colored.By.Factors=FALSE)
#'
#' ##------------------------------------------------------------------------##
#' resMFUZZcluster <- MFUZZclustersNumber(SEresNorm=resNorm,
#'                                        DATAnorm=FALSE,
#'                                        Method="hcpc",
#'                                        Max.clust=5,
#'                                        Plot.Cluster=TRUE,
#'                                        path.result=NULL)

MFUZZclustersNumber <- function(SEresNorm,
                                DATAnorm=TRUE,
                                Method="hcpc",
                                Max.clust=3,
                                Min.std=0.1,
                                Plot.Cluster=TRUE,
                                path.result=NULL) {
    ##-----------------------------------------------------------------------##
    ## To avoid "no visible binding for global variable" with devtools::check()
    Nb.clust <- value <- variable <- NULL

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check
    resErr <- ErrMFUZZclustersNumber(SEresNorm,
                                     DATAnorm=DATAnorm,
                                     Method=Method,
                                     Max.clust=Max.clust,
                                     Min.std=Min.std,
                                     Plot.Cluster=Plot.Cluster,
                                     path.result=path.result)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Preprocessing
    cSEdat <- SummarizedExperiment::colData(SEresNorm)

    if (c("Time")%in%colnames(cSEdat)) {
        Vect.time <- as.character(cSEdat$Time)
        timeLevels <- levels(as.factor(Vect.time))
        Nb.time <- length(timeLevels)
    } else {
        stop("Samples must belong to different times points.")
    }## if (c("Time")%in%colnames(cSEdat))

    if (c("Group")%in%colnames(cSEdat)) {
        Vect.group <- as.character(cSEdat$Group)
        groupLevels <- levels(as.factor(Vect.group))
        Nb.group <- length(groupLevels)

        colname.grp <- paste0(".", rep(groupLevels, each=Nb.time))
        ClustOptTitle <- paste0("Mfuzz_OptimalClusterNumber_",
                                paste0(groupLevels, collapse="_"), ".pdf")
    } else {
        Vect.group <- NULL
        groupLevels <- NULL
        Nb.group <- 1

        colname.grp <- ""
        ClustOptTitle <- "Mfuzz_OptimalClusterNumber.pdf"
    }## if (c("Group")%in%colnames(cSEdat))

    if (isTRUE(DATAnorm)) {
        aSE <- 2
    } else {
        aSE <- 1
    }## if (isTRUE(DATAnorm))

    ExprData <- SummarizedExperiment::assays(SEresNorm)[[aSE]]
    ExprData.f <- data.frame(ExprData)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Data for Mfuzz analysis and selection of the number of cluster
    mfuzzData <- matrix(NA, ncol=Nb.time*Nb.group, nrow=nrow(ExprData.f))
    row.names(mfuzzData) <- row.names(ExprData.f)

    Tps.info <- paste0("t", gsub("T", "", gsub("t", "", timeLevels)))

    colnames(mfuzzData) <- paste0("Mean_", rep(Tps.info, times=Nb.group),
                                  colname.grp)

    ##-----------------------------------------------------------------------##
    ## Filling the data
    for (g in seq_len(Nb.group)) {
        for (t in seq_len(Nb.time)) {
            Index.t <- which(Vect.time == timeLevels[t])

            if (is.null(Vect.group)) {
                mfuzzData[, t] <- apply(as.data.frame(ExprData.f[, Index.t]),
                                        1, mean)
            } else {
                Index.g <- which(Vect.group == groupLevels[g])
                Index.tg <- intersect(Index.t, Index.g)

                DataTG <- as.data.frame(ExprData.f[, Index.tg])
                mfuzzData[, Nb.time*(g-1) + t] <- apply(DataTG, 1, FUN=mean)
            }## if (is.null(Vect.group))
        }## for(t in seq_len(Nb.time))
    }## for(g in seq_len(Nb.group))

    ##-----------------------------------------------------------------------##
    ## Data which will contain the results of Kmeans
    Sum.nb.c <- data.frame(matrix(NA, nrow=Max.clust, ncol=Nb.group+1))

    if (Method == "hcpc") {
        Score <- "Tot.withinss.scaled"
    } else {
        Score <- "Scaled.height"
    }## if(Method == "hcpc")

    if (is.null(Vect.group) == TRUE) {
        colnames(Sum.nb.c) <- c("Nb.clust", as.character(Score))
    } else {
        colnames(Sum.nb.c) <- c("Nb.clust",
                                paste0(as.character(Score), "_", groupLevels))
    }## if(is.null(Vect.group)==TRUE)

    Sum.nb.c[, 1] <- c(1, seq(2, Max.clust))

    ##-----------------------------------------------------------------------##
    ## Kmeans
    clustOPT <- rep(NA, times=Nb.group)

    for (g in seq_len(Nb.group)) {
        ##-------------------------------------------------------------------##
        Std.g <- apply(mfuzzData[,seq_len(Nb.time) + Nb.time*(g-1)], 1, sd)
        GeneInf.Min.std.g <- which(Std.g < Min.std)

        if (length(GeneInf.Min.std.g) > 0) {
            GeneSelNbClust <- -GeneInf.Min.std.g
        } else {
            GeneSelNbClust <- seq_len(length(Std.g))
        }## if(length(GeneInf.Min.std.g)>0)

        mfuzzData.g <- data.frame(mfuzzData[GeneSelNbClust,
                                            seq_len(Nb.time)+Nb.time*(g-1)])

        ##-------------------------------------------------------------------##
        if (Method == "hcpc") {
            Nb.gene <- nrow(mfuzzData.g)

            if (Nb.gene <= 200) {
                res.pca <- FactoMineR::PCA(round(mfuzzData.g, digits=3),
                                           graph=FALSE)
                res.hcpc <- FactoMineR::HCPC(res.pca, graph=FALSE,
                                             nb.clust=-1, consol=TRUE, min=2)
            } else {
                res.kk <- NbClustKmeansHCPC(Nb.gene, 200, 50, 30000, 175)
                kkHCPC <- res.kk$Nkmeans

                options(warn=-1)

                cl <- stats::kmeans(round(data.frame(scale(mfuzzData.g)),
                                          digits=2),
                                    kkHCPC, iter.max=10)
                res.hcpc <- FactoMineR::HCPC(round(data.frame(cl$centers),
                                                   digits=2),
                                             graph=FALSE, nb.clust=-1,
                                             consol=FALSE, min=3)

                options(warn=0)
            }## if(nrow(mfuzzData.g)<=200)

            ## rev(res.hc$height)
            Clust.height <- rev(res.hcpc$call$t$tree$height)
            Sum.nb.c[, g+1] <- c(Clust.height/max(Clust.height),
                                 0)[seq_len(Max.clust)]

            inert.gain <- rev(res.hcpc$call$t$tree$height)
            intra <- rev(cumsum(rev(inert.gain)))
            quot <- intra[2:length(intra)]/intra[seq_len(length(intra)-1)]

            if (abs(which.min(quot) + 1 - res.hcpc$call$t$nb.clust) >2 ) {
                Index.nb.clust <- res.hcpc$call$t$nb.clust
            } else {
                Index.nb.clust <- max(which.min(quot) + 1,
                                      res.hcpc$call$t$nb.clust)
                ## Index.nb.clust<-res.hcpc$call$t$nb.clust
            }## if(abs(which.min(quot)+1-res.hcpc$call$t$nb.clust)>2)

        }## if(Method=="hcpc")

        ##-------------------------------------------------------------------##
        if (Method == "kmeans") {
            inertie.wihtin <- rep(0, times=Max.clust-1)
            cpt.clust <- 0

            options(warn = -1)

            for (k in seq(from=2, to=Max.clust, by=1)) {
                cpt.clust <- cpt.clust + 1
                clus <- stats::kmeans(round(data.frame(scale(mfuzzData.g)),
                                            digits=1),
                                      centers=k, nstart=5)
                inertie.wihtin[cpt.clust] <- clus$tot.withinss
            }## for(k in 2:Max.clust)

            options(warn = 0)

            ## DQ.within<-rev(cumsum(rev(inertie.wihtin)))
            ## quot.DQ.within<-DQ.within[-length(inertie.wihtin)]/DQ.within[-1]
            DQ.within <- rev(cumsum(rev(c(clus$totss, inertie.wihtin))))
            quot.DQ.within <- DQ.within[-1]/DQ.within[-length(inertie.wihtin)]
            Index.nb.clust <- which.min(quot.DQ.within) + 1

            Sum.nb.c[, g + 1] <- c(clus$totss, inertie.wihtin)/clus$totss
        }## if(Method=="kmeans")

        ##-------------------------------------------------------------------##
        clustOPT[g] <- Index.nb.clust
    }## for(g in seq_len(Nb.group))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## cluster graph preprocess
    NBsumCl <- Sum.nb.c

    if (Method == "hcpc") {
        Ylab <- "Scaled height (ward)"
    } else {
        Ylab <- "Scaled within-cluster inertia"
    }## if (Method == "hcpc")

    if (Nb.group > 1) {
        for (g in seq_len(Nb.group-1)) {
            NBsumCl[, g + 2] <- NBsumCl[, g + 2] + 0.1*g
        }## for (g in seq_len(Nb.group-1))

        varLevels <- groupLevels
        datOPT <- data.frame(Nb.clust=clustOPT,
                             value=diag(as.matrix(NBsumCl[clustOPT, -1])),
                             variable=rep("Optimal cluster", times=Nb.group))
    } else {
        varLevels <- "G1"
        datOPT <- data.frame(Nb.clust=clustOPT,
                             value=as.numeric(NBsumCl[clustOPT, 2]),
                             variable="Optimal cluster")
    }## if (Nb.group > 1)

    meltCL <- reshape2::melt(NBsumCl, id.vars=1)
    meltCL$variable <- as.factor(meltCL$variable)
    levels(meltCL$variable) <- varLevels

    ##-----------------------------------------------------------------------##
    ## cluster graph
    ggNBcl <- ggplot2::ggplot(meltCL,
                              ggplot2::aes(x=Nb.clust, y=value,
                                           group=variable)) +
        ggplot2::geom_line(ggplot2::aes(linetype=variable)) +
        ggplot2::geom_point() + ## ggplot2::ylim(0, 1 + 0.15*(Nb.group-1)) +
        ggplot2::geom_point(data=datOPT,
                            ggplot2::aes(x=Nb.clust, y=value,
                                         group=variable, color=variable),
                            shape=18, size=4) +
        ggplot2::scale_color_manual(values=c("blue"), guide="none") +
        ggplot2::ylab(Ylab) + ggplot2::xlab("Number of clusters")

    if (Nb.group > 1) {
        ggNBttl <- "Biological conditions"
        ggNBcl <- ggNBcl +
            ggplot2::guides(linetype=ggplot2::guide_legend(order=1,
                                                           title=ggNBttl),
                            color=ggplot2::guide_legend(order=2, title=""))
    } else {
        ggNBcl <- ggNBcl +
            ggplot2::guides(linetype="none",
                            color=ggplot2::guide_legend(order=2, title=""))
    }## if (Nb.group > 1)

    ggNBcl <- ggNBcl+
        ggplot2::theme(legend.position="bottom", legend.box="horizontal")

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Save and plot graph
    if (!is.null(path.result)) {
        grDevices::pdf(file=file.path(path.result, ClustOptTitle),
                       width=11, height=8)
        print(ggNBcl)
        grDevices::dev.off()
    }## if(!is.null(path.result))

    if (isTRUE(Plot.Cluster)) {
        print(ggNBcl)
    }## if (isTRUE(Plot.Cluster))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Data containing the number of cluster for each group
    if (is.null(Vect.group)) {
        clustData <- data.frame(Name="OneGroupOnly", ClusterKmeans=clustOPT)
    } else {
        clustData <- data.frame(Name=groupLevels, ClusterKmeans=clustOPT)
    }## if(is.null(Vect.group)==TRUE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## SE Mfuzz
    listMFUZZ <- list(Summary.Nb.Cluster=Sum.nb.c,
                      DataClustSel=clustData,
                      clustDATAplot=ggNBcl,
                      MfuzzData=mfuzzData)
    SEprepMFUZZ <- SEresNorm
    S4Vectors::metadata(SEprepMFUZZ)$Results[[1]][[4]] <- listMFUZZ

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    return(SEobj=SEprepMFUZZ)
}## MFUZZclustersNumber()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

NbClustKmeansHCPC <- function(NrowData, x1, y1, x2, y2){
    if (x2 <= x1 | y2 <= y1) {
        StopMfuzzMessage <- paste0("'x2' must be strictly greater than 'x1'",
                                   " and ",
                                   "'y2' must be strictly greater than 'y1'")
        stop(StopMfuzzMessage)
    }## if(x2 <= x1 | y2 <= y1)

    acoef <- log(y2/y1)/log(x2/x1)
    bcoef <- y2/(x2^acoef)
    NbClustKKhcpc <- ceiling(bcoef*NrowData^acoef)

    return(list(Nkmeans=NbClustKKhcpc,
                a=acoef,
                b=bcoef))
}## NbClustKmeansHCPC()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrMFUZZclustersNumber <- function(SEresNorm,
                                   DATAnorm=TRUE,
                                   Method="hcpc",
                                   Max.clust=3,
                                   Min.std=0.1,
                                   Plot.Cluster=TRUE,
                                   path.result=NULL) {
    ##-----------------------------------------------------------------------##
    ## Check SEresNorm (cf DATAplotExpressionGenes())
    Err1_MfuzzNcl <- ErrSEresNorm(SEresNorm=SEresNorm,
                                  DATAnorm=DATAnorm,
                                  path.result=path.result)

    Err2_MfuzzNcl <- ErrNNI(NNI=Max.clust, NNIname="Max.clust")

    Ngenes <- length(S4Vectors::rownames(SEresNorm))
    if(Max.clust<2 | floor(Max.clust)!=Max.clust | Max.clust>=Ngenes){
        ## length(S4Vectors::rownames(resDATAnormFission))
        Err_max <- paste0("'Max.clust' must be an integer ",
                          "greater or equal to 2 and ",
                          "lesser than the number of genes.")
        stop(Err_max)
    }## if(Max.clust<2 | floor(Max.clust)!=Max.clust | Max.clust>=Ngenes)

    ##-----------------------------------------------------------------------##
    if (!is.numeric(Min.std)) {
        stop("'Min.std' must be positive numeric values")
    }## if (!is.numeric(Min.std))

    if (Min.std < 0) {
        stop("'Min.std' must be positive numeric values")
    }## if (Min.std < 0)

    ##-----------------------------------------------------------------------##
    ## Plot.Cluster
    if (!isTRUE(Plot.Cluster) & !isFALSE(Plot.Cluster)) {
        stop("'Plot.Cluster' must be TRUE or FALSE.")
    }## if (!isTRUE(Plot.Cluster) & !isFALSE(Plot.Cluster))

    ## Different Method
    if (!Method%in%c("hcpc", "kmeans")) {
        stop("'Method' must be 'hcpc' or 'kmeans'.")
    }## if (!Method%in%c("hcpc", "kmeans"))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrMFUZZclustersNumber()

