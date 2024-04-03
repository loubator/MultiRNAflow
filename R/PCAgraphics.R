#' @title Multiple 2D and 3D PCA graphs.
#'
#' @description The function plots 2D and 3D PCA using the function
#' [PCArealization()]
#' which realizes a PCA analysis. This function is called repeatedly by the
#' function [PCAanalysis()]
#' if samples belong to different biological conditions and time points.
#'
#' @details All results are built from the results of our function
#' [DATAnormalization()].
#'
#' @param SEresNorm Results of the function
#' [DATAnormalization()].
#' @param DATAnorm \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' \code{TRUE} means the function uses the normalized data.
#' \code{FALSE} means the function uses the raw counts data.
#' @param gene.deletion \code{NULL} or a vector of characters or a vector of
#' integers. \code{NULL} as default.
#' If \code{gene.deletion} is a vector of characters, all genes with names in
#' \code{gene.deletion} will be deleted from the data set as input
#' \code{RawCounts} of our function
#' [DATAprepSE()].
#' If \code{gene.deletion} is a vector of integers,
#' all the corresponding row numbers will be deleted from the data set as input
#' \code{RawCounts} of our function
#' [DATAprepSE()].
#' If \code{gene.deletion=NULL} all genes will be used in the construction
#' of the PCA.
#' @param sample.deletion \code{NULL} or a vector of characters or
#' a vector of integers. \code{NULL} as default.
#' If \code{sample.deletion} is a vector of characters, all samples with names
#' in \code{sample.deletion} will not be used in the construction of the PCA.
#' If \code{sample.deletion} is a vector of integers,
#' all the corresponding column numbers will not be used in the construction
#' of the PCA from the data set as input \code{RawCounts} of our function
#' [DATAprepSE()].
#' If \code{sample.deletion=NULL} all samples will be used
#' in the construction of the PCA.
#' @param Plot.PCA \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, PCA graphs will be plotted.
#' Otherwise no graph will be plotted.
#' @param Mean.Accross.Time \code{TRUE} or \code{FALSE}.
#' \code{FALSE} as default.
#' If \code{FALSE} and if \code{Time.position}
#' (input of [DATAprepSE()])
#' is not set as \code{NULL}, consecutive time points within a sample
#' are linked to help visualization of temporal patterns.
#' If \code{TRUE} and if \code{Time.position} is not set as \code{NULL},
#' the mean per time of all genes is computed for each biological condition and
#' the means of consecutive time points within biological condition are
#' linked to help visualization of temporal patterns.
#' @param Color.Group \code{NULL} or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' If \code{Color.Group}
#' (input of [DATAprepSE()])
#' is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute
#' a color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param motion3D \code{TRUE} or \code{FALSE}.
#' If TRUE, the 3D PCA plots will also be plotted in a rgl window
#' (see [plot3Drgl::plotrgl()])
#' allowing to interactively rotate and zoom.
#' @param Phi Angle defining the colatitude direction for the 3D PCA plot
#' (see \code{Details} in
#' [graphics::persp()]).
#' @param Theta Angle defining the azimuthal direction for the 3D PCA plot
#' (see \code{Details} in
#' [graphics::persp()]).
#' @param Cex.point Non negative numeric value giving the size of points
#' in all PCA plots which are not automatically plotted by
#' [FactoMineR::PCA()].
#' @param Cex.label Non negative numeric value giving the size of
#' the labels associated to each point of the all PCA graphs
#' which are not automatically plotted by
#' [FactoMineR::PCA()].
#' @param epsilon Non negative numeric value giving the length between points
#' and their labels in all PCA plots which are not automatically plotted
#' by [FactoMineR::PCA()].
#' @param path.result Character or \code{NULL}.
#' Path to save the different PCA graphs.
#' If \code{NULL}, the different PCA graphs will not be saved in a folder.
#' \code{NULL} as default.
#' @param Name.file.pca Character or \code{NULL}.
#' If \code{Name.file.pca} is a character, \code{Name.file.pca} will be added
#' at the beginning of all names of the saved graphs.
#'
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEresNorm} with the outputs from the function
#' [FactoMineR::PCA()],
#' and plots several 2D and 3D PCA graphs depending
#' on the experimental design (if \code{Plot.PCA=TRUE}),
#' saved in the metadata \code{Results[[1]][[2]]} of \code{SEresNorm},
#' * When samples belong only to different biological conditions,
#' the function returns a 2D and two 3D PCA graphs.
#' In each graph, samples are colored with different colors for different
#' biological conditions. The two 3D PCA graphs are identical but
#' one of them will be opened in a rgl window
#' (see [plot3Drgl::plotrgl()])
#' and it allows to interactively rotate and zoom.
#' * When samples belong only to different time points, the function returns
#'   * One 2D PCA graph, one 3D PCA graph and the same 3D PCA graph
#'   in a rgl window where samples are colored with different colors
#'   for different time points.
#'   Furthermore, lines are drawn between each pair of consecutive points
#'   for each sample (if \code{Mean.Accross.Time=FALSE},
#'   otherwise it will be only between means).
#'   * The same graphs describe above but without lines.
#' * When samples belong to different time points and different
#' biological conditions, the function returns
#'   * One 2D PCA graph, one 3D PCA graph and the same 3D PCA graph
#'   in a rgl window where samples are colored with different colors
#'   for different time points.
#'   Furthermore, lines are drawn between each pair of consecutive points
#'   for each sample (if \code{Mean.Accross.Time=FALSE},
#'   otherwise it will be only between means).
#'   * The same graphs describe above but without lines.
#'   * The same six following graphs for each biological condition
#'   (one PCA analysis per biological condition).
#'   One 2D PCA graph, one 3D PCA graph and the same 3D PCA graph
#'   in a rgl window where samples belong to only one biological condition
#'   and are colored with different colors for different time points.
#'   Furthermore, lines are drawn between each pair of consecutive points
#'   for each sample (if \code{Mean.Accross.Time=FALSE},
#'   otherwise it will be only between means).
#'   The three others graphs are identical to the three previous ones
#'   but without lines.
#'
#' The interactive 3D graphs will be plotted only if \code{motion3D=TRUE}.
#'
#' @seealso This function is called by our function
#' [PCAanalysis()]
#' and calls our function
#' [PCArealization()].
#'
#' @importFrom S4Vectors metadata
#' @importFrom stats aggregate
#' @importFrom FactoMineR plot.PCA
#' @importFrom plot3D scatter3D text3D
#' @importFrom ggplotify as.ggplot
#' @importFrom graphics legend lines text
#' @importFrom grDevices dev.off pdf recordPlot
#' @importFrom plot3Drgl plotrgl
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
#' ## Normalization
#' resNorm <- DATAnormalization(SEres=resDATAprepSE,
#'                              Normalization="rle",
#'                              Plot.Boxplot=FALSE,
#'                              Colored.By.Factors=FALSE)
#' ## Color for each group
#' GROUPcolor <- data.frame(Name=c("G1", "G2"), Col=c("black", "red"))
#' ##------------------------------------------------------------------------##
#' resPCAgraph <- PCAgraphics(SEresNorm=resNorm,
#'                            DATAnorm=TRUE,
#'                            gene.deletion=c("Gene1", "Gene5"),
#'                            sample.deletion=c(2,6),
#'                            Plot.PCA=TRUE,
#'                            Mean.Accross.Time=FALSE,
#'                            Color.Group=GROUPcolor,
#'                            motion3D=FALSE,
#'                            Phi=25, Theta=140, Cex.label=0.7,
#'                            Cex.point=0.7, epsilon=0.2,
#'                            path.result=NULL, Name.file.pca=NULL)

PCAgraphics <- function(SEresNorm,
                        DATAnorm=TRUE,
                        gene.deletion=NULL,
                        sample.deletion=NULL,
                        Plot.PCA=TRUE,
                        Mean.Accross.Time=FALSE,
                        Color.Group=NULL,
                        motion3D=FALSE,
                        Phi=25, Theta=140, epsilon=0.2,
                        Cex.point=0.7, Cex.label=0.7,
                        path.result=NULL,
                        Name.file.pca=NULL) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check inputs
    resErr <- ErrPCAgraphics(Plot.PCA=Plot.PCA,
                             Mean.Accross.Time=Mean.Accross.Time,
                             motion3D=motion3D,
                             Phi=Phi, Theta=Theta, epsilon=epsilon,
                             Cex.point=Cex.point, Cex.label=Cex.label,
                             path.result=path.result,
                             Name.file.pca=Name.file.pca) ## Color.Group=NULL

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## PCA
    SEresPCA <- PCArealization(SEresNorm=SEresNorm,
                               DATAnorm=DATAnorm,
                               sample.deletion=sample.deletion,
                               gene.deletion=gene.deletion,
                               Supp.del.sample=FALSE)

    PCAlist <- S4Vectors::metadata(SEresPCA)$Results[[1]][[2]]
    listFCTRS <- PCAlist$List.Factors
    res.PCA <- PCAlist$PCAresults

    Vector.time <- listFCTRS$Vector.time
    Vector.group <- listFCTRS$Vector.group
    Vector.patient <- listFCTRS$Vector.patient
    LvlsPAT <- levels(factor(Vector.patient))

    PCAqualiSup <- res.PCA$call$quali.sup$quali.sup
    coordPCAind <- res.PCA$ind$coord

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Color assignment
    if (!is.null(Vector.time) & is.null(Vector.group)) {
        Tlevels <- as.character(levels(factor(Vector.time)))
        Tlevels <- paste0("t", gsub("t", "", gsub("T", "", Tlevels)))
        levels(PCAqualiSup$Quali.Sup.Time)<-Tlevels

        ##-------------------------------------------------------------------##
        TlevelsCOLOR <- myPaletteT(Nt=length(Tlevels))
        Color.Time <- data.frame(Name=Tlevels, Col=TlevelsCOLOR)

        data.color <- rbind(Color.Time, rep(NA, times=ncol(Color.Time)))
        data.color <- data.color[-which(is.na(data.color$Col)),]

        ind.color <- PCAqualiSup$Quali.Sup.Time
        levels(ind.color) <- data.color[, 2]
        LegendTitle <- "Time"

        ##-------------------------------------------------------------------##
        ## Color.Time <- NULL
        ## if (is.null(Color.Time)) {
        ##     TlevelsCOLOR <- c("#737373",## "#252525"
        ##                       scales::hue_pal()(length(Tlevels)-1))
        ##     Color.Time <- data.frame(Name=Tlevels, Col=TlevelsCOLOR)
        ##     data.color <- rbind(Color.Time,
        ##                         rep(NA, times=ncol(Color.Time)))
        ##     data.color <- data.color[-which(is.na(data.color$Col)),]
        ## } else {
        ##     Id.LevelColT <- order(Color.Time[, 1])
        ##     Color.Time <- data.frame(Name=Tlevels,
        ##                              Col=Color.Time[Id.LevelColT, 2])
        ##     data.color <- rbind(Color.Time, rep(NA,times=ncol(Color.Time)))
        ##     data.color <- data.color[-which(is.na(data.color$Col)),]
        ## }## if(is.null(Color.Time))

        ##-------------------------------------------------------------------##
    } else {
        Glevels <- levels(factor(Vector.group))

        if (is.null(Color.Group)) {
            MypaletteG <- myPaletteBC(Nbc=length(Glevels))
            Color.Group <- data.frame(Name=Glevels, Col=MypaletteG)

            data.color <- rbind(Color.Group,
                                rep(NA, times=ncol(Color.Group)))
            data.color <- data.color[-which(is.na(data.color$Col)),]
        } else {
            Id.LevelCol.G <- order(Color.Group[, 1])
            Color.Group <- data.frame(Name=Glevels,
                                      Col=Color.Group[Id.LevelCol.G, 2])

            data.color <- rbind(Color.Group,
                                rep(NA, times=ncol(Color.Group)))
            data.color <- data.color[-which(is.na(data.color$Col)),]
        }## if(is.null(Color.Group)==TRUE)

        ind.color <- PCAqualiSup$Quali.Sup.Group
        levels(ind.color) <- data.color[, 2]
        LegendTitle <- "Group"
    }## if(!is.null(Vector.time) & is.null(Vector.group))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (is.null(Name.file.pca)) {
        Name.file.pca.f <- "Graph"
    } else {
        Name.file.pca.f <- Name.file.pca
    }## if(is.null(Name.file.pca)==TRUE)

    TitlePCA2D <- paste0(Name.file.pca.f, "_PCA2D.pdf")
    TitlePCA3D <- paste0(Name.file.pca.f, "_PCA3D.pdf")

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## PCA plot from plot.PCA() ## Qualitative factor are in magenta
    options(ggrepel.max.overlaps=30)

    ggPCA2d <- FactoMineR::plot.PCA(res.PCA, axes=c(1, 2),
                                    choix="ind", habillage="ind",
                                    col.hab=as.character(ind.color),
                                    col.quali="magenta",
                                    shadowtext=TRUE, autoLab="y", title="",
                                    graph.type="ggplot", cex=Cex.label)

    options(ggrepel.max.overlaps=10)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 3D PCA colored by biological condition or time points
    ## par(mar=c(4.1, 4.1, 2.1, 2.1))#, xpd=TRUE)
    data.3D <- coordPCAind[, c(1, 2, 3)]

    ggPCA3d <- ggplotify::as.ggplot(
        function() PCAplot3D(PCAcoord3D=data.3D, matPCAeig=res.PCA$eig,
                             colorINDfactor=ind.color,
                             clabPCA=c("", data.color$Name), epsilon=0.2,
                             Phi=25, Theta=140, Cex.point=0.7, Cex.label=0.7)
    )

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (!is.null(path.result)) {
        grDevices::pdf(file=file.path(path.result, TitlePCA2D),
                       width=11, height=8)
        print(ggPCA2d)
        grDevices::dev.off()

        grDevices::pdf(file=file.path(path.result, TitlePCA3D),
                       width=11, height=8)
        print(ggPCA3d)
        grDevices::dev.off()
    }## if(!is.null(path.result))

    if (isTRUE(Plot.PCA)) {
        ## graphics::plot.new() ## clean up device
        print(ggPCA2d)

        ## graphics::plot.new() ## clean up device
        print(ggPCA3d)
        if (isTRUE(motion3D)) { ## 3D PCA in rgl windows
            plot3Drgl::plotrgl()
        }## if(isTRUE(motion3D))

    }## if(Plot.PCA==TRUE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (!is.null(Vector.time)) {
        ##-------------------------------------------------------------------##
        coord.t <- PCAqualiSup$Quali.Sup.Time
        ## n.row <- length(coord.t)
        ## index.order <- seq_len(n.row)
        ## nb.time <- length(unique(coord.t))

        data.name <- data.frame(rep("Mean", times=nrow(PCAqualiSup)),
                                PCAqualiSup)

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        resPCA2Dt <- PCAplot2Dt(PCAcoord2D=coordPCAind,
                                matPCAeig=res.PCA$eig,
                                coord.t=coord.t,
                                data.name=data.name,
                                Vector.group=Vector.group,
                                Vector.time=Vector.time,
                                Vector.patient=Vector.patient,
                                Mean.Accross.Time=Mean.Accross.Time,
                                colorINDfactor=ind.color,
                                epsilon=epsilon,
                                Cex.point=Cex.point,
                                Cex.label=Cex.label)

        resPCA3Dt <- ggplotify::as.ggplot(
            function() PCAplot3Dt(PCAcoord3D=coordPCAind,
                                  matPCAeig=res.PCA$eig,
                                  coord.t=coord.t,
                                  data.name=data.name,
                                  Vector.group=Vector.group,
                                  Vector.time=Vector.time,
                                  Vector.patient=Vector.patient,
                                  Mean.Accross.Time=Mean.Accross.Time,
                                  colorINDfactor=ind.color,
                                  epsilon=epsilon,
                                  Phi=Phi, Theta=Theta,
                                  Cex.point=Cex.point,
                                  Cex.label=Cex.label)
        )

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        TitlePCA2Dlk <- paste0(Name.file.pca.f, "_PCA2D_linked.pdf")
        TitlePCA3Dlk <- paste0(Name.file.pca.f, "_PCA3D_linked.pdf")

        if (!is.null(path.result)) {
            grDevices::pdf(file=file.path(path.result, TitlePCA2Dlk),
                           width=11, height=8)
            print(resPCA2Dt)
            grDevices::dev.off()

            grDevices::pdf(file=file.path(path.result, TitlePCA3Dlk),
                           width=11, height=8)
            print(resPCA3Dt)
            grDevices::dev.off()

        }## if(is.null(path.result)==FALSE)


        if (isTRUE(Plot.PCA)) {
            ## graphics::plot.new() ## clean up device
            print(resPCA2Dt)

            ##graphics::plot.new() ## clean up device
            print(resPCA3Dt)

            if (isTRUE(motion3D)) {
                plot3Drgl::plotrgl()
            }## if (isTRUE(motion3D))

        }## if(Plot.PCA==TRUE)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        listPCAplots <- vector(mode="list", length=4)
        names(listPCAplots) <- c("PCA_2D", "PCA_3D",
                                 "PCA_2DtemporalLinks",
                                 "PCA_3DtemporalLinks")
        listPCAplots[[1]] <- ggPCA2d
        listPCAplots[[2]] <- ggPCA3d
        listPCAplots[[3]] <- resPCA2Dt
        listPCAplots[[4]] <- resPCA3Dt
    } else {
        listPCAplots <- vector(mode="list", length=2)
        names(listPCAplots) <- c("PCA_2D", "PCA_3D")
        listPCAplots[[1]] <- ggPCA2d
        listPCAplots[[2]] <- ggPCA3d
    }## if(is.null(Vector.time) == FALSE)

    SEresPCAgraph <- SEresPCA
    PCAlist <- append(PCAlist, listPCAplots)
    S4Vectors::metadata(SEresPCAgraph)$Results[[1]][[2]] <- PCAlist

    ## grid::par(mar=c(5.1, 4.1, 4.1, 2.1))
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(SEobj=SEresPCAgraph)
}## PCAgraphics()

##---------------------------------------------------------------------------##
##--------------------ErrPCAgraphics-----------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrPCAgraphics <- function(Plot.PCA=TRUE,
                           Mean.Accross.Time=FALSE, ##Color.Group=NULL,
                           motion3D=FALSE,
                           Phi=25, Theta=140, epsilon=0.2,
                           Cex.point=0.7, Cex.label=0.7,
                           path.result=NULL,
                           Name.file.pca=NULL) {
    ##-----------------------------------------------------------------------##
    ## Check TRUE, FALSE
    if (!isTRUE(Plot.PCA) & !isFALSE(Plot.PCA)) {
        stop("'Plot.PCA' must be TRUE or FALSE.")
    }## if (!isTRUE(Plot.PCA) & !isFALSE(Plot.PCA))

    if (!isTRUE(Mean.Accross.Time) & !isFALSE(Mean.Accross.Time)) {
        stop("'Mean.Accross.Time' must be TRUE or FALSE.")
    }## if (!isTRUE(Mean.Accross.Time) & !isFALSE(Mean.Accross.Time))

    ##-----------------------------------------------------------------------##
    Err_PCAHCPCgraph <- ErrPCAHCPCgraphics(motion3D=motion3D, epsilon=epsilon,
                                           Phi=Phi, Theta=Theta,
                                           Cex.point=Cex.point,
                                           Cex.label=Cex.label,
                                           path.result=path.result)

    ##-----------------------------------------------------------------------##
    if (!is.null(Name.file.pca)) {
        if (!is.character(Name.file.pca)) {
            stop("'Name.file.pca' must be NULL or a character.")
        }## if (!is.character(Name.file.pca))
    }## if (!is.null(Name.file.pca))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrPCAgraphics()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##--------------------ErrPCAHCPCgraphics-------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrPCAHCPCgraphics <- function(motion3D=FALSE,
                               Phi=25, Theta=140, epsilon=0.2,
                               Cex.point=0.7, Cex.label=0.7,
                               path.result=NULL) {
    ##-----------------------------------------------------------------------##
    if (!isTRUE(motion3D) & !isFALSE(motion3D)) {
        stop("'motion3D' must be TRUE or FALSE.")
    }## if (!isTRUE(motion3D) & !isFALSE(motion3D))

    ##-----------------------------------------------------------------------##
    ## Check numeric
    Err_a <- "'Phi', 'Theta' and 'epsilon' must be positive numeric values"
    Err_cex <- "'Cex.point' and 'Cex.label' must be positive numeric values"

    if (!is.numeric(Phi) | !is.numeric(Theta) | !is.numeric(epsilon)) {
        stop(Err_a)
    }## if (!is.numeric(Phi) | !is.numeric(Theta) | !is.numeric(epsilon))

    if (min(c(Phi, Theta, epsilon))<0) {
        stop(Err_a)
    }## if (min(c(Phi, Theta, epsilon))<0)

    if (!is.numeric(Cex.point) | !is.numeric(Cex.label)) {
        stop(Err_cex)
    }## if (!is.numeric(Phi) | !is.numeric(Theta) | !is.numeric(epsilon))

    if (min(c(Cex.point, Cex.label))<0) {
        stop(Err_cex)
    }## if (min(c(Phi, Theta, epsilon))<0)

    ##-----------------------------------------------------------------------##
    ## Check path.result
    if (!is.null(path.result)) {
        if (!is.character(path.result)) {
            stop("'path.result' must be NULL or a character.")
        }## if (!is.character(path.result))
    }## if (!is.null(path.result))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrPCAHCPCgraphics()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##--------------------PCAplot3D----------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

PCAplot3D <- function(PCAcoord3D, matPCAeig,
                      colorINDfactor, clabPCA,
                      Phi=25, Theta=140, epsilon=0.2,
                      Cex.point=0.7, Cex.label=0.7) {
    ##-----------------------------------------------------------------------##
    ## ## clean up device : important to avoid problems when plotting all plots
    ## ## graphics::plot.new() ## put in main function
    ## cur_dev <- grDevices::dev.cur() ## store current device
    ## grDevices::pdf(file=NULL, width=11, height=8) ## open null device
    ## grDevices::dev.control("enable") ## turn on recording for the null device
    ## null_dev <- grDevices::dev.cur() ## store null device
    ## ## make sure we always clean up properly, even if smthg causes an error
    ## on.exit({
    ##     grDevices::dev.off(null_dev)
    ##     if (cur_dev > 1) grDevices::dev.set(cur_dev)
    ##     ## previous line: only set cur device if not null device
    ## })
    ##-----------------------------------------------------------------------##
    colorClust <- levels(colorINDfactor)
    colInd <- as.character(colorINDfactor)
    Nclust <- length(colorClust)

    matPCAeig <- round(matPCAeig, digits=2)

    ##-----------------------------------------------------------------------##
    plot3D::scatter3D(PCAcoord3D[, 1], PCAcoord3D[, 3], PCAcoord3D[, 2],
                      epsilon=epsilon, col=colInd,
                      clab=clabPCA, cex=Cex.point, theta=Theta, phi=Phi, d=2,
                      pch=20, colvar=NULL, bty="b2", ticktype="detailed",
                      xlab=paste0("dim1 (", matPCAeig[1, 2], "%)"),
                      ylab=paste0("dim3 (", matPCAeig[3, 2], "%)"),
                      zlab=paste0("dim2 (", matPCAeig[2, 2], "%)"))

    plot3D::text3D(PCAcoord3D[, 1] + epsilon,
                   PCAcoord3D[, 3] + epsilon,
                   PCAcoord3D[, 2] + epsilon,
                   labels=rownames(PCAcoord3D), add=TRUE, colkey=FALSE,
                   col=colInd, cex=Cex.label, font=2)

    ##-----------------------------------------------------------------------##
    ## ggPCA3d <- grDevices::recordPlot()
    ## return(ggPCA3d)
}## PCAplot3D()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##--------------------PCAplot2Dt---------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

PCAplot2Dt <- function(PCAcoord2D, matPCAeig, coord.t, data.name,
                       Vector.group, Vector.time, Vector.patient,
                       Mean.Accross.Time, colorINDfactor,
                       epsilon=0.2, Cex.point=0.7, Cex.label=0.7) {
    ##-----------------------------------------------------------------------##
    ## Plot preprocessing
    colorClust <- levels(colorINDfactor)
    colInd <- as.character(colorINDfactor)
    Nclust <- length(colorClust)

    LvlsPAT <- levels(factor(Vector.patient))
    NbPerCond <- as.numeric(table(Vector.patient))
    colorPoints <- as.character(colorINDfactor)

    if (!is.null(Vector.group)) {
        colorLinks <- as.character(colorINDfactor)
        vectGroup <- as.character(Vector.group)
    } else {
        colorLinks <- rep("grey", times=length(as.character(colorINDfactor)))
        vectGroup <- rep("G1", times=length(as.character(Vector.patient)))
    }## if(!is.null(Vector.group)) colorLinks <- "#999999"

    ##-----------------------------------------------------------------------##
    framePCA2t <- data.frame(PCAcoord2D[,c(1, 2)],
                             vGroup=vectGroup,
                             vPAT=Vector.patient,
                             vTime=Vector.time,
                             vColor=colorPoints)

    SMPLorder <- order(framePCA2t$vGroup, framePCA2t$vPAT, framePCA2t$vTime)
    colorLinks <- colorLinks[SMPLorder]
    framePCA2t <- framePCA2t[SMPLorder,]

    if (isFALSE(Mean.Accross.Time) & max(NbPerCond) > 1) {
        framePCA2t_2D <- framePCA2t
        framePCA2t_2D$vSizeP <-rep(Cex.point*3.5, times=nrow(framePCA2t_2D))
        framePCA2t_2D$vSizeL <-rep(Cex.label*3.5, times=nrow(framePCA2t_2D))
        framePCA2t_2D$vShape <-rep(16, times=nrow(framePCA2t_2D)) ## 16 <--> 19
    } else {
        name2D <- apply(unique(data.name[SMPLorder,]), 1, paste, collapse="_")
        name2D <- sort(as.character(name2D))

        Mean.t <- stats::aggregate(framePCA2t[,c(1, 2)],
                                   list(framePCA2t[,6],
                                        framePCA2t[,5]),
                                   mean)
        colnames(Mean.t)[c(1, 2)] <- c("vColor", "vTime")

        if (!is.null(Vector.group)) {
            Mean.t <- Mean.t[order(Mean.t$vColor, Mean.t$vTime),]
        } else {
            Mean.t <- Mean.t[order(Mean.t$vTime),]
        }## if (!is.null(Vector.group))

        Mean.t <- data.frame(Mean.t[, c(3, 4)],
                             vGroup=Mean.t$vColor,
                             vPAT=name2D,
                             Mean.t[,c(2, 1)])
        row.names(Mean.t) <- name2D

        if (!is.null(Vector.group)) {
            colorLinks <- as.character(Mean.t$vColor)
            MEANgroup <- 6
        } else {
            colorLinks <- rep("grey", times=nrow(Mean.t))
            Mean.t$vGroup <- rep("G1", times=nrow(Mean.t))
            MEANgroup <- 3
        }## if(!is.null(Vector.group)) colorLinks <- "#999999"

        framePCA2t_2D <- rbind(framePCA2t, Mean.t)
        framePCA2t_2D$vSizeP <-c(rep(Cex.point*2.5, times=nrow(framePCA2t)),
                                 rep(Cex.point*4, times=nrow(Mean.t)))
        framePCA2t_2D$vSizeL <-c(rep(Cex.label*2.5, times=nrow(framePCA2t)),
                                 rep(Cex.label*4, times=nrow(Mean.t)))
        framePCA2t_2D$vShape <-c(rep(16, times=nrow(framePCA2t)),
                                 rep(18, times=nrow(Mean.t))) ## 16 <--> 19
    }## if (isFALSE(Mean.Accross.Time) & max(NbPerCond) > 1)

    matPCAeig <- round(matPCAeig, digits=2)

    ##-----------------------------------------------------------------------##
    ggPCAt_2D <- ggplot2::ggplot(framePCA2t_2D,
                                 ggplot2::aes(x=framePCA2t_2D[, 1],
                                              y=framePCA2t_2D[, 2],
                                              group=framePCA2t_2D[, 4],
                                              label=row.names(framePCA2t_2D))) +
        ggplot2::theme_light() + ## ggplot2::theme_linedraw() +
        ggplot2::geom_hline(yintercept=0, linetype="dashed", linewidth=0.2) +
        ggplot2::geom_vline(xintercept=0, linetype="dashed", linewidth=0.2) +
        ggplot2::xlab(paste0("dim1 (", matPCAeig[1, 2], "%)")) +
        ggplot2::ylab(paste0("dim1 (", matPCAeig[2, 2], "%)")) +
        ggplot2::geom_point(colour=framePCA2t_2D[, 6],
                            size=framePCA2t_2D[, 7],
                            shape=framePCA2t_2D[, 9])

    if (isFALSE(Mean.Accross.Time) & max(NbPerCond) > 1) {
        ggPCAt_2D <- ggPCAt_2D +
            ggplot2::geom_path(colour=colorLinks) + ##framePCA2t_2D[, 5]
            ggrepel::geom_text_repel(data=framePCA2t_2D[, c(1, 2)],
                                     col=framePCA2t_2D[, 6],
                                     size=framePCA2t_2D[, 8],
                                     box.padding=0.5,
                                     max.overlaps=nrow(framePCA2t_2D)+2,
                                     segment.linetype=3,
                                     segment.colour=framePCA2t_2D[, 6],
                                     segment.size=0.5,
                                     min.segment.length=0)## "#7570b3",
    } else {
        ggPCAt_2D <- ggPCAt_2D +
            ggrepel::geom_text_repel(data=framePCA2t_2D[, c(1, 2)],
                                     col=framePCA2t_2D[, 6],
                                     size=framePCA2t_2D[, 8],
                                     box.padding=0.5,
                                     max.overlaps=nrow(framePCA2t_2D)+2,
                                     segment.linetype=3,
                                     segment.colour=framePCA2t_2D[, 6],
                                     segment.size=0.5,
                                     min.segment.length=0) +
            ggplot2::geom_path(data=Mean.t, inherit.aes=FALSE,
                               mapping=ggplot2::aes(x=Mean.t[, 1],
                                                    y=Mean.t[, 2],
                                                    group=Mean.t[, MEANgroup]),
                               colour=colorLinks)
    }## if (isFALSE(Mean.Accross.Time) & max(NbPerCond) > 1)

    ##-----------------------------------------------------------------------##
    return(ggPCAt_2D)
}## PCAplot2Dt()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##--------------------PCAplot3Dt---------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

PCAplot3Dt <- function(PCAcoord3D, matPCAeig, coord.t, data.name,
                       Vector.group, Vector.time, Vector.patient,
                       Mean.Accross.Time, colorINDfactor,
                       epsilon=0.2,
                       Cex.point=0.7, Cex.label=0.7,
                       Phi=25, Theta=140) {
    ##-----------------------------------------------------------------------##
    ## Plot preprocessing
    colorClust <- levels(colorINDfactor)
    colInd <- as.character(colorINDfactor)
    Nclust <- length(colorClust)

    LvlsPAT <- levels(factor(Vector.patient))
    NbPerCond <- as.numeric(table(Vector.patient))

    LvlsTime <- levels(factor(Vector.time))

    if (!is.null(Vector.group)) {
        colorLinks <- colInd
    } else {
        colorLinks <- rep("grey", times=length(colInd))
    }## if(!is.null(Vector.group))

    matPCAeig <- round(matPCAeig, digits=2)
    PCArgAlpha <- 0.05

    ##-----------------------------------------------------------------------##
    ## ## clean up device : important to avoid problems when plotting all plots
    ## ## graphics::plot.new() ## put in main function
    ## cur_dev <- grDevices::dev.cur() ## store current device
    ## grDevices::pdf(file=NULL, width=11, height=8) ## open null device
    ## grDevices::dev.control("enable") ## turn on recording for the null device
    ## null_dev <- grDevices::dev.cur() ## store null device
    ## ## make sure we always clean up properly, even if smthg causes an error
    ## on.exit({
    ##     grDevices::dev.off(null_dev)
    ##     if (cur_dev > 1) grDevices::dev.set(cur_dev)
    ##     ## previous line: only set cur device if not null device
    ## })
    ##-----------------------------------------------------------------------##
    if (isFALSE(Mean.Accross.Time) & max(NbPerCond) > 1) {
        ##-------------------------------------------------------------------##
        Smpl.sel <- LvlsPAT[1]
        IndexSmplSel <- which(Vector.patient == Smpl.sel)
        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

        if (length(coord.smpl) == 1) {
            coord.smplf <- c(coord.smpl, NA)
        } else {
            coord.smplf <- coord.smpl
        }## if(length(coord.smpl)==1)

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        plot3D::scatter3D(PCAcoord3D[coord.smplf, 1],
                          PCAcoord3D[coord.smplf, 3],
                          PCAcoord3D[coord.smplf, 2],
                          theta=Theta ,phi=Phi, colvar=NULL, cex=Cex.point,
                          type="b", pch=20, bty="b2", ticktype="detailed",
                          col=colorLinks[coord.smpl],
                          xlim=PCArange(PCAcoord3D[, 1], PCArgAlpha),
                          ylim=PCArange(PCAcoord3D[, 3], PCArgAlpha),
                          zlim=PCArange(PCAcoord3D[, 2], PCArgAlpha),
                          xlab=paste0("dim1 (", matPCAeig[1, 2], "%)"),
                          ylab=paste0("dim3 (", matPCAeig[3, 2], "%)"),
                          zlab=paste0("dim2 (", matPCAeig[2 ,2], "%)"))

        plot3D::text3D(PCAcoord3D[coord.smpl, 1],
                       PCAcoord3D[coord.smpl, 3],
                       PCAcoord3D[coord.smpl, 2],
                       labels=row.names(PCAcoord3D)[coord.smpl],
                       add=TRUE, colkey=FALSE, pch=20, cex=Cex.label,
                       col=colInd[coord.smpl])

        ## graphics::legend("right", title=LegendTitle, legend=data.color[,1],
        ##                  pch=20, horiz=FALSE, xpd=TRUE, # inset=c(-0.015),
        ##                  cex=0.8, col=data.color[,2 ]) ##cex=Cex.point*0.9

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        for (smpl in seq(from=2, to=length(LvlsPAT), by=1)) {
            ##---------------------------------------------------------------##
            Smpl.sel <- LvlsPAT[smpl]
            IndexSmplSel <- which(Vector.patient == Smpl.sel)
            times.smpl.sel <- factor(Vector.time[IndexSmplSel])
            coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

            if (length(coord.smpl) == 1) {
                coord.smplf <- c(coord.smpl, NA)
            } else {
                coord.smplf <- coord.smpl
            }## if (length(coord.smpl) == 1)

            ##---------------------------------------------------------------##
            ##---------------------------------------------------------------##
            plot3D::scatter3D(PCAcoord3D[coord.smplf, 1],
                              PCAcoord3D[coord.smplf, 3],
                              PCAcoord3D[coord.smplf, 2],
                              type="b", colvar=NULL, add=TRUE, cex=Cex.point,
                              bty="b2", pch=20, ticktype="detailed",
                              col=colorLinks[coord.smpl])

            plot3D::text3D(PCAcoord3D[coord.smpl, 1],
                           PCAcoord3D[coord.smpl, 3],
                           PCAcoord3D[coord.smpl, 2],
                           labels=row.names(PCAcoord3D)[coord.smpl],
                           add=TRUE, colkey=FALSE, cex=Cex.label, pch=20,
                           col=as.character(colorINDfactor)[coord.smpl])
        }## for(smpl in 2:length(LvlsPAT))

    } else {
        ##-------------------------------------------------------------------##
        Smpl.sel <- LvlsPAT[1]
        IndexSmplSel <- which(Vector.patient == Smpl.sel)
        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

        if (length(coord.smpl) == 1) {
            coord.smplf <- c(coord.smpl, NA)
        } else {
            coord.smplf <- coord.smpl
        }## if (length(coord.smpl) == 1)

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        plot3D::scatter3D(PCAcoord3D[coord.smplf, 1],
                          PCAcoord3D[coord.smplf, 3],
                          PCAcoord3D[coord.smplf, 2],
                          theta=Theta, phi=Phi, cex=Cex.point,
                          col=colorLinks[coord.smpl],
                          colvar=NULL, pch=20, bty="b2", type="p",
                          xlim=PCArange(PCAcoord3D[, 1], PCArgAlpha),
                          ylim=PCArange(PCAcoord3D[, 3], PCArgAlpha),
                          zlim=PCArange(PCAcoord3D[, 2], PCArgAlpha),
                          xlab=paste0("dim1 (", matPCAeig[1, 2], "%)"),
                          ylab=paste0("dim3 (", matPCAeig[3, 2], "%)"),
                          zlab=paste0("dim2 (", matPCAeig[2 ,2], "%)"))

        plot3D::text3D(PCAcoord3D[coord.smpl, 1],
                       PCAcoord3D[coord.smpl, 3],
                       PCAcoord3D[coord.smpl, 2],
                       labels=row.names(PCAcoord3D)[coord.smpl],
                       add=TRUE, colkey=FALSE, cex=Cex.label*0.67, pch=20,##0.6,
                       col=colInd[coord.smpl])

        ## graphics::legend("right", title=LegendTitle, legend=data.color[, 1],
        ##                  pch=20, horiz=FALSE, xpd=TRUE, cex=0.8,
        ##                  col=data.color[, 2]) ## inset=c(-0.015),

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        for (smpl in seq(from=2, to=length(LvlsPAT), by=1)) {
            Smpl.sel <- LvlsPAT[smpl]
            IndexSmplSel <- which(Vector.patient == Smpl.sel)
            times.smpl.sel <- factor(Vector.time[IndexSmplSel])
            coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

            if (length(coord.smpl) == 1) {
                coord.smplf <- c(coord.smpl, NA)
            } else {
                coord.smplf <- coord.smpl
            }## if(length(coord.smpl)==1)

            ##---------------------------------------------------------------##
            ##---------------------------------------------------------------##
            plot3D::scatter3D(PCAcoord3D[coord.smplf, 1],
                              PCAcoord3D[coord.smplf, 3],
                              PCAcoord3D[coord.smplf, 2],
                              add=TRUE, colvar=NULL, type="p", bty="b2",
                              cex=Cex.point, col=colorLinks[coord.smpl], pch=20)
            plot3D::text3D(PCAcoord3D[coord.smpl, 1],
                           PCAcoord3D[coord.smpl, 3],
                           PCAcoord3D[coord.smpl, 2],
                           labels=row.names(PCAcoord3D)[coord.smpl],
                           add=TRUE, colkey=FALSE, cex=Cex.label*0.67, pch=20,
                           col=colInd[coord.smpl])
        }## for(smpl in 2:length(LvlsPAT))

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        MEANseq <- seq_len(length(levels(factor(colorLinks))))

        for (col in MEANseq) {
            levelsLINKcol <- levels(factor(colorLinks))[col]
            Id.col <- which(colorLinks == levelsLINKcol)

            if (length(MEANseq) == 1) {
                MEANcolor <- colorClust ## myPaletteT(length(LvlsTime))
            } else {
                MEANcolor <- levelsLINKcol
            }## if (length(MEANseq) == 1)

            Mean.t <- stats::aggregate(PCAcoord3D[Id.col, c(1, 2, 3)],
                                       list(coord.t[Id.col]),
                                       mean)
            New.Name.3D <- apply(unique(data.name[Id.col,]), 1,
                                 paste,
                                 collapse="_")
            New.Name.3D <- as.character(New.Name.3D)

            plot3D::scatter3D(Mean.t[,2], Mean.t[,4], Mean.t[,3],
                              add=TRUE, colvar=NULL, cex=Cex.point, type="b",
                              pch=20, bty="b2", col=levelsLINKcol)

            plot3D::text3D(Mean.t[,2], Mean.t[,4], Mean.t[,3],
                           labels=New.Name.3D, cex=Cex.label,
                           col=MEANcolor, add=TRUE, colkey=FALSE, pch=20)
        }## for(col in 1:length(levels(factor(colorLinks))))
    }## if(Mean.Accross.Time==FALSE & max(NbPerCond)>1)

    ##-----------------------------------------------------------------------##
    ## ggPCA3Dt <- grDevices::recordPlot()
    ## return(ggPCA3Dt)
}## PCAplot3Dt()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##--------------------PCArange-----------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

PCArange <- function(v, alpha) {
    rg <- c(min(v) - alpha*(max(v) - min(v)), max(v) + alpha*(max(v) - min(v)))
    return(rg)
}## PCArange()



