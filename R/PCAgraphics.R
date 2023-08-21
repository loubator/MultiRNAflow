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
#' @param Supp.del.sample \code{TRUE} or \code{FALSE}.
#' If \code{FALSE}, the samples selected with \code{sample.deletion} will
#' be deleted.
#' If \code{TRUE}, the samples selected with \code{sample.deletion} will
#' be plotted.
#' These individuals are called supplementary individuals in
#' [FactoMineR::PCA()].
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
#' @param D3.mouvement \code{TRUE} or \code{FALSE}.
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
#' on the experimental design (if \code{Plot.PCA=TRUE})
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
#' The interactive 3D graphs will be plotted only if \code{D3.mouvement=TRUE}.
#'
#' @seealso This function is called by our function
#' [PCAanalysis()]
#' and calls our function
#' [PCArealization()].
#'
#' @importFrom S4Vectors metadata
#' @importFrom scales hue_pal
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats aggregate
#' @importFrom FactoMineR plot.PCA
#' @importFrom plot3D scatter3D text3D
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
#' ##-------------------------------------------------------------------------#
#' resPCAgraph <- PCAgraphics(SEresNorm=resNorm,
#'                            DATAnorm=TRUE,
#'                            gene.deletion=c("Gene1", "Gene5"),
#'                            sample.deletion=c(2,6),
#'                            Supp.del.sample=FALSE,
#'                            Plot.PCA=TRUE,
#'                            Mean.Accross.Time=FALSE,
#'                            Color.Group=GROUPcolor,
#'                            D3.mouvement=FALSE,
#'                            Phi=25, Theta=140, Cex.label=0.7,
#'                            Cex.point=0.7, epsilon=0.2,
#'                            path.result=NULL, Name.file.pca=NULL)

PCAgraphics <- function(SEresNorm,
                        DATAnorm=TRUE,
                        gene.deletion=NULL,
                        sample.deletion=NULL,
                        Supp.del.sample=FALSE,
                        Plot.PCA=TRUE,
                        Mean.Accross.Time=FALSE,
                        Color.Group=NULL,
                        D3.mouvement=FALSE,
                        Phi=25, Theta=140, epsilon=0.2,
                        Cex.point=0.7, Cex.label=0.7,
                        path.result=NULL,
                        Name.file.pca=NULL) {
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check TRUE, FALSE
    if (!isTRUE(Supp.del.sample) & !isFALSE(Supp.del.sample)) {
        stop("'Supp.del.sample' must be TRUE or FALSE.")
    }## if (!isTRUE(Supp.del.sample) & !isFALSE(Supp.del.sample))

    if (!isTRUE(Plot.PCA) & !isFALSE(Plot.PCA)) {
        stop("'Plot.PCA' must be TRUE or FALSE.")
    }## if (!isTRUE(Plot.PCA) & !isFALSE(Plot.PCA))

    if (!isTRUE(Mean.Accross.Time) & !isFALSE(Mean.Accross.Time)) {
        stop("'Mean.Accross.Time' must be TRUE or FALSE.")
    }## if (!isTRUE(Mean.Accross.Time) & !isFALSE(Mean.Accross.Time))

    if (!isTRUE(D3.mouvement) & !isFALSE(D3.mouvement)) {
        stop("'D3.mouvement' must be TRUE or FALSE.")
    }## if (!isTRUE(D3.mouvement) & !isFALSE(D3.mouvement))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check numeric
    if (!is.numeric(Phi) | !is.numeric(Theta) | !is.numeric(epsilon)) {
        Err_a <- "'Phi', 'Theta' and 'epsilon' must be positive numeric values"
        stop(Err_a)
        if (min(c(Phi, Theta, epsilon))<0) {
            stop(Err_a)
        }## if (min(c(Phi, Theta, epsilon))<0)
    }## if (!is.numeric(Phi) | !is.numeric(Theta) | !is.numeric(epsilon))

    if (!is.numeric(Cex.point) | !is.numeric(Cex.label)) {
        stop("'Cex.point' and 'Cex.label' must be positive numeric values")
        if (min(c(Cex.point, Cex.label))<0) {
            stop("'Cex.point' and 'Cex.label' must be positive numeric values")
        }## if (min(c(Phi, Theta, epsilon))<0)
    }## if (!is.numeric(Phi) | !is.numeric(Theta) | !is.numeric(epsilon))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check deletion
    if (!is.null(gene.deletion)) {
        if (!is.character(gene.deletion)) {
            if (!is.numeric(gene.deletion)) {
                Err_delint <- paste("'gene.deletion' must be either NULL",
                                    "either character or",
                                    "non negative integers.")
                stop(Err_delint)
            } else {
                if (sum(abs(floor(gene.deletion) - gene.deletion)) !=0) {
                    stop(Err_delint)
                }## if (floor(gene.deletion) != gene.deletion)
                if (min(gene.deletion) <= 0) {
                    stop(Err_delint)
                }## if (min(gene.deletion) <= 0)
            }## if (!is.numeric(gene.deletion))
        }## if (!is.character(gene.deletion))
    }## if (!is.null(gene.deletion))

    ## Check deletion
    if (!is.null(sample.deletion)) {
        if (!is.character(sample.deletion)) {
            if (!is.numeric(sample.deletion)) {
                Err_delint <- paste("'sample.deletion' must be either NULL",
                                    "either character or",
                                    "non negative integers.")
                stop(Err_delint)
            } else{
                if (sum(abs(floor(sample.deletion) - sample.deletion)) !=0) {
                    stop(Err_delint)
                }## if (floor(sample.deletion) != sample.deletion)
                if (min(sample.deletion) <= 0) {
                    stop(Err_delint)
                }## if (min(sample.deletion) <= 0)
            }## if (!is.numeric(sample.deletion))
        }## if (!is.character(sample.deletion))
    }## if (!is.null(sample.deletion))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check deletion
    if (!is.null(path.result)) {
        if (!is.character(path.result)) {
            stop("'path.result' must be NULL or a character.")
        }## if (!is.character(path.result))
    }## if (!is.null(path.result))

    if (!is.null(Name.file.pca)) {
        if (!is.character(Name.file.pca)) {
            stop("'Name.file.pca' must be NULL or a character.")
        }## if (!is.character(Name.file.pca))
    }## if (!is.null(Name.file.pca))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## PCA
    SEresPCA <- PCArealization(SEresNorm=SEresNorm,
                               DATAnorm=DATAnorm,
                               sample.deletion=sample.deletion,
                               gene.deletion=gene.deletion,
                               Supp.del.sample=Supp.del.sample)

    listFCTRS <- S4Vectors::metadata(SEresPCA)$PCA$List.Factors
    res.PCA <- S4Vectors::metadata(SEresPCA)$PCA$res.pca

    Vector.time <- listFCTRS$Vector.time
    Vector.group <- listFCTRS$Vector.group
    Vector.patient <- listFCTRS$Vector.patient
    LvlsPAT <- levels(factor(Vector.patient))

    PCAqualiSup <- res.PCA$call$quali.sup$quali.sup
    coordPCAind <- res.PCA$ind$coord

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if (isTRUE(Plot.PCA) | !is.null(path.result)) {
        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## Color assignment
        if (!is.null(Vector.time) & is.null(Vector.group)) {
            Tlevels <- as.character(levels(factor(Vector.time)))
            Tlevels <- paste0("t", gsub("t", "", gsub("T", "", Tlevels)))
            levels(PCAqualiSup$Quali.Sup.Time)<-Tlevels

            ##----------------------------------------------------------------#
            Color.Time <- NULL
            if (is.null(Color.Time)) {
                TlevelsCOLOR <- c("#737373",## "#252525"
                                  scales::hue_pal()(length(Tlevels)-1))
                Color.Time <- data.frame(Name=Tlevels, Col=TlevelsCOLOR)

                data.color <- rbind(Color.Time,
                                    rep(NA, times=ncol(Color.Time)))
                data.color <- data.color[-which(is.na(data.color$Col)),]
            } else {
                Id.LevelColT <- order(Color.Time[, 1])
                Color.Time <- data.frame(Name=Tlevels,
                                         Col=Color.Time[Id.LevelColT, 2])

                data.color <- rbind(Color.Time, rep(NA,times=ncol(Color.Time)))
                data.color <- data.color[-which(is.na(data.color$Col)),]
            }## if(is.null(Color.Time))

            ind.color <- PCAqualiSup$Quali.Sup.Time
            levels(ind.color) <- data.color[, 2]
            LegendTitle <- "Time"
            ##----------------------------------------------------------------#
        } else {
            Glevels <- levels(factor(Vector.group))

            if (is.null(Color.Group)) {
                MypaletteG <- c(RColorBrewer::brewer.pal(8, "Dark2"),
                                RColorBrewer::brewer.pal(8, "Set2"))

                if (length(Glevels) > 16) {
                    MypaG2 <- scales::hue_pal(l=90)(seq_len(length(Glevels)-1))
                    MypaletteG <- c(MypaletteG, MypaG2)
                }## if(length(Glevels)>16)

                MypaletteG <- MypaletteG[seq_len(length(Glevels))]
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

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        if (is.null(Name.file.pca) == TRUE) {
            Name.file.pca.f <- "Graph"
        } else {
            Name.file.pca.f <- Name.file.pca
        }## if(is.null(Name.file.pca)==TRUE)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## PCA plot from plot.PCA() ## Qualitative factor are in magenta
        options(ggrepel.max.overlaps=30)

        g.2DPCA <- FactoMineR::plot.PCA(res.PCA, axes=c(1,2),
                                        choix="ind", habillage="ind",
                                        col.hab=as.character(ind.color),
                                        col.quali="magenta",
                                        shadowtext=TRUE, autoLab="y",
                                        graph.type="ggplot", cex=Cex.label)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        if (!is.null(path.result)) {
            TitlePCA2D <- paste0(Name.file.pca.f, "_PCA2D.pdf")

            grDevices::pdf(file=file.path(path.result, TitlePCA2D),
                           width=11, height=8)
            print(g.2DPCA)
            grDevices::dev.off()

            if (isTRUE(Plot.PCA)) {
                print(g.2DPCA)
            }## if(Plot.PCA==TRUE)

        } else {
            if (isTRUE(Plot.PCA)) {
                print(g.2DPCA)
            }## if(isTRUE(Plot.PCA))
        }## if(!is.null(path.result))

        options(ggrepel.max.overlaps=10)

        ##--------------------------------------------------------------------#
        ## 3D PCA colored by biological condition or time points
        ## par(mar=c(4.1, 4.1, 2.1, 2.1))#, xpd=TRUE)
        data.3D <- coordPCAind[, c(1, 2, 3)]

        if (!is.null(path.result)) {
            TitlePCA3D <- paste0(Name.file.pca.f, "_PCA3D.pdf")

            grDevices::pdf(file=file.path(path.result, TitlePCA3D),
                           width=11, height=8)

            plot3D::scatter3D(data.3D[, 1], data.3D[, 3], data.3D[, 2],
                              epsilon=epsilon, pch=20, colvar=NULL,
                              col=as.character(ind.color), bty="b2",
                              cex=Cex.point, clab=c("", data.color$Name),
                              ticktype="detailed", theta=Theta, phi=Phi, d=2,
                              xlab=paste0("dim1 (",
                                         round(res.PCA$eig[, 2][1], digits=2),
                                         "%)"),
                              ylab=paste0("dim3 (",
                                          round(res.PCA$eig[, 2][3], digits=2),
                                         "%)"),
                              zlab=paste0("dim2 (",
                                         round(res.PCA$eig[, 2][2], digits=2),
                                         "%)"))

            plot3D::text3D(data.3D[, 1] + epsilon,
                           data.3D[, 3] + epsilon,
                           data.3D[, 2] + epsilon,
                           labels=rownames(data.3D), add=TRUE, colkey=FALSE,
                           col=as.character(ind.color), cex=Cex.label, font=2)

            grDevices::dev.off()


            if (isTRUE(Plot.PCA)) {
                plot3D::scatter3D(data.3D[, 1], data.3D[, 3], data.3D[, 2],
                                  pch=20, colvar=NULL, cex=Cex.point,
                                  col=as.character(ind.color), bty="b2", d=2,
                                  ticktype="detailed", theta=Theta, phi=Phi,
                                  clab=c("", data.color$Name), epsilon=epsilon,
                                  xlab=paste0("dim1 (",
                                              round(res.PCA$eig[,2][1],
                                                    digits=2), "%)"),
                                  ylab=paste0("dim3 (",
                                              round(res.PCA$eig[,2][3],
                                                    digits=2), "%)"),
                                  zlab=paste0("dim2 (",
                                              round(res.PCA$eig[,2][2],
                                                    digits=2), "%)"))

                plot3D::text3D(data.3D[, 1] + epsilon,
                               data.3D[, 3] + epsilon,
                               data.3D[, 2] + epsilon,
                               labels=rownames(data.3D), add=TRUE,
                               col=as.character(ind.color), colkey=FALSE,
                               cex=Cex.label, font=2)

                ##------------------------------------------------------------#
                ## 3D PCA in rgl windows
                if (D3.mouvement==TRUE) {
                    plot3Drgl::plotrgl()
                }## if(D3.mouvement==TRUE)
            }## Plot.PCA==TRUE

        } else {
            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            if (isTRUE(Plot.PCA)) {
                plot3D::scatter3D(data.3D[, 1], data.3D[, 3], data.3D[, 2],
                                  pch=20, colvar=NULL, theta=Theta, phi=Phi,
                                  col=as.character(ind.color), epsilon=epsilon,
                                  cex=Cex.point, bty="b2", ticktype="detailed",
                                  d=2, clab=c("", data.color$Name),
                                  xlab=paste0("dim1 (",
                                              round(res.PCA$eig[, 2][1],
                                                    digits=2), "%)"),
                                  ylab=paste0("dim3 (",
                                              round(res.PCA$eig[, 2][3],
                                                    digits=2), "%)"),
                                  zlab=paste0("dim2 (",
                                              round(res.PCA$eig[, 2][2],
                                                    digits=2), "%)"))

                plot3D::text3D(data.3D[, 1] + epsilon,
                               data.3D[, 3] + epsilon,
                               data.3D[, 2] + epsilon,
                               labels=rownames(data.3D), add=TRUE,
                               col=as.character(ind.color), cex=Cex.label,
                               colkey=FALSE, font=2)

                ##------------------------------------------------------------#
                ## 3D PCA in rgl windows
                if (isTRUE(D3.mouvement)) {
                    plot3Drgl::plotrgl()
                }## if(D3.mouvement==TRUE)

            }## Plot.PCA==TRUE
        }## if(is.null(path.result)==FALSE)
        ##--------------------------------------------------------------------#
        ## graphics::legend("right", title=LegendTitle,
        ##                  legend=data.color[,1],
        ##                  pch=20, horiz=FALSE, xpd=TRUE,
        ##                  cex=0.8, # inset=c(-0.015),
        ##                  col=data.color[,2])
        ##
        ##
        ##
        ## s3d.all <- grDevices::recordPlot()
        ## graphics::plot.new() ## clean up device
        ##--------------------------------------------------------------------#
        ### R plot without showing the graphic window
        ## dev.control('enable') # enable display list
        ## plot(rnorm(10))
        ## obj = recordPlot()
        ## dev.off()
        ## replayPlot(obj)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        if (!is.null(Vector.time)) {
            ##----------------------------------------------------------------#
            coord.t <- PCAqualiSup$Quali.Sup.Time
            n.row <- length(coord.t)
            index.order <- seq_len(n.row)
            nb.time <- length(unique(coord.t))

            ##----------------------------------------------------------------#
            if (!is.null(Vector.group)) {
                col.link <- as.character(ind.color)
            } else {
                col.link <- rep("grey", times=length(as.character(ind.color)))
            }## if(!is.null(Vector.group))

            ##----------------------------------------------------------------#
            NbPerCond <- as.numeric(table(Vector.patient))

            data.name <- data.frame(rep("Mean", times=nrow(PCAqualiSup)),
                                    PCAqualiSup)

            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            if (!is.null(path.result)) {
                TitlePCA2Dlk <- paste0(Name.file.pca.f, "_PCAlink2D.pdf")

                grDevices::pdf(file=file.path(path.result, TitlePCA2Dlk),
                               width=11, height=8)

                plot(NA, type="c", col="green", #main="2D PCA",
                     xlim=c(min(coordPCAind[, 1]) - epsilon,
                            max(coordPCAind[, 1]) + epsilon),
                     ylim=c(min(coordPCAind[, 2]) - epsilon,
                            max(coordPCAind[, 2]) + epsilon),
                     xlab=paste0("dim1 (",
                                 round(res.PCA$eig[, 2][1], digits=2), "%)"),
                     ylab=paste0("dim2 (",
                                round(res.PCA$eig[, 2][2], digits=2), "%)"))

                ##------------------------------------------------------------#
                ##------------------------------------------------------------#
                if (isFALSE(Mean.Accross.Time) & max(NbPerCond) > 1) {
                    ##--------------------------------------------------------#
                    for (smpl in seq_len(length(LvlsPAT))) {

                        Smpl.sel <- LvlsPAT[smpl]
                        IndexSmplSel <- which(Vector.patient == Smpl.sel)
                        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                        graphics::lines(coordPCAind[coord.smpl, 1],
                                        coordPCAind[coord.smpl, 2],
                                        pch=22, type="c", lty=1,
                                        col=col.link[coord.smpl])
                        graphics::text(coordPCAind[coord.smpl, 1] + epsilon,
                                       coordPCAind[coord.smpl, 2] + epsilon,
                                       labels=row.names(coordPCAind[coord.smpl,
                                                                    ]),
                                       cex=Cex.label,#0.6,
                                       col=as.character(ind.color)[coord.smpl])
                    }## for(smpl in 1:length(LvlsPAT))

                } else {

                    ##--------------------------------------------------------#
                    for (smpl in seq_len(length(LvlsPAT))) {

                        Smpl.sel <- LvlsPAT[smpl]
                        IndexSmplSel <- which(Vector.patient == Smpl.sel)
                        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                        graphics::points(coordPCAind[coord.smpl, 1],
                                         coordPCAind[coord.smpl, 2],
                                         pch=19,
                                         cex=Cex.point,
                                         col=col.link[coord.smpl])
                        graphics::text(coordPCAind[coord.smpl, 1] + epsilon,
                                       coordPCAind[coord.smpl, 2] + epsilon,
                                       labels=row.names(coordPCAind[coord.smpl,
                                                                    ]),
                                       cex=Cex.label,
                                       col=as.character(ind.color)[coord.smpl])
                    }## for(smpl in 1:length(LvlsPAT))

                    ##--------------------------------------------------------#
                    ##--------------------------------------------------------#
                    for (col in seq_len(length(levels(ind.color)))) {
                        Id.col <- which(as.character(ind.color)==levels(ind.color)[col])
                        Mean.t <- stats::aggregate(coordPCAind[Id.col, c(1,2)],
                                                   list(coord.t[Id.col]), mean)
                        New.Name.2D <- apply(unique(data.name[Id.col,]), 1,
                                             paste, collapse="_")

                        graphics::lines(Mean.t[, 2], Mean.t[, 3],
                                        type="c", pch=22, lwd=2,
                                        col=levels(ind.color)[col])
                        graphics::text(Mean.t[, 2], Mean.t[, 3],
                                       labels=New.Name.2D, cex=Cex.label*0.9,
                                       col=levels(ind.color)[col])
                    }## for(col in 1:length(levels(ind.color)))
                }## if(Mean.Accross.Time==FALSE & max(NbPerCond)>1)
                ##------------------------------------------------------------#
                grDevices::dev.off()
            }## if(is.null(path.result)==FALSE)
            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            if (isTRUE(Plot.PCA)) {

                plot(NA, type="c", col="green",
                     xlim=c(min(coordPCAind[, 1]) - epsilon,
                            max(coordPCAind[, 1]) + epsilon),
                     ylim=c(min(coordPCAind[, 2]) - epsilon,
                            max(coordPCAind[, 2]) + epsilon),
                     xlab=paste0("dim1 (",
                                round(res.PCA$eig[, 2][1], digits=2),"%)"),
                     ylab=paste0("dim2 (",
                                round(res.PCA$eig[, 2][2], digits=2),"%)"))
                ##
                ## xlegend<-max(coordPCAind[,1])
                ##         +epsilon+(max(coordPCAind[, 1])+epsilon)/8
                ##
                ## graphics::legend(x=xlegend, y=0,
                ##                  title=LegendTitle,
                ##                  legend = data.color[, 1],
                ##                  pch=20, horiz=FALSE, xpd = TRUE,
                ##                  cex=0.8, #inset=c(-0.03),
                ##                  col=data.color[, 2])

                ##------------------------------------------------------------#
                ##------------------------------------------------------------#
                if (isFALSE(Mean.Accross.Time) & max(NbPerCond) > 1) {
                    ##--------------------------------------------------------#
                    for (smpl in seq_len(length(LvlsPAT))) {

                        Smpl.sel <- LvlsPAT[smpl]
                        IndexSmplSel <- which(Vector.patient == Smpl.sel)
                        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                        graphics::lines(coordPCAind[coord.smpl, 1],
                                        coordPCAind[coord.smpl, 2],
                                        pch=22, type="c", lty=1,
                                        col=col.link[coord.smpl])
                        graphics::text(coordPCAind[coord.smpl, 1] + epsilon,
                                       coordPCAind[coord.smpl, 2] + epsilon,
                                       labels=row.names(coordPCAind[coord.smpl,
                                                                    ]),
                                       cex=Cex.label,
                                       col=as.character(ind.color)[coord.smpl])
                    }## for(smpl in 1:length(LvlsPAT))

                } else {
                    ##--------------------------------------------------------#
                    for (smpl in seq_len(length(LvlsPAT))) {
                        Smpl.sel <- LvlsPAT[smpl]
                        IndexSmplSel <- which(Vector.patient == Smpl.sel)
                        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                        graphics::points(coordPCAind[coord.smpl, 1],
                                         coordPCAind[coord.smpl, 2],
                                         pch=19, cex=Cex.point,
                                         col=col.link[coord.smpl])
                        graphics::text(coordPCAind[coord.smpl, 1] + epsilon,
                                       coordPCAind[coord.smpl, 2] + epsilon,
                                       labels=row.names(coordPCAind[coord.smpl,
                                                                    ]),
                                       cex=Cex.label,
                                       col=as.character(ind.color)[coord.smpl])
                    }## for(smpl in 1:length(LvlsPAT))

                    ##--------------------------------------------------------#
                    ##--------------------------------------------------------#
                    for (col in seq_len(length(levels(ind.color)))) {
                        Id.col <- which(as.character(ind.color)==levels(ind.color)[col])
                        Mean.t <- stats::aggregate(coordPCAind[Id.col, c(1,2)],
                                                   list(coord.t[Id.col]), mean)
                        New.Name.2D <- apply(unique(data.name[Id.col,]), 1,
                                             paste, collapse="_")

                        graphics::lines(Mean.t[, 2], Mean.t[, 3],
                                        type="c", pch=22, lwd=2,
                                        col=levels(ind.color)[col])
                        graphics::text(Mean.t[, 2], Mean.t[, 3],
                                       labels=New.Name.2D, cex=Cex.label*0.9,
                                       col=levels(ind.color)[col]) ##*0.6,
                    }## for(col in 1:length(levels(ind.color)))
                }## if(Mean.Accross.Time==FALSE & max(NbPerCond)>1)
            }## if(Plot.PCA==TRUE)
            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            ##
            ## g.2DPCA.t<-grDevices::recordPlot()
            ## graphics::plot.new() ## clean up device

            ##----------------------------------------------------------------#
            if (!is.null(path.result)) {
                TitlePCA3Dlk <- paste0(Name.file.pca.f, "_PCAlink3D.pdf")

                grDevices::pdf(file=file.path(path.result, TitlePCA3Dlk),
                               width=11, height=8)
                ##------------------------------------------------------------#
                if (isFALSE(Mean.Accross.Time) & max(NbPerCond) > 1) {
                    Smpl.sel <- LvlsPAT[1]
                    IndexSmplSel <- which(Vector.patient == Smpl.sel)
                    times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                    coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                    ##--------------------------------------------------------#
                    if (length(coord.smpl) == 1) {
                        coord.smplf <- c(coord.smpl, NA)
                    } else {
                        coord.smplf <- coord.smpl
                    }## if(length(coord.smpl)==1)

                    #---------------------------------------------------------#
                    ## s3d.2t<-
                    plot3D::scatter3D(coordPCAind[coord.smplf, 1],
                                      coordPCAind[coord.smplf, 3],
                                      coordPCAind[coord.smplf, 2],
                                      theta=Theta ,phi=Phi, colvar=NULL,
                                      type="b", pch=20, bty="b2",
                                      ticktype="detailed", cex=Cex.point,
                                      col=col.link[coord.smpl],
                                      xlim=c(min(coordPCAind[, 1]),
                                             max(coordPCAind[, 1])),
                                      ylim=c(min(coordPCAind[, 3]),
                                             max(coordPCAind[, 3])),
                                      zlim=c(min(coordPCAind[, 2]),
                                             max(coordPCAind[, 2])),
                                      xlab=paste0("dim1 (",
                                                 round(res.PCA$eig[,2][1],
                                                       digits=2),
                                                 "%)"),
                                      ylab=paste0("dim3 (",
                                                 round(res.PCA$eig[,2][3],
                                                       digits=2),
                                                 "%)"),
                                      zlab=paste0("dim2 (",
                                                 round(res.PCA$eig[,2][2],
                                                       digits=2),
                                                 "%)"))
                    #
                    plot3D::text3D(coordPCAind[coord.smpl, 1],
                                   coordPCAind[coord.smpl, 3],
                                   coordPCAind[coord.smpl, 2],
                                   labels=row.names(coordPCAind)[coord.smpl],
                                   add=TRUE,
                                   colkey=FALSE,
                                   cex=Cex.label, pch=20,
                                   col=as.character(ind.color)[coord.smpl])

                    ## graphics::legend("right", title=LegendTitle,
                    ##                  legend=data.color[,1],
                    ##                  pch=20, horiz=FALSE, xpd=TRUE,
                    ##                  cex=0.8, # inset=c(-0.015),
                    ##                  col=data.color[,2]) ##cex=Cex.point*0.9

                    ##--------------------------------------------------------#
                    for (smpl in seq(from=2, to=length(LvlsPAT), by=1)) {

                        Smpl.sel <- LvlsPAT[smpl]
                        IndexSmplSel <- which(Vector.patient == Smpl.sel)
                        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                        ##----------------------------------------------------#
                        if (length(coord.smpl) == 1) {
                            coord.smplf <- c(coord.smpl, NA)
                        } else {
                            coord.smplf <- coord.smpl
                        }

                        ##----------------------------------------------------#
                        plot3D::scatter3D(coordPCAind[coord.smplf, 1],
                                          coordPCAind[coord.smplf, 3],
                                          coordPCAind[coord.smplf, 2],
                                          type="b", colvar=NULL, add=TRUE,
                                          bty="b2", pch=20,
                                          cex=Cex.point, ticktype="detailed",
                                          col=col.link[coord.smpl])

                        plot3D::text3D(coordPCAind[coord.smpl, 1],
                                       coordPCAind[coord.smpl, 3],
                                       coordPCAind[coord.smpl, 2],
                                       labels=row.names(coordPCAind)[coord.smpl],
                                       add=TRUE, colkey=FALSE,
                                       cex=Cex.label, pch=20,
                                       col=as.character(ind.color)[coord.smpl])
                    }## for(smpl in 2:length(LvlsPAT))
                } else {
                    ##--------------------------------------------------------#
                    Smpl.sel <- LvlsPAT[1]
                    IndexSmplSel <- which(Vector.patient == Smpl.sel)
                    times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                    coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                    ##--------------------------------------------------------#
                    if (length(coord.smpl)==1) {
                        coord.smplf <- c(coord.smpl, NA)
                    } else {
                        coord.smplf <- coord.smpl
                    }## if(length(coord.smpl)==1)

                    ##--------------------------------------------------------#
                    ## s3d.2t<-
                    plot3D::scatter3D(coordPCAind[coord.smplf, 1],
                                      coordPCAind[coord.smplf, 3],
                                      coordPCAind[coord.smplf, 2],
                                      theta=Theta ,phi=Phi, pch=20, bty="b2",
                                      type="p", cex=Cex.point, colvar=NULL,
                                      col=col.link[coord.smpl],
                                      xlim=c(min(coordPCAind[, 1]),
                                             max(coordPCAind[, 1])),
                                      ylim=c(min(coordPCAind[, 3]),
                                             max(coordPCAind[, 3])),
                                      zlim=c(min(coordPCAind[, 2]),
                                             max(coordPCAind[, 2])),
                                      xlab=paste0("dim1 (",
                                                  round(res.PCA$eig[, 2][1],
                                                        digits=2),
                                                  "%)"),
                                      ylab=paste0("dim3 (",
                                                  round(res.PCA$eig[, 2][3],
                                                        digits=2),
                                                  "%)"),
                                      zlab=paste0("dim2 (",
                                                  round(res.PCA$eig[, 2][2],
                                                        digits=2),
                                                  "%)"))
                    ##
                    plot3D::text3D(coordPCAind[coord.smpl, 1],
                                   coordPCAind[coord.smpl, 3],
                                   coordPCAind[coord.smpl, 2],
                                   labels=row.names(coordPCAind)[coord.smpl],
                                   add=TRUE, colkey=FALSE,
                                   cex=Cex.label, pch=20, ##0.6,
                                   col=as.character(ind.color)[coord.smpl])

                    ## graphics::legend("right", title=LegendTitle,
                    ##                  legend = data.color[,1],
                    ##                  pch = 20, horiz=FALSE, xpd = TRUE,
                    ##                  cex=0.8, # inset=c(-0.015),
                    ##                  col = data.color[,2])

                    ##--------------------------------------------------------#
                    for (smpl in seq(from=2, to=length(LvlsPAT), by=1)) {
                        Smpl.sel <- LvlsPAT[smpl]
                        IndexSmplSel <- which(Vector.patient == Smpl.sel)
                        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                        ##----------------------------------------------------#
                        if (length(coord.smpl) == 1) {
                            coord.smplf<-c(coord.smpl, NA)
                        }else{
                            coord.smplf<-coord.smpl
                        }## if(length(coord.smpl)==1)

                        ##----------------------------------------------------#
                        plot3D::scatter3D(coordPCAind[coord.smplf, 1],
                                          coordPCAind[coord.smplf, 3],
                                          coordPCAind[coord.smplf, 2],
                                          type="p", colvar=NULL, bty="b2",
                                          add=TRUE, pch=20, cex=Cex.point,
                                          col=col.link[coord.smpl])
                        plot3D::text3D(coordPCAind[coord.smpl, 1],
                                       coordPCAind[coord.smpl, 3],
                                       coordPCAind[coord.smpl, 2],
                                       labels=row.names(coordPCAind)[coord.smpl],
                                       add=TRUE, colkey=FALSE,
                                       cex=Cex.label, pch=20,
                                       col=as.character(ind.color)[coord.smpl])
                    }## for(smpl in 2:length(LvlsPAT))

                    ##--------------------------------------------------------#
                    data.name <- data.frame(rep("Mean",
                                                times=nrow(PCAqualiSup)),
                                            PCAqualiSup)

                    ##--------------------------------------------------------#
                    for(col in seq_len(length(levels(factor(col.link))))){
                        Id.col <- which(col.link==levels(factor(col.link))[col])
                        Mean.t <- stats::aggregate(coordPCAind[Id.col,
                                                               c(1, 2, 3)],
                                                   list(coord.t[Id.col]), mean)
                        New.Name.3D <- apply(unique(data.name[Id.col,]), 1,
                                             paste,
                                             collapse="_")

                        plot3D::scatter3D(Mean.t[,2], Mean.t[,4], Mean.t[,3],
                                          type="b", colvar=NULL, cex=Cex.point,
                                          add=TRUE, pch=20, bty="b2",
                                          col=levels(factor(col.link))[col])

                        plot3D::text3D(Mean.t[,2], Mean.t[,4], Mean.t[,3],
                                       labels=New.Name.3D, cex=Cex.label,
                                       col=levels(factor(col.link))[col],
                                       add=TRUE, colkey=FALSE, pch=20)
                    }## for(col in 1:length(levels(factor(col.link))))
                }## if(Mean.Accross.Time==FALSE & max(NbPerCond)>1)
                ##------------------------------------------------------------#
                grDevices::dev.off()
            }## if(is.null(path.result)==FALSE)
            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            if (isTRUE(Plot.PCA)) {
                if (isFALSE(Mean.Accross.Time) & max(NbPerCond) > 1) {
                    Smpl.sel <- LvlsPAT[1]
                    IndexSmplSel <- which(Vector.patient == Smpl.sel)
                    times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                    coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                    ##--------------------------------------------------------#
                    if (length(coord.smpl) == 1) {
                        coord.smplf <- c(coord.smpl, NA)
                    } else {
                        coord.smplf <- coord.smpl
                    }## if(length(coord.smpl)==1)

                    ##--------------------------------------------------------#
                    ## s3d.2t<-
                    plot3D::scatter3D(coordPCAind[coord.smplf, 1],
                                      coordPCAind[coord.smplf, 3],
                                      coordPCAind[coord.smplf, 2],
                                      theta=Theta ,phi=Phi, colvar=NULL,
                                      type="b", pch=20, bty="b2",
                                      ticktype="detailed", cex=Cex.point,#0.6,
                                      col=col.link[coord.smpl],
                                      xlim=c(min(coordPCAind[, 1]),
                                             max(coordPCAind[, 1])),
                                      ylim=c(min(coordPCAind[, 3]),
                                             max(coordPCAind[, 3])),
                                      zlim=c(min(coordPCAind[, 2]),
                                             max(coordPCAind[, 2])),
                                      xlab=paste0("dim1 (",
                                                 round(res.PCA$eig[, 2][1],
                                                       digits=2),
                                                 "%)"),
                                      ylab=paste0("dim3 (",
                                                 round(res.PCA$eig[, 2][3],
                                                       digits=2),
                                                 "%)"),
                                      zlab=paste0("dim2 (",
                                                 round(res.PCA$eig[, 2][2],
                                                       digits=2),
                                                 "%)"))

                    plot3D::text3D(coordPCAind[coord.smpl, 1],
                                   coordPCAind[coord.smpl, 3],
                                   coordPCAind[coord.smpl, 2],
                                   labels=row.names(coordPCAind)[coord.smpl],
                                   add=TRUE, colkey=FALSE,
                                   cex=Cex.label, pch=20,
                                   col=as.character(ind.color)[coord.smpl])

                    ## graphics::legend("right", title=LegendTitle,
                    ##                  legend=data.color[,1],
                    ##                  pch=20, horiz=FALSE, xpd=TRUE,
                    ##                  cex=0.8, # inset=c(-0.015),
                    ##                  # cex = Cex.point*0.9
                    ##                  col=data.color[,2])

                    ##--------------------------------------------------------#
                    for(smpl in seq(from=2, to=length(LvlsPAT), by=1)){
                        Smpl.sel <- LvlsPAT[smpl]
                        IndexSmplSel <- which(Vector.patient == Smpl.sel)
                        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                        ##----------------------------------------------------#
                        if (length(coord.smpl) == 1) {
                            coord.smplf <- c(coord.smpl, NA)
                        } else {
                            coord.smplf <- coord.smpl
                        }## if (length(coord.smpl) == 1)

                        ##----------------------------------------------------#
                        plot3D::scatter3D(coordPCAind[coord.smplf, 1],
                                          coordPCAind[coord.smplf, 3],
                                          coordPCAind[coord.smplf, 2],
                                          type="b", colvar=NULL,
                                          add=TRUE, pch=20, cex=Cex.point,
                                          bty="b2", ticktype="detailed",
                                          col=col.link[coord.smpl])

                        plot3D::text3D(coordPCAind[coord.smpl, 1],
                                       coordPCAind[coord.smpl, 3],
                                       coordPCAind[coord.smpl, 2],
                                       labels=row.names(coordPCAind)[coord.smpl],
                                       add=TRUE, colkey=FALSE,
                                       cex=Cex.label, pch=20,
                                       col=as.character(ind.color)[coord.smpl])
                    }## for(smpl in 2:length(LvlsPAT))
                } else {
                    ##--------------------------------------------------------#
                    Smpl.sel <- LvlsPAT[1]
                    IndexSmplSel <- which(Vector.patient == Smpl.sel)
                    times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                    coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                    ##--------------------------------------------------------#
                    if (length(coord.smpl) == 1) {
                        coord.smplf <- c(coord.smpl, NA)
                    } else {
                        coord.smplf <- coord.smpl
                    }## if(length(coord.smpl)==1)

                    ##--------------------------------------------------------#
                    # s3d.2t<-
                    plot3D::scatter3D(coordPCAind[coord.smplf, 1],
                                      coordPCAind[coord.smplf, 3],
                                      coordPCAind[coord.smplf, 2],
                                      theta=Theta, phi=Phi, bty="b2", pch=20,
                                      type="p", cex=Cex.point, colvar=NULL,
                                      col=col.link[coord.smpl],
                                      xlim=c(min(coordPCAind[, 1]),
                                             max(coordPCAind[, 1])),
                                      ylim=c(min(coordPCAind[, 3]),
                                             max(coordPCAind[, 3])),
                                      zlim=c(min(coordPCAind[, 2]),
                                             max(coordPCAind[, 2])),
                                      xlab=paste0("dim1 (",
                                                 round(res.PCA$eig[, 2][1],
                                                       digits=2),
                                                 "%)"),
                                      ylab=paste0("dim3 (",
                                                 round(res.PCA$eig[, 2][3],
                                                       digits=2),
                                                 "%)"),
                                      zlab=paste0("dim2 (",
                                                 round(res.PCA$eig[, 2][2],
                                                       digits=2),
                                                 "%)"))
                    #
                    plot3D::text3D(coordPCAind[coord.smpl, 1],
                                   coordPCAind[coord.smpl, 3],
                                   coordPCAind[coord.smpl, 2],
                                   labels=row.names(coordPCAind)[coord.smpl],
                                   add=TRUE, colkey=FALSE,
                                   cex=Cex.label, pch=20,
                                   col=as.character(ind.color)[coord.smpl])

                    ## graphics::legend("right", title=LegendTitle,
                    ##                  legend = data.color[,1],
                    ##                  pch = 20, horiz=FALSE, xpd = TRUE,
                    ##                  cex=0.8, # inset=c(-0.015),
                    ##                  col = data.color[,2])

                    ##--------------------------------------------------------#
                    for(smpl in seq(from=2, to=length(LvlsPAT), by=1)){
                        Smpl.sel <- LvlsPAT[smpl]
                        IndexSmplSel <- which(Vector.patient == Smpl.sel)
                        times.smpl.sel <- factor(Vector.time[IndexSmplSel])
                        coord.smpl <- IndexSmplSel[order(times.smpl.sel)]

                        ##----------------------------------------------------#
                        if (length(coord.smpl) == 1) {
                            coord.smplf <- c(coord.smpl, NA)
                        } else {
                            coord.smplf <- coord.smpl
                        }## if(length(coord.smpl)==1)

                        ##----------------------------------------------------#
                        plot3D::scatter3D(coordPCAind[coord.smplf, 1],
                                          coordPCAind[coord.smplf, 3],
                                          coordPCAind[coord.smplf, 2],
                                          type="p", colvar=NULL, cex=Cex.point,
                                          add=TRUE, pch=20, bty="b2",
                                          col=col.link[coord.smpl])

                        plot3D::text3D(coordPCAind[coord.smpl, 1],
                                       coordPCAind[coord.smpl, 3],
                                       coordPCAind[coord.smpl, 2],
                                       labels=row.names(coordPCAind)[coord.smpl],
                                       add=TRUE, colkey=FALSE,
                                       cex=Cex.label, pch=20,
                                       col=as.character(ind.color)[coord.smpl])
                    }## for(smpl in 2:length(LvlsPAT))

                    ##--------------------------------------------------------#
                    data.name <- data.frame(rep("Mean",
                                                times=nrow(PCAqualiSup)),
                                            PCAqualiSup)

                    ##--------------------------------------------------------#
                    for (col in seq_len(length(levels(factor(col.link))))) {

                        Id.col <- which(col.link==levels(factor(col.link))[col])
                        Mean.t <- stats::aggregate(coordPCAind[Id.col,
                                                               c(1, 2, 3)],
                                                   list(coord.t[Id.col]), mean)
                        New.Name.3D <- apply(unique(data.name[Id.col,]), 1,
                                             paste, collapse="_")

                        plot3D::scatter3D(Mean.t[,2], Mean.t[,4], Mean.t[,3],
                                          type="b", colvar=NULL, cex=Cex.point,
                                          add=TRUE, pch=20, bty="b2",
                                          col=levels(factor(col.link))[col])
                        plot3D::text3D(Mean.t[, 2], Mean.t[, 4], Mean.t[, 3],
                                       labels=New.Name.3D, cex=Cex.label,
                                       col=levels(factor(col.link))[col],
                                       add=TRUE, colkey=FALSE, pch=20)
                    }## for(col in 1:length(levels(factor(col.link))))
                }## if(Mean.Accross.Time==FALSE & max(NbPerCond)>1)

                ##------------------------------------------------------------#
                ## 3D PCA in rgl window
                if (D3.mouvement == TRUE) {
                    plot3Drgl::plotrgl()
                }## if(D3.mouvement==TRUE)
            }## if(Plot.PCA==TRUE)

            ##----------------------------------------------------------------#
            ##----------------------------------------------------------------#
            ## s3d.2t<-grDevices::recordPlot()
            ## graphics::plot.new() ## clean up device
            ##
            ## List.plot.PCA<-vector(mode="list", length=4)
            ## names(List.plot.PCA)<-c("PCA.2D", "PCA.3D",
            ##                         "PCA.2D.temporal.links",
            ##                         "PCA.3D.temporal.links")
            ## List.plot.PCA[[1]]<-g.2DPCA
            ## List.plot.PCA[[2]]<-s3d.all
            ## List.plot.PCA[[3]]<-g.2DPCA.t
            ## List.plot.PCA[[4]]<-s3d.2t
        } else {
            ## List.plot.PCA<-vector(mode="list", length=2)
            ## names(List.plot.PCA)<-c("PCA.2D", "PCA.3D")
            ## List.plot.PCA[[1]]<-g.2DPCA
            ## List.plot.PCA[[2]]<-s3d.all
        }## if(is.null(Vector.time)==FALSE)

    } else {
        ## List.plot.PCA<-NULL
    }## if(Plot.PCA==TRUE | is.null(path.result)==FALSE)
    ##
    ## par(mar=c(5.1, 4.1, 4.1, 2.1))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output ## list(res.pca=res.PCA, List.plot.PCA=NULL)
    return(SEobj=SEresPCA)
}## PCAgraphics()
