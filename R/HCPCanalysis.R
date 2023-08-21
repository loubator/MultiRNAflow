#' @title Hierarchical clustering analysis with HCPC (Main function)
#'
#' @description The functions performs a hierarchical clustering on results
#' from a factor analysis with the R function
#' [FactoMineR::HCPC()].
#'
#' @details All results are built from the results of our function
#' [DATAnormalization()].
#'
#' The number of clusters is automatically selected by
#' [FactoMineR::HCPC()]
#' and is described in the section \code{Details} of
#' [FactoMineR::HCPC()].
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
#' @param Plot.HCPC \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, all graphs will be plotted.
#' Otherwise no graph will be plotted.
#' @param Color.Group \code{NULL} or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute
#' a color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param D3.mouvement \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the 3D PCA plots will also be plotted in a rgl window
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
#' @param Cex.label Non negative numeric value giving the size of the labels
#' associated to each point of the all PCA graphs which are not automatically
#' plotted by [FactoMineR::PCA()].
#' @param epsilon Non negative numeric value giving the length between points
#' and their labels in all PCA plots which are not automatically plotted
#' by [FactoMineR::PCA()].
#' @param path.result Character or \code{NULL}. Path to save all results.
#' If \code{path.result} contains a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.hcpc}" and a sub sub folder,
#' "1-3_HCPCanalysis_\code{Name.folder.hcpc}"
#' all results will be saved in the sub folder
#' "1_UnsupervisedAnalysis_\code{Name.folder.hcpc}/
#' 1-3_HCPCanalysis_\code{Name.folder.hcpc}".
#' Otherwise, a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.hcpc}" and/or a sub sub folder
#' "1-3_HCPCanalysis_\code{Name.folder.hcpc}"
#' will be created in \code{path.result} and all results will be saved in
#' "1_UnsupervisedAnalysis_\code{Name.folder.hcpc}/
#' 1-3_HCPCanalysis_\code{Name.folder.hcpc}".
#' If \code{NULL}, the results will not be saved in a folder.
#' \code{NULL} as default.
#' @param Name.folder.hcpc Character or \code{NULL}.
#' If \code{Name.folder.hcpc} is a character, the folder and sub folder names
#' which will contain the PCA graphs will respectively be
#' "1_UnsupervisedAnalysis_\code{Name.folder.hcpc}"
#' and "1-3_HCPCanalysis_\code{Name.folder.hcpc}".
#' Otherwise, the folder and sub folder names will respectively be
#' "1_UnsupervisedAnalysis" and "1-3_HCPCanalysis".
#'
#' @return The function returns
#' * a dendrogram (also called hierarchical tree) using the function
#' [factoextra::fviz_dend()]
#' * one 2D PCA and two 3D PCA produced by the function
#' [PCArealization()]
#' where samples are colored with different colors for different clusters.
#' The two 3D PCA graphs are identical but one of them will be opened
#' in a rgl window
#' (see [plot3Drgl::plotrgl()])
#' allowing to interactively rotate and zoom.
#' The interactive 3D graph will be plotted only if \code{D3.mouvement=TRUE}.
#' * A graph indicating for each sample, its cluster and
#' the time and/or biological condition associated to the sample.
#' * the outputs of
#' [FactoMineR::HCPC()].
#'
#' @seealso The function calls the functions
#' [PCArealization()] and
#' [FactoMineR::HCPC()].
#' The function
#' [FactoMineR::HCPC()]
#' will take as input the output of
#' [PCArealization()].
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom FactoMineR HCPC
#' @importFrom scales hue_pal
#' @importFrom ggsci pal_jco
#' @importFrom RColorBrewer brewer.pal
#' @importFrom factoextra fviz_dend fviz_cluster
#' @importFrom graphics legend
#' @importFrom plot3D scatter3D text3D
#' @importFrom plot3Drgl plotrgl
#' @importFrom grDevices pdf dev.off
#' @importFrom ggplot2 ggplot aes ylab ggtitle theme scale_y_continuous guides
#' scale_color_manual scale_fill_manual coord_flip geom_bar guide_legend
#' element_text guide_axis theme_minimal
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
#' ##-------------------------------------------------------------------------#
#' resHCPCanalysis <- HCPCanalysis(SEresNorm=resNorm,
#'                                 DATAnorm=TRUE,
#'                                 sample.deletion=NULL,
#'                                 Supp.del.sample=FALSE,
#'                                 gene.deletion=NULL,
#'                                 Plot.HCPC=TRUE,
#'                                 Color.Group=NULL,
#'                                 Phi=25, Theta=140,
#'                                 Cex.point=1, Cex.label=0.6, epsilon=0.4,
#'                                 D3.mouvement=FALSE,
#'                                 path.result=NULL,
#'                                 Name.folder.hcpc=NULL)

HCPCanalysis <- function(SEresNorm,
                         DATAnorm=TRUE,
                         gene.deletion=NULL,
                         sample.deletion=NULL,
                         Supp.del.sample=FALSE,
                         Plot.HCPC=TRUE,
                         Color.Group=NULL,
                         Phi=25, Theta=140, epsilon=0.2,
                         Cex.point=0.7, Cex.label=0.7,
                         D3.mouvement=FALSE,
                         path.result=NULL,
                         Name.folder.hcpc=NULL) {
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check deletion
    if (!is.null(path.result)) {
        if (!is.character(path.result)) {
            stop("'path.result' must be NULL or a character.")
        }## if (!is.character(path.result))
    }## if (!is.null(path.result))

    if (!is.null(Name.folder.hcpc)) {
        if (!is.character(Name.folder.hcpc)) {
            stop("'Name.folder.hcpc' must be NULL or a character.")
        }## if (!is.character(Name.folder.hcpc))
    }## if (!is.null(Name.folder.hcpc))

    if (!isTRUE(D3.mouvement) & !isFALSE(D3.mouvement)) {
        stop("'D3.mouvement' must be TRUE or FALSE.")
    }## if (!isTRUE(D3.mouvement) & !isFALSE(D3.mouvement))

    if (!isTRUE(Plot.HCPC) & !isFALSE(Plot.HCPC)) {
        stop("'Plot.HCPC' must be TRUE or FALSE.")
    }## if (!isTRUE(Plot.HCPC) & !isFALSE(Plot.HCPC))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## To avoid "no visible binding for global variable" with devtools::check()
    Samples <- len <- NULL

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## PCA preprocessing
    SEresPCA <- PCArealization(SEresNorm=SEresNorm,
                               DATAnorm=DATAnorm,
                               sample.deletion=sample.deletion,
                               gene.deletion=gene.deletion,
                               Supp.del.sample=Supp.del.sample)

    resPCA <- S4Vectors::metadata(SEresPCA)$PCA$res.pca
    listFCTRS <- S4Vectors::metadata(SEresPCA)$PCA$List.Factors

    Vector.time <- listFCTRS$Vector.time
    Vector.group <- listFCTRS$Vector.group
    Vector.patient <- listFCTRS$Vector.patient
    LvlsPAT <- levels(factor(Vector.patient))

    ## PCAqualiSup <- resPCA$call$quali.sup$quali.sup
    ## coordPCAind <- resPCA$ind$coord

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Folder creation if no existence
    if (is.null(Name.folder.hcpc)) {
        Name.folder.hcpc <- ""
        SubFolder.name <- "1_UnsupervisedAnalysis"
    } else {
        Name.folder.hcpc <- paste0("_", Name.folder.hcpc)
        SubFolder.name <- paste0("1_UnsupervisedAnalysis", Name.folder.hcpc)
    }## if(is.null(Name.folder.hcpc)==TRUE)

    if (!is.null(path.result)) {
        if (!SubFolder.name%in%dir(path=path.result)) {
            print("Folder creation")
            dir.create(path=file.path(path.result, SubFolder.name))
            path.result.f <- file.path(path.result, SubFolder.name)
        } else {
            path.result.f <- file.path(path.result, SubFolder.name)
        }## if(SubFolder.name%in%dir(path = path.result)==FALSE)
    } else {
        path.result.f <- NULL
    }## if(is.null(path.result)==FALSE)


    if (!is.null(path.result)) {
        SubFolder.name <- paste0("1-3_HCPCanalysis", Name.folder.hcpc)

        if (!SubFolder.name%in%dir(path=path.result.f)) {
            dir.create(path=file.path(path.result.f, SubFolder.name))
            path.result.new <- file.path(path.result.f, SubFolder.name)
        } else {
            path.result.new <- file.path(path.result.f, SubFolder.name)
        }## if(SubFolder.name%in%dir(path = path.result.f)==FALSE)

    } else {
        path.result.new <- NULL
    }## if(is.null(path.result)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## 2D PCA and clustering
    res.hcpc <- FactoMineR::HCPC(resPCA, nb.clust=-1, method="ward",
                                 graph=FALSE, consol=FALSE)

    if (isTRUE(Plot.HCPC) | !is.null(path.result)) {
        ##--------------------------------------------------------------------#
        ## Color for each cluster
        MyColclust <- c(ggsci::pal_jco("default")(4), "#35B600",
                        ggsci::pal_jco("default")(10)[5:10])
        NbClust <- res.hcpc$call$t$nb.clust

        if (NbClust > 10) {
            MyColclust2 <- scales::hue_pal(direction=-1)(NbClust-8)
            col.clust <- c(MyColclust, MyColclust2)[seq_len(NbClust)]
        } else {
            col.clust <- as.character(MyColclust)[seq_len(NbClust)]
        }## if(NbClust>10)

        colInd <- res.hcpc$data.clust$clust
        levels(colInd) <- col.clust

        TreeOrder <- res.hcpc$call$t$tree$order
        ColorderTree <- unique(as.numeric(res.hcpc$call$X$clust[TreeOrder]))

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## Shape for each cluster
        if (NbClust > 16) {
            ShapeFvizClust <- 16
        } else {
            ShapeFvizClust <- NULL
        }## if (NbClust > 16)

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## Dendogram and PCA graph with clustering results
        options(warn=-1, ggrepel.max.overlaps=35)

        g.2DPCA.clust <- factoextra::fviz_cluster(res.hcpc, repel=TRUE,
                                                  show.clust.cent=TRUE,
                                                  ggtheme=ggplot2::theme_minimal(),
                                                  main="PCA plot with HCPC clusters")+
            ggplot2::scale_colour_manual('Cluster',
                                         breaks=seq_len(NbClust),
                                         values=levels(colInd))+
            ggplot2::scale_fill_manual('Cluster', values=levels(colInd))+
            ## ggplot2::scale_shape_discrete(name='Cluster')+ ## fill/colour
            ggplot2::guides(colour=ggplot2::guide_legend(title="Cluster",
                                                         override.aes=ggplot2::aes(label="")))

        ##--------------------------------------------------------------------#
        ## Shape for each cluster
        if (NbClust > 6 & NbClust <= 16) {
            g.2DPCA.clust <- g.2DPCA.clust +
                ggplot2::scale_shape_manual(name='Cluster',
                                            values=c(16, 17, 15, 3, 7, 8, 18,
                                                     4, 10, 11, 12,
                                                     13, 14, 5,
                                                     6, 9)[seq_len(NbClust)])
        } else {
            g.2DPCA.clust <- g.2DPCA.clust +
                ggplot2::scale_shape_discrete(name='Cluster')
        }## if(NbClust>6)

        ##--------------------------------------------------------------------#
        OldLabel <- res.hcpc$call$t$tree$labels
        res.hcpc$call$t$tree$labels <- paste0(" ", OldLabel)
        InerGain <- res.hcpc$call$t$inert.gain[1]/5

        Dendo.clust <- factoextra::fviz_dend(res.hcpc,
                                             main="Dendrogram (Ward distance)",
                                             xlab="",
                                             ylab="Distance",
                                             cex=0.65,
                                             palette=levels(colInd)[ColorderTree],
                                             rect=TRUE, rect_fill=FALSE,
                                             rect_lty=7,
                                             rect_border="azure3",
                                             horiz=TRUE,
                                             labels_track_height=InerGain)

        res.hcpc$call$t$tree$labels <- OldLabel

        options(warn=0, ggrepel.max.overlaps=10)

        ##--------------------------------------------------------------------#
        ## Data describing the distribution of samples
        ## among cluster, times, groups
        NbGene <- nrow(resPCA$var$coord)
        RepSample <- rep(row.names(res.hcpc$data.clust),
                         times=ncol(res.hcpc$data.clust) - NbGene)
        RepSample <- factor(RepSample,
                            levels=row.names(res.hcpc$call$X)[TreeOrder])

        ClustCharac <- factor(paste0("Cluster.",
                                     as.numeric(res.hcpc$data.clust$clust)),
                              levels=paste0("Cluster.", seq_len(NbClust)))

        Index.Factor <- seq_len(ncol(res.hcpc$data.clust) - NbGene - 1)
        Attribute <- c(as.character(ClustCharac),
                       as.character(unlist(res.hcpc$data.clust[,
                                                               Index.Factor])))

        FactorLevels <- apply(data.frame(res.hcpc$data.clust[,Index.Factor]),
                              2, function(x) levels(factor(x)))

        Attribute <- factor(Attribute,
                            levels=c(as.character(unlist(FactorLevels)),
                                     levels(ClustCharac)))

        DataClustFactor <- data.frame(Samples=RepSample,
                                      Attribute=Attribute,
                                      len=rep(1, times=length(RepSample)))

        ##--------------------------------------------------------------------#
        ## Graph preprocessing
        DistriNameFile <- "SampleDistribution_Clusters_Times_Groups"
        gDistriTitle <- "Links between sample information and clusters"
        NbAttribute <- length(c(as.character(unlist(FactorLevels)),
                                levels(ClustCharac)))
        NbSample <- length(row.names(res.hcpc$data.clust))

        ##--------------------------------------------------------------------#
        if (NbAttribute > NbSample) {
            NcolLegend <- 2
        } else {
            NcolLegend <- 1
        }## if (NbAttribute > NbSample)

        Color.f <- c()

        if (!is.null(Vector.group)) {
            Glevels <- levels(factor(Vector.group))

            if (is.null(Color.Group)) {
                MypaletteG <- c(RColorBrewer::brewer.pal(8, "Dark2"),
                                RColorBrewer::brewer.pal(8, "Set2"))

                if (length(Glevels) > 16) {
                    MypalG <- scales::hue_pal(l=90)(seq_len(length(Glevels)-1))
                    MypaletteG <- c(MypaletteG, MypalG)
                }## if(length(Glevels)>16)

                MypaletteG <- MypaletteG[seq_len(length(Glevels))]
                Color.Group <- data.frame(Name=Glevels, Col=MypaletteG)

            } else {
                Id.LevelCol.G <- order(Color.Group[, 1])
                Color.Group <- data.frame(Name=Glevels,
                                          Col=Color.Group[Id.LevelCol.G, 2])
            }## if(is.null(Color.Group)==TRUE)

            Color.f <- c(Color.f, Color.Group$Col)
        } else {
            DistriNameFile <- gsub("_Groups", "", DistriNameFile, fixed=TRUE)
            gDistriTitle <- gsub("groups, ", "", gDistriTitle, fixed=TRUE)
        }## if(is.null(Vector.group)==FALSE)


        if (!is.null(Vector.time)) {
            Color.Time <- NULL
            Tlevels <- levels(factor(Vector.time))

            if (is.null(Color.Time)) {
                MypaletteT <- c("#737373",# "#252525"
                                scales::hue_pal()(length(Tlevels)-1))
                Color.Time <- data.frame(Name=Tlevels, Col=MypaletteT)
            } else {
                Id.LevelColT <- order(Color.Time[, 1])
                Color.Time <- data.frame(Name=Tlevels,
                                         Col=Color.Time[Id.LevelColT, 2])
            }## if(is.null(Color.Time)==TRUE)

            Color.f <- c(Color.f, Color.Time$Col)
        } else {
            DistriNameFile <- gsub("_Times", "", DistriNameFile, fixed=TRUE)
            gDistriTitle <- gsub(", times", "", gDistriTitle, fixed=TRUE)
        }## if(is.null(Vector.time) == FALSE)

        Color.f <- c(Color.f, col.clust)

        ##--------------------------------------------------------------------#
        gDistribution <- ggplot2::ggplot(data=DataClustFactor,
                                         aes(x=Samples, y=len, fill=Attribute,
                                             color=Attribute))+
            ggplot2::geom_bar(stat="identity")+
            ggplot2::coord_flip()+
            ggplot2::scale_fill_manual(values=Color.f)+
            ggplot2::scale_color_manual(values=Color.f)+
            ggplot2::scale_y_continuous(label=c("Cluster", "Time", "Group"),
                                        breaks=c(0.5, 1.5, 2.5),
                                        guide=ggplot2::guide_axis(angle=45))+
            ggplot2::ylab("")+
            ggplot2::ggtitle(gDistriTitle)+
            ggplot2::theme(axis.text.x=ggplot2::element_text(face="bold"),
                           legend.title=ggplot2::element_text(face="bold",
                                                              size=ggplot2::rel(0.8)),
                           legend.text=ggplot2::element_text(size=ggplot2::rel(0.6)))+
            ggplot2::guides(fill=ggplot2::guide_legend(ncol=NcolLegend))

        ##--------------------------------------------------------------------#
        ## Save graph
        if(!is.null(path.result)){

            DendoTitle<-paste0("Dendogram_HCPC", Name.folder.hcpc, ".pdf")
            grDevices::pdf(file=file.path(path.result.new, DendoTitle),
                           width=11, height=8)
            print(Dendo.clust)
            grDevices::dev.off()

            PCA2dTitle<-paste0("PCA2d_HCPC", Name.folder.hcpc, ".pdf")
            grDevices::pdf(file=file.path(path.result.new, PCA2dTitle),
                           width=11, height=8)
            print(g.2DPCA.clust) ## print(g.2DPCA.clust)
            grDevices::dev.off()

            DistiTitle<-paste0(DistriNameFile, Name.folder.hcpc, ".pdf")
            grDevices::pdf(file=file.path(path.result.new, DistiTitle),
                           width=11, height=8)
            print(gDistribution)
            grDevices::dev.off()

        } else {

            if (Plot.HCPC == TRUE) {
                print(Dendo.clust)
                print(g.2DPCA.clust)
                print(gDistribution)
            }## if(Plot.HCPC==TRUE)

        }## if(is.null(path.result)==FALSE)

        ##--------------------------------------------------------------------#
        ## 3D PCA colored by cluster
        data.3D <- resPCA$ind$coord[, c(1, 2, 3)]

        if (!is.null(path.result)) {
            PCA3dTitle<-paste0("PCA3d_HCPC", Name.folder.hcpc, ".pdf")
            grDevices::pdf(file=file.path(path.result.new, PCA3dTitle),
                           width=11, height=8)

            plot3D::scatter3D(data.3D[, 1], data.3D[, 3], data.3D[, 2],
                              col=as.character(colInd), ticktype="detailed",
                              pch=20, cex=Cex.point, colvar=NULL, bty = "b2",
                              theta=Theta, phi=Phi, d=2, epsilon=epsilon,
                              main = "3D PCA plot with HCPC clusters",
                              xlab=paste0("dim1 (",
                                          round(resPCA$eig[, 2][1], digits=2),
                                         "%)"),
                              ylab=paste0("dim3 (",
                                          round(resPCA$eig[, 2][3], digits=2),
                                         "%)"),
                              zlab=paste0("dim2 (",
                                          round(resPCA$eig[, 2][2], digits=2),
                                         "%)"))

            plot3D::text3D(data.3D[, 1] + epsilon,
                           data.3D[, 3] + epsilon,
                           data.3D[, 2] + epsilon,
                           labels=rownames(data.3D), cex=Cex.label, font=2,
                           col=as.character(colInd), add=TRUE, colkey=FALSE)

            graphics::legend("right", title="Cluster",
                             legend=seq_len(NbClust), pch=20, horiz=FALSE,
                             xpd=TRUE,  col=col.clust,
                             cex=Cex.point*0.8, inset=c(-0.015))## -0.15

            ## legend=#paste(rep("Cluster",times=NbClust), 1:NbClust),
            ## col=unique(col.clust)[order(unique(res.hcpc$data.clust$clust))],
            grDevices::dev.off()
        }## if(is.null(path.result)==FALSE)

        if (Plot.HCPC == TRUE) {
            plot3D::scatter3D(data.3D[, 1], data.3D[, 3], data.3D[, 2],
                              col=as.character(colInd),
                              epsilon=epsilon, cex=Cex.point, colvar=NULL, d=2,
                              ticktype="detailed", theta=Theta, phi=Phi, pch=20,
                              main="3D PCA plot with HCPC clusters", bty="b2",
                              xlab=paste0("dim1 (",
                                          round(resPCA$eig[, 2][1], digits=2),
                                          "%)"),
                              ylab=paste0("dim3 (",
                                          round(resPCA$eig[, 2][3], digits=2),
                                          "%)"),
                              zlab=paste0("dim2 (",
                                          round(resPCA$eig[, 2][2], digits=2),
                                          "%)"))

            plot3D::text3D(data.3D[, 1] + epsilon,
                           data.3D[, 3] + epsilon,
                           data.3D[, 2] + epsilon,
                           labels = rownames(data.3D), cex=Cex.label, font=2,
                           col=as.character(colInd), add=TRUE, colkey=FALSE)

            graphics::legend("right", title="Cluster",
                             legend=seq_len(NbClust),
                             pch=20, horiz=FALSE, xpd=TRUE,
                             col=col.clust,
                             cex=Cex.point*0.8, inset=c(-0.015))

            ## PCA.3D<-grDevices::recordPlot()
            ## graphics::plot.new() ## clean up device

            ##----------------------------------------------------------------#
            ## 3D PCA in rgl window
            if (D3.mouvement == TRUE) {
                plot3Drgl::plotrgl()
            }## if(D3.mouvement == TRUE)
        }## if(Plot.HCPC==TRUE)

        ##--------------------------------------------------------------------#
        List.plot.hcpc <- vector(mode="list", length=4)
        names(List.plot.hcpc) <- c("Dendrogram", "Sample.Information.Clusters",
                                   "PCA2D.clusters", "PCA3D.clusters")
        List.plot.hcpc[[1]] <- Dendo.clust
        List.plot.hcpc[[2]] <- gDistribution
        List.plot.hcpc[[3]] <- g.2DPCA.clust
        ## List.plot.hcpc[[4]] <- PCA.3D
    } else {
        List.plot.hcpc <- NULL
    }## if(Plot.HCPC==TRUE)

    ##------------------------------------------------------------------------#
    FactorCluster <- data.frame(Obs=row.names(resPCA$call$quali.sup$quali.sup),
                                resPCA$call$quali.sup$quali.sup,
                                Clust=res.hcpc$data.clust$clust)
    colnames(FactorCluster) <- gsub("Quali.Sup.", "",
                                    x=colnames(FactorCluster),
                                    fixed=TRUE)

    ##------------------------------------------------------------------------#
    if (!is.null(path.result)) {
        FileTitle <- paste0("Results", res.hcpc$call$t$nb.clust,
                            "Clusters", "_", "HCPC", Name.folder.hcpc, ".csv")
        utils::write.table(FactorCluster,
                           file=file.path(path.result.new, FileTitle),
                           sep=";", row.names=FALSE)
    }## if(is.null(path.result)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## SE HCPC
    listHCPC <- list(Res.hcpc=res.hcpc, Samples.FactorCluster=FactorCluster)
    NBcol <- length(unlist(S4Vectors::metadata(SEresPCA)$colDataINFO))

    SEprepHCPC <- SEresPCA
    SummarizedExperiment::colData(SEprepHCPC)$HCPC.clust <- FactorCluster$Clust
    S4Vectors::metadata(SEprepHCPC)$HCPC <- listHCPC
    S4Vectors::metadata(SEprepHCPC)$colDataINFO$colINFOclusterHCPC <- NBcol+1

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## List.plot.HCPC=List.plot.hcpc
    return(SEobj=SEprepHCPC)
}## HCPCanalysis()
