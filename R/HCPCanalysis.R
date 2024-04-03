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
#' @param Plot.HCPC \code{TRUE} or \code{FALSE}. \code{FALSE} as default.
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
#' @param motion3D \code{TRUE} or \code{FALSE}.
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
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEresNorm} with the outputs from the function
#' [FactoMineR::HCPC()],
#' (saved in the metadata \code{Results[[1]][[3]]} of \code{SEresNorm})
#' * a dendrogram (also called hierarchical tree) using the function
#' [factoextra::fviz_dend()]
#' * one 2D PCA and two 3D PCA produced by the function
#' [PCAgraphics()]
#' where samples are colored with different colors for different clusters.
#' The two 3D PCA graphs are identical but one of them will be opened
#' in a rgl window
#' (see [plot3Drgl::plotrgl()])
#' allowing to interactively rotate and zoom.
#' The interactive 3D graph will be plotted only if \code{motion3D=TRUE}.
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
#' @importFrom factoextra fviz_dend fviz_cluster
#' @importFrom plot3D scatter3D text3D
#' @importFrom plot3Drgl plotrgl
#' @importFrom ggplotify as.ggplot
#' @importFrom grDevices dev.cur pdf dev.control dev.off dev.set recordPlot
#' @importFrom graphics legend plot.new
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
#' ##------------------------------------------------------------------------##
#' resHCPCanalysis <- HCPCanalysis(SEresNorm=resNorm,
#'                                 DATAnorm=TRUE,
#'                                 sample.deletion=NULL,
#'                                 gene.deletion=NULL,
#'                                 Plot.HCPC=TRUE,
#'                                 Color.Group=NULL,
#'                                 Phi=25, Theta=140,
#'                                 Cex.point=1, Cex.label=0.6, epsilon=0.4,
#'                                 motion3D=FALSE,
#'                                 path.result=NULL,
#'                                 Name.folder.hcpc=NULL)

HCPCanalysis <- function(SEresNorm,
                         DATAnorm=TRUE,
                         gene.deletion=NULL,
                         sample.deletion=NULL,
                         Plot.HCPC=FALSE,
                         Color.Group=NULL,
                         Phi=25, Theta=140, epsilon=0.2,
                         Cex.point=0.7, Cex.label=0.7,
                         motion3D=FALSE,
                         path.result=NULL,
                         Name.folder.hcpc=NULL) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    resErr <- ErrHCPCanalysis(SEresNorm=SEresNorm,
                              DATAnorm=DATAnorm,
                              gene.deletion=gene.deletion,
                              sample.deletion=sample.deletion,
                              Plot.HCPC=Plot.HCPC,
                              Color.Group=Color.Group,
                              Phi=Phi, Theta=Theta, epsilon=epsilon,
                              Cex.point=Cex.point, Cex.label=Cex.label,
                              motion3D=motion3D,
                              path.result=path.result,
                              Name.folder.hcpc=Name.folder.hcpc)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## To avoid "no visible binding for global variable" with devtools::check()
    Samples <- len <- NULL

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## PCA preprocessing
    SEresPCA <- PCArealization(SEresNorm=SEresNorm,
                               DATAnorm=DATAnorm,
                               sample.deletion=sample.deletion,
                               gene.deletion=gene.deletion,
                               Supp.del.sample=FALSE)

    resPCA <- S4Vectors::metadata(SEresPCA)$Results[[1]][[2]]$PCAresults
    listFCTRS <- S4Vectors::metadata(SEresPCA)$Results[[1]][[2]]$List.Factors

    Vector.time <- listFCTRS$Vector.time
    Vector.group <- listFCTRS$Vector.group
    Vector.patient <- listFCTRS$Vector.patient
    LvlsPAT <- levels(factor(Vector.patient))
    NbGene <- nrow(resPCA$var$coord)

    ## PCAqualiSup <- resPCA$call$quali.sup$quali.sup
    ## coordPCAind <- resPCA$ind$coord
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Folder creation if no existence
    if (is.null(Name.folder.hcpc)) {
        Name.folder.hcpc <- ""
        SubFolder.name <- "1_UnsupervisedAnalysis"
    } else {
        Name.folder.hcpc <- paste0("_", Name.folder.hcpc)
        SubFolder.name <- paste0("1_UnsupervisedAnalysis", Name.folder.hcpc)
    }## if(is.null(Name.folder.hcpc)==TRUE)

    SubFolderHCPC <- paste0("1-3_HCPCanalysis", Name.folder.hcpc)

    if (!is.null(path.result)) {
        if (!SubFolder.name%in%dir(path=path.result)) {
            print("Folder creation")
            dir.create(path=file.path(path.result, SubFolder.name))
        }## if(SubFolder.name%in%dir(path = path.result)==FALSE)
        path.result.f <- file.path(path.result, SubFolder.name)
    } else {
        path.result.f <- NULL
    }## if(is.null(path.result)==FALSE)

    if (!is.null(path.result)) {
        if (!SubFolderHCPC%in%dir(path=path.result.f)) {
            dir.create(path=file.path(path.result.f, SubFolderHCPC))
        }## if(SubFolderHCPC%in%dir(path = path.result.f)==FALSE)
        path.result.new <- file.path(path.result.f, SubFolderHCPC)
    } else {
        path.result.new <- NULL
    }## if(is.null(path.result)==FALSE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## 2D PCA and clustering
    res.hcpc <- FactoMineR::HCPC(resPCA, nb.clust=-1, method="ward",
                                 graph=FALSE, consol=FALSE)

    dataHCPCclust <- res.hcpc$data.clust
    Index.Factor <- seq_len(ncol(dataHCPCclust) - NbGene - 1)

    dataHCPCclust_clust <- dataHCPCclust$clust
    HCPCsample <- row.names(dataHCPCclust)
    NbClust <- res.hcpc$call$t$nb.clust

    ##-----------------------------------------------------------------------##
    ## Color for each cluster
    col.clust <- myPaletteHCPC(Nclust=NbClust)
    colInd <- dataHCPCclust_clust
    levels(colInd) <- col.clust

    TreeOrder <- res.hcpc$call$t$tree$order
    ColorderTree <- unique(as.numeric(res.hcpc$call$X$clust[TreeOrder]))

    ##-----------------------------------------------------------------------##
    ## Shape for each cluster
    ## if(NbClust>16){ShapeFvizClust<-16}else{ShapeFvizClust<-NULL}
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Dendogram and PCA graph with clustering results
    options(warn=-1, ggrepel.max.overlaps=35)

    title2DHCPC <- "PCA plot with HCPC clusters"
    gl_overaes <- ggplot2::aes(label="")

    g2DPCAclust <- factoextra::fviz_cluster(res.hcpc, repel=TRUE,
                                            show.clust.cent=TRUE,
                                            ggtheme=ggplot2::theme_minimal(),
                                            main=title2DHCPC)+
        ggplot2::scale_colour_manual('Cluster',
                                     breaks=seq_len(NbClust),
                                     values=levels(colInd))+
        ggplot2::scale_fill_manual('Cluster', values=levels(colInd))+
        ## ggplot2::scale_shape_discrete(name='Cluster')+ ## fill/colour
        ggplot2::guides(colour=ggplot2::guide_legend(title="Cluster",
                                                     override.aes=gl_overaes))

    ##-----------------------------------------------------------------------##
    ## Shape for each cluster
    Shape_val <- c(16, 17, 15, 3, 7, 8, 18, 4, 10, 11, 12, 13, 14, 5, 6, 9)
    Shape_val <- Shape_val[seq_len(NbClust)]

    if (NbClust > 6 & NbClust <= 16) {
        g2DPCAclust <- g2DPCAclust +
            ggplot2::scale_shape_manual(name='Cluster', values=Shape_val)
    } else {
        g2DPCAclust <- g2DPCAclust +
            ggplot2::scale_shape_discrete(name='Cluster')
    }## if (NbClust > 6 & NbClust <= 16)

    ##-----------------------------------------------------------------------##
    OldLabel <- res.hcpc$call$t$tree$labels
    res.hcpc$call$t$tree$labels <- paste0(" ", OldLabel)
    InerGain <- res.hcpc$call$t$inert.gain[1]/5

    clustDendro <- factoextra::fviz_dend(res.hcpc,
                                         xlab="", ylab="Distance",
                                         main="Dendrogram (Ward distance)",
                                         palette=levels(colInd)[ColorderTree],
                                         rect=TRUE, rect_fill=FALSE, rect_lty=7,
                                         rect_border="azure3", horiz=TRUE,
                                         labels_track_height=InerGain)

    res.hcpc$call$t$tree$labels <- OldLabel

    options(warn=0, ggrepel.max.overlaps=10)

    ##-----------------------------------------------------------------------##
    ## Data describing the distribution of samples
    ## among cluster, times, groups
    RepSample <- rep(HCPCsample, times=ncol(dataHCPCclust) - NbGene)
    RepSample <- factor(RepSample, levels=row.names(res.hcpc$call$X)[TreeOrder])

    ClustCharac <- factor(paste0("Cluster.", as.numeric(dataHCPCclust_clust)),
                          levels=paste0("Cluster.", seq_len(NbClust)))

    Attribute <- c(as.character(ClustCharac),
                   as.character(unlist(dataHCPCclust[, Index.Factor])))

    FactorLevels <- apply(data.frame(dataHCPCclust[, Index.Factor]),
                          2, function(x) levels(factor(x)))

    LvlsAttribute <- c(as.character(unlist(FactorLevels)), levels(ClustCharac))
    Attribute <- factor(Attribute, levels=LvlsAttribute)

    dataClustFactor <- data.frame(Samples=RepSample,
                                  Attribute=Attribute,
                                  len=rep(1, times=length(RepSample)))

    ##-----------------------------------------------------------------------##
    ## Graph preprocessing
    DistriNameFile <- "SampleDistribution_Clusters_Times_Groups"
    gDistriTitle <- "Links between sample information and clusters"
    NbAttribute <- length(LvlsAttribute)
    NbSample <- length(HCPCsample)

    if (NbAttribute > NbSample) {
        NcolLegend <- 2
    } else {
        NcolLegend <- 1
    }## if (NbAttribute > NbSample)

    ##-----------------------------------------------------------------------##
    Color.f <- c()
    lb_gd <- c(1)

    if (!is.null(Vector.group)) {
        lb_gd <- c(lb_gd, 3)
        Glevels <- levels(factor(Vector.group))

        if (is.null(Color.Group)) {
            MypaletteG <- myPaletteBC(Nbc=length(Glevels))
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
    }## if(!is.null(Vector.group)==FALSE)


    if (!is.null(Vector.time)) {
        lb_gd <- sort(c(lb_gd, 2))

        Tlevels <- levels(factor(Vector.time))
        MypaletteT <- myPaletteT(Nt=length(Tlevels))
        Color.Time <- data.frame(Name=Tlevels, Col=MypaletteT)

        Color.f <- c(Color.f, Color.Time$Col)

        ## Color.Time <- NULL
        ## if (is.null(Color.Time)) { ## "#252525"
        ## } else {
        ##     Id.LevelColT <- order(Color.Time[, 1])
        ##     Color.Time <- data.frame(Name=Tlevels,
        ##                              Col=Color.Time[Id.LevelColT, 2])
        ## }## if (is.null(Color.Time))
    } else {
        DistriNameFile <- gsub("_Times", "", DistriNameFile, fixed=TRUE)
        gDistriTitle <- gsub(", times", "", gDistriTitle, fixed=TRUE)
    }## if(is.null(Vector.time) == FALSE)

    Color.f <- c(Color.f, col.clust)

    ##-----------------------------------------------------------------------##
    lgTitle_gd <- ggplot2::element_text(face="bold", size=ggplot2::rel(0.8))
    lgText_gd <- ggplot2::element_text(size=ggplot2::rel(0.6))

    name_gd <- c("Cluster", "Time", "Biological condition")[lb_gd]
    breaks_gd <- c(0.5, 1.5, 2.5)[seq_len(length(lb_gd))]

    gDistribution <- ggplot2::ggplot(data=dataClustFactor,
                                     aes(x=Samples, y=len,
                                         fill=Attribute, color=Attribute)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values=Color.f) +
        ggplot2::scale_color_manual(values=Color.f) +
        ggplot2::scale_y_continuous(label=name_gd,
                                    breaks=breaks_gd,
                                    guide=ggplot2::guide_axis(angle=45)) +
        ggplot2::ylab("") +
        ggplot2::ggtitle(gDistriTitle) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(face="bold"),
                       legend.title=lgTitle_gd, legend.text=lgText_gd) +
        ggplot2::guides(fill=ggplot2::guide_legend(ncol=NcolLegend))

    ##-----------------------------------------------------------------------##
    ## Save graph
    DendoTitle <- paste0("Dendogram_HCPC", Name.folder.hcpc, ".pdf")
    PCA2dTitle <- paste0("PCA2D_withHCPCclustering", Name.folder.hcpc, ".pdf")
    DistiTitle <- paste0("Cluster_SampleDistribution", Name.folder.hcpc, ".pdf")

    if (!is.null(path.result)) {
        grDevices::pdf(file=file.path(path.result.new, DendoTitle),
                       width=11, height=8)
        print(clustDendro)
        grDevices::dev.off()

        grDevices::pdf(file=file.path(path.result.new, PCA2dTitle),
                       width=11, height=8)
        print(g2DPCAclust)
        grDevices::dev.off()

        grDevices::pdf(file=file.path(path.result.new, DistiTitle),
                       width=11, height=8)
        print(gDistribution)
        grDevices::dev.off()

    }## if(is.null(path.result)==FALSE)

    if (isTRUE(Plot.HCPC)) {
        print(clustDendro)
        print(g2DPCAclust)
        print(gDistribution)
    }## if(isTRUE(Plot.HCPC))

    ##-----------------------------------------------------------------------##
    ## 3D PCA colored by cluster
    data3D <- resPCA$ind$coord[, c(1, 2, 3)]

    # res3Dhcpc <- HCPCpca3D(PCAcoord3D=data3D, matPCAeig=resPCA$eig,
    #                        colorINDfactor=colInd,
    #                        Phi=Phi, Theta=Theta, epsilon=epsilon,
    #                        Cex.point=Cex.point, Cex.label=Cex.label)

    res3Dhcpc <- ggplotify::as.ggplot(
        function() HCPCpca3D(PCAcoord3D=data3D,
                             matPCAeig=resPCA$eig,
                             colorINDfactor=colInd,
                             Phi=Phi, Theta=Theta,
                             epsilon=epsilon,
                             Cex.point=Cex.point,
                             Cex.label=Cex.label)
    )

    PCA3dTitle <- paste0("PCA3D_withHCPCclustering", Name.folder.hcpc, ".pdf")

    if (!is.null(path.result)) {
        grDevices::pdf(file=file.path(path.result.new, PCA3dTitle),
                       width=11, height=8)
        print(res3Dhcpc)
        grDevices::dev.off()
    }## if(is.null(path.result)==FALSE)

    if (isTRUE(Plot.HCPC)) {
        graphics::plot.new() ## clean up device
        print(res3Dhcpc)

        ## 3D PCA in rgl window
        if (isTRUE(motion3D)) {
            plot3Drgl::plotrgl()
        }## if (isTRUE(motion3D))
    }## if(isTRUE(Plot.HCPC))

    ##-----------------------------------------------------------------------##
    listHCPCplot <- vector(mode="list", length=4)
    names(listHCPCplot) <- c("Dendrogram",
                             "Cluster_SampleDistribution",
                             "PCA2DclustersHCPC", "PCA3DclustersHCPC")
    listHCPCplot[[1]] <- clustDendro
    listHCPCplot[[2]] <- gDistribution
    listHCPCplot[[3]] <- g2DPCAclust
    listHCPCplot[[4]] <- res3Dhcpc

    ##-----------------------------------------------------------------------##
    FactorCluster <- data.frame(Obs=row.names(resPCA$call$quali.sup$quali.sup),
                                resPCA$call$quali.sup$quali.sup,
                                Clust=dataHCPCclust$clust)
    colnames(FactorCluster) <- gsub("Quali.Sup.", "",
                                    x=colnames(FactorCluster),
                                    fixed=TRUE)

    ##-----------------------------------------------------------------------##
    FileTitle <- paste0("HCPCresults_", res.hcpc$call$t$nb.clust,
                        "Clusters", Name.folder.hcpc, ".csv")

    if (!is.null(path.result)) {
        utils::write.table(FactorCluster,
                           file=file.path(path.result.new, FileTitle),
                           sep=";", row.names=FALSE)
    }## if(is.null(path.result)==FALSE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## SE HCPC
    listHCPC <- list(resHCPC=res.hcpc, Samples.FactorCluster=FactorCluster)
    NBcol <- length(unlist(S4Vectors::metadata(SEresPCA)$colDataINFO))

    SEprepHCPC <- SEresPCA
    SummarizedExperiment::colData(SEprepHCPC)$HCPC.clust <- FactorCluster$Clust
    S4Vectors::metadata(SEprepHCPC)$colDataINFO$colINFOclusterHCPC <- NBcol+1
    S4Vectors::metadata(SEprepHCPC)$Results[[1]][[3]] <- append(listHCPC,
                                                                listHCPCplot)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## listHCPCplot=listHCPCplot
    return(SEobj=SEprepHCPC)
}## HCPCanalysis()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

myPaletteHCPC <- function(Nclust=4) {
    ## paletteClust <- c(ggsci::pal_jco("default")(4), "#35B600",
    ##                   ggsci::pal_jco("default")(10)[seq(5, 10, 1)])
    paletteClust <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
                      "#35B600", "#7AA6DCFF", "#003C67FF", "#8F7700FF",
                      "#3B3B3BFF", "#A73030FF", "#4A6990FF")

    if (Nclust > 10) {
        ## MyColclust2 <- scales::hue_pal(direction=-1)(Nclust-8)
        MyColclust2 <- rev(gg_color_hue(n=Nclust-8, l=65))
        paletteClust <- c(paletteClust, MyColclust2)
    }## if(Nclust>10)

    paletteClust <- as.character(paletteClust)[seq_len(Nclust)]

    return(paletteClust)
}## myPaletteHCPC()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrHCPCanalysis <- function(SEresNorm,
                            DATAnorm=TRUE,
                            gene.deletion=NULL,
                            sample.deletion=NULL,
                            Plot.HCPC=FALSE,
                            Color.Group=NULL,
                            Phi=25, Theta=140, epsilon=0.2,
                            Cex.point=0.7, Cex.label=0.7,
                            motion3D=FALSE,
                            path.result=NULL,
                            Name.folder.hcpc=NULL) {
    ##-----------------------------------------------------------------------##
    ## Check deletion
    ErrHCPC1 <- ErrPCAHCPCgraphics(motion3D=motion3D,
                                   Phi=Phi, Theta=Theta, epsilon=epsilon,
                                   Cex.point=Cex.point, Cex.label=Cex.label,
                                   path.result=path.result)

    ##-----------------------------------------------------------------------##
    if (!isTRUE(Plot.HCPC) & !isFALSE(Plot.HCPC)) {
        stop("'Plot.HCPC' must be TRUE or FALSE.")
    }## if (!isTRUE(Plot.HCPC) & !isFALSE(Plot.HCPC))

    if (!is.null(Name.folder.hcpc)) {
        if (!is.character(Name.folder.hcpc)) {
            stop("'Name.folder.hcpc' must be NULL or a character.")
        }## if (!is.character(Name.folder.hcpc))
    }## if (!is.null(Name.folder.hcpc))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrHCPCanalysis()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

HCPCpca3D <- function(PCAcoord3D, matPCAeig, colorINDfactor,
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
    CEXlgd <- Cex.point*0.8

    ##-----------------------------------------------------------------------##
    plot3D::scatter3D(PCAcoord3D[, 1], PCAcoord3D[, 3], PCAcoord3D[, 2],
                      col=colInd, epsilon=epsilon, cex=Cex.point, colvar=NULL,
                      ticktype="detailed", theta=Theta, phi=Phi, pch=20, d=2,
                      main="3D PCA plot with HCPC clusters", bty="b2",
                      xlab=paste0("dim1 (", round(matPCAeig[1, 2], digits=2),
                                  "%)"),
                      ylab=paste0("dim3 (", round(matPCAeig[3, 2], digits=2),
                                  "%)"),
                      zlab=paste0("dim2 (", round(matPCAeig[2, 2], digits=2),
                                  "%)"))

    plot3D::text3D(PCAcoord3D[, 1] + epsilon,
                   PCAcoord3D[, 3] + epsilon,
                   PCAcoord3D[, 2] + epsilon,
                   labels=rownames(PCAcoord3D), cex=Cex.label, col=colInd,
                   font=2, add=TRUE, colkey=FALSE)

    graphics::legend("right", title="Cluster", legend=seq_len(Nclust),
                     col=colorClust, cex=CEXlgd, text.width=CEXlgd/8,
                     pch=20, horiz=FALSE, xpd=TRUE, inset=c(-0.16))## c(-0.015)
    ##-----------------------------------------------------------------------##
    ## gHCPCpca <- grDevices::recordPlot()
    ## return(gHCPCpca)
}## HCPCpca3D()

