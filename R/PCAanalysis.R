#' @title Automatic PCA analysis (Main function)
#'
#' @description The functions performs an automatic principal component
#' analysis (PCA) from a gene expression dataset where samples can belong to
#' different biological conditions and/or time points.
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
#' @param path.result Character or \code{NULL}. Path to save all results.
#' If \code{path.result} contains a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}" and a sub sub folder,
#' "1-2_PCAanalysis_\code{Name.folder.pca}"
#' all results will be saved in the sub folder
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}/
#' 1-2_PCAanalysis_\code{Name.folder.pca}".
#' Otherwise, a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}" and/or a sub sub folder
#' "1-2_PCAanalysis_\code{Name.folder.pca}"
#' will be created in \code{path.result} and all results will be saved in
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}/
#' 1-2_PCAanalysis_\code{Name.folder.pca}".
#' If \code{NULL}, the results will not be saved in a folder.
#' \code{NULL} as default.
#' @param Name.folder.pca Character or \code{NULL}.
#' If \code{Name.folder.pca} is a character, the folder and sub folder names
#' which will contain the PCA graphs will respectively be
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}"
#' and "1-2_PCAanalysis_\code{Name.folder.pca}".
#' Otherwise, the folder and sub folder names will respectively be
#' "1_UnsupervisedAnalysis" and "1-2_PCAanalysis".
#'
#' @return The function returns the same SummarizedExperiment class object
#' \code{SEresNorm} with the outputs from the function
#' [FactoMineR::PCA()],
#' and several 2D and 3D PCA graphs depending on the experimental design
#' (if \code{Plot.PCA=TRUE}),
#' saved in the metadata \code{Results[[1]][[2]]} of \code{SEresNorm},
#' * When samples belong only to different biological conditions,
#' the function returns a 2D and two 3D PCA graphs.
#' In each graph, samples are colored with different colors for different
#' biological conditions. The two 3D PCA graphs are identical but one of them
#' will be opened in a rgl window
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
#' * When samples belong to different time points and
#' different biological conditions, the function returns
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
#' @seealso The function calls the R functions
#' [PCAgraphics()] and
#' [ColnamesToFactors()].
#'
#' @importFrom stats var as.formula
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors metadata
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
#' resPCAanalysis <- PCAanalysis(SEresNorm=resNorm,
#'                               DATAnorm=TRUE,
#'                               gene.deletion=c("Gene1", "Gene5"),
#'                               sample.deletion=c(2, 6),
#'                               Plot.PCA=TRUE,
#'                               Mean.Accross.Time=FALSE,
#'                               Color.Group=GROUPcolor,
#'                               motion3D=FALSE,
#'                               Phi=25, Theta=140,
#'                               Cex.label=0.7, Cex.point=0.7, epsilon=0.2,
#'                               path.result=NULL, Name.folder.pca=NULL)

PCAanalysis <- function(SEresNorm,
                        DATAnorm=TRUE,
                        gene.deletion=NULL,
                        sample.deletion=NULL,
                        Plot.PCA=TRUE,
                        Mean.Accross.Time=FALSE,
                        Color.Group=NULL,
                        Phi=25, Theta=140, epsilon=0.2,
                        Cex.point=0.7, Cex.label=0.7,
                        motion3D=FALSE,
                        path.result=NULL,
                        Name.folder.pca=NULL) {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Check inputs
    resErr <- ErrPCAanalysis(SEresNorm=SEresNorm, DATAnorm=DATAnorm,
                             gene.deletion=gene.deletion,
                             sample.deletion=sample.deletion,
                             Plot.PCA=Plot.PCA, motion3D=motion3D,
                             Mean.Accross.Time=Mean.Accross.Time,
                             Phi=Phi, Theta=Theta, epsilon=epsilon,
                             Cex.point=Cex.point, Cex.label=Cex.label,
                             Color.Group=Color.Group, path.result=path.result,
                             Name.folder.pca=Name.folder.pca)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Preprocessing
    if (isTRUE(DATAnorm)) {
        aSE <- 2
    } else {
        aSE <- 1
    }## if (isTRUE(DATAnorm))

    cSEdat <- SummarizedExperiment::colData(SEresNorm)

    if (c("Group")%in%colnames(cSEdat)) {
        Vector.group <- as.character(cSEdat$Group)
    } else {
        Vector.group <- NULL
    }## if (c("Group")%in%colnames(cSEdat))

    if (c("Time")%in%colnames(cSEdat)) {
        Vector.time <- as.character(cSEdat$Time)
    } else {
        Vector.time <- NULL
    }## if (c("Time")%in%colnames(cSEdat))

    tb.spinfoini <- as.numeric(table(cSEdat$ID))
    Var.sample <- stats::var(tb.spinfoini)
    max.tb <- max(tb.spinfoini)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Folder creation if no existence
    if (is.null(Name.folder.pca)) {
        Name.folder.pca <- ""
        SubFolder.name <- "1_UnsupervisedAnalysis"
    } else {
        Name.folder.pca <- paste0("_", Name.folder.pca)
        SubFolder.name <- paste0("1_UnsupervisedAnalysis", Name.folder.pca)
    }## if (is.null(Name.folder.pca))

    SubFolderPCA <- paste0("1-2_PCAanalysis", Name.folder.pca)

    if (!is.null(path.result)) {
        if (!SubFolder.name%in%dir(path=path.result)) {
            print("Folder creation")
            dir.create(path=file.path(path.result, SubFolder.name))
        }## if(!SubFolder.name%in%dir(path=path.result))
        path.result.f <- file.path(path.result, SubFolder.name)
    } else {
        path.result.f <- NULL
    }## if (!is.null(path.result))

    ##-----------------------------------------------------------------------##
    if (!is.null(Vector.group)) {
        if (!is.null(Vector.time)) {
            Name.file.pca <- paste0("allBiologicalConditions", Name.folder.pca)
        } else {
            Name.file.pca <- paste0("BiologicalCondition", Name.folder.pca)
        }## if(!is.null(Vector.time))
    } else {
        Name.file.pca <- paste0("Time", Name.folder.pca)
    }## if(!is.null(Vector.group))

    if (!is.null(path.result.f)) {
        if (!SubFolderPCA%in%dir(path=path.result.f)) {
            dir.create(path=file.path(path.result.f, SubFolderPCA))
        }## if(!SubFolderPCA%in%dir(path=path.result.f))
        path.result.new <- file.path(path.result.f, SubFolderPCA)
    } else {
        path.result.new <- NULL
    }## if(!is.null(path.result.f))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Main results ##tb.spinfoini[1]>1
    if (isFALSE(Mean.Accross.Time) & Var.sample==0 & max.tb>1) {
        SEresPCA <- PCAgraphics(SEresNorm=SEresNorm,
                                DATAnorm=DATAnorm,
                                sample.deletion=sample.deletion,
                                Plot.PCA=Plot.PCA,
                                Mean.Accross.Time=FALSE,
                                gene.deletion=gene.deletion,
                                Color.Group=Color.Group,
                                motion3D=motion3D,
                                Phi=Phi, Theta=Theta, epsilon=epsilon,
                                Cex.point=Cex.point, Cex.label=Cex.label,
                                path.result=path.result.new,
                                Name.file.pca=Name.file.pca)
    } else {
        SEresPCA <- PCAgraphics(SEresNorm=SEresNorm,
                                DATAnorm=DATAnorm,
                                sample.deletion=sample.deletion,
                                gene.deletion=gene.deletion,
                                Plot.PCA=Plot.PCA,
                                Mean.Accross.Time=TRUE,
                                Color.Group=Color.Group,
                                motion3D=motion3D,
                                Phi=Phi, Theta=Theta, epsilon=epsilon,
                                Cex.point=Cex.point, Cex.label=Cex.label,
                                path.result=path.result.new,
                                Name.file.pca=Name.file.pca)
    }## if(Mean.Accross.Time==FALSE & Var.sample==0 & max.tb>1)

    PCAlist <- S4Vectors::metadata(SEresPCA)$Results[[1]][[2]]

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    if (!is.null(Vector.group) & !is.null(Vector.time)) {
        ##-------------------------------------------------------------------##
        Data1 <- data.frame(SummarizedExperiment::assays(SEresNorm)[[1]])
        Data2 <- data.frame(SummarizedExperiment::assays(SEresNorm)[[2]])
        ## NameG <- as.character(SummarizedExperiment::rownames(SEresNorm))

        Tt.Del <- gsub("t", "", gsub("T", "", as.character(Vector.time)))
        Time.name <- levels(as.factor(paste0("T", Tt.Del)))

        Group.Levels <- levels(as.factor(Vector.group))

        if (!is.null(sample.deletion)) {
            colG <- S4Vectors::metadata(SEresNorm)$colGene
            RAWcolN <- S4Vectors::metadata(SEresNorm)$RAWcolnames
            delIDsample <- delFUNsample(sample.deletion=sample.deletion,
                                        RAWcolnames=RAWcolN, Column.gene=colG)
        }## if (is.null(sample.deletion))

        ##-------------------------------------------------------------------##
        ##-------------------------------------------------------------------##
        for (g in seq_len(length(Group.Levels))) {
            Index.g <- which(Vector.group == Group.Levels[g])

            if (is.null(sample.deletion)) {
                sample.deletion_g <- NULL
            } else {
                interIDs_gDEL <- intersect(sort(delIDsample), Index.g)
                if (length(interIDs_gDEL) == 0) {
                    sample.deletion_g <- NULL
                } else {
                    sample.deletion_g <- interIDs_gDEL
                }## if (length(interIDs_gDEL) == 0)
            }## if (is.null(sample.deletion))

            Data1g <- as.matrix(Data1[, Index.g])
            Data2g <- as.matrix(Data2[, Index.g])
            cSEdatg <- cSEdat[Index.g, -1]
            name2g <- names(SummarizedExperiment::assays(SEresNorm))[2]
            Meta1 <- S4Vectors::metadata(SEresNorm)
            formulag <- stats::as.formula(counts ~ Time)
            SEiden <- "SEresNormalization"

            SEresNormg <- SEobjFUN(Data1g, cSEdatg)
            SummarizedExperiment::assays(SEresNormg)[[2]] <- Data2g
            names(SummarizedExperiment::assays(SEresNormg))[2] <- name2g
            ## SummarizedExperiment::rownames(SEresNormg) <- NameG
            S4Vectors::metadata(SEresNormg) <- Meta1
            S4Vectors::metadata(SEresNormg)$formula <- formulag
            S4Vectors::metadata(SEresNormg)$SEidentification <- SEiden

            Name.file.pca.g <- paste0("BiologicalCondition_", Group.Levels[g])

            # if (isTRUE(Plot.PCA)) {
            #     graphics::plot.new() ## clean up device
            # }## if (isTRUE(Plot.PCA))

            ##---------------------------------------------------------------##
            ##---------------------------------------------------------------##
            if(isFALSE(Mean.Accross.Time) & Var.sample==0 & tb.spinfoini[1]>1){
                SEresPCAg <- PCAgraphics(SEresNorm=SEresNormg,
                                         DATAnorm=DATAnorm,
                                         sample.deletion=sample.deletion_g,
                                         gene.deletion=gene.deletion,
                                         Plot.PCA=Plot.PCA,
                                         Mean.Accross.Time=FALSE,
                                         Color.Group=Color.Group,
                                         motion3D=motion3D,
                                         Phi=Phi, Theta=Theta, epsilon=epsilon,
                                         Cex.point=Cex.point,
                                         Cex.label=Cex.label,
                                         path.result=path.result.new,
                                         Name.file.pca=Name.file.pca.g)
            } else {
                SEresPCAg <- PCAgraphics(SEresNorm=SEresNormg,
                                         DATAnorm=DATAnorm,
                                         sample.deletion=sample.deletion_g,
                                         gene.deletion=gene.deletion,
                                         Plot.PCA=Plot.PCA,
                                         Mean.Accross.Time=TRUE,
                                         Color.Group=Color.Group,
                                         motion3D=motion3D,
                                         Phi=Phi, Theta=Theta, epsilon=epsilon,
                                         Cex.point=Cex.point,
                                         Cex.label=Cex.label,
                                         path.result=path.result.new,
                                         Name.file.pca=Name.file.pca.g)
            }## if(isFALSE(Mean.Accross.Time)&Var.sample==0&tb.spinfoini[1]>1)

            PCAlist_BC <- S4Vectors::metadata(SEresPCAg)$Results[[1]][[2]]
            nPCAlist <- length(PCAlist)

            PCAlist <- append(PCAlist, list(PCAlist_BC[-c(1, 2, 3)]))
            names(PCAlist)[nPCAlist+1] <- paste0("PCA_BiologicalCondition_",
                                                 Group.Levels[g])
        }## for (g in seq_len(length(Group.Levels)))

    }## if (!is.null(Vector.group) & !is.null(Vector.time))

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## SE object
    SEresPCAall <- SEresPCA
    S4Vectors::metadata(SEresPCAall)$Results[[1]][[2]] <- PCAlist[-1]

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    ## Output
    return(SEobj=SEresPCAall)
}## PCAanalysis()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

delFUNsample <- function(sample.deletion=NULL,
                         RAWcolnames="sample",
                         Column.gene=NULL) {
    ##-----------------------------------------------------------------------##
    if (is.numeric(sample.deletion)) {
        if (is.null(Column.gene)) {
            delIDsample <- sample.deletion
        } else {
            ColSplDel <- RAWcolnames[sample.deletion]
            delIDsample <- which(RAWcolnames[-Column.gene]%in%ColSplDel)
            delIDsample <- delIDsample + 1
        }## if(is.null(Column.gene))
    } else {
        if (is.null(Column.gene)) {
            delIDsample <- which(RAWcolnames%in%sample.deletion)
        } else {
            delIDsample <- which(RAWcolnames[-Column.gene]%in%sample.deletion)
            delIDsample <- delIDsample + 1
        }## if(is.null(Column.gene))
    }## if(is.numeric(sample.deletion))

    ##-----------------------------------------------------------------------##
    return(sort(delIDsample))
}## delFUNsample()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

ErrPCAanalysis <- function(SEresNorm,
                           DATAnorm=TRUE,
                           gene.deletion=NULL,
                           sample.deletion=NULL,
                           Plot.PCA=TRUE,
                           Mean.Accross.Time=FALSE,
                           Color.Group=NULL,
                           Phi=25, Theta=140, epsilon=0.2,
                           Cex.point=0.7, Cex.label=0.7,
                           motion3D=FALSE,
                           path.result=NULL,
                           Name.folder.pca=NULL) {
    ##-----------------------------------------------------------------------##
    ## Check ErrPCArealization
    resErr1 <- ErrPCArealization(SEresNorm=SEresNorm,
                                 DATAnorm=DATAnorm,
                                 gene.deletion=gene.deletion,
                                 sample.deletion=sample.deletion,
                                 Supp.del.sample=FALSE)

    ##-----------------------------------------------------------------------##
    ## Check ErrPCAgraphics
    resErr2 <- ErrPCAgraphics(Plot.PCA=Plot.PCA,
                              Mean.Accross.Time=Mean.Accross.Time,
                              motion3D=motion3D,
                              Phi=Phi, Theta=Theta, epsilon=epsilon,
                              Cex.point=Cex.point, Cex.label=Cex.label,
                              path.result=path.result)

    ##-----------------------------------------------------------------------##
    ## Check Name.folder.pca
    if (!is.null(Name.folder.pca)) {
        if (!is.character(Name.folder.pca)) {
            stop("'Name.folder.pca' must be NULL or a character.")
        }## if (!is.character(Name.folder.pca))
    }## if (!is.null(Name.folder.pca))

    ##-----------------------------------------------------------------------##
    return(Message="No error")
}## ErrPCAanalysis()
