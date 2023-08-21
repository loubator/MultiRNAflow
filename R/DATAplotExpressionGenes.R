#' @title Plot expression of a subset of genes.
#'
#' @description The function allows to plot gene expression profiles
#' according to time and/or biological conditions.
#'
#' @details All results are built from the results of our function
#' [DATAnormalization()].
#'
#' @param SEresNorm Results of the function
#' [DATAnormalization()].
#' @param Vector.row.gene Vector of non negative integers indicating
#' the rows of the genes to be plotted.
#' @param DATAnorm \code{TRUE} or \code{FALSE}. \code{TRUE} by default.
#' \code{TRUE} means the function plots gene normalized expression profiles.
#' \code{FALSE} means the function plots gene raw expression profiles.
#' @param Color.Group \code{NULL} or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute
#' a color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used. \code{NULL} by default.
#' @param Plot.Expression \code{TRUE} or \code{FALSE}. \code{TRUE} by default.
#' If \code{TRUE}, the graph will be plotted.
#' Otherwise no graph will be plotted.
#' @param path.result Character or \code{NULL}. \code{NULL} by default.
#' Path to save all results.
#' If \code{path.result} contains a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}" and a sub sub folder,
#' "1-5_ProfileExpression_\code{Name.folder.profile}"
#' all results will be saved in the sub folder
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}/
#' 1-5_ProfileExpression_\code{Name.folder.profile}".
#' Otherwise, a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}" and/or a sub sub folder
#' "1-5_ProfileExpression_\code{Name.folder.profile}"
#' will be created in \code{path.result} and all results will be saved in
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}/
#' 1-5_ProfileExpression_\code{Name.folder.profile}".
#' If NULL, the results will not be saved in a folder. NULL as default.
#' @param Name.folder.profile Character or \code{NULL}. \code{NULL} by default.
#' If \code{Name.folder.profile} is a character, the folder and
#' sub folder names which will contain the PCA graphs will respectively be
#' "1_UnsupervisedAnalysis_\code{Name.folder.profile}" and
#' "1-5_ProfileExpression_\code{Name.folder.profile}".
#' Otherwise, the folder and sub folder names will respectively be
#' "1_UnsupervisedAnalysis" and "1-5_ProfileExpression".
#'
#' @return The function plots for each gene selected with
#' the input \code{Vector.row.gene}
#' * In the case where samples belong to different time points only :
#' the evolution of the expression of each replicate across time and
#' the evolution of the mean and the standard deviation of the expression
#' across time.
#' * In the case where samples belong to different biological conditions only:
#' a violin plot
#' (see [ggplot2::geom_violin()]),
#' and error bars (standard deviation)
#' (see [ggplot2::geom_errorbar()])
#' for each biological condition.
#' * In the case where samples belong to different time points and different
#' biological conditions : the evolution of the expression of each replicate
#' across time and the evolution of the mean and the standard deviation
#' of the expression across time for each biological condition.
#'
#' @seealso The function calls our R function
#' [DATAnormalization()]
#' fisrt, then
#' [DATAplotExpression1Gene()]
#' for each selected genes with \code{Vector.row.gene}.
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom SummarizedExperiment assays colData rownames
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
#' ##------------------------------------------------------------------------#
#' resEVOgenes <- DATAplotExpressionGenes(SEresNorm=resNorm,
#'                                        Vector.row.gene=c(1,3),
#'                                        DATAnorm=TRUE,
#'                                        Color.Group=NULL,
#'                                        Plot.Expression=TRUE,
#'                                        path.result=NULL,
#'                                        Name.folder.profile=NULL)

DATAplotExpressionGenes <- function(SEresNorm,
                                    Vector.row.gene,
                                    DATAnorm=TRUE,
                                    Color.Group=NULL,
                                    Plot.Expression=TRUE,
                                    path.result=NULL,
                                    Name.folder.profile=NULL) {
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check
    ## DATAprepSE
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    if (!is(SEresNorm, "SummarizedExperiment")) {
        stop(Err_SE)
    } else {
        codeDEres <- S4Vectors::metadata(SEresNorm)$SEidentification

        if (is.null(codeDEres)) {
            stop(Err_SE)
        }## if (is.null(codeDEres))

        if (codeDEres != "SEresNormalization") {
            stop(Err_SE)
        }## if (codeDEres != "SEresNormalization")
    }## if (!is(SEresNorm, "SummarizedExperiment"))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check 2
    Err_integers <- paste("'Vector.row.gene' must be a vector",
                          "of non negative integers.")
    if (!is.numeric(Vector.row.gene) & !is.integer(Vector.row.gene)) {
        stop(Err_integers)
    } else {
        if (sum(abs(floor(Vector.row.gene)-Vector.row.gene)) != 0) {
            stop(Err_integers)
        }## if(floor(Individual.position) != Individual.position)

        if (min(Vector.row.gene) <= 0) {
            stop(Err_integers)
        }## if (min(Vector.row.gene) <= 0)
    }## if(is.null(Individual.position))

    if (!isTRUE(DATAnorm) & !isFALSE(DATAnorm)) {
        stop("'DATAnorm' must be TRUE or FALSE.")
    }## if (!isTRUE(DATAnorm) & !isFALSE(DATAnorm))

    if (!is.null(Color.Group)) {
        if (!is.data.frame(Color.Group)) {
            stop("'Color.Group' must be NULL or a data.frame.")
        }## if (!is.data.frame(Color.Group))
    }## if (!is.null(Color.Group))

    if (!isTRUE(Plot.Expression) & !isFALSE(Plot.Expression)) {
        stop("'Plot.Expression' must be TRUE or FALSE.")
    }## if (!isTRUE(Plot.Expression) & !isFALSE(Plot.Expression))

    if (!is.null(path.result)) {
        if (!is.character(path.result)) {
            stop("'path.result' must be NULL or a character.")
        }## if (!is.character(path.result))
    }## if (!is.null(path.result))

    if (!is.null(Name.folder.profile)) {
        if (!is.character(Name.folder.profile)) {
            stop("'Name.folder.profile' must be NULL or a character.")
        }## if (!is.character(Name.folder.profile))
    }## if (!is.null(Name.folder.profile))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Folder creation if no existence
    if (is.null(Name.folder.profile)) {
        Name.folder.profile <- ""
        SubFolder.name <- "1_UnsupervisedAnalysis"
    } else {
        Name.folder.profile <- paste0("_", Name.folder.profile)
        SubFolder.name <- paste0("1_UnsupervisedAnalysis", Name.folder.profile)
    }## if(is.null(Name.folder.profile))

    if (!is.null(path.result)) {
        if(!SubFolder.name%in%dir(path=path.result)){
            print("Folder creation")
            dir.create(path=file.path(path.result, SubFolder.name))
            path.result.f <- file.path(path.result, SubFolder.name)
        } else {
            path.result.f <- file.path(path.result, SubFolder.name)
        }## if (!is.null(path.result))
    } else {
        path.result.f <- NULL
    }## if(!is.null(path.result)=)

    if (!is.null(path.result.f)) {
        nom.dossier.result <- paste0("1-5_ProfileExpressionAnalysis",
                                     Name.folder.profile)
        if (!nom.dossier.result%in%dir(path = path.result.f)) {
            dir.create(path=file.path(path.result.f, nom.dossier.result))
            path.result.new <- file.path(path.result.f, nom.dossier.result)
        } else {
            path.result.new <- file.path(path.result.f, nom.dossier.result)
        }## if(nom.dossier.result%in%dir(path = path.result.f)==FALSE)
    } else {
        path.result.new <- NULL
    }## if(is.null(path.result))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## SUB SE object
    if (DATAnorm == TRUE) {
        aSE <- 2
    } else {
        aSE <- 1
    }## if (DATAnorm == TRUE)

    NameG <- as.character(SummarizedExperiment::rownames(SEresNorm))
    assaySE <- data.frame(SummarizedExperiment::assays(SEresNorm)[[aSE]])
    cSEdat <- SummarizedExperiment::colData(SEresNorm)
    metaSelect <- S4Vectors::metadata(SEresNorm)[c("RAWcolnames", "colGene",
                                                   "colINFOfactors", "formula",
                                                   "SEidentification")]

    subSEnorm <- SEobjFUN(as.matrix(assaySE[Vector.row.gene,]), cSEdat)
    S4Vectors::metadata(subSEnorm) <- metaSelect
    ##------------------------------------------------------------------------#
    List.All.G <- vector(mode="list", length=length(Vector.row.gene))
    names(List.All.G) <- NameG[Vector.row.gene]

    cpt <- 0
    for (g.sel in Vector.row.gene) {
        cpt <- cpt+1
        PlotExpr1G <- DATAplotExpression1Gene(SEres=subSEnorm,
                                              row.gene=cpt,
                                              Color.Group=Color.Group)
        List.All.G[[cpt]] <- PlotExpr1G
    }## for(g.sel in Vector.row.gene)

    S4Vectors::metadata(subSEnorm)$List.plots <- List.All.G

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Save of all graph in a pdf file
    if (!is.null(path.result)) {
        grDevices::pdf(file.path(path.result.new,
                                 paste0("PlotsProfileGeneExpression",
                                        Name.folder.profile, ".pdf")),
                       width=11, height=8, onefile=TRUE)

        for (g.sel in seq_len(length(List.All.G))) {
            print(List.All.G[[g.sel]])
        }## for(g.sel in Vector.row.gene)

        grDevices::dev.off()
    }## if(is.null(path.result)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if (isTRUE(Plot.Expression)) {
        for (g.sel in seq_len(length(List.All.G))) {
            print(List.All.G[[g.sel]])
        }## for(g.sel in Vector.row.gene)
    }## if(isTRUE(Plot.Expression))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(SEobj=subSEnorm)
}## DATAplotExpressionGenes()
