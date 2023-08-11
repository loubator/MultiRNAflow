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
#' @param Vector.row.gene Vector of integer indicating the rows of the genes
#' to be plotted.
#' @param DATAnorm \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
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
#' \code{Color.Group} will not be used.
#' @param Plot.Expression \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, the graph will be plotted.
#' Otherwise no graph will be plotted.
#' @param path.result Character or \code{NULL}. Path to save all results.
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
#' @param Name.folder.profile Character or \code{NULL}.
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
    if (is.null(SEresNorm$SEidentification)) {
        stop(Err_SE)
    } else {
        if (SEresNorm$SEidentification != "SEresNormalization") {
            stop(Err_SE)
        }## if (SEresNorm$SEidentification != "SEresNormalization")
    }## if ((is.null(SEresNorm$SEidentification))

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

    assaySE <- data.frame(SummarizedExperiment::assays(SEresNorm$SEobj)[[aSE]])
    cSEdat <- SummarizedExperiment::colData(SEresNorm$SEobj)
    subSEnorm <- SEobjFUN(as.matrix(assaySE[Vector.row.gene,]), cSEdat)

    NameG <- as.character(SummarizedExperiment::rownames(subSEnorm))

    subSEresNorm <- SEresNorm
    subSEresNorm$SEobj <-subSEnorm

    ##------------------------------------------------------------------------#
    List.All.G <- vector(mode="list", length=length(Vector.row.gene))
    names(List.All.G) <- NameG

    cpt <- 0
    for (g.sel in Vector.row.gene) {
        cpt <- cpt+1
        PlotExpr1G <- DATAplotExpression1Gene(SEres=subSEresNorm,
                                              row.gene=cpt,
                                              Color.Group=Color.Group)
        List.All.G[[cpt]] <- PlotExpr1G
    }## for(g.sel in Vector.row.gene)

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
    return(list(SEobj=subSEnorm,
                List.plots=List.All.G))
}## DATAplotExpressionGenes()
