#' @title Barplot of DE genes from a contingency table.
#'
#' @description From a contingency table between two variables,
#' the function plots a barplot of the frequency distribution of one variable
#' against the other (see \code{Details}).
#'
#' @details
#' A contingency table (or cross-tabulation) is a table that displays the
#' frequency distribution of two variables (each containing several levels),
#' i.e. the number of observation recorded per pair of levels.
#' The function plots a single barplot from \code{ContingencyTable}.
#'
#' This function is called by [DEanalysisGroup()] and
#' [DEanalysisTimeAndGroup()].
#' These two functions produce several contingency tables,
#' giving information about specific and particular
#' DE genes, as described below.
#'
#' First, we look for all genes that are DE between at least two biological
#' conditions.
#' A gene will be called specific to a given biological condition BC1,
#' if the gene is DE between BC1 and any other biological conditions,
#' but not DE between any pair of other biological conditions.
#' Then each DE gene will be categorized as follow:
#' * If a gene is not specific, the gene will be categorized as 'Other'.
#' The category 'Other' does not exist when there are only two biological
#' conditions.
#' * If a gene is specific to a given biological condition BC1 and expressions
#' in BC1 are higher than in the other biological conditions,
#' the gene will be categorized as 'Upregulated'.
#' * If a gene is specific to a given biological condition BC1 and expressions
#' in BC1 are lower than in the other biological conditions,
#' the gene will be categorized as 'Downregulated'.
#'
#' The functions [DEanalysisGroup()] and [DEanalysisTimeAndGroup()] produce two
#' contingency table that allow to plot both
#' * the number of genes categorized as 'Other', 'Upregulated'
#' and 'Downregulated'
#' (only when there are strictly more than two biological conditions).
#' * the number of genes categorized 'Upregulated' and 'Downregulated'.
#'
#' Second, we look for all genes that are DE between at least one time point
#' (except t0) and t0 for each biological condition.
#' A gene will be categorized as 'particular' to a given biological condition
#' BC1 for a given time point ti (except t0), if the gene is DE between ti and
#' t0 for the biological condition BC1, but not DE between ti and t0 for
#' the other biological conditions.
#' A gene will be categorized as 'common' to all biological conditions,
#' if the gene is DE between ti and t0 for all biological conditions.
#' Otherwise, a gene will categorized as 'Other'.
#'
#' The function [DEanalysisTimeAndGroup()] produces a contingency table that
#' allow to plot the number of 'specific', 'common' and 'other' genes for
#' each ti (except t0).
#'
#' @param ContingencyTable A numeric data.frame,
#' corresponding to a contingency table, of dimension N1*N2, with N1 and N2,
#' respectively the number of levels in the first and second variable
#' (see examples and details).
#' @param dodge \code{TRUE} or \code{FALSE}.
#' \code{FALSE} means multiple bars in the barplot
#' (one per level of the first variable)
#' one for each fixed level of the other variable.
#' \code{TRUE} means multiple bars will be dodged side-to-side
#' (see [ggplot2::geom_bar()]).
#'
#' @return A barplot using [ggplot2] (see details).
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_bar position_dodge xlab ylab guide_axis
#' scale_x_discrete scale_fill_manual
#'
#' @seealso The [DEplotBarplot()] function
#' * is used by the following functions of our package: [DEanalysisGroup()]
#' and [DEanalysisTimeAndGroup()].
#' * calls the R package [ggplot2] in order to plot the barplot.
#'
#' @export
#'
#' @examples
#' ## Data simulation
#' CrossTabulation <- matrix(c(75,30,10,5, 5,35,5,20, 220,235,285,275),
#'                           ncol=4, byrow=TRUE)
#' colnames(CrossTabulation) <- c("A", "B", "C", "D")
#' row.names(CrossTabulation) <- c("Spe.Pos", "Spe.Neg", "Other")
#'
#' ##------------------------------------------------------------------------##
#' res.dodgeTRUE <- DEplotBarplot(ContingencyTable=CrossTabulation,dodge=FALSE)
#' res.dodgeTRUE
#'
#' res.dodgeFALSE <- DEplotBarplot(ContingencyTable=CrossTabulation,dodge=TRUE)
#' res.dodgeFALSE

DEplotBarplot <- function(ContingencyTable,
                          dodge=TRUE) {
    ##-----------------------------------------------------------------------##
    ## To avoid "no visible binding for global variable" with devtools::check()
    variable <- value <- Attribute <- NULL

    ##-----------------------------------------------------------------------##
    ## Data preprocessing for graph
    barplotData_ini <- data.frame(Attribute=row.names(ContingencyTable),
                                  ContingencyTable)
    barplotData <- reshape2::melt(barplotData_ini, id.vars=c("Attribute"))

    ColLegend <- barplotData$Attribute
    levels_f <- c(levels(factor(ColLegend))[-1], levels(factor(ColLegend))[1])

    barplotData$Attribute <- factor(barplotData$Attribute, levels=levels_f)

    ##-----------------------------------------------------------------------##
    ## Graph
    q.contGT <- ggplot2::ggplot(data=barplotData,
                                ggplot2::aes(x=variable, y=value,
                                             fill=Attribute)) +
        ggplot2::ylab("Number of Genes") + ggplot2::xlab("") +
        ggplot2::scale_x_discrete(guide=ggplot2::guide_axis(angle=45))

    if (isTRUE(dodge)) {
        q.contGT <- q.contGT +
            ggplot2::geom_bar(stat="identity", color="black",
                              position=ggplot2::position_dodge())
    } else {
        q.contGT <- q.contGT +
            ggplot2::geom_bar(stat="identity", color="black")
    }## if(isTRUE(dodge))

    ##-----------------------------------------------------------------------##
    if (length(unique(barplotData$Attribute)) == 3) {
        q.contGT <- q.contGT +
            ggplot2::scale_fill_manual(values=c("#999999", "#E41A1C",
                                                "steelblue"))
    }## if(length(unique(barplotData$Attribute))==3)

    if (length(unique(barplotData$Attribute)) == 2) {
        q.contGT <- q.contGT +
            ggplot2::scale_fill_manual(values=c("#E41A1C", "steelblue"))
    }## if(length(unique(barplotData$Attribute))==2)

    ##-----------------------------------------------------------------------##
    ## Output
    return(graph.cont.g.t=q.contGT)
}## DEplotBarplot()
