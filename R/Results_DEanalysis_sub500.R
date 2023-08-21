#' @title DE results of three dataset
#'
#' @description The list Results_DEanalysis_sub500 contains the results of
#' [DEanalysisGlobal()]
#' for each of the following raw counts :
#' RawCounts_Weger2021_MOUSEsub500,
#' RawCounts_Leong2014_FISSIONsub500wt and
#' RawCounts_Schleiss2021_CLLsub500
#'
#' @format A list of 3 SummarizedExperiment class object
#'
#' @details
#' Each list in Results_DEanalysis_sub500 contains only the necessary outputs
#' of [DEanalysisGlobal()],
#' needed for the functions:
#' [DEplotVolcanoMA()],
#' [DEplotHeatmaps()],
#' [GSEApreprocessing()],
#' and
#' [GSEAQuickAnalysis()],
#' for each of the following raw counts :
#' RawCounts_Weger2021_MOUSEsub500,
#' RawCounts_Leong2014_FISSIONsub500wt and
#' RawCounts_Schleiss2021_CLLsub500
#'
#' @return \code{Results_DEanalysis_sub500} contains the outputs of
#' [DEanalysisGlobal()] of:
#' RawCounts_Weger2021_MOUSEsub500,
#' RawCounts_Leong2014_FISSIONsub500wt and
#' RawCounts_Schleiss2021_CLLsub500
#'
#' @usage data(Results_DEanalysis_sub500)
#'
#' @examples
#' data(Results_DEanalysis_sub500)
"Results_DEanalysis_sub500"
