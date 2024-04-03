#' @title Homo sapiens transcript database
#'
#' @description The database is a data.frame which contains transcript length
#' of homo sapiens genes (40452 genes).
#'
#' @format A data frame with 500 rows (genes) and 13 columns (samples).
#' The column names are as follow
#'  \describe{
#'            \item{symbol}{ENSEMBL gene names.}
#'            \item{Median.length.RNA}{The sample is the first replica (r1)
#'            of the biological condition N1wt and T1wt.}
#'  }
#'
#' @details
#' The first column contains genes symbol of the homo sapiens organism and
#' the second column contains the median of transcript length for each gene
#' of the first column.
#'
#' @source {HGNC, ENSEMBL and NCBI database.}
#'
#' @return Mouse dataset with four biological conditions.
#'
#' @usage data(Transcript_HomoSapiens_Database)
#'
#' @examples
#' data(Transcript_HomoSapiens_Database)
"Transcript_HomoSapiens_Database"
