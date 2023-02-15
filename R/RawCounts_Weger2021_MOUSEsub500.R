#' @title Mouse count data with four biological conditions,
#' six time measurements and 500 genes.
#'
#' @description This time series count data (read counts) represents
#' the temporal transcriptional response
#' (six time measurements across the course of a day) of Bmal1 wild-type (WT)
#' and Cry1/2 WT, Bmal1 KO and Cry1/2 WT, Bmal1 (WT) and Cry1/2 KO,
#' and Bmal1 KO and Cry1/2 KO mice under an ad libitum (AL) or
#' night restricted feeding (RF) regimen.
#' Therefore, there are eight biological conditions.
#' As there are only two mice per biological condition,
#' we will not consider the effect of the regimen.
#' The original dataset has 40327 genes but we kept only 500 genes
#' in order to increase the speed of each function in our algorithm.
#'
#' @format A data frame with 500 rows (genes) and 97 columns (samples).
#' The column names are as follow
#'  \describe{
#'            \item{Gene}{ENSEMBL gene names.}
#'            \item{BmKo_t0_r1}{The sample is the first replica (r1) of the biological condition Bmal1 and KO at time t0 (00h).}
#'            \item{BmKo_t1_r1}{The sample is the first replica (r1) of the biological condition Bmal1 and KO at time t1 (04h).}
#'            \item{BmKo_t2_r1}{The sample is the first replica (r1) of the biological condition Bmal1 and KO at time t2 (08h).}
#'            \item{BmKo_t3_r1}{The sample is the first replica (r1) of the biological condition Bmal1 and KO at time t3 (12h).}
#'            \item{BmKo_t4_r1}{The sample is the first replica (r1) of the biological condition Bmal1 and KO at time t4 (16h).}
#'            \item{BmKo_t5_r1}{The sample is the first replica (r1) of the biological condition Bmal1 and KO at time t5 (20h).}
#'            \item{BmKo_t0_r2}{The sample is the first replica (r2) of the biological condition Bmal1 and KO at time t0 (00h).}
#'            \item{BmKo_t1_r2}{The sample is the first replica (r2) of the biological condition Bmal1 and KO at time t1 (04h).}
#'            \item{BmKo_t2_r2}{The sample is the first replica (r2) of the biological condition Bmal1 and KO at time t2 (08h).}
#'            \item{BmKo_t3_r2}{The sample is the first replica (r2) of the biological condition Bmal1 and KO at time t3 (12h).}
#'            \item{BmKo_t4_r2}{The sample is the first replica (r2) of the biological condition Bmal1 and KO at time t4 (16h).}
#'            \item{BmKo_t5_r2}{The sample is the first replica (r2) of the biological condition Bmal1 and KO at time t5 (20h).}
#'            \item{BmKo_t0_r3}{The sample is the first replica (r3) of the biological condition Bmal1 and KO at time t0 (00h).}
#'            \item{BmKo_t1_r3}{The sample is the first replica (r3) of the biological condition Bmal1 and KO at time t1 (04h).}
#'            \item{BmKo_t2_r3}{The sample is the first replica (r3) of the biological condition Bmal1 and KO at time t2 (08h).}
#'            \item{BmKo_t3_r3}{The sample is the first replica (r3) of the biological condition Bmal1 and KO at time t3 (12h).}
#'            \item{BmKo_t4_r3}{The sample is the first replica (r3) of the biological condition Bmal1 and KO at time t4 (16h).}
#'            \item{BmKo_t5_r3}{The sample is the first replica (r3) of the biological condition Bmal1 and KO at time t5 (20h).}
#'            \item{BmKo_t0_r4}{The sample is the first replica (r4) of the biological condition Bmal1 and KO at time t0 (00h).}
#'            \item{BmKo_t1_r4}{The sample is the first replica (r4) of the biological condition Bmal1 and KO at time t1 (04h).}
#'            \item{BmKo_t2_r4}{The sample is the first replica (r4) of the biological condition Bmal1 and KO at time t2 (08h).}
#'            \item{BmKo_t3_r4}{The sample is the first replica (r4) of the biological condition Bmal1 and KO at time t3 (12h).}
#'            \item{BmKo_t4_r4}{The sample is the first replica (r4) of the biological condition Bmal1 and KO at time t4 (16h).}
#'            \item{BmKo_t5_r4}{The sample is the first replica (r4) of the biological condition Bmal1 and KO at time t5 (20h).}
#'            \item{BmWt_t0_r5}{The sample is the first replica (r5) of the biological condition Bmal1 and wild-type at time t0 (00h).}
#'            \item{BmWt_t1_r5}{The sample is the first replica (r5) of the biological condition Bmal1 and wild-type at time t1 (04h).}
#'            \item{BmWt_t2_r5}{The sample is the first replica (r5) of the biological condition Bmal1 and wild-type at time t2 (08h).}
#'            \item{BmWt_t3_r5}{The sample is the first replica (r5) of the biological condition Bmal1 and wild-type at time t3 (12h).}
#'            \item{BmWt_t4_r5}{The sample is the first replica (r5) of the biological condition Bmal1 and wild-type at time t4 (16h).}
#'            \item{BmWt_t5_r5}{The sample is the first replica (r5) of the biological condition Bmal1 and wild-type at time t5 (20h).}
#'            \item{BmWt_t0_r6}{The sample is the first replica (r6) of the biological condition Bmal1 and wild-type at time t0 (00h).}
#'            \item{BmWt_t1_r6}{The sample is the first replica (r6) of the biological condition Bmal1 and wild-type at time t1 (04h).}
#'            \item{BmWt_t2_r6}{The sample is the first replica (r6) of the biological condition Bmal1 and wild-type at time t2 (08h).}
#'            \item{BmWt_t3_r6}{The sample is the first replica (r6) of the biological condition Bmal1 and wild-type at time t3 (12h).}
#'            \item{BmWt_t4_r6}{The sample is the first replica (r6) of the biological condition Bmal1 and wild-type at time t4 (16h).}
#'            \item{BmWt_t5_r6}{The sample is the first replica (r6) of the biological condition Bmal1 and wild-type at time t5 (20h).}
#'            \item{BmWt_t0_r7}{The sample is the first replica (r7) of the biological condition Bmal1 and wild-type at time t0 (00h).}
#'            \item{BmWt_t1_r7}{The sample is the first replica (r7) of the biological condition Bmal1 and wild-type at time t1 (04h).}
#'            \item{BmWt_t2_r7}{The sample is the first replica (r7) of the biological condition Bmal1 and wild-type at time t2 (08h).}
#'            \item{BmWt_t3_r7}{The sample is the first replica (r7) of the biological condition Bmal1 and wild-type at time t3 (12h).}
#'            \item{BmWt_t4_r7}{The sample is the first replica (r7) of the biological condition Bmal1 and wild-type at time t4 (16h).}
#'            \item{BmWt_t5_r7}{The sample is the first replica (r7) of the biological condition Bmal1 and wild-type at time t5 (20h).}
#'            \item{BmWt_t0_r8}{The sample is the first replica (r8) of the biological condition Bmal1 and wild-type at time t0 (00h).}
#'            \item{BmWt_t1_r8}{The sample is the first replica (r8) of the biological condition Bmal1 and wild-type at time t1 (04h).}
#'            \item{BmWt_t2_r8}{The sample is the first replica (r8) of the biological condition Bmal1 and wild-type at time t2 (08h).}
#'            \item{BmWt_t3_r8}{The sample is the first replica (r8) of the biological condition Bmal1 and wild-type at time t3 (12h).}
#'            \item{BmWt_t4_r8}{The sample is the first replica (r8) of the biological condition Bmal1 and wild-type at time t4 (16h).}
#'            \item{BmWt_t5_r8}{The sample is the first replica (r8) of the biological condition Bmal1 and wild-type at time t5 (20h).}
#'            \item{CrKo_t0_r9}{The sample is the first replica (r9) of the biological condition Cry1/2 and KO at time t0 (00h).}
#'            \item{CrKo_t1_r9}{The sample is the first replica (r9) of the biological condition Cry1/2 and KO at time t1 (04h).}
#'            \item{CrKo_t2_r9}{The sample is the first replica (r9) of the biological condition Cry1/2 and KO at time t2 (08h).}
#'            \item{CrKo_t3_r9}{The sample is the first replica (r9) of the biological condition Cry1/2 and KO at time t3 (12h).}
#'            \item{CrKo_t4_r9}{The sample is the first replica (r9) of the biological condition Cry1/2 and KO at time t4 (16h).}
#'            \item{CrKo_t5_r9}{The sample is the first replica (r9) of the biological condition Cry1/2 and KO at time t5 (20h).}
#'            \item{CrKo_t0_r10}{The sample is the first replica (r10) of the biological condition Cry1/2 and KO at time t0 (00h).}
#'            \item{CrKo_t1_r10}{The sample is the first replica (r10) of the biological condition Cry1/2 and KO at time t1 (04h).}
#'            \item{CrKo_t2_r10}{The sample is the first replica (r10) of the biological condition Cry1/2 and KO at time t2 (08h).}
#'            \item{CrKo_t3_r10}{The sample is the first replica (r10) of the biological condition Cry1/2 and KO at time t3 (12h).}
#'            \item{CrKo_t4_r10}{The sample is the first replica (r10) of the biological condition Cry1/2 and KO at time t4 (16h).}
#'            \item{CrKo_t5_r10}{The sample is the first replica (r10) of the biological condition Cry1/2 and KO at time t5 (20h).}
#'            \item{CrKo_t0_r11}{The sample is the first replica (r11) of the biological condition Cry1/2 and KO at time t0 (00h).}
#'            \item{CrKo_t1_r11}{The sample is the first replica (r11) of the biological condition Cry1/2 and KO at time t1 (04h).}
#'            \item{CrKo_t2_r11}{The sample is the first replica (r11) of the biological condition Cry1/2 and KO at time t2 (08h).}
#'            \item{CrKo_t3_r11}{The sample is the first replica (r11) of the biological condition Cry1/2 and KO at time t3 (12h).}
#'            \item{CrKo_t4_r11}{The sample is the first replica (r11) of the biological condition Cry1/2 and KO at time t4 (16h).}
#'            \item{CrKo_t5_r11}{The sample is the first replica (r11) of the biological condition Cry1/2 and KO at time t5 (20h).}
#'            \item{CrKo_t0_r12}{The sample is the first replica (r12) of the biological condition Cry1/2 and KO at time t0 (00h).}
#'            \item{CrKo_t1_r12}{The sample is the first replica (r12) of the biological condition Cry1/2 and KO at time t1 (04h).}
#'            \item{CrKo_t2_r12}{The sample is the first replica (r12) of the biological condition Cry1/2 and KO at time t2 (08h).}
#'            \item{CrKo_t3_r12}{The sample is the first replica (r12) of the biological condition Cry1/2 and KO at time t3 (12h).}
#'            \item{CrKo_t4_r12}{The sample is the first replica (r12) of the biological condition Cry1/2 and KO at time t4 (16h).}
#'            \item{CrKo_t5_r12}{The sample is the first replica (r12) of the biological condition Cry1/2 and KO at time t5 (20h).}
#'            \item{CrWt_t0_r13}{The sample is the first replica (r13) of the biological condition Cry1/2 and wild-type at time t0 (00h).}
#'            \item{CrWt_t1_r13}{The sample is the first replica (r13) of the biological condition Cry1/2 and wild-type at time t1 (04h).}
#'            \item{CrWt_t2_r13}{The sample is the first replica (r13) of the biological condition Cry1/2 and wild-type at time t2 (08h).}
#'            \item{CrWt_t3_r13}{The sample is the first replica (r13) of the biological condition Cry1/2 and wild-type at time t3 (12h).}
#'            \item{CrWt_t4_r13}{The sample is the first replica (r13) of the biological condition Cry1/2 and wild-type at time t4 (16h).}
#'            \item{CrWt_t5_r13}{The sample is the first replica (r13) of the biological condition Cry1/2 and wild-type at time t5 (20h).}
#'            \item{CrWt_t0_r14}{The sample is the first replica (r14) of the biological condition Cry1/2 and wild-type at time t0 (00h).}
#'            \item{CrWt_t1_r14}{The sample is the first replica (r14) of the biological condition Cry1/2 and wild-type at time t1 (04h).}
#'            \item{CrWt_t2_r14}{The sample is the first replica (r14) of the biological condition Cry1/2 and wild-type at time t2 (08h).}
#'            \item{CrWt_t3_r14}{The sample is the first replica (r14) of the biological condition Cry1/2 and wild-type at time t3 (12h).}
#'            \item{CrWt_t4_r14}{The sample is the first replica (r14) of the biological condition Cry1/2 and wild-type at time t4 (16h).}
#'            \item{CrWt_t5_r14}{The sample is the first replica (r14) of the biological condition Cry1/2 and wild-type at time t5 (20h).}
#'            \item{CrWt_t0_r15}{The sample is the first replica (r15) of the biological condition Cry1/2 and wild-type at time t0 (00h).}
#'            \item{CrWt_t1_r15}{The sample is the first replica (r15) of the biological condition Cry1/2 and wild-type at time t1 (04h).}
#'            \item{CrWt_t2_r15}{The sample is the first replica (r15) of the biological condition Cry1/2 and wild-type at time t2 (08h).}
#'            \item{CrWt_t3_r15}{The sample is the first replica (r15) of the biological condition Cry1/2 and wild-type at time t3 (12h).}
#'            \item{CrWt_t4_r15}{The sample is the first replica (r15) of the biological condition Cry1/2 and wild-type at time t4 (16h).}
#'            \item{CrWt_t5_r15}{The sample is the first replica (r15) of the biological condition Cry1/2 and wild-type at time t5 (20h).}
#'            \item{CrWt_t0_r16}{The sample is the first replica (r16) of the biological condition Cry1/2 and wild-type at time t0 (00h).}
#'            \item{CrWt_t1_r16}{The sample is the first replica (r16) of the biological condition Cry1/2 and wild-type at time t1 (04h).}
#'            \item{CrWt_t2_r16}{The sample is the first replica (r16) of the biological condition Cry1/2 and wild-type at time t2 (08h).}
#'            \item{CrWt_t3_r16}{The sample is the first replica (r16) of the biological condition Cry1/2 and wild-type at time t3 (12h).}
#'            \item{CrWt_t4_r16}{The sample is the first replica (r16) of the biological condition Cry1/2 and wild-type at time t4 (16h).}
#'            \item{CrWt_t5_r16}{The sample is the first replica (r16) of the biological condition Cry1/2 and wild-type at time t5 (20h).}
#'  }
#'
#' @details The data is used in order to describe our algorithm in the case
#' where samples belong to different time points.
#'
#' We kept 500 genes only in order to increase the speed for each example.
#'
#' @source {This dataset comes from Gene Expression Omnibus (GEO)
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135898}.
#' The name of the samples was renamed in order to be used with our package.}
#'
#' @usage data(RawCounts_Weger2021_MOUSEsub500)
#'
#' @references
#' Weger BD, Gobet C, David FPA, Atger F et al.
#' 'Systematic analysis of differential rhythmic liver gene expression mediated
#' by the circadian clock and feeding rhythms'.
#' Proc Natl Acad Sci USA 2021 Jan 19;118(3). PMID:33452134. GEO:GSE135898.
#'
#' @examples
#' data(RawCounts_Weger2021_MOUSEsub500)
"RawCounts_Weger2021_MOUSEsub500"
