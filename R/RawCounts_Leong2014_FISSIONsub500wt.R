#' @title Yeast times series raw counts data after stimulation
#' with and without silencing
#'
#' @description Raw counts data for fission yeast RNA-Seq experiment with
#' two groups (wt and mut), 6 times (0, 15min, 30min, 60min, 120min, 180min)
#' and 3 replicates for each group and time.
#' The original dataset has 7039 genes but we kept only 500 genes
#' in order to increase the speed of each function in our algorithm.
#'
#' @format A data frame with 500 rows (genes) and 37 columns (samples).
#' The column names are as follow
#'  \describe{
#'            \item{Gene}{Gene name}
#'            \item{wt_t0_r1}{The sample is the first replica (r1) of the biological condition control (wt) at time t0 (0 min)}
#'            \item{wt_t0_r2}{The sample is the second replica (r2) of the biological condition control (wt) at time t0 (0 min)}
#'            \item{wt_t0_r3}{The sample is the third replica (r3) of the biological condition control (wt) at time t0 (0 min)}
#'            \item{wt_t1_r1}{The sample is the first replica (r1) of the biological condition control (wt) at time t1 (15 min)}
#'            \item{wt_t1_r2}{The sample is the second replica (r2) of the biological condition control (wt) at time t1 (15 min)}
#'            \item{wt_t1_r3}{The sample is the third replica (r3) of the biological condition control (wt) at time t1 (15 min)}
#'            \item{wt_t2_r1}{The sample is the first replica (r1) of the biological condition control (wt) at time t2 (30 min)}
#'            \item{wt_t2_r2}{The sample is the second replica (r2) of the biological condition control (wt) at time t2 (30 min)}
#'            \item{wt_t2_r3}{The sample is the third replica (r3) of the biological condition control (wt) at time t2 (30 min)}
#'            \item{wt_t3_r1}{The sample is the first replica (r1) of the biological condition control (wt) at time t3 (60 min)}
#'            \item{wt_t3_r2}{The sample is the second replica (r2) of the biological condition control (wt) at time t3 (60 min)}
#'            \item{wt_t3_r3}{The sample is the third replica (r3) of the biological condition control (wt) at time t3 (60 min)}
#'            \item{wt_t4_r1}{The sample is the first replica (r1) of the biological condition control (wt) at time t4 (120 min)}
#'            \item{wt_t4_r2}{The sample is the second replica (r2) of the biological condition control (wt) at time t4 (120 min)}
#'            \item{wt_t4_r3}{The sample is the third replica (r3) of the biological condition control (wt) at time t4 (120 min)}
#'            \item{wt_t5_r1}{The sample is the first replica (r1) of the biological condition control (wt) at time t5 (180 min)}
#'            \item{wt_t5_r2}{The sample is the second replica (r2) of the biological condition control (wt) at time t5 (180 min)}
#'            \item{wt_t5_r3}{The sample is the third replica (r3) of the biological condition control (wt) at time t5 (180 min)}
#'  }
#'
#' @details
#' The following is quoted from the GEO series GSE56761 (link in source):
#'
#' Summary:
#' "Mitogen Activated Protein Kinase (MAPK) signaling cascades transduce
#' information arising from events external to the cell,
#' such as environmental stresses, to a variety of downstream effectors and
#' transcription factors.
#' The fission yeast stress activated MAP kinase (SAPK) pathway is conserved
#' with the p38 and JNK pathways in humans, and comprises the MAPKKKs Win1,
#' Wis4, the MAPKK Wis1, and the MAPK, Sty1. Sty1 and its main downstream
#' effector Atf1 regulate a large set of core environmental stress
#' response genes.
#' The fission yeast genome encodes three other ATF proteins: Atf21, Atf31
#' and Pcr1. Among these, atf21 is specifically induced under conditions of
#' high osmolarity.
#' We have therefore instigated a programme to investigate the role played by
#' non coding RNAs (ncRNAs) in response to osmotic stress challenge in wild
#' type and atf21Delta cells.
#' By integrating global proteomics and RNA sequencing data, we identified
#' a systematic program in which elevated antisense RNAs arising both
#' from ncRNAs and from 3'-overlapping convergent gene-pairs is directly
#' associated with substantial reductions in protein levels throughout
#' the fission yeast genome. We also found an xtensive array of ncRNAs with
#' trans associations that have the potential to influence different
#' biological processes and stress responses in fission yeast,
#' suggesting ncRNAs comprise additional components of the SAPK regulatory
#' system".
#'
#' Overall design:
#' "Global transcription profiles of fission yeast wild type (WT) and
#' atf21del strains over an osmotic stress time course following treatment
#' with 1M sorbitol at 0, 15, 30, 60, 120 and 180 mins.
#' Strand-specific single end sequencing of total RNA was performed in
#' biological triplicates on the Applied Biosystems
#' SOLiD 5500xl Genetic Analyzer System".
#'
#' We kept 500 genes only in order to increase the speed for each example.
#'
#'
#' @source {This dataset can be found in the R Package fission.
#' \url{https://bioconductor.org/packages/release/data/experiment/html/fission.html}
#' Link of GEO series GSE56761:
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56761}.
#' The name of the samples was renamed in order to be used with our package.}
#'
#' @return Yeast dataset with 6 time measurements.
#'
#' @usage data(RawCounts_Leong2014_FISSIONsub500wt)
#'
#' @references
#' Leong HS, Dawson K, Wirth C, Li Y et al.
#' 'A global non-coding RNA system modulates fission yeast protein levels
#' in response to stress'.
#' Nat Commun 2014 May 23;5:3947. PMID:24853205. GEO:GSE56761.
#'
#' @examples
#' data(RawCounts_Leong2014_FISSIONsub500wt)
"RawCounts_Leong2014_FISSIONsub500wt"
