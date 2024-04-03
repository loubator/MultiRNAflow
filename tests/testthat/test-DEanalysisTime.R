testthat::test_that("Test DEanalysisTime", {
    ##-----------------------------------------------------------------------##
    data("RawCounts_Leong2014_FISSIONsub500wt")
    ## We take only the first three time for the speed of the example
    RawCounts_Fission_3t <- RawCounts_Leong2014_FISSIONsub500wt[seq_len(200),
                                                                seq_len(10)]

    ## Preprocessing step
    resDATAprepSEfission <- DATAprepSE(RawCounts=RawCounts_Fission_3t,
                                       Column.gene=1,
                                       Group.position=NULL,
                                       Time.position=2,
                                       Individual.position=3)

    DESeq2preprocess <- S4Vectors::metadata(resDATAprepSEfission)$DESeq2obj
    DESeq2obj <- DESeq2preprocess$DESeq2preproceesing
    ddsDEt <- DESeq2::DESeq(DESeq2obj, quiet=TRUE, test="LRT", reduced=~1)

    ##-----------------------------------------------------------------------##
    RawCounts_Fission_2t <- RawCounts_Leong2014_FISSIONsub500wt[seq_len(200),
                                                                seq_len(7)]

    ## Preprocessing step
    resDATAprepSEfission2t <- DATAprepSE(RawCounts=RawCounts_Fission_2t,
                                         Column.gene=1,
                                         Group.position=NULL,
                                         Time.position=2,
                                         Individual.position=3)

    DESeq2preprocess2t <- S4Vectors::metadata(resDATAprepSEfission2t)$DESeq2obj
    DESeq2obj2t <- DESeq2preprocess2t$DESeq2preproceesing

    ddsDE2t <- DESeq2::DESeq(DESeq2obj2t, quiet=TRUE, test="LRT", reduced=~1)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DEanalysisTime(DESeq.result=list(1, 2),
                                          pval.min=0.05,
                                          pval.vect.t=NULL,
                                          log.FC.min=1,
                                          LRT.supp.info = FALSE,
                                          Plot.DE.graph=TRUE,
                                          path.result=NULL,
                                          SubFile.name=NULL),
                           "Res.DE.analysis must be a 'DESeqDataSet' object",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_s4_class(DEanalysisTime(DESeq.result=ddsDEt,
                                             pval.min=0.05,
                                             pval.vect.t=c(0.01, 0.05, 0.05),
                                             log.FC.min=1,
                                             LRT.supp.info=FALSE,
                                             Plot.DE.graph=TRUE,
                                             path.result=NULL,
                                             SubFile.name="test"),
                              "DESeqDataSet")

    testthat::expect_s4_class(DEanalysisTime(DESeq.result=ddsDEt,
                                             pval.min=0.05,
                                             pval.vect.t=c(0.01),
                                             log.FC.min=1,
                                             LRT.supp.info=TRUE,
                                             Plot.DE.graph=FALSE,
                                             path.result=NULL,
                                             SubFile.name=NULL),
                              "DESeqDataSet")

    testthat::expect_s4_class(DEanalysisTime(DESeq.result=ddsDE2t,
                                             pval.min=0.05,
                                             pval.vect.t=NULL,
                                             log.FC.min=1,
                                             LRT.supp.info=TRUE,
                                             Plot.DE.graph=TRUE,
                                             path.result=NULL,
                                             SubFile.name=NULL),
                              "DESeqDataSet")

})
