testthat::test_that("Test DEresultGroup", {
    ##-----------------------------------------------------------------------##
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    ## No time points. We take only two groups for the speed of the example
    RawCounts_T1Wt <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200),
                                                            seq_len(10)]

    ## Preprocessing step
    resDATAprepSEmus1 <- DATAprepSE(RawCounts=RawCounts_T1Wt,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    DESeq2preprocess <- S4Vectors::metadata(resDATAprepSEmus1)$DESeq2obj
    DESeq2obj <- DESeq2preprocess$DESeq2preproceesing
    dds.DE.G <- DESeq2::DESeq(DESeq2obj, quiet=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DEresultGroup(DESeq.result=list(1,2),
                                         LRT.supp.info=TRUE,
                                         pval.min=0.05,
                                         log.FC.min=1),
                           "Res.DE.analysis must be a 'DESeqDataSet' object",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_s4_class(DEresultGroup(DESeq.result=dds.DE.G,
                                            LRT.supp.info=FALSE,
                                            pval.min=0.05,
                                            log.FC.min=1),
                              "DESeqDataSet")## "SummarizedExperiment
})
