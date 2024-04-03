testthat::test_that("Test DEanalysisGroup", {
    ##-----------------------------------------------------------------------##
    ## Data
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
    dds.DE.G <- DESeq2::DESeq(DESeq2obj, quiet=TRUE, betaPrior=FALSE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DEanalysisGroup(DESeq.result=list(1,2),
                                           pval.min=0.05,
                                           log.FC.min=1,
                                           LRT.supp.info=TRUE,
                                           Plot.DE.graph=TRUE,
                                           path.result=NULL,
                                           SubFile.name=NULL),
                           "Res.DE.analysis must be a 'DESeqDataSet' object",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_s4_class(DEanalysisGroup(DESeq.result=dds.DE.G,
                                              pval.min=0.01,
                                              log.FC.min=1,
                                              LRT.supp.info=FALSE,
                                              Plot.DE.graph=TRUE,
                                              path.result=NULL,
                                              SubFile.name="test"),
                              "DESeqDataSet")

})
