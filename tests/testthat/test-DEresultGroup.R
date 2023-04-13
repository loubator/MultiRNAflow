test_that("Test DEresultGroup", {
    expect_error(DEresultGroup(DESeq.result=list(1,2),
                               LRT.supp.info=TRUE,
                               pval.min=0.05,
                               log.FC.min=1),
                 "Res.DE.analysis must be a 'DESeqDataSet' object",
                 fixed=TRUE)
})
