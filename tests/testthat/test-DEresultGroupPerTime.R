test_that("Test DEresultGroupPerTime", {
    expect_error(DEresultGroupPerTime(DESeq.result=list(1,2),
                                      LRT.supp.info=TRUE,
                                      pval.min=0.05,
                                      log.FC.min=1),
                 "Res.DE.analysis must a 'DESeqDataSet' object",
                 fixed=TRUE)
})
