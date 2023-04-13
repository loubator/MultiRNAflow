test_that("Test DEanalysisTimeAndGroup", {
    expect_error(DEanalysisTimeAndGroup(DESeq.result=list(1, 2),
                                        pval.min=0.05,
                                        pval.vect.t=NULL,
                                        log.FC.min=1,
                                        LRT.supp.info=FALSE,
                                        Plot.DE.graph=TRUE,
                                        path.result=NULL,
                                        SubFile.name=NULL),
                 "Res.DE.analysis must be a 'DESeqDataSet' object",
                 fixed=TRUE)
})
