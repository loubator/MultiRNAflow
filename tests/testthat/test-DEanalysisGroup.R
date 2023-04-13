test_that("Test DEanalysisGroup", {
    expect_error(DEanalysisGroup(DESeq.result=list(1,2),
                                 pval.min=0.05,
                                 log.FC.min=1,
                                 LRT.supp.info=TRUE,
                                 Plot.DE.graph=TRUE,
                                 path.result=NULL,
                                 SubFile.name=NULL),
                 "Res.DE.analysis must a 'DESeqDataSet' object",
                 fixed=TRUE)
})
