test_that("Test DEplotVolcanoMA", {
    expect_error(DEplotVolcanoMA(Res.DE.analysis=matrix(0, nrow=3, ncol=2),
                                 NbGene.plotted=2,
                                 SizeLabel=3,
                                 Display.plots=FALSE,
                                 Save.plots=FALSE),
                 paste("Res.DE.analysis must be a list or",
                       "a 'DESeqDataSet' object"),
                 fixed=TRUE)
})
