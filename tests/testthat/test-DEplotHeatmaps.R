test_that("Test DEplotHeatmaps", {
    expect_error(DEplotHeatmaps(Res.DE.analysis=matrix(0, ncol=2, nrow=3),
                                ColumnsCriteria=2),
                 paste("Res.DE.analysis must be a list or",
                       "a 'DESeqDataSet' object"),
                 fixed=TRUE)
})
