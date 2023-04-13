test_that("Test GSEApreprocessing", {
    expect_error(GSEApreprocessing(Res.DE.analysis=matrix(0,ncol=2,nrow=3),
                                   ColumnsCriteria=2,
                                   Set.Operation="union",
                                   Rnk.files=TRUE,
                                   Save.files=FALSE),
                 paste("Res.DE.analysis must be a list or",
                       "a 'DESeqDataSet' object"),
                 fixed=TRUE)
})
