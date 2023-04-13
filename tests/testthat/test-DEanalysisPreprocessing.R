test_that("Test DEanalysisPreprocessing", {
    #
    expect_error(DEanalysisPreprocessing(RawCounts=matrix(0, ncol=3,
                                                          nrow=2),
                                         Column.gene=1,
                                         Group.position=NULL,
                                         Time.position=NULL,
                                         Individual.position=2),
                 "'Time.position' and 'Group.position' can not be both NULL",
                 fixed=TRUE)
})
