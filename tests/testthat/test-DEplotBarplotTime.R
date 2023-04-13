test_that("Test DEplotBarplotTime", {
    Dat1.FTP<-matrix(0,ncol=3, nrow=4)

    expect_error(DEplotBarplotTime(table.DE.time=Dat1.FTP,
                                   Log2.FC.matrix=NULL),
                 "No DE genes", fixed=TRUE)
})
