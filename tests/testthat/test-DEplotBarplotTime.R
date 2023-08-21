test_that("Test DEplotBarplotTime", {
    ##------------------------------------------------------------------------#
    Dat1.FTP<-matrix(0,ncol=3, nrow=4)

    expect_error(DEplotBarplotTime(table.DE.time=Dat1.FTP,
                                   Log2.FC.matrix=NULL),
                 "No DE genes", fixed=TRUE)

    ##------------------------------------------------------------------------#
    set.seed(1994)
    Dat1.FTP<-matrix(sample(c(0,1), replace=TRUE, size=120, c(0.3,0.7)), ncol=3)
    Dat2.FTP<-matrix(round(rnorm(n=120, mean=0, sd=1),digits=2), ncol=3)
    colnames(Dat1.FTP)=paste0("t", 1:3)
    colnames(Dat2.FTP)=paste0("t", 1:3)

    expect_s3_class(DEplotBarplotTime(table.DE.time=Dat1.FTP,
                                      Log2.FC.matrix=Dat2.FTP)$g.nb.DEPerTime,
                    "ggplot")

})
