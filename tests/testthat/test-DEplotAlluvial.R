test_that("multiplication DEplotAlluvial", {
    ##------------------------------------------------------------------------#
    set.seed(1994)
    NbTime.vst0<-4
    BinTable<-matrix(sample(c(0,1),replace=TRUE,
                            size=NbTime.vst0*120,c(0.60,0.40)),
                     ncol=NbTime.vst0)
    colnames(BinTable)<-paste0("t", 1:NbTime.vst0)

    ##------------------------------------------------------------------------#
    expect_s3_class(DEplotAlluvial(table.DE.time=BinTable,
                                   Temporal.Group=TRUE)$g.alluvial,
                    "ggplot")

    expect_s3_class(DEplotAlluvial(table.DE.time=BinTable,
                                   Temporal.Group=FALSE),
                    "ggplot")
})
