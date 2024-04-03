testthat::test_that("Test DEplotAlluvial", {
    ##-----------------------------------------------------------------------##
    set.seed(1994)
    NbTime.vst0 <- 4
    BinTable <- matrix(sample(c(0, 1), replace=TRUE,
                              size=NbTime.vst0*120,
                              prob=c(0.60, 0.40)),
                       ncol=NbTime.vst0)
    colnames(BinTable) <- paste0("t", seq_len(NbTime.vst0))

    BinTable2 <- BinTable
    row.names(BinTable2) <- paste0("Gene", seq_len(nrow(BinTable)))
    colnames(BinTable2) <- NULL

    ##-----------------------------------------------------------------------##
    testthat::expect_s3_class(DEplotAlluvial(table.DE.time=BinTable,
                                             Temporal.Group=TRUE)$g.alluvial,
                              "ggplot")

    testthat::expect_s3_class(DEplotAlluvial(table.DE.time=BinTable2,
                                             Temporal.Group=TRUE,
                                             title.alluvial="TitleTest",
                                             title.evolution="tt2")$g.alluvial,
                              "ggplot")

    testthat::expect_s3_class(DEplotAlluvial(table.DE.time=BinTable,
                                             Temporal.Group=FALSE),
                              "ggplot")

    testthat::expect_s3_class(DEplotAlluvial(table.DE.time=BinTable,
                                             Temporal.Group=FALSE,
                                             title.alluvial="TitleTest",
                                             title.evolution="TitleTest2"),
                              "ggplot")
})
