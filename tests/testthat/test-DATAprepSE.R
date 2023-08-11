test_that("Test DATAprepSE", {
    ##
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2<-RawCounts_Antoszewski2022_MOUSEsub500

    resCTFmus2 <- ColnamesToFactors(ExprData=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    mus2cDat1 <- data.frame(Group=resCTFmus2$Group.Info)
    mus2cDat2 <- data.frame(BC=resCTFmus2$Group.Info,
                            ID=resCTFmus2$Individual.info)

    ##------------------------------------------------------------------------#
    Err_NULL_mus2 <- paste0("'Time.position', 'Group.position' and colData ",
                            "can not be all 'NULL'")

    expect_error(DATAnormalization(DATAprepSE(RawCounts=datamus2,
                                              Column.gene=1,
                                              Group.position=NULL,
                                              Time.position=NULL,
                                              Individual.position=2,
                                              colData=NULL)),
                 Err_NULL_mus2,
                 fixed=TRUE)

    expect_error(DATAnormalization(DATAprepSE(RawCounts=datamus2,
                                              Column.gene=1,
                                              Group.position=2,
                                              Time.position=NULL,
                                              Individual.position=2,
                                              colData=mus2cDat1)),
                 "'colData' must have two or three coulumns",
                 fixed=TRUE)

    ##
    Err_colname_mus2 <- paste0("The column names of 'colData' must be ",
                               "either 'Group', 'ID'",
                               "either 'Time', 'ID'")

    expect_error(DATAnormalization(DATAprepSE(RawCounts=datamus2,
                                              Column.gene=1,
                                              Group.position=2,
                                              Time.position=NULL,
                                              Individual.position=2,
                                              colData=mus2cDat2)),
                 Err_colname_mus2,
                 fixed=TRUE)
})
