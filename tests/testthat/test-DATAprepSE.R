test_that("Test DATAprepSE", {
    ##
    ##------------------------------------------------------------------------#
    data("RawCounts_Leong2014_FISSIONsub500wt")
    dataFission <- RawCounts_Leong2014_FISSIONsub500wt

    resCTFfission <- ColnamesToFactors(ExprData=dataFission,
                                       Column.gene=1,
                                       Group.position=NULL,
                                       Time.position=2,
                                       Individual.position=3)

    fissioncDat <- data.frame(Time=resCTFfission$Time.Info,
                              ID=resCTFfission$Individual.info)

    ##------------------------------------------------------------------------#
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500

    resCTFmus2 <- ColnamesToFactors(ExprData=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    mus2cDat1 <- data.frame(Group=resCTFmus2$Group.Info)
    mus2cDat2 <- data.frame(BC=resCTFmus2$Group.Info,
                            ID=resCTFmus2$Individual.info)
    mus2cDat3 <- data.frame(BC=resCTFmus2$Group.Info,
                            Tps=rep(1, times=length(resCTFmus2$Group.Info)),
                            ID=resCTFmus2$Individual.info)

    mus2cDat4 <- mus2cDat2
    colnames(mus2cDat4)[1] <- "Group"

    ##------------------------------------------------------------------------#
    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1.5,
                            Group.position=1,
                            Time.position=NULL,
                            Individual.position=2,
                            colData=NULL),
                 "'Column.gene' must be an integer.",
                 fixed=TRUE)

    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=1.5,
                            Time.position=NULL,
                            Individual.position=2,
                            colData=NULL),
                 "'Group.position' must be an integer.",
                 fixed=TRUE)

    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=1,
                            Time.position=1.5,
                            Individual.position=2,
                            colData=NULL),
                 "'Time.position' must be an integer.",
                 fixed=TRUE)

    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=1,
                            Time.position=NULL,
                            Individual.position=2.5,
                            colData=NULL),
                 "'Individual.position' must be an integer.",
                 fixed=TRUE)

    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=1,
                            Time.position=NULL,
                            Individual.position=NULL,
                            colData=NULL),
                 "Every sample must have an indidual name (name or number).",
                 fixed=TRUE)

    expect_error(DATAprepSE(RawCounts=datamus2[,2],
                            Column.gene=1,
                            Group.position=1,
                            Time.position=NULL,
                            Individual.position=2,
                            colData=NULL),
                 "'RawCounts' must be a matrix of class data.frame.",
                 fixed=TRUE)


    ##------------------------------------------------------------------------#
    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=NULL,
                            Time.position=NULL,
                            Individual.position=2,
                            colData=NULL),
                 paste0("'Time.position', 'Group.position' and colData ",
                        "can not be all 'NULL'"),
                 fixed=TRUE)

    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=2,
                            Time.position=NULL,
                            Individual.position=2,
                            colData=rep(1, times=nrow(mus2cDat1))),
                 "'colData' must be a matrix of class data.frame.",
                 fixed=TRUE)

    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=2,
                            Time.position=NULL,
                            Individual.position=2,
                            colData=mus2cDat1),
                 "'colData' must have two or three coulumns",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=2,
                            Time.position=NULL,
                            Individual.position=2,
                            colData=mus2cDat4[-1,]),
                 paste("The number of rows of 'colData'",
                       "must be equal to the number of",
                       "samples (numbers of column (Nc) of",
                       "'RawCounts' if 'Column.gene==NULL',",
                       "Nc - 1 otherwise)."),
                 fixed=TRUE)

    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=2,
                            Time.position=NULL,
                            Individual.position=2,
                            colData=mus2cDat2),
                 paste0("The column names of 'colData' must be ",
                        "either 'Group', 'ID'",
                        "either 'Time', 'ID'"),
                 fixed=TRUE)

    expect_error(DATAprepSE(RawCounts=datamus2,
                            Column.gene=1,
                            Group.position=2,
                            Time.position=NULL,
                            Individual.position=2,
                            colData=mus2cDat3),
                 paste0("The column names of 'colData' must be ",
                        "'Group' ,'Time', 'ID'"),
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_s4_class(DATAprepSE(RawCounts=datamus2,
                               Column.gene=1,
                               Group.position=1,
                               Time.position=NULL,
                               Individual.position=2,
                               colData=NULL),
                    "SummarizedExperiment")

    expect_s4_class(DATAprepSE(RawCounts=dataFission[, -1],
                               Column.gene=NULL,
                               Group.position=NULL,
                               Time.position=2,
                               Individual.position=3,
                               colData=NULL),
                    "SummarizedExperiment")

    expect_s4_class(DATAprepSE(RawCounts=dataFission[, -1],
                               Column.gene=NULL,
                               Group.position=NULL,
                               Time.position=NULL,
                               Individual.position=NULL,
                               colData=fissioncDat),
                    "SummarizedExperiment")
})
