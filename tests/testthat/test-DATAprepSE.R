testthat::test_that("Test DATAprepSE", {
    ##
    ##-----------------------------------------------------------------------##
    data("Transcript_HomoSapiens_Database")
    data("RawCounts_Schleiss2021_CLLsub500")

    col0Leuk <- seq(2, 55, 9)
    col1Leuk <- seq(3, 55, 9)
    col2Leuk <- seq(4, 55, 9)
    col3tLeuk <- sort(c(1, col0Leuk, col1Leuk, col2Leuk))

    dataLeuk <- RawCounts_Schleiss2021_CLLsub500[seq_len(27), col3tLeuk]

    dataLeuk2 <- RawCounts_Schleiss2021_CLLsub500[seq_len(27)[-c(2, 10, 20)],
                                                  c(col3tLeuk)]

    dataLeuk2dupliGenes <- dataLeuk2
    dataLeuk2dupliGenes[2, 1] <- dataLeuk2dupliGenes[1, 1]

    BgCdEx <- rep(c("P", "NP"), each=9)
    TimeEx <- rep(paste0("t", seq_len(3) - 1), times=6)
    IndvEx <- rep(paste0("pcl", seq_len(6)), each=3)
    colDataEx <- data.frame(Group=BgCdEx, Time=TimeEx, ID=IndvEx)

    ##-----------------------------------------------------------------------##
    data("RawCounts_Leong2014_FISSIONsub500wt")
    dataFission <- RawCounts_Leong2014_FISSIONsub500wt[seq_len(20),]

    resCTFfission <- ColnamesToFactors(ExprData=dataFission,
                                       Column.gene=1,
                                       Group.position=NULL,
                                       Time.position=2,
                                       Individual.position=3)

    fissioncDat <- data.frame(Time=resCTFfission$Time.Info,
                              ID=resCTFfission$Individual.info)

    ##-----------------------------------------------------------------------##
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(20),]

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

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1.5,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL),
                           "'Column.gene' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene="column1",
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL),
                           "'Column.gene' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=dataLeuk2dupliGenes,
                                      Column.gene=1,
                                      Group.position=2,
                                      Time.position=4,
                                      Individual.position=3,
                                      colData=NULL),
                           "There are 1 duplicated genes: ABCA7",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1.5,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL),
                           "'Group.position' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=1.5,
                                      Individual.position=2,
                                      colData=NULL),
                           "'Time.position' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2.5,
                                      colData=NULL),
                           paste("'Individual.position' must be",
                                 "a non-negative integer."),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=NULL,
                                      colData=NULL),
                           paste0("Every sample must have an indidual name ",
                                  "(name or number)."),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2[,2],
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL),
                           "'RawCounts' must be a matrix of class data.frame.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=NULL,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL),
                           paste0("'Time.position', 'Group.position' and ",
                                  "'colData' can not be all 'NULL'"),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=2,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=rep(1, times=nrow(mus2cDat1))),
                           "'colData' must be a matrix of class data.frame.",
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=2,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=mus2cDat1),
                           "'colData' must have two or three coulumns",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
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

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=2,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=mus2cDat2),
                           paste0("The column names of 'colData' must be ",
                                  "either 'Group', 'ID'",
                                  "either 'Time', 'ID'"),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=2,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=mus2cDat3),
                           paste0("The column names of 'colData' must be ",
                                  "'Group' ,'Time', 'ID'"),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL,
                                      RNAlength="chr"),
                           paste0("'RNAlength' must be either NULL, ",
                                  "'hsapiens' or a two columns data.frame."),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=NULL,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL,
                                      RNAlength="hsapiens"),
                           paste0("If 'RNAlength' is not NULL, ",
                                  "'Column.gene' can not be NULL too."),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL,
                                      RNAlength=data.frame(c(1 ,2),
                                                           c(1, 2),
                                                           c(1, 2))),
                           "'RNAlength' must be a two columns data.frame.",
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL,
                                      RNAlength=data.frame(c(1 ,2), c(1, 2))),
                           paste0("The first column of 'RNAlength' must be ",
                                  "character and the second numeric."),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL,
                                      RNAlength=data.frame(c("Gene1", "Gene2"),
                                                           c(1, 2))),
                           "No corresponding genes.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL,
                                      SUMfilter="chr"),
                           paste0("'SUMfilter' must be a positive numeric ",
                                  "value"),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL,
                                      SUMfilter=-1),
                           paste0("'SUMfilter' must be a positive numeric ",
                                  "value"),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL,
                                      VARfilter="chr"),
                           paste0("'VARfilter' must be a positive numeric ",
                                  "value"),
                           fixed=TRUE)

    testthat::expect_error(DATAprepSE(RawCounts=datamus2,
                                      Column.gene=1,
                                      Group.position=1,
                                      Time.position=NULL,
                                      Individual.position=2,
                                      colData=NULL,
                                      VARfilter=-1),
                           paste0("'VARfilter' must be a positive numeric ",
                                  "value"),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    res1_prepSE <- DATAprepSE(RawCounts=datamus2,
                              Column.gene=1, Group.position=1,
                              Time.position=NULL, Individual.position=2,
                              colData=NULL, SUMfilter=200)

    res2_prepSE <- DATAprepSE(RawCounts=dataFission[, -1],
                              Column.gene=NULL, Group.position=NULL,
                              Time.position=2, Individual.position=3,
                              colData=NULL, VARfilter=500)

    res3_prepSE <- DATAprepSE(RawCounts=dataFission[, -1],
                              Column.gene=NULL, Group.position=NULL,
                              Time.position=NULL, Individual.position=NULL,
                              colData=fissioncDat)

    res4_prepSE <- DATAprepSE(RawCounts=dataLeuk,
                              Column.gene=1, Group.position=2,
                              Time.position=4, Individual.position=3,
                              colData=NULL, RNAlength="hsapiens")

    res5_prepSE <- DATAprepSE(RawCounts=dataLeuk2,
                              Column.gene=1, Group.position=2,
                              Time.position=4, Individual.position=3,
                              colData=colDataEx, RNAlength="hsapiens")

    res6_prepSE <- DATAprepSE(RawCounts=datamus2,
                              Column.gene=1, Group.position=NULL,
                              Time.position=NULL, Individual.position=NULL,
                              colData=mus2cDat4)

    testthat::expect_s4_class(res1_prepSE, "SummarizedExperiment")
    testthat::expect_s4_class(res2_prepSE, "SummarizedExperiment")
    testthat::expect_s4_class(res3_prepSE, "SummarizedExperiment")
    testthat::expect_s4_class(res4_prepSE, "SummarizedExperiment")
    testthat::expect_s4_class(res5_prepSE, "SummarizedExperiment")
    testthat::expect_s4_class(res6_prepSE, "SummarizedExperiment")
})
