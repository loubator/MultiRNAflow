testthat::test_that("Test PCApreprocessing", {
    ##
    ##-----------------------------------------------------------------------##
    data("RawCounts_Leong2014_FISSIONsub500wt")
    dataFission <- RawCounts_Leong2014_FISSIONsub500wt[seq_len(100),]

    ## preprocessing
    SEresPrepFission <- DATAprepSE(RawCounts=dataFission,
                                   Column.gene=1,
                                   Group.position=NULL,
                                   Time.position=2,
                                   Individual.position=3)

    ## normalization
    SEresNormFission <- DATAnormalization(SEres=SEresPrepFission,
                                          Normalization="rle",
                                          Plot.Boxplot=FALSE,
                                          Colored.By.Factors=FALSE)

    ##-----------------------------------------------------------------------##
    ## data importation
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(100),]

    ## preprocessing
    SEresPrepMus2 <- DATAprepSE(RawCounts=datamus2,
                                Column.gene=1,
                                Group.position=1,
                                Time.position=NULL,
                                Individual.position=2)

    ## normalization
    SEresNormMus2 <- DATAnormalization(SEres=SEresPrepMus2,
                                       Normalization="rle",
                                       Plot.Boxplot=TRUE,
                                       Colored.By.Factors=TRUE)

    SEresNormMus2Id <- SEresNormMus2
    S4Vectors::metadata(SEresNormMus2Id)$SEidentification <- NULL

    ##-----------------------------------------------------------------------##
    data("RawCounts_Schleiss2021_CLLsub500")

    col0Leuk <- seq(2, 55, 9)
    col1Leuk <- seq(3, 55, 9)
    ## col2Leuk <- seq(4, 55, 9)
    col3tLeuk <- sort(c(1, col0Leuk, col1Leuk))## , col2Leuk

    dataLeuk <- RawCounts_Schleiss2021_CLLsub500[seq_len(27), col3tLeuk]

    ## preprocessing
    SEresPrepLeuk <- DATAprepSE(RawCounts=dataLeuk,
                                Column.gene=1,
                                Group.position=2,
                                Time.position=4,
                                Individual.position=3)

    ## normalization
    SEresNormLeuk <- DATAnormalization(SEres=SEresPrepLeuk,
                                       Normalization="rle",
                                       Plot.Boxplot=TRUE,
                                       Colored.By.Factors=TRUE)

    ##-----------------------------------------------------------------------##
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    testthat::expect_error(PCApreprocessing(SEresNorm=datamus2,
                                            DATAnorm=TRUE),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(PCApreprocessing(SEresNorm=SEresPrepMus2,
                                            DATAnorm=TRUE),
                           Err_SE,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(PCApreprocessing(SEresNorm=SEresNormMus2,
                                            DATAnorm=1.1),
                           "'DATAnorm' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(PCApreprocessing(SEresNorm=SEresNormMus2Id,
                                            DATAnorm=TRUE),
                           paste0("'SEresNorm' mut be the results of ",
                                  "the function 'DATAnormalization().'"),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    res1_PCAprep <- PCApreprocessing(SEresNorm=SEresNormMus2, DATAnorm=FALSE)
    res2_PCAprep <- PCApreprocessing(SEresNorm=SEresNormFission, DATAnorm=TRUE)
    res3_PCAprep <- PCApreprocessing(SEresNorm=SEresNormLeuk, DATAnorm=TRUE)

    testthat::expect_s4_class(res1_PCAprep, "SummarizedExperiment")
    testthat::expect_s4_class(res2_PCAprep, "SummarizedExperiment")
    testthat::expect_s4_class(res3_PCAprep, "SummarizedExperiment")

})
