test_that("Test PCApreprocessing", {
    ##
    ##------------------------------------------------------------------------#
    data("RawCounts_Leong2014_FISSIONsub500wt")
    dataFission <- RawCounts_Leong2014_FISSIONsub500wt[seq_len(100),]

    ## preprocessing
    resDATAprepSEfission <- DATAprepSE(RawCounts=dataFission,
                                       Column.gene=1,
                                       Group.position=NULL,
                                       Time.position=2,
                                       Individual.position=3)

    ## normalization
    resDATAnormFission <- DATAnormalization(SEres=resDATAprepSEfission,
                                            Normalization="rle",
                                            Plot.Boxplot=FALSE,
                                            Colored.By.Factors=FALSE)

    ##------------------------------------------------------------------------#
    ## data importation
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2<-RawCounts_Antoszewski2022_MOUSEsub500[seq_len(100),]

    ## preprocessing
    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    ## normalization
    resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
                                         Normalization="rle",
                                         Plot.Boxplot=TRUE,
                                         Colored.By.Factors=TRUE)

    ##------------------------------------------------------------------------#
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    expect_error(PCApreprocessing(SEresNorm=datamus2,
                                  DATAnorm=TRUE),
                 Err_SE,
                 fixed=TRUE)

    expect_error(PCApreprocessing(SEresNorm=resDATAprepSEmus2,
                                  DATAnorm=TRUE),
                 Err_SE,
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(PCApreprocessing(SEresNorm=resDATAnormMus2,
                                  DATAnorm=1.1),
                 "'DATAnorm' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_s4_class(PCApreprocessing(SEresNorm=resDATAnormMus2,
                                     DATAnorm=FALSE),
                    "SummarizedExperiment")

    expect_s4_class(PCApreprocessing(SEresNorm=resDATAnormFission,
                                     DATAnorm=TRUE),
                    "SummarizedExperiment")
})
