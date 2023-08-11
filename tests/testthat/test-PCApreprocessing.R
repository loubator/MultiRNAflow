test_that("Test PCApreprocessing", {
    ##
    ##------------------------------------------------------------------------#
    ## data importation
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2<-RawCounts_Antoszewski2022_MOUSEsub500

    ## preprocessing
    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    ## normalization
    # resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
    #                                      Normalization="rle",
    #                                      Plot.Boxplot=TRUE,
    #                                      Colored.By.Factors=TRUE)

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
})
