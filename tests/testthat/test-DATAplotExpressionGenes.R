test_that("multiplication works", {
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

    expect_error(DATAplotExpressionGenes(SEresNorm=datamus2,
                                         Vector.row.gene=c(1,3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=TRUE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAprepSEmus2,
                                         Vector.row.gene=c(1,3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=TRUE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 Err_SE,
                 fixed=TRUE)
})
