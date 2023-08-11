test_that("Test DATAplotBoxplotSamples", {
    ##
    ##------------------------------------------------------------------------#
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2<-RawCounts_Antoszewski2022_MOUSEsub500

    # resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
    #                                 Column.gene=1,
    #                                 Group.position=1,
    #                                 Time.position=NULL,
    #                                 Individual.position=2)

    ##------------------------------------------------------------------------#
    Err_SEmus2 <- paste0("'SEres' mut be the results of either the function ",
                         "'DATAprepSE()' or 'DATAnormalization()'.")

    expect_error(DATAplotBoxplotSamples(SEres=datamus2,
                                        Log2.transformation=TRUE,
                                        Colored.By.Factors=TRUE,
                                        Color.Group=NULL,
                                        Plot.genes=FALSE,
                                        y.label=NULL),
                 Err_SEmus2,
                 fixed=TRUE)
})
