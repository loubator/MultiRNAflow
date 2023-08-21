test_that("Test DATAplotExpression1Gene", {
    ## Simulation raw counts
    resSIMcount <- RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
                                       Nb.Gene=10)
    # Preprocessing step
    resDATAprepSE <- DATAprepSE(RawCounts=resSIMcount$Sim.dat,
                                Column.gene=1,
                                Group.position=1,
                                Time.position=2,
                                Individual.position=3)

    resDATAprepNOid <- resDATAprepSE
    S4Vectors::metadata(resDATAprepNOid)$SEidentification <- "Test"

    resDATAprepNULLid <- resDATAprepSE
    NewMETA <- S4Vectors::metadata(resDATAprepNULLid)[seq_len(4)]
    S4Vectors::metadata(resDATAprepNULLid) <- NewMETA

    ##------------------------------------------------------------------------#
    Err_SE <- paste0("'SEres' mut be the results of either the function ",
                     "'DATAprepSE()' or 'DATAnormalization()'.")

    expect_error(DATAplotExpression1Gene(SEres=resSIMcount$Sim.dat,
                                         row.gene=1,
                                         Color.Group=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(DATAplotExpression1Gene(SEres=resDATAprepNOid,
                                         row.gene=1,
                                         Color.Group=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(DATAplotExpression1Gene(SEres=resDATAprepNULLid,
                                         row.gene=1,
                                         Color.Group=NULL),
                 Err_SE,
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(DATAplotExpression1Gene(SEres=resDATAprepSE,
                                         row.gene=1.5,
                                         Color.Group=NULL),
                 "'row.gene' must be a non negative integer.",
                 fixed=TRUE)

    expect_error(DATAplotExpression1Gene(SEres=resDATAprepSE,
                                         row.gene=-1,
                                         Color.Group=NULL),
                 "'row.gene' must be a non negative integer.",
                 fixed=TRUE)

    expect_error(DATAplotExpression1Gene(SEres=resDATAprepSE,
                                         row.gene="Try",
                                         Color.Group=NULL),
                 "'row.gene' must be a non negative integer.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(DATAplotExpression1Gene(SEres=resDATAprepSE,
                                         row.gene=1,
                                         Color.Group=1.5),
                 "'Color.Group' must be NULL or a data.frame.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_s3_class(DATAplotExpression1Gene(SEres=resDATAprepSE,
                                            row.gene=1,
                                            Color.Group=NULL),
                    "ggplot")
})
