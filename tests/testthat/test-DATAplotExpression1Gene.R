test_that("Test DATAplotExpression1Gene", {
    ## Simulation raw counts
    resSIMcount <- RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
                                       Nb.Gene=10)
    ## Preprocessing step
    # resDATAprepSE <- DATAprepSE(RawCounts=resSIMcount$Sim.dat,
    #                             Column.gene=1,
    #                             Group.position=1,
    #                             Time.position=2,
    #                             Individual.position=3)

    ##------------------------------------------------------------------------#
    Err_SE <- paste0("'SEres' mut be the results of either the function ",
                     "'DATAprepSE()' or 'DATAnormalization()'.")

    expect_error(DATAplotExpression1Gene(SEres=resSIMcount$Sim.dat,
                                         row.gene=1,
                                         Color.Group=NULL),
                 Err_SE,
                 fixed=TRUE)
})
