test_that("Test DATAplotExpression1Gene", {
    #
    res.sim.count<-RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
                                       Nb.Gene=10)

    #
    expect_error(DATAplotExpression1Gene(ExprData=res.sim.count$Sim.dat,
                                         row.gene=1,
                                         Column.gene=1,
                                         Group.position=NULL,
                                         Time.position=NULL,
                                         Individual.position=3,
                                         Color.Group=NULL),
                 "Samples must belong to at least one time or one group",
                 fixed=TRUE)
})
