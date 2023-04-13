test_that("Test PCApreprocessing", {
    #
    res.sim.count<-RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
                                       Nb.Gene=10)
    #
    expect_error(PCApreprocessing(ExprData=res.sim.count$Sim.dat,
                                  Column.gene=1,
                                  Group.position=NULL,
                                  Time.position=NULL,
                                  Individual.position=3),
                 "Samples must belong to at least one time or one group",
                 fixed=TRUE)
})
