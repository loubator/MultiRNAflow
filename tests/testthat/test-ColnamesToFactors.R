test_that("Test ColnamesToFactors", {
    #
    Data.sim<-RawCountsSimulation(Nb.Group=3,
                                  Nb.Time=2,
                                  Nb.per.GT=3,
                                  Nb.Gene=10)
    #
    subdata<-Data.sim$Sim.dat[, -ncol(Data.sim$Sim.dat)]
    #
    expect_error(ColnamesToFactors(ExprData=subdata,
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=2,
                                   Individual.position=3),
                 paste("Every individual must have a unique name,",
                       "must be associated to a unique group and",
                       "must be associated only once to each",
                       "of the same 2 time measurements"),
                 fixed=TRUE)
})
