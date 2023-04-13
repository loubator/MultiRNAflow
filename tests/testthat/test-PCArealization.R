test_that("Test PCArealization", {
    #
    res.sim.count<-RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
                                       Nb.Gene=10)
    #
    expect_error(PCArealization(ExprData=res.sim.count$Sim.dat,
                                Column.gene=1,
                                Group.position=NULL,
                                Time.position=NULL,
                                Individual.position=3,
                                gene.deletion=NULL,
                                sample.deletion=NULL,
                                Supp.del.sample=FALSE),
                 "Samples must belong to at least one time or one group",
                 fixed=TRUE)
})
