test_that("Test PCAanalysis", {
    #
    res.sim.count<-RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
                                       Nb.Gene=10)
    #
    expect_error(PCAanalysis(ExprData=res.sim.count$Sim.dat,
                             Column.gene=1,
                             Group.position=NULL,
                             Time.position=NULL,
                             Individual.position=3,
                             gene.deletion=NULL,
                             sample.deletion=NULL,
                             Supp.del.sample=FALSE,
                             Plot.PCA=TRUE,
                             Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL,
                             Name.folder.pca=NULL),
                 "Samples must belong to at least one time or one group",
                 fixed=TRUE)
})
