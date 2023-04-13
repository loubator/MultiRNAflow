test_that("Test MFUZZclustersNumber", {
    data.clust.sim<-matrix(0, nrow=30, ncol=12)
    #
    expect_error(MFUZZclustersNumber(ExprData=data.clust.sim,
                                     Column.gene=NULL,
                                     Group.position=1,
                                     Time.position=2,
                                     Individual.position=3,
                                     Method="hcpc",
                                     Max.clust=1,
                                     Plot.Cluster=TRUE,
                                     path.result=NULL),
                 "'Max.clust' must be an integer greater or equal to 2.",
                 fixed=TRUE)
    #
    expect_error(MFUZZclustersNumber(ExprData=data.clust.sim,
                                     Column.gene=NULL,
                                     Group.position=1,
                                     Time.position=2,
                                     Individual.position=3,
                                     Method="hcpc",
                                     Max.clust=40,
                                     Plot.Cluster=TRUE,
                                     path.result=NULL),
                 "'Max.clust' must be an integer greater or equal to 2.",
                 fixed=TRUE)
})
