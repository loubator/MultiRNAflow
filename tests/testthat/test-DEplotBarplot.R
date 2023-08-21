test_that("multiplication DEplotBarplot", {
    ## Data simulation
    CrossTabulation<-matrix(c(75,30,10,5, 5,35,5,20, 220,235,285,275),
                            ncol=4, byrow=TRUE)
    colnames(CrossTabulation)<-c("A", "B", "C", "D")
    row.names(CrossTabulation)<-c("Spe.Pos", "Spe.Neg", "Other")

    ##------------------------------------------------------------------------#
    expect_s3_class(DEplotBarplot(ContingencyTable=CrossTabulation,
                                  dodge=TRUE),
                    "ggplot")
    expect_s3_class(DEplotBarplot(ContingencyTable=CrossTabulation,
                                  dodge=FALSE),
                    "ggplot")
})
