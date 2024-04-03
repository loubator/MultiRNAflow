testthat::test_that("Test DEplotBarplot", {
    ## Data simulation
    CrossTabulation <- matrix(c(75,30,10,5, 5,35,5,20, 220,235,285,275),
                              ncol=4, byrow=TRUE)
    colnames(CrossTabulation) <- c("A", "B", "C", "D")
    row.names(CrossTabulation) <- c("Spe.Pos", "Spe.Neg", "Other")

    CrossTabulation2 <- CrossTabulation[-3,]

    ##------------------------------------------------------------------------#
    testthat::expect_s3_class(DEplotBarplot(ContingencyTable=CrossTabulation,
                                            dodge=TRUE),
                              "ggplot")
    testthat::expect_s3_class(DEplotBarplot(ContingencyTable=CrossTabulation,
                                            dodge=FALSE),
                              "ggplot")
    testthat::expect_s3_class(DEplotBarplot(ContingencyTable=CrossTabulation2,
                                            dodge=TRUE),
                              "ggplot")
})
