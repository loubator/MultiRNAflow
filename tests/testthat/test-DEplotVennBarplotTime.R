testthat::test_that("Test DEplotVennBarplotTime", {
    set.seed(1994)
    Nb.Time <- 4 ## Number of time measurement
    ##-----------------------------------------------------------------------##
    table.DE.time.ex <- matrix(sample(c(0,1), replace=TRUE,
                                      size=40*(Nb.Time-1), c(0.2, 0.8)),
                               ncol=Nb.Time-1)
    colnames(table.DE.time.ex) <- paste0("t", seq_len(Nb.Time-1))
    ##-----------------------------------------------------------------------##
    Log2.FC.matrix.ex <- matrix(round(rnorm(n=40*(Nb.Time-1), mean=0, sd=1),
                                      digits=2),
                                ncol=(Nb.Time-1))
    colnames(Log2.FC.matrix.ex) <- paste0("t", seq_len(Nb.Time-1))
    Log2.FC.matrix.ex2 <- abs(Log2.FC.matrix.ex)
    ##-----------------------------------------------------------------------##
    res.VennBarplot <- DEplotVennBarplotTime(table.DE.time=table.DE.time.ex,
                                             Log2.FC.matrix=Log2.FC.matrix.ex)
    res.VennBarplot2 <- DEplotVennBarplotTime(table.DE.time=table.DE.time.ex,
                                              Log2.FC.matrix=Log2.FC.matrix.ex2)

    testthat::expect_s3_class(res.VennBarplot$Upset.graph, "upset")
    testthat::expect_s3_class(res.VennBarplot2$Upset.graph.with.nb.over,
                              "upset")
})
