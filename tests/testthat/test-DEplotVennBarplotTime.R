test_that("multiplication DEplotVennBarplotTime", {
    set.seed(1994)
    Nb.Time<-4 ## Number of time measurement
    ##-------------------------------------------------------------------------#
    table.DE.time.ex=matrix(sample(c(0,1), replace=TRUE,
                                   size=40*(Nb.Time-1), c(0.2, 0.8)),
                            ncol=Nb.Time-1)
    colnames(table.DE.time.ex)=paste0("t", 1:(Nb.Time-1))
    ##-------------------------------------------------------------------------#
    Log2.FC.matrix.ex=matrix(round(rnorm(n=40*(Nb.Time-1), mean=0, sd=1),
                                   digits=2),
                             ncol=(Nb.Time-1))
    colnames(Log2.FC.matrix.ex)=paste0("t", 1:(Nb.Time-1))
    ##-------------------------------------------------------------------------#
    res.VennBarplot <- DEplotVennBarplotTime(table.DE.time=table.DE.time.ex,
                                             Log2.FC.matrix=Log2.FC.matrix.ex)

    expect_s3_class(res.VennBarplot$Upset.graph,
                    "upset")
})
