test_that("multiplication DEplotVennBarplotGroup", {
    set.seed(1994)
    ##-------------------------------------------------------------------------#
    ## Binary matrix
    Bin.Table.G<-matrix(c(sample(c(0,1), replace=TRUE, size=240,c(0.75,0.35)),
                          sample(c(0,1), replace=TRUE, size=240,c(0.3,0.7)),
                          rep(0,18)),
                        ncol=6, byrow=TRUE)
    colnames(Bin.Table.G)=c(".A..B.",".A..C.",".A..D.",
                            ".B..C.",".B..D.",".C..D.")
    ##-------------------------------------------------------------------------#
    ## Results
    rVenn <- DEplotVennBarplotGroup(Mat.DE.pair.group=Bin.Table.G)
    expect_s3_class(rVenn$Upset.global,
                    "upset")
})
