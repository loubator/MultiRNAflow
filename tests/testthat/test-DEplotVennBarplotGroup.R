testthat::test_that("Test DEplotVennBarplotGroup", {
    set.seed(1994)
    ##-----------------------------------------------------------------------##
    ## Binary matrix
    Bin.Table.G <- matrix(c(sample(c(0,1), replace=TRUE, size=240,
                                   c(0.75,0.35)),
                            sample(c(0,1), replace=TRUE, size=240,
                                   c(0.3,0.7)),
                            rep(0, 18)),
                          ncol=6, byrow=TRUE)
    colnames(Bin.Table.G) <- c(".A..B.",".A..C.",".A..D.",
                               ".B..C.",".B..D.",".C..D.")

    Bin.Table.G2 <- data.frame(Bin.Table.G[,1])
    colnames(Bin.Table.G2) <- c(".A..B.")

    Bin.Table.G3 <- Bin.Table.G
    Bin.Table.G3 <- Bin.Table.G3[-which(apply(Bin.Table.G, 1, sum) == 0),]

    res3apply <- apply(Bin.Table.G3, 1, function(x) paste0(x, collapse=""))
    res3table <- table(res3apply)
    res3names <- names(res3table)[which(as.numeric(res3table) > 3)]

    Bin.Table.G3 <- Bin.Table.G3[which(res3apply%in%res3names),]

    ##-----------------------------------------------------------------------##
    ## Results
    rVenn <- DEplotVennBarplotGroup(Mat.DE.pair.group=Bin.Table.G)
    rVenn2 <- DEplotVennBarplotGroup(Mat.DE.pair.group=Bin.Table.G2)
    rVenn3 <- DEplotVennBarplotGroup(Mat.DE.pair.group=Bin.Table.G3)

    testthat::expect_s3_class(rVenn$Upset.global, "upset")
    testthat::expect_null(rVenn2$Upset.global)
    testthat::expect_s3_class(rVenn3$Upset.global, "upset")
    testthat::expect_null(rVenn3$Upset.threshold)
})
