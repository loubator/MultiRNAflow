testthat::test_that("Test DEplotBarplotFacetGrid", {
    Group.ex <- c('G1', 'G2',' G3')
    Time.ex <- c('t1', 't2', 't3', 't4')
    Spe.sign.ex <- c("Pos","Neg")

    GtimesT <- length(Group.ex)*length(Time.ex)

    Nb.Spe <- sample(3:60, GtimesT, replace=FALSE)
    Nb.Spe.sign <- sample(3:60, 2*GtimesT, replace=FALSE)

    colorBC <- data.frame(BC=Group.ex, Color=c("blue", "red", "black"))

    ##-----------------------------------------------------------------------##
    Melt.Dat.1 <- data.frame(Group=rep(Group.ex, times=length(Time.ex)),
                             Time=rep(Time.ex, each=length(Group.ex)),
                             Nb.Spe.DE=Nb.Spe)

    testthat::expect_s3_class(DEplotBarplotFacetGrid(Data=Melt.Dat.1,
                                                     Abs.col=2, Legend.col=2,
                                                     Facet.col=1, Value.col=3,
                                                     Color.Legend=NULL),
                              "ggplot")

    ##-----------------------------------------------------------------------##
    Melt.Dat.2 <- data.frame(Group=rep(Group.ex, times=length(Time.ex)*2),
                             Time=rep(Time.ex, each=length(Group.ex)*2),
                             Spe.sign=rep(Spe.sign.ex, times=2*GtimesT),
                             Nb.Spe.DE=Nb.Spe.sign)

    testthat::expect_s3_class(DEplotBarplotFacetGrid(Data=Melt.Dat.2,
                                                     Abs.col=1,
                                                     Legend.col=3,
                                                     Facet.col=2,
                                                     Value.col=4,
                                                     Color.Legend=NULL),
                              "ggplot")

    testthat::expect_s3_class(DEplotBarplotFacetGrid(Data=Melt.Dat.2,
                                                     Abs.col=1,
                                                     Legend.col=3,
                                                     Facet.col=2,
                                                     Value.col=4,
                                                     Color.Legend=colorBC),
                              "ggplot")
})
