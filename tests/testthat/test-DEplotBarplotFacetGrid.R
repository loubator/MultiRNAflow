test_that("multiplication works", {
    Group.ex<-c('G1', 'G2',' G3')
    Time.ex<-c('t1', 't2', 't3', 't4')
    Spe.sign.ex<-c("Pos","Neg")
    Nb.Spe<-sample(3:60, length(Group.ex)*length(Time.ex), replace=FALSE)
    Nb.Spe.sign<-sample(3:60, length(Group.ex)*length(Time.ex)*2,replace=FALSE)
    ##------------------------------------------------------------------------#
    Melt.Dat.1<-data.frame(Group=rep(Group.ex,times=length(Time.ex)),
                           Time=rep(Time.ex,each=length(Group.ex)),
                           Nb.Spe.DE=Nb.Spe)

    expect_s3_class(DEplotBarplotFacetGrid(Data=Melt.Dat.1,
                                           Abs.col=2, Legend.col=2,
                                           Facet.col=1, Value.col=3,
                                           Color.Legend=NULL),
                    "ggplot")
    ##------------------------------------------------------------------------#
    Melt.Dat.2=data.frame(Group=rep(Group.ex,times=length(Time.ex)*2),
                          Time=rep(Time.ex,each=length(Group.ex)*2),
                          Spe.sign=rep(Spe.sign.ex,
                                       times=length(Time.ex)*length(Group.ex)*2),
                          Nb.Spe.DE=Nb.Spe.sign)

    expect_s3_class(DEplotBarplotFacetGrid(Data=Melt.Dat.2,
                                           Abs.col=1,
                                           Legend.col=3,
                                           Facet.col=2,
                                           Value.col=4,
                                           Color.Legend=NULL),
                    "ggplot")
})
