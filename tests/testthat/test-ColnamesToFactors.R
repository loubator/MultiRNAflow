test_that("Test ColnamesToFactors", {
    ##------------------------------------------------------------------------#
    Data.sim <- RawCountsSimulation(Nb.Group=3,
                                    Nb.Time=2,
                                    Nb.per.GT=3,
                                    Nb.Gene=10)
    subdata <- Data.sim$Sim.dat[, -ncol(Data.sim$Sim.dat)]

    Data.sim2 <- Data.sim
    colnames(Data.sim2$Sim.dat) <- gsub("Ind", "", colnames(Data.sim$Sim.dat),
                                        fixed=TRUE)
    Data.sim2$Vect.Sample <- factor(gsub("Ind", "", Data.sim2$Vect.Sample,
                                         fixed=TRUE))

    DataTonly.sim <- RawCountsSimulation(Nb.Group=1,
                                         Nb.Time=4,
                                         Nb.per.GT=3,
                                         Nb.Gene=10)
    subdataTonly.sim <- DataTonly.sim$Sim.dat[, -ncol(DataTonly.sim$Sim.dat)]

    DataGonly.sim <- RawCountsSimulation(Nb.Group=3,
                                         Nb.Time=1,
                                         Nb.per.GT=3,
                                         Nb.Gene=10)
    subdataGonly.sim <- DataTonly.sim$Sim.dat[, -ncol(DataGonly.sim$Sim.dat)]

    Data1smpl.sim <- RawCountsSimulation(Nb.Group=3,
                                         Nb.Time=2,
                                         Nb.per.GT=1,
                                         Nb.Gene=10)

    DataT1smpl.sim <- RawCountsSimulation(Nb.Group=1,
                                          Nb.Time=2,
                                          Nb.per.GT=1,
                                          Nb.Gene=10)

    DataG1smpl.sim <- RawCountsSimulation(Nb.Group=2,
                                          Nb.Time=1,
                                          Nb.per.GT=1,
                                          Nb.Gene=10)

    ##------------------------------------------------------------------------#
    expect_error(ColnamesToFactors(ExprData=subdata,
                                   Column.gene=1,
                                   Group.position=NULL,
                                   Time.position=NULL,
                                   Individual.position=3),
                 "Samples must belong to at least one time or one group.",
                 fixed=TRUE)

    expect_error(ColnamesToFactors(ExprData=subdata,
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=2,
                                   Individual.position=NULL),
                 "Every sample must have an indidual name (name or number).",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error( ColnamesToFactors(ExprData=Data.sim$Sim.dat,
                                    Column.gene=1.5,
                                    Group.position=1,
                                    Time.position=2,
                                    Individual.position=3),
                  "'Column.gene' must be an integer.",
                  fixed=TRUE)

    expect_error(ColnamesToFactors(ExprData=Data.sim$Sim.dat,
                                   Column.gene=1,
                                   Group.position=1.5,
                                   Time.position=2,
                                   Individual.position=3),
                 "'Group.position' must be an integer.",
                 fixed=TRUE)

    expect_error(ColnamesToFactors(ExprData=Data.sim$Sim.dat,
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=2.5,
                                   Individual.position=3),
                 "'Time.position' must be an integer.",
                 fixed=TRUE)

    expect_error(ColnamesToFactors(ExprData=Data.sim$Sim.dat,
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=2,
                                   Individual.position=3.5),
                 "'Individual.position' must be an integer.",
                 fixed=TRUE)

    expect_error(ColnamesToFactors(ExprData=Data.sim$Sim.dat[,2],
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=2,
                                   Individual.position=3),
                 "'ExprData' must be a matrix of class data.frame.",
                 fixed=TRUE)


    ##------------------------------------------------------------------------#
    expect_error(ColnamesToFactors(ExprData=subdata,
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=2,
                                   Individual.position=3),
                 paste("Every individual must have a unique name,",
                       "must be associated to a unique group and",
                       "must be associated only once to each",
                       "of the same 2 time measurements."),
                 fixed=TRUE)

    expect_error(ColnamesToFactors(ExprData=subdataTonly.sim,
                                   Column.gene=1,
                                   Group.position=NULL,
                                   Time.position=1,
                                   Individual.position=2),
                 paste("Every individual must have a unique name and",
                       "must be associated only once to each",
                       "of the same 4 time measurements."),
                 fixed=TRUE)

    expect_error(ColnamesToFactors(ExprData=subdataGonly.sim,
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=NULL,
                                   Individual.position=2),
                 "Every individual must be associated to only one group.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(ColnamesToFactors(ExprData=Data1smpl.sim$Sim.dat,
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=2,
                                   Individual.position=3),
                 "Each group must have at least two individuals.",
                 fixed=TRUE)

    expect_error(ColnamesToFactors(ExprData=DataT1smpl.sim$Sim.dat,
                                   Column.gene=1,
                                   Group.position=NULL,
                                   Time.position=1,
                                   Individual.position=2),
                 paste("The data must contain the temporal",
                       "expression of at least two individuals."),
                 fixed=TRUE)

    expect_error(ColnamesToFactors(ExprData=DataG1smpl.sim$Sim.dat,
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=NULL,
                                   Individual.position=2),
                 "Each group must have at least two individuals.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_vector(ColnamesToFactors(ExprData=Data.sim2$Sim.dat,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=2,
                                    Individual.position=3)$Final.Name)

})
