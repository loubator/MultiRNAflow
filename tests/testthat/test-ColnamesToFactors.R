testthat::test_that("Test ColnamesToFactors", {
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
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

    ##-----------------------------------------------------------------------##
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

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    testthat::expect_error(ErrNNI(NNI="NNItest", NNIname="NNItest"),
                           "'NNItest' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(ErrNNI(NNI=-1, NNIname="NNItest"),
                           "'NNItest' must be a non-negative integer.",
                           fixed=TRUE)

    res_ErrNNI <- ErrNNI(NNI=1, NNIname="NNItest")
    testthat::expect_equal(res_ErrNNI, "No error")

    ##-----------------------------------------------------------------------##
    res_ErrPosition <- ErrPosition(Column.gene=1, Group.position=2,
                                   Time.position=3, Individual.position=4)
    testthat::expect_equal(res_ErrPosition, "No error")

    ##-----------------------------------------------------------------------##
    testthat::expect_error(ColnamesToFactors(ExprData=subdata,
                                             Column.gene=1,
                                             Group.position=NULL,
                                             Time.position=NULL,
                                             Individual.position=3),
                           paste0("Samples must belong to at least ",
                                  "one time or one group."),
                           fixed=TRUE)

    testthat::expect_error(ColnamesToFactors(ExprData=subdata,
                                             Column.gene=1,
                                             Group.position=1,
                                             Time.position=2,
                                             Individual.position=NULL),
                           paste0("Every sample must have an indidual name ",
                                  "(name or number)."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error( ColnamesToFactors(ExprData=Data.sim$Sim.dat,
                                              Column.gene=1.5,
                                              Group.position=1,
                                              Time.position=2,
                                              Individual.position=3),
                            "'Column.gene' must be a non-negative integer.",
                            fixed=TRUE)

    testthat::expect_error(ColnamesToFactors(ExprData=Data.sim$Sim.dat,
                                             Column.gene=1,
                                             Group.position=1.5,
                                             Time.position=2,
                                             Individual.position=3),
                           "'Group.position' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(ColnamesToFactors(ExprData=Data.sim$Sim.dat,
                                             Column.gene=1,
                                             Group.position=1,
                                             Time.position=2.5,
                                             Individual.position=3),
                           "'Time.position' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(ColnamesToFactors(ExprData=Data.sim$Sim.dat,
                                             Column.gene=1,
                                             Group.position=1,
                                             Time.position=2,
                                             Individual.position=3.5),
                           paste("'Individual.position' must be",
                                 "a non-negative integer."),
                           fixed=TRUE)

    testthat::expect_error(ColnamesToFactors(ExprData=Data.sim$Sim.dat[,2],
                                             Column.gene=1,
                                             Group.position=1,
                                             Time.position=2,
                                             Individual.position=3),
                           "'ExprData' must be a matrix of class data.frame.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(ColnamesToFactors(ExprData=subdata,
                                             Column.gene=1,
                                             Group.position=1,
                                             Time.position=2,
                                             Individual.position=3),
                           paste("Every individual must have a unique name,",
                                 "must be associated to a unique group and",
                                 "must be associated only once to each",
                                 "of the same 2 time measurements."),
                           fixed=TRUE)

    testthat::expect_error(ColnamesToFactors(ExprData=subdataTonly.sim,
                                             Column.gene=1,
                                             Group.position=NULL,
                                             Time.position=1,
                                             Individual.position=2),
                           paste("Every individual must have a unique name and",
                                 "must be associated only once to each",
                                 "of the same 4 time measurements."),
                           fixed=TRUE)

    testthat::expect_error(ColnamesToFactors(ExprData=subdataGonly.sim,
                                             Column.gene=1,
                                             Group.position=1,
                                             Time.position=NULL,
                                             Individual.position=2),
                           paste0("Every individual must be associated to ",
                                  "only one group."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(ColnamesToFactors(ExprData=Data1smpl.sim$Sim.dat,
                                             Column.gene=1,
                                             Group.position=1,
                                             Time.position=2,
                                             Individual.position=3),
                           "Each group must have at least two individuals.",
                           fixed=TRUE)

    testthat::expect_error(ColnamesToFactors(ExprData=DataT1smpl.sim$Sim.dat,
                                             Column.gene=1,
                                             Group.position=NULL,
                                             Time.position=1,
                                             Individual.position=2),
                           paste("The data must contain the temporal",
                                 "expression of at least two individuals."),
                           fixed=TRUE)

    testthat::expect_error(ColnamesToFactors(ExprData=DataG1smpl.sim$Sim.dat,
                                             Column.gene=1,
                                             Group.position=1,
                                             Time.position=NULL,
                                             Individual.position=2),
                           "Each group must have at least two individuals.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    res1_ColTOfct <- ColnamesToFactors(ExprData=Data.sim2$Sim.dat,
                                       Column.gene=1,
                                       Group.position=1,
                                       Time.position=2,
                                       Individual.position=3)

    res2_ColTOfct <- ColnamesToFactors(ExprData=Data.sim2$Sim.dat[,-1],
                                       Column.gene=NULL,
                                       Group.position=1,
                                       Time.position=2,
                                       Individual.position=3)

    res3_ColTOfct <- ColnamesToFactors(ExprData=DataTonly.sim$Sim.dat,
                                       Column.gene=1,
                                       Group.position=NULL,
                                       Time.position=1,
                                       Individual.position=2)

    res4_ColTOfct <- ColnamesToFactors(ExprData=DataGonly.sim$Sim.dat,
                                       Column.gene=1,
                                       Group.position=1,
                                       Time.position=NULL,
                                       Individual.position=2)

    testthat::expect_vector(res1_ColTOfct$Final.Name)
    testthat::expect_vector(res2_ColTOfct$Final.Name)
    testthat::expect_vector(res3_ColTOfct$Final.Name)
    testthat::expect_vector(res4_ColTOfct$Final.Name)

})
