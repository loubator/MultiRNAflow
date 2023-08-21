test_that("Test RawCountsSimulation", {
    expect_error(RawCountsSimulation(Nb.Group=1,
                                     Nb.Time=1,
                                     Nb.per.GT=3,
                                     Nb.Gene=10),
                 "At least two groups or two times are demanded.",
                 fixed=TRUE)

    expect_error(RawCountsSimulation(Nb.Group=2.5,
                                     Nb.Time=1,
                                     Nb.per.GT=3,
                                     Nb.Gene=10),
                 "'Nb.Group' must be an integer.",
                 fixed=TRUE)

    expect_error(RawCountsSimulation(Nb.Group=1,
                                     Nb.Time=2.5,
                                     Nb.per.GT=3,
                                     Nb.Gene=10),
                 "'Nb.Time' must be an integer.",
                 fixed=TRUE)

    expect_error(RawCountsSimulation(Nb.Group=1,
                                     Nb.Time=3,
                                     Nb.per.GT=4.5,
                                     Nb.Gene=10),
                 "'Nb.per.GT' must be an integer.",
                 fixed=TRUE)

    expect_error(RawCountsSimulation(Nb.Group=2,
                                     Nb.Time=3,
                                     Nb.per.GT=3,
                                     Nb.Gene=10.5),
                 "'Nb.Gene' must be an integer.",
                 fixed=TRUE)
})
