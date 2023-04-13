test_that("Test DEanalysisSubData", {
    expect_error(DEanalysisSubData(Data=matrix(1, ncol=10, nrow=10),
                                   Res.DE.analysis=matrix(1, ncol=10, nrow=10),
                                   ColumnsCriteria=1,
                                   Set.Operation="OTHER",
                                   Save.SubData=FALSE),
                 "Set.Operation mut be 'union', 'intersect' or 'setdiff'",
                 fixed=TRUE)
    #
    expect_error(DEanalysisSubData(Data=matrix(1, ncol=10, nrow=10),
                                   Res.DE.analysis=matrix(1, ncol=10, nrow=10),
                                   ColumnsCriteria=1.5,
                                   Set.Operation="union",
                                   Save.SubData=FALSE),
                 "'ColumnsCriteria' must be integers or characters",
                 fixed=TRUE)
    #

    Stop.nrow<-paste("The number of rows of 'Data' and, ",
                     "the length or the number of rows of 'Res.DE.analysis', ",
                     "must be identical.", sep= "")

    expect_error(DEanalysisSubData(Data=matrix(1, ncol=10, nrow=10),
                                   Res.DE.analysis=matrix(1, ncol=10, nrow=20),
                                   ColumnsCriteria=2,
                                   Set.Operation="union",
                                   Save.SubData=FALSE),
                 Stop.nrow, fixed=TRUE)
})
