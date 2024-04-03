testthat::test_that("Test DATAplotExpression1Gene", {
    ##
    ##-----------------------------------------------------------------------##
    ## Simulation raw counts
    resSIMcountTG <- RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=4,
                                         Nb.Gene=10)
    # Preprocessing step
    SEprepTGsim <- DATAprepSE(RawCounts=resSIMcountTG$Sim.dat,
                              Column.gene=1, Group.position=1,
                              Time.position=2, Individual.position=3)

    resDATAprepNOid <- SEprepTGsim
    S4Vectors::metadata(resDATAprepNOid)$SEidentification <- "Test"

    resDATAprepNULLid <- SEprepTGsim
    NewMETA <- S4Vectors::metadata(resDATAprepNULLid)[seq_len(4)]
    S4Vectors::metadata(resDATAprepNULLid) <- NewMETA

    ##-----------------------------------------------------------------------##
    ## Simulation raw counts
    resSIMcountT <- RawCountsSimulation(Nb.Group=1, Nb.Time=3, Nb.per.GT=4,
                                        Nb.Gene=10)
    # Preprocessing step
    SEprepTsim <- DATAprepSE(RawCounts=resSIMcountT$Sim.dat,
                             Column.gene=1, Group.position=NULL,
                             Time.position=1, Individual.position=2)

    ## Simulation raw counts
    resSIMcountG <- RawCountsSimulation(Nb.Group=3, Nb.Time=1, Nb.per.GT=4,
                                        Nb.Gene=10)
    # Preprocessing step
    SEprepGsim <- DATAprepSE(RawCounts=resSIMcountG$Sim.dat,
                             Column.gene=1, Group.position=1,
                             Time.position=NULL, Individual.position=2)

    ## Simulation raw counts
    resSIMcountG2 <- RawCountsSimulation(Nb.Group=3, Nb.Time=1, Nb.per.GT=51,
                                         Nb.Gene=10)
    # Preprocessing step
    SEprepGsim2 <- DATAprepSE(RawCounts=resSIMcountG2$Sim.dat,
                              Column.gene=1, Group.position=1,
                              Time.position=NULL, Individual.position=2)

    dataCOL <- data.frame(BC=c("G1", "G2", "G3"), COL=c("red", "blue", "black"))

    ##-----------------------------------------------------------------------##
    Err_SE <- paste0("'SEres' mut be the results of either the function ",
                     "'DATAprepSE()' or 'DATAnormalization()'.")

    testthat::expect_error(DATAplotExpression1Gene(SEres=resSIMcountTG$Sim.dat,
                                                   row.gene=1,
                                                   Color.Group=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpression1Gene(SEres=resDATAprepNOid,
                                                   row.gene=1,
                                                   Color.Group=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpression1Gene(SEres=resDATAprepNULLid,
                                                   row.gene=1,
                                                   Color.Group=NULL),
                           Err_SE,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAplotExpression1Gene(SEres=SEprepTGsim,
                                                   row.gene=1.5,
                                                   Color.Group=NULL),
                           "'row.gene' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpression1Gene(SEres=SEprepTGsim,
                                                   row.gene=-1,
                                                   Color.Group=NULL),
                           "'row.gene' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpression1Gene(SEres=SEprepTGsim,
                                                   row.gene="Try",
                                                   Color.Group=NULL),
                           "'row.gene' must be a non-negative integer.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAplotExpression1Gene(SEres=SEprepTGsim,
                                                   row.gene=1,
                                                   Color.Group=1.5),
                           "'Color.Group' must be NULL or a data.frame.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    res1_profile1gene <- DATAplotExpression1Gene(SEres=SEprepTGsim,
                                                 row.gene=1, Color.Group=NULL)

    res2_profile1gene <- DATAplotExpression1Gene(SEres=SEprepTsim,
                                                 row.gene=1, Color.Group=NULL)

    res3_profile1gene <- DATAplotExpression1Gene(SEres=SEprepGsim,
                                                 row.gene=1,
                                                 Color.Group=dataCOL)

    res4_profile1gene <- DATAplotExpression1Gene(SEres=SEprepGsim2,
                                                 row.gene=1,
                                                 Color.Group=NULL)

    testthat::expect_s3_class(res1_profile1gene, "ggplot")
    testthat::expect_s3_class(res2_profile1gene, "ggplot")
    testthat::expect_s3_class(res3_profile1gene, "ggplot")
    testthat::expect_s3_class(res4_profile1gene, "ggplot")
})
