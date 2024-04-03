testthat::test_that("Test DEanalysisGlobal", {
    ##-----------------------------------------------------------------------##
    ## Data
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    ## No time points. We take only two groups for the speed of the example
    RawCounts_T1Wt <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200),
                                                            seq_len(7)]

    ## Preprocessing step
    resDATAprepSEmus1 <- DATAprepSE(RawCounts=RawCounts_T1Wt,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    ##-----------------------------------------------------------------------##
    data("RawCounts_Leong2014_FISSIONsub500wt")
    ## We take only the first three time for the speed of the example
    RawCounts_Fission_3t <- RawCounts_Leong2014_FISSIONsub500wt[seq_len(200),
                                                                seq_len(7)]

    ## Preprocessing step
    SEprepFission <- DATAprepSE(RawCounts=RawCounts_Fission_3t,
                                Column.gene=1,
                                Group.position=NULL,
                                Time.position=2,
                                Individual.position=3)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DEanalysisGlobal(SEres=list(1,2),
                                            pval.min=0.05,
                                            pval.vect.t=NULL,
                                            log.FC.min=1,
                                            LRT.supp.info=FALSE,
                                            Plot.DE.graph=TRUE,
                                            path.result=NULL,
                                            Name.folder.DE=NULL),
                           paste0("'SEres' mut be the results of either the ",
                                  "function 'DATAprepSE()' or ",
                                  "'DATAnormalization()'."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    resDEmus1 <- DEanalysisGlobal(SEres=resDATAprepSEmus1,
                                  pval.min=0.05,
                                  pval.vect.t=NULL,
                                  log.FC.min=1,
                                  LRT.supp.info=FALSE,
                                  Plot.DE.graph=TRUE,
                                  path.result=NULL,
                                  Name.folder.DE="test")

    resDEfission <- DEanalysisGlobal(SEres=SEprepFission,
                                     pval.min=0.05,
                                     pval.vect.t=c(0.05, 0.05, 0.05),
                                     log.FC.min=1,
                                     LRT.supp.info=TRUE,
                                     Plot.DE.graph=TRUE,
                                     path.result=NULL,
                                     Name.folder.DE=NULL)

    res1 <- Glossary(path.result=NULL, Case=1)
    res2 <- Glossary(path.result=NULL, Case=2)
    res3 <- Glossary(path.result=NULL, Case=3)

    testthat::expect_s4_class(resDEmus1, "SummarizedExperiment")
    testthat::expect_s4_class(resDEfission, "SummarizedExperiment")
    testthat::expect_type(res1, "list")
    testthat::expect_type(res2, "list")
    testthat::expect_type(res3, "list")

})
