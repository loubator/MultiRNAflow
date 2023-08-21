test_that("Test DEanalysisSubData", {
    ##------------------------------------------------------------------------#
    expect_error(DEanalysisSubData(SEresDE=matrix(1, ncol=10, nrow=10),
                                   ColumnsCriteria=1,
                                   Set.Operation="OTHER",
                                   Save.SubData=FALSE),
                 paste0("'SEresDE' mut be the results of the function ",
                        "'DEanalysisGlobal()'."),
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    ## data importation
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    ## No time points. We take only two groups for the speed of the example
    dataT1wt <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200), seq_len(7)]

    ## Preprocessing with Results of DEanalysisGlobal()
    resDATAprepSE <- DATAprepSE(RawCounts=dataT1wt,
                                Column.gene=1,
                                Group.position=1,
                                Time.position=NULL,
                                Individual.position=2)

    ##------------------------------------------------------------------------#
    ## DE analysis
    resDET1wt <- DEanalysisGlobal(SEres=resDATAprepSE,
                                  pval.min=0.05,
                                  pval.vect.t=NULL,
                                  log.FC.min=1,
                                  LRT.supp.info=FALSE,
                                  Plot.DE.graph=FALSE,
                                  path.result=NULL,
                                  Name.folder.DE=NULL)

    ##------------------------------------------------------------------------#
    expect_error(DEanalysisSubData(SEresDE=resDET1wt,
                                   ColumnsCriteria=1.5,
                                   Set.Operation="union",
                                   Save.SubData=FALSE),
                 "'ColumnsCriteria' must be integers or characters",
                 fixed=TRUE)

    expect_error(DEanalysisSubData(SEresDE=resDET1wt,
                                   ColumnsCriteria=1,
                                   Set.Operation="OTHER",
                                   Save.SubData=FALSE),
                 "Set.Operation mut be 'union', 'intersect' or 'setdiff'",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_s4_class(DEanalysisSubData(SEresDE=resDET1wt,
                                      ColumnsCriteria=2,
                                      Set.Operation="union",
                                      Save.SubData=FALSE),
                    "SummarizedExperiment")

})
