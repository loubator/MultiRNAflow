testthat::test_that("Test DEplotVolcanoMA", {
    ##-----------------------------------------------------------------------##
    data("Results_DEanalysis_sub500")
    resDEmus1 <- Results_DEanalysis_sub500$DE_Antoszewski2022_MOUSEsub500
    resDEleuk <- Results_DEanalysis_sub500$DE_Schleiss2021_CLLsub500
    resDET1wt <- Results_DEanalysis_sub500$DE_Leong2014_FISSIONsub500wt

    DEsummaryMus1 <- SummarizedExperiment::rowData(resDEmus1)
    DEsummaryMus1$Pvalue.adjusted_N1haT1wt.versus.N1haT1ko.[c(3, 4)] <- 0
    SummarizedExperiment::rowData(resDEmus1) <- DEsummaryMus1

    resDEmus1_2 <- resDEmus1
    resDEmus1_3 <- resDEmus1
    S4Vectors::metadata(resDEmus1_2)$DESeq2obj$SEidentification <- NULL
    S4Vectors::metadata(resDEmus1_3)$DESeq2obj$SEidentification <- "test"

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DEplotVolcanoMA(SEresDE=matrix(0, nrow=3, ncol=2),
                                           NbGene.plotted=2,
                                           SizeLabel=3,
                                           Display.plots=FALSE,
                                           Save.plots=FALSE),
                           paste0("'SEresDE' mut be the results of ",
                                  "the function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(DEplotVolcanoMA(SEresDE=resDEmus1_2,
                                           NbGene.plotted=2,
                                           SizeLabel=3,
                                           Display.plots=FALSE,
                                           Save.plots=FALSE),
                           paste0("'SEresDE' mut be the results of ",
                                  "the function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(DEplotVolcanoMA(SEresDE=resDEmus1_3,
                                           NbGene.plotted=2,
                                           SizeLabel=3,
                                           Display.plots=FALSE,
                                           Save.plots=FALSE),
                           paste0("'SEresDE' mut be the results of ",
                                  "the function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    ## ##the two previous lines is equivalent to the following uncomment lines
    ## ## data importation
    ## data(RawCounts_Antoszewski2022_MOUSEsub500)
    ## ## No time points. We take only two groups for the speed of the example
    ## dataMus1<-RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200),seq_len(7)]
    ##
    ## ## Preprocessing with Results of DEanalysisGlobal()
    ## resDATAprepSE <- DATAprepSE(RawCounts=dataMus1,
    ##                             Column.gene=1,
    ##                             Group.position=1,
    ##                             Time.position=NULL,
    ##                             Individual.position=2)
    ## ## DE analysis
    ## resDEmus1 <- DEanalysisGlobal(SEres=resDATAprepSE,
    ##                               pval.min=0.05,
    ##                               pval.vect.t=NULL,
    ##                               log.FC.min=1,
    ##                               LRT.supp.info=FALSE,
    ##                               Plot.DE.graph=FALSE,
    ##                               path.result=NULL,
    ##                               Name.folder.DE=NULL)
    ##-----------------------------------------------------------------------##
    ## ## Volcano MA
    ## resVolcanoMA <- DEplotVolcanoMA(SEresDE=resDEmus1,
    ##                                 NbGene.plotted=5,
    ##                                 Display.plots=TRUE,
    ##                                 Save.plots=FALSE)
    ## expect_s3_class(resVolcanoMA$List.plot.Volcano[[1]], "ggplot")
    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    testthat::expect_error(DEplotVolcanoMA(SEresDE=resDEmus1,
                                           NbGene.plotted=1.1,
                                           SizeLabel=3,
                                           Display.plots=FALSE,
                                           Save.plots=FALSE),
                           "'NbGene.plotted' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(DEplotVolcanoMA(SEresDE=resDEmus1,
                                           NbGene.plotted=5,
                                           SizeLabel=0,
                                           Display.plots=FALSE,
                                           Save.plots=FALSE),
                           paste("'SizeLabel' must be strictly positive",
                                 "numeric value"),
                           fixed=TRUE)

    testthat::expect_error(DEplotVolcanoMA(SEresDE=resDEmus1,
                                           NbGene.plotted=5,
                                           SizeLabel="test",
                                           Display.plots=FALSE,
                                           Save.plots=FALSE),
                           paste("'SizeLabel' must be strictly positive",
                                 "numeric value"),
                           fixed=TRUE)

    testthat::expect_error(DEplotVolcanoMA(SEresDE=resDEmus1,
                                           NbGene.plotted=5,
                                           SizeLabel=3,
                                           Display.plots="test",
                                           Save.plots=FALSE),
                           "'Display.plots' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DEplotVolcanoMA(SEresDE=resDEmus1,
                                           NbGene.plotted=5,
                                           SizeLabel=3,
                                           Display.plots=FALSE,
                                           Save.plots="test"),
                           "'Save.plots' must be TRUE or FALSE.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    res1_MAvolcano <- DEplotVolcanoMA(SEresDE=resDEmus1, NbGene.plotted=5,
                                      Display.plots=TRUE, Save.plots=FALSE)

    res2_MAvolcano <- DEplotVolcanoMA(SEresDE=resDEleuk, NbGene.plotted=5,
                                      Display.plots=FALSE, Save.plots=FALSE)

    res3_MAvolcano <- DEplotVolcanoMA(SEresDE=resDEleuk, NbGene.plotted=5,
                                      Display.plots=FALSE, Save.plots=TRUE)

    res4_MAvolcano <- DEplotVolcanoMA(SEresDE=resDET1wt, NbGene.plotted=5,
                                      Display.plots=FALSE, Save.plots=FALSE)

    testthat::expect_s4_class(res1_MAvolcano, "SummarizedExperiment")
    testthat::expect_s4_class(res2_MAvolcano, "SummarizedExperiment")
    testthat::expect_s4_class(res3_MAvolcano, "SummarizedExperiment")
    testthat::expect_s4_class(res4_MAvolcano, "SummarizedExperiment")

})
