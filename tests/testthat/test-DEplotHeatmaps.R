testthat::test_that("Test DEplotHeatmaps", {
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
    resDEmus1_4 <- resDEmus1
    S4Vectors::metadata(resDEmus1_2)$DESeq2obj$SEidentification <- NULL
    S4Vectors::metadata(resDEmus1_3)$DESeq2obj$SEidentification <- "test"

    DEsummaryMus1_4 <- SummarizedExperiment::rowData(resDEmus1_4)
    DEsummaryMus1_4$Specific.genes_N1haT1ko <- 0
    DEsummaryMus1_4$Specific.genes_N1haT1ko[c(1, 2)] <- 1
    SummarizedExperiment::rowData(resDEmus1_4) <- DEsummaryMus1_4

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
    # resHeatmap <- DEplotHeatmaps(SEresDE=resDET1wt,
    #                              ColumnsCriteria=3, ##Specific genes N1haT1ko
    #                              Set.Operation="union",
    #                              NbGene.analysis=20,
    #                              Color.Group=NULL,
    #                              SizeLabelRows=5,
    #                              SizeLabelCols=5,
    #                              Display.plots=TRUE,
    #                              Save.plots=FALSE)
    # testthat::expect_s4_class(resHeatmap$Heatmap.Correlation, "Heatmap")
    ##-----------------------------------------------------------------------##

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DEplotHeatmaps(SEresDE=matrix(0, ncol=2, nrow=3),
                                          ColumnsCriteria=2,
                                          Display.plots=TRUE),
                           paste0("'SEresDE' mut be the results of ",
                                  "the function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1_2,
                                          ColumnsCriteria=2,
                                          Display.plots=FALSE),
                           paste0("'SEresDE' mut be the results of ",
                                  "the function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1_3,
                                          ColumnsCriteria=2,
                                          Display.plots=TRUE),
                           paste0("'SEresDE' mut be the results of ",
                                  "the function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1_4,
                                          ColumnsCriteria=3,
                                          Display.plots=TRUE),
                           paste("Only 2 genes have been selected.",
                                 "The function will work",
                                 "with at least three genes."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    vNNImessage <- paste("'ColumnsCriteria' must be a non-negative integer,",
                         "a vector of non-negative integers or 'NULL'")

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1,
                                          ColumnsCriteria=c(-1, 0),
                                          Set.Operation="union",
                                          NbGene.analysis=10,
                                          Color.Group=dataCOL,
                                          SizeLabelRows=5,
                                          SizeLabelCols=5,
                                          Display.plots=TRUE,
                                          Save.plots=TRUE),
                           vNNImessage,
                           fixed=TRUE)

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1,
                                          ColumnsCriteria="columnCr",
                                          Set.Operation="union",
                                          NbGene.analysis=10,
                                          Color.Group=dataCOL,
                                          SizeLabelRows=5,
                                          SizeLabelCols=5,
                                          Display.plots=TRUE,
                                          Save.plots=TRUE),
                           vNNImessage,
                           fixed=TRUE)

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1,
                                          ColumnsCriteria=3,
                                          Set.Operation="test",
                                          NbGene.analysis=10,
                                          Color.Group=dataCOL,
                                          SizeLabelRows=5,
                                          SizeLabelCols=5,
                                          Display.plots=TRUE,
                                          Save.plots=TRUE),
                           paste("'Set.Operation' mut be 'union', 'intersect'",
                                 "or 'setdiff'"),
                           fixed=TRUE)

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1,
                                          ColumnsCriteria=3,
                                          Set.Operation="union",
                                          NbGene.analysis=1.1,
                                          Color.Group=dataCOL,
                                          SizeLabelRows=5,
                                          SizeLabelCols=5,
                                          Display.plots=TRUE,
                                          Save.plots=TRUE),
                           "'NbGene.analysis' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1,
                                          ColumnsCriteria=3,
                                          Set.Operation="union",
                                          NbGene.analysis=5,
                                          Color.Group="test",
                                          SizeLabelRows=5,
                                          SizeLabelCols=5,
                                          Display.plots=TRUE,
                                          Save.plots=TRUE),
                           "'Color.Group' must be NULL or a data.frame.",
                           fixed=TRUE)

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1,
                                          ColumnsCriteria=3,
                                          Set.Operation="union",
                                          NbGene.analysis=5,
                                          Color.Group=NULL,
                                          SizeLabelRows="test",
                                          SizeLabelCols=5,
                                          Display.plots=TRUE,
                                          Save.plots=TRUE),
                           paste("'SizeLabelRows' and 'SizeLabelCols' must be",
                                 "strictly positive numeric values"),
                           fixed=TRUE)

    testthat::expect_error(DEplotHeatmaps(SEresDE=resDEmus1,
                                          ColumnsCriteria=3,
                                          Set.Operation="union",
                                          NbGene.analysis=5,
                                          Color.Group=NULL,
                                          SizeLabelRows=5,
                                          SizeLabelCols=-1,
                                          Display.plots=TRUE,
                                          Save.plots=TRUE),
                           paste("'SizeLabelRows' and 'SizeLabelCols' must be",
                                 "strictly positive numeric values"),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    ## Specific genes N1haT1ko
    dataCOL <- data.frame(BC=c("N1wtT1wt", "N1haT1wt", "N1haT1ko", "N1wtT1ko"),
                          COL=c("red", "blue", "green", "black"))
    dataCOLleuk <- data.frame(BC=c("NP", "P"), COL=c("black", "red"))

    res1_Heatmaps <- DEplotHeatmaps(SEresDE=resDEmus1,
                                    ColumnsCriteria=seq(3, 6, 1),
                                    Set.Operation="union",
                                    NbGene.analysis=10,
                                    Color.Group=dataCOL,
                                    SizeLabelRows=5,
                                    SizeLabelCols=5,
                                    Display.plots=TRUE,
                                    Save.plots=TRUE)

    res2_Heatmaps <- DEplotHeatmaps(SEresDE=resDEmus1,
                                    ColumnsCriteria=NULL,
                                    Set.Operation="union",
                                    NbGene.analysis=10,
                                    Color.Group=NULL,
                                    SizeLabelRows=5,
                                    SizeLabelCols=5,
                                    Display.plots=FALSE,
                                    Save.plots=FALSE)

    res3_Heatmaps <- DEplotHeatmaps(SEresDE=resDET1wt,
                                    ColumnsCriteria=2,
                                    Set.Operation="union",
                                    NbGene.analysis=10,
                                    Color.Group=NULL,
                                    SizeLabelRows=5,
                                    SizeLabelCols=5,
                                    Display.plots=FALSE,
                                    Save.plots=FALSE)

    res4_Heatmaps <- DEplotHeatmaps(SEresDE=resDEleuk,
                                    ColumnsCriteria=c(2, 3),
                                    Set.Operation="union",
                                    NbGene.analysis=10,
                                    Color.Group=NULL,
                                    SizeLabelRows=5,
                                    SizeLabelCols=5,
                                    Display.plots=FALSE,
                                    Save.plots=FALSE)

    res5_Heatmaps <- DEplotHeatmaps(SEresDE=resDEleuk,
                                    ColumnsCriteria=c(2, 3),
                                    Set.Operation="union",
                                    NbGene.analysis=10,
                                    Color.Group=dataCOLleuk,
                                    SizeLabelRows=5,
                                    SizeLabelCols=5,
                                    Display.plots=FALSE,
                                    Save.plots=FALSE)

    testthat::expect_s4_class(res1_Heatmaps, "SummarizedExperiment")
    testthat::expect_s4_class(res2_Heatmaps, "SummarizedExperiment")
    testthat::expect_s4_class(res3_Heatmaps, "SummarizedExperiment")
    testthat::expect_s4_class(res4_Heatmaps, "SummarizedExperiment")
    testthat::expect_s4_class(res5_Heatmaps, "SummarizedExperiment")

})

