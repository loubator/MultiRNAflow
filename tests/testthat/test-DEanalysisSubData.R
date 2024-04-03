testthat::test_that("Test DEanalysisSubData", {
    ##-----------------------------------------------------------------------##
    data("Results_DEanalysis_sub500")
    resDEmus1 <- Results_DEanalysis_sub500$DE_Antoszewski2022_MOUSEsub500

    resDEmus1_2 <- resDEmus1
    resDEmus1_3 <- resDEmus1
    S4Vectors::metadata(resDEmus1_2)$DESeq2obj$SEidentification <- NULL
    S4Vectors::metadata(resDEmus1_3)$DESeq2obj$SEidentification <- "test"

    SummarizedExperiment::rowData(resDEmus1)[,6] <- 1

    ## resDEleuk <- Results_DEanalysis_sub500$DE_Schleiss2021_CLLsub500
    ## resDET1wt <- Results_DEanalysis_sub500$DE_Leong2014_FISSIONsub500wt


    ##-----------------------------------------------------------------------##
    ## ## the  previous lines is equivalent to the following uncomment lines
    ## ## data importation
    ## data("RawCounts_Antoszewski2022_MOUSEsub500")
    ## ## No time points. We take only two groups for the speed of the example
    ## dataMus1<-RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200),seq_len(7)]
    ##
    ## ## Preprocessing with Results of DEanalysisGlobal()
    ## SEprepMus1 <- DATAprepSE(RawCounts=dataMus1,
    ##                          Column.gene=1,
    ##                          Group.position=1,
    ##                          Time.position=NULL,
    ##                          Individual.position=2)
    ## ## DE analysis
    ## SEmus1DE <- DEanalysisGlobal(SEres=SEprepMus1,
    ##                              pval.min=0.05,
    ##                              pval.vect.t=NULL,
    ##                              log.FC.min=1,
    ##                              LRT.supp.info=FALSE,
    ##                              Plot.DE.graph=FALSE,
    ##                              path.result=NULL,
    ##                              Name.folder.DE=NULL)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DEanalysisSubData(SEresDE=matrix(1, ncol=9, nrow=9),
                                             ColumnsCriteria=1,
                                             Set.Operation="union",
                                             Save.SubData=FALSE),
                           paste0("'SEresDE' mut be the results of the ",
                                  "function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(DEanalysisSubData(SEresDE=resDEmus1_2,
                                             ColumnsCriteria=1,
                                             Set.Operation="union",
                                             Save.SubData=FALSE),
                           paste0("'SEresDE' mut be the results of the ",
                                  "function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(DEanalysisSubData(SEresDE=resDEmus1_3,
                                             ColumnsCriteria=1,
                                             Set.Operation="union",
                                             Save.SubData=FALSE),
                           paste0("'SEresDE' mut be the results of the ",
                                  "function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    colDataframe <- data.frame(matrix(1, ncol=2, nrow=2))
    testthat::expect_error(DEanalysisSubData(SEresDE=resDEmus1,
                                             ColumnsCriteria=colDataframe,
                                             Set.Operation="union",
                                             Save.SubData=FALSE),
                           "'ColumnsCriteria' must be integers or characters",
                           fixed=TRUE)

    testthat::expect_error(DEanalysisSubData(SEresDE=resDEmus1,
                                             ColumnsCriteria=1.5,
                                             Set.Operation="union",
                                             Save.SubData=FALSE),
                           "'ColumnsCriteria' must be integers or characters",
                           fixed=TRUE)

    testthat::expect_error(DEanalysisSubData(SEresDE=resDEmus1,
                                             ColumnsCriteria="test",
                                             Set.Operation="union",
                                             Save.SubData=FALSE),
                           "The element 1 of 'ColumnsCriteria' is not correct.",
                           fixed=TRUE)

    testthat::expect_error(DEanalysisSubData(SEresDE=resDEmus1,
                                             ColumnsCriteria=50,
                                             Set.Operation="union",
                                             Save.SubData=FALSE),
                           paste("Integers of 'ColumnsCriteria' must be",
                                 "between 1 and 28"),
                           fixed=TRUE)

    testthat::expect_error(DEanalysisSubData(SEresDE=resDEmus1,
                                             ColumnsCriteria=1,
                                             Set.Operation="OTHER",
                                             Save.SubData=FALSE),
                           paste("Set.Operation mut be 'union', 'intersect'",
                                 "or 'setdiff'"),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    selColnames <- c("Specific.genes_N1haT1ko", "Specific.genes_N1haT1wt")

    testthat::expect_s4_class(DEanalysisSubData(SEresDE=resDEmus1,
                                                ColumnsCriteria=2,
                                                Set.Operation="union",
                                                Save.SubData=FALSE),
                              "SummarizedExperiment")

    testthat::expect_output(DEanalysisSubData(SEresDE=resDEmus1,
                                              ColumnsCriteria=6,
                                              Set.Operation="union",
                                              Save.SubData=FALSE),
                            paste("All rows are selected because there is",
                                  "no 0 in the column selected.",
                                  "The original 'SEresDE' is returned."))

    testthat::expect_output(DEanalysisSubData(SEresDE=resDEmus1,
                                              ColumnsCriteria=c(3, 4),
                                              Set.Operation="intersect",
                                              Save.SubData=FALSE),
                            paste("No selection because the column selected",
                                  "is full of 0.",
                                  "The original 'SEresDE' is returned."))

    testthat::expect_s4_class(DEanalysisSubData(SEresDE=resDEmus1,
                                                ColumnsCriteria=seq(3, 5, 1),
                                                Set.Operation="union",
                                                Save.SubData=FALSE),
                              "SummarizedExperiment")

    testthat::expect_s4_class(DEanalysisSubData(SEresDE=resDEmus1,
                                                ColumnsCriteria=c(2, 3),
                                                Set.Operation="intersect",
                                                Save.SubData=TRUE),
                              "SummarizedExperiment")

    testthat::expect_s4_class(DEanalysisSubData(SEresDE=resDEmus1,
                                                ColumnsCriteria=selColnames,
                                                Set.Operation="setdiff",
                                                Save.SubData=FALSE),
                              "SummarizedExperiment")

})
