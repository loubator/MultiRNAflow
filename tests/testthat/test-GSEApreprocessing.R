testthat::test_that("Test GSEApreprocessing", {
    ##-----------------------------------------------------------------------##
    data("Results_DEanalysis_sub500")
    resDEmus1 <- Results_DEanalysis_sub500$DE_Antoszewski2022_MOUSEsub500

    resDEmus1_2 <- resDEmus1
    resDEmus1_3 <- resDEmus1
    S4Vectors::metadata(resDEmus1_2)$DESeq2obj$SEidentification <- NULL
    S4Vectors::metadata(resDEmus1_3)$DESeq2obj$SEidentification <- "test"

    resDEleuk <- Results_DEanalysis_sub500$DE_Schleiss2021_CLLsub500
    ## resDET1wt <- Results_DEanalysis_sub500$DE_Leong2014_FISSIONsub500wt


    ##-----------------------------------------------------------------------##
    ## ##the two previous lines is equivalent to the following uncomment lines
    ## ## data importation
    ## data(RawCounts_Antoszewski2022_MOUSEsub500)
    ## ## No time points. We take only two groups for the speed of the example
    ## dataT1wt<-RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200),seq_len(7)]
    ##
    ## ## Preprocessing with Results of DEanalysisGlobal()
    ## resDATAprepSE <- DATAprepSE(RawCounts=dataT1wt,
    ##                             Column.gene=1,
    ##                             Group.position=1,
    ##                             Time.position=NULL,
    ##                             Individual.position=2)
    ## ## DE analysis
    ## resDET1wt <- DEanalysisGlobal(SEres=resDATAprepSE,
    ##                               pval.min=0.05,
    ##                               pval.vect.t=NULL,
    ##                               log.FC.min=1,
    ##                               LRT.supp.info=FALSE,
    ##                               Plot.DE.graph=FALSE,
    ##                               path.result=NULL,
    ##                               Name.folder.DE=NULL)
    ##-----------------------------------------------------------------------##

    ##-----------------------------------------------------------------------##
    testthat::expect_error(GSEApreprocessing(SEresDE=matrix(0,ncol=2,nrow=3),
                                             ColumnsCriteria=2,
                                             Set.Operation="union",
                                             Rnk.files=TRUE,
                                             Save.files=FALSE),
                           paste0("'SEresDE' mut be the results of the ",
                                  "function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(GSEApreprocessing(SEresDE=resDEmus1_2,
                                             ColumnsCriteria=2,
                                             Set.Operation="union",
                                             Rnk.files=TRUE,
                                             Save.files=FALSE),
                           paste0("'SEresDE' mut be the results of the ",
                                  "function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(GSEApreprocessing(SEresDE=resDEmus1_3,
                                             ColumnsCriteria=2,
                                             Set.Operation="union",
                                             Rnk.files=TRUE,
                                             Save.files=FALSE),
                           paste0("'SEresDE' mut be the results of the ",
                                  "function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(GSEApreprocessing(SEresDE=resDEmus1,
                                             ColumnsCriteria=seq(3, 6, 1),
                                             Set.Operation="union",
                                             Rnk.files="Rnk.files",
                                             Save.files=FALSE),
                           "'Rnk.files' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(GSEApreprocessing(SEresDE=resDEmus1,
                                             ColumnsCriteria=seq(3, 6, 1),
                                             Set.Operation="union",
                                             Rnk.files=TRUE,
                                             Save.files="save"),
                           "'Save.files' must be TRUE or FALSE.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_type(GSEApreprocessing(SEresDE=resDEmus1,
                                            ColumnsCriteria=seq(3, 6, 1),
                                            Set.Operation="union",
                                            Rnk.files=TRUE,
                                            Save.files=FALSE)$TargetDEGene,
                          "character")

    testthat::expect_type(GSEApreprocessing(SEresDE=resDEmus1,
                                            ColumnsCriteria=seq(3, 6, 1),
                                            Set.Operation="union",
                                            Rnk.files=FALSE,
                                            Save.files=TRUE)$TargetDEGene,
                          "character")

    testthat::expect_type(GSEApreprocessing(SEresDE=resDEleuk,
                                            ColumnsCriteria=c(2, 3),
                                            Set.Operation="union",
                                            Rnk.files=TRUE,
                                            Save.files=FALSE)$TargetDEGene,
                          "character")

    ## testthat::expect_type(GSEApreprocessing(SEresDE=resDET1wt,
    ##                                         ColumnsCriteria=2,
    ##                                         Set.Operation="union",
    ##                                         Rnk.files=TRUE,
    ##                                         Save.files=FALSE)$TargetDEGene,
    ##                       "character")
    ##

})
