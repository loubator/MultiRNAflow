testthat::test_that("Test DATAplotExpressionGenes", {
    ##
    ##-----------------------------------------------------------------------##
    ## data importation
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(50), seq_len(7)]

    dataCOL <- data.frame(BC=c("N1wtT1wt", "N1haT1wt"), COL=c("red", "blue"))

    ## preprocessing
    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    ## normalization
    SEresNormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
                                       Normalization="rle",
                                       Plot.Boxplot=FALSE,
                                       Colored.By.Factors=FALSE)

    resDATAprepNULLid <- SEresNormMus2
    NewMETA <- S4Vectors::metadata(resDATAprepNULLid)[seq_len(4)]
    S4Vectors::metadata(resDATAprepNULLid) <- NewMETA

    ##-----------------------------------------------------------------------##
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=datamus2,
                                                   Vector.row.gene=c(1, 3),
                                                   DATAnorm=TRUE,
                                                   Color.Group=NULL,
                                                   Plot.Expression=FALSE,
                                                   path.result=NULL,
                                                   Name.folder.profile=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=resDATAprepSEmus2,
                                                   Vector.row.gene=c(1, 3),
                                                   DATAnorm=TRUE,
                                                   Color.Group=NULL,
                                                   Plot.Expression=FALSE,
                                                   path.result=NULL,
                                                   Name.folder.profile=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=resDATAprepNULLid,
                                                   Vector.row.gene=c(1, 3),
                                                   DATAnorm=TRUE,
                                                   Color.Group=NULL,
                                                   Plot.Expression=FALSE,
                                                   path.result=NULL,
                                                   Name.folder.profile=NULL),
                           Err_SE,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    Err_integers <- paste("'Vector.row.gene' must be a vector",
                          "of non negative integers.")

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                   Vector.row.gene=c("test",
                                                                     "try"),
                                                   DATAnorm=TRUE,
                                                   Color.Group=NULL,
                                                   Plot.Expression=FALSE,
                                                   path.result=NULL,
                                                   Name.folder.profile=NULL),
                           Err_integers,
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                   Vector.row.gene=c(1.1, 3),
                                                   DATAnorm=TRUE,
                                                   Color.Group=NULL,
                                                   Plot.Expression=FALSE,
                                                   path.result=NULL,
                                                   Name.folder.profile=NULL),
                           Err_integers,
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                   Vector.row.gene=c(1, -3),
                                                   DATAnorm=TRUE,
                                                   Color.Group=NULL,
                                                   Plot.Expression=FALSE,
                                                   path.result=NULL,
                                                   Name.folder.profile=NULL),
                           Err_integers,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                   Vector.row.gene=c(1, 3),
                                                   DATAnorm=1,
                                                   Color.Group=NULL,
                                                   Plot.Expression=FALSE,
                                                   path.result=NULL,
                                                   Name.folder.profile=NULL),
                           "'DATAnorm' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                   Vector.row.gene=c(1, 3),
                                                   DATAnorm=TRUE,
                                                   Color.Group=1,
                                                   Plot.Expression=FALSE,
                                                   path.result=NULL,
                                                   Name.folder.profile=NULL),
                           "'Color.Group' must be NULL or a data.frame.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                   Vector.row.gene=c(1, 3),
                                                   DATAnorm=TRUE,
                                                   Color.Group=NULL,
                                                   Plot.Expression=1,
                                                   path.result=NULL,
                                                   Name.folder.profile=NULL),
                           "'Plot.Expression' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                   Vector.row.gene=c(1, 3),
                                                   DATAnorm=TRUE,
                                                   Color.Group=NULL,
                                                   Plot.Expression=FALSE,
                                                   path.result=2,
                                                   Name.folder.profile=NULL),
                           "'path.result' must be NULL or a character.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                   Vector.row.gene=c(1, 3),
                                                   DATAnorm=TRUE,
                                                   Color.Group=NULL,
                                                   Plot.Expression=FALSE,
                                                   path.result=NULL,
                                                   Name.folder.profile=1),
                           "'Name.folder.profile' must be NULL or a character.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    res1_Errprofile <- ErrSEresNorm(SEresNorm=SEresNormMus2,
                                    DATAnorm=TRUE, path.result=NULL)
    res2_Errprofile <- ErrProfileGenes(Vector.row.gene=c(1, 3),
                                       Color.Group=dataCOL,
                                       Plot.Expression=FALSE,
                                       Name.folder.profile=NULL)

    testthat::expect_equal(res1_Errprofile, "No error")
    testthat::expect_equal(res2_Errprofile, "No error")


    res1_Ynorm <- myYlabelNorm("vst")
    res2_Ynorm <- myYlabelNorm("rlog")
    res3_Ynorm <- myYlabelNorm("rle")
    res4_Ynorm <- myYlabelNorm("rleRPKM")

    testthat::expect_equal(res1_Ynorm, "vst normalized counts")
    testthat::expect_equal(res2_Ynorm, "rlog normalized counts")
    testthat::expect_equal(res3_Ynorm, "rle normalized counts")
    testthat::expect_equal(res4_Ynorm, "rle and RPKM normalized counts")

    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    res1_profileGenes <- DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                 Vector.row.gene=c(1, 3),
                                                 DATAnorm=TRUE,
                                                 Color.Group=dataCOL,
                                                 Plot.Expression=FALSE,
                                                 path.result=NULL,
                                                 Name.folder.profile=NULL)

    res2_profileGenes <- DATAplotExpressionGenes(SEresNorm=SEresNormMus2,
                                                 Vector.row.gene=c(3),
                                                 DATAnorm=FALSE,
                                                 Color.Group=NULL,
                                                 Plot.Expression=TRUE,
                                                 path.result=NULL,
                                                 Name.folder.profile="tt")

    testthat::expect_s4_class(res1_profileGenes, "SummarizedExperiment")
    testthat::expect_s4_class(res2_profileGenes, "SummarizedExperiment")

})
