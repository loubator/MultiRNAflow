test_that("multiplication works", {
    ##
    ##------------------------------------------------------------------------#
    ## data importation
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(100),]

    ## preprocessing
    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    ## normalization
    resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
                                         Normalization="rle",
                                         Plot.Boxplot=FALSE,
                                         Colored.By.Factors=FALSE)

    resDATAprepNULLid <- resDATAnormMus2
    NewMETA <- S4Vectors::metadata(resDATAprepNULLid)[seq_len(4)]
    S4Vectors::metadata(resDATAprepNULLid) <- NewMETA

    ##------------------------------------------------------------------------#
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    expect_error(DATAplotExpressionGenes(SEresNorm=datamus2,
                                         Vector.row.gene=c(1, 3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=FALSE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAprepSEmus2,
                                         Vector.row.gene=c(1, 3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=FALSE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAprepNULLid,
                                         Vector.row.gene=c(1, 3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=FALSE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 Err_SE,
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    Err_integers <- paste("'Vector.row.gene' must be a vector",
                          "of non negative integers.")

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                         Vector.row.gene=c("test", "try"),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=FALSE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 Err_integers,
                 fixed=TRUE)

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                         Vector.row.gene=c(1.1, 3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=FALSE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 Err_integers,
                 fixed=TRUE)

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                         Vector.row.gene=c(1, -3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=FALSE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 Err_integers,
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                         Vector.row.gene=c(1, 3),
                                         DATAnorm=1,
                                         Color.Group=NULL,
                                         Plot.Expression=FALSE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 "'DATAnorm' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                         Vector.row.gene=c(1, 3),
                                         DATAnorm=TRUE,
                                         Color.Group=1,
                                         Plot.Expression=FALSE,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 "'Color.Group' must be NULL or a data.frame.",
                 fixed=TRUE)

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                         Vector.row.gene=c(1, 3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=1,
                                         path.result=NULL,
                                         Name.folder.profile=NULL),
                 "'Plot.Expression' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                         Vector.row.gene=c(1, 3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=FALSE,
                                         path.result=2,
                                         Name.folder.profile=NULL),
                 "'path.result' must be NULL or a character.",
                 fixed=TRUE)

    expect_error(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                         Vector.row.gene=c(1, 3),
                                         DATAnorm=TRUE,
                                         Color.Group=NULL,
                                         Plot.Expression=FALSE,
                                         path.result=NULL,
                                         Name.folder.profile=1),
                 "'Name.folder.profile' must be NULL or a character.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_s4_class(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                            Vector.row.gene=c(1, 3),
                                            DATAnorm=TRUE,
                                            Color.Group=NULL,
                                            Plot.Expression=FALSE,
                                            path.result=NULL,
                                            Name.folder.profile=NULL),
                    "SummarizedExperiment")

    expect_s4_class(DATAplotExpressionGenes(SEresNorm=resDATAnormMus2,
                                            Vector.row.gene=c(1, 3),
                                            DATAnorm=FALSE,
                                            Color.Group=NULL,
                                            Plot.Expression=TRUE,
                                            path.result=NULL,
                                            Name.folder.profile="test"),
                    "SummarizedExperiment")

})
