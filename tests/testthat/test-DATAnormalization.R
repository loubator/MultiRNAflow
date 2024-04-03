testthat::test_that("Test DATAnormalization", {
    ##
    ##-----------------------------------------------------------------------##
    data("Transcript_HomoSapiens_Database")
    data("RawCounts_Schleiss2021_CLLsub500")

    col0Leuk <- seq(2, 55, 9)
    col1Leuk <- seq(3, 55, 9)
    col2Leuk <- seq(4, 55, 9)
    col3tLeuk <- sort(c(1, col0Leuk, col1Leuk, col2Leuk))

    dataLeuk <- RawCounts_Schleiss2021_CLLsub500[seq_len(27), col3tLeuk]
    SEresDATAprepLeuk <- DATAprepSE(RawCounts=dataLeuk,
                                    Column.gene=1,
                                    Group.position=2,
                                    Time.position=4,
                                    Individual.position=3,
                                    RNAlength="hsapiens")

    SEresDATAprepLeuk2 <- DATAprepSE(RawCounts=dataLeuk,
                                     Column.gene=1,
                                     Group.position=2,
                                     Time.position=4,
                                     Individual.position=3,
                                     RNAlength=NULL)

    ##-----------------------------------------------------------------------##
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    datamus500 <- RawCounts_Antoszewski2022_MOUSEsub500
    datamus2 <- datamus500[seq_len(100),]

    SEprepMus2 <- DATAprepSE(RawCounts=datamus2[,seq_len(10)],
                             Column.gene=1,
                             Group.position=1,
                             Time.position=NULL,
                             Individual.position=2)

    resDATAprepNOid <- SEprepMus2
    S4Vectors::metadata(resDATAprepNOid)$SEidentification <- "Test"

    resDATAprepNULLid <- SEprepMus2
    NewMETA <- S4Vectors::metadata(resDATAprepNULLid)[seq_len(4)]
    S4Vectors::metadata(resDATAprepNULLid) <- NewMETA

    ##-----------------------------------------------------------------------##
    datamus500_1 <- datamus500[, seq(1, 7)]
    datamus500_2 <- datamus500[, c(1, seq(8, 13))]
    datamus500_3 <- datamus500[seq_len(351), c(1, 2, 3, 4, 11, 12, 13)]
    colnames(datamus500_3) <- colnames(datamus500_2) <- colnames(datamus500_1)

    datamus3 <- rbind(datamus500_1, datamus500_2, datamus500_3)
    datamus3$Gene <- paste0("Gene", seq_len(nrow(datamus3)))

    SEresDATAprepMus3 <- DATAprepSE(RawCounts=datamus3,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAnormalization(SEres=datamus2,
                                             Normalization="other",
                                             Plot.Boxplot=TRUE,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE,
                                             path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           paste0("'SEres' mut be the results of the function",
                                  " 'DATAprepSE()'"),
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=resDATAprepNOid,
                                             Normalization="vst",
                                             Plot.Boxplot=TRUE,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE,
                                             path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           paste0("'SEres' mut be the results of the function",
                                  " 'DATAprepSE()'"),
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=resDATAprepNULLid,
                                             Normalization="vst",
                                             Plot.Boxplot=TRUE,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE,
                                             path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           paste0("'SEres' mut be the results of the function",
                                  " 'DATAprepSE()'"),
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=SEprepMus2,
                                             Normalization="other",
                                             Plot.Boxplot=TRUE,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE, path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           paste0("'Normalization' mut be 'vst', 'rlog', 'rle'",
                                  " or 'rpkm'"),
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=SEresDATAprepLeuk2,
                                             Normalization="rpkm",
                                             Plot.Boxplot=FALSE,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE,
                                             path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           paste0("The input 'RNAlength' of the function ",
                                  "'DATAprepSE' can not be NULL."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAnormalization(SEres=SEprepMus2,
                                             Normalization="vst",
                                             Plot.Boxplot=1,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE,
                                             path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=1.1),
                           "'Blind.rlog.vst' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=SEprepMus2,
                                             Normalization="vst",
                                             Plot.Boxplot=1,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE,
                                             path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           "'Plot.Boxplot' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=SEprepMus2,
                                             Normalization="vst",
                                             Plot.Boxplot=TRUE,
                                             Colored.By.Factors=2,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE,
                                             path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           "'Colored.By.Factors' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=SEprepMus2,
                                             Normalization="vst",
                                             Plot.Boxplot=TRUE,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=3,
                                             Plot.genes=FALSE,
                                             path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           "'Color.Group' must be NULL or a data.frame.",
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=SEprepMus2,
                                             Normalization="vst",
                                             Plot.Boxplot=TRUE,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=4, path.result=NULL,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           "'Plot.genes' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=SEprepMus2,
                                             Normalization="vst",
                                             Plot.Boxplot=FALSE,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE, path.result=5,
                                             Name.folder.norm=NULL,
                                             Blind.rlog.vst=FALSE),
                           "'path.result' must be NULL or a character.",
                           fixed=TRUE)

    testthat::expect_error(DATAnormalization(SEres=SEprepMus2,
                                             Normalization="vst",
                                             Plot.Boxplot=FALSE,
                                             Colored.By.Factors=FALSE,
                                             Color.Group=NULL,
                                             Plot.genes=FALSE,
                                             path.result=NULL,
                                             Name.folder.norm=6,
                                             Blind.rlog.vst=FALSE),
                           "'Name.folder.norm' must be NULL or a character.",
                           fixed=TRUE)


    ##-----------------------------------------------------------------------##
    ##-----------------------------------------------------------------------##
    res_ErrNorm1 <- ErrDATAnormalization1(SEres=SEprepMus2,
                                          Normalization="vst",
                                          Blind.rlog.vst=FALSE)

    res_ErrNorm2 <- ErrDATAnormalization2(Plot.Boxplot=TRUE,
                                          Colored.By.Factors=FALSE,
                                          Color.Group=NULL,
                                          Plot.genes=FALSE,
                                          path.result=NULL,
                                          Name.folder.norm=NULL)

    testthat::expect_equal(res_ErrNorm1, "No error")
    testthat::expect_equal(res_ErrNorm2, "No error")

    ##-----------------------------------------------------------------------##
    res1_DATAnorm <- DATAnormalization(SEres=SEresDATAprepLeuk,
                                       Normalization="rpkm",
                                       Blind.rlog.vst=FALSE, Plot.Boxplot=FALSE,
                                       Colored.By.Factors=FALSE,
                                       Color.Group=NULL, Plot.genes=FALSE,
                                       path.result=NULL, Name.folder.norm=NULL)

    res2_DATAnorm <- DATAnormalization(SEres=SEprepMus2,
                                       Normalization="rle",
                                       Blind.rlog.vst=FALSE, Plot.Boxplot=TRUE,
                                       Colored.By.Factors=FALSE,
                                       Color.Group=NULL, Plot.genes=FALSE,
                                       path.result=NULL, Name.folder.norm=NULL)

    res3_DATAnorm <- DATAnormalization(SEres=SEprepMus2,
                                       Normalization="vst",
                                       Blind.rlog.vst=FALSE, Plot.Boxplot=TRUE,
                                       Colored.By.Factors=TRUE,
                                       Color.Group=NULL, Plot.genes=FALSE,
                                       path.result=NULL, Name.folder.norm=NULL)

    res4_DATAnorm <- DATAnormalization(SEres=SEresDATAprepMus3,
                                       Normalization="vst",
                                       Blind.rlog.vst=FALSE, Plot.Boxplot=FALSE,
                                       Colored.By.Factors=FALSE,
                                       Color.Group=NULL, Plot.genes=TRUE,
                                       path.result=NULL, Name.folder.norm=NULL)

    res5_DATAnorm <- DATAnormalization(SEres=SEprepMus2,
                                       Normalization="rlog",
                                       Blind.rlog.vst=TRUE, Plot.Boxplot=FALSE,
                                       Colored.By.Factors=FALSE,
                                       Color.Group=NULL, Plot.genes=FALSE,
                                       path.result=NULL, Name.folder.norm="Nm")

    testthat::expect_s4_class(res1_DATAnorm, "SummarizedExperiment")
    testthat::expect_s4_class(res2_DATAnorm, "SummarizedExperiment")
    testthat::expect_s4_class(res3_DATAnorm, "SummarizedExperiment")
    testthat::expect_s4_class(res4_DATAnorm, "SummarizedExperiment")
    testthat::expect_s4_class(res5_DATAnorm, "SummarizedExperiment")
})
