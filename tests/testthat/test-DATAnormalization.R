test_that("Test DATAnormalization", {
    ##
    ##------------------------------------------------------------------------#
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2<-RawCounts_Antoszewski2022_MOUSEsub500

    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    resDATAprepNOid <- resDATAprepSEmus2
    S4Vectors::metadata(resDATAprepNOid)$SEidentification <- "Test"

    resDATAprepNULLid <- resDATAprepSEmus2
    NewMETA <- S4Vectors::metadata(resDATAprepNULLid)[seq_len(4)]
    S4Vectors::metadata(resDATAprepNULLid) <- NewMETA

    ##------------------------------------------------------------------------#
    expect_error(DATAnormalization(SEres=datamus2,
                                   Normalization="other", Plot.Boxplot=TRUE,
                                   Colored.By.Factors=FALSE, Color.Group=NULL,
                                   Plot.genes=FALSE, path.result=NULL,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=FALSE),
                 "'SEres' mut be the results of the function 'DATAprepSE()'",
                 fixed=TRUE)

    expect_error(DATAnormalization(SEres=resDATAprepNOid,
                                   Normalization="vst", Plot.Boxplot=TRUE,
                                   Colored.By.Factors=FALSE, Color.Group=NULL,
                                   Plot.genes=FALSE, path.result=NULL,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=FALSE),
                 "'SEres' mut be the results of the function 'DATAprepSE()'",
                 fixed=TRUE)

    expect_error(DATAnormalization(SEres=resDATAprepNULLid,
                                   Normalization="vst", Plot.Boxplot=TRUE,
                                   Colored.By.Factors=FALSE, Color.Group=NULL,
                                   Plot.genes=FALSE, path.result=NULL,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=FALSE),
                 "'SEres' mut be the results of the function 'DATAprepSE()'",
                 fixed=TRUE)

    expect_error(DATAnormalization(SEres=resDATAprepSEmus2,
                                   Normalization="other", Plot.Boxplot=TRUE,
                                   Colored.By.Factors=FALSE, Color.Group=NULL,
                                   Plot.genes=FALSE, path.result=NULL,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=FALSE),
                 "Normalization mut be 'vst', 'rlog' or 'rle'",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(DATAnormalization(SEres=resDATAprepSEmus2,
                                   Normalization="vst",
                                   Plot.Boxplot=1,
                                   Colored.By.Factors=FALSE, Color.Group=NULL,
                                   Plot.genes=FALSE, path.result=NULL,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=1.1),
                 "'Blind.rlog.vst' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(DATAnormalization(SEres=resDATAprepSEmus2,
                                   Normalization="vst", Plot.Boxplot=1,
                                   Colored.By.Factors=FALSE, Color.Group=NULL,
                                   Plot.genes=FALSE, path.result=NULL,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=FALSE),
                 "'Plot.Boxplot' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(DATAnormalization(SEres=resDATAprepSEmus2,
                                   Normalization="vst", Plot.Boxplot=TRUE,
                                   Colored.By.Factors=2, Color.Group=NULL,
                                   Plot.genes=FALSE, path.result=NULL,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=FALSE),
                 "'Colored.By.Factors' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(DATAnormalization(SEres=resDATAprepSEmus2,
                                   Normalization="vst", Plot.Boxplot=TRUE,
                                   Colored.By.Factors=FALSE, Color.Group=3,
                                   Plot.genes=FALSE, path.result=NULL,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=FALSE),
                 "'Color.Group' must be NULL or a data.frame.",
                 fixed=TRUE)

    expect_error(DATAnormalization(SEres=resDATAprepSEmus2,
                                   Normalization="vst", Plot.Boxplot=TRUE,
                                   Colored.By.Factors=FALSE, Color.Group=NULL,
                                   Plot.genes=4, path.result=NULL,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=FALSE),
                 "'Plot.genes' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(DATAnormalization(SEres=resDATAprepSEmus2,
                                   Normalization="vst", Plot.Boxplot=FALSE,
                                   Colored.By.Factors=FALSE, Color.Group=NULL,
                                   Plot.genes=FALSE, path.result=5,
                                   Name.folder.norm=NULL,
                                   Blind.rlog.vst=FALSE),
                 "'path.result' must be NULL or a character.",
                 fixed=TRUE)

    expect_error(DATAnormalization(SEres=resDATAprepSEmus2,
                                   Normalization="vst", Plot.Boxplot=FALSE,
                                   Colored.By.Factors=FALSE, Color.Group=NULL,
                                   Plot.genes=FALSE, path.result=NULL,
                                   Name.folder.norm=6,
                                   Blind.rlog.vst=FALSE),
                 "'Name.folder.norm' must be NULL or a character.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_s4_class(DATAnormalization(SEres=resDATAprepSEmus2,
                                      Normalization="vst", Plot.Boxplot=FALSE,
                                      Colored.By.Factors=FALSE,
                                      Color.Group=NULL,
                                      Plot.genes=FALSE, path.result=NULL,
                                      Name.folder.norm=NULL,
                                      Blind.rlog.vst=FALSE),
                    "SummarizedExperiment")

    expect_s4_class(DATAnormalization(SEres=resDATAprepSEmus2,
                                      Normalization="vst",
                                      Plot.Boxplot=TRUE,
                                      Colored.By.Factors=FALSE,
                                      Color.Group=NULL,
                                      Plot.genes=FALSE, path.result=NULL,
                                      Name.folder.norm=NULL,
                                      Blind.rlog.vst=FALSE),
                    "SummarizedExperiment")

    expect_s4_class(DATAnormalization(SEres=resDATAprepSEmus2,
                                      Normalization="vst",
                                      Plot.Boxplot=TRUE,
                                      Colored.By.Factors=TRUE,
                                      Color.Group=NULL,
                                      Plot.genes=FALSE, path.result=NULL,
                                      Name.folder.norm=NULL,
                                      Blind.rlog.vst=FALSE),
                    "SummarizedExperiment")

    expect_s4_class(DATAnormalization(SEres=resDATAprepSEmus2,
                                      Normalization="vst",
                                      Plot.Boxplot=TRUE,
                                      Colored.By.Factors=FALSE,
                                      Color.Group=NULL,
                                      Plot.genes=TRUE, path.result=NULL,
                                      Name.folder.norm=NULL,
                                      Blind.rlog.vst=FALSE),
                    "SummarizedExperiment")

    expect_s4_class(DATAnormalization(SEres=resDATAprepSEmus2,
                                      Normalization="rlog",
                                      Plot.Boxplot=FALSE,
                                      Colored.By.Factors=FALSE,
                                      Color.Group=NULL,
                                      Plot.genes=FALSE, path.result=NULL,
                                      Name.folder.norm="Name",
                                      Blind.rlog.vst=TRUE),
                    "SummarizedExperiment")

})
