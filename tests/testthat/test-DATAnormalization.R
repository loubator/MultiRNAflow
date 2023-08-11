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

    ##------------------------------------------------------------------------#
    expect_error(DATAnormalization(SEres=datamus2,
                                   Normalization="other", Plot.Boxplot=TRUE,
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
})
