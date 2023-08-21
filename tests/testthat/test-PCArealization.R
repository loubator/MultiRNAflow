test_that("Test PCArealization", {
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

    ##------------------------------------------------------------------------#
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    expect_error(PCArealization(SEresNorm=datamus2,
                                DATAnorm=TRUE,
                                gene.deletion=c("Gene3", "Gene5"),
                                sample.deletion=c(3, 8),
                                Supp.del.sample=TRUE),
                 Err_SE,
                 fixed=TRUE)

    expect_error(PCArealization(SEresNorm=resDATAprepSEmus2,
                                DATAnorm=TRUE,
                                gene.deletion=c("Gene3", "Gene5"),
                                sample.deletion=c(3, 8),
                                Supp.del.sample=TRUE),
                 Err_SE,
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(PCArealization(SEresNorm=resDATAnormMus2,
                                DATAnorm=TRUE,
                                gene.deletion=c("Gene3", "Gene5"),
                                sample.deletion=c(3, 8),
                                Supp.del.sample=1.1),
                 "'Supp.del.sample' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(PCArealization(SEresNorm=resDATAnormMus2,
                                DATAnorm=TRUE,
                                gene.deletion=c("Gene3", "Gene5"),
                                sample.deletion=c(3, 8),
                                Supp.del.sample=FALSE),
                 "The genes 1,2 selected are not present in your dataset.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_s4_class(PCArealization(SEresNorm=resDATAnormMus2,
                                   DATAnorm=TRUE,
                                   gene.deletion=NULL,
                                   sample.deletion=NULL,
                                   Supp.del.sample=FALSE),
                    "SummarizedExperiment")

    expect_s4_class(PCArealization(SEresNorm=resDATAnormMus2,
                                   DATAnorm=TRUE,
                                   gene.deletion=c(1,2),
                                   sample.deletion=c(3, 8),
                                   Supp.del.sample=FALSE),
                    "SummarizedExperiment")

    expect_s4_class(PCArealization(SEresNorm=resDATAnormMus2,
                                   DATAnorm=TRUE,
                                   gene.deletion=c("ENSMUSG00000025777",
                                                   "ENSMUSG00000066877"),
                                   sample.deletion=c(3, 8),
                                   Supp.del.sample=TRUE),
                    "SummarizedExperiment")
})
