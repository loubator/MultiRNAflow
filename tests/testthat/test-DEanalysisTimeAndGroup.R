testthat::test_that("Test DEanalysisTimeAndGroup", {
    ##------------------------------------------------------------------------#
    data("RawCounts_Schleiss2021_CLLsub500")
    ## We take only the first three times (both group) for the speed of
    ## the example
    Index3t <- c(2:4,11:13,20:22, 29:31,38:40,47:49)
    RawCounts_3t <- RawCounts_Schleiss2021_CLLsub500[seq_len(200), c(1,Index3t)]

    ## Preprocessing step
    resDATAprepSEleuk <- DATAprepSE(RawCounts=RawCounts_3t,
                                    Column.gene=1,
                                    Group.position=2,
                                    Time.position=4,
                                    Individual.position=3)

    DESeq2preprocess <- S4Vectors::metadata(resDATAprepSEleuk)$DESeq2obj
    DESeq2obj <- DESeq2preprocess$DESeq2preproceesing
    dds.DE <- DESeq2::DESeq(DESeq2obj, test="LRT", reduced=~1)

    ##------------------------------------------------------------------------#
    data("RawCounts_Weger2021_MOUSEsub500")
    ## We take only the first three times (both group) for the speed of
    ## the example
    colnamesMus2 <- colnames(RawCounts_Weger2021_MOUSEsub500)
    BC3 <- c(grep("BmKo", colnamesMus2, fixed=TRUE),
             grep("BmWt", colnamesMus2, fixed=TRUE),
             grep("CrKo", colnamesMus2, fixed=TRUE))
    Time3 <- c(grep("t0", colnamesMus2, fixed=TRUE),
               grep("t1", colnamesMus2, fixed=TRUE),
               grep("t2", colnamesMus2, fixed=TRUE))
    Index3tMus2 <- intersect(BC3, Time3)

    Mus2_3t <- RawCounts_Weger2021_MOUSEsub500[seq_len(200), c(1, Index3tMus2)]

    ## Preprocessing step
    resDATAprepSEmus2 <- DATAprepSE(RawCounts=Mus2_3t,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=2,
                                    Individual.position=3)

    DESeq2preprocessMus2 <- S4Vectors::metadata(resDATAprepSEmus2)$DESeq2obj
    DESeq2objMus2 <- DESeq2preprocessMus2$DESeq2preproceesing
    dds.DE.mus2 <- DESeq2::DESeq(DESeq2objMus2)

    ##------------------------------------------------------------------------#
    testthat::expect_error(DEanalysisTimeAndGroup(DESeq.result=list(1, 2),
                                                  pval.min=0.05,
                                                  pval.vect.t=NULL,
                                                  log.FC.min=1,
                                                  LRT.supp.info=FALSE,
                                                  Plot.DE.graph=TRUE,
                                                  path.result=NULL,
                                                  SubFile.name=NULL),
                           "Res.DE.analysis must be a 'DESeqDataSet' object",
                           fixed=TRUE)

    ##------------------------------------------------------------------------#
    testthat::expect_s4_class(DEanalysisTimeAndGroup(DESeq.result=dds.DE,
                                                     LRT.supp.info=TRUE,
                                                     pval.min=0.05,
                                                     pval.vect.t=NULL,
                                                     log.FC.min=0.1,
                                                     Plot.DE.graph=TRUE,
                                                     path.result=NULL,
                                                     SubFile.name="test"),
                              "DESeqDataSet")

    testthat::expect_s4_class(DEanalysisTimeAndGroup(DESeq.result=dds.DE.mus2,
                                                     LRT.supp.info=FALSE,
                                                     pval.min=0.2,
                                                     pval.vect.t=c(0.3),
                                                     log.FC.min=0.1,
                                                     Plot.DE.graph=TRUE,
                                                     path.result=NULL,
                                                     SubFile.name=NULL),
                              "DESeqDataSet")
})
