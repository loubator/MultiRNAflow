test_that("Test DEanalysisTimeAndGroup", {
    expect_error(DEanalysisTimeAndGroup(DESeq.result=list(1, 2),
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
    ##------------------------------------------------------------------------#
    data(RawCounts_Schleiss2021_CLLsub500)
    ## We take only the first three times (both group) for the speed of
    ## the example
    Index3t<-c(2:4,11:13,20:22, 29:31,38:40,47:49)
    RawCounts_3t<-RawCounts_Schleiss2021_CLLsub500[seq_len(200), c(1,Index3t)]

    ## Preprocessing step
    resDATAprepSEleuk <- DATAprepSE(RawCounts=RawCounts_3t,
                                    Column.gene=1,
                                    Group.position=2,
                                    Time.position=4,
                                    Individual.position=3)

    DESeq2preprocess <- S4Vectors::metadata(resDATAprepSEleuk)$DESeq2obj
    DESeq2obj <- DESeq2preprocess$DESeq2preproceesing

    ##------------------------------------------------------------------------#
    dds.DE<-DESeq2::DESeq(DESeq2obj)

    expect_s4_class(DEanalysisTimeAndGroup(DESeq.result=dds.DE,
                                           LRT.supp.info=FALSE,
                                           pval.min=0.05,
                                           pval.vect.t=NULL,
                                           log.FC.min=0.1,
                                           Plot.DE.graph=TRUE,
                                           path.result=NULL,
                                           SubFile.name=NULL),
                    "DESeqDataSet")
})
