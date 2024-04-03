testthat::test_that("Test PCArealization", {
    ##
    ##-----------------------------------------------------------------------##
    ## data importation
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(50), seq_len(7)]

    ## preprocessing
    SEresPrepMus2 <- DATAprepSE(RawCounts=datamus2,
                                Column.gene=1,
                                Group.position=1,
                                Time.position=NULL,
                                Individual.position=2)

    ## normalization
    SEresNormMus2 <- DATAnormalization(SEres=SEresPrepMus2,
                                       Normalization="rle",
                                       Plot.Boxplot=FALSE,
                                       Colored.By.Factors=FALSE)

    ## preprocessing
    SEresPrepMus2_c <- DATAprepSE(RawCounts=datamus2[,-1],
                                  Column.gene=NULL,
                                  Group.position=1,
                                  Time.position=NULL,
                                  Individual.position=2)

    ## normalization
    SEresNormMus2_c <- DATAnormalization(SEres=SEresPrepMus2_c,
                                         Normalization="rle",
                                         Plot.Boxplot=FALSE,
                                         Colored.By.Factors=FALSE)

    ##-----------------------------------------------------------------------##
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    testthat::expect_error(PCArealization(SEresNorm=datamus2,
                                          DATAnorm=TRUE,
                                          gene.deletion=c("Gene3", "Gene5"),
                                          sample.deletion=c(3, 8),
                                          Supp.del.sample=TRUE),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(PCArealization(SEresNorm=SEresPrepMus2,
                                          DATAnorm=TRUE,
                                          gene.deletion=c("Gene3", "Gene5"),
                                          sample.deletion=c(3, 8),
                                          Supp.del.sample=TRUE),
                           Err_SE,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=1.1,
                                          gene.deletion=NULL,
                                          sample.deletion=NULL,
                                          Supp.del.sample=FALSE),
                           "'DATAnorm' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=TRUE,
                                          gene.deletion=NULL,
                                          sample.deletion=NULL,
                                          Supp.del.sample=1.1),
                           "'Supp.del.sample' must be TRUE or FALSE.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    Err_sdel <- paste("'sample.deletion' must be either NULL",
                      "either character or non negative integers.")
    Err_smpl <- paste0("The elements 'smpl1', 'smpl2' are not present ",
                       "among the column names of the raw count data.")

    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=TRUE,
                                          gene.deletion=NULL,
                                          sample.deletion=-1,
                                          Supp.del.sample=FALSE),
                           Err_sdel,
                           fixed=TRUE)

    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=FALSE,
                                          gene.deletion=NULL,
                                          sample.deletion=c(1.1),
                                          Supp.del.sample=FALSE),
                           Err_sdel,
                           fixed=TRUE)

    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=TRUE,
                                          gene.deletion=c(1, 2),
                                          sample.deletion=data.frame(c(1)),
                                          Supp.del.sample=TRUE),
                           Err_sdel,
                           fixed=TRUE)

    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=TRUE,
                                          gene.deletion=c(1, 2),
                                          sample.deletion=c("smpl1", "smpl2"),
                                          Supp.del.sample=TRUE),
                           Err_smpl,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    Err_gdel <- paste("'gene.deletion' must be either NULL",
                      "either character or non negative integers.")

    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=FALSE,
                                          gene.deletion=-1,
                                          sample.deletion=NULL,
                                          Supp.del.sample=FALSE),
                           Err_gdel,
                           fixed=TRUE)

    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=FALSE,
                                          gene.deletion=c(1.1),
                                          sample.deletion=NULL,
                                          Supp.del.sample=FALSE),
                           Err_gdel,
                           fixed=TRUE)

    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=FALSE,
                                          gene.deletion=data.frame(c(1,1)),
                                          sample.deletion=NULL,
                                          Supp.del.sample=FALSE),
                           Err_gdel,
                           fixed=TRUE)

    testthat::expect_error(PCArealization(SEresNorm=SEresNormMus2,
                                          DATAnorm=TRUE,
                                          gene.deletion=c("Gene3", "Gene5"),
                                          sample.deletion=c(3, 8),
                                          Supp.del.sample=FALSE),
                           paste0("The genes 1, 2 selected are not present ",
                                  "in your dataset."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    genDel <- c("ENSMUSG00000025777", "ENSMUSG00000066877")

    res1_PCAres <- PCArealization(SEresNorm=SEresNormMus2,
                                  DATAnorm=TRUE,
                                  gene.deletion=c(1, 2),
                                  sample.deletion=c(3, 8),
                                  Supp.del.sample=FALSE)

    res2_PCAres <- PCArealization(SEresNorm=SEresNormMus2,
                                  DATAnorm=TRUE,
                                  gene.deletion=genDel,##c(1, 2),
                                  sample.deletion="N1wtT1wt_r2",
                                  Supp.del.sample=TRUE)

    res3_PCAres <- PCArealization(SEresNorm=SEresNormMus2_c,
                                  DATAnorm=TRUE,
                                  gene.deletion=NULL,
                                  sample.deletion="N1wtT1wt_r2",
                                  Supp.del.sample=TRUE)

    res4_PCAres <- PCArealization(SEresNorm=SEresNormMus2_c,
                                  DATAnorm=TRUE,
                                  gene.deletion=NULL,
                                  sample.deletion=c(3, 8),
                                  Supp.del.sample=TRUE)

    testthat::expect_s4_class(res1_PCAres, "SummarizedExperiment")
    testthat::expect_s4_class(res2_PCAres, "SummarizedExperiment")
    testthat::expect_s4_class(res3_PCAres, "SummarizedExperiment")
    testthat::expect_s4_class(res4_PCAres, "SummarizedExperiment")

})
