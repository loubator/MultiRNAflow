testthat::test_that("Test PCAanalysis", {
    ##
    ##-----------------------------------------------------------------------##
    ## data importation
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(50), seq_len(10)]

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

    resDATAnormMus2Id <- resDATAnormMus2
    S4Vectors::metadata(resDATAnormMus2Id)$SEidentification <- "test"

    resDATAnormMus2Id2 <- resDATAnormMus2
    S4Vectors::metadata(resDATAnormMus2Id2)$SEidentification <- NULL
    ##-----------------------------------------------------------------------##
    ## data importation
    dataSIM <- RawCountsSimulation(Nb.Group=2,
                                   Nb.Time=3,
                                   Nb.per.GT=5,
                                   Nb.Gene=30)$Sim.dat

    ## preprocessing
    resDATAprepSEsim <- DATAprepSE(RawCounts=dataSIM,
                                   Column.gene=1,
                                   Group.position=1,
                                   Time.position=2,
                                   Individual.position=3)

    ## normalization
    resDATAnormSim <- DATAnormalization(SEres=resDATAprepSEsim,
                                        Normalization="rle",
                                        Plot.Boxplot=FALSE,
                                        Colored.By.Factors=TRUE)

    ##------------------------------------------------------------------------#
    ## ## data importation
    ## dataSIMt <- RawCountsSimulation(Nb.Group=1,
    ##                                 Nb.Time=3,
    ##                                 Nb.per.GT=5,
    ##                                 Nb.Gene=30)$Sim.dat
    ##
    ## ## preprocessing
    ## resDATAprepSEsimt <- DATAprepSE(RawCounts=dataSIMt,
    ##                                 Column.gene=1,
    ##                                 Group.position=NULL,
    ##                                 Time.position=1,
    ##                                 Individual.position=2)
    ##
    ## ## normalization
    ## resDATAnormSimt <- DATAnormalization(SEres=resDATAprepSEsimt,
    ##                                      Normalization="rle",
    ##                                      Plot.Boxplot=FALSE,
    ##                                      Colored.By.Factors=TRUE)

    ##-----------------------------------------------------------------------##
    COLnamesSim <- colnames(dataSIM)
    COLnamesSim2 <- COLnamesSim[-1]

    res1_delFUNsample <- delFUNsample(sample.deletion=2,
                                      RAWcolnames=COLnamesSim,
                                      Column.gene=1)
    res2_delFUNsample <- delFUNsample(sample.deletion="G1_t0_Ind1",
                                      RAWcolnames=COLnamesSim,
                                      Column.gene=1)
    res3_delFUNsample <- delFUNsample(sample.deletion=1,
                                      RAWcolnames=COLnamesSim2,
                                      Column.gene=NULL)
    res4_delFUNsample <- delFUNsample(sample.deletion="G1_t0_Ind1",
                                      RAWcolnames=COLnamesSim2,
                                      Column.gene=NULL)

    testthat::expect_equal(res1_delFUNsample, 2)
    testthat::expect_equal(res2_delFUNsample, 2)
    testthat::expect_equal(res3_delFUNsample, 1)
    testthat::expect_equal(res4_delFUNsample, 1)

    ##-----------------------------------------------------------------------##
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    testthat::expect_error(PCAanalysis(SEresNorm=datamus2,
                                       DATAnorm=TRUE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL,
                                       Name.folder.pca=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(PCAanalysis(SEresNorm=resDATAprepSEmus2,
                                       DATAnorm=TRUE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL,
                                       Name.folder.pca=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(PCAanalysis(SEresNorm=resDATAnormMus2Id,
                                       DATAnorm=TRUE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL,
                                       Name.folder.pca=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(PCAanalysis(SEresNorm=resDATAnormMus2Id2,
                                       DATAnorm=TRUE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL,
                                       Name.folder.pca=NULL),
                           Err_SE,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                                       DATAnorm=1.1,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL,
                                       Name.folder.pca=NULL),
                           "'DATAnorm' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=1.1,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL,
                                       Name.folder.pca=NULL),
                           "'Plot.PCA' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=1.1,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL,
                                       Name.folder.pca=NULL),
                           "'Mean.Accross.Time' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=1.1,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL,
                                       Name.folder.pca=NULL),
                           "'motion3D' must be TRUE or FALSE.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=1.1,
                                       Name.folder.pca=NULL),
                           "'path.result' must be NULL or a character.",
                           fixed=TRUE)

    testthat::expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.folder.pca=1),
                           "'Name.folder.pca' must be NULL or a character.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    res1_PCAanalysis <- PCAanalysis(SEresNorm=resDATAnormMus2,
                                    DATAnorm=FALSE,
                                    gene.deletion=c(2, 4),
                                    sample.deletion=c(2, 6),
                                    Plot.PCA=FALSE,
                                    Mean.Accross.Time=FALSE,
                                    Color.Group=NULL,
                                    motion3D=FALSE,
                                    Phi=25, Theta=140, Cex.label=0.7,
                                    Cex.point=0.7, epsilon=0.2,
                                    path.result=NULL,
                                    Name.folder.pca=NULL)

    res2_PCAanalysis <- PCAanalysis(SEresNorm=resDATAnormSim,
                                    DATAnorm=TRUE,
                                    gene.deletion=NULL,
                                    sample.deletion=c("G1_t0_Ind2",
                                                      "G1_t1_Ind3"),
                                    Plot.PCA=FALSE,
                                    Mean.Accross.Time=FALSE,
                                    Color.Group=NULL,
                                    motion3D=FALSE,
                                    Phi=25, Theta=140, Cex.label=0.7,
                                    Cex.point=0.7, epsilon=0.2,
                                    path.result=NULL,
                                    Name.folder.pca="test")

    res3_PCAanalysis <- PCAanalysis(SEresNorm=resDATAnormSim,
                                    DATAnorm=TRUE,
                                    gene.deletion=NULL,
                                    sample.deletion=NULL,
                                    Plot.PCA=FALSE,
                                    Mean.Accross.Time=TRUE,
                                    Color.Group=NULL,
                                    motion3D=FALSE,
                                    Phi=25, Theta=140, Cex.label=0.7,
                                    Cex.point=0.7, epsilon=0.2,
                                    path.result=NULL,
                                    Name.folder.pca="test")

    testthat::expect_s4_class(res1_PCAanalysis, "SummarizedExperiment")
    testthat::expect_s4_class(res2_PCAanalysis, "SummarizedExperiment")
    testthat::expect_s4_class(res3_PCAanalysis, "SummarizedExperiment")

    ## res4_PCAanalysis <- PCAanalysis(SEresNorm=resDATAnormSimt,
    ##                                 DATAnorm=FALSE,
    ##                                 gene.deletion=NULL,
    ##                                 sample.deletion=NULL,
    ##                                 Plot.PCA=FALSE,
    ##                                 Mean.Accross.Time=TRUE,
    ##                                 Color.Group=NULL,
    ##                                 motion3D=FALSE,
    ##                                 Phi=25, Theta=140, Cex.label=0.7,
    ##                                 Cex.point=0.7, epsilon=0.2,
    ##                                 path.result=NULL,
    ##                                 Name.folder.pca=NULL)
    ## testthat::expect_s4_class(res4_PCAanalysis, "SummarizedExperiment")

})
