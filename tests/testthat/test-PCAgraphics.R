testthat::test_that("Test PCAgraphics", {
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
                                       Colored.By.Factors=TRUE)

    dataCOL <- data.frame(BC=c("N1wtT1wt", "N1haT1wt"), COL=c("red", "blue"))

    ##-----------------------------------------------------------------------##
    ## data importation
    dataSIM <- RawCountsSimulation(Nb.Group=2,
                                   Nb.Time=3,
                                   Nb.per.GT=5,
                                   Nb.Gene=50)$Sim.dat

    ## preprocessing
    SEresPrepSIM <- DATAprepSE(RawCounts=dataSIM,
                               Column.gene=1,
                               Group.position=1,
                               Time.position=2,
                               Individual.position=3)

    ## normalization
    SEresNormSIM <- DATAnormalization(SEres=SEresPrepSIM,
                                      Normalization="rle",
                                      Plot.Boxplot=FALSE,
                                      Colored.By.Factors=TRUE)

    ##-----------------------------------------------------------------------##
    ## data importation
    dataSIMt <- RawCountsSimulation(Nb.Group=1,
                                    Nb.Time=3,
                                    Nb.per.GT=5,
                                    Nb.Gene=50)$Sim.dat

    ## preprocessing
    SEresPrepSIMt <- DATAprepSE(RawCounts=dataSIMt,
                                Column.gene=1,
                                Group.position=NULL,
                                Time.position=1,
                                Individual.position=2)

    ## normalization
    SEresNormSIMt <- DATAnormalization(SEres=SEresPrepSIMt,
                                       Normalization="rle",
                                       Plot.Boxplot=FALSE,
                                       Colored.By.Factors=TRUE)

    ##-----------------------------------------------------------------------##
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    testthat::expect_error(PCAgraphics(SEresNorm=datamus2,
                                       DATAnorm=TRUE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresPrepMus2,
                                       DATAnorm=TRUE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_SE,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=1.1,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           "'DATAnorm' must be TRUE or FALSE.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    Err_sdel <- paste("'sample.deletion' must be either NULL",
                      "either character or non negative integers.")
    Err_smpl <- paste0("The elements 'smpl1' are not present ",
                       "among the column names of the raw count data.")

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=-1,
                                       Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_sdel,
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=c(1.1),
                                       Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_sdel,
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion="smpl1",
                                       Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_smpl,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    Err_gdel <- paste("'gene.deletion' must be either NULL",
                      "either character or non negative integers.")
    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=-1,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_gdel,
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=c(1.1),
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_gdel,
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=data.frame(c(1,1)),
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_gdel,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=1.1,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           "'Plot.PCA' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=1.1,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           "'Mean.Accross.Time' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=1.1,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           "'motion3D' must be TRUE or FALSE.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=1.1, Name.file.pca=NULL),
                           "'path.result' must be NULL or a character.",
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=1),
                           "'Name.file.pca' must be NULL or a character.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    Err_a <- "'Phi', 'Theta' and 'epsilon' must be positive numeric values"
    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi="test", Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_a,
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=-1, Theta=140, Cex.label=0.7,
                                       Cex.point=0.7, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_a,
                           fixed=TRUE)

    Err_cex <- "'Cex.point' and 'Cex.label' must be positive numeric values"
    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point="test", epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_cex,
                           fixed=TRUE)

    testthat::expect_error(PCAgraphics(SEresNorm=SEresNormMus2,
                                       DATAnorm=FALSE,
                                       gene.deletion=NULL,
                                       sample.deletion=NULL,
                                       Plot.PCA=FALSE,
                                       Mean.Accross.Time=FALSE,
                                       Color.Group=NULL,
                                       motion3D=FALSE,
                                       Phi=25, Theta=140, Cex.label=0.7,
                                       Cex.point=-1, epsilon=0.2,
                                       path.result=NULL, Name.file.pca=NULL),
                           Err_cex,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    res1_PCAgraph <- PCAgraphics(SEresNorm=SEresNormMus2,
                                 DATAnorm=FALSE,
                                 gene.deletion=NULL,
                                 sample.deletion=NULL,
                                 Plot.PCA=TRUE,
                                 Mean.Accross.Time=FALSE,
                                 Color.Group=dataCOL,
                                 motion3D=FALSE,
                                 Phi=25, Theta=140, Cex.label=0.7,
                                 Cex.point=0.7, epsilon=0.2,
                                 path.result=NULL,
                                 Name.file.pca="Test")

    res2_PCAgraph <- PCAgraphics(SEresNorm=SEresNormSIM,
                                 DATAnorm=FALSE,
                                 gene.deletion=NULL,
                                 sample.deletion=NULL,
                                 Plot.PCA=FALSE,
                                 Mean.Accross.Time=FALSE,
                                 Color.Group=NULL,
                                 motion3D=FALSE,
                                 Phi=25, Theta=140, Cex.label=0.7,
                                 Cex.point=0.7, epsilon=0.2,
                                 path.result=NULL,
                                 Name.file.pca=NULL)

    res3_PCAgraph <- PCAgraphics(SEresNorm=SEresNormSIM,
                                 DATAnorm=FALSE,
                                 gene.deletion=NULL,
                                 sample.deletion=NULL,
                                 Plot.PCA=TRUE,
                                 Mean.Accross.Time=TRUE,
                                 Color.Group=NULL,
                                 motion3D=FALSE,
                                 Phi=25, Theta=140, Cex.label=0.7,
                                 Cex.point=0.7, epsilon=0.2,
                                 path.result=NULL,
                                 Name.file.pca=NULL)

    res4_PCAgraph <- PCAgraphics(SEresNorm=SEresNormSIMt,
                                 DATAnorm=FALSE,
                                 gene.deletion=NULL,
                                 sample.deletion=NULL,
                                 Plot.PCA=FALSE,
                                 Mean.Accross.Time=TRUE,
                                 Color.Group=NULL,
                                 motion3D=FALSE,
                                 Phi=25, Theta=140, Cex.label=0.7,
                                 Cex.point=0.7, epsilon=0.2,
                                 path.result=NULL,
                                 Name.file.pca=NULL)

    testthat::expect_s4_class(res1_PCAgraph, "SummarizedExperiment")
    testthat::expect_s4_class(res2_PCAgraph, "SummarizedExperiment")
    testthat::expect_s4_class(res3_PCAgraph, "SummarizedExperiment")
    testthat::expect_s4_class(res4_PCAgraph, "SummarizedExperiment")

})
