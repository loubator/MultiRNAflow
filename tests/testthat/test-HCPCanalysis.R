testthat::test_that("Test HCPCanalysis", {
    ##
    ##-----------------------------------------------------------------------##
    ## data importation
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(50), seq_len(7)]
    dataCOL <- data.frame(BC=c("N1wtT1wt", "N1haT1wt"), COL=c("red", "blue"))

    ## preprocessing
    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    ## normalization
    resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
                                         Normalization="rle",
                                         Plot.Boxplot=TRUE,
                                         Colored.By.Factors=TRUE)

    ##------------------------------------------------------------------------#
    data("RawCounts_Leong2014_FISSIONsub500wt")
    ## We take only the first three time for the speed of the example
    RawCounts_Fission_3t <- RawCounts_Leong2014_FISSIONsub500wt[seq_len(30),
                                                                seq_len(10)]

    ## Preprocessing step
    resDATAprepSEfission <- DATAprepSE(RawCounts=RawCounts_Fission_3t,
                                       Column.gene=1,
                                       Group.position=NULL,
                                       Time.position=2,
                                       Individual.position=3)

    ## normalization
    resDATAnormSEfission <- DATAnormalization(SEres=resDATAprepSEfission,
                                              Normalization="rle",
                                              Plot.Boxplot=TRUE,
                                              Colored.By.Factors=TRUE)

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

    ##-----------------------------------------------------------------------##
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    testthat::expect_error(HCPCanalysis(SEresNorm=datamus2,
                                        DATAnorm=TRUE,
                                        gene.deletion=NULL,
                                        sample.deletion=NULL,
                                        Plot.HCPC=FALSE,
                                        Color.Group=NULL,
                                        motion3D=FALSE,
                                        Phi=25, Theta=140, Cex.label=0.7,
                                        Cex.point=0.7, epsilon=0.2,
                                        path.result=NULL,
                                        Name.folder.hcpc=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(HCPCanalysis(SEresNorm=resDATAprepSEmus2,
                                        DATAnorm=TRUE,
                                        gene.deletion=NULL,
                                        sample.deletion=NULL,
                                        Plot.HCPC=FALSE,
                                        Color.Group=NULL,
                                        motion3D=FALSE,
                                        Phi=25, Theta=140, Cex.label=0.7,
                                        Cex.point=0.7, epsilon=0.2,
                                        path.result=NULL,
                                        Name.folder.hcpc=NULL),
                           Err_SE,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                                        DATAnorm=1.1,
                                        gene.deletion=NULL,
                                        sample.deletion=NULL,
                                        Plot.HCPC=FALSE,
                                        Color.Group=NULL,
                                        motion3D=FALSE,
                                        Phi=25, Theta=140, Cex.label=0.7,
                                        Cex.point=0.7, epsilon=0.2,
                                        path.result=NULL,
                                        Name.folder.hcpc=NULL),
                           "'DATAnorm' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                                        DATAnorm=FALSE,
                                        gene.deletion=NULL,
                                        sample.deletion=NULL,
                                        Plot.HCPC=1.1,
                                        Color.Group=NULL,
                                        motion3D=FALSE,
                                        Phi=25, Theta=140, Cex.label=0.7,
                                        Cex.point=0.7, epsilon=0.2,
                                        path.result=NULL,
                                        Name.folder.hcpc=NULL),
                           "'Plot.HCPC' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                                        DATAnorm=FALSE,
                                        gene.deletion=NULL,
                                        sample.deletion=NULL,
                                        Plot.HCPC=FALSE,
                                        Color.Group=NULL,
                                        motion3D=1.1,
                                        Phi=25, Theta=140, Cex.label=0.7,
                                        Cex.point=0.7, epsilon=0.2,
                                        path.result=NULL,
                                        Name.folder.hcpc=NULL),
                           "'motion3D' must be TRUE or FALSE.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                                        DATAnorm=FALSE,
                                        gene.deletion=NULL,
                                        sample.deletion=NULL,
                                        Plot.HCPC=FALSE,
                                        Color.Group=NULL,
                                        motion3D=FALSE,
                                        Phi=25, Theta=140, Cex.label=0.7,
                                        Cex.point=0.7, epsilon=0.2,
                                        path.result=1.1,
                                        Name.folder.hcpc=NULL),
                           "'path.result' must be NULL or a character.",
                           fixed=TRUE)

    testthat::expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                                        DATAnorm=FALSE,
                                        gene.deletion=NULL,
                                        sample.deletion=NULL,
                                        Plot.HCPC=FALSE,
                                        Color.Group=NULL,
                                        motion3D=FALSE,
                                        Phi=25, Theta=140, Cex.label=0.7,
                                        Cex.point=0.7, epsilon=0.2,
                                        path.result=NULL, Name.folder.hcpc=1),
                           "'Name.folder.hcpc' must be NULL or a character.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    color1_HCPC <- myPaletteHCPC(Nclust=3)
    color2_HCPC <- myPaletteHCPC(Nclust=18)

    testthat::expect_vector(color1_HCPC)
    testthat::expect_vector(color2_HCPC)

    ##-----------------------------------------------------------------------##
    res1_HCPC <- HCPCanalysis(SEresNorm=resDATAnormMus2, DATAnorm=FALSE,
                              gene.deletion=NULL, sample.deletion=NULL,
                              motion3D=FALSE, Plot.HCPC=FALSE,
                              Color.Group=dataCOL, Phi=25, Theta=140,
                              Cex.label=0.7, Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc=NULL)

    res2_HCPC <- HCPCanalysis(SEresNorm=resDATAnormSim, DATAnorm=FALSE,
                              gene.deletion=NULL, sample.deletion=NULL,
                              motion3D=FALSE, Plot.HCPC=TRUE, Color.Group=NULL,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc="test")

    res3_HCPC <- HCPCanalysis(SEresNorm=resDATAnormSEfission, DATAnorm=FALSE,
                              gene.deletion=NULL, sample.deletion=NULL,
                              motion3D=FALSE, Plot.HCPC=TRUE, Color.Group=NULL,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc="test")

    testthat::expect_s4_class(res1_HCPC, "SummarizedExperiment")
    testthat::expect_s4_class(res2_HCPC, "SummarizedExperiment")
    testthat::expect_s4_class(res3_HCPC, "SummarizedExperiment")

})
