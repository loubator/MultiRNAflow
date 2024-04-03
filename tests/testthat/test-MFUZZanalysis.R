testthat::test_that("Test MFUZZanalysis", {
    ##
    ##-----------------------------------------------------------------------##
    ## data importation
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(100),]
    ## preprocessing
    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)
    resDATAnormMus2 <- resDATAprepSEmus2
    resDATAnormMus2$SEidentification <- "SEresNormalization"
    ## normalization
    resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
                                         Normalization="rle",
                                         Plot.Boxplot=FALSE,
                                         Colored.By.Factors=FALSE)

    ##-----------------------------------------------------------------------##
    ## data importation ## preprocessing
    data("RawCounts_Leong2014_FISSIONsub500wt")
    datafission <- RawCounts_Leong2014_FISSIONsub500wt[seq_len(50),]

    resDATAprepSEfission <- DATAprepSE(RawCounts=datafission,
                                       Column.gene=1,
                                       Group.position=NULL,
                                       Time.position=2,
                                       Individual.position=3)

    ## normalization
    resDATAnormFission <- DATAnormalization(SEres=resDATAprepSEfission,
                                            Normalization="rle",
                                            Plot.Boxplot=FALSE,
                                            Colored.By.Factors=FALSE)

    resDATAnormFissionId <- resDATAnormFission
    S4Vectors::metadata(resDATAnormFissionId)$SEidentification <- "test"

    resDATAnormFissionId2 <- resDATAnormFission
    S4Vectors::metadata(resDATAnormFissionId2)$SEidentification <- NULL
    ##-----------------------------------------------------------------------##
    ## data importation
    dataSIM <- RawCountsSimulation(Nb.Group=2,
                                   Nb.Time=3,
                                   Nb.per.GT=5,
                                   Nb.Gene=50)$Sim.dat

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

    testthat::expect_error(MFUZZanalysis(SEresNorm=datafission,
                                         DATAnorm=TRUE,
                                         Method="hcpc",
                                         Max.clust=4,
                                         path.result=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAprepSEfission,
                                         DATAnorm=TRUE,
                                         Method="hcpc",
                                         Max.clust=4,
                                         path.result=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormFissionId,
                                         DATAnorm=TRUE,
                                         Method="hcpc",
                                         Max.clust=4,
                                         path.result=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormFissionId2,
                                         DATAnorm=TRUE,
                                         Method="hcpc",
                                         Max.clust=4,
                                         path.result=NULL),
                           Err_SE,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    Err_max <- paste("'Max.clust' must be an integer greater or equal to 2",
                      "and lesser than the number of genes.")
    ErrNcl_BC <- paste("The number of clusters selected for each biological",
                       "condition  must be an integer greater or equal to 2.")

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormMus2,
                                         DATAnorm=FALSE,
                                         DataNumberCluster=NULL,
                                         Method="hcpc",
                                         Max.clust=4,
                                         Membership=0.5,
                                         Min.std=0.1,
                                         Plot.Mfuzz=FALSE,
                                         path.result=NULL,
                                         Name.folder.mfuzz=NULL),
                           "Samples must belong to different times points.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormFission,
                                         DATAnorm=1.1,
                                         DataNumberCluster=NULL,
                                         Method="hcpc",
                                         Max.clust=100,
                                         Membership=0.5,
                                         Min.std=0.1,
                                         Plot.Mfuzz=FALSE,
                                         path.result=NULL,
                                         Name.folder.mfuzz=NULL),
                           "'DATAnorm' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormFission,
                                         DATAnorm=FALSE,
                                         DataNumberCluster=NULL,
                                         Method="1.1",
                                         Max.clust=4,
                                         Membership=0.5,
                                         Min.std=0.1,
                                         Plot.Mfuzz=FALSE,
                                         path.result=NULL,
                                         Name.folder.mfuzz=NULL),
                           "'Method' must be 'hcpc' or 'kmeans'.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormFission,
                                         DATAnorm=FALSE,
                                         DataNumberCluster=NULL,
                                         Method="hcpc",
                                         Max.clust=4,
                                         Membership=0.5,
                                         Min.std=0.1,
                                         Plot.Mfuzz=1.1,
                                         path.result=NULL,
                                         Name.folder.mfuzz=NULL),
                           "'Plot.Mfuzz' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormFission,
                                         DATAnorm=FALSE,
                                         DataNumberCluster=NULL,
                                         Method="hcpc",
                                         Max.clust=1.1,
                                         Membership=0.5,
                                         Min.std=0.1,
                                         Plot.Mfuzz=FALSE,
                                         path.result=NULL,
                                         Name.folder.mfuzz=NULL),
                           "'Max.clust' must be a non-negative integer.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormFission,
                                         DATAnorm=FALSE,
                                         DataNumberCluster=NULL,
                                         Method="hcpc",
                                         Max.clust=100,
                                         Membership=0.5,
                                         Min.std=0.1,
                                         Plot.Mfuzz=FALSE,
                                         path.result=1,
                                         Name.folder.mfuzz=NULL),
                           "'path.result' must be NULL or a character.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormFission,
                                         DATAnorm=FALSE,
                                         DataNumberCluster=NULL,
                                         Method="hcpc",
                                         Max.clust=100,
                                         Membership=0.5,
                                         Min.std=0.1,
                                         Plot.Mfuzz=FALSE,
                                         path.result=NULL,
                                         Name.folder.mfuzz=1),
                           "'Name.folder.mfuzz' must be NULL or a character.",
                           fixed=TRUE)

    clMfuzz <- data.frame(c("G1", "G2"), c(1, 2.3))
    testthat::expect_error(MFUZZanalysis(SEresNorm=resDATAnormSim,
                                         DATAnorm=TRUE,
                                         DataNumberCluster=clMfuzz,
                                         Method="hcpc",
                                         Max.clust=10,
                                         Membership=0.5,
                                         Min.std=0.1,
                                         Plot.Mfuzz=FALSE,
                                         path.result=NULL,
                                         Name.folder.mfuzz="Test"),
                           ErrNcl_BC,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    clFission <- data.frame(BC=c("G1"), CL=c(7))
    clMfuzz2 <- data.frame(BC=c("G1", "G2"), CL=c(3, 4))

    res1_Mfuzz <- MFUZZanalysis(SEresNorm=resDATAnormFission,
                                DATAnorm=FALSE,
                                DataNumberCluster=clFission,
                                Method="hcpc",
                                Max.clust=10,
                                Membership=0.5,
                                Min.std=0.1,
                                Plot.Mfuzz=FALSE,
                                path.result=NULL,
                                Name.folder.mfuzz=NULL)

    res2_Mfuzz <- MFUZZanalysis(SEresNorm=resDATAnormSim,
                                DATAnorm=TRUE,
                                DataNumberCluster=NULL,
                                Method="hcpc",
                                Max.clust=10,
                                Membership=0.5,
                                Min.std=0.1,
                                Plot.Mfuzz=TRUE,
                                path.result=NULL,
                                Name.folder.mfuzz="Test")

    res3_Mfuzz <- MFUZZanalysis(SEresNorm=resDATAnormSim,
                                DATAnorm=TRUE,
                                DataNumberCluster=clMfuzz2,
                                Method="hcpc",
                                Max.clust=10,
                                Membership=0.5,
                                Min.std=0.1,
                                Plot.Mfuzz=FALSE,
                                path.result=NULL,
                                Name.folder.mfuzz=NULL)

    testthat::expect_s4_class(res1_Mfuzz, "SummarizedExperiment")
    testthat::expect_s4_class(res2_Mfuzz, "SummarizedExperiment")
    testthat::expect_s4_class(res3_Mfuzz, "SummarizedExperiment")

})
