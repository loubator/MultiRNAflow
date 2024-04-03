testthat::test_that("Test MFUZZclustersNumber", {
    ##
    ##-----------------------------------------------------------------------##
    ## data importation
    data("RawCounts_Antoszewski2022_MOUSEsub500")
    datamus2<-RawCounts_Antoszewski2022_MOUSEsub500[seq_len(100), seq_len(10)]

    ## preprocessing
    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)
    resDATAnormMus2 <- resDATAprepSEmus2
    ## normalization
    resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
                                         Normalization="rle",
                                         Plot.Boxplot=FALSE,
                                         Colored.By.Factors=FALSE)


    #------------------------------------------------------------------------#
    data("RawCounts_Leong2014_FISSIONsub500wt")
    datafission <- RawCounts_Leong2014_FISSIONsub500wt[seq_len(250),]

    ## preprocessing
    resDATAprepSEfission <- DATAprepSE(RawCounts=datafission,
                                       Column.gene=1,
                                       Group.position=NULL,
                                       Time.position=2,
                                       Individual.position=3)
    SEprepSEfissionIdNULL <- resDATAprepSEfission
    S4Vectors::metadata(SEprepSEfissionIdNULL)$SEidentification <- NULL

    ## normalization
    resDATAnormFission <- DATAnormalization(SEres=resDATAprepSEfission,
                                            Normalization="rle",
                                            Plot.Boxplot=FALSE,
                                            Colored.By.Factors=FALSE)

    #------------------------------------------------------------------------#
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

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=datafission,
                                               DATAnorm=TRUE,
                                               Method="hcpc",
                                               Max.clust=4,
                                               Plot.Cluster=TRUE,
                                               path.result=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAprepSEfission,
                                               DATAnorm=TRUE,
                                               Method="hcpc",
                                               Max.clust=4,
                                               Plot.Cluster=TRUE,
                                               path.result=NULL),
                           Err_SE,
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=SEprepSEfissionIdNULL,
                                               DATAnorm=TRUE,
                                               Method="hcpc",
                                               Max.clust=4,
                                               Plot.Cluster=TRUE,
                                               path.result=NULL),
                           Err_SE,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    Err_max1 <- paste0("'Max.clust' must be an integer ",
                       "greater or equal to 2 and ",
                       "lesser than the number of genes.")
    Err_max2 <- "'Max.clust' must be a non-negative integer."

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAnormMus2,
                                               DATAnorm=FALSE,
                                               Method="hcpc",
                                               Max.clust=4,
                                               Plot.Cluster=FALSE,
                                               path.result=NULL),
                           "Samples must belong to different times points.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                               DATAnorm=1.1,
                                               Method="hcpc",
                                               Max.clust=4,
                                               Plot.Cluster=FALSE,
                                               path.result=NULL),
                           "'DATAnorm' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                               DATAnorm=TRUE,
                                               Method="1.1",
                                               Max.clust=4,
                                               Plot.Cluster=FALSE,
                                               path.result=NULL),
                           "'Method' must be 'hcpc' or 'kmeans'.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                               DATAnorm=TRUE,
                                               Method="hcpc",
                                               Max.clust=4,
                                               Plot.Cluster=1.2,
                                               path.result=NULL),
                           "'Plot.Cluster' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                               DATAnorm=TRUE,
                                               Method="hcpc",
                                               Max.clust=1,
                                               Plot.Cluster=FALSE,
                                               path.result=NULL),
                           Err_max1,
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                               DATAnorm=TRUE,
                                               Method="hcpc",
                                               Max.clust=-1,
                                               Plot.Cluster=FALSE,
                                               path.result=NULL),
                           Err_max2,
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                               DATAnorm=TRUE,
                                               Method="hcpc",
                                               Max.clust=4,
                                               Min.std="coucou",
                                               Plot.Cluster=FALSE,
                                               path.result=NULL),
                           "'Min.std' must be positive numeric values",
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                               DATAnorm=TRUE,
                                               Method="hcpc",
                                               Max.clust=4,
                                               Min.std=-2,
                                               Plot.Cluster=FALSE,
                                               path.result=NULL),
                           "'Min.std' must be positive numeric values",
                           fixed=TRUE)

    testthat::expect_error(MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                               DATAnorm=TRUE,
                                               Method="hcpc",
                                               Max.clust=4,
                                               Plot.Cluster=FALSE,
                                               path.result=1),
                           "'path.result' must be NULL or a character.",
                           fixed=TRUE)

    testthat::expect_error(NbClustKmeansHCPC(10, x1=9, y1=8, x2=3, y2=4),
                           paste0("'x2' must be strictly greater than 'x1'",
                                  " and ",
                                  "'y2' must be strictly greater than 'y1'"),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    res1_MfuzzNcl <- MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                         DATAnorm=TRUE, Method="kmeans",
                                         Max.clust=4,
                                         Plot.Cluster=FALSE, path.result=NULL)

    res2_MfuzzNcl <- MFUZZclustersNumber(SEresNorm=resDATAnormSim,
                                         DATAnorm=TRUE, Method="hcpc",
                                         Max.clust=4,
                                         Plot.Cluster=TRUE, path.result=NULL)

    res3_MfuzzNcl <- MFUZZclustersNumber(SEresNorm=resDATAnormFission,
                                         DATAnorm=FALSE, Method="hcpc",
                                         Max.clust=4,
                                         Plot.Cluster=TRUE, path.result=NULL)

    testthat::expect_s4_class(res1_MfuzzNcl, "SummarizedExperiment")
    testthat::expect_s4_class(res2_MfuzzNcl, "SummarizedExperiment")
    testthat::expect_s4_class(res3_MfuzzNcl, "SummarizedExperiment")

})
