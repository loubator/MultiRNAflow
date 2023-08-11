test_that("Test MFUZZanalysis", {
    ##
    ##------------------------------------------------------------------------#
    ## data importation
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2<-RawCounts_Antoszewski2022_MOUSEsub500
    ## preprocessing
    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)
    resDATAnormMus2 <- resDATAprepSEmus2
    resDATAnormMus2$SEidentification <- "SEresNormalization"
    ## normalization
    # resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
    #                                      Normalization="rle",
    #                                      Plot.Boxplot=TRUE,
    #                                      Colored.By.Factors=TRUE)


    data(RawCounts_Leong2014_FISSIONsub500wt)
    datafission<-RawCounts_Leong2014_FISSIONsub500wt
    ## preprocessing
    resDATAprepSEfission <- DATAprepSE(RawCounts=datafission,
                                       Column.gene=1,
                                       Group.position=NULL,
                                       Time.position=2,
                                       Individual.position=3)
    resDATAnormFission <- resDATAprepSEfission
    resDATAnormFission$SEidentification <- "SEresNormalization"

    ##------------------------------------------------------------------------#
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")
    Err_max <- paste0("'Max.clust' must be an integer ",
                      "greater or equal to 2 and ",
                      "lesser than the number of genes.")

    expect_error(MFUZZanalysis(SEresNorm=datafission,
                               DATAnorm=TRUE,
                               DataNumberCluster=NULL,
                               Method="hcpc",
                               Max.clust=6,
                               Membership=0.5,
                               Min.std=0.1,
                               Plot.Mfuzz=TRUE,
                               path.result=NULL,
                               Name.folder.mfuzz=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(MFUZZanalysis(SEresNorm=resDATAprepSEfission,
                               DATAnorm=TRUE,
                               DataNumberCluster=NULL,
                               Method="hcpc",
                               Max.clust=6,
                               Membership=0.5,
                               Min.std=0.1,
                               Plot.Mfuzz=TRUE,
                               path.result=NULL,
                               Name.folder.mfuzz=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(MFUZZanalysis(SEresNorm=resDATAnormFission,
                               DATAnorm=FALSE,
                               DataNumberCluster=NULL,
                               Method="hcpc",
                               Max.clust=10000,
                               Membership=0.5,
                               Min.std=0.1,
                               Plot.Mfuzz=TRUE,
                               path.result=NULL,
                               Name.folder.mfuzz=NULL),
                 Err_max,
                 fixed=TRUE)

    expect_error(MFUZZanalysis(SEresNorm=resDATAnormMus2,
                               DATAnorm=FALSE,
                               DataNumberCluster=NULL,
                               Method="hcpc",
                               Max.clust=100,
                               Membership=0.5,
                               Min.std=0.1,
                               Plot.Mfuzz=TRUE,
                               path.result=NULL,
                               Name.folder.mfuzz=NULL),
                 "Samples must belong to different times points.",
                 fixed=TRUE)
})
