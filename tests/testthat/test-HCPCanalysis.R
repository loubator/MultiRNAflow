test_that("Test HCPCanalysis", {
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

    ## normalization
    resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
                                         Normalization="rle",
                                         Plot.Boxplot=TRUE,
                                         Colored.By.Factors=TRUE)

    ##------------------------------------------------------------------------#
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    expect_error(HCPCanalysis(SEresNorm=datamus2,
                              DATAnorm=TRUE,
                              gene.deletion=NULL,
                              sample.deletion=NULL, Supp.del.sample=FALSE,
                              Plot.HCPC=FALSE,
                              Color.Group=NULL,
                              D3.mouvement=FALSE,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(HCPCanalysis(SEresNorm=resDATAprepSEmus2,
                              DATAnorm=TRUE,
                              gene.deletion=NULL,
                              sample.deletion=NULL, Supp.del.sample=FALSE,
                              Plot.HCPC=FALSE,
                              Color.Group=NULL,
                              D3.mouvement=FALSE,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc=NULL),
                 Err_SE,
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                              DATAnorm=1.1,
                              gene.deletion=NULL,
                              sample.deletion=NULL, Supp.del.sample=FALSE,
                              Plot.HCPC=FALSE,
                              Color.Group=NULL,
                              D3.mouvement=FALSE,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc=NULL),
                 "'DATAnorm' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                              DATAnorm=FALSE,
                              gene.deletion=NULL,
                              sample.deletion=NULL, Supp.del.sample=1.1,
                              Plot.HCPC=FALSE,
                              Color.Group=NULL,
                              D3.mouvement=FALSE,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc=NULL),
                 "'Supp.del.sample' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                              DATAnorm=FALSE,
                              gene.deletion=NULL,
                              sample.deletion=NULL, Supp.del.sample=FALSE,
                              Plot.HCPC=1.1,
                              Color.Group=NULL,
                              D3.mouvement=FALSE,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc=NULL),
                 "'Plot.HCPC' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                              DATAnorm=FALSE,
                              gene.deletion=NULL,
                              sample.deletion=NULL, Supp.del.sample=FALSE,
                              Plot.HCPC=FALSE,
                              Color.Group=NULL,
                              D3.mouvement=1.1,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc=NULL),
                 "'D3.mouvement' must be TRUE or FALSE.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                              DATAnorm=FALSE,
                              gene.deletion=NULL,
                              sample.deletion=NULL, Supp.del.sample=FALSE,
                              Plot.HCPC=FALSE,
                              Color.Group=NULL,
                              D3.mouvement=FALSE,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=1.1, Name.folder.hcpc=NULL),
                 "'path.result' must be NULL or a character.",
                 fixed=TRUE)


    expect_error(HCPCanalysis(SEresNorm=resDATAnormMus2,
                              DATAnorm=FALSE,
                              gene.deletion=NULL,
                              sample.deletion=NULL, Supp.del.sample=FALSE,
                              Plot.HCPC=FALSE,
                              Color.Group=NULL,
                              D3.mouvement=FALSE,
                              Phi=25, Theta=140, Cex.label=0.7,
                              Cex.point=0.7, epsilon=0.2,
                              path.result=NULL, Name.folder.hcpc=1),
                 "'Name.folder.hcpc' must be NULL or a character.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_s4_class(HCPCanalysis(SEresNorm=resDATAnormMus2,
                                 DATAnorm=FALSE,
                                 gene.deletion=NULL,
                                 sample.deletion=NULL, Supp.del.sample=FALSE,
                                 Plot.HCPC=FALSE,
                                 Color.Group=NULL,
                                 D3.mouvement=FALSE,
                                 Phi=25, Theta=140, Cex.label=0.7,
                                 Cex.point=0.7, epsilon=0.2,
                                 path.result=NULL, Name.folder.hcpc=NULL),
                    "SummarizedExperiment")

    expect_s4_class(HCPCanalysis(SEresNorm=resDATAnormMus2,
                                 DATAnorm=FALSE,
                                 gene.deletion=NULL,
                                 sample.deletion=NULL, Supp.del.sample=FALSE,
                                 Plot.HCPC=TRUE,
                                 Color.Group=NULL,
                                 D3.mouvement=FALSE,
                                 Phi=25, Theta=140, Cex.label=0.7,
                                 Cex.point=0.7, epsilon=0.2,
                                 path.result=NULL, Name.folder.hcpc="test"),
                    "SummarizedExperiment")
})
