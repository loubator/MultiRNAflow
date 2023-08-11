test_that("Test PCAgraphics", {
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
    # resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
    #                                      Normalization="rle",
    #                                      Plot.Boxplot=TRUE,
    #                                      Colored.By.Factors=TRUE)

    ##------------------------------------------------------------------------#
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    expect_error(PCAgraphics(SEresNorm=datamus2,
                             DATAnorm=TRUE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=TRUE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(PCAgraphics(SEresNorm=resDATAprepSEmus2,
                             DATAnorm=TRUE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=TRUE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 Err_SE,
                 fixed=TRUE)
})
