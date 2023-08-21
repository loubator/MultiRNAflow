test_that("Test PCAanalysis", {
    ##
    ##------------------------------------------------------------------------#
    ## data importation
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2<-RawCounts_Antoszewski2022_MOUSEsub500[seq_len(100),]

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
    ##------------------------------------------------------------------------#
    ## data importation
    dataSIM <- RawCountsSimulation(Nb.Group=2,
                                   Nb.Time=3,
                                   Nb.per.GT=5,
                                   Nb.Gene=100)$Sim.dat

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
    Err_SE <- paste0("'SEresNorm' mut be the results of the function ",
                     "'DATAnormalization().'")

    expect_error(PCAanalysis(SEresNorm=datamus2,
                             DATAnorm=TRUE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.folder.pca=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(PCAanalysis(SEresNorm=resDATAprepSEmus2,
                             DATAnorm=TRUE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.folder.pca=NULL),
                 Err_SE,
                 fixed=TRUE)

    expect_error(PCAanalysis(SEresNorm=resDATAnormMus2Id,
                             DATAnorm=TRUE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.folder.pca=NULL),
                 Err_SE,
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                             DATAnorm=1.1,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.folder.pca=NULL),
                 "'DATAnorm' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=1.1,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.folder.pca=NULL),
                 "'Supp.del.sample' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=1.1, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.folder.pca=NULL),
                 "'Plot.PCA' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=1.1,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.folder.pca=NULL),
                 "'Mean.Accross.Time' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=1.1,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.folder.pca=NULL),
                 "'D3.mouvement' must be TRUE or FALSE.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=1.1, Name.folder.pca=NULL),
                 "'path.result' must be NULL or a character.",
                 fixed=TRUE)


    expect_error(PCAanalysis(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.folder.pca=1),
                 "'Name.folder.pca' must be NULL or a character.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_s4_class(PCAanalysis(SEresNorm=resDATAnormMus2,
                                DATAnorm=FALSE,
                                gene.deletion=NULL,
                                sample.deletion=NULL, Supp.del.sample=FALSE,
                                Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                                Color.Group=NULL,
                                D3.mouvement=FALSE,
                                Phi=25, Theta=140, Cex.label=0.7,
                                Cex.point=0.7, epsilon=0.2,
                                path.result=NULL, Name.folder.pca=NULL),
                    "SummarizedExperiment")

    expect_s4_class(PCAanalysis(SEresNorm=resDATAnormSim,
                                DATAnorm=TRUE,
                                gene.deletion=NULL,
                                sample.deletion=NULL, Supp.del.sample=FALSE,
                                Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                                Color.Group=NULL,
                                D3.mouvement=FALSE,
                                Phi=25, Theta=140, Cex.label=0.7,
                                Cex.point=0.7, epsilon=0.2,
                                path.result=NULL, Name.folder.pca="test"),
                    "SummarizedExperiment")

    expect_s4_class(PCAanalysis(SEresNorm=resDATAnormSim,
                                DATAnorm=FALSE,
                                gene.deletion=NULL,
                                sample.deletion=NULL, Supp.del.sample=FALSE,
                                Plot.PCA=FALSE, Mean.Accross.Time=TRUE,
                                Color.Group=NULL,
                                D3.mouvement=FALSE,
                                Phi=25, Theta=140, Cex.label=0.7,
                                Cex.point=0.7, epsilon=0.2,
                                path.result=NULL, Name.folder.pca=NULL),
                    "SummarizedExperiment")
})
