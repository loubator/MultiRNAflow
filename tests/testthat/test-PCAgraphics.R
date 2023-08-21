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
    resDATAnormMus2 <- DATAnormalization(SEres=resDATAprepSEmus2,
                                         Normalization="rle",
                                         Plot.Boxplot=TRUE,
                                         Colored.By.Factors=TRUE)

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

    expect_error(PCAgraphics(SEresNorm=datamus2,
                             DATAnorm=TRUE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
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
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 Err_SE,
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=1.1,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 "'DATAnorm' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=1.1,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 "'Supp.del.sample' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=1.1, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 "'Plot.PCA' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=1.1,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 "'Mean.Accross.Time' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=1.1,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 "'D3.mouvement' must be TRUE or FALSE.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=1.1, Name.file.pca=NULL),
                 "'path.result' must be NULL or a character.",
                 fixed=TRUE)


    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=1),
                 "'Name.file.pca' must be NULL or a character.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    Err_a <- "'Phi', 'Theta' and 'epsilon' must be positive numeric values"
    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi="test", Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 Err_a,
                 fixed=TRUE)

    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point="test", epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 "'Cex.point' and 'Cex.label' must be positive numeric values",
                 fixed=TRUE)

    Err_gdel <- paste("'gene.deletion' must be either NULL",
                      "either character or non negative integers.")
    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=data.frame(c(1,1)),
                             sample.deletion=NULL, Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 Err_gdel,
                 fixed=TRUE)

    Err_sdel <- paste("'sample.deletion' must be either NULL",
                      "either character or non negative integers.")
    expect_error(PCAgraphics(SEresNorm=resDATAnormMus2,
                             DATAnorm=FALSE,
                             gene.deletion=NULL,
                             sample.deletion=data.frame(c(1,1)),
                             Supp.del.sample=FALSE,
                             Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                             Color.Group=NULL,
                             D3.mouvement=FALSE,
                             Phi=25, Theta=140, Cex.label=0.7,
                             Cex.point=0.7, epsilon=0.2,
                             path.result=NULL, Name.file.pca=NULL),
                 Err_sdel,
                 fixed=TRUE)


    ##------------------------------------------------------------------------#
    expect_s4_class(PCAgraphics(SEresNorm=resDATAnormMus2,
                                DATAnorm=FALSE,
                                gene.deletion=NULL,
                                sample.deletion=NULL, Supp.del.sample=FALSE,
                                Plot.PCA=TRUE, Mean.Accross.Time=FALSE,
                                Color.Group=NULL,
                                D3.mouvement=FALSE,
                                Phi=25, Theta=140, Cex.label=0.7,
                                Cex.point=0.7, epsilon=0.2,
                                path.result=NULL, Name.file.pca="Test"),
                    "SummarizedExperiment")

    expect_s4_class(PCAgraphics(SEresNorm=resDATAnormSim,
                                DATAnorm=FALSE,
                                gene.deletion=NULL,
                                sample.deletion=NULL, Supp.del.sample=FALSE,
                                Plot.PCA=FALSE, Mean.Accross.Time=FALSE,
                                Color.Group=NULL,
                                D3.mouvement=FALSE,
                                Phi=25, Theta=140, Cex.label=0.7,
                                Cex.point=0.7, epsilon=0.2,
                                path.result=NULL, Name.file.pca=NULL),
                    "SummarizedExperiment")

    expect_s4_class(PCAgraphics(SEresNorm=resDATAnormSim,
                                DATAnorm=FALSE,
                                gene.deletion=NULL,
                                sample.deletion=NULL, Supp.del.sample=FALSE,
                                Plot.PCA=TRUE, Mean.Accross.Time=FALSE,
                                Color.Group=NULL,
                                D3.mouvement=FALSE,
                                Phi=25, Theta=140, Cex.label=0.7,
                                Cex.point=0.7, epsilon=0.2,
                                path.result=NULL, Name.file.pca=NULL),
                    "SummarizedExperiment")

    expect_s4_class(PCAgraphics(SEresNorm=resDATAnormSim,
                                DATAnorm=FALSE,
                                gene.deletion=NULL,
                                sample.deletion=NULL, Supp.del.sample=FALSE,
                                Plot.PCA=FALSE, Mean.Accross.Time=TRUE,
                                Color.Group=NULL,
                                D3.mouvement=FALSE,
                                Phi=25, Theta=140, Cex.label=0.7,
                                Cex.point=0.7, epsilon=0.2,
                                path.result=NULL, Name.file.pca=NULL),
                    "SummarizedExperiment")

    expect_s4_class(PCAgraphics(SEresNorm=resDATAnormSim,
                                DATAnorm=FALSE,
                                gene.deletion=NULL,
                                sample.deletion=NULL, Supp.del.sample=FALSE,
                                Plot.PCA=TRUE, Mean.Accross.Time=TRUE,
                                Color.Group=NULL,
                                D3.mouvement=FALSE,
                                Phi=25, Theta=140, Cex.label=0.7,
                                Cex.point=0.7, epsilon=0.2,
                                path.result=NULL, Name.file.pca=NULL),
                    "SummarizedExperiment")

})
