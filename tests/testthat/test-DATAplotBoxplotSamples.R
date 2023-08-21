test_that("Test DATAplotBoxplotSamples", {
    ##
    ##------------------------------------------------------------------------#
    data(RawCounts_Antoszewski2022_MOUSEsub500)
    datamus2 <- RawCounts_Antoszewski2022_MOUSEsub500[seq_len(100),]

    resDATAprepSEmus2 <- DATAprepSE(RawCounts=datamus2,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=NULL,
                                    Individual.position=2)

    resDATAprepSEmus2Id <- resDATAprepSEmus2
    S4Vectors::metadata(resDATAprepSEmus2Id)$SEidentification <- "test"

    dataCOL <- data.frame(BC=c("N1wtT1wt", "N1haT1wt", "N1haT1ko", "N1wtT1ko"),
                          COL=c("red", "blue", "green", "black"))

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
    ## data importation
    dataSIMt <- RawCountsSimulation(Nb.Group=1,
                                    Nb.Time=3,
                                    Nb.per.GT=5,
                                    Nb.Gene=100)$Sim.dat

    ## preprocessing
    resDATAprepSEsimt <- DATAprepSE(RawCounts=dataSIMt,
                                    Column.gene=1,
                                    Group.position=NULL,
                                    Time.position=1,
                                    Individual.position=2)

    ## normalization
    resDATAnormSimt <- DATAnormalization(SEres=resDATAprepSEsimt,
                                         Normalization="rle",
                                         Plot.Boxplot=FALSE,
                                         Colored.By.Factors=TRUE)

    ##------------------------------------------------------------------------#
    ## data importation
    dataSIM2 <- RawCountsSimulation(Nb.Group=1,
                                    Nb.Time=3,
                                    Nb.per.GT=5,
                                    Nb.Gene=100)$Sim.dat

    ## preprocessing
    resDATAprepSEsim2 <- DATAprepSE(RawCounts=dataSIM,
                                    Column.gene=1,
                                    Group.position=1,
                                    Time.position=2,
                                    Individual.position=3)

    ## normalization
    resDATAnormSim2 <- DATAnormalization(SEres=resDATAprepSEsim,
                                         Normalization="rle",
                                         Plot.Boxplot=FALSE,
                                         Colored.By.Factors=TRUE)

    ##------------------------------------------------------------------------#
    Err_SEmus2 <- paste0("'SEres' mut be the results of either the function ",
                         "'DATAprepSE()' or 'DATAnormalization()'.")

    expect_error(DATAplotBoxplotSamples(SEres=datamus2,
                                        Log2.transformation=TRUE,
                                        Colored.By.Factors=TRUE,
                                        Color.Group=NULL,
                                        Plot.genes=FALSE,
                                        y.label=NULL),
                 Err_SEmus2,
                 fixed=TRUE)

    expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2Id,
                                        Log2.transformation=TRUE,
                                        Colored.By.Factors=TRUE,
                                        Color.Group=NULL,
                                        Plot.genes=FALSE,
                                        y.label=NULL),
                 Err_SEmus2,
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                        Log2.transformation=1.1,
                                        Colored.By.Factors=TRUE,
                                        Color.Group=NULL,
                                        Plot.genes=FALSE,
                                        y.label=NULL),
                 "'Log2.transformation' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                        Log2.transformation=FALSE,
                                        Colored.By.Factors=1.1,
                                        Color.Group=NULL,
                                        Plot.genes=FALSE,
                                        y.label=NULL),
                 "'Colored.By.Factors' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                        Log2.transformation=FALSE,
                                        Colored.By.Factors=FALSE,
                                        Color.Group=NULL,
                                        Plot.genes=1.1,
                                        y.label=NULL),
                 "'Plot.genes' must be TRUE or FALSE.",
                 fixed=TRUE)

    expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                        Log2.transformation=FALSE,
                                        Colored.By.Factors=FALSE,
                                        Color.Group=NULL,
                                        Plot.genes=TRUE,
                                        y.label=1.1),
                 "'y.label' must be NULL or a character.",
                 fixed=TRUE)

    expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                        Log2.transformation=FALSE,
                                        Colored.By.Factors=FALSE,
                                        Color.Group="Test",
                                        Plot.genes=TRUE,
                                        y.label=NULL),
                 "'Color.Group' must be NULL or a data.frame.",
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                           Log2.transformation=FALSE,
                                           Colored.By.Factors=FALSE,
                                           Color.Group=NULL,
                                           Plot.genes=TRUE,
                                           y.label=NULL),
                    "ggplot")

    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                           Log2.transformation=FALSE,
                                           Colored.By.Factors=TRUE,
                                           Color.Group=NULL,
                                           Plot.genes=FALSE,
                                           y.label=NULL),
                    "ggplot")

    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                           Log2.transformation=FALSE,
                                           Colored.By.Factors=TRUE,
                                           Color.Group=NULL,
                                           Plot.genes=FALSE,
                                           y.label="test"),
                    "ggplot")

    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                           Log2.transformation=TRUE,
                                           Colored.By.Factors=TRUE,
                                           Color.Group=dataCOL,
                                           Plot.genes=FALSE,
                                           y.label=NULL),
                    "ggplot")

    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAnormSim,
                                           Log2.transformation=TRUE,
                                           Colored.By.Factors=TRUE,
                                           Color.Group=NULL,
                                           Plot.genes=FALSE,
                                           y.label=NULL),
                    "ggplot")

    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAnormSim,
                                           Log2.transformation=FALSE,
                                           Colored.By.Factors=TRUE,
                                           Color.Group=NULL,
                                           Plot.genes=TRUE,
                                           y.label=NULL),
                    "ggplot")

    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAnormSim2,
                                           Log2.transformation=TRUE,
                                           Colored.By.Factors=TRUE,
                                           Color.Group=NULL,
                                           Plot.genes=FALSE,
                                           y.label=NULL),
                    "ggplot")

    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAnormSim2,
                                           Log2.transformation=FALSE,
                                           Colored.By.Factors=TRUE,
                                           Color.Group=NULL,
                                           Plot.genes=TRUE,
                                           y.label=NULL),
                    "ggplot")

    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAnormSimt,
                                           Log2.transformation=FALSE,
                                           Colored.By.Factors=TRUE,
                                           Color.Group=NULL,
                                           Plot.genes=FALSE,
                                           y.label=NULL),
                    "ggplot")

    expect_s3_class(DATAplotBoxplotSamples(SEres=resDATAnormSimt,
                                           Log2.transformation=FALSE,
                                           Colored.By.Factors=TRUE,
                                           Color.Group=NULL,
                                           Plot.genes=TRUE,
                                           y.label=NULL),
                    "ggplot")

})
