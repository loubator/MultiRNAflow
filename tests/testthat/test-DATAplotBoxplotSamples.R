testthat::test_that("Test DATAplotBoxplotSamples", {
    ##
    ##-----------------------------------------------------------------------##
    data("RawCounts_Antoszewski2022_MOUSEsub500")
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

    ##-----------------------------------------------------------------------##
    ## data importation
    dataSIM <- RawCountsSimulation(Nb.Group=2,
                                   Nb.Time=3,
                                   Nb.per.GT=4,
                                   Nb.Gene=100)$Sim.dat

    ## preprocessing
    SEprepSIM <- DATAprepSE(RawCounts=dataSIM,
                            Column.gene=1,
                            Group.position=1,
                            Time.position=2,
                            Individual.position=3)

    ## normalization
    SEnormSIM <- DATAnormalization(SEres=SEprepSIM,
                                   Normalization="rle",
                                   Plot.Boxplot=FALSE,
                                   Colored.By.Factors=TRUE)

    ##-----------------------------------------------------------------------##
    ## data importation
    dataSIMt <- RawCountsSimulation(Nb.Group=1,
                                    Nb.Time=3,
                                    Nb.per.GT=4,
                                    Nb.Gene=100)$Sim.dat

    ## preprocessing
    SEprepSIMt <- DATAprepSE(RawCounts=dataSIMt,
                             Column.gene=1,
                             Group.position=NULL,
                             Time.position=1,
                             Individual.position=2)

    ## normalization
    SEnormSIMt <- DATAnormalization(SEres=SEprepSIMt,
                                    Normalization="rle",
                                    Plot.Boxplot=FALSE,
                                    Colored.By.Factors=TRUE)

    ##-----------------------------------------------------------------------##
    Err_SEmus2 <- paste0("'SEres' mut be the results of either the function ",
                         "'DATAprepSE()' or 'DATAnormalization()'.")

    testthat::expect_error(DATAplotBoxplotSamples(SEres=datamus2,
                                                  Log2.transformation=TRUE,
                                                  Colored.By.Factors=TRUE,
                                                  Color.Group=NULL,
                                                  Plot.genes=FALSE,
                                                  y.label=NULL),
                           Err_SEmus2,
                           fixed=TRUE)

    testthat::expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2Id,
                                                  Log2.transformation=TRUE,
                                                  Colored.By.Factors=TRUE,
                                                  Color.Group=NULL,
                                                  Plot.genes=FALSE,
                                                  y.label=NULL),
                           Err_SEmus2,
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                                  Log2.transformation=1.1,
                                                  Colored.By.Factors=TRUE,
                                                  Color.Group=NULL,
                                                  Plot.genes=FALSE,
                                                  y.label=NULL),
                           "'Log2.transformation' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                                  Log2.transformation=FALSE,
                                                  Colored.By.Factors=1.1,
                                                  Color.Group=NULL,
                                                  Plot.genes=FALSE,
                                                  y.label=NULL),
                           "'Colored.By.Factors' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                                  Log2.transformation=FALSE,
                                                  Colored.By.Factors=FALSE,
                                                  Color.Group=NULL,
                                                  Plot.genes=1.1,
                                                  y.label=NULL),
                           "'Plot.genes' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                                  Log2.transformation=FALSE,
                                                  Colored.By.Factors=FALSE,
                                                  Color.Group=NULL,
                                                  Plot.genes=TRUE,
                                                  y.label=1.1),
                           "'y.label' must be NULL or a character.",
                           fixed=TRUE)

    testthat::expect_error(DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                                  Log2.transformation=FALSE,
                                                  Colored.By.Factors=FALSE,
                                                  Color.Group="Test",
                                                  Plot.genes=TRUE,
                                                  y.label=NULL),
                           "'Color.Group' must be NULL or a data.frame.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    ## res1_palG <- myPaletteBC(Nbc=4)
    ## testthat::expect_vector(res1_palG)
    res2_palG <- myPaletteBC(Nbc=20)
    testthat::expect_vector(res2_palG)

    ##-----------------------------------------------------------------------##
    res1_NORMboxplot <- DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                               Log2.transformation=FALSE,
                                               Colored.By.Factors=FALSE,
                                               Color.Group=NULL,
                                               Plot.genes=TRUE, y.label="test")

    res2_NORMboxplot <- DATAplotBoxplotSamples(SEres=resDATAprepSEmus2,
                                               Log2.transformation=TRUE,
                                               Colored.By.Factors=TRUE,
                                               Color.Group=dataCOL,
                                               Plot.genes=FALSE, y.label=NULL)

    res3_NORMboxplot <- DATAplotBoxplotSamples(SEres=SEnormSIM,
                                               Log2.transformation=FALSE,
                                               Colored.By.Factors=TRUE,
                                               Color.Group=NULL,
                                               Plot.genes=FALSE,
                                               y.label=NULL)

    res4_NORMboxplot <- DATAplotBoxplotSamples(SEres=SEnormSIMt,
                                               Log2.transformation=TRUE,
                                               Colored.By.Factors=TRUE,
                                               Color.Group=NULL,
                                               Plot.genes=TRUE,
                                               y.label=NULL)

    testthat::expect_s3_class(res1_NORMboxplot, "ggplot")
    testthat::expect_s3_class(res2_NORMboxplot, "ggplot")
    testthat::expect_s3_class(res3_NORMboxplot, "ggplot")
    testthat::expect_s3_class(res4_NORMboxplot, "ggplot")

})
