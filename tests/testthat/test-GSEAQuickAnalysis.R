testthat::test_that("Test GSEAQuickAnalysis", {
    ##-----------------------------------------------------------------------##
    data("Results_DEanalysis_sub500")
    resDEmus1 <- Results_DEanalysis_sub500$DE_Antoszewski2022_MOUSEsub500

    resDEmus1_2 <- resDEmus1
    resDEmus1_3 <- resDEmus1
    S4Vectors::metadata(resDEmus1_2)$DESeq2obj$SEidentification <- NULL
    S4Vectors::metadata(resDEmus1_3)$DESeq2obj$SEidentification <- "test"

    ## resDEleuk <- Results_DEanalysis_sub500$DE_Schleiss2021_CLLsub500
    ## resDET1wt <- Results_DEanalysis_sub500$DE_Leong2014_FISSIONsub500wt

    ##-----------------------------------------------------------------------##
    ## ##the two previous lines is equivalent to the following uncomment lines
    ## ## data importation
    ## data(RawCounts_Antoszewski2022_MOUSEsub500)
    ## ## No time points. We take only two groups for the speed of the example
    ## dataT1wt<-RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200),seq_len(7)]
    ##
    ## ## Preprocessing with Results of DEanalysisGlobal()
    ## resDATAprepSE <- DATAprepSE(RawCounts=dataT1wt,
    ##                             Column.gene=1,
    ##                             Group.position=1,
    ##                             Time.position=NULL,
    ##                             Individual.position=2)
    ## ## DE analysis
    ## resDET1wt <- DEanalysisGlobal(SEres=resDATAprepSE,
    ##                               pval.min=0.05,
    ##                               pval.vect.t=NULL,
    ##                               log.FC.min=1,
    ##                               LRT.supp.info=FALSE,
    ##                               Plot.DE.graph=FALSE,
    ##                               path.result=NULL,
    ##                               Name.folder.DE=NULL)
    ## GOmat[c(11, 15, 17, 18, 40, 46, 47, 49),], t7,t8 leuk P
    ##-----------------------------------------------------------------------##
    GOid <- c("GO:0000280", "GO:1903047", "GO:0030071", "GO:1902099",
              "GO:0000278", "GO:0010948", "GO:0001772", "KEGG:04110", "GO:---")
    GOsource <- c(rep("GO:BP", times=6), "GO:CC", "KEGG", "GO:MF")
    GOname <- c("nuclear division", "mitotic cell cycle process",
                "regulation of mitotic metaphase/anaphase transition",
                "regulation of metaphase/anaphase transition of cell cycle",
                "mitotic cell cycle",
                "negative regulation of cell cycle process",
                "immunological synapse", "Cell cycle",
                "regulation of metaphase/anaphase transition of cel")
    GOpvalue <- c(0.005206831, 0.014458134, 0.026515973, 0.031107582,
                  0.088544011, 0.144214195, 0.150970732, 0.199608750,
                  0.5)
    GOsignificant <- c(rep(TRUE, times=4), rep(FALSE, times=4), FALSE)
    GOtermsize <- c(447, 745, 90, 93, 895, 302, 46, 157, 100)
    GOquerysize <- c(rep(72, times=6), 75, 32, 50)
    GOintersection <- c(10, 12, 5, 5, 12, 7, 3, 4, 20)
    GOprecision <- c(0.13888889, 0.16666667, 0.06944444, 0.06944444,
                     0.16666667, 0.09722222, 0.04000000, 0.12500000, 0.5)
    GOrecall <- c(0.02237136, 0.01610738, 0.05555556, 0.05376344,
                  0.01340782, 0.02317881, 0.06521739, 0.02547771, 0.5)
    GOparents <- c("(GO:0048285)", "(GO:0000278)_(GO:0022402)",
                   "(GO:0007091)_(GO:1901990)_(GO:1902099)",
                   "(GO:0033045)_(GO:0044784)_(GO:1901987)", "(GO:0007049)",
                   "(GO:0010564)_(GO:0022402)_(GO:0045786)",
                   "(GO:0005886)_(GO:0110165)", "(KEGG:00000)", "GO:-1001-")
    GOgeneid <- c(paste("(AURKB)_(BUB1B)_(ESPL1)_(KIFC1)_(NCAPG)",
                        "(SPC24)_(SMC2)_(TTK)_(TPX2)_(TTN)", sep="_"),
                  paste("(ACVR1)_(AURKB)_(BUB1B)_(DTL)_(ESPL1)_(KIFC1)",
                        "(NCAPG)_(SPC24)_(SMC2)_(TTK)_(TPX2)_(TTN)", sep="_"),
                  "(AURKB)_(BUB1B)_(ESPL1)_(SPC24)_(TTK)",
                  "(AURKB)_(BUB1B)_(ESPL1)_(SPC24)_(TTK)",
                  paste("(ACVR1)_(AURKB)_(BUB1B)_(DTL)_(ESPL1)_(KIFC1)",
                        "(NCAPG)_(SPC24)_(SMC2)_(TTK)_(TPX2)_(TTN)", sep="_"),
                  "(ACVR1)_(AURKB)_(BUB1B)_(DTL)_(ESPL1)_(SPC24)_(TTK)",
                  "(CD6)_(LGALS3)_(VAV3)", "(AURKB)_(BUB1B)_(ESPL1)_(TTK)",
                  "(0001)_(0010)_(0100)_(1000)")

    GOmat_test <- data.frame(term_id=GOid, source=GOsource, term_name=GOname,
                             p_value=GOpvalue, significant=GOsignificant,
                             term_size=GOtermsize, query_size=GOquerysize,
                             intersection_size=GOintersection,
                             precision=GOprecision, recall=GOrecall,
                             GOparents=GOparents, Gene_id=GOgeneid)

    GOmat_test2 <- GOmat_test
    GOmat_test2$p_value[9] <- 1

    RESgpro_test <- RESgprofiler2(RESgost=GOmat_test, MaxNumberGO=10)
    RESgpro_test2 <- RESgprofiler2(RESgost=GOmat_test2, MaxNumberGO=10)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(GSEAQuickAnalysis(Internet.Connection=FALSE,
                                             SEresDE=matrix(0, ncol=2, nrow=3),
                                             ColumnsCriteria=1,
                                             ColumnsLog2ordered=NULL,
                                             Set.Operation="union",
                                             Organism="hsapiens",
                                             MaxNumberGO=20,
                                             Background=FALSE,
                                             Display.plots=FALSE,
                                             Save.plots=FALSE),
                           paste0("'SEresDE' mut be the results of the ",
                                  "function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(GSEAQuickAnalysis(Internet.Connection=FALSE,
                                             SEresDE=resDEmus1_2,
                                             ColumnsCriteria=1,
                                             ColumnsLog2ordered=NULL,
                                             Set.Operation="union",
                                             Organism="hsapiens",
                                             MaxNumberGO=20,
                                             Background=FALSE,
                                             Display.plots=FALSE,
                                             Save.plots=FALSE),
                           paste0("'SEresDE' mut be the results of the ",
                                  "function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    testthat::expect_error(GSEAQuickAnalysis(Internet.Connection=FALSE,
                                             SEresDE=resDEmus1_3,
                                             ColumnsCriteria=1,
                                             ColumnsLog2ordered=NULL,
                                             Set.Operation="union",
                                             Organism="hsapiens",
                                             MaxNumberGO=20,
                                             Background=FALSE,
                                             Display.plots=FALSE,
                                             Save.plots=FALSE),
                           paste0("'SEresDE' mut be the results of the ",
                                  "function 'DEanalysisGlobal()'."),
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    testthat::expect_error(GSEAQuickAnalysis(Internet.Connection=FALSE,
                                             SEresDE=resDEmus1,
                                             ColumnsCriteria=2,
                                             ColumnsLog2ordered="colselcted",
                                             Set.Operation="union",
                                             Organism="hsapiens",
                                             MaxNumberGO=20,
                                             Background=FALSE,
                                             Display.plots=TRUE,
                                             Save.plots=FALSE),
                           paste("'ColumnsLog2ordered' must be a non-negative",
                                 "integer, a vector of non-negative integers",
                                 "or 'NULL'"),
                           fixed=TRUE)

    testthat::expect_error(GSEAQuickAnalysis(Internet.Connection=FALSE,
                                             SEresDE=resDEmus1,
                                             ColumnsCriteria=2,
                                             ColumnsLog2ordered=NULL,
                                             Set.Operation="union",
                                             Organism="hsapiens",
                                             MaxNumberGO=20,
                                             Background="Background",
                                             Display.plots=TRUE,
                                             Save.plots=FALSE),
                           "'Background' must be TRUE or FALSE.",
                           fixed=TRUE)

    testthat::expect_error(GSEAQuickAnalysis(Internet.Connection="internet",
                                             SEresDE=resDEmus1,
                                             ColumnsCriteria=2,
                                             ColumnsLog2ordered=NULL,
                                             Set.Operation="union",
                                             Organism="hsapiens",
                                             MaxNumberGO=20,
                                             Background=FALSE,
                                             Display.plots=TRUE,
                                             Save.plots=FALSE),
                           "'Internet.Connection' must be TRUE or FALSE.",
                           fixed=TRUE)

    ##-----------------------------------------------------------------------##
    MessageInternet <- paste("Once the user is sure to have",
                             "an internet connection, the user must set",
                             "'Internet.Connection=TRUE' in order to realize",
                             "the enrichment analysis", sep=" ")

    testthat::expect_equal(myCollapse(c("test1", "test2")), "(test1)_(test2)")

    testthat::expect_message(GSEAQuickAnalysis(Internet.Connection=FALSE,
                                               SEresDE=resDEmus1,
                                               ColumnsCriteria=2,
                                               ColumnsLog2ordered=NULL,
                                               Set.Operation="union",
                                               Organism="hsapiens",
                                               MaxNumberGO=20,
                                               Background=FALSE,
                                               Display.plots=TRUE,
                                               Save.plots=FALSE),
                             MessageInternet,
                             fixed=TRUE)

    testthat::expect_message(GSEAQuickAnalysis(Internet.Connection=FALSE,
                                               SEresDE=resDEmus1,
                                               ColumnsCriteria=2,
                                               ColumnsLog2ordered=11,
                                               Set.Operation="union",
                                               Organism="hsapiens",
                                               MaxNumberGO=20,
                                               Background=TRUE,
                                               Display.plots=FALSE,
                                               Save.plots=TRUE),
                             MessageInternet,
                             fixed=TRUE)

    testthat::expect_s3_class(RESgpro_test$glolipop, "ggplot")
    testthat::expect_s3_class(RESgpro_test2$gManhattan, "ggplot")
})
