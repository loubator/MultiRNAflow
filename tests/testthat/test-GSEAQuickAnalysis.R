test_that("Test GSEAQuickAnalysis", {
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    expect_error(GSEAQuickAnalysis(Internect.Connection=TRUE,
                                   SEresDE=matrix(0, ncol=2, nrow=3),
                                   ColumnsCriteria=1,
                                   ColumnsLog2ordered=NULL,
                                   Set.Operation="union",
                                   Organism="hsapiens",
                                   MaxNumberGO=20,
                                   Background=FALSE,
                                   Display.plots=TRUE,
                                   Save.plots=FALSE),
                 paste0("'SEresDE' mut be the results of the function ",
                        "'DEanalysisGlobal()'."),
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    MessageInternet<-paste("Once the user is sure to have",
                           "an internet connection, the user must set",
                           "'Internect.Connection=TRUE' in order to realize",
                           "the enrichment analysis", sep=" ")

    expect_message(GSEAQuickAnalysis(Internect.Connection=FALSE,
                                     SEresDE=list(matrix(0,ncol=2, nrow=3)),
                                     ColumnsCriteria=1,
                                     ColumnsLog2ordered=NULL,
                                     Set.Operation="union",
                                     Organism="hsapiens",
                                     MaxNumberGO=20,
                                     Background=FALSE,
                                     Display.plots=TRUE,
                                     Save.plots=FALSE),
                   MessageInternet,
                   fixed=TRUE)
})
