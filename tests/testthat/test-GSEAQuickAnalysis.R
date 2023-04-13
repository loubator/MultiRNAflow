test_that("Test GSEAQuickAnalysis", {
    expect_error(GSEAQuickAnalysis(Internect.Connection=TRUE,
                                   Res.DE.analysis=matrix(0,ncol=2,
                                                          nrow=3),
                                   ColumnsCriteria=1,
                                   ColumnsLog2ordered=NULL,
                                   Set.Operation="union",
                                   Organism="hsapiens",
                                   MaxNumberGO=20,
                                   Background=FALSE,
                                   Display.plots=TRUE,
                                   Save.plots=FALSE),
                 paste("Res.DE.analysis must be a list or",
                       "a 'DESeqDataSet' object"),
                 fixed=TRUE)

    MessageInternet<-paste("Once the user is sure to have",
                           "an internet connection, the user must set",
                           "'Internect.Connection=TRUE' in order to realize",
                           "the enrichment analysis", sep=" ")

    expect_message(GSEAQuickAnalysis(Internect.Connection=FALSE,
                                     Res.DE.analysis=list(matrix(0,ncol=2,
                                                                 nrow=3)),
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
