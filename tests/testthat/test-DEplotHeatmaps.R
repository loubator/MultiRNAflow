test_that("Test DEplotHeatmaps", {
    ##------------------------------------------------------------------------#
    expect_error(DEplotHeatmaps(SEresDE=matrix(0, ncol=2, nrow=3),
                                ColumnsCriteria=2),
                 paste0("'SEresDE' mut be the results of the function ",
                        "'DEanalysisGlobal()'."),
                 fixed=TRUE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    data(Results_DEanalysis_sub500)
    resDET1wt <- Results_DEanalysis_sub500$DE_Antoszewski2022_MOUSEsub500

    ## ##the two previous lines is equivalent to the following uncomment lines
    # ## data importation
    # data(RawCounts_Antoszewski2022_MOUSEsub500)
    # ## No time points. We take only two groups for the speed of the example
    # dataT1wt<-RawCounts_Antoszewski2022_MOUSEsub500[seq_len(200),seq_len(7)]
    #
    # ## Preprocessing with Results of DEanalysisGlobal()
    # resDATAprepSE <- DATAprepSE(RawCounts=dataT1wt,
    #                             Column.gene=1,
    #                             Group.position=1,
    #                             Time.position=NULL,
    #                             Individual.position=2)
    # ## DE analysis
    # resDET1wt <- DEanalysisGlobal(SEres=resDATAprepSE,
    #                               pval.min=0.05,
    #                               pval.vect.t=NULL,
    #                               log.FC.min=1,
    #                               LRT.supp.info=FALSE,
    #                               Plot.DE.graph=FALSE,
    #                               path.result=NULL,
    #                               Name.folder.DE=NULL)

    ##------------------------------------------------------------------------#
    resHeatmap <- DEplotHeatmaps(SEresDE=resDET1wt,
                                 ColumnsCriteria=3, ## Specific genes N1haT1ko
                                 Set.Operation="union",
                                 NbGene.analysis=20,
                                 Color.Group=NULL,
                                 SizeLabelRows=5,
                                 SizeLabelCols=5,
                                 Display.plots=TRUE,
                                 Save.plots=FALSE)

    ##------------------------------------------------------------------------#
    expect_s4_class(resHeatmap$Heatmap.Correlation,
                    "Heatmap")
})
