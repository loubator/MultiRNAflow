test_that("Test DATAnormalization", {
  data(RawCounts_Antoszewski2022_MOUSEsub500)
  #
  expect_error(DATAnormalization(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
                                 Column.gene=1, Individual.position=2,
                                 Group.position=NULL, Time.position=NULL,
                                 Normalization="rle", Plot.Boxplot=TRUE,
                                 Colored.By.Factors=FALSE, Color.Group=NULL,
                                 Plot.genes=FALSE, path.result=NULL,
                                 Name.folder.norm=NULL,Blind.rlog.vst=FALSE),
               "'Time.position' and 'Group.position' can not be both NULL",
               fixed=TRUE)
})
