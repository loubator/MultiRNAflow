## code to prepare `DATASET` dataset goes here

# usethis::use_data_raw()

# usethis::use_data(DATASET, overwrite = TRUE)

##===========================================================================##
# Useful function
##===========================================================================##
RdSel<-function(NumericVector, Nselected="max", SortedCriteriaVector=NULL){
  if(Nselected=="max"){
    SubVector<-NumericVector
  }else{
    if(Nselected>0){
      NselectedF<-min(Nselected,length(NumericVector))
      if(is.null(SortedCriteriaVector)==TRUE){
        RdIdSel<-sort(sample(x=1:length(NumericVector),size=NselectedF))
        SubVector<-NumericVector[as.numeric(RdIdSel)]
      }else{
        SortIdSel<-sort(order(SortedCriteriaVector,
                              decreasing=TRUE)[1:NselectedF])
        SubVector<-NumericVector[as.numeric(SortIdSel)]
      }# if(is.null(SortedCriteriaVector)==TRUE)
    }else{
      SubVector<-NULL
    }# if(Nselected>0)
  }# if(Nselected=="max")
  return(SubVector)
}# (RdSel)

##===========================================================================##
## GEO DataSet accession GSE169116 (4 Biocond (NOTCH1,TCF), 1 time, 3 replicats)
##===========================================================================##
#(20220121)_(Gene expression profiling of LSK cells with hyperactivated Notch1
# and/or Tcf1 knock-out)_(Antoszewski.M et al)_(Mus musculus)

DirDatRpackage<-"/Users/rodolphe/BUL/UL_THESE/CreationPackage/Data_Rpackage"
name.GSE169116<-"GSE169116_annotated_raw_counts"
DataGenes.3<-read.table(file=paste(DirDatRpackage,"/", name.GSE169116, ".txt",
                                   sep=""),
                        header=TRUE, sep="\t")
# head(DataGenes.3)
CountDat_Mus_Antoszewski2022<-DataGenes.3
colnames(CountDat_Mus_Antoszewski2022)<-c("Gene",
                                          paste(rep(c("N1wtT1wt","N1haT1wt",
                                                      "N1haT1ko","N1wtT1ko"),
                                                    each=3),
                                                paste("r",1:12,sep=""),
                                                sep="_"))
# usethis::use_data(CountDat_Mus_Antoszewski2022, overwrite = TRUE)
# c("NOTCHWTandTCFWT","NICDandTCFWT","NICDandTCFKO","NOTCHWTandTCFKO")

##---------------------------------------------------------------------------##
## Preprocessing step
resDATAprepSEmus1 <- DATAprepSE(RawCounts=CountDat_Mus_Antoszewski2022,
                                Column.gene=1,
                                Group.position=1,
                                Time.position=NULL,
                                Individual.position=2)

DEresMus<-DEanalysisGlobal(SEres=resDATAprepSEmus1,
                           pval.min=0.05,
                           pval.vect.t=NULL,
                           log.FC.min=1,
                           LRT.supp.info=FALSE,
                           path.result=NULL,
                           Name.folder.DE=NULL)

mus1DEres <- SummarizedExperiment::rowData(DEresMus)
##---------------------------------------------------------------------------##
Id.Spe.Mus<-sort(unique(c(RdSel(which(mus1DEres$Specific.genes_N1haT1ko>0),111),
                          RdSel(which(mus1DEres$Specific.genes_N1haT1wt>0),20),
                          RdSel(which(mus1DEres$Specific.genes_N1wtT1ko>0),80),
                          RdSel(which(mus1DEres$Specific.genes_N1wtT1wt>0),90),
                          RdSel(which(mus1DEres$DE.1pair.of.Group.minimum>0),
                                60))))# 357
#
IdNoDE.Mus<-RdSel(setdiff(c(1:39017),
                          which(mus1DEres$DE.1pair.of.Group.minimum>0)),
                  500-length(Id.Spe.Mus))
IdMus500<-sort(c(Id.Spe.Mus, IdNoDE.Mus))
RawCounts_Antoszewski2022_MOUSEsub500<-CountDat_Mus_Antoszewski2022[IdMus500,]
#
# usethis::use_data(RawCounts_Antoszewski2022_MOUSEsub500, overwrite=TRUE)
# usethis::use_r("RawCounts_Antoszewski2022_MOUSEsub500")
#

# Id.T1.wtVSko<-sort(c(which(mus1DEres$DE.per.pair.group..N1wtT1wt..N1wtT1ko.>0)[1:30],
#                      IdNoDE.Mus[1:20]))
# colname.T1.wtVSko<-c(1:4, 11:13)
# CountDat_Mus_Antoszewski2022_sub50Tcf1<-CountDat_Mus_Antoszewski2022[Id.T1.wtVSko,
#                                                                      colname.T1.wtVSko]
# colnames(CountDat_Mus_Antoszewski2022_sub50Tcf1)<-gsub("N1wt", "",
#                                                        colnames(CountDat_Mus_Antoszewski2022_sub50Tcf1),
#                                                        fixed=TRUE)
#
# usethis::use_data(CountDat_Mus_Antoszewski2022_sub50Tcf1, overwrite = TRUE)
# usethis::use_r("CountDat_Mus_Antoszewski2022_sub50Tcf1")
##---------------------------------------------------------------------------##
#
##===========================================================================##
## GEO DataSet accession GSE56761 (2 Biocond (wt,mut), 9 times, 3 replicats)
##===========================================================================##
# (2014)_(A global non coding RNA system modulates fission yeast protein levels
# in response to stress)_(Leong et al)_(Schizosaccharomyces pombe)
#
# name.GSE56761="GSE56761_RPKM_expressions"
# DataGenes.2<-read.table(file=paste(DirDatRpackage,"/", name.GSE56761,".",
#                                    "txt",sep=""),
#                         header = TRUE, sep = "\t")
# head(DataGenes.2)
# library(fission)
# data("fission")
load(file.path(DirDatRpackage, "fission.RData"))
CountFission<-fission@assays$data$counts
# head(CountFission) # fission@colData@listData$id
RawCounts_Leong2014_FISSION<-data.frame(Gene=row.names(CountFission),
                                        CountFission)

colnames(RawCounts_Leong2014_FISSION)[-1]<-paste(rep(c("wt","mut"),each=18),
                                                 paste("t",c(rep(0:5,each=3),
                                                             rep(0:5,each=3)),
                                                       sep=""),
                                                 paste("r",c(rep(1:3,times=6),
                                                             rep(4:6,times=6)),
                                                       sep=""),
                                                 sep="_")
# usethis::use_data(RawCounts_Leong2014_FISSION, overwrite = TRUE)
RawCounts_Leong2014_FISSIONwt<-RawCounts_Leong2014_FISSION[,1:19]
##---------------------------------------------------------------------------##
resDATAprepSEfission <- DATAprepSE(RawCounts=RawCounts_Leong2014_FISSIONwt,
                                   Column.gene=1,
                                   Group.position=NULL,
                                   Time.position=2,
                                   Individual.position=3)

DEresYeastWt<-DEanalysisGlobal(SEres=resDATAprepSEfission,
                               pval.min=0.05,
                               pval.vect.t=NULL,
                               log.FC.min=1,
                               LRT.supp.info=FALSE,
                               path.result=NULL,
                               Name.folder.DE=NULL)

fissionDEres <- SummarizedExperiment::rowData(DEresYeastWt)

##---------------------------------------------------------------------------##
DEtimePattern<-fissionDEres$DE.Temporal.Pattern
WTfissionPeak<-sort(c(RdSel(grep(".00001",DEtimePattern,fixed=TRUE), "max"),
                      RdSel(grep(".0001",DEtimePattern,fixed=TRUE), "max"),
                      RdSel(grep(".001",DEtimePattern,fixed=TRUE),90),
                      RdSel(grep(".01",DEtimePattern,fixed=TRUE),100),
                      RdSel(grep(".01",DEtimePattern,fixed=TRUE),95)))# 313
#
WTfission1tmin<-which(fissionDEres$DE.1time.minimum==1)
# WTfission1tmin.NoPeak<-unique(WTfission1tmin[!WTfission1tmin%in%WTfissionPeak])# union(setdiff(a,b),setdiff(a,b))
WTfission1tmin.NoPeak.f<-RdSel(WTfission1tmin[-intersect(WTfission1tmin,WTfissionPeak)],
                               400-length(WTfissionPeak))
#
DEfissionWt<-sort(unique(c(WTfissionPeak, WTfission1tmin.NoPeak.f)))
noDEfissionWt<-RdSel(which(fissionDEres$DE.1time.minimum==0),
                     500-length(DEfissionWt))
##---------------------------------------------------------------------------##
fission500genes<-sort(unique(c(DEfissionWt, noDEfissionWt)))
# c(1:36,38:109), c(37,110,273,327,345,355,404,432,545)
RawCounts_Leong2014_FISSIONsub500wt<-RawCounts_Leong2014_FISSIONwt[fission500genes,]
#
# usethis::use_data(RawCounts_Leong2014_FISSIONsub500wt, overwrite=TRUE)
# usethis::use_r("RawCounts_Leong2014_FISSIONsub500wt") # creer une fonction
##---------------------------------------------------------------------------##
# Id.Yeast.mut3t<-sort(c(grep(".01",DEresYeast$DE.results$Pattern.DE_mut,fixed=TRUE)[1:11],
#                        grep(".1",DEresYeast$DE.results$Pattern.DE_mut,fixed=TRUE)[1:12]))
# Id.Yeast.wt3t<-sort(c(grep(".01",DEresYeast$DE.results$Pattern.DE_wt,fixed=TRUE)[1:14],
#                       grep(".1",DEresYeast$DE.results$Pattern.DE_wt,fixed=TRUE)[1:13]))
# #
# I.sel.DEyeast10<-unique(which(DEresYeast$DE.results$Intersect_mut_T>0),
#                         which(DEresYeast$DE.results$Intersect_wt_T>0))
# #
# IdYeast75.3t<-sort(unique(c(Id.Yeast.mut3t,Id.Yeast.wt3t,I.sel.DEyeast10,
#                             which(DEresYeast$DE.results$DE.1tmin_mut>0)[1:16],
#                             which(DEresYeast$DE.results$DE.1tmin_wt>0)[1:20],
#                             noDEfission[1:10])))
# #
# Id.Mus.3T<-c(1:10,20:28)
# #
# CountDat_FissionPombe_Leong2014_3Tsub75<-RawCounts_Leong2014_FISSION[IdYeast75.3t,Id.Mus.3T]
#
# usethis::use_data(CountDat_FissionPombe_Leong2014_3Tsub75, overwrite = TRUE)
# usethis::use_r("CountDat_FissionPombe_Leong2014_3Tsub75") # pour creer une fonction
##---------------------------------------------------------------------------##
#
##===========================================================================##
## GEO DataSet accession GSE130385 (2 Biocond (P,NP), 9 times, 3 replicats)
##===========================================================================##
# (2021)_(Temporal transcriptional response from primary human chronic lymphocytic
# leukemia (CLL)-cells after B-cell receptor stimulation.)_(Schleiss.C, Vallat.L
# et al)_(Homo sapiens)

DirDatRpackage<-"/Users/rodolphe/BUL/UL_THESE/CreationPackage/Data_Rpackage"
name.GSE130385<-"GSE130385_Raw_counts"
DataGenes<-read.table(file=paste(DirDatRpackage,"/",
                                 name.GSE130385, ".", "csv", sep=""),
                      header=TRUE, sep=";")
# head(DataGenes)
RawCounts_Schleiss2021_CLL<-DataGenes
colnames(RawCounts_Schleiss2021_CLL)[2:55]<-paste("CLL",
                                                  rep(c("P","NP"),each=27),
                                                  paste("r",rep(1:6,each=9),
                                                        sep=""),
                                                  paste("t",rep(0:8,times=6),
                                                        sep=""),
                                                  sep="_")
# usethis::use_data(RawCounts_Schleiss2021_CLL, overwrite = TRUE)
# usethis::use_r("RawCounts_Schleiss2021_CLL") # pour creer une fonction
##---------------------------------------------------------------------------##
resDATAprepSEleuk <- DATAprepSE(RawCounts=RawCounts_Schleiss2021_CLL,
                                Column.gene=1,
                                Group.position=2,
                                Time.position=4,
                                Individual.position=3)

DEresCLL<-DEanalysisGlobal(SEres=resDATAprepSEleuk,
                           pval.min=0.05,
                           pval.vect.t=NULL,
                           log.FC.min=1,
                           LRT.supp.info=FALSE,
                           path.result=NULL,
                           Name.folder.DE=NULL)

leukDEres <- SummarizedExperiment::rowData(DEresCLL)

##---------------------------------------------------------------------------##
PatDE_NP<-leukDEres$Pattern.DE_NP
PatDE_P<-leukDEres$Pattern.DE_P

set.seed(1994)

NPcll<-sort(c(RdSel(grep(".00000001",PatDE_NP,fixed=TRUE),30),
              RdSel(grep(".0000001",PatDE_NP,fixed=TRUE),25),
              RdSel(grep(".000001",PatDE_NP,fixed=TRUE),30),
              RdSel(grep(".00001",PatDE_NP,fixed=TRUE),20),
              RdSel(grep(".0001",PatDE_NP,fixed=TRUE),20),
              RdSel(grep(".001",PatDE_NP,fixed=TRUE),70),
              RdSel(grep(".01",PatDE_NP,fixed=TRUE),60),
              RdSel(grep(".1",PatDE_NP,fixed=TRUE),20))) # 275
Pcll<-sort(c(RdSel(grep(".00000001",PatDE_P,fixed=TRUE),30),
             RdSel(grep(".0000001",PatDE_P,fixed=TRUE),40),
             RdSel(grep(".000001",PatDE_P,fixed=TRUE),45),
             RdSel(grep(".00001",PatDE_P,fixed=TRUE),20),
             RdSel(grep(".0001",PatDE_P,fixed=TRUE),25),
             RdSel(grep(".001",PatDE_P,fixed=TRUE),90),
             RdSel(grep(".01",PatDE_P,fixed=TRUE),55),
             RdSel(grep(".1",PatDE_P,fixed=TRUE),20)))# 325
#
SignatureNP<-which(apply(leukDEres[,2:9],1,sum)>0)# 609
SignatureP<-which(apply(leukDEres[,10:17],1,sum)>0)# 568
#
Spe.P.NP.1tmin<-which(leukDEres$Specific.genes_NP_1t.min>0)# 2121
##setdiff(NPcll,Pcll),
DEgeneCLL<-sort(unique(c(NPcll[-which(NPcll%in%intersect(NPcll,Pcll))],
                         RdSel(SignatureNP, 50), RdSel(SignatureP, 70),
                         RdSel(Spe.P.NP.1tmin, 80))))#452
#
DE1tmin_NP<-which(leukDEres$DE.1time.minimum_NP>0)
DE1tmin_P<-which(leukDEres$DE.1time.minimum_P>0)
noDEgeneCLL<-setdiff(c(1:25369),
                     sort(unique(c(DE1tmin_NP, DE1tmin_P, Spe.P.NP.1tmin))))
noDEgeneCLL<-RdSel(noDEgeneCLL,500-length(DEgeneCLL))

##---------------------------------------------------------------------------##
cll500genes<-sort(c(DEgeneCLL, noDEgeneCLL))
RawCounts_Schleiss2021_CLLsub500<-RawCounts_Schleiss2021_CLL[cll500genes,]
#
# usethis::use_data(RawCounts_Schleiss2021_CLLsub500, overwrite=TRUE)
# usethis::use_r("RawCounts_Schleiss2021_CLLsub500") # pour creer une fonction
##---------------------------------------------------------------------------##
#
##===========================================================================##
## GEO DataSet accession GSE135898
## (4 Biocond (BmKo,BmWt,CrKo,CrWt), 6 times, 4 replicats)
##===========================================================================##
# (2021)_(Temporal profiles of gene expression in Cry1/2 KO, Bmal1 KO under
# night restricted feeding and ad libitum feeding regimen.)_(Weger BD, Gobet C,
# Martin E, Gachon F, Naef F)_(Mus musculus)
#
dat.geo<-"/Users/rodolphe/BUL/UL_THESE/CreationPackage/GSE135898_Data/GSE135898_counts_v2.txt"
Datageo<-read.table(file=dat.geo, header=TRUE, sep="\t", dec=",", row.names=1)
# sort(colnames(Datageo))
# Beware, the first column corresponds to "Bmal1_KO1_AL_RNA_zt04"

Datageo2 <- Datageo[-which(is.na(Datageo), arr.ind = TRUE)[, 1],]

cl0<-c("Bmal1_KO1_AL_RNA_zt04", colnames(Datageo)[-96])
cl1<-gsub("CRY_KO_1", "CRY_KO1", x=cl0, fixed = TRUE)
cl2<-gsub("_KO", "_KO_", x=cl1, fixed = TRUE)
cl3<-gsub("_WT", "_WT_", x=cl2, fixed = TRUE)
cl4<-gsub("_RNA_", "_", x=cl3, fixed = TRUE)
cl5<-gsub("Bmal1", "Bm", x=cl4, fixed = TRUE)
cl6<-gsub("CRY", "Cr", x=cl5, fixed = TRUE)
cl7<-gsub("_AL_", "_Al_", x=cl6, fixed = TRUE)
cl8<-gsub("_RF_", "_Rf_", x=cl7, fixed = TRUE)
cl9<-gsub("_KO_", "Ko_", x=cl8, fixed = TRUE)
cl10<-gsub("_WT_", "Wt_", x=cl9, fixed = TRUE)
cl11<-gsub("_zt20", "_t5", x=cl10, fixed = TRUE)
cl12<-gsub("_zt16", "_t4", x=cl11, fixed = TRUE)
cl13<-gsub("_zt12", "_t3", x=cl12, fixed = TRUE)
cl14<-gsub("_zt08", "_t2", x=cl13, fixed = TRUE)
cl15<-gsub("_zt04", "_t1", x=cl14, fixed = TRUE)
cl16<-gsub("_zt0", "_t0", x=cl15, fixed = TRUE)
#
Colnames.matMus<-matrix(unlist(strsplit(cl16, split="_", fixed=TRUE)),
                        ncol=length(cl16))
#
Alim.delete<-TRUE
if(Alim.delete==TRUE){
  colnames(Datageo2)<-paste(Colnames.matMus[1,], Colnames.matMus[4,],
                            paste("r",rep(1:16,each=6),sep=""), sep="_")
}else{
  colnames(Datageo2)<-paste(paste(Colnames.matMus[1,], Colnames.matMus[3,],
                                  sep=""),
                            Colnames.matMus[4,],
                            paste("r",rep(1:16,each=6),sep=""),sep="_")
}# if(Alim.delete==TRUE)
# Colnames.matMus[2,]

SumDatageo2<-apply(Datageo2[,-1], MARGIN=1, sum)
varDatageo2<-apply(Datageo2[,-1], MARGIN=1, var)
# which(is.na(Datageo2),arr.ind = TRUE)
# Union.inter=which(SumDatageo2>100)#14986
# Union.inter=which(SumDatageo2>0)#26452
Union.inter<-intersect(which(SumDatageo2>20), which(varDatageo2>10))#12495
Datageo3<-data.frame(Gene=row.names(Datageo2)[Union.inter],
                     Datageo2[Union.inter,])
#
Path.save.GSE135898<-"/Users/rodolphe/BUL/UL_THESE/CreationPackage/GSE135898_Data"
#

##---------------------------------------------------------------------------##
resDATAprepSEmus2 <- DATAprepSE(RawCounts=Datageo3,#Datageo.f,# Datageo.Ncoln.3
                                Column.gene=1,
                                Group.position=1,
                                Time.position=2,
                                Individual.position=3)

ResGSE135898 <- DEanalysisGlobal(SEres=resDATAprepSEmus2,
                                 log.FC.min=1,
                                 pval.min=0.05,
                                 pval.vect.t=NULL,
                                 LRT.supp.info=FALSE,
                                 path.result=NULL,#Path.save.GSE135898,
                                 Name.folder.DE=NULL)

mus2DEres <- SummarizedExperiment::rowData(ResGSE135898)

##---------------------------------------------------------------------------##
# The original dataset has 25369 genes but we kept only 500 genes
# in order to increase the speed of each function in our algorithm.
PatDE_BmKo <- mus2DEres$Pattern.DE_BmKo
PatDE_BmWt <- mus2DEres$Pattern.DE_BmWt
PatDE_CrKo <- mus2DEres$Pattern.DE_CrKo
PatDE_CrWt <- mus2DEres$Pattern.DE_CrWt

set.seed(1994)

MusBmKo <- sort(c(grep(".00001",PatDE_BmKo,fixed=TRUE),
                  grep(".0001",PatDE_BmKo,fixed=TRUE),
                  RdSel(grep(".001",PatDE_BmKo,fixed=TRUE),30),
                  RdSel(grep(".01",PatDE_BmKo,fixed=TRUE),15),
                  RdSel(grep(".1",PatDE_BmKo,fixed=TRUE),8))) # 63

MusBmWt <- sort(c(grep(".00001",PatDE_BmWt,fixed=TRUE),
                  RdSel(grep(".0001",PatDE_BmWt,fixed=TRUE),60),
                  RdSel(grep(".001",PatDE_BmWt,fixed=TRUE),130),
                  RdSel(grep(".01",PatDE_BmWt,fixed=TRUE),155),
                  grep(".1",PatDE_BmWt,fixed=TRUE))) # 355

MusCrKo <- sort(c(grep(".00001",PatDE_CrKo,fixed=TRUE),
                  grep(".0001",PatDE_CrKo,fixed=TRUE),
                  grep(".001",PatDE_CrKo,fixed=TRUE),
                  grep(".01",PatDE_CrKo,fixed=TRUE),
                  grep(".1",PatDE_CrKo,fixed=TRUE))) # 48

MusCrWt <- sort(c(grep(".00001",PatDE_CrWt,fixed=TRUE),
                  RdSel(grep(".0001",PatDE_CrWt,fixed=TRUE),35),
                  RdSel(grep(".001",PatDE_CrWt,fixed=TRUE),145),
                  RdSel(grep(".01",PatDE_CrWt,fixed=TRUE),150),
                  grep(".1",PatDE_CrWt,fixed=TRUE))) # 367

SignatureBmKo <- which(apply(mus2DEres[,2:6],1,sum)>0)# 39
SignatureBmWt <- which(apply(mus2DEres[,7:11],1,sum)>0)# 4
SignatureCrKo <- which(apply(mus2DEres[,12:16],1,sum)>0)# 11
SignatureCrWt <- which(apply(mus2DEres[,17:21],1,sum)>0)# 20

Spe.BmKo.1tmin <- which(mus2DEres$Specific.genes_BmKo_1t.min>0)# 611
Spe.BmWt.1tmin <- which(mus2DEres$Specific.genes_BmWt_1t.min>0)# 7
Spe.CrKo.1tmin <- which(mus2DEres$Specific.genes_CrKo_1t.min>0)# 203
Spe.CrWt.1tmin <- which(mus2DEres$Specific.genes_CrWt_1t.min>0)# 33

DEMus.Sel <- sort(unique(c(setdiff(setdiff(setdiff(MusBmKo,MusBmWt), MusCrKo),
                                   MusCrWt),
                           RdSel(SignatureBmKo, 20),
                           SignatureBmWt, SignatureCrKo, SignatureCrWt,
                           RdSel(Spe.BmKo.1tmin, 230),
                           Spe.BmWt.1tmin,
                           RdSel(Spe.CrKo.1tmin, 100),
                           RdSel(Spe.CrWt.1tmin, 15))))# 413

DE1tmin_BmKo <- which(mus2DEres$DE.1time.minimum_BmKo>0)
DE1tmin_BmWt <- which(mus2DEres$DE.1time.minimum_BmWt>0)
DE1tmin_CrKo <- which(mus2DEres$DE.1time.minimum_CrKo>0)
DE1tmin_CrWt <- which(mus2DEres$DE.1time.minimum_CrWt>0)

NoDEMus <- setdiff(seq_len(nrow(mus2DEres)),
                   sort(unique(c(DE1tmin_BmKo, DE1tmin_BmWt,
                                 DE1tmin_CrKo, DE1tmin_CrWt,
                                 Spe.BmKo.1tmin, Spe.BmWt.1tmin,
                                 Spe.CrKo.1tmin, Spe.CrWt.1tmin))))
DEMus.NoSel <- NoDEMus[sort(sample(x=seq_len(length(NoDEMus)),
                                   size=500-length(DEMus.Sel)))]
##---------------------------------------------------------------------------##
MusBmCrKoWt <- sort(c(DEMus.Sel, DEMus.NoSel))
# id.Colnames.Mus<-c(1, 1+order(rep(seq(0,90,6),each=6)
# +rep(c(2,3,1,4,5,6),times=16)))
id.Colnames.Mus <- seq_len(ncol(Datageo3))
RawCounts_Weger2021_MOUSEsub500 <- Datageo3[MusBmCrKoWt, id.Colnames.Mus]
#
# usethis::use_data(RawCounts_Weger2021_MOUSEsub500, overwrite = TRUE)
# usethis::use_r("RawCounts_Weger2021_MOUSEsub500")
#
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
## DE results

resSEmus1500 <- DATAprepSE(RawCounts=RawCounts_Antoszewski2022_MOUSEsub500,
                           Column.gene=1,
                           Group.position=1,
                           Time.position=NULL,
                           Individual.position=2)

resDEMus1<-DEanalysisGlobal(SEres=resSEmus1500,
                            pval.min=0.05,
                            pval.vect.t=NULL,
                            log.FC.min=1,
                            LRT.supp.info=FALSE,
                            Plot.DE.graph=FALSE,
                            path.result=NULL,
                            Name.folder.DE=NULL)

resSEfission500 <- DATAprepSE(RawCounts=RawCounts_Leong2014_FISSIONsub500wt,
                              Column.gene=1,
                              Group.position=NULL,
                              Time.position=2,
                              Individual.position=3)

resDEFission<-DEanalysisGlobal(SEres=resSEfission500,
                               pval.min=0.05,
                               pval.vect.t=NULL,
                               log.FC.min=1,
                               LRT.supp.info=FALSE,
                               Plot.DE.graph =FALSE,
                               path.result=NULL,#FolderResultsFISSION
                               Name.folder.DE=NULL)

resSEleuk500 <- DATAprepSE(RawCounts=RawCounts_Schleiss2021_CLLsub500,
                           Column.gene=1,
                           Group.position=2,
                           Time.position=4,
                           Individual.position=3)

resDELeuk<-DEanalysisGlobal(SEres=resSEleuk500,
                            pval.min=0.05,
                            pval.vect.t=NULL,
                            log.FC.min=1,
                            LRT.supp.info=FALSE,
                            Plot.DE.graph =FALSE,
                            path.result=NULL,
                            Name.folder.DE=NULL)

resSEmus2500 <- DATAprepSE(RawCounts=RawCounts_Weger2021_MOUSEsub500,
                           Column.gene=1,
                           Group.position=1,
                           Time.position=2,
                           Individual.position=3)

resDEMUS2<-DEanalysisGlobal(SEres=resSEmus2500,
                            pval.min=0.05,
                            pval.vect.t=NULL,
                            log.FC.min=1,
                            LRT.supp.info=FALSE,
                            Plot.DE.graph =FALSE,
                            path.result=NULL,
                            Name.folder.DE=NULL)

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
## names(metadata(resDEMus1)$DESeq2obj)
resDEmus1SUB <- resDEMus1
colDatamus1 <- SummarizedExperiment::colData(resDEMus1)[, c(1, 2)]
RLEmus1 <- SummarizedExperiment::assays(resDEMus1)##$rle
DESeq2mus1 <- S4Vectors::metadata(resDEMus1)$DESeq2obj[c(4, 5, 6, 7)] ## 2
metaDEmus1 <- S4Vectors::metadata(resDEMus1)[-c(1, 2)]

SummarizedExperiment::colData(resDEmus1SUB) <- colDatamus1
SummarizedExperiment::assays(resDEmus1SUB) <- RLEmus1 ##list(rle=RLEmus1)
S4Vectors::metadata(resDEmus1SUB) <- metaDEmus1
S4Vectors::metadata(resDEmus1SUB)$DESeq2obj <- DESeq2mus1

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
## names(metadata(resDEFission)$DESeq2obj)
resDEfissionSUB <- resDEFission
colDatafission <- SummarizedExperiment::colData(resDEFission)[, c(1, 2)]
RLEfission <- SummarizedExperiment::assays(resDEFission)##$rle
DESeq2fission <- S4Vectors::metadata(resDEFission)$DESeq2obj[c(4, 5, 6, 7)]
metaDEfission <- S4Vectors::metadata(resDEFission)[-c(1, 2)]

SummarizedExperiment::colData(resDEfissionSUB) <- colDatafission
SummarizedExperiment::assays(resDEfissionSUB) <- RLEfission##list(rle=RLEfission)
S4Vectors::metadata(resDEfissionSUB) <- metaDEfission
S4Vectors::metadata(resDEfissionSUB)$DESeq2obj <- DESeq2fission

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
## names(metadata(resDELeuk)$DESeq2obj)
resDELeukSUB <- resDELeuk
colDataLeuk <- SummarizedExperiment::colData(resDELeuk)[, c(1,2,3)]
RLEleuk <- SummarizedExperiment::assays(resDELeuk)##$rle
DESeq2leuk <- S4Vectors::metadata(resDELeuk)$DESeq2obj[c(4, 5, 6, 7)]
metaDEleuk <- S4Vectors::metadata(resDELeuk)[-c(1, 2)]

SummarizedExperiment::colData(resDELeukSUB) <- colDataLeuk
SummarizedExperiment::assays(resDELeukSUB) <- RLEleuk##list(rle=RLEleuk)
S4Vectors::metadata(resDELeukSUB) <- metaDEleuk
S4Vectors::metadata(resDELeukSUB)$DESeq2obj <- DESeq2leuk

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
## names(metadata(resDELeukSUB)$DESeq2obj)
# resDEmus2SUB <- resDEMUS2
# colDatamus2 <- SummarizedExperiment::colData(resDEMUS2)[, c(1,2,3)]
# RLEmus2 <- SummarizedExperiment::assays(resDEMUS2)$rle
# DESeq2mus2 <- S4Vectors::metadata(resDEMUS2)$DESeq2obj[c(5, 6, 7)]
# metaDEmus2 <- S4Vectors::metadata(resDEMUS2)[-c(1, 2, 3)]
#
# SummarizedExperiment::colData(resDEmus2SUB) <- colDatamus2
# SummarizedExperiment::assays(resDEmus2SUB) <- list(rle=RLEmus2)
# S4Vectors::metadata(resDEmus2SUB) <- metaDEmus2
# S4Vectors::metadata(resDEmus2SUB)$DESeq2obj <- DESeq2mus2

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
Results_DEanalysis_sub500<-list(DE_Antoszewski2022_MOUSEsub500=resDEmus1SUB,
                                DE_Leong2014_FISSIONsub500wt=resDEfissionSUB,
                                DE_Schleiss2021_CLLsub500=resDELeukSUB)
## DE_Weger2021_MOUSEsub500=resDEmus2SUB

# usethis::use_data(Results_DEanalysis_sub500, overwrite=TRUE)
# usethis::use_r("Results_DEanalysis_sub500") # pour creer une fonction
##---------------------------------------------------------------------------##


##===========================================================================##
## Database
##===========================================================================##

DirDataExpr <- file.path("/Users/rodolphe/Documents/THESE-ATER",
                         "20181001-20230831_THESE/THESE_AUTRE/THESE_HELP_VALLAT",
                         "GeneInfo_HGNC.NCBI.ENSEMBL/Final_Pattern_DATA")
NameDatabase <- paste("231102_MergeDatabase_selectionLV2", ".csv", sep="")

Database <- read.csv2(file=file.path(DirDataExpr, NameDatabase),
                      header=TRUE, sep=";")
## str(Database)
##---------------------------------------------------------------------------##


Database2 <- Database[, c(2, 13, 5)]
## str(Database2)
## length(intersect(which(Database2$locus_group == "protein-coding gene"),
##                  which(is.na(Database2[,2]))))
## table(Database2$locus_group[which(!is.na(Database2[,2]))])
## SELrow <- intersect(which(Database2$locus_group == "protein-coding gene"),
##                     which(!is.na(Database2[,2])))

SELrow <- which(!is.na(Database2[, 2]))
Transcript_HomoSapiens_Database <- Database2[SELrow, seq_len(2)]

## object.size(Transcript_HomoSapiens_Database)/(8*1000)
## usethis::use_data(Transcript_HomoSapiens_Database, overwrite=TRUE)
## usethis::use_r("Transcript_HomoSapiens_Database") ## pour creer une fonction

