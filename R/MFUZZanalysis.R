#' @title Clustering of temporal patterns (Main function).
#'
#' @description The function performs a soft clustering of temporal patterns
#' based on the fuzzy c-means algorithm using the R package \code{Mfuzz}.
#'
#' @details The column names of \code{ExprData} must be a vector of strings
#' of characters containing
#' * a string of characters (if \eqn{k=1}) which is the label of the column
#' containing gene names.
#' * \eqn{N_s} sample names which must be strings of characters containing
#' at least : the name of the individual (e.g patient, mouse, yeasts culture),
#' its biological condition (if there is at least two) and
#' the time where data have been collected if there is at least two;
#' (must be either 't0', 'T0' or '0' for time 0,
#' 't1', 'T1' or '1' for time 1, ...).
#'
#' All these sample information must be separated by underscores
#' in the sample name. For instance 'CLL_P_t0_r1',
#' corresponds to the patient 'r1' belonging to the biological condition 'P'
#' and where data were collected at time 't0'.
#' I this example, 'CLL' describe the type of cells
#' (here chronic lymphocytic leukemia) and is not used in our analysis.
#'
#' In the string of characters 'CLL_P_t0_r1',
#' 'r1' is localized after the third underscore, so \code{Individual.position=4},
#' 'P' is localized after the first underscore, so \code{Group.position=2} and
#' 't0' is localized after the second underscore, so \code{Time.position=3}.
#'
#' The \code{Mfuzz} package works with datasets where rows correspond to genes
#' and columns correspond to times.
#' If \code{ExprData} contains several replicates per time,
#' the algorithm computes the mean of replicates for each gene
#' before using [Mfuzz::mfuzz()].
#' When there are several biological conditions, the algorithm realizes
#' the [Mfuzz::mfuzz()] analysis for each biological condition.
#'
#' @param ExprData Data.frame with \eqn{N_g} rows and (\eqn{N_{s+k}}) columns,
#' where \eqn{N_g} is the number of genes,
#' \eqn{N_s} is the number of samples and
#' \eqn{k=1} if a column is used to specify gene names, or \eqn{k=0} otherwise.
#' If \eqn{k=1}, the position of the column containing gene names is given
#' by \code{Column.gene}.
#' The data.frame contains numeric values giving gene expressions of each gene
#' in each sample.
#' Gene expressions can be raw counts or normalized raw counts.
#' Column names of the data.frame must describe each sample's information
#' (individual, biological condition and time) and have the structure described
#' in the section \code{Details}.
#' @param Column.gene Integer indicating the column where gene names are given.
#' Set \code{Column.gene=NULL} if there is no such column.
#' @param Group.position Integer indicating the position of group information
#' in the string of characters in each sample names (see \code{Details}).
#' Set \code{Group.position=NULL} if there is only one or no biological
#' information in the string of character in each sample name.
#' @param Time.position Integer indicating the position of time measurement
#' information in the string of characters in each sample names
#' (see \code{Details}).
#' Set \code{Time.position=NULL} if there is only one or no time measurement
#' information in the string of character in each sample name.
#' @param Individual.position Integer indicating the position of the name of
#' the individual (e.g patient, replicate, mouse, yeasts culture ...)
#' in the string of characters in each sample names (see \code{Details}).
#' The names of different individuals must be all different.
#' Furthermore, if individual names are just numbers, they will be transform in
#' a vector of class "character" by [CharacterNumbers()] and
#' a "r" will be added to each individual name ("r" for replicate).
#' @param DataNumberCluster Data.frame or \code{NULL}. \code{NULL} as default.
#' If \code{DataNumberCluster} is a data.frame where the first column contains
#' the name of the biological conditions and the second the number of cluster
#' selected for each biological condition.
#' If \code{DataNumberCluster=NULL}, a number of clusters will be automatically
#' computed for each biological condition (see [MFUZZclustersNumber()]).
#' @param Method "kmeans" or "hcpc". The method used for selecting the number
#' of cluster to be used for the temporal cluster analysis (see \code{Details}).
#' Only used if \code{DataNumberCluster} is not \code{NULL}.
#' @param Max.clust Integer strictly superior to 1 indicating the maximum
#' number of clusters.
#' \code{Max.clust} will be used only if \code{DataNumberCluster=NULL}
#' @param Membership Numeric value between 0 and 1.
#' For each cluster, genes with membership values below the threshold
#' \code{Membership} will not be displayed.
#' The membership values correspond to the probability of gene to belong
#' to each cluster.
#' @param Min.std Numeric positive value.
#' All genes where their standard deviations are smaller than the threshold
#' \code{Min.std} will be excluded.
#' @param path.result Character or \code{NULL}. Path to save all results.
#' If \code{path.result} contains a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.mfuzz}" and a sub sub folder,
#' "1-4_MFUZZanalysis_\code{Name.folder.mfuzz}"
#' all results will be saved in the sub folder
#' "1_UnsupervisedAnalysis_\code{Name.folder.mfuzz}/1-4_MFUZZanalysis_\code{Name.folder.mfuzz}.
#' Otherwise, a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.mfuzz}" and/or a sub sub folder
#' "1-4_MFUZZanalysis_\code{Name.folder.mfuzz}"
#' will be created in \code{path.result} and all results will be saved in
#' "1_UnsupervisedAnalysis_\code{Name.folder.mfuzz}/1-4_MFUZZanalysis_\code{Name.folder.mfuzz}".
#' If \code{NULL}, the results will not be saved in a folder.
#' \code{NULL} as default.
#' @param Name.folder.mfuzz Character or \code{NULL}.
#' If \code{Name.folder.mfuzz} is a character, the folder and sub folder names
#' which will contain the PCA graphs will respectively be
#' "1_UnsupervisedAnalysis_\code{Name.folder.mfuzz}" and
#' "1-4_MFUZZanalysis_\code{Name.folder.mfuzz}".
#' Otherwise, the folder and sub folder names will respectively be
#' "1_UnsupervisedAnalysis" and "1-4_MFUZZanalysis".
#'
#' @return The function returns
#' * the final data used for the Mfuzz analysis (see \code{Details}).
#' * the cluster associated to each gene.
#' * plots generated by [Mfuzz::mfuzz.plot2()] for each biological condition.
#'
#' @seealso The function uses the function [MFUZZclustersNumber()] in order to
#' compute the optimal number of cluster for each biological condition
#' with the kmeans method.
#'
#' @importFrom methods new
#' @importFrom Mfuzz filter.NA fill.NA filter.std standardise mestimate
#' mfuzz mfuzz.plot2
#' @importFrom utils write.table
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics par
#'
#' @export
#'
#' @examples
#' # Data simulation
#' data.clust.sim<-matrix(rnorm(12*10*3,
#'                              mean=rep(c(rep(c(1,6,9,4,3,1,6.5,0.7,10),
#'                                             times=2),
#'                                         rep(c(2,3.6,3.7,5,7.9,8,7.5,3.5,3.4),
#'                                             times=2)),
#'                              each=10), sd=0.2),
#'                        nrow=30, ncol=12)
#' #
#' colnames(data.clust.sim)<-c("G1_t0_r1","G1_t1_r1","G1_t2_r1",
#'                             "G1_t0_r2","G1_t1_r2","G1_t2_r2",
#'                             "G2_t0_r3","G2_t1_r3","G2_t2_r3",
#'                             "G2_t0_r4","G2_t1_r4","G2_t2_r4")
#' #--------------------------------------------------------------------------#
#' # Plot the temporal expression of each individual
#' graphics::matplot(t(rbind(data.clust.sim[,1:3], data.clust.sim[,4:6],
#'                           data.clust.sim[,7:9], data.clust.sim[,10:12])),
#'                           type=c("b"), pch=19, col=rep(c("black","red"),
#'                           each=6*10), xlab="Time", ylab="Gene expression")
#' #--------------------------------------------------------------------------#
#' MFUZZanalysis(ExprData=data.clust.sim,
#'               Column.gene=NULL,
#'               Group.position=1,
#'               Time.position=2,
#'               Individual.position=3,
#'               DataNumberCluster=NULL,
#'               Membership=0.5,
#'               Min.std=0.1,
#'               path.result=NULL)

MFUZZanalysis<-function(ExprData,
                        Column.gene,
                        Group.position,
                        Time.position,
                        Individual.position,
                        DataNumberCluster,
                        Method="hcpc",
                        Max.clust=10,
                        Membership,
                        Min.std=0.1,
                        path.result=NULL,
                        Name.folder.mfuzz=NULL){
  #---------------------------------------------------------------------------#
  # Folder creation if no existence
  #---------------------------------------------------------------------------#
  if(is.null(Name.folder.mfuzz)==TRUE){
    Name.folder.mfuzz<-""
    SubFolder.name<-"1_UnsupervisedAnalysis"
  }else{
    Name.folder.mfuzz<-paste("_", Name.folder.mfuzz, sep="")
    SubFolder.name<-paste("1_UnsupervisedAnalysis", Name.folder.mfuzz, sep="")
  }# if(is.null(Name.folder.mfuzz)==TRUE)
  #
  if(is.null(path.result)==FALSE){
    if(SubFolder.name%in%dir(path = path.result)==FALSE){
      print("Folder creation")
      dir.create(path=paste(path.result,"/",SubFolder.name,sep=""))
      path.result.f<-paste(path.result,"/",SubFolder.name,sep="")
    }else{
      path.result.f<-paste(path.result,"/",SubFolder.name,sep="")
    }# if(SubFolder.name%in%dir(path = path.result)==FALSE)
  }else{
    path.result.f<-NULL
  }# if(is.null(path.result)==FALSE)
  #
  if(is.null(path.result.f)==FALSE){
    nom.dossier.result<-paste("1-4_MFUZZanalysis",Name.folder.mfuzz,sep="")
    if(nom.dossier.result%in%dir(path = path.result.f)==FALSE){
      dir.create(path=paste(path.result.f,"/",nom.dossier.result,sep=""))
      path.result.new<-paste(path.result.f,"/",nom.dossier.result,sep="")
    }else{
      path.result.new<-paste(path.result.f,"/",nom.dossier.result,sep="")
    }# if(nom.dossier.result%in%dir(path = path.result.f)==FALSE)
  }else{
    path.result.new<-NULL
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  # Preprocessing
  res.factors<-ColnamesToFactors(ExprData=ExprData,
                                 Column.gene=Column.gene,
                                 Group.position=Group.position,
                                 Time.position=Time.position,
                                 Individual.position=Individual.position)
  #
  Vect.time<-res.factors$Time.Info
  Vect.group<-res.factors$Group.Info
  Nb.time<-length(levels(as.factor(Vect.time)))
  if(is.null(Vect.group)==TRUE){
    Nb.group<-1
    Name.plot.Mfuzz<-c("ClustersNumbers", "Mfuzz.Plots")
  }else{
    Group.Lvls<-levels(as.factor(Vect.group))
    Nb.group<-length(Group.Lvls)
    Name.plot.Mfuzz<-c("ClustersNumbers",
                       paste("Mfuzz.Plots.Group_", Group.Lvls, sep=""))
  }# if(is.null(Vect.group)==TRUE)
  #---------------------------------------------------------------------------#
  List.plot.Mfuzz<-vector(mode="list", length=Nb.group+1)
  names(List.plot.Mfuzz)<-Name.plot.Mfuzz
  #---------------------------------------------------------------------------#
  # Selection the number of clusters
  if(is.null(DataNumberCluster)==TRUE){
    resClustNumber<-MFUZZclustersNumber(ExprData=ExprData,
                                        Column.gene=Column.gene,
                                        Group.position=Group.position,
                                        Time.position=Time.position,
                                        Individual.position=Individual.position,
                                        Method=Method,
                                        Max.clust=Max.clust,
                                        Min.std=Min.std,
                                        path.result=path.result.new)
    #-------------------------------------------------------------------------#
    DataNumberCluster<-resClustNumber$DataClustSel
    List.plot.Mfuzz[[1]]<-resClustNumber$Plot.Number.Cluster
    #-------------------------------------------------------------------------#
  }else{
    Sum.clust<-sum(DataNumberCluster[,2])
    if(Sum.clust<2 | floor(Sum.clust)!=Sum.clust){
      stop("Max.clust must be an integer greater or equal to 2")
    }# if(Sum.clust<2 | floor(Sum.clust)!=Sum.clust)
  }# if(is.null(DataNumberCluster)==TRUE)
  #---------------------------------------------------------------------------#
  # Name of all genes
  if(is.null(Column.gene)==TRUE){
    if(is.null(row.names(ExprData))==TRUE){
      Name.G<-as.character(seq_len(nrow(ExprData)))
    }else{
      Name.G<-row.names(ExprData)
    }# if(is.null(row.names(ExprData))==TRUE)
  }else{
    Name.G<-ExprData[,Column.gene]
  }# if(is.null(Column.gene)==TRUE)
  #---------------------------------------------------------------------------#
  # data with only expression
  if(is.null(Column.gene)==TRUE){
    ind.col.expr<-seq_len(ncol(ExprData))
  }else{
    ind.col.expr<-seq_len(ncol(ExprData))[-Column.gene]
  }# if(is.null(Column.gene)==TRUE)
  ExprData.f<-ExprData[,ind.col.expr]
  #---------------------------------------------------------------------------#
  Data.mfuzz<-matrix(NA, ncol=Nb.time*Nb.group, nrow=nrow(ExprData))
  row.names(Data.mfuzz)<-Name.G # row.names(ExprData)
  #
  Tps.info<-paste("t",gsub("T","",
                           gsub("t","",levels(as.factor(Vect.time)))),sep="")
  if(is.null(Vect.group)==TRUE){
    colname.grp<-""
  }else{
    colname.grp<-paste(".", rep(Group.Lvls,each=Nb.time), sep="")
  }# if(is.null(Vect.group)==TRUE)
  colnames(Data.mfuzz)<-paste("Mean_", rep(Tps.info,times=Nb.group),
                              colname.grp, sep="")
  #---------------------------------------------------------------------------#
  for(g in seq_len(Nb.group)){
    for(t in seq_len(Nb.time)){
      if(is.null(Vect.group)==TRUE){
        Index.t<-which(Vect.time==levels(as.factor(Vect.time))[t])
        Data.mfuzz[,t]<-apply(as.data.frame(ExprData.f[,Index.t]), 1, mean)
      }else{
        Index.g<-which(Vect.group==Group.Lvls[g])
        Index.t<-which(Vect.time==levels(as.factor(Vect.time))[t])
        ExprData.f.TG<-as.data.frame(ExprData.f[,intersect(Index.t,Index.g)])
        Data.mfuzz[,Nb.time*(g-1)+t]<-apply(ExprData.f.TG, 1, FUN=mean)
      }# if(is.null(Vect.group)==TRUE)
    }# for(t in 1:Nb.time)
  }# for(g in 1:Nb.group)
  #---------------------------------------------------------------------------#
  # data which will contain the Mfuzz results
  dat.mfuzz<-data.frame(matrix(NA, ncol=1+3*Nb.group, nrow=nrow(ExprData)))
  row.names(dat.mfuzz)<-row.names(Data.mfuzz)
  #
  if(is.null(Vect.group)==TRUE){
    colnames(dat.mfuzz)<-c("Gene.Name", "Cluster", "Membership",
                           paste("Cluster.alpha_", Membership, sep=""))
  }else{
    Pref.name<-c("Gene.Name",
                 rep(c("Cluster", "Membership",
                       paste("Cluster.alpha_", Membership, sep="")),
                     times=Nb.group))
    Suff.name<-c("", rep(Group.Lvls,each=3))
    colnames(dat.mfuzz)<-paste(Pref.name ,c("", rep("_",times=3*Nb.group)),
                               Suff.name, sep="")
  }# if(is.null(Vect.group)==TRUE){
  #---------------------------------------------------------------------------#
  # Mfuzz analysis
  Data.mfuzz.excl<-Data.mfuzz
  row.names(Data.mfuzz.excl)<-as.character(seq_len(nrow(Data.mfuzz)))
  #
  for(g in seq_len(Nb.group)){# 1:Nb.group
    eset<-methods::new('ExpressionSet',
                       exprs=Data.mfuzz.excl[,seq_len(Nb.time)+Nb.time*(g-1)])
    eset.r<-Mfuzz::filter.NA(eset, thres=0.25)
    eset.f<-Mfuzz::fill.NA(eset.r, mode="mean")
    eset.fstd<-Mfuzz::filter.std(eset.f, min.std=Min.std, visu=FALSE)
    eset.s<-Mfuzz::standardise(eset.fstd)
    # Estimation of the parameter m
    m.choose<-Mfuzz::mestimate(eset.s)
    # c.choose=cselection(eset,m.choose,crange=seq(4,32,4),repeats=5,visu=TRUE)
    ## cselection impossible when the data is too big
    if(Nb.group>1){
      Id.c.g<-which(DataNumberCluster[,1]==Group.Lvls[g])
    }else{
      Id.c.g<-1
    }# if(Nb.group>1)
    Nb.c.g<-DataNumberCluster[Id.c.g,2]
    cl<-Mfuzz::mfuzz(eset.s, centers=Nb.c.g, m=m.choose)
    #-------------------------------------------------------------------------#
    # Filling the data which will contain the Mfuzz results
    Cluster.alpha<-rep(0, times=length(cl$cluster))
    #
    for(c in seq_len(length(cl$cluster))){# 1:length(cl$cluster)
      if(cl$membership[c,cl$cluster[c]]>=Membership){
        Cluster.alpha[c]<-cl$cluster[c]
      }# if(cl$membership[c,cl$cluster[c]]>=Membership)
    }# for(c in 1:length(cl$cluster))
    #
    gene.not.excluded<-as.numeric(row.names(eset.fstd@assayData$exprs))
    dat.mfuzz[gene.not.excluded,3*(g-1)+1+1]<-as.numeric(cl$cluster)
    dat.mfuzz[gene.not.excluded,3*(g-1)+2+1]<-round(apply(cl$membership,1,max),
                                                    digits=2)
    dat.mfuzz[gene.not.excluded,3*(g-1)+3+1]<-Cluster.alpha
    #-------------------------------------------------------------------------#
    if(is.null(Vect.group)==TRUE){
      Suff.save<-c("")
    }else{
      Suff.save<-paste(Group.Lvls,"_",sep="")
    }# if(is.null(Vect.group)==TRUE)
    #-------------------------------------------------------------------------#
    nrow.plot.f<-ceiling(sqrt(Nb.c.g))
    ncol.plot.f<-ceiling(Nb.c.g/ceiling(sqrt(Nb.c.g)))
    #-------------------------------------------------------------------------#
    # Mfuzz plot
    if(is.null(path.result)==FALSE){
      grDevices::pdf(file=paste(path.result.new,"/Clusters_",
                                Suff.save[g], "alpha", Membership, "_",
                                "Plot_MFUZZ", Name.folder.mfuzz,
                                ".pdf",sep=""),
                     width=11, height=8)# width = 8, height = 11
      plot.mfuzz<-Mfuzz::mfuzz.plot2(eset.s, cl=cl, min.mem=Membership,
                                     centre=TRUE,
                                     time.labels=levels(as.factor(Vect.time)),
                                     mfrow=c(nrow.plot.f, ncol.plot.f),
                                     x11=FALSE)
      grDevices::dev.off()
    }# if(is.null(path.result)==FALSE)
    # graphics::par(mfrow=c(nrow.plot.f,ncol.plot.f))
    Mfuzz::mfuzz.plot2(eset.s, cl=cl, min.mem=Membership, centre=TRUE,
                       time.labels=levels(as.factor(Vect.time)),
                       mfrow=c(nrow.plot.f,ncol.plot.f), x11=FALSE)
    # The following lines ensure good ouput graphics
    if(Nb.c.g>nrow.plot.f*ncol.plot.f){
      mod.gr<-Nb.c.g%%nrow.plot.f*ncol.plot.f
      for(p in seq_len(mod.gr)){# 1:mod.gr
        plot(0, type='n', axes=FALSE, ann=FALSE)
      }# for(p in 1:mod.gr)
    }# if(Nb.c.g >nrow.plot.f*ncol.plot.f)
    #
    if(Nb.c.g<nrow.plot.f*ncol.plot.f){
      diff.gr<-nrow.plot.f*ncol.plot.f-Nb.c.g
      for(p in seq_len(diff.gr)){# 1:diff.gr
        plot(0,type='n',axes=FALSE,ann=FALSE)
      }# for(p in 1:diff.gr)
    }# if(Nb.c.g <nrow.plot.f*ncol.plot.f)
    graphics::par(mfrow=c(1,1))
    #
    Mfuzz.Plot.g<-grDevices::recordPlot()
    graphics::plot.new() ## clean up device
    List.plot.Mfuzz[[g+1]]<-Mfuzz.Plot.g
    #
  }# for(g in 1:Nb.group)
  dat.mfuzz[,1]<-Name.G
  #---------------------------------------------------------------------------#
  # Mfuff.plot<-grDevices::recordPlot()
  # graphics::plot.new() ## clean up device
  #---------------------------------------------------------------------------#
  # csv file with Mfuzz results
  if(is.null(path.result)==FALSE){
    utils::write.table(dat.mfuzz,
                       file=paste(path.result.new, "/Clustering_Results_",
                                  paste0(Group.Lvls, collapse="_"),
                                  "_alpha", Membership, "_MFUZZ",
                                  Name.folder.mfuzz, ".csv",sep=""),
                       sep=";", row.names=FALSE)
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  return(list(Data.Mfuzz=Data.mfuzz,
              Result.Mfuzz=dat.mfuzz,
              Plot.Mfuzz=List.plot.Mfuzz))
}# MFUZZanalysis()
