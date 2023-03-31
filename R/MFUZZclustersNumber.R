#' @title Automatic choice of the number of clusters to use for
#' the Mfuzz analysis
#'
#' @description The function uses [stats::kmeans()] or [FactoMineR::HCPC()]
#' in order to compute the number of cluster for the [Mfuzz::mfuzz()] analysis.
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
#' The kmeans method or the hierarchical clustering method,
#' respectively included in [stats::kmeans()] and [FactoMineR::HCPC()],
#' is used in order to compute the optimal number of clusters.
#' If there are several biological conditions, the algorithm computes
#' one optimal number of clusters per biological condition.
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
#' @param Individual.position Integer indicating the position of the name
#' of the individual (e.g patient, replicate, mouse, yeasts culture ...)
#' in the string of characters in each sample names (see \code{Details}).
#' The names of different individuals must be all different.
#' Furthermore, if individual names are just numbers, they will be transform in
#' a vector of class "character" by [CharacterNumbers()] and
#' a "r" will be added to each individual name ("r" for replicate).
#' @param Method "kmeans" or "hcpc". The method used for selecting the number
#' of cluster to be used for the temporal cluster analysis (see \code{Details}).
#' \code{Method="kmeans"} is advised for large number of genes.
#' @param Max.clust Integer strictly superior to 1 indicating
#' the maximum number of clusters. The default is \code{Max.clust=10}.
#' @param Min.std Numeric positive value. All genes where their
#' standard deviations are smaller than the threshold Min.std will be excluded.
#' @param Plot.Cluster \code{TRUE} or \code{FALSE}. \code{TRUE} as default.
#' If \code{TRUE}, the output graph will be plotted.
#' Otherwise the graph will be plotted.
#' @param path.result Character or \code{NULL}.
#' Path to save the plot described in the section \code{Value}.
#' If \code{NULL}, the graph will not be saved in a folder.
#' \code{NULL} as default.
#'
#' @return The function returns
#' * the optimal number of clusters for each biological condition
#' (between 2 and \code{Max.clust}).
#' * a data.frame with (\eqn{N_{bc}+1}) columns and \code{Max.clust} rows with
#' \eqn{N_{bc}} the number of biological conditions.
#'   * If \code{Method="kmeans"}, the ith rows and the jth column correspond to
#'     the within-cluster intertia (see \code{tot.withinss} from
#'     [stats::kmeans()]) dividing by the sum of the variance of each row of
#'     \code{ExprData} of the (j-1)th biological condition computed by
#'     [stats::kmeans()] with i clusters.
#'     When there is only one cluster, the within-cluster intertia corresponds
#'     to the sum of the variance of each row of \code{ExprData}
#'     (see \code{Details}).
#'     The first column contains integers between 1 and \code{Max.clust} which
#'     corresponds to the number of clusters selected for the [stats::kmeans()]
#'     analysis.
#'   * If \code{Method="hcpc"}, the jth column correspond to the clustering
#'   heights (see the output \code{height} from [FactoMineR::HCPC()])
#'   dividing by the maximum value of \code{height}.
#'   The first column contains integers between 1 and \code{Max.clust} which
#'   corresponds to the number of clusters selected for the [stats::kmeans()]
#'   analysis.
#' * a plot which gives
#'   * If \code{Method="kmeans"}, the evolution of the weighted within-cluster
#'   intertia per number of clusters (from 1 to \code{Max.clust})
#'   for each biological condition.
#'   The optimal number of cluster for each biological condition
#'   will be colored in blue.
#'   * If \code{Method="hcpc"}, the evolution of the scaled height per
#'   number of clusters (from 1 to \code{Max.clust})
#'   for each biological condition.
#'   The optimal number of cluster for each biological condition will be
#'   colored in blue.
#'
#' @seealso The function is called by [MFUZZanalysis()].
#'
#' @importFrom stats kmeans
#' @importFrom FactoMineR HCPC
#' @importFrom graphics lines legend
#' @importFrom grDevices pdf dev.off
#'
#' @export
#'
#' @examples
#' # Data simulation
#' set.seed(33)
#' data.clust.sim<-matrix(rnorm(12*10*3, sd=0.2,
#'                              mean=rep(c(rep(c(1,6,9,4,3,1,6.5,0.7,10),
#'                                             times=2),
#'                                         rep(c(2,3.6,3.7,5,7.9,8,7.5,3.5,3.4),
#'                                         times=2)),
#'                              each=10)),
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
#'                   type=c("b"), pch=19, col=rep(c("black","red"), each=6*10),
#'                   xlab="Time", ylab="Gene expression")
#' #--------------------------------------------------------------------------#
#' MFUZZclustersNumber(ExprData=data.clust.sim,
#'                     Column.gene=NULL,
#'                     Group.position=1,
#'                     Time.position=2,
#'                     Individual.position=3,
#'                     Method="hcpc",
#'                     Max.clust=5,
#'                     Plot.Cluster=TRUE,
#'                     path.result=NULL)

MFUZZclustersNumber<-function(ExprData,
                              Column.gene,
                              Group.position,
                              Time.position,
                              Individual.position,
                              Method="hcpc",
                              Max.clust,
                              Min.std=0.1,
                              Plot.Cluster=TRUE,
                              path.result=NULL){
  #---------------------------------------------------------------------------#
  # Condition on Max.clust
  if(Max.clust<2 | floor(Max.clust)!=Max.clust | Max.clust>=nrow(ExprData)){
    stop("'Max.clust' must be an integer greater or equal to 2.")
  }# if(Max.clust<2 | floor(Max.clust)!=Max.clust)
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
  #
  if(is.null(Vect.group)==TRUE){
    Nb.group<-1
  }else{
    Nb.group<-length(levels(as.factor(Vect.group)))
  }# if(is.null(Vect.group)==TRUE)
  #---------------------------------------------------------------------------#
  # data with only expression
  if(is.null(Column.gene)==TRUE){
    ind.col.expr<-seq_len(ncol(ExprData))
  }else{
    ind.col.expr<-seq_len(ncol(ExprData))[-Column.gene]
  }# if(is.null(Column.gene)==TRUE)
  ExprData.f<-ExprData[,ind.col.expr]
  #---------------------------------------------------------------------------#
  # Data for Mfuzz analysis and selection of the number of cluster
  Data.mfuzz<-matrix(NA, ncol=Nb.time*Nb.group, nrow=nrow(ExprData))
  row.names(Data.mfuzz)<-row.names(ExprData)
  #
  Tps.info<-paste("t", gsub("T","", gsub("t","",levels(as.factor(Vect.time)))),
                  sep="")
  Grp.info<-levels(as.factor(Vect.group))
  if(is.null(Vect.group)==TRUE){
    colname.grp<-""
  }else{
    colname.grp<-paste(".", rep(Grp.info,each=Nb.time), sep="")
  }# if(is.null(Vect.group)==TRUE)
  colnames(Data.mfuzz)<-paste("Mean_",rep(Tps.info,times=Nb.group),colname.grp,
                              sep="")
  #---------------------------------------------------------------------------#
  # Filling the data
  for(g in seq_len(Nb.group)){# 1:Nb.group
    for(t in seq_len(Nb.time)){# 1:Nb.time
      if(is.null(Vect.group)==TRUE){
        Index.t<-which(Vect.time==levels(as.factor(Vect.time))[t])
        Data.mfuzz[,t]<-apply(as.data.frame(ExprData.f[,Index.t]),1,mean)
      }else{
        Index.g<-which(Vect.group==levels(as.factor(Vect.group))[g])
        Index.t<-which(Vect.time==levels(as.factor(Vect.time))[t])
        ExprData.f.TG<-as.data.frame(ExprData.f[,intersect(Index.t,Index.g)])
        Data.mfuzz[,Nb.time*(g-1)+t]<-apply(ExprData.f.TG, 1, FUN=mean)
      }# if(is.null(Vect.group)==TRUE)
    }# for(t in 1:Nb.time)
  }# for(g in 1:Nb.group)
  #---------------------------------------------------------------------------#
  # Data which will contain the results of Kmeans
  Sum.nb.c<-data.frame(matrix(NA, nrow=Max.clust, ncol=Nb.group+1))
  Nb.time<-length(levels(as.factor(Vect.time)))
  if(Method=="hcpc"){Score<-"Tot.withinss.scaled"}else{Score<-"Scaled.height"}
  if(is.null(Vect.group)==TRUE){
    colnames(Sum.nb.c)<-c("Nb.clust", as.character(Score))
  }else{
    colnames(Sum.nb.c)<-c("Nb.clust",
                          paste(as.character(Score), "_",
                                levels(as.factor(Vect.group)), sep=""))
  }# if(is.null(Vect.group)==TRUE)
  Sum.nb.c[,1]<-c(1, 2:Max.clust)
  #---------------------------------------------------------------------------#
  # Kmeans
  Opti.clust<-rep(NA, times=Nb.group)
  col.c<-rep("black", times=Max.clust*Nb.group)
  pch.c<-rep(3, times=Max.clust*Nb.group)
  #
  for(g in seq_len(Nb.group)){
    Std.g<-apply(Data.mfuzz[,seq_len(Nb.time) +Nb.time*(g-1)], 1, sd)
    GeneInf.Min.std.g<-which(Std.g<Min.std)
    if(length(GeneInf.Min.std.g)>0){
      GeneSelNbClust<--GeneInf.Min.std.g
    }else{
      GeneSelNbClust<-seq_len(length(Std.g))
    }# if(length(GeneInf.Min.std.g)>0)
    #
    Data.mfuzz.g<-data.frame(Data.mfuzz[GeneSelNbClust,
                                        seq_len(Nb.time)+Nb.time*(g-1)])
    if(Method=="hcpc"){
      Nb.gene<-nrow(Data.mfuzz.g)
      if(Nb.gene<=200){
        res.pca<-FactoMineR::PCA(round(Data.mfuzz.g,digits=3), graph=FALSE)
        res.hcpc<-FactoMineR::HCPC(res.pca, graph=FALSE,
                                   nb.clust=-1, consol=TRUE, min=2)
      }else{
        res.kk<-NbClustKmeansHCPC(Nb.gene,200,50,30000,175)
        kkHCPC<-res.kk$Nkmeans
        #
        options(warn = -1)
        #
        cl<-stats::kmeans(round(data.frame(scale(Data.mfuzz.g)), digits=2),
                          kkHCPC, iter.max=10)
        res.hcpc<-FactoMineR::HCPC(round(data.frame(cl$centers), digits=2),
                                   graph=FALSE, nb.clust=-1,consol=FALSE,min=3)
        options(warn = 0)
      }# if(nrow(Data.mfuzz.g)<=200)
      #
      Clust.height<-rev(res.hcpc$call$t$tree$height)# rev(res.hc$height)
      Sum.nb.c[,g+1]<-c(Clust.height/max(Clust.height),0)[seq_len(Max.clust)]
      #
      inert.gain<-rev(res.hcpc$call$t$tree$height)
      intra<-rev(cumsum(rev(inert.gain)))
      quot<-intra[2:length(intra)]/intra[seq_len(length(intra)-1)]
      #
      if(abs(which.min(quot)+1-res.hcpc$call$t$nb.clust)>2){
        Index.nb.clust<-res.hcpc$call$t$nb.clust
      }else{
        Index.nb.clust<-max(which.min(quot)+1,res.hcpc$call$t$nb.clust)
        # Index.nb.clust<-res.hcpc$call$t$nb.clust
      }# if(abs(which.min(quot)+1-res.hcpc$call$t$nb.clust)>2)
    }# if(Method=="hcpc")
    #
    if(Method=="kmeans"){
      inertie.wihtin<-rep(0, times=Max.clust-1)
      cpt.clust<-0
      options(warn = -1)
      for(k in seq(from=2, to=Max.clust, by=1)){
        cpt.clust<-cpt.clust+1
        clus<-stats::kmeans(round(data.frame(scale(Data.mfuzz.g)), digits=1),
                            centers=k, nstart=5)
        inertie.wihtin[cpt.clust]<-clus$tot.withinss
      }# for(k in 2:Max.clust)
      options(warn = 0)
      # DQ.within<-rev(cumsum(rev(inertie.wihtin)))
      # quot.DQ.within<-DQ.within[-length(inertie.wihtin)]/DQ.within[-1]
      DQ.within<-rev(cumsum(rev(c(clus$totss, inertie.wihtin))))
      quot.DQ.within<-DQ.within[-1]/DQ.within[-length(inertie.wihtin)]
      Index.nb.clust<-which.min(quot.DQ.within)+1
      #
      Sum.nb.c[,g+1]<-c(clus$totss,inertie.wihtin)/clus$totss
    }# if(Method=="kmeans")
    #
    col.c[Index.nb.clust+Max.clust*(g-1)]<-"blue"
    pch.c[Index.nb.clust+Max.clust*(g-1)]<-19
    Opti.clust[g]<-Index.nb.clust
  }# for(g in 1:Nb.group)
  #---------------------------------------------------------------------------#
  # Save graph
  if(Method=="hcpc"){Ylab<-"Scaled height (ward)"}else{
    Ylab<-"Scaled within-cluster inertia"}
  #
  if(is.null(path.result)==FALSE){
    grDevices::pdf(file=paste(path.result,"/Clustering_OptimalClusterNumber_",
                              paste0(levels(as.factor(Vect.group)),
                                     collapse="_"), ".pdf", sep=""),
                   width=11, height=8)
    #
    plot(Sum.nb.c[,1], Sum.nb.c[,2], type="b", ylim=c(0,1 +0.2*(Nb.group-1)),
         pch=pch.c[seq_len(Max.clust)], col=col.c[seq_len(Max.clust)],
         xlab="Number of cluster", ylab=Ylab)
    #
    if(Nb.group>1){
      for(g in seq_len(Nb.group-1)){
        # Add a second line
        graphics::lines(Sum.nb.c[,1], Sum.nb.c[,g+2]+0.2*g, type="b", lty=g+1,
                        pch=pch.c[seq_len(Max.clust)+ Max.clust*g],
                        col=col.c[seq_len(Max.clust)+ Max.clust*g])
      }# for(g in 1:(Nb.group-1))
      # Add a legend to the plot
      graphics::legend("topright", legend=levels(as.factor(Vect.group)),
                       col=c("black"), lty=seq_len(Nb.group), cex=0.8)
      #,inset=.02
    }# if(Nb.group>1)
    graphics::legend("top", legend=c("Cluster optimal"),
                     col=c("blue"),pch=19,cex=0.7)#,inset=c(0.5,0.02)
    #
    grDevices::dev.off()
  }# if(is.null(path.result)==FALSE)
  #---------------------------------------------------------------------------#
  if(Plot.Cluster==TRUE){
    #
    plot(Sum.nb.c[,1], Sum.nb.c[,2], type="b", ylim=c(0,1 +0.2*(Nb.group-1)),
         pch=pch.c[seq_len(Max.clust)], col=col.c[seq_len(Max.clust)],
         xlab="Number of cluster", ylab=Ylab)
    #
    if(Nb.group>1){
      for(g in seq_len(Nb.group-1)){
        # Add a second line
        graphics::lines(Sum.nb.c[,1], Sum.nb.c[,g+2]+0.2*g, type="b", lty=g+1,
                        pch=pch.c[seq_len(Max.clust)+ Max.clust*g],
                        col=col.c[seq_len(Max.clust)+ Max.clust*g])
      }# for(g in 1:(Nb.group-1))
      # Add a legend to the plot
      graphics::legend("topright", legend=levels(as.factor(Vect.group)),
                       col=c("black"), lty=seq_len(Nb.group), cex=0.8)
      #,inset=.02
    }# if(Nb.group>1)
    graphics::legend("top", legend=c("Cluster optimal"),
                     col=c("blue"),pch=19,cex=0.7)#,inset=c(0.5,0.02)
    #
    # Number.Cluster.plot<-grDevices::recordPlot()
    # graphics::plot.new() ## clean up device
  }# if(Plot.Cluster==TRUE)
  #---------------------------------------------------------------------------#
  # Data containing the number of cluster for each group
  if(is.null(Vect.group)==TRUE){
    data.clust.kmeans<-data.frame(Name="OneGroupOnly",ClusterKmeans=Opti.clust)
  }else{
    data.clust.kmeans<-data.frame(Name=levels(as.factor(Vect.group)),
                                  ClusterKmeans=Opti.clust)
  }# if(is.null(Vect.group)==TRUE)
  #---------------------------------------------------------------------------#
  #Plot.Number.Cluster=Number.Cluster.plot
  return(list(Summary.Nb.Cluster=Sum.nb.c,
              DataClustSel=data.clust.kmeans))
}# MFUZZclustersNumber()

NbClustKmeansHCPC<-function(NrowData,x1,y1,x2,y2){
  if(x2<=x1 | y2<=y1){
    StopMfuzzMessage<-paste("'x2' must be strictly greater than 'x1' and",
                            "'y2' must be strictly greater than 'y1'", sep="")
    stop(StopMfuzzMessage)
  }# if(x2<=x1 | y2<=y1)
  acoef<-log(y2/y1)/log(x2/x1)
  bcoef<-y2/(x2^acoef)
  NbClustKKhcpc<-ceiling(bcoef*NrowData^acoef)
  return(list(Nkmeans=NbClustKKhcpc, a=acoef, b=bcoef))
}# NbClustKmeansHCPC()
