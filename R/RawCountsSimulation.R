#' @title RNA-seq raw counts data simulation
#'
#' @description The function simulates an in silico RNA-seq raw counts data
#' inspired from the model used in the \code{DESeq2} package.
#' It is used in some examples of other functions.
#'
#' @param Nb.Group Non negative integer.
#' Number of biological condition (minimum 1).
#' @param Nb.Time Non negative integer. Number of time points (minimum 1).
#' @param Nb.per.GT Non negative integer.
#' Number of sample for each condition and time (minimum 1).
#' @param Nb.Gene Non negative integer. Number of genes (minimum 1)
#'
#' @return A simulated RNA-seq raw counts data.
#'
#' @importFrom stats model.matrix rnorm rnbinom rexp runif
#'
#' @export
#'
#' @examples
#' RawCountsSimulation(Nb.Group=3,Nb.Time=5,Nb.per.GT=7,Nb.Gene=50)
#' # RawCountsSimulation(Nb.Group=1,Nb.Time=5,Nb.per.GT=7,Nb.Gene=50)
#' # RawCountsSimulation(Nb.Group=3,Nb.Time=1,Nb.per.GT=7,Nb.Gene=50)

RawCountsSimulation<-function(Nb.Group,Nb.Time,Nb.per.GT,Nb.Gene){
  if(max(Nb.Group,Nb.Time)==1){
    stop("At least two groups or two times are demanded.")
  }# if(max(Nb.Group,Nb.Time)==1)
  #
  Tsim<-factor(rep(rep(paste("t", seq(from=0,to=Nb.Time-1,by=1), sep=""),
                       times=Nb.Group),
                   each=Nb.per.GT))
  Gsim<-factor(rep(paste("G", seq_len(Nb.Group), sep=""),# 1:Nb.Group
                   each=Nb.Time*Nb.per.GT))
  #
  S.pg<-rep(seq_len(Nb.per.GT), times=Nb.Time*Nb.Group)# 1:Nb.per.GT
  S.pg<-S.pg+rep(Nb.per.GT*seq(from=0, to=Nb.Group-1, by=1),# (0:(Nb.Group-1))
                 each=Nb.Time*Nb.per.GT)
  #
  Patient<-factor(paste("Ind", S.pg, sep=""))
  Data.sim<-as.data.frame(matrix(0, nrow=Nb.Gene,
                                 ncol=Nb.per.GT*Nb.Group*Nb.Time))
  row.names(Data.sim)<-NULL#paste("Gene",1:Nb.Gene,sep="")
  #---------------------------------------------------------------------------#
  if(min(Nb.Group,Nb.Time)>1){
    colnames(Data.sim)<-paste(Gsim,"_",Tsim,"_",Patient,sep="")
    Design.M<-matrix(stats::model.matrix(~Tsim+Gsim+Gsim:Tsim),
                     ncol=length(colnames(stats::model.matrix(~Tsim+Gsim+Gsim:Tsim))))
    colnames(Design.M)<-colnames(stats::model.matrix(~Tsim+Gsim+Gsim:Tsim))
  }# if(min(Nb.Group,Nb.Time)>1)
  #
  if(Nb.Group==1){
    colnames(Data.sim)<-paste(Tsim,"_",Patient,sep="")
    Design.M<-matrix(stats::model.matrix(~Tsim),
                     ncol=length(colnames(stats::model.matrix(~Tsim))))
    colnames(Design.M)<-colnames(stats::model.matrix(~Tsim))
  }# if(Nb.Group==1)
  #
  if(Nb.Time==1){
    colnames(Data.sim)<-paste(Gsim,"_",Patient,sep="")
    Design.M<-matrix(stats::model.matrix(~Gsim),
                     ncol=length(colnames(stats::model.matrix(~Gsim))))
    colnames(Design.M)<-colnames(stats::model.matrix(~Gsim))
  }# if(Nb.Time==1)
  #---------------------------------------------------------------------------#
  bio.noise<-stats::runif(n=Nb.Gene, min=0,max=3)
  Beta0<-stats::rnorm(n=Nb.Gene,mean=5,sd=4)
  for(g in seq_len(Nb.Gene)){# 1:Nb.Gene
    Beta.g<-stats::rnorm(n=ncol(Design.M)-1,mean=0,
                         sd=stats::rexp(n=ncol(Design.M)-1, rate=0.3))
    mu.s.tg<-2^as.numeric(Design.M%*%c(Beta0[g],Beta.g))
    Data.sim[g,]<-stats::rnbinom(n=ncol(Data.sim),size=bio.noise[g],mu=mu.s.tg)
  }# for(g in 1:Nb.Gene)
  Data.sim.f<-cbind(Gene=paste("Gene", seq_len(Nb.Gene), sep=""), Data.sim)
  row.names(Data.sim.f)<-paste("Gene", seq_len(Nb.Gene), sep="")
  #---------------------------------------------------------------------------#
  return(list(Sim.dat=Data.sim.f,
              Vect.Group=Gsim,
              Vect.Time=Tsim,
              Vect.Sample=Patient))
}# RawCountsSimulation
