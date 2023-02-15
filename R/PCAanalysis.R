#' @title Automatic PCA analysis (Main function)
#'
#' @description The functions performs an automatic principal component
#' analysis (PCA) from a gene expression dataset where samples can belong to
#' different biological conditions and/or time points.
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
#' @param Mean.Accross.Time TRUE or FALSE.
#' If FALSE and if \code{Time.position} is not set as NULL,
#' consecutive time points within a sample are linked to help visualization of
#' temporal patterns.
#' If TRUE and if \code{Time.position} is not set as NULL,
#' the mean per time of all genes is computed for each biological condition and
#' the means of consecutive time points within biological condition are linked
#' to help visualization of temporal patterns.
#' @param gene.deletion NULL or a vector of characters or a vector of integers.
#' If \code{gene.deletion} is a vector of characters, all genes with names in
#' \code{gene.deletion} will be deleted from \code{ExprData}.
#' If \code{gene.deletion} is a vector of integers,
#' all the corresponding row numbers of \code{ExprData} will be deleted.
#' If \code{gene.deletion=NULL} all genes of \code{ExprData} will be used in
#' the construction of the PCA.
#' @param sample.deletion \code{NULL} or a vector of characters or
#' a vector of integers.
#' If \code{sample.deletion} is a vector of characters, all samples with names
#' in \code{sample.deletion} will not be used in the construction of the PCA.
#' If \code{sample.deletion} is a vector of integers, all the corresponding
#' column numbers of \code{ExprData} will not be used in the construction
#' of the PCA.
#' If \code{sample.deletion=NULL} all samples will be used in the construction
#' of the PCA.
#' @param Supp.del.sample \code{TRUE} or \code{FALSE}.
#' If \code{FALSE}, the samples selected with \code{sample.deletion}
#' will be deleted.
#' If \code{TRUE}, the samples selected with \code{sample.deletion}
#' will be plotted.
#' These individuals are called supplementary individuals in [FactoMineR::PCA()].
#' @param Color.Group NULL or a data.frame with \eqn{N_{bc}} rows and
#' two columns where \eqn{N_{bc}} is the number of biological conditions.
#' If \code{Color.Group} is a data.frame, the first column must contain
#' the name of each biological condition and the second column must contain
#' the colors associated to each biological condition.
#' If \code{Color.Group=NULL}, the function will automatically attribute
#' a color for each biological condition.
#' If samples belong to different time points only,
#' \code{Color.Group} will not be used.
#' @param D3.mouvement \code{TRUE} or \code{FALSE}.
#' If \code{TRUE}, the 3D PCA plots will also be plotted in a rgl window
#' (see [plot3Drgl::plotrgl()]) allowing to interactively rotate and zoom.
#' @param Phi Angle defining the colatitude direction for the 3D PCA plot
#' (see \code{Details} in [graphics::persp()]).
#' @param Theta Angle defining the azimuthal direction for the 3D PCA plot
#' (see \code{Details} in [graphics::persp()]).
#' @param Cex.point Non negative numeric value giving the size of points
#' in all PCA plots which are not automatically plotted by [FactoMineR::PCA()].
#' @param Cex.label Non negative numeric value giving the size of
#' the labels associated to each point of the all PCA graphs which are not
#' automatically plotted by [FactoMineR::PCA()].
#' @param epsilon Non negative numeric value giving the length between points
#' and their labels in all PCA plots which are not automatically plotted
#' by [FactoMineR::PCA()].
#' @param path.result Character or \code{NULL}. Path to save all results.
#' If \code{path.result} contains a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}" and a sub sub folder,
#' "1-2_PCAanalysis_\code{Name.folder.pca}"
#' all results will be saved in the sub folder
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}/1-2_PCAanalysis_\code{Name.folder.pca}".
#' Otherwise, a sub folder entitled
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}" and/or a sub sub folder
#' "1-2_PCAanalysis_\code{Name.folder.pca}"
#' will be created in \code{path.result} and all results will be saved in
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}/1-2_PCAanalysis_\code{Name.folder.pca}".
#' If \code{NULL}, the results will not be saved in a folder.
#' \code{NULL} as default.
#' @param Name.folder.pca Character or \code{NULL}.
#' If \code{Name.folder.pca} is a character, the folder and sub folder names
#' which will contain the PCA graphs will respectively be
#' "1_UnsupervisedAnalysis_\code{Name.folder.pca}"
#' and "1-2_PCAanalysis_\code{Name.folder.pca}".
#' Otherwise, the folder and sub folder names will respectively be
#' "1_UnsupervisedAnalysis" and "1-2_PCAanalysis".
#'
#' @return The function returns several 2D and 3D PCA graphs depending on
#' the experimental design
#' * When samples belong only to different biological conditions,
#' the function returns a 2D and two 3D PCA graphs.
#' In each graph, samples are colored with different colors for different
#' biological conditions. The two 3D PCA graphs are identical but one of them
#' will be opened in a rgl window (see [plot3Drgl::plotrgl()]) and it allows to
#' interactively rotate and zoom.
#' * When samples belong only to different time points, the function returns
#'   * One 2D PCA graph, one 3D PCA graph and the same 3D PCA graph
#'   in a rgl window where samples are colored with different colors
#'   for different time points.
#'   Furthermore, lines are drawn between each pair of consecutive points
#'   for each sample (if \code{Mean.Accross.Time=FALSE},
#'   otherwise it will be only between means).
#'   * The same graphs describe above but without lines.
#' * When samples belong to different time points and
#' different biological conditions, the function returns
#'   * One 2D PCA graph, one 3D PCA graph and the same 3D PCA graph
#'   in a rgl window where samples are colored with different colors
#'   for different time points.
#'   Furthermore, lines are drawn between each pair of consecutive points
#'   for each sample (if \code{Mean.Accross.Time=FALSE},
#'   otherwise it will be only between means).
#'   * The same graphs describe above but without lines.
#'   * The same six following graphs for each biological condition
#'   (one PCA analysis per biological condition).
#'   One 2D PCA graph, one 3D PCA graph and the same 3D PCA graph
#'   in a rgl window where samples belong to only one biological condition and
#'   are colored with different colors for different time points.
#'   Furthermore, lines are drawn between each pair of consecutive points
#'   for each sample (if \code{Mean.Accross.Time=FALSE},
#'   otherwise it will be only between means).
#'   The three others graphs are identical to the three previous ones
#'   but without lines.
#'
#' The interactive 3D graphs will be plotted only if \code{D3.mouvement=TRUE}.
#'
#' @seealso The function calls the R functions [PCAgraphics()]
#' and [ColnamesToFactors()].
#'
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#' res.sim.count<-RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=2,
#'                                    Nb.Gene=20)
#' #--------------------------------------------------------------------------#
#' Res.acp.analysis<-PCAanalysis(ExprData=res.sim.count$Sim.dat,
#'                               Column.gene=1,
#'                               Group.position=1,
#'                               Time.position=2,
#'                               Individual.position=3,
#'                               Mean.Accross.Time=FALSE,
#'                               gene.deletion=NULL,
#'                               sample.deletion=NULL,
#'                               Supp.del.sample=FALSE,
#'                               Color.Group=NULL,
#'                               Phi=25, Theta=140, Cex.label=0.7,Cex.point=0.7,
#'                               epsilon=0.2,
#'                               D3.mouvement=FALSE,
#'                               path.result=NULL,
#'                               Name.folder.pca=NULL)

PCAanalysis<-function(ExprData,
                      Column.gene,
                      Group.position,
                      Time.position,
                      Individual.position,
                      Mean.Accross.Time,
                      gene.deletion,
                      sample.deletion,
                      Supp.del.sample=FALSE,
                      Color.Group=NULL,
                      Phi,Theta, Cex.point, Cex.label, epsilon,
                      D3.mouvement=FALSE,
                      path.result=NULL,
                      Name.folder.pca=NULL){
  #---------------------------------------------------------------------------#
  # Preprocessing
  res.Factors<-ColnamesToFactors(ExprData=ExprData,
                                 Column.gene=Column.gene,
                                 Group.position=Group.position,
                                 Time.position=Time.position,
                                 Individual.position=Individual.position)
  Vector.group<-res.Factors$Group.Info
  Vector.time<-res.Factors$Time.Info
  #---------------------------------------------------------------------------#
  # Folder creation if no existence
  #---------------------------------------------------------------------------#
  if(is.null(Name.folder.pca)==TRUE){
    Name.folder.pca<-""
    SubFolder.name<-"1_UnsupervisedAnalysis"
  }else{
    Name.folder.pca<-paste("_",Name.folder.pca,sep="")
    SubFolder.name<-paste("1_UnsupervisedAnalysis",Name.folder.pca,sep="")
  }# if(is.null(Name.folder.pca)==TRUE)
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
  #---------------------------------------------------------------------------#
  if(is.null(Vector.group)==FALSE){
    if(is.null(Vector.time)==FALSE){
      Name.file.pca<-paste("BothGroupTime",Name.folder.pca,sep="")
    }else{
      Name.file.pca<-paste("Group",Name.folder.pca,sep="")
    }# if(is.null(Vector.time)==FALSE)
  }else{
    Name.file.pca<-paste("Time",Name.folder.pca,sep="")
  }# if(is.null(Vector.group)==FALSE)
  #
  if(is.null(path.result.f)==FALSE){
    nom.dossier.result<-paste("1-2_PCAanalysis",Name.folder.pca,sep="")
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
  # Main results
  #---------------------------------------------------------------------------#
  tb.spinfoini<-as.numeric(table(res.Factors$Individual.info))
  max.tb<-max(tb.spinfoini)
  Var.sample<-stats::var(tb.spinfoini)
  if(Mean.Accross.Time==FALSE & Var.sample==0 & max.tb>1){#tb.spinfoini[1]>1
    res.PCA<-PCAgraphics(ExprData=ExprData,
                         Column.gene=Column.gene,
                         Group.position=Group.position,
                         Time.position=Time.position,
                         Individual.position=Individual.position,
                         Mean.Accross.Time=FALSE,
                         sample.deletion=sample.deletion,
                         Supp.del.sample=Supp.del.sample,
                         gene.deletion=gene.deletion,
                         Color.Group=Color.Group,
                         D3.mouvement=D3.mouvement,
                         Phi=Phi,Theta=Theta, epsilon=epsilon,
                         Cex.point=Cex.point, Cex.label=Cex.label,
                         path.result=path.result.new,
                         Name.file.pca=Name.file.pca)
  }else{
    res.PCA<-PCAgraphics(ExprData=ExprData,
                         Column.gene=Column.gene,
                         Group.position=Group.position,
                         Time.position=Time.position,
                         Individual.position=Individual.position,
                         Mean.Accross.Time=TRUE,
                         sample.deletion=sample.deletion,
                         Supp.del.sample=Supp.del.sample,
                         gene.deletion=gene.deletion,
                         Color.Group=Color.Group,
                         D3.mouvement=D3.mouvement,
                         Phi=Phi,Theta=Theta, epsilon=epsilon,
                         Cex.point=Cex.point, Cex.label=Cex.label,
                         path.result=path.result.new,
                         Name.file.pca=Name.file.pca)
  }# if(Mean.Accross.Time==FALSE & Var.sample==0 & max.tb>1)
  #---------------------------------------------------------------------------#
  if(is.null(Vector.group)==FALSE & is.null(Vector.time)==FALSE){
    Tt.Del<-gsub("t","", gsub("T","",as.character(Vector.time)))
    Time.name<-levels(as.factor(paste("T",Tt.Del,sep="")))
    #-------------------------------------------------------------------------#
    Group.Levels<-levels(as.factor(Vector.group))
    res.PCA.per.g<-vector(mode="list", length=length(Group.Levels))
    names(res.PCA.per.g)<-paste("PCA.results.Group_", Group.Levels, sep="")
    #
    for(g in seq_len(length(Group.Levels))){
      Index.g<-which(Vector.group==Group.Levels[g])
      #
      if(is.null(Column.gene)==TRUE){
        Sub.data.g<-ExprData[,Index.g]
        Index.g.f<-Index.g
      }else{
        Sub.data.g<-cbind(ExprData[,Column.gene],
                          ExprData[,-Column.gene][,Index.g])
      }# if(is.null(Column.gene)==TRUE)
      #-----------------------------------------------------------------------#
      Name.file.pca.g<-paste("Group_", Group.Levels[g], sep="")
      #-----------------------------------------------------------------------#
      tb.spinfoini<-as.numeric(table(res.Factors$Individual.info))
      Var.sample<-stats::var(tb.spinfoini)
      #-----------------------------------------------------------------------#
      if(Mean.Accross.Time==FALSE & Var.sample==0 & tb.spinfoini[1]>1){
        res.PCA.g<-PCAgraphics(ExprData=Sub.data.g,
                               Column.gene=Column.gene,
                               Group.position=NULL,
                               Time.position=Time.position,
                               Individual.position=Individual.position,
                               Mean.Accross.Time=FALSE,
                               sample.deletion=sample.deletion,
                               Supp.del.sample=Supp.del.sample,
                               gene.deletion=gene.deletion,
                               Color.Group=Color.Group,
                               D3.mouvement=D3.mouvement,
                               Phi=Phi, Theta=Theta, epsilon=epsilon,
                               Cex.point=Cex.point, Cex.label=Cex.label,
                               path.result=path.result.new,
                               Name.file.pca=Name.file.pca.g)
      }else{
        res.PCA.g<-PCAgraphics(ExprData=Sub.data.g,
                               Column.gene=Column.gene,
                               Group.position=NULL,
                               Time.position=Time.position,
                               Individual.position=Individual.position,
                               Mean.Accross.Time=TRUE,
                               sample.deletion=sample.deletion,
                               Supp.del.sample=Supp.del.sample,
                               gene.deletion=gene.deletion,
                               Color.Group=Color.Group,
                               D3.mouvement=D3.mouvement,
                               Phi=Phi, Theta=Theta, epsilon=epsilon,
                               Cex.point=Cex.point, Cex.label=Cex.label,
                               path.result=path.result.new,
                               Name.file.pca=Name.file.pca.g)
      }# if(Mean.Accross.Time==FALSE & Var.sample==0 & tb.spinfoini[1]>1)
      names(res.PCA.g)<-paste(names(res.PCA.g),".Group_",Group.Levels[g],sep="")
      res.PCA.per.g[[g]]<-res.PCA.g
    }# for(g in 1:length(Group.Levels))
    return(list(res.pca=res.PCA$res.pca,
                List.plot.PCA=res.PCA$List.plot.PCA,
                PCA.results.per.Group=res.PCA.per.g))
  }else{
    return(list(res.pca=res.PCA$res.pca,
                List.plot.PCA=res.PCA$List.plot.PCA))
  }# if(is.null(Vector.group)==FALSE & is.null(Vector.time)==FALSE)
}# PCAanalysis()
