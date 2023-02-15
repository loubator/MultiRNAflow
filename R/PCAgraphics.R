#' @title Multiple 2D and 3D PCA graphs.
#'
#' @description The function plots 2D and 3D PCA using the function
#' [PCArealization()] which realizes a PCA analysis.
#' This function is called repeatedly by the function [PCAanalysis()]
#' if samples belong to different biological conditions and time points.
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
#' in each sample. Gene expressions can be raw counts or normalized raw counts.
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
#' @param Mean.Accross.Time \code{TRUE} or \code{FALSE}.
#' If \code{FALSE} and if \code{Time.position} is not set as \code{NULL},
#' consecutive time points within a sample are linked to help visualization of
#' temporal patterns.
#' If \code{TRUE} and if \code{Time.position} is not set as \code{NULL},
#' the mean per time of all genes is computed for each biological condition and
#' the means of consecutive time points within biological condition are
#' linked to help visualization of temporal patterns.
#' @param gene.deletion \code{NULL} or a vector of characters or
#' a vector of integers.
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
#' If \code{sample.deletion} is a vector of integers,
#' all the corresponding column numbers of \code{ExprData} will not be used
#' in the construction of the PCA.
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
#' If TRUE, the 3D PCA plots will also be plotted in a rgl window
#' (see [plot3Drgl::plotrgl()]) allowing to interactively rotate and zoom.
#' @param Phi Angle defining the colatitude direction for the 3D PCA plot
#' (see \code{Details} in [graphics::persp()]).
#' @param Theta Angle defining the azimuthal direction for the 3D PCA plot
#' (see \code{Details} in [graphics::persp()]).
#' @param Cex.point Non negative numeric value giving the size of points
#' in all PCA plots which are not automatically plotted by [FactoMineR::PCA()].
#' @param Cex.label Non negative numeric value giving the size of
#' the labels associated to each point of the all PCA graphs
#' which are not automatically plotted by [FactoMineR::PCA()].
#' @param epsilon Non negative numeric value giving the length between points
#' and their labels in all PCA plots which are not automatically plotted
#' by [FactoMineR::PCA()].
#' @param path.result Character or \code{NULL}.
#' Path to save the different PCA graphs.
#' If \code{NULL}, the different PCA graphs will not be saved in a folder.
#' \code{NULL} as default.
#' @param Name.file.pca Character or \code{NULL}.
#' If \code{Name.file.pca} is a character, \code{Name.file.pca} will be added
#' at the beginning of all names of the saved graphs.
#'
#' @return The function returns the outputs from the function
#' [FactoMineR::PCA()] and several 2D and 3D PCA graphs depending
#' on the experimental design
#' * When samples belong only to different biological conditions,
#' the function returns a 2D and two 3D PCA graphs.
#' In each graph, samples are colored with different colors for different
#' biological conditions. The two 3D PCA graphs are identical but
#' one of them will be opened in a rgl window (see [plot3Drgl::plotrgl()]) and
#' it allows to interactively rotate and zoom.
#' * When samples belong only to different time points, the function returns
#'   * One 2D PCA graph, one 3D PCA graph and the same 3D PCA graph
#'   in a rgl window where samples are colored with different colors
#'   for different time points.
#'   Furthermore, lines are drawn between each pair of consecutive points
#'   for each sample (if \code{Mean.Accross.Time=FALSE},
#'   otherwise it will be only between means).
#'   * The same graphs describe above but without lines.
#' * When samples belong to different time points and different
#' biological conditions, the function returns
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
#'   in a rgl window where samples belong to only one biological condition
#'   and are colored with different colors for different time points.
#'   Furthermore, lines are drawn between each pair of consecutive points
#'   for each sample (if \code{Mean.Accross.Time=FALSE},
#'   otherwise it will be only between means).
#'   The three others graphs are identical to the three previous ones
#'   but without lines.
#'
#' The interactive 3D graphs will be plotted only if \code{D3.mouvement=TRUE}.
#'
#' @seealso This function is called by the function [PCAanalysis()].
#'
#' @importFrom scales hue_pal
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats aggregate
#' @importFrom FactoMineR plot.PCA
#' @importFrom plot3D scatter3D text3D
#' @importFrom graphics legend lines text
#' @importFrom grDevices dev.off pdf
#' @importFrom plot3Drgl plotrgl
#'
#' @export
#'
#' @examples
#' res.sim.count<-RawCountsSimulation(Nb.Group=2, Nb.Time=3, Nb.per.GT=2,
#'                                    Nb.Gene=20)
#' #
#' data.color.group=data.frame(Name=c("G1","G2"), Col=c("black","red"))
#' #--------------------------------------------------------------------------#
#' PCAgraphics(ExprData=res.sim.count$Sim.dat,
#'             Column.gene=1,
#'             Group.position=1,
#'             Time.position=2,
#'             Individual.position=3,
#'             Mean.Accross.Time=FALSE,
#'             gene.deletion=c("Gene1","Gene5"),
#'             sample.deletion=c(2,6),
#'             Supp.del.sample=FALSE,
#'             Color.Group=data.color.group,
#'             D3.mouvement=FALSE,
#'             Phi=25, Theta=140, Cex.label=0.7, Cex.point=0.7, epsilon=0.2,
#'             path.result=NULL,
#'             Name.file.pca=NULL)

PCAgraphics<-function(ExprData,
                      Column.gene,
                      Group.position,
                      Time.position,
                      Individual.position,
                      Mean.Accross.Time,
                      gene.deletion,
                      sample.deletion,
                      Supp.del.sample,
                      Color.Group=NULL,
                      D3.mouvement=FALSE,
                      Phi, Theta, Cex.point, Cex.label, epsilon,
                      path.result=NULL,
                      Name.file.pca=NULL){
  #---------------------------------------------------------------------------#
  # PCA
  res.pca.analysis<-PCArealization(ExprData=ExprData,
                                   Column.gene=Column.gene,
                                   Group.position=Group.position,
                                   Time.position=Time.position,
                                   Individual.position=Individual.position,
                                   sample.deletion=sample.deletion,
                                   gene.deletion=gene.deletion,
                                   Supp.del.sample=Supp.del.sample)
  #
  Vector.time<-res.pca.analysis$List.Factors$Vector.time
  Vector.group<-res.pca.analysis$List.Factors$Vector.group
  Vector.patient<-res.pca.analysis$List.Factors$Vector.patient
  res.PCA<-res.pca.analysis$res.pca
  #---------------------------------------------------------------------------#
  # Color assignment
  if(is.null(Time.position)==FALSE & is.null(Group.position)==TRUE){
    Tlevels<-paste("t",
                   gsub("t","",
                        gsub("T","",as.character(levels(factor(Vector.time))))),
                   sep="")
    levels(res.PCA$call$quali.sup$quali.sup$Quali.Sup.Time)<-Tlevels
    #
    Color.Time<-NULL
    if(is.null(Color.Time)==TRUE){
      Color.Time<-data.frame(Name=Tlevels,
                             Col=c("#737373",# "#252525"
                                   scales::hue_pal()(length(Tlevels)-1)))
      #
      data.color<-rbind(Color.Time, rep(NA,times=ncol(Color.Time)))
      data.color<-data.color[-which(is.na(data.color$Col)),]
    }else{
      Id.LevelColT<-order(Color.Time[,1])
      Color.Time<-data.frame(Name=Tlevels, Col=Color.Time[Id.LevelColT,2])
      #
      data.color<-rbind(Color.Time,rep(NA,times=ncol(Color.Time)))
      data.color<-data.color[-which(is.na(data.color$Col)),]
    }# if(is.null(Color.Time)==TRUE)
    #
    ind.color<-res.PCA$call$quali.sup$quali.sup$Quali.Sup.Time
    levels(ind.color)<-data.color[,2]
    LegendTitle<-"Time"
  }else{
    Glevels<-levels(factor(Vector.group))
    #
    if(is.null(Color.Group)==TRUE){
      MypaletteG<-c(RColorBrewer::brewer.pal(8,"Dark2"),
                    RColorBrewer::brewer.pal(8,"Set2"))
      #
      if(length(Glevels)>16){
        MypaletteG<-c(MypaletteG,
                      scales::hue_pal(l=90)(seq_len(length(Glevels)-1)))
      }# if(length(Glevels)>16)
      #
      Color.Group<-data.frame(Name=Glevels,
                              Col=MypaletteG[seq_len(length(Glevels))])
      #
      data.color<-rbind(Color.Group, rep(NA,times=ncol(Color.Group)))
      data.color<-data.color[-which(is.na(data.color$Col)),]
    }else{
      Id.LevelCol.G<-order(Color.Group[,1])
      Color.Group<-data.frame(Name=Glevels, Col=Color.Group[Id.LevelCol.G,2])
      #
      data.color<-rbind(Color.Group, rep(NA,times=ncol(Color.Group)))
      data.color<-data.color[-which(is.na(data.color$Col)),]
    }# if(is.null(Color.Group)==TRUE)
    #
    ind.color<-res.PCA$call$quali.sup$quali.sup$Quali.Sup.Group
    levels(ind.color)<-data.color[,2]
    LegendTitle<-"Group"
  }# if(is.null(Time.position)==FALSE & is.null(Group.position)==TRUE)
  #---------------------------------------------------------------------------#
  if(is.null(Name.file.pca)==TRUE){
    Name.file.pca.f<-"Graph"
  }else{
    Name.file.pca.f<-Name.file.pca
  }# if(is.null(Name.file.pca)==TRUE)
  #---------------------------------------------------------------------------#
  # PCA plot from plot.PCA()
  if(is.null(path.result)==FALSE){
    grDevices::pdf(file = paste(path.result,"/",Name.file.pca.f,"_PCA2D.pdf",
                                sep=""),
                   width = 11, height = 8)#width = 8, height = 11
  }
  options(ggrepel.max.overlaps=30)
  #
  g.2DPCA<-FactoMineR::plot.PCA(res.PCA, axes=c(1,2),
                                choix="ind", habillage="ind",
                                col.hab=as.character(ind.color),
                                col.quali="magenta",
                                shadowtext=TRUE, autoLab="y",
                                graph.type="ggplot", cex=Cex.label)
  # Qualitative factor are in magenta
  print(g.2DPCA)
  options(ggrepel.max.overlaps=10)
  #
  if(is.null(path.result)==FALSE){
    grDevices::dev.off()
  }
  #---------------------------------------------------------------------------#
  # PCAs
  par(mar=c(4.1, 4.1, 2.1, 2.1))#, xpd=TRUE)
  # 3D PCA colored by biological condition or time points
  data.3D<-res.PCA$ind$coord[,c(1,2,3)]
  if(is.null(path.result)==FALSE){
    grDevices::pdf(file=paste(path.result, "/", Name.file.pca.f, "_PCA3D.pdf",
                              sep=""),
                   width=11, height=8)#width = 8, height = 11
  }# if(is.null(path.result)==FALSE)
  #
  # s3d.all<-
  plot3D::scatter3D(data.3D[,1], data.3D[,3], data.3D[,2], pch=20, colvar=NULL,
                    col=as.character(ind.color), cex=Cex.point, bty="b2",
                    ticktype = "detailed", theta=Theta, phi=Phi, d=2,
                    clab=c("",data.color$Name), epsilon=epsilon,
                    xlab=paste("dim1 (",round(res.PCA$eig[,2][1],digits=2),
                               "%)",sep=""),
                    ylab=paste("dim3 (",round(res.PCA$eig[,2][3],digits=2),
                               "%)",sep=""),
                    zlab=paste("dim2 (",round(res.PCA$eig[,2][2],digits=2),
                               "%)",sep=""))
  #
  plot3D::text3D(data.3D[,1]+epsilon, data.3D[,3]+epsilon, data.3D[,2]+epsilon,
                 labels=rownames(data.3D),  add=TRUE,
                 col=as.character(ind.color),
                 colkey = FALSE, cex = Cex.label,#0.5,
                 font=2)
  # graphics::legend("right", title=LegendTitle,
  #                  legend=data.color[,1],
  #                  pch=20, horiz=FALSE, xpd=TRUE,
  #                  cex=0.8, # inset=c(-0.015),
  #                  col=data.color[,2])
  if(is.null(path.result)==FALSE){
    grDevices::dev.off()
  }# if(is.null(path.result)==FALSE)
  #
  s3d.all <- grDevices::recordPlot()
  graphics::plot.new() ## clean up device
  #---------------------------------------------------------------------------#
  # 3D PCA in rgl windows
  if(D3.mouvement==TRUE){
    plot3Drgl::plotrgl()
  }# if(D3.mouvement==TRUE)
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  if(is.null(Vector.time)==FALSE){
    n.row<-length(res.PCA$call$quali.sup$quali.sup$Quali.Sup.Time)
    index.order<-seq_len(n.row)# 1:n.row
    nb.time<-length(unique(res.PCA$call$quali.sup$quali.sup$Quali.Sup.Time))
    #-------------------------------------------------------------------------#
    if(is.null(Vector.group)==FALSE){
      col.link<-as.character(ind.color)
    }else{
      col.link<-rep("grey", times=length(as.character(ind.color)))
    }# if(is.null(Vector.group)==FALSE)
    #-------------------------------------------------------------------------#
    if(is.null(path.result)==FALSE){
      grDevices::pdf(file=paste(path.result, "/", Name.file.pca.f,
                                "_PCAlink2D.pdf",sep=""),
                     width=11, height = 8)#width = 8, height = 11
    }
    # g.2DPCA.t<-
    plot(NA, type="c", col="green", #main="2D PCA",
         xlim=c(min(res.PCA$ind$coord[,1])-epsilon,
                max(res.PCA$ind$coord[,1])+epsilon),
         ylim=c(min(res.PCA$ind$coord[,2])-epsilon,
                max(res.PCA$ind$coord[,2])+epsilon),
         xlab=paste("dim1 (",round(res.PCA$eig[,2][1],digits=2),"%)",sep=""),
         ylab=paste("dim2 (",round(res.PCA$eig[,2][2],digits=2),"%)",sep=""))
    #
    # xlegend<-max(res.PCA$ind$coord[,1])
    #         +epsilon+(max(res.PCA$ind$coord[,1])+epsilon)/8
    #
    # graphics::legend(x=xlegend, y=0,
    #                  title=LegendTitle,
    #                  legend = data.color[,1],
    #                  pch=20, horiz=FALSE, xpd = TRUE,
    #                  cex=0.8, #inset=c(-0.03),
    #                  col=data.color[,2])
    #-------------------------------------------------------------------------#
    NbPerCond<-as.numeric(table(Vector.patient))
    #-------------------------------------------------------------------------#
    #-------------------------------------------------------------------------#
    if(Mean.Accross.Time==FALSE & max(NbPerCond)>1){
      #-----------------------------------------------------------------------#
      for(smpl in seq_len(length(levels(factor(Vector.patient))))){
        # 1:length(levels(factor(Vector.patient)))
        Smpl.sel<-levels(factor(Vector.patient))[smpl]
        times.smpl.sel<-factor(Vector.time[which(Vector.patient==Smpl.sel)])
        coord.smpl<-which(Vector.patient==Smpl.sel)[order(times.smpl.sel)]
        graphics::lines(res.PCA$ind$coord[coord.smpl,1],
                        res.PCA$ind$coord[coord.smpl,2],
                        pch=22,type="c",lty=1,#lty=2,
                        col=col.link[coord.smpl])
        graphics::text(res.PCA$ind$coord[coord.smpl,1]+epsilon,
                       res.PCA$ind$coord[coord.smpl,2]+epsilon,
                       labels=row.names(res.PCA$ind$coord[coord.smpl,]),
                       cex=Cex.label,#0.6,
                       col=as.character(ind.color)[coord.smpl])
      }# for(smpl in 1:length(levels(factor(Vector.patient))))
    }else{
      #-----------------------------------------------------------------------#
      for(smpl in seq_len(length(levels(factor(Vector.patient))))){
        # 1:length(levels(factor(Vector.patient)))
        Smpl.sel<-levels(factor(Vector.patient))[smpl]
        times.smpl.sel<-factor(Vector.time[which(Vector.patient==Smpl.sel)])
        coord.smpl<-which(Vector.patient==Smpl.sel)[order(times.smpl.sel)]
        graphics::points(res.PCA$ind$coord[coord.smpl,1],
                         res.PCA$ind$coord[coord.smpl,2],
                         pch=19,
                         cex=Cex.point,
                         col=col.link[coord.smpl])
        graphics::text(res.PCA$ind$coord[coord.smpl,1]+epsilon,
                       res.PCA$ind$coord[coord.smpl,2]+epsilon,
                       labels=row.names(res.PCA$ind$coord[coord.smpl,]),
                       cex=Cex.label,#0.6,
                       col=as.character(ind.color)[coord.smpl])
      }# for(smpl in 1:length(levels(factor(Vector.patient))))
      #-----------------------------------------------------------------------#
      #-----------------------------------------------------------------------#
      coord.t<-res.PCA$call$quali.sup$quali.sup$Quali.Sup.Time
      data.name<-data.frame(rep("Mean",
                                times=nrow(res.PCA$call$quali.sup$quali.sup)),
                            res.PCA$call$quali.sup$quali.sup)
      #-----------------------------------------------------------------------#
      for(col in seq_len(length(levels(ind.color)))){
        Id.col<-which(as.character(ind.color)==levels(ind.color)[col])
        Mean.t<-stats::aggregate(res.PCA$ind$coord[Id.col,c(1,2)],
                                 list(coord.t[Id.col]), mean)
        New.Name.2D<-apply(unique(data.name[Id.col,]), 1, paste, collapse="_")
        graphics::lines(Mean.t[,2], Mean.t[,3],
                        type="c", pch=22, lwd=2,
                        col=levels(ind.color)[col])
        graphics::text(Mean.t[,2], Mean.t[,3],
                       labels=New.Name.2D,
                       cex=Cex.label*0.9,#0.6,
                       col=levels(ind.color)[col])
      }# for(col in 1:length(levels(ind.color)))
    }# if(Mean.Accross.Time==FALSE & max(NbPerCond)>1)
    #-------------------------------------------------------------------------#
    #-------------------------------------------------------------------------#
    if(is.null(path.result)==FALSE){
      grDevices::dev.off()
    }# if(is.null(path.result)==FALSE)
    #
    g.2DPCA.t<-grDevices::recordPlot()
    graphics::plot.new() ## clean up device
    #-------------------------------------------------------------------------#
    if(is.null(path.result)==FALSE){
      grDevices::pdf(file=paste(path.result, "/", Name.file.pca.f,
                                "_PCAlink3D.pdf", sep=""),
                     width=11, height = 8)#width = 8, height = 11
    }# if(is.null(path.result)==FALSE)
    #-------------------------------------------------------------------------#
    if(Mean.Accross.Time==FALSE & max(NbPerCond)>1){
      Smpl.sel<-levels(factor(Vector.patient))[1]
      times.smpl.sel<-factor(Vector.time[which(Vector.patient==Smpl.sel)])
      coord.smpl<-which(Vector.patient==Smpl.sel)[order(times.smpl.sel)]
      #-----------------------------------------------------------------------#
      if(length(coord.smpl)==1){
        coord.smplf<-c(coord.smpl,NA)
      }else{coord.smplf<-coord.smpl}# if(length(coord.smpl)==1)
      #-----------------------------------------------------------------------#
      # s3d.2t<-
      plot3D::scatter3D(res.PCA$ind$coord[coord.smplf,1],
                        res.PCA$ind$coord[coord.smplf,3],
                        res.PCA$ind$coord[coord.smplf,2],
                        theta=Theta ,phi=Phi, colvar=NULL,
                        type="b", pch=20, bty="b2",#bty = "g"
                        ticktype="detailed",
                        cex=Cex.point,#0.6,
                        col=col.link[coord.smpl], # main = "3D PCA",
                        xlim=c(min(res.PCA$ind$coord[,1]),
                               max(res.PCA$ind$coord[,1])),
                        ylim=c(min(res.PCA$ind$coord[,3]),
                               max(res.PCA$ind$coord[,3])),
                        zlim=c(min(res.PCA$ind$coord[,2]),
                               max(res.PCA$ind$coord[,2])),
                        xlab=paste("dim1 (",round(res.PCA$eig[,2][1],digits=2),
                                   "%)",sep=""),
                        ylab=paste("dim3 (",round(res.PCA$eig[,2][3],digits=2),
                                   "%)",sep=""),
                        zlab=paste("dim2 (",round(res.PCA$eig[,2][2],digits=2),
                                   "%)",sep=""))
      #
      plot3D::text3D(res.PCA$ind$coord[coord.smpl,1],
                     res.PCA$ind$coord[coord.smpl,3],
                     res.PCA$ind$coord[coord.smpl,2],
                     labels=row.names(res.PCA$ind$coord)[coord.smpl],
                     add=TRUE, colkey=FALSE, cex=Cex.label, pch=20,#0.6,
                     col=as.character(ind.color)[coord.smpl])
      #
      # graphics::legend("right", title=LegendTitle,
      #                  legend=data.color[,1],
      #                  pch=20, horiz=FALSE, xpd=TRUE,
      #                  cex=0.8, # inset=c(-0.015),# cex = Cex.point*0.9
      #                  col=data.color[,2])
      #-----------------------------------------------------------------------#
      for(smpl in seq(from=2, to=length(levels(factor(Vector.patient))),by=1)){
        # 2:length(levels(factor(Vector.patient)))
        Smpl.sel<-levels(factor(Vector.patient))[smpl]
        times.smpl.sel<-factor(Vector.time[which(Vector.patient==Smpl.sel)])
        coord.smpl<-which(Vector.patient==Smpl.sel)[order(times.smpl.sel)]
        #---------------------------------------------------------------------#
        if(length(coord.smpl)==1){coord.smplf<-c(coord.smpl,NA)}else{
          coord.smplf<-coord.smpl}
        #---------------------------------------------------------------------#
        plot3D::scatter3D(res.PCA$ind$coord[coord.smplf,1],
                          res.PCA$ind$coord[coord.smplf,3],
                          res.PCA$ind$coord[coord.smplf,2],
                          type="b", colvar=NULL, add=TRUE, pch=20,
                          bty="b2", ticktype="detailed", cex=Cex.point,
                          col=col.link[coord.smpl])
        #0.6, #bty = "g"
        plot3D::text3D(res.PCA$ind$coord[coord.smpl,1],
                       res.PCA$ind$coord[coord.smpl,3],
                       res.PCA$ind$coord[coord.smpl,2],
                       labels=row.names(res.PCA$ind$coord)[coord.smpl],
                       add=TRUE, colkey=FALSE, cex=Cex.label, pch=20, #0.6,
                       col=as.character(ind.color)[coord.smpl])
      }# for(smpl in 2:length(levels(factor(Vector.patient))))
    }else{
      #-----------------------------------------------------------------------#
      Smpl.sel<-levels(factor(Vector.patient))[1]
      times.smpl.sel<-factor(Vector.time[which(Vector.patient==Smpl.sel)])
      coord.smpl<-which(Vector.patient==Smpl.sel)[order(times.smpl.sel)]
      #-----------------------------------------------------------------------#
      if(length(coord.smpl)==1){
        coord.smplf<-c(coord.smpl,NA)
      }else{
        coord.smplf<-coord.smpl}# if(length(coord.smpl)==1)
      #-----------------------------------------------------------------------#
      # s3d.2t<-
      plot3D::scatter3D(res.PCA$ind$coord[coord.smplf,1],
                        res.PCA$ind$coord[coord.smplf,3],
                        res.PCA$ind$coord[coord.smplf,2],
                        theta=Theta ,phi=Phi, colvar=NULL,
                        type="p", pch=20, cex=Cex.point, bty="b2",
                        col=col.link[coord.smpl],# main="3D PCA",#bty="g",#0.6,
                        xlim=c(min(res.PCA$ind$coord[,1]),
                               max(res.PCA$ind$coord[,1])),
                        ylim=c(min(res.PCA$ind$coord[,3]),
                               max(res.PCA$ind$coord[,3])),
                        zlim=c(min(res.PCA$ind$coord[,2]),
                               max(res.PCA$ind$coord[,2])),
                        xlab=paste("dim1 (",round(res.PCA$eig[,2][1],digits=2),
                                   "%)",sep=""),
                        ylab=paste("dim3 (",round(res.PCA$eig[,2][3],digits=2),
                                   "%)",sep=""),
                        zlab=paste("dim2 (",round(res.PCA$eig[,2][2],digits=2),
                                   "%)",sep=""))
      #
      plot3D::text3D(res.PCA$ind$coord[coord.smpl,1],
                     res.PCA$ind$coord[coord.smpl,3],
                     res.PCA$ind$coord[coord.smpl,2],
                     labels=row.names(res.PCA$ind$coord)[coord.smpl],
                     add=TRUE, colkey=FALSE, cex=Cex.label, pch=20, #0.6,
                     col=as.character(ind.color)[coord.smpl])
      #
      # graphics::legend("right", title=LegendTitle,
      #                  legend = data.color[,1],
      #                  pch = 20, horiz=FALSE, xpd = TRUE,
      #                  cex=0.8, # inset=c(-0.015),
      #                  col = data.color[,2])
      #-----------------------------------------------------------------------#
      for(smpl in seq(from=2, to=length(levels(factor(Vector.patient))),by=1)){
        # 2:length(levels(factor(Vector.patient)))
        Smpl.sel<-levels(factor(Vector.patient))[smpl]
        times.smpl.sel<-factor(Vector.time[which(Vector.patient==Smpl.sel)])
        coord.smpl<-which(Vector.patient==Smpl.sel)[order(times.smpl.sel)]
        #---------------------------------------------------------------------#
        if(length(coord.smpl)==1){
          coord.smplf<-c(coord.smpl,NA)
        }else{coord.smplf<-coord.smpl}# if(length(coord.smpl)==1)
        #---------------------------------------------------------------------#
        plot3D::scatter3D(res.PCA$ind$coord[coord.smplf,1],
                          res.PCA$ind$coord[coord.smplf,3],
                          res.PCA$ind$coord[coord.smplf,2],
                          type="p", colvar=NULL, add=TRUE, pch=20,
                          cex=Cex.point, bty="b2",#bty = "g", #0.6,
                          col=col.link[coord.smpl])
        plot3D::text3D(res.PCA$ind$coord[coord.smpl,1],
                       res.PCA$ind$coord[coord.smpl,3],
                       res.PCA$ind$coord[coord.smpl,2],
                       labels=row.names(res.PCA$ind$coord)[coord.smpl],
                       add=TRUE, colkey=FALSE, cex=Cex.label, pch=20,#0.6,
                       col=as.character(ind.color)[coord.smpl])
      }# for(smpl in 2:length(levels(factor(Vector.patient))))
      #-----------------------------------------------------------------------#
      coord.t<-res.PCA$call$quali.sup$quali.sup$Quali.Sup.Time
      data.name<-data.frame(rep("Mean",
                                times=nrow(res.PCA$call$quali.sup$quali.sup)),
                            res.PCA$call$quali.sup$quali.sup)
      #-----------------------------------------------------------------------#
      for(col in seq_len(length(levels(factor(col.link))))){
        # 1:length(levels(factor(col.link)))
        Id.col<-which(col.link==levels(factor(col.link))[col])
        Mean.t<-stats::aggregate(res.PCA$ind$coord[Id.col,c(1,2,3)],
                                 list(coord.t[Id.col]), mean)
        New.Name.3D<-apply(unique(data.name[Id.col,]), 1, paste, collapse="_" )
        #
        plot3D::scatter3D(Mean.t[,2], Mean.t[,4], Mean.t[,3],
                          type="b", colvar=NULL, add=TRUE, pch=20, bty="b2",
                          cex=Cex.point,#0.6, #bty = "g"
                          col=levels(factor(col.link))[col])
        plot3D::text3D(Mean.t[,2], Mean.t[,4], Mean.t[,3],
                       labels=New.Name.3D, col=levels(factor(col.link))[col],
                       add=TRUE, colkey=FALSE, cex=Cex.label, pch=20) #0.6,
      }# for(col in 1:length(levels(factor(col.link))))
    }# if(Mean.Accross.Time==FALSE & max(NbPerCond)>1)
    #-------------------------------------------------------------------------#
    if(is.null(path.result)==FALSE){
      grDevices::dev.off()
    }# if(is.null(path.result)==FALSE)
    s3d.2t<-grDevices::recordPlot()
    graphics::plot.new() ## clean up device
    #-------------------------------------------------------------------------#
    # 3D PCA in rgl window
    if(D3.mouvement==TRUE){
      plot3Drgl::plotrgl()
    }# if(D3.mouvement==TRUE)
    #
    List.plot.PCA<-vector(mode="list", length=4)
    names(List.plot.PCA)<-c("PCA.2D", "PCA.3D",
                            "PCA.2D.temporal.links", "PCA.3D.temporal.links")
    List.plot.PCA[[1]]<-g.2DPCA
    List.plot.PCA[[2]]<-s3d.all
    List.plot.PCA[[3]]<-g.2DPCA.t
    List.plot.PCA[[4]]<-s3d.2t
  }else{
    List.plot.PCA<-vector(mode="list", length=2)
    names(List.plot.PCA)<-c("PCA.2D", "PCA.3D")
    List.plot.PCA[[1]]<-g.2DPCA
    List.plot.PCA[[2]]<-s3d.all
  }# if(is.null(Vector.time)==FALSE)
  #
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  #
  return(list(res.pca=res.PCA,
              List.plot.PCA=List.plot.PCA))
}# PCAgraphics()
