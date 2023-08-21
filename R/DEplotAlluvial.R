#' @title Alluvial graphs of differentially expressed (DE) genes
#'
#' @description
#' The function takes as input a binary table with \eqn{N_g} lines
#' corresponding to genes and
#' * if \code{Temporal.Group=TRUE} : \eqn{T-1} columns corresponding to times
#' (with \eqn{T} the number of time points).
#' A '1' in the n-th row and t-th column means that the n-th gene is
#' differentially expressed (DE) at time t, compared with
#' the reference time t0.
#' * if \code{Temporal.Group=FALSE} :
#' \eqn{G} columns corresponding to the number of group.
#' A '1' in the \eqn{n}-th row and \eqn{g}-th column means
#' that the n-th gene is
#'   * DE at least one time ti, compared with the reference time t0,
#'   for the group \eqn{g}.
#'   * specific at least one time ti, compared with the reference time t0,
#'   for the group \eqn{g} (see [DEanalysisTimeAndGroup()]
#'   for the notion 'specific').
#'   * a signature gene at least one time ti, compared with the reference time
#'   t0, for the group \eqn{g} (see [DEanalysisTimeAndGroup()]
#'   for the notion 'signature').
#'
#' The function plots
#' * if \code{Temporal.Group=TRUE}, two graphs: an alluvial graph and
#' a plot showing the time evolution of the number of DE genes within
#' each temporal group. By temporal group, we mean the sets of genes which
#' are first DE at the same time.
#' * if \code{Temporal.Group=FALSE} : an alluvial graph.
#'
#' @param table.DE.time Binary matrix (table filled with 0 and 1) with
#' \eqn{N_g} rows and \eqn{T-1} columns with \eqn{N_g} the number of genes and
#' \eqn{T-1} the number of time points.
#' @param Temporal.Group \code{TRUE} or \code{FALSE},
#' \code{FALSE} as default (see \code{Description}).
#' @param title.alluvial String of characters or \code{NULL},
#' \code{NULL} as default. The input \code{title.allluvial} corresponds
#' to the title of the alluvial graph.
#' If \code{title} is a string of characters,
#' \code{title} will be the title of the alluvial graph.
#' If \code{title=NULL}, the title of the alluvial graph will be
#' 'Alluvial graph'.
#' @param title.evolution String of characters or \code{NULL},
#' \code{NULL} as default. Only applied if \code{Temporal.Group=TRUE}.
#' The input \code{title.evolution} corresponds to the title of
#' the second graph (see \code{Description}).
#' If \code{title} is a string of characters, it will be to the title of
#' the second graph.
#' If \code{title=NULL}, the title of the second graph will be
#' 'Time evolution of the number of DE genes within each temporal group'.
#'
#' @details The names of the columns of the table will be the axis labels
#' in the plots.
#' If the table has no column names, the function will automatically create
#' column names (t1,t2,...).
#'
#' @return The function returns, as described in \code{description}
#' * if \code{Temporal.Group=TRUE}, two graphs: an alluvial graph and
#' a plot showing the time evolution of the number of DE genes within
#' each temporal group.
#' By temporal group, we mean the sets of genes which are first DE
#' at the same time.
#' * if \code{Temporal.Group=FALSE} : an alluvial graph.
#'
#' @seealso The [DEplotAlluvial()] function
#' * is used by the following functions of our package : [DEanalysisTime()]
#' and [DEanalysisTimeAndGroup()].
#' * calls the R package [ggplot2] in order to plot the two graphs.
#'
#' @importFrom ggalluvial stat_stratum stat_alluvium
#' @importFrom ggplot2 ggplot aes theme scale_x_continuous ggtitle labs guides
#' geom_line geom_point guide_legend geom_area scale_size_manual xlab ylab
#' scale_fill_manual
#' @importFrom stats aggregate
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#' set.seed(1994)
#'
#' NbTime.vst0<-4
#' BinTable<-matrix(sample(c(0,1),replace=TRUE,
#'                         size=NbTime.vst0*120,c(0.60,0.40)),
#'                  ncol=NbTime.vst0)
#' colnames(BinTable)<-paste0("t", 1:NbTime.vst0)
#'
#' ##-------------------------------------------------------------------------#
#' res.alluvial<-DEplotAlluvial(table.DE.time=BinTable)
#' print(res.alluvial$g.alluvial)
#' print(res.alluvial$g.alluvial.freq)

DEplotAlluvial<-function(table.DE.time,
                         Temporal.Group=TRUE,
                         title.alluvial=NULL,
                         title.evolution=NULL){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if(isTRUE(Temporal.Group)){
        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 1) Parameters : column names
        if(is.null(colnames(table.DE.time))){
            Dat.col.name<-paste0("t", seq_len(ncol(table.DE.time)))
        }else{
            Dat.col.name<-colnames(table.DE.time)
        }## if(is.null(colnames(table.DE.time)))
        ## To avoid "no visible binding for global variable" with
        ## devtools::check()
        variable<-Freq<-Gene<-Time_group<-value<-Category<-NULL

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 2) Data for graphics
        ## Only DE genes at least at one times versus t0 are kept
        gene.DE.number<-which(apply(table.DE.time, MARGIN=1,
                                    function(x) sum(x))!=0)

        if(length(gene.DE.number)==1){
            table.DE.time.DE.gene<-t(as.matrix(table.DE.time[gene.DE.number,]))
        }else{
            table.DE.time.DE.gene<-table.DE.time[gene.DE.number,]
        }
        colnames(table.DE.time.DE.gene)<-Dat.col.name

        ##--------------------------------------------------------------------#
        ## Position of DE genes in table.group.gene.DE.peak. Temporal group is
        ## calculated and corresponds to the first time a gene is DE.
        Time.first.peak<-apply(as.matrix(table.DE.time.DE.gene),
                               MARGIN=1,
                               FUN=function(x) which(cumsum(x)==1)[1])

        if(is.null(row.names(table.DE.time.DE.gene))){
            name.g<-as.character(seq_len(length(gene.DE.number)))
        }else{
            name.g<-row.names(table.DE.time.DE.gene)
        }# if(is.null(row.names(table.DE.time.DE.gene)))

        DE.genes.info<-data.frame(Index=gene.DE.number,
                                  Name=name.g,
                                  T.group=as.numeric(Time.first.peak))

        ##--------------------------------------------------------------------#
        ## Data for alluvial graph per gene :
        ## Patterns (1 for DE vs t0), time group & Freq
        data.alluvial<-cbind(table.DE.time.DE.gene,
                             data.frame(Time_group=DE.genes.info$T.group,
                                        Freq=rep(1,times=nrow(table.DE.time.DE.gene))))

        data.alluvial2.t.ini<-stats::aggregate(x=list(Freq=rep(1,nrow(data.alluvial))),
                                               by=data.alluvial[,seq_len(ncol(table.DE.time)+1)],
                                               FUN=length)
        #
        Id.T.0DE<-which(apply(data.alluvial2.t.ini[,seq_len(ncol(table.DE.time))],
                              2, FUN=function(x) 0%in%x))
        Id.T.1DE<-which(apply(data.alluvial2.t.ini[,seq_len(ncol(table.DE.time))],
                              2, FUN=function(x) 1%in%x))
        Fill.col.bar.ini<-rep(c('black','grey'),times=ncol(table.DE.time))
        #
        if(length(Id.T.1DE)==0 | length(Id.T.0DE)==0){
            if(length(Id.T.1DE)==0){
                Keep.col<-2*(c(as.numeric(Id.T.0DE))-1)+1
            }else{
                Keep.col<-2*c(as.numeric(Id.T.1DE))
            }
        }else{
            Keep.col<-sort(c(2*(c(as.numeric(Id.T.0DE))-1)+1,2*c(as.numeric(Id.T.1DE))))
        }# if(length(Id.T.1DE)==0 | length(Id.T.0DE)==0)
        #
        Fill.col.bar<-Fill.col.bar.ini[Keep.col]
        #
        if(length(unique(data.alluvial2.t.ini$Time_group))<=ncol(table.DE.time)){
            add.TG<-cbind(diag(rep(1,times=ncol(table.DE.time))),
                          seq_len(ncol(table.DE.time)),
                          rep(0,times=ncol(table.DE.time)))
            colnames(add.TG)<-colnames(data.alluvial2.t.ini)
            data.alluvial2.t<-rbind(data.alluvial2.t.ini,
                                    add.TG[-sort(unique(data.alluvial2.t.ini$Time_group)),])
        }else{
            data.alluvial2.t<-data.alluvial2.t.ini
        }# if(length(unique(data.alluvial2.t.ini$Time_group))<=ncol(table.DE.time))
        data.alluvial2.tf<-cbind(reshape2::melt(data=cbind(Gene=seq_len(nrow(data.alluvial2.t)),
                                                           data.alluvial2.t),
                                                id.vars=c("Gene","Time_group",
                                                          "Freq")))

        data.alluvial2.tf$Time_group<-as.factor(data.alluvial2.tf$Time_group)
        data.alluvial2.tf$variable<-as.numeric(data.alluvial2.tf$variable)*50
        data.alluvial2.tf$value<-as.factor(data.alluvial2.tf$value)
        data.alluvial2.tf$Gene<-as.factor(data.alluvial2.tf$Gene)
        #
        if(is.null(colnames(table.DE.time))){
            label.test.allu<-paste0("t", seq_len(ncol(table.DE.time)))
        }else{
            label.test.allu<-colnames(table.DE.time)
        }# if(is.null(colnames(table.DE.time)))
        #-------------------------------------------------------------------------#
        # Data for alluvial graph : number DE gene per temporal group at each time
        data.alluvial2.prepare<-stats::aggregate(data.alluvial[,seq_len(ncol(table.DE.time))],
                                                 by=list(Category=data.alluvial$Time_group),
                                                 FUN=sum)
        # 1:ncol(table.DE.time)
        if(length(unique(data.alluvial2.t$Time_group))<=ncol(table.DE.time)){
            add.TG.2<-cbind(seq_len(ncol(table.DE.time)),
                            diag(rep(0,times=ncol(table.DE.time))))
            colnames(add.TG.2)<-colnames(data.alluvial2.prepare)
            #
            data.alluvial2.prepare.f<-rbind(data.alluvial2.prepare,
                                            add.TG.2[-sort(unique(data.alluvial2.t.ini$Time_group)),])
        }else{
            data.alluvial2.prepare.f<-data.alluvial2.prepare
        }# if(length(unique(data.alluvial2.t$Time_group))<=ncol(table.DE.time))
        #
        data.alluvial2<-reshape2::melt(data=data.alluvial2.prepare.f,
                                       id.vars=c("Category"))

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 3) Graphs
        ## Graph alluvial 1
        Guide.Legend.col<-ggplot2::guide_legend(override.aes=list(colour=c("black"),
                                                                  fill=c('black',
                                                                         'grey')))

        levels(data.alluvial2.tf$value)<-c("no", "yes")[as.numeric(levels(data.alluvial2.tf$value))+1]

        q.alluvial<-ggplot2::ggplot(data.alluvial2.tf,
                                    ggplot2::aes(x=variable,
                                                 y=Freq,
                                                 alluvium=Gene,
                                                 fill=Time_group,
                                                 stratum=value)) +
            ggalluvial::stat_alluvium(geom="flow", lode.guidance="forward",
                                      width=5, reverse=FALSE) +
            ggalluvial::stat_stratum(width=5,
                                     mapping=ggplot2::aes(size=value),
                                     fill=Fill.col.bar, reverse=FALSE,
                                     alpha=0.7) +
            ggplot2::labs(x="Time", y="Number of genes",
                          fill="Temporal group") +
            ggplot2::scale_x_continuous(breaks=seq_len(ncol(table.DE.time))*50,
                                        labels=label.test.allu)+
            ggplot2::scale_size_manual("DE versus t0", values=c(0.5 ,0.5),
                                       guide=Guide.Legend.col)
        #
        if(is.null(title.alluvial)){
            q.alluvial<-q.alluvial+ ggplot2::ggtitle("Alluvial graph")
        }else{
            q.alluvial<-q.alluvial+ ggplot2::ggtitle(title.alluvial)
        }# if(is.null(title.alluvial))

        ##--------------------------------------------------------------------#
        ## Graph alluvial 2
        q.freq.per.time.per.group<-ggplot2::ggplot(data=data.alluvial2,
                                                   ggplot2::aes(x=variable,
                                                                y=value,
                                                                group=as.character(Category)))+
            ggplot2::guides(color=ggplot2::guide_legend(title="Temporal group",
                                                        title.position="left"),
                            fill=ggplot2::guide_legend(title="Temporal group",
                                                       title.position="left"))+
            ggplot2::geom_line(ggplot2::aes(color=as.character(Category)),
                               linewidth=1) +
            ggplot2::geom_point(ggplot2::aes(color=as.character(Category)),
                                size=2) +
            ggplot2::geom_area(aes(fill=as.character(Category),
                                   group=as.character(Category)),
                               alpha=0.2, position='identity')+
            ggplot2::theme(axis.text.x=ggplot2::element_text(angle=0,
                                                             hjust=0.5),
                           legend.position="bottom") +
            ggplot2::xlab("Time") + ggplot2::ylab("Number of genes")
        #
        if(is.null(title.evolution)){
            q.freq.per.time.per.group<-q.freq.per.time.per.group+
                ggplot2::ggtitle("Time evolution of the number of DE genes within each temporal group.")
        }else{
            q.freq.per.time.per.group<-q.freq.per.time.per.group+
                ggplot2::ggtitle(title.evolution)
        }# if(is.null(title.evolution))

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 4) Output
        return(list(g.alluvial=q.alluvial,
                    g.alluvial.freq=q.freq.per.time.per.group))
    }else{
        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 1) Parameters : column names
        if(is.null(colnames(table.DE.time))){
            Dat.col.name<-paste0("t", seq_len(ncol(table.DE.time)))
        }else{
            Dat.col.name<-colnames(table.DE.time)
        }# if(is.null(colnames(table.DE.time)))

        ## To avoid "no visible binding for global variable" with
        ## devtools::check()
        variable<-Freq<-Gene<-Attribute<-value<-Category<-NULL

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 2) Data for graphs
        ## Only DE genes at least at one times versus t0 are kept
        gene.DE.number<-which(apply(table.DE.time, MARGIN=1,
                                    function(x) sum(x))!=0)

        if(length(gene.DE.number)==1){
            table.DE.time.DE.gene<-t(as.matrix(table.DE.time[gene.DE.number,]))
        }else{
            table.DE.time.DE.gene<-table.DE.time[gene.DE.number,]
        }
        colnames(table.DE.time.DE.gene)<-Dat.col.name

        ##--------------------------------------------------------------------#
        ## Position of DE genes in table.group.gene.DE.peak
        ## Temporal group is calculated and corresponds to the first time
        ## a gene is DE.
        DE.nb.G<-apply(as.data.frame(table.DE.time.DE.gene),
                       MARGIN=1,
                       FUN=sum)
        DE.nb.G[-which(DE.nb.G%in%c(1,ncol(table.DE.time.DE.gene)))]<-ncol(table.DE.time.DE.gene)+10
        DE.nb.G<-as.factor(as.numeric(DE.nb.G))
        IdAttribute<-which(c(1,ncol(table.DE.time.DE.gene),
                             ncol(table.DE.time.DE.gene)+10)%in%levels(DE.nb.G))

        levels(DE.nb.G)<-c("Specific", "Common", "Other")[IdAttribute]

        if(is.null(row.names(table.DE.time.DE.gene))){
            name.g<-as.character(seq_len(length(gene.DE.number)))
        }else{
            name.g<-row.names(table.DE.time.DE.gene)
        }# if(is.null(row.names(table.DE.time.DE.gene)))
        #
        DE.genes.info<-data.frame(Index=gene.DE.number,
                                  Name=name.g,
                                  Attribute=DE.nb.G)

        ##--------------------------------------------------------------------#
        ## Data for alluvial graph : Patterns (1 for DE vs t0),
        ## time group & Freq
        data.alluvial<-cbind(table.DE.time.DE.gene,
                             data.frame(Attribute=DE.genes.info$Attribute,
                                        Freq=rep(1,times=nrow(table.DE.time.DE.gene))))

        data.alluvial2.t.ini<-stats::aggregate(x=list(Freq=rep(1,nrow(data.alluvial))),
                                               by=data.alluvial[,seq_len(ncol(table.DE.time)+1)],
                                               FUN=length)

        Id.T.0DE<-which(apply(data.alluvial2.t.ini[,seq_len(ncol(table.DE.time))],
                              2, FUN=function(x) 0%in%x))
        Id.T.1DE<-which(apply(data.alluvial2.t.ini[,seq_len(ncol(table.DE.time))],
                              2, FUN=function(x) 1%in%x))
        Fill.col.bar.ini<-rep(c('black','grey'),times=ncol(table.DE.time))


        if(length(Id.T.1DE)==0 | length(Id.T.0DE)==0){
            if(length(Id.T.1DE)==0){
                Keep.col<-2*(c(as.numeric(Id.T.0DE))-1)+1
            }else{
                Keep.col<-2*c(as.numeric(Id.T.1DE))
            }# if(length(Id.T.1DE)==0)
        }else{
            Keep.col<-sort(c(2*(c(as.numeric(Id.T.0DE))-1)+1,
                             2*c(as.numeric(Id.T.1DE))))
        }# if(length(Id.T.1DE)==0 | length(Id.T.0DE)==0)
        #
        Fill.col.bar<-Fill.col.bar.ini[Keep.col]
        #
        data.alluvial2.t<-data.alluvial2.t.ini
        data.alluvial2.tf<-cbind(reshape2::melt(data=cbind(Gene=seq_len(nrow(data.alluvial2.t)),
                                                           data.alluvial2.t),
                                                id.vars=c("Gene", "Attribute",
                                                          "Freq")))
        data.alluvial2.tf[!duplicated(data.alluvial2.tf), ]
        # data.alluvial2.t$variable=as.numeric(data.alluvial2.t$variable)
        data.alluvial2.tf$Attribute<-as.factor(data.alluvial2.tf$Attribute)
        data.alluvial2.tf$variable<-as.numeric(data.alluvial2.tf$variable)*50
        data.alluvial2.tf$value<-as.factor(data.alluvial2.tf$value)
        # levels(data.alluvial2.tf$value)<-c("grey","black")
        data.alluvial2.tf$Gene<-as.factor(data.alluvial2.tf$Gene)
        #
        if(is.null(colnames(table.DE.time))){
            label.test.allu<-paste0("t", seq_len(ncol(table.DE.time)))
        }else{# 1:ncol(table.DE.time) c('t1','t2','t3')
            label.test.allu<-colnames(table.DE.time)
        }# if(is.null(colnames(table.DE.time)))

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## 3) Graphs
        # Graph alluvial 1
        Guide.Legend.col<-ggplot2::guide_legend(override.aes=list(colour=c("black"),
                                                                  fill=c('black','grey')))
        #
        ColorAllu<-data.alluvial2.tf$Attribute
        levels(ColorAllu)<-c("#0099B4FF", "#42B540FF",
                             "#F39B7FB2")[IdAttribute]
        ## 4DBBD5B2, #00A087B2, #42B54099
        alphaAllu<-data.alluvial2.tf$Attribute
        levels(alphaAllu)<-c(0.7, 1, 0.6)[IdAttribute]
        #
        levels(data.alluvial2.tf$value)<-c("no", "yes")[as.numeric(levels(data.alluvial2.tf$value))+1]
        #
        g.alluvial<-ggplot2::ggplot(data.alluvial2.tf,
                                    ggplot2::aes(x=variable, y=Freq,
                                                 alluvium=Gene,
                                                 fill=Attribute,
                                                 stratum=value)) +
            ggalluvial::stat_alluvium(geom="flow", lode.guidance="forward",
                                      width=5, reverse=FALSE,#color="grey",
                                      alpha=as.numeric(as.character(alphaAllu))) +
            ggplot2::scale_fill_manual(values=levels(ColorAllu))+
            ggalluvial::stat_stratum(width=5, mapping=ggplot2::aes(size=value),
                                     fill=Fill.col.bar,
                                     reverse=FALSE, alpha=0.7) +
            ggplot2::labs(x="", y="Number of genes")+
            ggplot2::scale_x_continuous(breaks=seq_len(ncol(table.DE.time))*50,
                                        labels=label.test.allu,
                                        guide=guide_axis(angle=45))+
            ggplot2::scale_size_manual("DE versus t0", values=c(0.5, 0.5),
                                       guide=Guide.Legend.col)
        #
        if(is.null(title.alluvial)){
            g.alluvial<-g.alluvial+ ggplot2::ggtitle("Alluvial graph")
        }else{
            g.alluvial<-g.alluvial+ ggplot2::ggtitle(title.alluvial)
        }# if(is.null(title.alluvial))

        ##--------------------------------------------------------------------#
        ##--------------------------------------------------------------------#
        ## Output
        return(g.alluvial)
    }## if(isTRUE(Temporal.Group))
}## DEplotAlluvial()
