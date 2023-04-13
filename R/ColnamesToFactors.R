#' @title Extraction of factors information and suitable column names creation
#' from the column names of a dataset.
#'
#' @description
#' This function generates new reduced column names according to the presence
#' of biological conditions and/or time points, and extract the different
#' factors (individual's names, time measurements, biological conditions)
#' from the column names of the dataset (see \code{Details}).
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
#' 'r1' is localized after the third underscore,
#' so \code{Individual.position=4},
#' 'P' is localized after the first underscore, so \code{Group.position=2} and
#' 't0' is localized after the second underscore, so \code{Time.position=3}.
#'
#' @seealso The [ColnamesToFactors()] function is used by the following
#' functions of our package : [DEanalysisPreprocessing()],
#' [PCApreprocessing()], [MFUZZclustersNumber()] and [MFUZZanalysis()].
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
#'
#' @return The function returns new column names of the dataset,
#' a vector indicating the name of the individual for each sample,
#' a vector indicating the time for each sample and/or
#' a vector indicating the biological condition for each sample.
#'
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#' ## Data simulated with our function RawCountsSimulation()
#' Data.sim<-RawCountsSimulation(Nb.Group=3, Nb.Time=2, Nb.per.GT=3,
#'                               Nb.Gene=10)
#' ##-------------------------------------------------------------------------#
#' res.test.colnames=ColnamesToFactors(ExprData=Data.sim$Sim.dat,
#'                                     Column.gene=1,
#'                                     Group.position=1,
#'                                     Time.position=2,
#'                                     Individual.position=3)
#' print(res.test.colnames)

ColnamesToFactors<-function(ExprData,
                            Column.gene,
                            Group.position,
                            Time.position,
                            Individual.position){
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Every sample must have an indidual name
    if(is.null(Individual.position)){
        stop("Every sample must have an indidual name (name or number).")
    }## if(is.null(Individual.position))

    ## Biological condition & Times point absent
    if(is.null(Group.position) & is.null(Time.position)){
        stop("Samples must belong to at least one time or one group.")
    }## if(is.null(Group.position) & is.null(Time.position))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Column names with underscore
    colnames.with.underscore<-colnames(ExprData)
    ## Index of each sample
    if(is.null(Column.gene)){
        ind.col.expr<-seq_len(length(colnames.with.underscore))
    }else{
        ind.col.expr<-seq_len(length(colnames.with.underscore))[-Column.gene]
    }## if(is.null(Column.gene))
    ##
    Vect.colnames<-colnames.with.underscore[ind.col.expr]

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Transform colnames into a matrix containing in each row
    ## Individual information, Group information and or Time information
    Colnames.matrix.info<-matrix(unlist(strsplit(Vect.colnames,
                                                 split="_",
                                                 fixed=TRUE)),
                                 ncol=length(Vect.colnames))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Case when individual names are numbers
    Individual.Names<-Colnames.matrix.info[Individual.position,]
    if(!is.numeric(Individual.Names)){
        Individual.info<-Individual.Names
    }else{
        Individual.info<-paste("r", CharacterNumbers(Individual.Names), sep="")
    }## if(!is.numeric(Individual.Names))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Biological conditions & Times points present
    if(!is.null(Group.position) & !is.null(Time.position)){
        ##--------------------------------------------------------------------#
        ## Setting
        Tps.info.ini<-as.character(Colnames.matrix.info[Time.position,])
        Tps.info<-gsub("T", "", gsub("t", "", Tps.info.ini))
        Time.info.f<-paste0("t", CharacterNumbers(as.numeric(Tps.info)))
        Group.info<-as.character(Colnames.matrix.info[Group.position,])

        ##--------------------------------------------------------------------#
        ## Cases when algorithm must stop
        Contingency.GT<-matrix(table(Group.info, Time.info.f),
                               ncol=ncol(table(Group.info, Time.info.f)),
                               dimnames=dimnames(table(Group.info,
                                                       Time.info.f)))
        Var.GT<-sum(apply(Contingency.GT, 1, stats::var))
        ##
        Contingency.IT<-matrix(table(Individual.info,Time.info.f),
                               ncol=ncol(table(Individual.info, Time.info.f)),
                               dimnames=dimnames(table(Individual.info,
                                                       Time.info.f)))
        Var.IT<-sum(apply(Contingency.IT, 1, stats::var))

        ##--------------------------------------------------------------------#
        ## Check, stop
        if(Var.GT + Var.IT + max(Contingency.IT) - 1 > 0){
            Stop.BC.T<-paste("Every individual must have a unique name,",
                             "must be associated to a unique group and must",
                             "be associated only once to each of the same",
                             ncol(Contingency.IT), "time measurements.")
            stop(Stop.BC.T)
        }## if(Var.GT + Var.IT + max(Contingency.IT)-1>0)
        ##
        if(min(c(table(Group.info,Time.info.f))) < 2){
            stop("Each group must have at least two individuals.")
        }## if(min(c(table(Group.info,Time.info.f)))<2)

        ##--------------------------------------------------------------------#
        ## Final name
        Final.Name<-paste0(Group.info, ".", Individual.info, ".", Time.info.f)
    }## if(is.null(Group.position)==FALSE & is.null(Time.position)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Biological condition present & Times points absent
    if(is.null(Group.position)==FALSE & is.null(Time.position)==TRUE){
        ##--------------------------------------------------------------------#
        ## Setting
        Time.info.f<-Tps.info<-NULL
        Group.info<-as.character(Colnames.matrix.info[Group.position,])

        ##--------------------------------------------------------------------#
        ## Cases when algorithm must stop
        Contingency.IG<-matrix(table(Individual.info,Group.info),
                               ncol=ncol(table(Individual.info,Group.info)),
                               dimnames=dimnames(table(Individual.info,
                                                       Group.info)))
        ##
        Var.IG<-stats::var(apply(Contingency.IG, 1, sum))

        ##--------------------------------------------------------------------#
        ## Check, stop
        if(Var.IG + max(apply(Contingency.IG,1,sum)) - 1 > 0){
            stop("Every individual must be associated to only one group.")
        }## if(Var.IG + max(apply(Contingency.IG,1,sum))-1 >0)
        ##
        if(min(apply(Contingency.IG,2,sum)) < 2){
            stop("Each group must have at least two individuals.")
        }## if(min(apply(Contingency.IG,2,sum))<2)

        ##--------------------------------------------------------------------#
        ## Final name
        Final.Name<-paste(Group.info, ".", Individual.info, sep="")
    }## if(is.null(Group.position)==FALSE & is.null(Time.position)==TRUE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Biological condition absent & Times points present
    if(is.null(Group.position) & !is.null(Time.position)){
        ##--------------------------------------------------------------------#
        ## Setting
        Tps.info.ini<-as.character(Colnames.matrix.info[Time.position,])
        Tps.info<-gsub("T", "", gsub("t", "", Tps.info.ini))
        Time.info.f<-paste0("t", CharacterNumbers(as.numeric(Tps.info)))
        Group.info<-NULL

        ##--------------------------------------------------------------------#
        ## Cases when algorithm must stop
        if(length(Time.info.f) == length(unique(Time.info.f))){
            Stop.Tinfo<-paste("The data must contain the temporal expression",
                              "of at least two individuals.", sep=" ")
            stop(Stop.Tinfo)
        }## if(length(Time.info.f)==length(unique(Time.info.f)))
        ##
        Contingency.IT<-matrix(table(Individual.info, Time.info.f),
                               ncol=ncol(table(Individual.info, Time.info.f)),
                               dimnames=dimnames(table(Individual.info,
                                                       Time.info.f)))
        Var.IT<-sum(apply(Contingency.IT, 1, stats::var))

        ##--------------------------------------------------------------------#
        ## Check, stop
        if(Var.IT + max(Contingency.IT) - 1 > 0){
            Stop.T<-paste("Every individual must have a unique name and",
                          "must be associated only once to each of the same",
                          ncol(Contingency.IT), "time measurements.")
            stop(Stop.T)
        }## if(Var.IT + max(Contingency.IT)-1>0)

        ##--------------------------------------------------------------------#
        Final.Name<-paste0(Individual.info, ".", Time.info.f)
    }## if(is.null(Group.position)==TRUE & is.null(Time.position)==FALSE)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Data with initial and final sample names
    data.code.names<-data.frame(Initial.name=Vect.colnames,
                                Final.Name=Final.Name)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(list(Final.Name=Final.Name,
                Group.Info=Group.info,
                Time.Info=Time.info.f,
                Individual.info=Individual.info,
                Data.code.names=data.code.names))
}## ColnamesToFactors()
