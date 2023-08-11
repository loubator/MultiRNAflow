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
#' RawCountsSimulation(Nb.Group=3, Nb.Time=5, Nb.per.GT=7, Nb.Gene=50)
#' ## RawCountsSimulation(Nb.Group=1, Nb.Time=5, Nb.per.GT=7, Nb.Gene=50)
#' ## RawCountsSimulation(Nb.Group=3, Nb.Time=1, Nb.per.GT=7, Nb.Gene=50)

RawCountsSimulation <- function(Nb.Group,
                                Nb.Time,
                                Nb.per.GT,
                                Nb.Gene) {
    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Check
    if (max(Nb.Group,Nb.Time) == 1) {
        stop("At least two groups or two times are demanded.")
    }## if(max(Nb.Group,Nb.Time)==1)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Setting
    Tsim <- factor(rep(rep(paste0("t", seq(from=0, to=Nb.Time-1, by=1)),
                           times=Nb.Group),
                       each=Nb.per.GT))
    Gsim <- factor(rep(paste0("G", seq_len(Nb.Group)),
                       each=Nb.Time*Nb.per.GT))

    S.pg <- rep(seq_len(Nb.per.GT), times=(Nb.Time * Nb.Group))
    S.pg <- S.pg + rep(Nb.per.GT * seq(from=0, to=Nb.Group-1, by=1),
                       each=(Nb.Time * Nb.per.GT))

    Patient <- factor(paste0("Ind", S.pg))
    Data.sim <- as.data.frame(matrix(0, nrow=Nb.Gene,
                                   ncol=(Nb.per.GT * Nb.Group * Nb.Time)))
    row.names(Data.sim) <- NULL ##paste0("Gene", seq_len(Nb.Gene))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if (min(Nb.Group,Nb.Time) > 1) {
        colnames(Data.sim) <- paste0(Gsim, "_", Tsim, "_", Patient)
        Mmodel <- stats::model.matrix(~ Tsim + Gsim + Gsim:Tsim)
        Design.M <- matrix(Mmodel, ncol=length(colnames(Mmodel)))
        colnames(Design.M) <- colnames(Mmodel)
    }## if(min(Nb.Group,Nb.Time) >1 )

    if (Nb.Group == 1) {
        colnames(Data.sim) <- paste0(Tsim, "_", Patient)
        Design.M <- matrix(stats::model.matrix(~ Tsim),
                           ncol=length(colnames(stats::model.matrix(~ Tsim))))
        colnames(Design.M) <- colnames(stats::model.matrix(~ Tsim))
    }## if(Nb.Group == 1)

    if (Nb.Time == 1) {
        colnames(Data.sim) <- paste0(Gsim, "_", Patient)
        Design.M <- matrix(stats::model.matrix(~ Gsim),
                           ncol=length(colnames(stats::model.matrix(~ Gsim))))
        colnames(Design.M) <- colnames(stats::model.matrix(~ Gsim))
    }## if(Nb.Time == 1)

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Simulation
    bio.noise <- stats::runif(n=Nb.Gene, min=0.5, max=3)
    Beta0 <- stats::rnorm(n=Nb.Gene, mean=5, sd=4)
    for (g in seq_len(Nb.Gene)) {
        Beta.g <- stats::rnorm(n=ncol(Design.M) - 1, mean=0,
                               sd=stats::rexp(n=ncol(Design.M) - 1, rate=0.3))
        mu.s.tg <- 2^min(abs(as.numeric(Design.M%*%c(Beta0[g], Beta.g))),
                         18)
        Data.sim[g,] <- stats::rnbinom(n=ncol(Data.sim),
                                       size=bio.noise[g],
                                       mu=mu.s.tg) + 1
    }## for(g in 1:Nb.Gene)

    Data.sim.f <- cbind(Gene=paste0("Gene", seq_len(Nb.Gene)), Data.sim)
    row.names(Data.sim.f) <- paste0("Gene", seq_len(Nb.Gene))

    ##------------------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    ## Output
    return(list(Sim.dat=Data.sim.f,
                Vect.Group=Gsim,
                Vect.Time=Tsim,
                Vect.Sample=Patient))
}## RawCountsSimulation()
