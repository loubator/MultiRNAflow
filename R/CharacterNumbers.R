#' @title Transformation of a vector of integers into
#' a vector of class "character".
#'
#' @description Transformation of a vector of integers into
#' a vector of class "character" so that lexicographic order of characters
#' corresponds to the numerical order of time measurements.
#'
#' @param Vect.number Vector of integers.
#'
#' @details An appropriate number of character "0" is added in front of
#' the numerical characters corresponding to the decimal writing of
#' each integer in \code{Vect.number} so that the order of elements of
#' the vector is preserved. For example, "9">"11", but "09"<"11".
#'
#' @return A vector where each integer is transformed in class "character".
#'
#' @seealso The function is called by
#' [ColnamesToFactors()].
#'
#' @export
#'
#' @examples
#' CharacterNumbers(Vect.number=c(0,1,9,11,90,99,100,101))
#' CharacterNumbers(Vect.number=0:11)
#' CharacterNumbers(Vect.number=1:8)


CharacterNumbers <- function(Vect.number) {
    ##-----------------------------------------------------------------------##
    ## Check
    if (!is.numeric(Vect.number)) {
        stop("'Vect.number' must be a numeric vector")
    }## if(!is.numeric(Vect.number))

    ##-----------------------------------------------------------------------##
    ## Setting
    vNBno0 <- Vect.number
    vNBno0[which(vNBno0 == 0)] <- 1.5
    ## because of log10 which is used to computed the maximum number of digits
    Max.digit <- floor(log10(abs(max(vNBno0)))) + 1

    ##-----------------------------------------------------------------------##
    ## Vector which will containsthe same numbers reshaped as character
    Ch.numb <- rep(NA, times=length(vNBno0))
    if (Max.digit == 1) {
        ## Numbers belong to 0 and 9 so no need to add '0' in front of numbers
        Ch.numb <- as.character(vNBno0)
    } else {
        for (dg in seq_len(Max.digit)) {##1:Max.digit
            I.dig <- which(floor(log10(abs(vNBno0))) + 1 == dg)
            if (dg < Max.digit) {
                Ch.numb[I.dig] <- paste0(paste0(rep(0, times=Max.digit-dg),
                                                collapse=""),
                                         vNBno0[I.dig])
            } else {
                Ch.numb[I.dig] <- as.character(vNBno0[I.dig])
            }## if(dg<Max.digit)
        }## for (dg in seq_len(Max.digit))
    }## if(Max.digit==1) ## paste(..., sep = "") is equivalent to paste0(...)

    ##-----------------------------------------------------------------------##
    ## return
    return(Ch.Number=gsub("1.5", "0", Ch.numb, fixed=TRUE))
}## CharacterNumbers()
