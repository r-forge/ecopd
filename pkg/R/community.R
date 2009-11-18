# Generate site-by-species matrix
siteBySpecies <- function(phy, presence.only=FALSE,
    na.zero=FALSE, transpose=FALSE) {

    if (presence.only) {
        dat <- presence(phy, na.zero=na.zero)
    } else {
        dat <- abundance(phy, na.zero=na.zero)
    }
    mat <- as.matrix(dat)
    if (!transpose) mat <- t(mat)
    return(mat)

}
