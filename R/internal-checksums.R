#' Generate MD5 checksum from file
#'
#' @note Updated 2021-02-02.
#' @noRd
.md5 <- function(file) {
    md5(file = .cacheIt(file))
}



#' Generate SHA256 checksum from file
#'
#' @note Updated 2021-02-02.
#' @noRd
.sha256 <- function(file) {
    sha256(file = .cacheIt(file))
}
