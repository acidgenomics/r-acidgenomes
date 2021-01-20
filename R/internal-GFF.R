#' Get genome metadata from a GFF file
#'
#' @note Updated 2021-01-18.
#' @noRd
#'
#' @examples
#' url <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' genomeBuild <- .detectGenomeBuildFromGFF(url)
#' print(genomeBuild)

## Ensembl genomes

.detectGenomeBuildFromGFF <- function(file) {
    df <- getGFFMetadata(file, nMax = 500L)
    if (!is(df, "DataFrame")) return(NULL)
    ## GENCODE files have a description key that contains the genome build.
    if (isTRUE("description" %in% df[["key"]])) {
        string <- df[df[["key"]] == "description", "value", drop = TRUE]
        x <- str_match(
            string = string,
            pattern = "genome \\(([^\\)]+)\\)"
        )[1L, 2L]
        if (isString(x)) return(x)
    }
    ## Otherwise we can parse for standard "genome-build" key, which is
    ## supported by Ensembl and RefSeq.
    .getValue <- function(key) {
        x <- df[match(x = key, table = df[["key"]]), "value", drop = TRUE]
        if (is.na(x)) return(NULL)
        x
    }
    x <- .getValue("genome-build")
    if (isString(x)) return(x)
    NULL
}
