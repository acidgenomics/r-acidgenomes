#' Make EnsDb object from a GFF/GTF file
#'
#' Wrapper for ensembldb importer functions.
#'
#' @export
#' @note Updated 2023-04-26.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `EnsDb`.
#'
#' @seealso
#' - `ensembldb::ensDbFromGff()`.
#' - `ensembldb::ensDbFromGtf()`.
#'
#' @examples
#' ## > file <- AcidBase::pasteURL(
#' ## >     "ftp.ensembl.org",
#' ## >     "pub",
#' ## >     "release-102",
#' ## >     "gtf",
#' ## >     "homo_sapiens",
#' ## >     "Homo_sapiens.GRCh38.102.gtf.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > edb <- makeEnsDbFromGFF(file)
#' ## > print(edb)
makeEnsDbFromGFF <- function(file) {
    assert(
        requireNamespaces("ensembldb"),
        .isSupportedGFF(file)
    )
    if (isAFile(file)) {
        file <- realpath(file)
    }
    pattern <- .gffPatterns[["ensembl"]]
    assert(
        grepl(pattern = pattern, x = basename(file)),
        msg = "Failed to detect Ensembl GFF file."
    )
    alert(sprintf("Making {.cls %s} from {.file %s}.", "EnsDb", file))
    args <- list()
    x <- stri_match_first_regex(
        str = basename(file),
        pattern = pattern
    )[1L, , drop = TRUE]
    args[["genomeBuild"]] <- x[[4L]]
    args[["organism"]] <- gsub(pattern = "_", replacement = " ", x = x[[3L]])
    args[["version"]] <- as.integer(x[[5L]])
    args[["outfile"]] <- tempfile()
    ext <- fileExt(file)
    tmpfile <- .cacheIt(file)
    meta <- .getGFFMetadata(file)
    if (grepl(pattern = "gff", x = ext, ignore.case = TRUE)) {
        what <- ensembldb::ensDbFromGff
        args <- append(x = args, values = list("gff" = tmpfile))
    } else if (grepl(pattern = "gtf", x = ext, ignore.case = TRUE)) {
        what <- ensembldb::ensDbFromGtf
        args <- append(x = args, values = list("gtf" = tmpfile))
    } else {
        ## nocov start
        abort(sprintf(
            "Unsupported file: {.file %s}.", basename(file)
        ))
        ## nocov end
    }
    quietly({
        sqlite <- do.call(what = what, args = args)
        edb <- ensembldb::EnsDb(x = sqlite)
        file.remove(sqlite)
    })
    assert(is(edb, "EnsDb"))
    attr(x = edb, which = "args") <- args
    attr(x = edb, which = "gffMetadata") <- meta
    assert(validObject(edb))
    edb
}
