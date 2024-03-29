#' Make EnsDb object from a GFF/GTF file
#'
#' Wrapper for ensembldb importer functions.
#'
#' @export
#' @note Updated 2023-11-29.
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
#' ## > file <- AcidBase::pasteUrl(
#' ## >     "ftp.ensembl.org",
#' ## >     "pub",
#' ## >     "release-102",
#' ## >     "gtf",
#' ## >     "homo_sapiens",
#' ## >     "Homo_sapiens.GRCh38.102.gtf.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > edb <- makeEnsDbFromGff(file)
#' ## > print(edb)
makeEnsDbFromGff <- function(file) {
    assert(
        requireNamespaces("ensembldb"),
        .isSupportedGff(file)
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
    x <- strMatch(
        x = basename(file),
        pattern = pattern
    )[1L, , drop = TRUE]
    args[["genomeBuild"]] <- x[[4L]]
    args[["organism"]] <- gsub(pattern = "_", replacement = " ", x = x[[3L]])
    args[["version"]] <- as.integer(x[[5L]])
    args[["outfile"]] <- tempfile()
    ext <- fileExt(file)
    tmpfile <- .cacheIt(file)
    meta <- .getGffMetadata(file)
    if (grepl(pattern = "gff", x = ext, ignore.case = TRUE)) {
        what <- ensembldb::ensDbFromGff
        args <- append(x = args, values = list("gff" = tmpfile))
    } else if (grepl(pattern = "gtf", x = ext, ignore.case = TRUE)) {
        what <- ensembldb::ensDbFromGtf
        args <- append(x = args, values = list("gtf" = tmpfile))
    } else {
        abort(sprintf(
            "Unsupported file: {.file %s}.", basename(file)
        ))
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
