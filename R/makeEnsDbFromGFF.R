#' Make EnsDb object from a GFF/GTF file
#'
#' Wrapper for ensembldb importer functions.
#'
#' @name makeEnsDbFromGFF
#' @note Updated 2021-01-23.
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
#' file <- pasteURL(AcidGenomesTestsURL, "ensembl.gtf")
#' edb <- makeEnsDbFromGFF(file)
makeEnsDbFromGFF <- function(file) {
    assert(isString(file))
    requireNamespaces("ensembldb")
    alert(sprintf("Making {.var %s} from {.file %s}.", "EnsDb", file))
    args <- list()
    meta <- getGFFMetadata(file)
    args[["outfile"]] <- tempfile()
    args[["genomeBuild"]] <- meta[["genomeBuild"]]
    args[["organism"]] <- meta[["organism"]]
    args[["version"]] <- meta[["release"]]
    pattern <- .gffPatterns[["ensembl"]]
    if (isTRUE(grepl(pattern = pattern, x = basename(file)))) {
        x <- str_match(
            string = basename(file),
            pattern = pattern
        )[1L, , drop = TRUE]
        if (!isString(args[["organism"]])) {
            args[["organism"]] <-
                gsub(pattern = "_", replacement = " ", x = x[[3L]])
        }
        if (!isString(args[["genomeBuild"]])) {
            args[["genomeBuild"]] <- x[[4L]]
        }
        if (!isInt(args[["version"]])) {
            args[["version"]] <- as.integer(x[[5L]])
        }
    }
    ext <- fileExt(file)
    tmpfile <- .cacheIt(file)
    if (grepl(pattern = "gff", x = ext, ignore.case = TRUE)) {
        what <- ensembldb::ensDbFromGff
        args <- append(x = args, values = list("gff" = tmpfile))
    } else if (grepl(pattern = "gtf", x = ext, ignore.case = TRUE)) {
        what <- ensembldb::ensDbFromGtf
        args <- append(x = args, values = list("gtf" = tmpfile))
    } else {
        stop("Unsupported file extension.")  # nocov
    }
    suppressWarnings({
        suppressMessages({
            sqlite <- do.call(what = what, args = args)
        })
    })
    edb <- ensembldb::EnsDb(x = sqlite)
    assert(is(edb, "EnsDb"))
    attr(x = edb, which = "args") <- args
    attr(x = edb, which = "gffMetadata") <- meta
    validObject(edb)
    edb
}
