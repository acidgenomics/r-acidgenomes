#' Make EnsDb object from a GFF/GTF file
#'
#' @export
#' @note Updated 2021-01-14.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `EnsDb`.
#'
#' @examples
#' file <- pasteURL(AcidGenomesTestsURL, "ensembl.gtf")
#' edb <- makeEnsDbFromGFF(
#'     file = file,
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     release = 100L
#' )
makeEnsDbFromGFF <- function(
    file,
    organism = NULL,
    genomeBuild = NULL,
    release = NULL
) {
    assert(
        isString(file),
        isString(organism, nullOK = TRUE),
        isString(genomeBuild, nullOK = TRUE),
        isInt(release, nullOK = TRUE)
    )
    requireNamespaces("ensembldb")
    alert(sprintf("Making {.var %s} from {.file %s}.", "EnsDb", file))
    tmpfile <- .cacheIt(file)
    ext <- fileExt(file)
    args <- list("outfile" = tempfile())
    ## e.g. 11f4650926be_Homo_sapiens.GRCh38.102.gtf.gz
    pattern <- paste0(
        "^([a-z0-9]+_)?",        # temp prefix from BiocFileCache.
        "([A-Z][a-z]+_[a-z]+)",  # organism (e.g. "Homo_sapiens").
        "\\.([A-Za-z0-9]+)",     # genomeVersion (e.g. "GRCh38").
        "\\.([0-9]+)",           # (Ensembl release) version (e.g. "102").
        "\\.g[ft]f",
        "(\\.gz)?$"
    )
    if (isTRUE(grepl(pattern = pattern, x = basename(file)))) {
        alert("Detecting genome metadata from file name.")
        match <- str_match(string = basename(file), pattern = pattern)
        organism <- gsub(pattern = "_", replacement = " ", x = match[1L, 3L])
        genomeBuild <- match[1L, 4L]
        release <- as.integer(match[1L, 5L])
    } else if (any(is.null(organism), is.null(genomeBuild), is.null(release))) {
        ## nocov start
        stop(sprintf(
            paste(
                "Failed to match genome metadata from file name ('%s').",
                "Define values manually: %s."
            ),
            file,
            toString(c("organism", "genomeBuild", "release"))
        ))
        ## nocov end
    }
    args <- append(
        x = args,
        values = list(
            "organism" = organism,
            "genomeVersion" = genomeBuild,
            "version" = release
        )
    )
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
    validObject(edb)
    edb
}
