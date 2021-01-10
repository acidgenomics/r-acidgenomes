#' Make EnsDb object from a GFF/GTF file
#'
#' @note Updated 2021-01-09.
#' @noRd
.makeEnsDbFromGFF <- function(file) {
    alert(sprintf("Making {.var %s} from {.file %s}.", "EnsDb", file))
    requireNamespaces("ensembldb")
    assert(isAFile(file))
    ext <- fileExt(file)
    args <- list(
        ## Alternatively, can consider caching using BiocFileCache.
        outfile = tempfile()
    )
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
        args[["organism"]] <-
            gsub(pattern = "_", replacement = " ", x = match[1L, 3L])
        args[["genomeVersion"]] <- match[1L, 4L]
        args[["version"]] <- as.integer(match[1L, 5L])
    } else {
        alertWarning("Failed to detect genome metadata from file name.")
    }
    if (grepl(pattern = "gff", x = ext, ignore.case = TRUE)) {
        what <- ensembldb::ensDbFromGff
        args <- c("gff" = file, args)
    } else if (grepl(pattern = "gtf", x = ext, ignore.case = TRUE)) {
        what <- ensembldb::ensDbFromGtf
        args <- c("gtf" = file, args)
    } else {
        stop("Unsupported file extension.")  # nocov
    }
    suppressWarnings({
        sqlite <- do.call(what = what, args = args)
    })
    edb <- ensembldb::EnsDb(x = sqlite)
    assert(is(edb, "EnsDb"))
    validObject(edb)
    edb
}
