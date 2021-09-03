#' Download multiple genome files in a single call
#'
#' @note Updated 2021-08-03.
#' @noRd
#'
#' @return `character`
#'   Local file paths.
.downloadURLs <- function(
    urls,
    outputDir,
    cache
) {
    assert(
        allAreURLs(urls),
        isString(outputDir),
        isFlag(cache)
    )
    outputDir <- initDir(outputDir)
    destFiles <- vapply(
        X = urls,
        FUN = function(url) {
            file.path(outputDir, basename(url))
        },
        FUN.VALUE = character(1L),
        USE.NAMES = TRUE
    )
    assert(identical(names(urls), names(destFiles)))
    if (isTRUE(cache)) {
        cacheFiles <- .cacheIt(file = urls)
        file.copy(
            from = cacheFiles,
            to = destFiles,
            overwrite = TRUE,
            recursive = FALSE
        )
    } else {
        mapply(
            url = urls,
            destfile = destFiles,
            FUN = download,
            SIMPLIFY = TRUE,
            USE.NAMES = FALSE
        )
    }
    out <- destFiles
    names(out) <- names(urls)
    assert(allAreFiles(out))
    invisible(out)
}



#' Get an internal function from the package
#'
#' @note Updated 2021-03-03.
#' @noRd
.getFun <- function(x) {
    fun <- get(
        x = x,
        envir = asNamespace(.pkgName),
        inherits = FALSE
    )
    assert(is.function(fun))
    fun
}



#' Is the input GFF file supported in the package?
#'
## See `.gffPatterns` for pattern matching details.
#'
#' @note Updated 2021-08-06.
#' @noRd
.isSupportedGFF <- function(file) {
    ok <- isString(file)
    if (!ok) { return(FALSE) }
    denylist <- c(
        "flybase_gff" = paste0(
            "^([a-z0-9]+_)?",
            "^([^-]+)",
            "-([^-]+)",
            "-(r[0-9]+\\.[0-9]+)",
            "\\.gff",
            "(\\.gz)?$"
        ),
        "wormbase_gff" = paste0(
            "^([a-z0-9]+_)?",
            "^([a-z]_[a-z]+)",
            "\\.([A-Z0-9]+)",
            "\\.(WS[0-9]+)",
            "\\.([a-z_]+)",
            "\\.gff3",
            "(\\.gz)?$"
        )
    )
    if (isMatchingRegex(
        pattern = denylist[["flybase_gff"]],
        x = basename(file)
    )) {
        alertWarning("Use FlyBase GTF instead of GFF.")
        return(FALSE)
    } else if (isMatchingRegex(
        pattern = denylist[["wormbase_gff"]],
        x = basename(file)
    )) {
        alertWarning("Use WormBase GTF instead of GFF.")
        return(FALSE)
    }
    TRUE
}



#' Map a genome build to UCSC
#'
#' @details
#' Currently used internally for handoff to `Seqinfo` (e.g. GENCODE genome).
#'
#' @note Updated 2021-08-04.
#' @noRd
#'
#' @examples
#' .mapGenomeBuildToUCSC("GRCh38")
.mapGenomeBuildToUCSC <- function(x) {
    assert(isString(x))
    switch(
        EXPR = x,
        "BDGP5"             = "dm3",       # Drosophila melanogaster
        "BDGP6"             = "dm6",       # Drosophila melanogaster
        "CanFam3.1"         = "canFam3",   # Canis familiaris
        "GRCg6a"            = "galGal6",   # Gallus gallus
        "GRCh37"            = "hg19",      # Homo sapiens
        "GRCh38"            = "hg38",      # Homo sapiens
        "GRCm37"            = "mm9",       # Mus musculus
        "GRCm38"            = "mm10",      # Mus musculus
        "GRCm39"            = "mm39",      # Mus musculus
        "GRCz10"            = "danRer10",  # Danio rerio
        "GRCz11"            = "danRer11",  # Danio rerio
        "Galgal4"           = "galGal4",   # Gallus gallus
        "Gallus_gallus-5.0" = "galGal5",   # Gallus gallus
        "JGI_4.1"           = "xenTro2",   # Xenopus tropicalis
        "JGI_4.2"           = "xenTro3",   # Xenopus tropicalis
        "R64-1-1"           = "sacCer3",   # Saccharomyces cerevisiae
        "Rnor_5.0"          = "rn5",       # Rattus norvegicus
        "Rnor_6.0"          = "rn6",       # Rattus norvegicus
        "Sscrofa11.1"       = "susScr11",  # Sus scrofa
        "TAIR10"            = "araTha1",   # Arabidopsis thaliana
        "WBcel235"          = "ce11",      # Caenorhabditis elegans
        "WS220"             = "ce10",      # Caenorhabditis elegans
        abort(sprintf("Unsupported genome build: {.val %s}.", x))
    )
}
