#' Make EnsDb object from a GFF/GTF file
#'
#' Wrapper for ensembldb importer functions.
#'
#' @name makeEnsDbFromGFF
#' @note Updated 2021-01-14.
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
#' edb <- makeEnsDbFromGFF(
#'     file = file,
#'     organism = "Homo sapiens",
#'     genomeBuild = "GRCh38",
#'     release = 100L
#' )
NULL



#' @describeIn makeEnsDbFromGFF Primary function.
#' @export
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
    ## Attempt to detect the genome metadata from the file name, if possible.
    if (
        is.null(organism) &&
        is.null(genomeBuild) &&
        is.null(release)
    ) {
        ## Ensembl file name metadata pattern.
        ## e.g. 11f4650926be_Homo_sapiens.GRCh38.102.gtf.gz
        ensemblPattern <- paste0(
            "^([a-z0-9]+_)?",        # temp prefix from BiocFileCache.
            "([A-Z][a-z]+_[a-z]+)",  # organism (e.g. "Homo_sapiens").
            "\\.([A-Za-z0-9]+)",     # genomeVersion (e.g. "GRCh38").
            "\\.([0-9]+)",           # (Ensembl release) version (e.g. "102").
            "\\.g[ft]f",
            "(\\.gz)?$"
        )
        ## GENCODE file name metadata pattern.
        ## - Human: gencode.v32.annotation.gtf.gz
        ## - Mouse: gencode.vM25.annotation.gtf.gz
        gencodePattern <- paste0(
            "gencode",
            "\\.v([M0-9]+)",
            "(lift37)?",
            "\\.annotation",
            "\\.g[ft]f",
            "(\\.gz)?$"
        )
        if (isTRUE(grepl(
            pattern = ensemblPattern,
            x = basename(file)
        ))) {
            alert("Detecting Ensembl genome metadata from file name.")
            match <- str_match(
                string = basename(file),
                pattern = ensemblPattern
            )
            organism <- gsub(
                pattern = "_",
                replacement = " ",
                x = match[1L, 3L]
            )
            genomeBuild <- match[1L, 4L]
            release <- as.integer(match[1L, 5L])
        } else if (isTRUE(grepl(
            pattern = gencodePattern,
            x = basename(file)
        ))) {
            alert("Detecting GENCODE genome metadata from file name.")
            match <- str_match(
                string = basename(file),
                pattern = gencodePattern
            )
            release <- match[1L, 2L]
            if (grepl("^M", release)) {
                organism <- "Mus musculus"
            } else {
                organism <- "Homo sapiens"
                release <- as.integer(release)
                ## GRCh38
                if (match[1L, 3L] == "lift37") {
                    genomeBuild <- "GRCh37"
                }
            }
            ## Assembly (genome build) is documented in the first commented
            ## lines of the file (line 1 for GTF; line 2 for GFF3).
            if (is.null(genomeBuild)) {
                x <- import(file, format = "lines", nMax = 2L)
                x <- grep(pattern = "description:", x = x, value = TRUE)
                match <- str_match(
                    string = x,
                    pattern = "\\sgenome\\s\\(([^\\)]+)\\),"
                )
                genomeBuild <- match[1L, 2L]
            }
        }
    }
    if (any(is.null(organism), is.null(genomeBuild), is.null(release))) {
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



#' @describeIn makeEnsDbFromGFF Alias for GTF files.
#' @export
makeEnsDbFromGTF <- makeEnsDbFromGFF
