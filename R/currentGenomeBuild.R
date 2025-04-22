#' Current genome build
#'
#' Fetch the current genome build (assembly) version from online resources.
#'
#' @name currentGenomeBuild
#' @note Updated 2023-10-04.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `character(1)`.
#' Genome assembly build version.
#'
#' @seealso
#' - [Ensembl REST API](https://rest.ensembl.org).
#' - [UCSC REST API](https://genome.ucsc.edu/goldenPath/help/api.html).
#' - [RefSeq genomes](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/).
#'
#' @examples
#' ## Ensembl.
#' try({
#'     currentEnsemblGenomeBuild("Homo sapiens")
#'     currentEnsemblGenomeBuild("Mus musculus")
#' })
#'
#' ## GENCODE.
#' try({
#'     currentGencodeGenomeBuild("Homo sapiens")
#'     currentGencodeGenomeBuild("Mus musculus")
#' })
#'
#' ## RefSeq.
#' try({
#'     currentRefseqGenomeBuild(
#'         organism = "Homo sapiens",
#'         taxonomicGroup = "vertebrate_mammalian"
#'     )
#'     currentRefseqGenomeBuild(
#'         organism = "Mus musculus",
#'         taxonomicGroup = "vertebrate_mammalian"
#'     )
#' })
#'
#' ## UCSC.
#' try({
#'     currentUcscGenomeBuild("Homo sapiens")
#'     currentUcscGenomeBuild("Mus musculus")
#' })
NULL


#' Return a simple (minimal) genome build version
#'
#' @details
#' For example, sanitize "GRCh38.p13" to simply "GRCh38".
#'
#' @note Updated 2021-01-07.
#' @noRd
.simpleGenomeBuild <- function(x) {
    assert(isString(x))
    x <- sub(pattern = "\\.[^\\.]+$", replacement = "", x = x)
    x
}


## Updated 2023-07-27.
#' @rdname currentGenomeBuild
#' @export
currentEnsemblGenomeBuild <-
    function(organism) {
        assert(isOrganism(organism))
        organism <- snakeCase(organism)
        json <- getJson(pasteUrl(
            "rest.ensembl.org",
            "info",
            "assembly",
            paste0(organism, "?", "content-type=application/json"),
            protocol = "https"
        ))
        assert(isSubset("assembly_name", names(json)))
        out <- json[["assembly_name"]]
        assert(isString(out))
        out
    }


## Updated 2021-01-31.
#' @rdname currentGenomeBuild
#' @export
currentGencodeGenomeBuild <-
    function(organism) {
        organism <- match.arg(
            arg = organism,
            choices = c("Homo sapiens", "Mus musculus")
        )
        currentEnsemblGenomeBuild(organism)
    }


## Alternate approach using URL only:
## https://ftp.ncbi.nlm.nih.gov/genomes/refseq/<taxonomic_group>/<organism>/
## latest_assembly_versions/

## Updated 2021-01-14.
#' @rdname currentGenomeBuild
#' @export
#'
#' @param taxonomicGroup `character(1)`.
#' *Only applies to RefSeq*.
#' FTP server taxonomic group subdirectory path (e.g. "vertebrate_mammalian").
#' Defining this manually avoids having to query the FTP server.
currentRefseqGenomeBuild <-
    function(organism, taxonomicGroup = NULL) {
        assert(
            isOrganism(organism),
            isString(taxonomicGroup, nullOk = TRUE)
        )
        baseUrl <- .getRefSeqGenomeUrl(
            organism = organism,
            taxonomicGroup = taxonomicGroup,
            quiet = TRUE
        )
        summary <- .getRefSeqAssemblySummary(
            file = pasteUrl(baseUrl, "assembly_summary.txt")
        )
        assert(isSubset("ftp_path", names(summary)))
        out <- basename(summary[["ftp_path"]])
        out
    }


## Updated 2023-04-14.
#' @rdname currentGenomeBuild
#' @export
currentUcscGenomeBuild <-
    function(organism) {
        assert(isOrganism(organism))
        json <- getJson("https://api.genome.ucsc.edu/list/ucscGenomes")
        assert(isSubset("ucscGenomes", names(json)))
        json <- json[["ucscGenomes"]]
        lst <- Map(
            name = names(json),
            x = json,
            f = function(name, x) {
                ## Other useful keys: description, sourceName.
                c(
                    "build" = name,
                    "active" = x[["active"]],
                    "orderKey" = x[["orderKey"]],
                    "scientificName" = x[["scientificName"]]
                )
            },
            USE.NAMES = FALSE
        )
        df <- rbindToDataFrame(lst)
        df <- df[df[["active"]] == 1L, , drop = FALSE]
        if (!isSubset(organism, unique(df[["scientificName"]]))) {
            abort(sprintf("Invalid organism: {.val %s}.", organism))
        }
        df <- df[df[["scientificName"]] == organism, , drop = FALSE]
        ## The latest genome build has the lowest "orderKey" value.
        df <- df[order(df[["orderKey"]]), , drop = FALSE]
        out <- df[["build"]][[1L]]
        assert(isString(out))
        out
    }
