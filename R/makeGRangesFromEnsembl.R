## FIXME ADD SUPPORT FOR EXONS AND CDS.

## FIXME NEED TO ENSURE CALL IS SLOTTED INTO OBJECT.

## The tximeta is thinking about this approach very similarly.
##
## FIXME CONSIDER Returning useful Seqinfo, similar to tximeta.
##
## Useful functions to reference:
## - tximeta:::gtf2RefSeq
##
## See also:
## - https://github.com/mikelove/tximeta/blob/master/R/tximeta.R
## - https://github.com/mikelove/tximeta/blob/master/tests/testthat/test_tximeta.R



#' Make GRanges from Ensembl
#'
#' Quickly obtain gene and transcript annotations from
#' [Ensembl](https://www.ensembl.org/) using
#' [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/) and
#' [ensembldb](https://bioconductor.org/packages/ensembldb/).
#'
#' Simply specify the desired organism, using the full latin name. For example,
#' we can obtain human annotations with `Homo sapiens`. Optionally, specific
#' Ensembl genome builds (e.g. `GRCh38`) and release versions (e.g. `87`) are
#' supported.
#'
#' Under the hood, this function fetches annotations from AnnotationHub using
#' the ensembldb package. AnnotationHub supports versioned Ensembl releases,
#' back to version 87.
#'
#' Genome build: use `"GRCh38"` instead of `"hg38"` for the genome build, since
#' we're querying Ensembl and not UCSC.
#'
#' @section Broad class definitions:
#'
#' For gene and transcript annotations, a `broadClass` column is added, which
#' generalizes the gene types into a smaller number of semantically-meaningful
#' groups:
#'
#'   - `coding`.
#'   - `noncoding`.
#'   - `pseudo`.
#'   - `small`.
#'   - `decaying`.
#'   - `ig` (immunoglobulin).
#'   - `tcr` (T cell receptor).
#'   - `other`.
#'
#' @section GRCh37 (hg19) legacy annotations:
#'
#' [makeGRangesFromEnsembl()] supports the legacy *Homo sapiens* GRCh37 (release
#' 75) build by internally querying the [EnsDb.Hsapiens.v75][] package.
#' Alternatively, the corresponding GTF/GFF file can be loaded directly from
#' GENCODE or Ensembl.
#'
#' [EnsDb.Hsapiens.v75]: https://bioconductor.org/packages/EnsDb.Hsapiens.v75/
#'
#' @section AnnotationHub queries:
#'
#' Here's how to perform manual, customized AnnotationHub queries.
#'
#' ```
#' library(AnnotationHub)
#' library(ensembldb)
#' ah <- AnnotationHub()
#'
#' # Human ensembldb (EnsDb) records.
#' ahs <- query(
#'     x = ah,
#'     pattern = c(
#'         "Homo sapiens",
#'         "GRCh38",
#'         "Ensembl",
#'         "EnsDb"
#'     )
#' )
#' mcols(ahs)
#' print(ahs)
#' # EnsDb (Ensembl GRCh38 94; 2018-10-11)
#' ah[["AH64923"]]
#'
#' # Human UCSC TxDb records.
#' ahs <- query(
#'     x = ah,
#'     pattern = c(
#'         "Homo sapiens",
#'         "UCSC",
#'         "TxDb",
#'         "knownGene"
#'     )
#' )
#' mcols(ahs)
#' print(ahs)
#' # TxDb (UCSC hg38 GENCODE 24; 2016-12-22)
#' ah[["AH52260"]]
#' ```
#'
#' @name makeGRangesFromEnsembl
#' @note Updated 2021-01-14.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `GRanges`.
#'
#' @seealso
#' - [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/).
#' - [GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures/).
#' - [ensembldb](https://bioconductor.org/packages/ensembldb/).
#' - `ensembldb::ensDbFromGtf()`.
#' - `GenomicFeatures::makeTxDbFromGFF()`.
#'
#' @examples
#' ## Genes
#' x <- makeGRangesFromEnsembl("Homo sapiens", level = "genes")
#' summary(x)
#'
#' ## Transcripts
#' x <- makeGRangesFromEnsembl("Homo sapiens", level = "transcripts")
#' summary(x)
NULL



#' Internal GRanges from Ensembl generator with additional options
#'
#' @note Updated 2021-01-14.
#' @noRd
.makeGRangesFromEnsembl <- function(
    organism,
    level,
    genomeBuild,
    release,
    ignoreVersion,
    synonyms,
    ## Internal-only arguments:
    broadClass = TRUE
) {
    assert(
        isFlag(ignoreVersion),
        isFlag(broadClass),
        isFlag(synonyms)
    )
    level <- match.arg(level)
    alert("Making {.var GRanges} from Ensembl.")
    edb <- getEnsDb(
        organism = organism,
        genomeBuild = genomeBuild,
        release = release
    )
    assert(
        is(edb, "EnsDb"),
        isString(attr(edb, "id"))
    )
    .makeGRangesFromEnsDb(
        object = edb,
        level = level,
        ignoreVersion = ignoreVersion,
        broadClass = broadClass,
        synonyms = synonyms
    )
}

formals(.makeGRangesFromEnsembl)[["level"]] <-
    formals(.makeGRangesFromEnsDb)[["level"]]



#' @rdname makeGRangesFromEnsembl
#' @export
makeGRangesFromEnsembl <- function(
    organism,
    level,
    genomeBuild = NULL,
    release = NULL,
    ignoreVersion = TRUE,
    synonyms = FALSE
) {
    .makeGRangesFromEnsembl(
        organism = organism,
        level = match.arg(level),
        genomeBuild = genomeBuild,
        release = release,
        ignoreVersion = ignoreVersion,
        synonyms = synonyms
    )
}

formals(makeGRangesFromEnsembl)[["level"]] <-
    formals(.makeGRangesFromEnsembl)[["level"]]



#' @describeIn makeGRangesFromEnsembl
#' Legacy convenience function that calls [makeGRangesFromEnsembl()] and returns
#' a `tbl_df` (tibble) instead of `GRanges`. Note that `GRanges` can also be
#' coerced using [`as.data.frame()`][base::as.data.frame].
#' @export
annotable <- function(...) {
        gr <- makeGRangesFromEnsembl(...)
        mcols(gr) <- decode(mcols(gr))
        as_tibble(gr, rownames = NULL)
    }
