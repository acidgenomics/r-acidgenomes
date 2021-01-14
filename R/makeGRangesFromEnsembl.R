## FIXME RETHINK METADATA RETURN STRUCTURE.
## FIXME NEED TO ENSURE CALL IS SLOTTED INTO OBJECT.



#' Make GRanges from Ensembl
#'
#' Quickly obtain gene and transcript annotations from
#' [Ensembl](https://www.ensembl.org/) using
#' [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/) and
#' [ensembldb](https://bioconductor.org/packages/ensembldb/).
#'
#' Simply specify the desired organism, using the full Latin name. For example,
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
#' - [ensembldb](https://bioconductor.org/packages/ensembldb/).
#' - `ensembldb::ensDbFromGtf()`.
#'
#' @examples
#' ## Get annotations from Ensembl via AnnotationHub query.
#' genes <- makeGRangesFromEnsembl(
#'     organism = "Homo sapiens",
#'     level = "genes"
#' )
#' summary(genes)
#' transcripts <- makeGRangesFromEnsembl(
#'     organism = "Homo sapiens",
#'     level = "transcripts"
#' )
#' summary(transcripts)
#'
#' ## Get annotations from specific EnsDb object or package.
#' if ("EnsDb.Hsapiens.v75" %in% rownames(installed.packages())) {
#'     genes <- makeGRangesFromEnsDb(
#'         object = "EnsDb.Hsapiens.v75",
#'         level = "genes"
#'     )
#'     summary(genes)
#' }
NULL



#' Make GRanges from EnsDb object
#'
#' Internal variant with more options that we don't want to expose to user.
#'
#' @note Updated 2021-01-14.
#' @noRd
.makeGRangesFromEnsDb <- function(
    object,
    level,
    ignoreVersion,
    synonyms,
    ## Internal-only arguments:
    broadClass = TRUE
) {
    assert(
        isFlag(ignoreVersion),
        isFlag(synonyms),
        isFlag(broadClass)
    )
    level <- match.arg(level)
    alert("Making {.var GRanges} from {.var EnsDb}.")
    if (isString(object)) {
        package <- object
        requireNamespaces(package)
        object <- get(
            x = package,
            envir = asNamespace(package),
            inherits = FALSE
        )
    }
    assert(is(object, "EnsDb"))
    metadata <- .getEnsDbMetadata(object, level = level)
    args <- list(
        "x" = object,
        "order.type" = "asc",
        "return.type" = "GRanges"
    )
    ## FIXME USE LIST COLUMNS HERE INSTEAD.
    geneCols <- listColumns(object, "gene")
    geneCols <- c(
        "gene_id",
        "gene_name",
        "gene_biotype",
        "seq_coord_system",
        "entrezid"
    )
    switch(
        EXPR = level,
        "genes" = {
            fun <- ensembldb::genes
            args <- append(
                x = args,
                values = list(
                    "columns" = geneCols,
                    "order.by" = "gene_id"
                )
            )
        },
        "transcripts" = {
            fun <- ensembldb::transcripts
            args <- append(
                x = args,
                values = list(
                    "columns" = c(
                        "tx_id",
                        "tx_name",
                        "tx_biotype",
                        "tx_cds_seq_start",
                        "tx_cds_seq_end",
                        geneCols
                    ),
                    "order.by" = "tx_id"
                )
            )
        }
    )
    ## This step can warn about out-of-bound ranges that need to be trimmed.
    ## We're taking care of trimming on the `.makeGRanges` call below.
    suppressWarnings({
        gr <- do.call(what = fun, args = args)
    })
    assert(is(gr, "GRanges"))
    metadata(gr) <- metadata
    .makeGRanges(
        object = gr,
        ignoreVersion = ignoreVersion,
        broadClass = broadClass,
        synonyms = synonyms
    )
}



#' Make GRanges from Ensembl via AnnotationHub query
#'
#' Internal variant with more options that we don't want to expose to user.
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



#' @describeIn makeGRangesFromEnsembl Obtain annotations from Ensembl by
#'   querying AnnotationHub.
#' @export
makeGRangesFromEnsembl <- function(
    organism,
    level = c("genes", "transcripts"),
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



#' @describeIn makeGRangesFromEnsembl Use a specific `EnsDb` object as the
#'   annotation source. Alternatively, can pass in an EnsDb package name as
#'   a `character(1)`.
#' @export
#'
#' @param object `EnsDb` or `character(1)`.
#'   `EnsDb` object or name of specific annotation package containing a
#'   versioned EnsDb object (e.g. "EnsDb.Hsapiens.v75").
makeGRangesFromEnsDb <- function(
    object,
    level = c("genes", "transcripts"),
    ignoreVersion = TRUE,
    synonyms = FALSE
) {
    .makeGRangesFromEnsDb(
        object = object,
        level = match.arg(level),
        ignoreVersion = ignoreVersion,
        synonyms = synonyms
    )
}



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
