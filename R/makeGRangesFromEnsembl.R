## FIXME This isn't detaching Bioconductor packages correctly:
## AnnotationDbi, AnnotationFilter, Biobase, BiocGenerics, ensembldb,
## GenomeInfoDb, GenomicFeatures, GenomicRanges, IRanges, S4Vectors, stats4



#' Make GenomicRanges from Ensembl
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
#' - `coding`.
#' - `noncoding`.
#' - `pseudo`.
#' - `small`.
#' - `decaying`.
#' - `ig` (immunoglobulin).
#' - `tcr` (T cell receptor).
#' - `other`.
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
#' @name makeGRangesFromEnsembl
#' @note Updated 2023-04-12.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `EnsemblGenes` or `EnsemblTranscripts`.
#'
#' @seealso
#' - [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/).
#' - [ensembldb](https://bioconductor.org/packages/ensembldb/).
#' - `ensembldb::ensDbFromGff()`, `ensembldb::ensDbFromGtf()`.
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
#' if (goalie::isInstalled("EnsDb.Hsapiens.v75")) {
#'     genes <- makeGRangesFromEnsDb(
#'         object = "EnsDb.Hsapiens.v75",
#'         level = "genes"
#'     )
#'     summary(genes)
#' }
NULL



#' @describeIn makeGRangesFromEnsembl Obtain annotations from Ensembl by
#' querying AnnotationHub.
#' @export
makeGRangesFromEnsembl <-
    function(organism,
             level = c("genes", "transcripts"),
             genomeBuild = NULL,
             release = NULL,
             ignoreVersion = TRUE) {
        assert(isFlag(ignoreVersion))
        level <- match.arg(level)
        alert(sprintf("Making {.cls %s} from Ensembl.", "GenomicRanges"))
        edb <- .getEnsDb(
            organism = organism,
            genomeBuild = genomeBuild,
            release = release
        )
        gr <- makeGRangesFromEnsDb(
            object = edb,
            level = level,
            ignoreVersion = ignoreVersion
        )
        metadata(gr)[["call"]] <- tryCatch(
            expr = standardizeCall(),
            error = function(e) NULL
        )
        gr
    }



#' @describeIn makeGRangesFromEnsembl Use a specific `EnsDb` object as the
#' annotation source. Alternatively, can pass in an EnsDb package name as
#' a `character(1)`.
#' @export
#'
#' @param object `EnsDb` or `character(1)`.
#' `EnsDb` object or name of specific annotation package containing a
#' versioned EnsDb object (e.g. "EnsDb.Hsapiens.v75").
makeGRangesFromEnsDb <-
    function(object,
             level = c("genes", "transcripts"),
             ignoreVersion = TRUE) {
        pkgs <- .packages()
        assert(
            requireNamespaces("ensembldb"),
            isFlag(ignoreVersion)
        )
        level <- match.arg(level)
        alert(sprintf(
            "Making {.cls %s} from {.cls %s}.",
            "GenomicRanges", "EnsDb"
        ))
        if (isString(object)) {
            package <- object
            assert(requireNamespaces(package))
            object <- get(
                x = package,
                envir = asNamespace(package),
                inherits = FALSE
            )
        }
        assert(is(object, "EnsDb"))
        args <- list(
            "x" = object,
            "order.type" = "asc",
            "return.type" = "GRanges"
        )
        suppressWarnings({
            geneCols <- ensembldb::listColumns(object, "gene")
        })
        geneCols <- sort(unique(c(geneCols, "entrezid")))
        suppressWarnings({
            txCols <- ensembldb::listColumns(object, "tx")
        })
        txCols <- sort(unique(c(txCols, geneCols)))
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
                        "columns" = txCols,
                        "order.by" = "tx_id"
                    )
                )
            }
        )
        suppressWarnings({
            gr <- do.call(what = fun, args = args)
        })
        assert(is(gr, "GenomicRanges"))
        metadata(gr) <- .getEnsDbMetadata(object = object, level = level)
        gr <- .makeGRanges(
            object = gr,
            ignoreVersion = ignoreVersion
        )
        metadata(gr)[["call"]] <- tryCatch(
            expr = standardizeCall(),
            error = function(e) NULL
        )
        forceDetach(keep = pkgs)
        gr
    }
