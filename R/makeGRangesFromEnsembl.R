## FIXME Now returning GRCh37 112 by default instead of expected GRCh38, which is behind at 111.
## FIXME Need to think of approach for dealing with this edge case.
## FIXME With GRCh37 112 update, we may be able to sunset legacy EnsDb.Hsapiens.v75 usage.

## nolint start
#' Make genomic ranges (`GRanges`) from Ensembl
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
#' @section Ensembl canonical transcript:
#'
#' Fetching the canonical transcripts currently requires connecting to the
#' [Ensembl Perl API][]. This information is currently processed by the
#' ensembldb package using the `fetchTablesFromEnsembl` function.
#'
#' See also:
#' - http://useast.ensembl.org/info/genome/genebuild/canonical.html
#' - https://www.ensembl.info/2021/04/16/update-to-the-ensembl-canonical-transcript-set/
#' - https://github.com/jorainer/ensembldb/blob/devel/inst/perl/get_gene_transcript_exon_tables.pl
#' - https://github.com/jorainer/ensembldb/issues/123
#' - https://github.com/jorainer/ensembldb/blob/devel/R/functions-create-EnsDb.R#L36
#' - `ensembldb::fetchTablesFromEnsembl()`.
#'
#' [Ensembl Perl API]: http://useast.ensembl.org/info/docs/api/index.html
#'
#' @name makeGRangesFromEnsembl
#' @note Updated 2025-04-12.
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
#' - [Ensembl biotypes](https://useast.ensembl.org/info/genome/genebuild/biotypes.html).
#' - [Gene/transcript biotypes in GENCODE and Ensembl](https://www.gencodegenes.org/pages/biotypes.html).
#' - [Locus reference genomic](http://www.lrg-sequence.org/).
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
#' ## > if (goalie::isInstalled("EnsDb.Hsapiens.v75")) {
#' ## >     genes <- makeGRangesFromEnsDb(
#' ## >         object = "EnsDb.Hsapiens.v75",
#' ## >         level = "genes"
#' ## >     )
#' ## >     summary(genes)
#' ## > }
NULL
## nolint end

#' @describeIn makeGRangesFromEnsembl Obtain annotations from Ensembl by
#' querying AnnotationHub.
#' @export
makeGRangesFromEnsembl <-
    function(
        organism,
        level = c("genes", "transcripts", "exons"),
        genomeBuild = NULL,
        release = NULL,
        ignoreVersion = FALSE,
        extraMcols = FALSE
    ) {
        assert(
            isFlag(ignoreVersion),
            isFlag(extraMcols)
        )
        level <- match.arg(level)
        alert(sprintf("Making {.cls %s} from Ensembl.", "GRanges"))
        edb <- getEnsDb(
            organism = organism,
            genomeBuild = genomeBuild,
            release = release
        )
        gr <- makeGRangesFromEnsDb(
            object = edb,
            level = level,
            ignoreVersion = ignoreVersion,
            extraMcols = extraMcols
        )
        metadata(gr)[["call"]] <- tryCatch(
            expr = standardizeCall(),
            error = function(e) {
                NULL
            }
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
    function(
        object,
        level = c("genes", "transcripts", "exons"),
        ignoreVersion = FALSE,
        extraMcols = FALSE
    ) {
        assert(
            requireNamespaces("ensembldb"),
            isFlag(ignoreVersion),
            isFlag(extraMcols)
        )
        level <- match.arg(level)
        alert(sprintf("Making {.cls %s} from {.cls %s}.", "GRanges", "EnsDb"))
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
        switch(
            EXPR = level,
            "exons" = {
                fun <- ensembldb::exons
                colKeys <- c("exon", "gene", "tx")
                idCol <- "exon_id"
            },
            "genes" = {
                fun <- ensembldb::genes
                colKeys <- "gene"
                idCol <- "gene_id"
            },
            "transcripts" = {
                fun <- ensembldb::transcripts
                colKeys <- c("gene", "tx")
                idCol <- "tx_id"
            }
        )
        quietly({
            cols <- ensembldb::listColumns(object, colKeys)
        })
        cols <- c(cols, "entrezid")
        cols <- sort(unique(cols))
        if (isFALSE(ignoreVersion)) {
            idVerCol <- paste(idCol, "version", sep = "_")
            if (isSubset(idVerCol, cols)) {
                idCol <- idVerCol
            }
        }
        args <- list(
            "x" = object,
            "columns" = cols,
            "order.by" = idCol,
            "order.type" = "asc",
            "return.type" = "GRanges"
        )
        quietly({
            gr <- do.call(what = fun, args = args)
        })
        assert(is(gr, "GRanges"))
        metadata(gr) <- .getEnsDbMetadata(object = object, level = level)
        gr <- .makeGRanges(
            object = gr,
            ignoreVersion = ignoreVersion,
            extraMcols = extraMcols
        )
        metadata(gr)[["call"]] <- tryCatch(
            expr = standardizeCall(),
            error = function(e) {
                NULL
            }
        )
        gr
    }
