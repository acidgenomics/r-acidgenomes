## FIXME Can we support keeping track of transcripts at exon level?
## This appears to work correctly for our GFF3 file parsing.
##
## FIXME Ensure we filter "LRG_" genes, as this doesn't match up with the GFF
## file conventions.



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
#' @note Updated 2023-12-20.
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
    function(organism,
             level = c("genes", "transcripts", "exons"),
             genomeBuild = NULL,
             release = NULL,
             ignoreVersion = FALSE,
             extraMcols = TRUE) {
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
    function(object,
             level = c("genes", "transcripts", "exons"),
             ignoreVersion = FALSE,
             extraMcols = TRUE) {
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
        if (isFALSE(ignoreVersion)) {
            idCol <- paste(idCol, "version", sep = "_")
        }
        quietly({
            cols <- ensembldb::listColumns(object, colKeys)
        })
        cols <- c(cols, "entrezid")
        cols <- sort(unique(cols))
        assert(
            isSubset(idCol, cols),
            msg = sprintf(
                paste(
                    "Unsupported {.pkg %s} key: {.var %s}.",
                    "Likely need to set {.val %s}."
                ),
                "ensembldb", idCol, "ignoreVersion = TRUE"
            )
        )
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
        ## Ensure we always remove LRG gene annotations, which are defined in
        ## ensembldb output but not in primary GFF3/GTF files.
        isLrg <- mcols(gr)[["gene_biotype"]] == "LRG_gene"
        if (any(isLrg)) {
            alertInfo(sprintf(
                "Removing %d LRG %s.",
                sum(isLrg),
                ngettext(
                    n = sum(isLrg),
                    msg1 = substr(
                        x = level,
                        start = 1L,
                        stop = nchar(level) - 1L
                    ),
                    msg2 = level
                )
            ))
            gr <- gr[!isLrg]
        }
        if (
            identical(level, "exons") &&
            hasDuplicates(mcols(gr)[["exon_id"]])
        ) {
            ## e.g. "ENSE00000000021" multimaps to multiple transcripts here
            ## but not in the GFF file definition...need to resolve.
            ## e.g. Homo sapiens exon "ENSE00001132905".
            ## Keep: "ENSG00000291317" (TMEM276).
            ## Drop: "ENSG00000291316" (no gene name; novel protein).
            dupes <- dupes(mcols(gr)[["exon_id"]])
            assert(length(dupes) <= 10L)
            alert(sprintf(
                "Resolving %d duplicate exon-to-gene %s: %s.",
                length(dupes),
                ngettext(
                    n = length(dupes),
                    msg1 = "mapping",
                    msg2 = "mappings"
                ),
                toInlineString(dupes)
            ))
            ## ensembldb currently returns empty columns instead of setting NA.
            keep <- !{
                mcols(gr)[["exon_id"]] %in% dupes &
                    nzchar(mcols(gr)[["gene_name"]])
            }
            gr <- gr[keep]
            assert(hasNoDuplicates(mcols(gr)[["exon_id"]]))
        }
        assert(hasNoDuplicates(mcols(gr)[[idCol]]))
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
