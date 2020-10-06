## nolint start

#' Make GRanges from a GFF/GTF file
#'
#' @details
#' Remote URLs and compressed files are supported.
#'
#' @section Recommendations:
#'
#' - **Use GTF over GFF3.** We recommend using a GTF file instead of a GFF3
#'   file, when possible. The file format is more compact and easier to parse.
#' - **Use Ensembl over RefSeq.** We generally recommend using Ensembl over
#'   RefSeq, if possible. It's better supported in R and generally used by most
#'   NGS vendors.
#'
#' @section GFF/GTF specification:
#'
#' The GFF (General Feature Format) format consists of one line per feature,
#' each containing 9 columns of data, plus optional track definition lines.
#'
#' The GTF (General Transfer Format) format is identical to GFF version 2.
#'
#' The UCSC website has detailed conventions on the GFF3 format, including
#' the metadata columns.
#'
#' @section Feature type:
#'
#' - `CDS`: **C**o**D**ing **S**sequence. A contiguous sequence that contains a
#'   genomic interval bounded by start and stop codons. CDS refers to the
#'   portion of a genomic DNA sequence that is translated, from the start codon
#'   to the stop codon.
#' - `exon`: Genomic interval containing 5' UTR (`five_prime_UTR`), CDS, and
#'   3' UTR (`three_prime_UTR`).
#' - `mRNA`: Processed (spliced) mRNA transcript.
#'
#' See also:
#'
#' - [gffutils documentation](https://pythonhosted.org/gffutils/)
#' - [GenBank GFF documentation](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/)
#' - [stringtie GFF documentation](https://ccb.jhu.edu/software/stringtie/gff.shtml)
#' - [gmod.org GFF wiki](http://gmod.org/wiki/GFF)
#' - [Brent Lab GTF2 spec notes](https://mblab.wustl.edu/GTF2.html)
#' - [Sequence Ontology GFF3 spec notes](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
#'
#' @section Supported sources:
#'
#' Currently [makeGRangesFromGFF()] supports genomes from these sources:
#'
#' - Ensembl (GTF, GFF3).
#' - GENCODE (GTF, GFF3).
#' - RefSeq (GTF, GFF3).
#' - FlyBase (GTF).
#' - WormBase (GTF).
#'
#' @section Ensembl:
#'
#' Note that [makeGRangesFromEnsembl()] offers native support for Ensembl genome
#' builds and returns additional useful metadata that isn't defined inside a
#' GFF/GTF file.
#'
#' If you must load a GFF/GTF file directly, then use [makeGRangesFromGFF()].
#'
#' @section GENCODE vs. Ensembl:
#'
#' Annotations available from Ensembl and GENCODE are very similar.
#'
#' The GENCODE annotation is made by merging the manual gene annotation produced
#' by the Ensembl-Havana team and the Ensembl-genebuild automated gene
#' annotation. The GENCODE annotation is the default gene annotation displayed
#' in the Ensembl browser. The GENCODE releases coincide with the Ensembl
#' releases, although GENCODE can skip an Ensembl release if there is no update
#' to the annotation with respect to the previous release. In practical terms,
#' the GENCODE annotation is essentially identical to the Ensembl annotation.
#'
#' However, GENCODE handles pseudoautosomal regions (PAR) differently than
#' Ensembl. The Ensembl GTF file only includes this annotation once, for
#' chromosome X. However, GENCODE GTF/GFF3 files include the annotation in the
#' PAR regions of both chromosomes. You'll see these genes contain a "_PAR_Y"
#' suffix.
#'
#' Additionally, GENCODE GFF/GTF files import with a gene identifier containing
#' a suffix, which differs slightly from the Ensembl GFF/GTF spec
#' (e.g. GENCODE: `ENSG00000000003.14`; Ensembl: `ENSG00000000003`).
#'
#' The [GENCODE FAQ](https://www.gencodegenes.org/pages/faq.html) has additional
#' details.
#'
#' @section RefSeq:
#'
#' Refer to the
#' [current RefSeq spec](ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt)
#' for details.
#'
#' See also:
#'
#' - [RefSeq FAQ](https://www.ncbi.nlm.nih.gov/books/NBK50679/)
#' - ftp://ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz
#'
#' @section UCSC:
#'
#' Loading UCSC genome annotations from a GFF/GTF file are
#' *intentionally not supported* by this function.
#'
#' We recommend using a pre-built `TxDb` package from Bioconductor instead.
#' For example, load `TxDb.Hsapiens.UCSC.hg38.knownGene` for hg38.
#'
#' For reference, note that UCSC doesn't provide direct GFF/GTF file downloads.
#' Use of the [hgTables](https://genome.ucsc.edu/cgi-bin/hgTables) table
#' browser is required in a web browser.
#'
#' Select the following options to download hg38:
#'
#' - clade: `Mammal`
#' - genome: `Human`
#' - assembly: `Dec. 2013 (GRCh38/hg38)`
#' - group: `Genes and Gene Predictions`
#' - track: `GENCODE v29`
#' - table: `knownGene`
#' - region: `genome`
#' - output format: `GTF - gene transfer format`
#' - output file: `<Enter a file name>`
#'
#' Related URLs:
#'
#' - [UCSC hgTables](http://genome.ucsc.edu/cgi-bin/hgTables)
#' - [UCSC downloads](http://hgdownload.soe.ucsc.edu/downloads.html)
#' - [UCSC hg38 FTP](ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/)
#'
#' @section Example URLs:
#'
#' - Ensembl *Homo sapiens* 95
#'   [GTF](ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz),
#'   [GFF3](ftp://ftp.ensembl.org/pub/release-95/gff3/homo_sapiens/Homo_sapiens.GRCh38.95.gff3.gz)
#' - GENCODE *Homo sapiens* 29
#'   [GTF](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz),
#'   [GFF3](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz)
#' - RefSeq *Homo sapiens*
#'   [GTF and GFF3](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/)
#' - FlyBase *Drosophila melanogaster* r6.24
#'   [GTF](ftp://ftp.flybase.net/releases/FB2018_05/dmel_r6.24/gtf/dmel-all-r6.24.gtf.gz)
#' - WormBase *Caenorhabditis elegans* WS267
#'   [GTF](ftp://ftp.wormbase.org/pub/wormbase/releases/WS267/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS267.canonical_geneset.gtf.gz)
#'
#' @export
#' @note Updated 2020-10-05.
#'
#' @inheritParams AcidRoxygen::params
#' @param .checkAgainstTxDb `logical(1)`.
#'   Enable strict mode, intended for development and unit testing only.
#'   Generate an internal `TxDb` using [GenomicFeatures::makeTxDbFromGRanges()]
#'   and check that the [`ranges()`][IRanges::ranges],
#'   [`seqnames()`][GenomeInfoDb::seqnames], and identifiers defined in
#'   [`names()`][base::names] are identical. Doesn't work for all GFF/GTF files
#'   due to some current limitations in the GenomicFeatures package, so this is
#'   disabled by default. Generally, GenomicFeatures parses GTF files better
#'   than GFF files. However, it's a useful sanity check and should be enabled
#'   if possible.
#'
#' @return `GRanges`.
#'
#' @seealso
#' - [rtracklayer::import()].
#' - [GenomicFeatures::makeTxDbFromGRanges()].
#' - [GenomicFeatures::makeTxDbFromGFF()].
#'
#' @examples
#' file <- pasteURL(basejumpTestsURL, "ensembl.gtf", protocol = "none")
#'
#' ## Genes
#' x <- makeGRangesFromGFF(file = file, level = "genes")
#' summary(x)
#'
#' ## Transcripts
#' x <- makeGRangesFromGFF(file = file, level = "transcripts")
#' summary(x)

## nolint end

makeGRangesFromGFF <- function(
    file,
    level = c("genes", "transcripts"),
    ignoreTxVersion = TRUE,
    synonyms = TRUE,
    .checkAgainstTxDb = FALSE
) {
    assert(
        isString(file),
        isFlag(ignoreTxVersion),
        isFlag(synonyms),
        isFlag(.checkAgainstTxDb)
    )
    level <- match.arg(level)
    cli_alert(sprintf(
        fmt = "Making {.var GRanges} from GFF file ({.file %s}).",
        basename(file)
    ))

    ## Import ------------------------------------------------------------------
    ## This step uses `rtracklayer::import()` internally.
    ## Note that this can generate very large objects from GFF3.
    object <- import(file)
    ## Slot the source (e.g. Ensembl) and type (e.g. GTF) into `metadata()`.
    object <- .slotGFFDetectInfo(object)
    ## Note that this metadata is used by TxDb caller for GFF3 (see below).
    metadata(object)[["file"]] <- file
    ## Pull detection strings from GRanges `metadata()`.
    detect <- metadata(object)[["detect"]]
    assert(is.character(detect))
    source <- detect[["source"]]
    type <- detect[["type"]]
    assert(isString(source), isString(type))

    ## Pre-flight checks -------------------------------------------------------
    ## Not currently allowing FlyBase or WormBase GFF files. They're too
    ## complicated to parse, and the sites offer GTF files instead anyway.
    if (source %in% c("FlyBase", "WormBase") && type != "GTF") {
        stop(sprintf("Only GTF files from %s are supported.", source))  # nocov
    }

    ## TxDb (GenomicFeatures) --------------------------------------------------
    ## Run this step prior to any GRanges sanitization steps.
    ## This check may be removed in a future update.
    if (isTRUE(.checkAgainstTxDb)) {
        ## nocov start
        cli_alert_warning(
            "Strict mode enabled. Checking against {.var TxDb}."
        )
        txdb <- .makeTxDbFromGFF(object)
        ## nocov end
    } else {
        txdb <- NULL
    }

    ## Standardize -------------------------------------------------------------
    ## Standardize FlyBase, GENCODE, and RefSeq files to follow expected
    ## Ensembl-like naming conventions. This step must be run after
    ## `.makeTxDbFromGFF()` but prior to `.makeGenesFromGFF()` and/or
    ## `.makeTranscriptsFromGFF()` calls (see below).
    object <- switch(
        EXPR     = source,
        FlyBase  = .standardizeFlyBaseToEnsembl(object),
        GENCODE  = .standardizeGencodeToEnsembl(object),
        RefSeq   = .standardizeRefSeqToEnsembl(object),
        object
    )
    mcolnames <- colnames(mcols(object))
    assert(
        isSubset(x = c("gene_id", "transcript_id"), y = mcolnames),
        ## `gene_type` needs to be renamed to `gene_biotype`, if defined.
        areDisjointSets(x = c("gene", "gene_type"), y = mcolnames)
    )
    rm(mcolnames)
    genes <- object
    if (level == "transcripts") {
        transcripts <- object
    }
    rm(object)

    ## Genes -------------------------------------------------------------------
    ## `makeGRangesFromGFF()` attempts to always returns gene-level metadata,
    ## even when transcripts are requested. We'll merge this object into the
    ## transcript-level GRanges below, when possible.
    genes <- genes[!is.na(mcols(genes)[["gene_id"]])]
    genes <- genes[is.na(mcols(genes)[["transcript_id"]])]
    assert(hasLength(genes))
    if (source == "Ensembl" && type == "GFF3") {
        genes <- .makeGenesFromEnsemblGFF3(genes)
    } else if (source == "Ensembl" && type == "GTF") {
        genes <- .makeGenesFromEnsemblGTF(genes)
    } else if (source == "FlyBase" && type == "GTF") {
        genes <- .makeGenesFromFlyBaseGTF(genes)
    } else if (source == "GENCODE" && type == "GFF3") {
        genes <- .makeGenesFromGencodeGFF3(genes)
    } else if (source == "GENCODE" && type == "GTF") {
        genes <- .makeGenesFromGencodeGTF(genes)
    } else if (source == "RefSeq" && type == "GFF3") {
        genes <- .makeGenesFromRefSeqGFF3(genes)
    } else if (source == "RefSeq" && type == "GTF") {
        genes <- .makeGenesFromRefSeqGTF(genes)
    } else if (source == "WormBase" && type == "GTF") {
        genes <- .makeGenesFromWormBaseGTF(genes)
    } else {
        ## nocov start
        stop(
            "Failed to make gene-level GRanges.\n",
            "Unsupported GFF source file."
        )
        ## nocov end
    }
    ## Remove GFF-specific parent columns, etc.
    if (type == "GFF3") {
        genes <- .minimizeGFF3(genes)
    }
    ## Drop columns that contain all `NA`.
    mcols(genes) <- removeNA(mcols(genes))
    ## Set names and stash metadata.
    names(genes) <- mcols(genes)[["gene_id"]]
    metadata(genes)[["level"]] <- "genes"
    if (level == "genes") {
        out <- genes
    }

    ## Transcripts -------------------------------------------------------------
    if (level == "transcripts") {
        transcripts <-
            transcripts[!is.na(mcols(transcripts)[["transcript_id"]])]
        assert(hasLength(transcripts))
        if (source == "Ensembl" && type == "GFF3") {
            transcripts <- .makeTranscriptsFromEnsemblGFF3(transcripts)
        } else if (source == "Ensembl" && type == "GTF") {
            transcripts <- .makeTranscriptsFromEnsemblGTF(transcripts)
        } else if (source == "FlyBase" && type == "GTF") {
            transcripts <- .makeTranscriptsFromFlyBaseGTF(transcripts)
        } else if (source == "GENCODE" && type == "GFF3") {
            transcripts <- .makeTranscriptsFromGencodeGFF3(transcripts)
        } else if (source == "GENCODE" && type == "GTF") {
            transcripts <- .makeTranscriptsFromGencodeGTF(transcripts)
        } else if (source == "RefSeq" && type == "GFF3") {
            transcripts <- .makeTranscriptsFromRefSeqGFF3(transcripts)
        } else if (source == "RefSeq" && type == "GTF") {
            transcripts <- .makeTranscriptsFromRefSeqGTF(transcripts)
        } else if (source == "WormBase" && type == "GTF") {
            transcripts <- .makeTranscriptsFromWormBaseGTF(transcripts)
        } else {
            ## nocov start
            stop(
                "Failed to make transcript-level GRanges.\n",
                "Unsupported GFF file format."
            )
            ## nocov end
        }
        ## Remove GFF-specific parent columns, etc.
        if (type == "GFF3") {
            transcripts <- .minimizeGFF3(transcripts)
        }
        ## Drop columns that contain all NA.
        ## This step is necessary for gene-level merge step to work.
        ## Relatively slow for large datasets, and could use a speed boost.
        mcols(transcripts) <- removeNA(mcols(transcripts))
        ## Set names and stash metadata.
        names(transcripts) <- mcols(transcripts)[["transcript_id"]]
        metadata(transcripts)[["level"]] <- "transcripts"
        ## Skip gene-level metadata merge for GRanges that have been split
        ## into GRangesList.
        if (
            is(genes, "GRanges") &&
            ## This step is necessary for RefSeq GFF3.
            !anyDuplicated(mcols(genes)[["gene_id"]])
        ) {
            ## By default, merge the gene-level annotations into the
            ## transcript-level ones, for objects that have ranges 1:1 with the
            ## identifiers.
            out <- .mergeGenesIntoTranscripts(transcripts, genes)
        } else {
            cli_alert_warning("Skipping gene metadata merge.")
            out <- transcripts
        }
    }

    ## Return ------------------------------------------------------------------
    out <- .makeGRanges(
        object = out,
        ignoreTxVersion = ignoreTxVersion,
        synonyms = synonyms
    )
    ## Ensure that the ranges match GenomicFeatures output, if desired.
    ## Note that this step will error out for WormBase currently. Need to think
    ## of a reworked approach that avoids this.
    if (is(out, "GRanges") && is(txdb, "TxDb")) {
        ## nocov start
        assert(isTRUE(ignoreTxVersion))
        .checkGRangesAgainstTxDb(gr = out, txdb = txdb)
        ## nocov end
    }
    out
}



## Aliases =====================================================================
#' @describeIn makeGRangesFromGFF GTF file extension alias.
#'   Runs the same internal code as [makeGRangesFromGFF()].
#' @export
makeGRangesFromGTF <- makeGRangesFromGFF
