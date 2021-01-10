## FIXME ENSURE REFSEQ TRANSCRIPTS RETURN AS FLAT GRANGES OBJECT.
## FIXME TEST FLYBASE GFF AND WORMBASE GFF.

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
#' - Ensembl *Homo sapiens* GRCh38 102
#'   [GTF](ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz),
#'   [GFF3](ftp://ftp.ensembl.org/pub/release-102/gff3/homo_sapiens/Homo_sapiens.GRCh38.102.gff3.gz)
#' - GENCODE *Homo sapiens* GRCh38 32
#'   [GTF](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz),
#'   [GFF3](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gff3.gz)
#' - RefSeq *Homo sapiens* GCF_000001405.39 GRCh38
#'   [GTF](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz)
#'   [GFF3](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz)
#' - FlyBase *Drosophila melanogaster* r6.24
#'   [GTF](ftp://ftp.flybase.net/releases/FB2020_06/dmel_r6.37/gtf/dmel-all-r6.37.gtf.gz)
#' - WormBase *Caenorhabditis elegans* WS267
#'   [GTF](ftp://ftp.wormbase.org/pub/wormbase/releases/WS279/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.279.canonical_geneset.gtf.gz)
#'
#' @export
#' @note Updated 2021-01-10.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
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
#' - `tximeta:::getRanges()`.
#' - `columns()`.
#'
#' @examples
#' file <- pasteURL(AcidGenomesTestsURL, "ensembl.gtf")
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
    broadClass = TRUE,
    synonyms = FALSE
) {
    assert(
        isString(file),
        isFlag(ignoreTxVersion),
        isFlag(broadClass),
        isFlag(synonyms)
    )
    level <- match.arg(level)
    alert(sprintf(
        fmt = "Making {.var GRanges} from GFF file ({.file %s}).",
        basename(file)
    ))
    if (isAURL(file)) {
        tmpfile <- cacheURL(url = file)
    } else {
        tmpfile <- file
    }
    detect <- .detectGFF(tmpfile)
    meta[["detect"]] <- detect
    source <- detect[["source"]]
    type <- detect[["type"]]
    assert(isString(source), isString(type))
    ## Not currently allowing FlyBase or WormBase GFF files. They're too
    ## complicated to parse, and the sites offer GTF files instead anyway.
    ## FIXME CONSIDER TAKING THIS CHECK OUT...
    if (isSubset(source, c("FlyBase", "WormBase")) && type != "GTF") {
        stop(sprintf("Only GTF files from %s are supported.", source))  # nocov
    }
    ## Generate EnsDb or TxDb from input GFF/GTF file.
    ## Use `columns()` on EnsDb or TxDb to check for available metadata.
    args <- list()
    if (isSubset(source, c("Ensembl", "GENCODE"))) {
        db <- .makeEnsDbFromGFF(tmpfile)
        switch(
            EXPR = level,
            "genes" = {
                fun <- ensembldb::genes
                orderBy <- "gene_id"
                columns <- c(
                    "gene_id",
                    "gene_name",
                    "gene_biotype"
                )
            },
            "transcripts" = {
                fun <- ensembldb::transcripts
                orderBy <- "transcript_id"
                columns <- c(
                    "tx_id",
                    "tx_name",
                    "tx_biotype",
                    "gene_id",
                    "gene_name",
                    "gene_biotype"
                )
            }
        )
        args[["order.by"]] <- orderBy
        args[["return.type"]] <- "GRanges"
    } else {
        ## FIXME CHECK AGAINST COLUMNS AND DETERMINE WHICH ONES ARE OK TO USE.
        ## FIXME SOME ANNOTATION TYPES MAY SUPPORT GENE_NAME AND BIOTYPE...
        db <- .makeTxDbFromGFF(tmpfile)
        switch(
            EXPR = level,
            "genes" = {
                fun <- GenomicFeatures::genes
                columns <- "gene_id"
            },
            "transcripts" = {
                fun <- GenomicFeatures::transcripts
                columns <- c("tx_id", "tx_name", "gene_id")
            }
        )
    }
    args[["x"]] <- db
    args[["columns"]] <- columns
    gr <- do.call(what = fun, args = args)
    assert(is(gr, "GRanges"))
    out <- .makeGRanges(
        object = gr,
        ignoreTxVersion = ignoreTxVersion,
        broadClass = broadClass,
        synonyms = synonyms
    )
    ## FIXME WE NEED TO DECLARE WHICH PACKAGE GENERATED THIS RANGES.
    ## FIXME THIS NEEDS TO INCLUDE ORGANISM.
    metadata(out)[["file"]] <- file
    metadata(out)[["call"]] <- match.call()
    out
}



## Aliases =====================================================================
#' @describeIn makeGRangesFromGFF GTF file extension alias.
#'   Runs the same internal code as [makeGRangesFromGFF()].
#' @export
makeGRangesFromGTF <- makeGRangesFromGFF
