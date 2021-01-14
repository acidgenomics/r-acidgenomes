## FIXME REWORK, CALLING makeTxDbFromGFF internally.
## FIXME ENSURE REFSEQ TRANSCRIPTS RETURN AS FLAT GRANGES OBJECT.
## FIXME TEST FLYBASE GFF AND WORMBASE GFF.
## FIXME INCLUDE GENEVERSION HERE IF POSSIBLE WHEN `IGNOREVERSION` = FALSE
## FIXME RETHINK ALLOWING BROADCLASS AND SYNONYMS HERE.
## FIXME CURRENT RELEASE VERSION DOESNT SLOT ORGANISM HERE CORRECTLY.
##       RETHINK THAT FOR GFF.
## FIXME MAKE SURE FILE IS CORRECT URL, NOT TMPFILE BEFORE RELEASING.
## FIXME synonyms only works with Ensembl identifiers, consider making that more
##       clear in documentation.
## FIXME INDICATE TO THE USER MORE CLEARLY THAT THIS STEP IS SLOW.



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
#' - RefSeq *Homo sapiens* GRCh38.p12
#'   [GFF3](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gff.gz)
#' - FlyBase *Drosophila melanogaster* r6.24
#'   [GTF](ftp://ftp.flybase.net/releases/FB2020_06/dmel_r6.37/gtf/dmel-all-r6.37.gtf.gz)
#' - WormBase *Caenorhabditis elegans* WS267
#'   [GTF](ftp://ftp.wormbase.org/pub/wormbase/releases/WS279/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.279.canonical_geneset.gtf.gz)
#'
#' @section Exons vs. CDS:
#' - Exons: `gene - introns`.
#' - CDS: `exons - UTRs`.
#'
#' @name makeGRangesFromGFF
#' @note Updated 2021-01-14.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `GRanges`.
#'
#' @seealso
#' - `rtracklayer::import()`.
#' - `GenomicFeatures::makeTxDbFromGRanges()`.
#' - `GenomicFeatures::makeTxDbFromGFF()`.
#' - `GenomicFeatures:::.make_splicings()`.
#' - `tximeta:::getRanges()`.
#' - `AnnotationDbi::columns()`.
#' - `GenomeInfoDb::GenomeDescription-class`, which describes useful `organism`,
#'   `commonName`, `providerVersion`, `provider`, and `releaseDate` accessors.
#'
#' @examples
#' file <- pasteURL(AcidGenomesTestsURL, "ensembl.gtf")
#'
#' ## Genes.
#' genes <- makeGRangesFromGFF(file = file, level = "genes")
#' summary(genes)
#'
#' ## Transcripts.
#' transcripts <- makeGRangesFromGFF(file = file, level = "transcripts")
#' summary(transcripts)
NULL

## nolint end



#' Make GRanges from a GFF file
#'
#' Internal variant with more options that we don't want to expose to user.
#'
#' @note Updated 2021-01-14.
#' @noRd
.makeGRangesFromGFF <- function(
    file,
    level = c("genes", "transcripts"),
    ignoreVersion = TRUE,
    synonyms = FALSE,
    ## Internal-only arguments:
    broadClass = TRUE
) {
    assert(
        isString(file),
        isFlag(ignoreVersion),
        isFlag(synonyms),
        isFlag(broadClass)
    )
    level <- match.arg(level)
    alert(sprintf(
        fmt = "Making {.var GRanges} from GFF file ({.file %s}).",
        basename(file)
    ))
    tmpfile <- .cacheIt(file)
    ## Load raw GFF/GTF ranges into memory using `rtracklayer::import()`. We're
    ## using this downstream for file source detection and extra metadata that
    ## currently isn't supported in GenomicFeatures TxDb generation.
    rawGRanges <- import(tmpfile)
    detect <- .detectGFF(rawGRanges)
    source <- detect[["source"]]
    type <- detect[["type"]]
    assert(isString(source), isString(type))
    if (isSubset(source, c("FlyBase", "WormBase")) && type != "GTF") {
        stop(sprintf("Only GTF files from %s are supported.", source))  # nocov
    } else if (source == "RefSeq" && type != "GFF") {
        ## https://github.com/Bioconductor/GenomicFeatures/issues/26
        stop(sprintf("Only GFF files from %s are supported.", source))  # nocov
    }
    ## Use ensembldb for Ensembl and GENCODE files, otherwise handoff to
    ## GenomicFeatures and generate a TxDb object.
    if (isSubset(source, c("Ensembl", "GENCODE"))) {
        db <- makeEnsDbFromGFF(tmpfile)
        gr <- .makeGRangesFromEnsDb(
            object = db,
            level = level,
            ignoreVersion = ignoreVersion,
            broadClass = broadClass,
            synonyms = synonyms
        )
    } else {
        if (isTRUE(synonyms)) {
            stop(paste(
                "Synonyms are only supported for genomes from",
                "Ensembl and GENCODE."
            ))
        }
        db <- .makeTxDbFromGFF(tmpfile)
        ## FIXME 47 sequences (1 circular) from an unspecified genomes; no seqlengths...argh
        ## FIXME WE NEED TO IMPROVE THIS METADATA...
        gr1 <- .makeGRangesFromTxDb(object = db, level = level)
        gr2

        ## FIXME FOR REFSEQ USE ASSEMBLY_REPORT
        ## tximeta:::gtf2RefSeq

        ## FIXME Need to define seqlengths here.
        ## FIXME NEED TO ADD RICHER METADATA HERE...
        ## FIXME NEED TO SLOT SEQLENGTHS HERE.

        out <- .makeGRanges(
            object = gr,
            ignoreVersion = ignoreVersion,
            broadClass = broadClass,
            synonyms = synonyms
        )
    }
    assert(is(gr, "GRanges"))

    ## FIXME ATTEMPT TO SLOT THE GENOME BUILD FROM THE FILE NAME HERE.
    ## FIXME WE NEED TO DECLARE WHICH PACKAGE GENERATED THIS RANGES.
    ## FIXME THIS NEEDS TO INCLUDE ORGANISM.
    metadata(out)[["detect"]] <- detect
    metadata(out)[["file"]] <- file
    metadata(out)[["call"]] <- match.call()

    ## FIXME RETHINK THIS.
    ## Metadata assert checks before return.
    if (!isSubset(source, "FlyBase")) {
        seqinfo(gr)
        seqlengths(gr)
        genome(gr)
    }

    out
}



#' @describeIn makeGRangesFromGFF Primary function.
#' @export
makeGRangesFromGFF <- function(
    file,
    level = c("genes", "transcripts"),
    ignoreVersion = TRUE,
    synonyms = FALSE
) {
    .makeGRangesFromGFF(
        file = file,
        level = match.arg(level),
        ignoreVersion = ignoreVersion,
        synonyms = synonyms
    )
}



#' @describeIn makeGRangesFromGFF Alias for GTF files.
#' @export
makeGRangesFromGTF <- makeGRangesFromGFF
