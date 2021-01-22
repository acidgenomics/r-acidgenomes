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
#' @name makeGRangesFromGFF
#' @note Updated 2021-01-20.
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
#' - Ensembl
#' - GENCODE
#' - RefSeq
#' - UCSC
#' - FlyBase
#' - WormBase
#'
#' @section Ensembl:
#'
#' Note that [makeGRangesFromEnsembl()] offers native support for Ensembl genome
#' builds and returns additional useful metadata that isn't defined inside a
#' GFF/GTF file.
#'
#' If you must load a GFF/GTF file directly, then use [makeGRangesFromGFF()].
#'
#' Example URLs:
#'
#' - Ensembl *Homo sapiens* GRCh38.p13, release 102
#'   [GTF](ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz),
#'   [GFF3](ftp://ftp.ensembl.org/pub/release-102/gff3/homo_sapiens/Homo_sapiens.GRCh38.102.gff3.gz)
#' - Ensembl *Homo sapiens* GRCh37, release 102 (87)
#'   [GTF](ftp://ftp.ensembl.org/pub/grch37/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz),
#'   [GFF3](ftp://ftp.ensembl.org/pub/grch37/release-102/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz)
#'
#' @section GENCODE:
#'
#' Example URLs:
#'
#' - GENCODE *Homo sapiens* GRCh38.p13, release 36
#'   [GTF](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz),
#'   [GFF3](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gff3.gz)
#' - GENCODE *Homo sapiens* GRCh37, release 36
#'   [GTF](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh37_mapping/gencode.v36lift37.annotation.gtf.gz),
#'   [GFF3](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh37_mapping/gencode.v36lift37.annotation.gff3.gz)
#' - GENCODE *Mus musculus* GRCm38.p6, release M25
#'   [GTF](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz),
#'   [GFF3](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gff3.gz)
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
#' Example URLs:
#'
#' - RefSeq *Homo sapiens* GRCh38.p12
#'   [GTF](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gtf.gz),
#'   [GFF3](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gff.gz)
#'
#' See also:
#'
#' - [RefSeq FAQ](https://www.ncbi.nlm.nih.gov/books/NBK50679/)
#' - ftp://ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz
#'
#' @section UCSC:
#'
#' Example URLs:
#'
#' - UCSC *Homo sapiens* hg38 GTF files:
#'   [hg38.ensGene.gtf.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz),
#'   [hg38.knownGene.gtf.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz),
#'   [hg38.ncbiRefSeq.gtf.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz),
#'   [hg38.refGene.gtf.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz)
#'
#' Related URLs:
#'
#' - [UCSC downloads](http://hgdownload.soe.ucsc.edu/downloads.html)
#' - [UCSC hg38 bigZips FTP](ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
#' - [UCSC hgTables](http://genome.ucsc.edu/cgi-bin/hgTables)
#'
#' @section FlyBase:
#'
#' Example URLs:
#'
#' - FlyBase *Drosophila melanogaster* r6.24
#'   [GTF](ftp://ftp.flybase.net/releases/FB2020_06/dmel_r6.37/gtf/dmel-all-r6.37.gtf.gz),
#'   [GFF3](ftp://ftp.flybase.net/releases/FB2020_06/dmel_r6.37/gff/dmel-all-r6.37.gff.gz)
#'
#' @section WormBase:
#'
#' Example URLs:
#'
#' - WormBase *Caenorhabditis elegans* WS267
#'   [GTF](ftp://ftp.wormbase.org/pub/wormbase/releases/WS279/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz),
#'   [GFF3](ftp://ftp.wormbase.org/pub/wormbase/releases/WS279/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS279.annotations.gff3.gz)
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
#' @note Updated 2021-01-22.
#' @noRd
.makeGRangesFromGFF <- function(
    file,
    level = c("genes", "transcripts"),
    organism = NULL,
    genomeBuild = NULL,
    release = NULL,
    ignoreVersion = TRUE,
    synonyms = FALSE,
    seqinfo = NULL,
    ## Internal-only arguments:
    broadClass = TRUE
) {
    assert(
        isString(file),
        isOrganism(organism, nullOK = TRUE),
        isString(genomeBuild, nullOK = TRUE),
        isString(release, nullOK = TRUE) || isInt(release, nullOK = TRUE),
        isFlag(ignoreVersion),
        isFlag(synonyms),
        is(seqinfo, "Seqinfo") || is.null(seqinfo),
        isFlag(broadClass)
    )
    level <- match.arg(level)
    alert(sprintf(
        fmt = "Making {.var GRanges} from GFF file ({.file %s}).",
        basename(file)
    ))
    txdb <- makeTxDbFromGFF(file = tmpfile, seqinfo = seqinfo)
    gffMeta <- attr(x = txdb, which = "gffMetadata", exact = TRUE)



    ## FIXME SAFE TO TAKE THIS OUT AND JUST CALL DIRECTLY BELOW IN
    ##       RTRACKLAYER IMPORT STEP INSTEAD?
    rawRanges <- import(tmpfile)

    ## FIXME RETHINK THIS....
    ## FIXME PARSE THE GFF FILE DIRECTLY AND GET THE METADATA THAT WAY\
    ##       INSTEAD. REWORK THE INTERNAL TXDB FUNCTION TO ALSO DETECT IF
    ##       GTF OR GFF. CAN USE FILE NAME....
    ## FIXME CALL THESE IN THE RTRACKLAYER IMPORT STEP INSTEAD....
    source <- .grangesSource(rawRanges)
    format <- .grangesFormat(rawRanges)
    assert(isString(source), isString(format))
    if (isTRUE(synonyms) && !isSubset(source, c("Ensembl", "GENCODE"))) {
        ## nocov start
        stop(sprintf(
            "Synonyms only supported for genomes from: %s.",
            toString(c("Ensembl", "GENCODE"))
        ))
        ## nocov end
    }
    meta <- list(
        "call" = match.call(),
        "date" = Sys.Date(),
        "file" = file,
        "format" = format,
        "level" = level,
        "md5" = .md5(file = tmpfile),
        "sha256" = .sha256(file = tmpfile),
        "source" = source
    )



    ## FIXME HAND OFF TO ENSEMBLDB FOR ENSEMBL FILES...




    gr1 <- .makeGRangesFromTxDb(object = txdb, level = level)
    metadata(gr1) <- meta
    gr2 <- .makeGRangesFromRtracklayer(
        object = rawRanges,
        level = level,
        format = format,
        source = source
    )
    metadata(gr2) <- meta

    ## FIXME SPLIT THIS OUT TO SEPARATE FUNCTION.
    ## FIXME NEED TO CHECK THAT ALL IDS DEFINED IN LEFT SIDE ARE IN RIGHT SIDE.
    ##       USE ISSUBSET FOR THIS...

    idCol1 <- .matchGRangesNamesColumn(gr1)
    idCol2 <- .matchGRangesNamesColumn(gr2)
    assert(identical(idCol1, idCol2))
    idCol <- idCol1
    mcols1 <- mcols(gr1)
    mcols2 <- mcols(gr2)
    extra <- setdiff(colnames(mcols2), colnames(mcols1))
    mcols <- leftJoin(x = mcols1, y = mcols2[c(idCol, extra)], by = idCol)
    assert(identical(mcols[[idCol]], mcols1[[idCol]]))
    gr <- gr1
    mcols(gr) <- mcols
    ## FIXME Need to handle genomeBuild, organism, release for Ensembl here.
    ## FIXME Need to call detectOrganism here internally to slot?
    .makeGRanges(
        object = gr,
        ignoreVersion = ignoreVersion,
        broadClass = broadClass,
        synonyms = synonyms
    )
}



#' @describeIn makeGRangesFromGFF Primary function.
#' @export
makeGRangesFromGFF <- function(
    file,
    level = c("genes", "transcripts"),
    organism = NULL,
    genomeBuild = NULL,
    release = NULL,
    ignoreVersion = TRUE,
    synonyms = FALSE,
    seqinfo = NULL
) {
    .makeGRangesFromGFF(
        file = file,
        level = match.arg(level),
        organism = organism,
        genomeBuild = genomeBuild,
        release = release,
        ignoreVersion = ignoreVersion,
        synonyms = synonyms,
        seqinfo = seqinfo
    )
}
