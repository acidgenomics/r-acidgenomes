## NOTE Can consider using ensembldb to parse GFF files once v2.15.2 is
##      available on Bioconductor. Current stable version has issues parsing
##      compressed GFF3 files.
##      See issue:
##      - https://github.com/jorainer/ensembldb/issues/114



## nolint start

#' Make GenomicRanges from a GFF/GTF file
#'
#' @export
#' @note Updated 2021-08-06.
#'
#' @details
#' Remote URLs and compressed files are supported.
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
#'   [hg38.ensGene.gtf.gz](ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz),
#'   [hg38.knownGene.gtf.gz](ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz),
#'   [hg38.ncbiRefSeq.gtf.gz](ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz),
#'   [hg38.refGene.gtf.gz](ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz)
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
#' @return `GenomicRanges`.
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
#' ## Some examples here are commented because they are CPU-intensive and
#' ## can cause CI timeouts.
#'
#' ## Ensembl ====
#' file <- pasteURL(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-102",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.102.gtf.gz",
#'     protocol = "ftp"
#' )
#' genes <- makeGRangesFromGFF(
#'     file = file,
#'     level = "genes",
#'     ignoreVersion = FALSE
#' )
#' summary(genes)
#' ## > transcripts <- makeGRangesFromGFF(
#' ## >     file = file,
#' ## >     level = "transcripts",
#' ## >     ignoreVersion = FALSE
#' ## > )
#' ## > summary(transcripts)
#'
#' ## GENCODE ====
#' ## > file <- pasteURL(
#' ## >     "ftp.ebi.ac.uk",
#' ## >     "pub",
#' ## >     "databases",
#' ## >     "gencode",
#' ## >     "Gencode_human",
#' ## >     "release_36",
#' ## >     "gencode.v36.annotation.gtf.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > genes <- makeGRangesFromGFF(file = file, level = "genes")
#' ## > summary(genes)
#' ## > transcripts <- makeGRangesFromGFF(file = file, level = "transcripts")
#' ## > summary(transcripts)
#'
#' ## RefSeq ====
#' ## > file <- pasteURL(
#' ## >     "ftp.ncbi.nlm.nih.gov",
#' ## >     "genomes",
#' ## >     "refseq",
#' ## >     "vertebrate_mammalian",
#' ## >     "Homo_sapiens",
#' ## >     "all_assembly_versions",
#' ## >     "GCF_000001405.39_GRCh38.p13",
#' ## >     "GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > genes <- makeGRangesFromGFF(file = file, level = "genes")
#' ## > summary(genes)
#' ## > transcripts <- makeGRangesFromGFF(file = file, level = "transcripts")
#' ## > summary(transcripts)
#'
#' ## UCSC ====
#' ## > file <- pasteURL(
#' ## >     "hgdownload.soe.ucsc.edu",
#' ## >     "goldenPath",
#' ## >     "hg38",
#' ## >     "bigZips",
#' ## >     "genes",
#' ## >     "hg38.ensGene.gtf.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > genes <- makeGRangesFromGFF(file = file, level = "genes")
#' ## > summary(genes)
#' ## > transcripts <- makeGRangesFromGFF(file = file, level = "transcripts")
#' ## > summary(transcripts)

## nolint end

makeGRangesFromGFF <- function(
    file,
    level = c("genes", "transcripts"),
    ignoreVersion = TRUE,
    synonyms = FALSE
) {
    assert(
        .isSupportedGFF(file),
        isFlag(ignoreVersion),
        isFlag(synonyms)
    )
    level <- match.arg(level)
    alert(sprintf(
        fmt = "Making {.cls %s} from GFF file ({.file %s}).",
        "GenomicRanges", basename(file)
    ))
    if (
        isMatchingRegex(
            pattern = .gffPatterns[["ucsc"]],
            x = basename(file)
        )
    ) {
        alertInfo("UCSC genome annotation file detected.")
        txdb <- makeTxDbFromGFF(file)
        gr <- makeGRangesFromTxDb(
            object = txdb,
            level = level,
            ignoreVersion = ignoreVersion,
            synonyms = synonyms
        )
    } else {
        gr <- .makeGRangesFromRtracklayer(
            file = file,
            level = level,
            ignoreVersion = ignoreVersion,
            synonyms = synonyms
        )
    }
    metadata(gr)[["call"]] <- tryCatch(
        expr = standardizeCall(),
        error = function(e) NULL
    )
    gr
}
