## FIXME Need to debug transcript parsing mismatch with ensembldb:
## > head(setdiff(names(ensdb), names(gff)))
## [1] "ENST00000680009.1" "ENST00000630627.1" "ENST00000630624.1"
## [4] "ENST00000628275.2" "ENST00000628424.1" "ENST00000630311.1"



## nolint start
#' Make genomic ranges (`GRanges`) from a GFF/GTF file
#'
#' @export
#' @note Updated 2023-12-21.
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
#' genomic interval bounded by start and stop codons. CDS refers to the
#' portion of a genomic DNA sequence that is translated, from the start codon
#' to the stop codon.
#' - `exon`: Genomic interval containing 5' UTR (`five_prime_UTR`), CDS, and
#' 3' UTR (`three_prime_UTR`).
#' - `mRNA`: Processed (spliced) mRNA transcript.
#'
#' See also:
#'
#' - [gffutils documentation](https://daler.github.io/gffutils/)
#' - [GenBank GFF documentation](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/)
#' - [stringtie GFF documentation](https://ccb.jhu.edu/software/stringtie/gff.shtml)
#' - [Sequence Ontology GFF3 spec notes](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
#'
#' @section Supported sources:
#'
#' Currently [makeGRangesFromGff()] supports genomes from these sources:
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
#' If you must load a GFF/GTF file directly, then use [makeGRangesFromGff()].
#'
#' Example URLs:
#'
#' - Ensembl *Homo sapiens* GRCh38.p13, release 108:
#' [GTF](https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz),
#' [GFF3](https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/Homo_sapiens.GRCh38.108.gff3.gz)
#' - Ensembl *Homo sapiens* GRCh37, release 108 (87):
#' [GTF](https://ftp.ensembl.org/pub/grch37/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz),
#' [GFF3](https://ftp.ensembl.org/pub/grch37/release-108/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz)
#'
#' @section GENCODE:
#'
#' Example URLs:
#'
#' - GENCODE *Homo sapiens* GRCh38.p13, release 42:
#' [GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz),
#' [GFF3](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gff3.gz)
#' - GENCODE *Homo sapiens* GRCh37, release 42:
#' [GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh37_mapping/gencode.v42lift37.annotation.gtf.gz),
#' [GFF3](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh37_mapping/gencode.v42lift37.annotation.gff3.gz)
#' - GENCODE *Mus musculus* GRCm39, release M31:
#' [GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.annotation.gtf.gz),
#' [GFF3](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.annotation.gff3.gz)
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
#' [current RefSeq spec](https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt)
#' for details.
#'
#' Example URLs:
#'
#' - RefSeq *Homo sapiens* GRCh38.p14
#' [GTF](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz),
#' [GFF3](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz)
#'
#' See also:
#'
#' - [RefSeq FAQ](https://www.ncbi.nlm.nih.gov/books/NBK50679/)
#' - https://ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz
#'
#' @section UCSC:
#'
#' Example URLs:
#'
#' - UCSC *Homo sapiens* hg38 GTF files:
#' [hg38.knownGene.gtf.gz](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz),
#' [hg38.ncbiRefSeq.gtf.gz](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz),
#'
#' Related URLs:
#'
#' - [UCSC downloads](https://hgdownload.soe.ucsc.edu/downloads.html)
#' - [UCSC hg38 bigZips](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/)
#' - [UCSC hgTables](https://genome.ucsc.edu/cgi-bin/hgTables)
#'
#' @section FlyBase:
#'
#' Example URLs:
#'
#' - FlyBase *Drosophila melanogaster* r6.49
#' [GTF](https://ftp.flybase.net/releases/FB2022_06/dmel_r6.49/gtf/dmel-all-r6.49.gtf.gz),
#' [GFF3](https://ftp.flybase.net/releases/FB2022_06/dmel_r6.49/gff/dmel-all-r6.49.gff.gz)
#'
#' @section WormBase:
#'
#' Example URLs:
#'
#' - WormBase *Caenorhabditis elegans* WS287
#' [GTF](ftp://ftp.wormbase.org/pub/wormbase/releases/WS287/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS287.canonical_geneset.gtf.gz),
#' [GFF3](ftp://ftp.wormbase.org/pub/wormbase/releases/WS287/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS287.annotations.gff3.gz)
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
#' `commonName`, `providerVersion`, `provider`, and `releaseDate` accessors.
#'
#' @examples
#' ## Ensembl ====
#' file <- AcidBase::pasteUrl(
#'     "ftp.ensembl.org",
#'     "pub",
#'     "release-108",
#'     "gtf",
#'     "homo_sapiens",
#'     "Homo_sapiens.GRCh38.108.gtf.gz",
#'     protocol = "ftp"
#' )
#' genes <- makeGRangesFromGff(file = file, level = "genes")
#' summary(genes)
#' ## > transcripts <- makeGRangesFromGff(file = file, level = "transcripts")
#' ## > summary(transcripts)
#'
#' ## GENCODE ====
#' ## > file <- AcidBase::pasteUrl(
#' ## >     "ftp.ebi.ac.uk",
#' ## >     "pub",
#' ## >     "databases",
#' ## >     "gencode",
#' ## >     "Gencode_human",
#' ## >     "release_42",
#' ## >     "gencode.v42.annotation.gtf.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > genes <- makeGRangesFromGff(file = file, level = "genes")
#' ## > summary(genes)
#' ## > transcripts <- makeGRangesFromGff(file = file, level = "transcripts")
#' ## > summary(transcripts)
#'
#' ## RefSeq ====
#' ## > file <- AcidBase::pasteUrl(
#' ## >     "ftp.ncbi.nlm.nih.gov",
#' ## >     "genomes",
#' ## >     "refseq",
#' ## >     "vertebrate_mammalian",
#' ## >     "Homo_sapiens",
#' ## >     "all_assembly_versions",
#' ## >     "GCF_000001405.40_GRCh38.p14",
#' ## >     "GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > genes <- makeGRangesFromGff(file = file, level = "genes")
#' ## > summary(genes)
#' ## > transcripts <- makeGRangesFromGff(file = file, level = "transcripts")
#' ## > summary(transcripts)
#'
#' ## UCSC ====
#' ## > file <- AcidBase::pasteUrl(
#' ## >     "hgdownload.soe.ucsc.edu",
#' ## >     "goldenPath",
#' ## >     "hg38",
#' ## >     "bigZips",
#' ## >     "genes",
#' ## >     "hg38.knownGene.gtf.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > genes <- makeGRangesFromGff(file = file, level = "genes")
#' ## > summary(genes)
#' ## > transcripts <- makeGRangesFromGff(file = file, level = "transcripts")
#' ## > summary(transcripts)
## nolint end
makeGRangesFromGff <-
    function(file,
             level = c("genes", "transcripts", "exons"),
             ignoreVersion = FALSE,
             extraMcols = TRUE) {
        assert(
            .isSupportedGff(file),
            isFlag(ignoreVersion),
            isFlag(extraMcols)
        )
        level <- match.arg(level)
        if (isAFile(file)) {
            file <- realpath(file)
        }
        alert(sprintf(
            fmt = "Making {.cls %s} from GFF file ({.file %s}).",
            "GRanges", file
        ))
        tmpfile <- .cacheIt(file)
        meta <- .getGffMetadata(tmpfile)
        if (isAUrl(file)) {
            meta[["url"]] <- file
        }
        meta[["call"]] <- tryCatch(
            expr = standardizeCall(),
            error = function(e) {
                NULL
            }
        )
        if (identical(meta[["provider"]], "UCSC")) {
            alertInfo("UCSC genome annotation file detected.")
            txdb <- .makeTxDbFromGff(file = tmpfile, meta = meta)
            gr <- .makeGRangesFromTxDb(
                object = txdb,
                level = level,
                ignoreVersion = ignoreVersion,
                extraMcols = extraMcols
            )
        } else {
            gr <- .makeGRangesFromRtracklayer(
                file = tmpfile,
                level = level,
                ignoreVersion = ignoreVersion,
                extraMcols = extraMcols,
                meta = meta
            )
        }
        gr
    }
