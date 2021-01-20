data(
    DFrame,
    GRanges,
    RangedSummarizedExperiment,
    SingleCellExperiment,
    SummarizedExperiment_transcripts,
    matrix,
    matrix_lfc,
    sparseMatrix,
    package = "AcidTest",
    envir = environment()
)

df <- DFrame
gr <- GRanges
lfc <- matrix_lfc
mat <- matrix
rse <- RangedSummarizedExperiment
sce <- SingleCellExperiment
sparse <- sparseMatrix
txse <- SummarizedExperiment_transcripts

organism <- "Homo sapiens"
release <- 97L

## nolint start
DataFrame <- S4Vectors::DataFrame
GRanges <- GenomicRanges::GRanges
IRanges <- IRanges::IRanges
SummarizedExperiment <- SummarizedExperiment::SummarizedExperiment
assay <- SummarizedExperiment::assay
`assay<-` <- SummarizedExperiment::`assay<-`
cause <- goalie::cause
hasInternet <- goalie::hasInternet
rowRanges <- SummarizedExperiment::rowRanges
`rowRanges<-` <- SummarizedExperiment::`rowRanges<-`
skip_on_docker <- goalie::skip_on_docker
str_pad <- stringr::str_pad
tibble <- tibble::tibble
## nolint end

gffUrls <- c(
    pasteURL(
        "ftp.ensembl.org",
        "pub",
        "release-102",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.102.gtf.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ensembl.org",
        "pub",
        "release-102",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.102.gff3.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-102",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gtf.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-102",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gff3.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_36",
        "gencode.v36.annotation.gtf.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_36",
        "gencode.v36.annotation.gff3.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_36",
        "GRCh37_mapping",
        "gencode.v36lift37.annotation.gtf.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_36",
        "GRCh37_mapping",
        "gencode.v36lift37.annotation.gff3.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M25",
        "gencode.vM25.annotation.gtf.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M25",
        "gencode.vM25.annotation.gff3.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "all_assembly_versions",
        "GCF_000001405.38_GRCh38.p12",
        "GCF_000001405.38_GRCh38.p12_genomic.gtf.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "all_assembly_versions",
        "GCF_000001405.38_GRCh38.p12",
        "GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
        protocol = "ftp"
    ),
    ## Note that this file doesn't contain any metadata comments.
    pasteURL(
        "ftp.flybase.net",
        "releases",
        "FB2020_06",
        "dmel_r6.37",
        "gtf",
        "dmel-all-r6.37.gtf.gz",
        protocol = "ftp"
    ),
    ## This file is very large and slow to parse.
    pasteURL(
        "ftp.flybase.net",
        "releases",
        "FB2020_06",
        "dmel_r6.37",
        "gff",
        "dmel-all-r6.37.gff.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS279",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz",
        protocol = "ftp"
    ),
    pasteURL(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS279",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS279.annotations.gff3.gz",
        protocol = "ftp"
    )
)


options(acid.test = TRUE)
