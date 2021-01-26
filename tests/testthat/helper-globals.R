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

gffs <- c(
    "ensembl_grch37_gff3" = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-102",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch37_gtf" = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-102",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gtf.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gff3" = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "release-102",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.102.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gtf" = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "release-102",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.102.gtf.gz",
        protocol = "ftp"
    ),
    "flybase_gff3" = pasteURL(
        "ftp.flybase.net",
        "releases",
        "FB2020_06",
        "dmel_r6.37",
        "gff",
        "dmel-all-r6.37.gff.gz",
        protocol = "ftp"
    ),
    "flybase_gtf" = pasteURL(
        "ftp.flybase.net",
        "releases",
        "FB2020_06",
        "dmel_r6.37",
        "gtf",
        "dmel-all-r6.37.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grch37_gff3" = pasteURL(
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
    "gencode_grch37_gtf" = pasteURL(
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
    "gencode_grch38_gff3" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_36",
        "gencode.v36.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grch38_gtf" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_36",
        "gencode.v36.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grcm38_gff3" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M25",
        "gencode.vM25.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grcm38_gtf" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M25",
        "gencode.vM25.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "refseq_grch38_gff3" = pasteURL(
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
    "refseq_grch38_gtf" = pasteURL(
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
    "refseq_grch38_pipeline_gff3" = pasteURL(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "all",
        "GCA",
        "000",
        "001",
        "405",
        "GCA_000001405.15_GRCh38",
        "seqs_for_alignment_pipelines.ucsc_ids",
        "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz",
        protocol = "ftp"
    ),
    "refseq_grch38_pipeline_gtf" = pasteURL(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "all",
        "GCA",
        "000",
        "001",
        "405",
        "GCA_000001405.15_GRCh38",
        "seqs_for_alignment_pipelines.ucsc_ids",
        "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz",
        protocol = "ftp"
    ),
    "ucsc_hg38_ensgene_gtf" = pasteURL(
        "hgdownload.soe.ucsc.edu",
        "goldenPath",
        "hg38",
        "bigZips",
        "genes",
        "hg38.ensGene.gtf.gz",
        protocol = "ftp"
    ),
    "ucsc_hg38_knowngene_gtf" = pasteURL(
        "hgdownload.soe.ucsc.edu",
        "goldenPath",
        "hg38",
        "bigZips",
        "genes",
        "hg38.knownGene.gtf.gz",
        protocol = "ftp"
    ),
    "ucsc_hg38_ncbirefseq_gtf" = pasteURL(
        "hgdownload.soe.ucsc.edu",
        "goldenPath",
        "hg38",
        "bigZips",
        "genes",
        "hg38.ncbiRefSeq.gtf.gz",
        protocol = "ftp"
    ),
    "ucsc_hg38_refgene_gtf" = pasteURL(
        "hgdownload.soe.ucsc.edu",
        "goldenPath",
        "hg38",
        "bigZips",
        "genes",
        "hg38.refGene.gtf.gz",
        protocol = "ftp"
    ),
    "wormbase_gff3" = pasteURL(
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
    ),
    "wormbase_gtf" = pasteURL(
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
    )
)

options(acid.test = TRUE)
