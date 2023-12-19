## nolint start

DataFrame <- S4Vectors::DataFrame
data <- utils::data
hasInternet <- goalie::hasInternet
pasteUrl <- AcidBase::pasteUrl
seqinfo <- GenomeInfoDb::seqinfo
seqlengths <- GenomeInfoDb::seqlengths
seqnames <- GenomeInfoDb::seqnames
simpleClass <- AcidBase::simpleClass
tempdir2 <- AcidBase::tempdir2
unlink2 <- AcidBase::unlink2

CompressedIntegerList <- structure(
    .Data = "CompressedIntegerList",
    package = "IRanges"
)
## Rle encoding disabled in 0.7.3 update.
## > Rle <- structure(
## >     .Data = "Rle",
## >     package = "S4Vectors"
## > )

data(
    GRanges,
    package = "AcidTest",
    envir = environment()
)
gr <- GRanges

## nolint end

## Updated 2023-09-26.
txFastas <- c(
    "ensembl" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-108",
        "fasta",
        "homo_sapiens",
        "cdna",
        "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        protocol = "ftp"
    ),
    "flybase" = pasteUrl(
        "ftp.flybase.net",
        "releases",
        "FB2022_06",
        "dmel_r6.49",
        "fasta",
        "dmel-all-transcript-r6.49.fasta.gz",
        protocol = "ftp"
    ),
    "gencode" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_42",
        "gencode.v42.transcripts.fa.gz",
        protocol = "ftp"
    ),
    "wormbase" = pasteUrl(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS287",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS287.mRNA_transcripts.fa.gz",
        protocol = "ftp"
    )
)

## Updated 2023-01-30.
gffs <- c(
    "ensembl_grch37_gff3" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-108",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch37_gtf" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-108",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gtf.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gff3" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-108",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.108.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gtf" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-108",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.108.gtf.gz",
        protocol = "ftp"
    ),
    ## NOTE Not supported yet in the package.
    "flybase_gff3" = pasteUrl(
        "ftp.flybase.net",
        "releases",
        "FB2022_06",
        "dmel_r6.49",
        "gff",
        "dmel-all-r6.49.gff.gz",
        protocol = "ftp"
    ),
    "flybase_gtf" = pasteUrl(
        "ftp.flybase.net",
        "releases",
        "FB2022_06",
        "dmel_r6.49",
        "gtf",
        "dmel-all-r6.49.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grch37_gff3" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_42",
        "GRCh37_mapping",
        "gencode.v42lift37.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grch37_gtf" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_42",
        "GRCh37_mapping",
        "gencode.v42lift37.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grch38_gff3" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_42",
        "gencode.v42.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grch38_gtf" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_42",
        "gencode.v42.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grcm38_gff3" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M31",
        "gencode.vM31.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grcm38_gtf" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M31",
        "gencode.vM31.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "refseq_grch38_gff3" = pasteUrl(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "all_assembly_versions",
        "GCF_000001405.40_GRCh38.p14",
        "GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
        protocol = "ftp"
    ),
    "refseq_grch38_gtf" = pasteUrl(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "all_assembly_versions",
        "GCF_000001405.40_GRCh38.p14",
        "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
        protocol = "ftp"
    ),
    "refseq_grch38_pipeline_gff3" = pasteUrl(
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
    "refseq_grch38_pipeline_gtf" = pasteUrl(
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
    "ucsc_hg38_knowngene_gtf" = pasteUrl(
        "hgdownload.soe.ucsc.edu",
        "goldenPath",
        "hg38",
        "bigZips",
        "genes",
        "hg38.knownGene.gtf.gz",
        protocol = "ftp"
    ),
    "ucsc_hg38_ncbirefseq_gtf" = pasteUrl(
        "hgdownload.soe.ucsc.edu",
        "goldenPath",
        "hg38",
        "bigZips",
        "genes",
        "hg38.ncbiRefSeq.gtf.gz",
        protocol = "ftp"
    ),
    ## NOTE Not supported yet in the package.
    "wormbase_gff3" = pasteUrl(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS287",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS287.annotations.gff3.gz",
        protocol = "ftp"
    ),
    "wormbase_gtf" = pasteUrl(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS287",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS287.canonical_geneset.gtf.gz",
        protocol = "ftp"
    )
)
