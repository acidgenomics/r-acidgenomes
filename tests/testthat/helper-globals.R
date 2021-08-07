data(
    GRanges,
    package = "AcidTest",
    envir = environment()
)

gr <- GRanges

## nolint start
CompressedIntegerList <- structure(
    .Data = "CompressedIntegerList",
    package = "IRanges"
)
Rle <- structure(
    .Data = "Rle",
    package = "S4Vectors"
)
DataFrame <- S4Vectors::DataFrame
hasInternet <- goalie::hasInternet
seqinfo <- GenomeInfoDb::seqinfo
seqlengths <- GenomeInfoDb::seqlengths
seqnames <- GenomeInfoDb::seqnames
simpleClass <- AcidBase::simpleClass
## nolint end

## NOTE These aren't currently used in unit tests.
## Updated 2021-08-05.
fastas <- c(
    "ensembl" = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "release-104",
        "fasta",
        "homo_sapiens",
        "cdna",
        "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        protocol = "ftp"
    ),
    "flybase" = pasteURL(
        "ftp.flybase.net",
        "releases",
        "FB2021_03",
        "dmel_r6.40",
        "fasta",
        "dmel-all-transcript-r6.40.fasta.gz",
        protocol = "ftp"
    ),
    "gencode" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_38",
        "gencode.v38.transcripts.fa.gz",
        protocol = "ftp"
    ),
    "wormbase" = pasteURL(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS281",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS281.mRNA_transcripts.fa.gz",
        protocol = "ftp"
    )
)

## Updated 2021-08-05.
gffs <- c(
    "ensembl_grch37_gff3" = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-104",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch37_gtf" = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-104",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gtf.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gff3" = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "release-104",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.104.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gtf" = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "release-104",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.104.gtf.gz",
        protocol = "ftp"
    ),
    ## NOTE Not supported yet in the package.
    "flybase_gff3" = pasteURL(
        "ftp.flybase.net",
        "releases",
        "FB2021_03",
        "dmel_r6.40",
        "gff",
        "dmel-all-r6.40.gff.gz",
        protocol = "ftp"
    ),
    "flybase_gtf" = pasteURL(
        "ftp.flybase.net",
        "releases",
        "FB2021_03",
        "dmel_r6.40",
        "gtf",
        "dmel-all-r6.40.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grch37_gff3" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_38",
        "GRCh37_mapping",
        "gencode.v38lift37.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grch37_gtf" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_38",
        "GRCh37_mapping",
        "gencode.v38lift37.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grch38_gff3" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_38",
        "gencode.v38.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grch38_gtf" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_38",
        "gencode.v38.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grcm38_gff3" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M27",
        "gencode.vM27.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grcm38_gtf" = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M27",
        "gencode.vM27.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "refseq_grch38_gff3" = pasteURL(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "all_assembly_versions",
        "GCF_000001405.39_GRCh38.p13",
        "GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
        protocol = "ftp"
    ),
    ## NOTE TxDb fails on this: "some CDS cannot be mapped to an exon".
    "refseq_grch38_gtf" = pasteURL(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "all_assembly_versions",
        "GCF_000001405.39_GRCh38.p13",
        "GCF_000001405.39_GRCh38.p13_genomic.gtf.gz",
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
    ## NOTE TxDb fails on this: "some CDS cannot be mapped to an exon".
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
    ## NOTE Not supported yet in the package.
    "wormbase_gff3" = pasteURL(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS281",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS281.annotations.gff3.gz",
        protocol = "ftp"
    ),
    "wormbase_gtf" = pasteURL(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS281",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS281.canonical_geneset.gtf.gz",
        protocol = "ftp"
    )
)
