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

data(
    GRanges,
    package = "AcidTest",
    envir = environment()
)
gr <- GRanges

## nolint end

n <- list(
    "celegans" = list(
        "wormbase" = c(
            "exons" = NA, # FIXME
            "genes" = NA, # FIXME
            "transcripts" = 32004L
        )
    ),
    "dmelanogaster" = list(
        "flybase" = c(
            "exons" = NA, # FIXME
            "genes" = 17873L,
            "transcripts" = 35703L
        )
    ),
    "hsapiens" = list(
        "ensembl" = c(
            ## GRCh38 release 110.
            ## Ensembl preserves "TEC" and "artifact" biotypes in FASTA and
            ## GFF files, which differs from GENCODE.
            "exons" = 733353L,
            "genes" = 68974L,
            "transcripts" = 275741L
        ),
        "gencode" = c(
            ## GRCh38 release 44.
            ## GENCODE removes "TEC" and "artifact" biotypes in FASTA and GFF
            ## files, which differs from Ensembl.
            "exons" = NA, # FIXME
            "genes" = NA, # FIXME
            "transcripts" = 252835L
        )
    )
)

## Updated 2023-12-21.
txFastas <- c(
    "ensembl_cdna" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-110",
        "fasta",
        "homo_sapiens",
        "cdna",
        "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        protocol = "ftp"
    ),
    "ensembl_ncrna" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-110",
        "fasta",
        "homo_sapiens",
        "ncrna",
        "Homo_sapiens.GRCh38.ncrna.fa.gz",
        protocol = "ftp"
    ),
    "flybase" = pasteUrl(
        "ftp.flybase.net",
        "releases",
        "FB2023_06",
        "dmel_r6.55",
        "fasta",
        "dmel-all-transcript-r6.55.fasta.gz",
        protocol = "ftp"
    ),
    "gencode" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_44",
        "gencode.v44.transcripts.fa.gz",
        protocol = "ftp"
    ),
    "wormbase" = pasteUrl(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS290",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS290.mRNA_transcripts.fa.gz",
        protocol = "ftp"
    )
)

## Updated 2023-12-21.
gffs <- c(
    "ensembl_grch37_gff3" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-110",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch37_gff3_scaff" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-110",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch37_gtf" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-110",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.gtf.gz",
        protocol = "ftp"
    ),
    "ensembl_grch37_gtf_scaff" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "grch37",
        "release-110",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gtf.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gff3" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-110",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.110.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gff3_scaff" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-110",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.110.chr_patch_hapl_scaff.gff3.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gtf" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-110",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.110.gtf.gz",
        protocol = "ftp"
    ),
    "ensembl_grch38_gtf_scaff" = pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-110",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.110.chr_patch_hapl_scaff.gtf.gz",
        protocol = "ftp"
    ),
    ## NOTE Not supported yet in the package.
    "flybase_gff3" = pasteUrl(
        "ftp.flybase.net",
        "releases",
        "FB2023_06",
        "dmel_r6.55",
        "gff",
        "dmel-all-r6.55.gff.gz",
        protocol = "ftp"
    ),
    "flybase_gtf" = pasteUrl(
        "ftp.flybase.net",
        "releases",
        "FB2023_06",
        "dmel_r6.55",
        "gtf",
        "dmel-all-r6.55.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grch37_gff3" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_44",
        "GRCh37_mapping",
        "gencode.v44lift37.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grch37_gtf" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_44",
        "GRCh37_mapping",
        "gencode.v44lift37.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grch38_gff3" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_44",
        "gencode.v44.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grch38_gtf" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_44",
        "gencode.v44.annotation.gtf.gz",
        protocol = "ftp"
    ),
    "gencode_grcm38_gff3" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M33",
        "gencode.vM33.annotation.gff3.gz",
        protocol = "ftp"
    ),
    "gencode_grcm38_gtf" = pasteUrl(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_mouse",
        "release_M33",
        "gencode.vM33.annotation.gtf.gz",
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
    ## FIXME Add coverage of this in longtests.
    "refseq_grch38_pipeline_gff3" = pasteUrl(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "all",
        "GCA", "000", "001", "405",
        "GCA_000001405.15_GRCh38",
        "seqs_for_alignment_pipelines.ucsc_ids",
        "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz",
        protocol = "ftp"
    ),
    ## FIXME Add coverage of this in longtests.
    "refseq_grch38_pipeline_gtf" = pasteUrl(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "all",
        "GCA", "000", "001", "405",
        "GCA_000001405.15_GRCh38",
        "seqs_for_alignment_pipelines.ucsc_ids",
        "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz",
        protocol = "ftp"
    ),
    ## FIXME Add coverage of this in longtests.
    "refseq_t2t_gff3" = pasteUrl(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "all_assembly_versions",
        "GCF_009914755.1_T2T-CHM13v2.0",
        "GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
        protocol = "ftp"
    ),
    ## FIXME Add coverage of this in longtests.
    "refseq_t2t_gtf" = pasteUrl(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "all_assembly_versions",
        "GCF_009914755.1_T2T-CHM13v2.0",
        "GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz",
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
        "WS290",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS290.annotations.gff3.gz",
        protocol = "ftp"
    ),
    "wormbase_gtf" = pasteUrl(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS290",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS290.canonical_geneset.gtf.gz",
        protocol = "ftp"
    )
)
