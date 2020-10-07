## Updated 2020-07-26.

localOrRemoteFile <- pipette::localOrRemoteFile
pasteURL <- acidbase::pasteURL

## Enable this for even more verbose checks.
checkAgainstTxDb <- TRUE

levels <- c("genes", "transcripts")

files <- list(
    ensemblGFF3 = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "release-100",
        "gff3",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.100.gff3.gz",
        protocol = "ftp"
    ),
    ensemblGTF = pasteURL(
        "ftp.ensembl.org",
        "pub",
        "release-100",
        "gtf",
        "homo_sapiens",
        "Homo_sapiens.GRCh38.100.gtf.gz",
        protocol = "ftp"
    ),
    flybaseGTF = pasteURL(
        "ftp.flybase.net",
        "releases",
        "FB2020_03",
        "dmel_r6.34",
        "gtf",
        "dmel-all-r6.34.gtf.gz",
        protocol = "ftp"
    ),
    gencodeGFF3 = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_34",
        "gencode.v34.annotation.gff3.gz",
        protocol = "ftp"
    ),
    gencodeGTF = pasteURL(
        "ftp.ebi.ac.uk",
        "pub",
        "databases",
        "gencode",
        "Gencode_human",
        "release_34",
        "gencode.v34.annotation.gtf.gz",
        protocol = "ftp"
    ),
    refseqGFF3 = pasteURL(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "reference",
        "GCF_000001405.39_GRCh38.p13",
        "GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
        protocol = "ftp"
    ),
    refseqGTF = pasteURL(
        "ftp.ncbi.nlm.nih.gov",
        "genomes",
        "refseq",
        "vertebrate_mammalian",
        "Homo_sapiens",
        "reference",
        "GCF_000001405.39_GRCh38.p13",
        "GCF_000001405.39_GRCh38.p13_genomic.gtf.gz",
        protocol = "ftp"
    ),
    wormbaseGTF = pasteURL(
        "ftp.wormbase.org",
        "pub",
        "wormbase",
        "releases",
        "WS277",
        "species",
        "c_elegans",
        "PRJNA13758",
        "c_elegans.PRJNA13758.WS277.canonical_geneset.gtf.gz",
        protocol = "ftp"
    )
)



## Ensembl =====================================================================
context("Ensembl")

skip_if_not(hasInternet(url = "ftp://ftp.ensembl.org/"))
length <- c(genes = 60683L, transcripts = 227954L)
mapply(
    testName = c(
        "GTF",
        "GFF3"
    ),
    file = c(
        localOrRemoteFile(files[["ensemblGTF"]]),
        localOrRemoteFile(files[["ensemblGFF3"]])
    ),
    FUN = function(testName, file) {
        test_that(testName, {
            mapply(
                level = levels,
                length = length,
                FUN = function(level, length) {
                    object <- makeGRangesFromGFF(
                        file = file,
                        level = level,
                        .checkAgainstTxDb = checkAgainstTxDb
                    )
                    expect_s4_class(object, "GRanges")
                    expect_identical(length(object), length)
                },
                SIMPLIFY = FALSE
            )
        })
    },
    SIMPLIFY = FALSE
)



## FlyBase =====================================================================
context("FlyBase")

## Expecting warnings about broad class and dropped transcripts.

skip_if_not(hasInternet(url = "ftp://ftp.flybase.net/"))
length <- c(genes = 17875L, transcripts = 35647L)
file <- localOrRemoteFile(files[["flybaseGTF"]])
test_that("GTF", {
    mapply(
        level = levels,
        length = length,
        FUN = function(level, length) {
            suppressWarnings({
                object <- makeGRangesFromGFF(
                    file = file,
                    level = level,
                    .checkAgainstTxDb = checkAgainstTxDb
                )
            })
            expect_s4_class(object, "GRanges")
            expect_identical(length(object), length)
        },
        SIMPLIFY = FALSE
    )
})



## GENCODE =====================================================================
context("GENCODE")

## Expecting warning about PAR genes mismatch in strict mode.
## Note that we're dropping PAR Y chromosome dupes on genes and transcripts.
## 45 PAR dupe genes; 161 PAR dupe transcripts.

skip_if_not(hasInternet(url = "ftp://ftp.ebi.ac.uk/"))
mapply(
    testName = c(
        "GTF",
        "GFF3"
    ),
    file = c(
        localOrRemoteFile(files[["gencodeGTF"]]),
        localOrRemoteFile(files[["gencodeGFF3"]])
    ),
    length = list(
        c(genes = 60669L, transcripts = 228048L),
        c(genes = 60624L, transcripts = 227887L)
    ),
    FUN = function(testName, file, length) {
        test_that(testName, {
            mapply(
                level = levels,
                length = length,
                FUN = function(level, length) {
                    object <- makeGRangesFromGFF(
                        file = file,
                        level = level,
                        .checkAgainstTxDb = checkAgainstTxDb
                    )
                    expect_s4_class(object, "GRanges")
                    expect_identical(length(object), length)
                },
                SIMPLIFY = FALSE
            )
        })
    },
    SIMPLIFY = FALSE
)



## RefSeq ======================================================================
context("RefSeq")

# GenomicFeatures currently chokes on RefSeq files.

skip_if_not(hasInternet(url = "ftp://ftp.ncbi.nlm.nih.gov/"))
mapply(
    testName = c(
        "GTF",
        "GFF3"
    ),
    file = c(
        localOrRemoteFile(files[["refseqGTF"]]),
        localOrRemoteFile(files[["refseqGFF3"]])
    ),
    length = list(
        c(genes = 61197L, transcripts = 174165L),
        c(genes = 54606L, transcripts = 160987L)
    ),
    FUN = function(testName, file, class, length) {
        test_that(testName, {
            mapply(
                level = levels,
                length = length,
                FUN = function(level, length) {
                    object <- makeGRangesFromGFF(
                        file = file,
                        level = level,
                        .checkAgainstTxDb = FALSE
                    )
                    ## > expect_s4_class(object, class)
                    expect_identical(length(object), length)
                },
                SIMPLIFY = FALSE
            )
        })
    },
    SIMPLIFY = FALSE
)



## WormBase ====================================================================
context("WormBase")

skip_if_not(hasInternet(url = "ftp://ftp.wormbase.org/"))
length <- c(46932L, 66847L)
file <- localOrRemoteFile(files[["wormbaseGTF"]])
test_that("GTF", {
    mapply(
        level = levels,
        length = length,
        FUN = function(level, length) {
            object <- makeGRangesFromGFF(
                file = file,
                level = level,
                .checkAgainstTxDb = checkAgainstTxDb
            )
            expect_s4_class(object, "GRanges")
            expect_identical(length(object), length)
        },
        SIMPLIFY = FALSE
    )
})
