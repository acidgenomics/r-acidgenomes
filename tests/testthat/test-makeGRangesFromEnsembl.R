context("makeGRangesFromEnsembl")

skip_if_not(hasInternet())

Rle <- structure("Rle", package = "S4Vectors")  # nolint

test_that("Genes", {
    object <- makeGRangesFromEnsembl(
        organism = organism,
        level = "genes",
        release = release
    )
    expect_s4_class(object, "GRanges")
    expect_identical(
        object = metadata(object)[
            c("ensemblRelease", "genomeBuild", "id", "level", "organism")],
        expected = list(
            ensemblRelease = 97L,
            genomeBuild = "GRCh38",
            id = "AH73881",
            level = "genes",
            organism = "Homo sapiens"
        )
    )
    expect_identical(length(object), 67667L)
    expect_identical(
        object = head(names(object), 3L),
        expected = c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419")
    )
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            broadClass = Rle,
            description = Rle,
            entrezID = "list",
            geneBiotype = Rle,
            geneID = Rle,
            geneIDVersion = Rle,
            geneName = Rle,
            seqCoordSystem = Rle
        )
    )
})

## Transcript verion metadata isn't saved in older EnsDb releases, such as v87.
test_that("Transcripts", {
    object <- makeGRangesFromEnsembl(
        organism = organism,
        level = "transcripts",
        release = release,
        ignoreTxVersion = FALSE
    )
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 248916L)
    expect_identical(
        object = head(names(object), n = 2L),
        expected = c("ENST00000000233.10", "ENST00000000412.8")
    )
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            broadClass = Rle,
            description = Rle,
            entrezID = "list",
            geneBiotype = Rle,
            geneID = Rle,
            geneIDVersion = Rle,
            geneName = Rle,
            seqCoordSystem = Rle,
            transcriptBiotype = Rle,
            transcriptCdsSeqEnd = Rle,
            transcriptCdsSeqStart = Rle,
            transcriptID = Rle,
            transcriptIDVersion = Rle,
            transcriptName = Rle,
            transcriptSupportLevel = Rle
        )
    )
})

test_that("GRCh37 (hg19)", {
    ## Conditionally test if optional EnsDb.Hsapiens.v75 package is installed.
    skip_if_not("EnsDb.Hsapiens.v75" %in% rownames(installed.packages()))
    ## Genes
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        genomeBuild = "GRCh37"
    )
    expect_is(object, "GRanges")
    expect_identical(length(object), 64102L)
    expect_identical(head(names(object), 1L), "ENSG00000000003")
    ## Transcripts
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        genomeBuild = "GRCh37"
    )
    expect_is(object, "GRanges")
    expect_identical(length(object), 215647L)
    expect_identical(head(names(object), 1L), "ENST00000000233")
})

test_that("UCSC identifier matching (hg38)", {
    x <- makeGRangesFromEnsembl(organism = "Homo sapiens", genomeBuild = "hg38")
    expect_s4_class(x, "GRanges")
})

test_that("Organism with 3 words", {
    x <- makeGRangesFromEnsembl(organism = "Canis lupus familiaris")
    expect_s4_class(x, "GRanges")
})

test_that("Invalid parameters", {
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "Homo sapiens",
            release = 86L
        ),
        regexp = ">= 87"
    )
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "AAA",
            genomeBuild = "BBB"
        ),
        regexp = "No ID matched on AnnotationHub"
    )
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = c("Homo sapiens", "Mus musculus")
        ),
        regexp = "isString"
    )
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "Homo sapiens",
            level = "XXX"
        ),
        regexp = "'arg' should be one of \"genes\", \"transcripts\""
    )
})
