context("makeGRangesFromEnsembl")

skip_if_not(hasInternet())

test_that("Genes", {
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        release = 102L,
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_identical(
        object = metadata(object)[c(
            "annotationHubId",
            "genomeBuild",
            "organism",
            "release"
        )],
        expected = list(
            "annotationHubId" = "AH89180",
            "genomeBuild" = "GRCh38",
            "organism" = "Homo sapiens",
            "release" = 102L
        )
    )
    expect_identical(length(object), 68001L)
    expect_identical(
        object = head(names(object), 3L),
        expected = c(
            "ENSG00000228572",
            "ENSG00000182378",
            "ENSG00000226179"
        )
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "canonicalTranscript" = "Rle",
            "description" = "Rle",
            "entrezId" = "CompressedIntegerList",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "seqCoordSystem" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENSG00000223972"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "1",
            "start" = "11869",
            "end" = "14409",
            "width" = "2541",
            "strand" = "+",
            "broadClass" = "pseudo",
            "canonicalTranscript" = "ENST00000456328",
            "description" = paste(
                "DEAD/H-box helicase 11 like 1 (pseudogene)",
                "[Source:HGNC Symbol;Acc:HGNC:37102]"
            ),
            "entrezId" = "c(84771, 727856, 100287102, 100287596, 102725121)",
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "seqCoordSystem" = "chromosome"
        )
    )
    expect_true(isSubset(
        x = c(
            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
            "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"
        ),
        y = levels(seqnames(object))
    ))
    expect_identical(
        object = as.data.frame(seqinfo(object))["1", , drop = TRUE],
        expected = list(
            "seqlengths" = 248956422L,
            "isCircular" = FALSE,
            "genome" = "GRCh38"
        )
    )
})

## Transcript verion metadata isn't saved in older EnsDb releases, such as v87.
test_that("Transcripts", {
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        release = 102L,
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_identical(length(object), 254853L)
    expect_identical(
        object = head(names(object), n = 2L),
        expected = c("ENST00000635602.1", "ENST00000635506.1")
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "canonicalTranscript" = "Rle",
            "description" = "Rle",
            "entrezId" = "CompressedIntegerList",
            "gcContent" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "geneSeqEnd" = "Rle",
            "geneSeqStart" = "Rle",
            "seqCoordSystem" = "Rle",
            "txBiotype" = "Rle",
            "txCdsSeqEnd" = "Rle",
            "txCdsSeqStart" = "Rle",
            "txId" = "Rle",
            "txIdNoVersion" = "Rle",
            "txIdVersion" = "Rle",
            "txName" = "Rle",
            "txSupportLevel" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000635602.1"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "7",
            "start" = "12704",
            "end" = "27199",
            "width" = "14496",
            "strand" = "+",
            "broadClass" = "noncoding",
            "canonicalTranscript" = "ENST00000635602",
            "description" = "novel transcript",
            "entrezId" = "NA",
            "gcContent" = "42.090395480226",
            "geneBiotype" = "lncRNA",
            "geneId" = "ENSG00000283061.1",
            "geneIdNoVersion" = "ENSG00000283061",
            "geneIdVersion" = "ENSG00000283061.1",
            "geneName" = "AC215522.3",
            "geneSeqEnd" = "27234",
            "geneSeqStart" = "12704",
            "seqCoordSystem" = "chromosome",
            "txBiotype" = "lncRNA",
            "txCdsSeqEnd" = NA_character_,
            "txCdsSeqStart" = NA_character_,
            "txId" = "ENST00000635602.1",
            "txIdNoVersion" = "ENST00000635602",
            "txIdVersion" = "ENST00000635602.1",
            "txName" = "ENST00000635602",
            "txSupportLevel" = "5"
        )
    )
})

test_that("GRCh37", {
    ## Conditionally test if optional EnsDb.Hsapiens.v75 package is installed.
    skip_if_not_installed("EnsDb.Hsapiens.v75")
    ## Genes
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        genomeBuild = "GRCh37",
        ignoreVersion = TRUE
    )
    expect_is(object, "EnsemblGenes")
    expect_identical(length(object), 64102L)
    expect_identical(head(names(object), 1L), "ENSG00000228572")
    ## Transcripts
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        genomeBuild = "GRCh37",
        ignoreVersion = TRUE
    )
    expect_is(object, "EnsemblTranscripts")
    expect_identical(length(object), 215647L)
    expect_identical(head(names(object), 1L), "ENST00000478759")
})

test_that("Organism with 3 words", {
    x <- makeGRangesFromEnsembl(
        organism = "Canis lupus familiaris",
        level = "genes"
    )
    expect_s4_class(x, "EnsemblGenes")
})

test_that("Invalid parameters", {
    ## Currently only supports releases back to Ensembl 87.
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "Homo sapiens",
            release = 86L
        ),
        regexp = "No entry matched on AnnotationHub"
    )
    expect_error(
        object = makeGRangesFromEnsembl(
            organism = "Homo sapiens",
            genomeBuild = "BBB"
        ),
        regexp = "No entry matched on AnnotationHub"
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
