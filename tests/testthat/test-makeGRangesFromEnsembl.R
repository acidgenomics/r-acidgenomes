test_that("Genes", {
    expect_false(any(biocPkgs %in% (.packages())))
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        release = 108L,
        ignoreVersion = TRUE
    )
    expect_false(any(biocPkgs %in% (.packages())))
    expect_s4_class(object, "EnsemblGenes")
    expect_identical(
        object = metadata(object)[c(
            "annotationHubId",
            "genomeBuild",
            "organism",
            "release"
        )],
        expected = list(
            "annotationHubId" = "AH109336",
            "genomeBuild" = "GRCh38",
            "organism" = "Homo sapiens",
            "release" = 108L
        )
    )
    expect_length(object, 70616L)
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
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "ncbiGeneId" = "CompressedIntegerList",
            "seqCoordSystem" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENSG00000116044"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "2",
            "start" = "177218667",
            "end" = "177392756",
            "width" = "174090",
            "strand" = "-",
            "broadClass" = "coding",
            "canonicalTranscript" = "ENST00000397062",
            "description" = paste(
                "NFE2 like bZIP transcription factor 2",
                "[Source:HGNC Symbol;Acc:HGNC:7782]"
            ),
            "geneBiotype" = "protein_coding",
            "geneId" = "ENSG00000116044",
            "geneIdVersion" = "ENSG00000116044.17",
            "geneName" = "NFE2L2",
            "ncbiGeneId" = "4780",
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

test_that("Transcripts", {
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        release = 108L,
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, 275721L)
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
            "gcContent" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "geneSeqEnd" = "Rle",
            "geneSeqStart" = "Rle",
            "ncbiGeneId" = "CompressedIntegerList",
            "seqCoordSystem" = "Rle",
            "txBiotype" = "Rle",
            "txCdsSeqEnd" = "Rle",
            "txCdsSeqStart" = "Rle",
            "txExternalName" = "Rle",
            "txId" = "Rle",
            "txIdNoVersion" = "Rle",
            "txIdVersion" = "Rle",
            "txIsCanonical" = "Rle",
            "txName" = "Rle",
            "txSupportLevel" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000397062.8"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "2",
            "start" = "177230308",
            "end" = "177264727",
            "width" = "34420",
            "strand" = "-",
            "broadClass" = "coding",
            "canonicalTranscript" = "ENST00000397062",
            "description" = paste(
                "NFE2 like bZIP transcription factor 2",
                "[Source:HGNC Symbol;Acc:HGNC:7782]"
            ),
            "gcContent" = "42.1095666394113",
            "geneBiotype" = "protein_coding",
            "geneId" = "ENSG00000116044.17",
            "geneIdNoVersion" = "ENSG00000116044",
            "geneIdVersion" = "ENSG00000116044.17",
            "geneName" = "NFE2L2",
            "geneSeqEnd" = "177392756",
            "geneSeqStart" = "177218667",
            "ncbiGeneId" = "4780",
            "seqCoordSystem" = "chromosome",
            "txBiotype" = "protein_coding",
            "txCdsSeqEnd" = "177264576",
            "txCdsSeqStart" = "177230785",
            "txExternalName" = "NFE2L2-201",
            "txId" = "ENST00000397062.8",
            "txIdNoVersion" = "ENST00000397062",
            "txIdVersion" = "ENST00000397062.8",
            "txIsCanonical" = "1",
            "txName" = "ENST00000397062",
            "txSupportLevel" = "1"
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
    expect_s4_class(object, "EnsemblGenes")
    expect_length(object, 64102L)
    expect_identical(head(names(object), 1L), "ENSG00000228572")
    ## Transcripts
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        genomeBuild = "GRCh37",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, 215647L)
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
