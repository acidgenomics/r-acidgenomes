## Exons are covered in longtests.

test_that("Homo sapiens : genes", {
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "genes",
        release = 110L,
        ignoreVersion = FALSE
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
            "annotationHubId" = "AH113665",
            "genomeBuild" = "GRCh38",
            "organism" = "Homo sapiens",
            "release" = 110L
        )
    )
    expect_length(
        object = object,
        n = n[["hsapiens"]][["ensembl"]][["genes"]]
    )
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["geneId"]])
    )
    expect_true(allAreMatchingRegex(
        x = names(object),
        pattern = "^ENSG[0-9]{11}.[0-9]+$"
    ))
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "canonicalTranscript" = "character",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSynonyms" = "CompressedCharacterList",
            "ncbiGeneId" = "CompressedIntegerList",
            "seqCoordSystem" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENSG00000116044.17"]), # nolint
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
            "geneId" = "ENSG00000116044.17",
            "geneIdNoVersion" = "ENSG00000116044",
            "geneIdVersion" = "ENSG00000116044.17",
            "geneName" = "NFE2L2",
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
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

test_that("Homo sapiens : transcripts", {
    object <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        level = "transcripts",
        release = 110L,
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(
        object = object,
        n = n[["hsapiens"]][["ensembl"]][["transcripts"]]
    )
    ## FIXME Check that genes match here.
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["txId"]])
    )
    expect_true(allAreMatchingRegex(
        x = names(object),
        pattern = "^ENST[0-9]{11}.[0-9]+$"
    ))
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "canonicalTranscript" = "character",
            "description" = "character",
            "gcContent" = "numeric",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSeqEnd" = "integer",
            "geneSeqStart" = "integer",
            "geneSynonyms" = "CompressedCharacterList",
            "ncbiGeneId" = "CompressedIntegerList",
            "seqCoordSystem" = "factor",
            "txBiotype" = "factor",
            "txCdsSeqEnd" = "integer",
            "txCdsSeqStart" = "integer",
            "txExternalName" = "character",
            "txId" = "character",
            "txIdNoVersion" = "character",
            "txIdVersion" = "character",
            "txIsCanonical" = "logical",
            "txName" = "character",
            "txSupportLevel" = "factor"
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
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
            "ncbiGeneId" = "4780",
            "seqCoordSystem" = "chromosome",
            "txBiotype" = "protein_coding",
            "txCdsSeqEnd" = "177264576",
            "txCdsSeqStart" = "177230785",
            "txExternalName" = "NFE2L2-201",
            "txId" = "ENST00000397062.8",
            "txIdNoVersion" = "ENST00000397062",
            "txIdVersion" = "ENST00000397062.8",
            "txIsCanonical" = "TRUE",
            "txName" = "ENST00000397062",
            "txSupportLevel" = "1"
        )
    )
})
