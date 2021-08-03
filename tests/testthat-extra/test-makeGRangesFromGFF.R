context("makeGRangesFromGFF : FlyBase")

skip_if_not(hasInternet())
file <- gffs[["flybase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "FlyBaseGenes")
    expect_identical(length(object), 17875L)
    expect_identical(names(object)[[1L]], "FBgn0031208")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "source" = Rle,
            "type" = Rle
        )
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "FlyBaseTranscripts")
    expect_identical(length(object), 35643L)
    expect_identical(names(object)[[1L]], "FBtr0475186")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "source" = Rle,
            "txId" = Rle,
            "txName" = Rle,
            "type" = Rle
        )
    )
})



context("makeGRangesFromGFF : GENCODE")

skip_if_not(hasInternet())
file <- gffs[["gencode_grch38_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_identical(length(object), 60660L)
    expect_identical(names(object)[[1L]], "ENSG00000223972.5")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "hgncId" = Rle,
            "level" = Rle,
            "source" = Rle,
            "tag" = Rle,
            "type" = Rle
        )
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_identical(length(object), 232117L)
    expect_identical(
        object = names(object)[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["txId"]])[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["txIdNoVersion"]])[[1L]],
        expected = "ENST00000456328"
    )
    expect_identical(
        object = as.character(mcols(object)[["txIdVersion"]])[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneId"]])[[1L]],
        expected = "ENSG00000223972.5"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneIdNoVersion"]])[[1L]],
        expected = "ENSG00000223972"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneIdVersion"]])[[1L]],
        expected = "ENSG00000223972.5"
    )
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "ccdsId" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "havanaTranscript" = Rle,
            "hgncId" = Rle,
            "level" = Rle,
            "ont" = Rle,
            "proteinId" = Rle,
            "source" = Rle,
            "tag" = Rle,
            "txBiotype" = Rle,
            "txId" = Rle,
            "txIdNoVersion" = Rle,
            "txIdVersion" = Rle,
            "txName" = Rle,
            "txSupportLevel" = Rle,
            "type" = Rle
        )
    )
})

file <- gffs[["gencode_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_identical(length(object), 60660L)
    expect_identical(names(object)[[1L]], "ENSG00000223972.5")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "hgncId" = Rle,
            "level" = Rle,
            "ont" = "CompressedCharacterList",
            "source" = Rle,
            "tag" = "CompressedCharacterList",
            "type" = Rle
        )
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_identical(length(object), 232117L)
    expect_identical(
        object = names(object)[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["txId"]])[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["txIdNoVersion"]])[[1L]],
        expected = "ENST00000456328"
    )
    expect_identical(
        object = as.character(mcols(object)[["txIdVersion"]])[[1L]],
        expected = "ENST00000456328.2"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneId"]])[[1L]],
        expected = "ENSG00000223972.5"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneIdNoVersion"]])[[1L]],
        expected = "ENSG00000223972"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneIdVersion"]])[[1L]],
        expected = "ENSG00000223972.5"
    )
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "ccdsId" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneIdNoVersion" = Rle,
            "geneIdVersion" = Rle,
            "geneName" = Rle,
            "havanaGene" = Rle,
            "havanaTranscript" = Rle,
            "hgncId" = Rle,
            "level" = Rle,
            "ont" = "CompressedCharacterList",
            "proteinId" = Rle,
            "source" = Rle,
            "tag" = "CompressedCharacterList",
            "txBiotype" = Rle,
            "txId" = Rle,
            "txIdNoVersion" = Rle,
            "txIdVersion" = Rle,
            "txName" = Rle,
            "txSupportLevel" = Rle,
            "type" = Rle
        )
    )
})



context("makeGRangesFromGFF : RefSeq")

skip_if_not(hasInternet())
file <- gffs[["refseq_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "RefSeqGenes")
    expect_identical(length(object), 54651L)
    expect_identical(names(object)[[1L]], "A1BG")
    expect_identical(
        object = lapply(mcols(object[[1L]]), class),
        expected = list(
            "broadClass" = Rle,
            "description" = Rle,
            "endRange" = "CompressedCharacterList",
            "exception" = Rle,
            "gbkey" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSynonym" = "CompressedCharacterList",
            "partial" = Rle,
            "pseudo" = Rle,
            "source" = Rle,
            "startRange" = "CompressedCharacterList",
            "translExcept" = "CompressedCharacterList",
            "type" = Rle
        )
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "RefSeqTranscripts")
    expect_identical(length(object), 163975L)
    expect_identical(
        object = names(object)[[1L]],
        expected = "NM_000014.6"
    )
    expect_identical(
        object = as.character(mcols(object[[1L]])[["txId"]])[[1L]],
        expected = "NM_000014.6"
    )
    expect_identical(
        object = as.character(mcols(object[[1L]])[["geneId"]])[[1L]],
        expected = "A2M"
    )
    AsIs <- "list"  # nolint
    expect_identical(
        object = lapply(mcols(object[[1L]]), class),
        expected = list(
            "broadClass" = Rle,
            "description" = Rle,
            "endRange" = AsIs,
            "exception" = Rle,
            "gbkey" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSynonym" = AsIs,
            "inference" = Rle,
            "modelEvidence" = Rle,
            "partial" = Rle,
            "product" = Rle,
            "pseudo" = Rle,
            "source" = Rle,
            "startRange" = AsIs,
            "tag" = Rle,
            "translExcept" = AsIs,
            "txId" = Rle,
            "txName" = Rle,
            "type" = Rle
        )
    )
})



context("makeGRangesFromGFF : UCSC")

skip_if_not(hasInternet())
file <- gffs[["ucsc_hg38_ensgene_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "UCSCGenes")
    expect_identical(length(object), 64252L)
    expect_identical(names(object)[[1L]], "ENSG00000223972")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "geneId" = Rle,
            "geneName" = Rle
        )
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "UCSCTranscripts")
    expect_identical(length(object), 208239L)
    expect_identical(names(object)[[1L]], "ENST00000456328")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "geneId" = Rle,
            "geneName" = Rle,
            "txBiotype" = Rle,
            "txChrom" = Rle,
            "txEnd" = Rle,
            "txId" = Rle,
            "txName" = Rle,
            "txNumber" = Rle,
            "txStart" = Rle,
            "txStrand" = Rle
        )
    )
})



context("makeGRangesFromGFF : WormBase")

skip_if_not(hasInternet())
file <- gffs[["wormbase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 46934L)
    expect_identical(names(object)[[1L]], "WBGene00022276")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSource" = Rle,
            "source" = Rle,
            "type" = Rle
        )
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "GRanges")
    expect_identical(length(object), 59897L)
    expect_identical(names(object)[[1L]], "Y74C9A.2a.3")
    expect_identical(
        object = lapply(mcols(object), class),
        expected = list(
            "broadClass" = Rle,
            "geneBiotype" = Rle,
            "geneId" = Rle,
            "geneName" = Rle,
            "geneSource" = Rle,
            "source" = Rle,
            "txBiotype" = Rle,
            "txId" = Rle,
            "txName" = Rle,
            "txSource" = Rle,
            "type" = Rle
        )
    )
})
