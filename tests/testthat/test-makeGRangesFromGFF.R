context("makeGRangesFromGFF : Ensembl")

skip_if_not(hasInternet())
file <- gffs[["ensembl_grch38_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_identical(length(object), 60664L)
    expect_identical(names(object)[[1L]], "ENSG00000223972.5")
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "geneSource" = "Rle",
            "source" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38.p13"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 104L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "2c1be3a18925eb033052f4019512eab4cb3227b2a87e7766787452fa6ebbb79c"  # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_identical(length(object), 236920L)
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
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "ccdsId" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "geneSource" = "Rle",
            "source" = "Rle",
            "tag" = "Rle",
            "txBiotype" = "Rle",
            "txId" = "Rle",
            "txIdNoVersion" = "Rle",
            "txIdVersion" = "Rle",
            "txName" = "Rle",
            "txSource" = "Rle",
            "txSupportLevel" = "Rle",
            "type" = "Rle"
        )
    )
})

file <- gffs[["ensembl_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_identical(length(object), 60664L)
    expect_identical(names(object)[[1L]], "ENSG00000223972.5")
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "description" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "logicName" = "Rle",
            "source" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38.p13"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 104L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "313ad46bd4af78b45b9f5d8407bbcbd3f87f4be0747060e84b3b5eb931530ec1"  # nolint
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_identical(length(object), 236918L)
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
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "ccdsId" = "Rle",
            "description" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "logicName" = "Rle",
            "source" = "Rle",
            "tag" = "Rle",
            "txBiotype" = "Rle",
            "txId" = "Rle",
            "txIdNoVersion" = "Rle",
            "txIdVersion" = "Rle",
            "txName" = "Rle",
            "txSupportLevel" = "Rle",
            "type" = "Rle"
        )
    )
})



context("makeGRangesFromGFF : FlyBase")

skip_if_not(hasInternet())
file <- gffs[["flybase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "FlyBaseGenes")
    expect_identical(length(object), 17874L)
    expect_identical(names(object)[[1L]], "FBgn0031208")
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "source" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "r6.40"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Drosophila melanogaster"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = "r6.40"
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "6530ec492803869afec01a5c76a9e47e44e50955fe1cfee2b965687d4b01d99d"  # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "FlyBaseTranscripts")
    expect_identical(length(object), 35642L)
    expect_identical(names(object)[[1L]], "FBtr0475186")
    expect_identical(
        object = as.character(mcols(object)[["txId"]])[[1L]],
        expected = "FBtr0475186"
    )
    expect_identical(
        object = as.character(mcols(object)[["geneId"]])[[1L]],
        expected = "FBgn0031208"
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "source" = "Rle",
            "txId" = "Rle",
            "txName" = "Rle",
            "type" = "Rle"
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
    expect_identical(length(object), 60649L)
    expect_identical(names(object)[[1L]], "ENSG00000223972.5")
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "havanaGene" = "Rle",
            "hgncId" = "Rle",
            "level" = "Rle",
            "source" = "Rle",
            "tag" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 38L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "22020df0d3356e965868f4b193e89fa13e838b950a574349f7fcd461ac01c050"  # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_identical(length(object), 237012L)
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
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "ccdsId" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "havanaGene" = "Rle",
            "havanaTranscript" = "Rle",
            "hgncId" = "Rle",
            "level" = "Rle",
            "ont" = "Rle",
            "proteinId" = "Rle",
            "source" = "Rle",
            "tag" = "Rle",
            "txBiotype" = "Rle",
            "txId" = "Rle",
            "txIdNoVersion" = "Rle",
            "txIdVersion" = "Rle",
            "txName" = "Rle",
            "txSupportLevel" = "Rle",
            "type" = "Rle"
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
    expect_identical(length(object), 60649L)
    expect_identical(names(object)[[1L]], "ENSG00000223972.5")
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "havanaGene" = "Rle",
            "hgncId" = "Rle",
            "level" = "Rle",
            "ont" = "CompressedCharacterList",
            "source" = "Rle",
            "tag" = "CompressedCharacterList",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 38L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "177d45a232e299918785ab78d01b65ef18fac1e3017d96d7b06b9964e9bc997b"  # nolint
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_identical(length(object), 237012L)
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
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "ccdsId" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdNoVersion" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "havanaGene" = "Rle",
            "havanaTranscript" = "Rle",
            "hgncId" = "Rle",
            "level" = "Rle",
            "ont" = "CompressedCharacterList",
            "proteinId" = "Rle",
            "source" = "Rle",
            "tag" = "CompressedCharacterList",
            "txBiotype" = "Rle",
            "txId" = "Rle",
            "txIdNoVersion" = "Rle",
            "txIdVersion" = "Rle",
            "txName" = "Rle",
            "txSupportLevel" = "Rle",
            "type" = "Rle"
        )
    )
})



context("makeGRangesFromGFF : RefSeq")

skip_if_not(hasInternet())
file <- gffs[["refseq_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "RefSeqGenes")
    ## This changes over time, so don't hard-code (2021-08-05).
    ## > expect_identical(length(object), 54583L)
    expect_true(isSubset("A1BG", names(object)))
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "description" = "Rle",
            "endRange" = "CompressedCharacterList",
            "exception" = "Rle",
            "experiment" = "CompressedCharacterList",
            "gbkey" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "geneSynonym" = "CompressedCharacterList",
            "partial" = "Rle",
            "pseudo" = "Rle",
            "source" = "Rle",
            "startRange" = "CompressedCharacterList",
            "translExcept" = "CompressedCharacterList",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38.p13"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "RefSeqTranscripts")
    ## This changes over time, so don't hard-code (2021-08-05).
    ## > expect_identical(length(object), 163975L)
    expect_true(isSubset("NM_000014.6", names(object)))
    expect_identical(
        object = as.character(mcols(object[["NM_000014.6"]])[["txId"]])[[1L]],
        expected = "NM_000014.6"
    )
    expect_identical(
        object = as.character(mcols(object[["NM_000014.6"]])[["geneId"]])[[1L]],
        expected = "A2M"
    )
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "description" = "Rle",
            "endRange" = "CompressedCharacterList",
            "exception" = "Rle",
            "experiment" = "CompressedCharacterList",
            "gbkey" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "geneSynonym" = "CompressedCharacterList",
            "inference" = "Rle",
            "modelEvidence" = "Rle",
            "partial" = "Rle",
            "product" = "Rle",
            "pseudo" = "Rle",
            "source" = "Rle",
            "startRange" = "CompressedCharacterList",
            "tag" = "Rle",
            "translExcept" = "CompressedCharacterList",
            "txId" = "Rle",
            "txName" = "Rle",
            "type" = "Rle"
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
    ## This changes over time, so don't hard-code (2021-08-05).
    ## > expect_identical(length(object), 64252L)
    expect_true(isSubset("ENSG00000223972", names(object)))
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "geneId" = "Rle",
            "geneName" = "Rle"
        )
    )
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "hg38"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "UCSCTranscripts")
    ## This changes over time, so don't hard-code (2021-08-05).
    ## > expect_identical(length(object), 208239L)
    expect_true(isSubset("ENST00000456328", names(object)))
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "geneId" = "Rle",
            "geneName" = "Rle",
            "txBiotype" = "Rle",
            "txChrom" = "Rle",
            "txEnd" = "Rle",
            "txId" = "Rle",
            "txName" = "Rle",
            "txNumber" = "Rle",
            "txStart" = "Rle",
            "txStrand" = "Rle"
        )
    )
})



context("makeGRangesFromGFF : WormBase")

## FIXME Drop back to checking against WS280?
## FIXME WS280 is dropping "geneName" now as well...need to fix.

skip_if_not(hasInternet())
file <- gffs[["wormbase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "WormBaseGenes")
    expect_identical(length(object), 46925L)
    expect_identical(names(object)[[1L]], "WBGene00022276")
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "geneSource" = "Rle",
            "geneVersion" = "Rle",
            "source" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "WS281"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Caenorhabditis elegans"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = "WS281"
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "ea539b392680ce78db96d9d2f3ed2bb0e64f875685bc57a95022d5914a16a403"  # nolint
    )
})

## FIXME Now this is dropping gene name argh...
## FIXME How to get WormBase 281 to keep track of gene name???
test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "WormBaseTranscripts")
    expect_identical(length(object), 59961L)
    expect_identical(names(object)[[1L]], "Y74C9A.2a.3")
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "geneSource" = "Rle",
            "geneVersion" = "Rle",
            "source" = "Rle",
            "txBiotype" = "Rle",
            "txId" = "Rle",
            "txName" = "Rle",
            "txSource" = "Rle",
            "type" = "Rle"
        )
    )
})



context("makeGRangesFromGFF : bcbio")

test_that("bcbio `ref-transcripts.gtf` file", {
    file <- file.path("cache", "ref-transcripts.gtf")
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "EnsemblGenes")
    expect_length(object, n = 60L)
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, n = 167L)
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
})
