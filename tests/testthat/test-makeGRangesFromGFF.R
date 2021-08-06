## FIXME Need to check seqinfo here.
## FIXME RefSeq seqinfo is currently failing...what's up with that?
## FIXME Need to improve consistency of "geneId" and "geneName" checks.
## FIXME Check seqinfo, seqlevels, seqlengths.
## FIXME Automatic assignment of Seqinfo failed for RefSeq...need to resolve?

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
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["geneId"]])
    )
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
        object = vapply(
            X = as.data.frame(object["ENSG00000223972.5"]),
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
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972.5",
            "geneIdNoVersion" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "geneSource" = "havana",
            "source" = "havana",
            "type" = "gene"
        )
    )
    expect_identical(
        object = levels(seqnames(object))[seq_len(25L)],
        expected = c(
            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
            "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"
        )
    )
    expect_identical(
        object = as.data.frame(seqinfo(object))["1", , drop = TRUE],
        expected = list(
            "seqlengths" = 248956422L,
            "isCircular" = NA,
            "genome" = "GRCh38.p13"
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
        object = metadata(object)[["md5"]],
        expected = "1615e866df72f3396ae29714469d5bb4"
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
        object = names(object),
        expected = as.character(mcols(object)[["txId"]])
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
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000456328.2"]),
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
            "ccdsId" = NA_character_,
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972.5",
            "geneIdNoVersion" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "geneSource" = "havana",
            "source" = "havana",
            "tag" = "basic",
            "txBiotype" = "processed_transcript",
            "txId" = "ENST00000456328.2",
            "txIdNoVersion" = "ENST00000456328",
            "txIdVersion" = "ENST00000456328.2",
            "txName" = "DDX11L1-202",
            "txSource" = "havana",
            "txSupportLevel" = "1",
            "type" = "transcript"
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
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["geneId"]])
    )
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
        object = vapply(
            X = as.data.frame(object["ENSG00000223972.5"]),
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
            "description" = paste(
                "DEAD/H-box helicase 11 like 1 (pseudogene)",
                "[Source:HGNC Symbol;Acc:HGNC:37102]"
            ),
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972.5",
            "geneIdNoVersion" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "logicName" = "havana_homo_sapiens",
            "source" = "havana",
            "type" = "pseudogene"
        )
    )
    expect_identical(
        object = levels(seqnames(object))[seq_len(25L)],
        expected = c(
            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
            "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"
        )
    )
    expect_identical(
        object = as.data.frame(seqinfo(object))["1", , drop = TRUE],
        expected = list(
            "seqlengths" = 248956422L,
            "isCircular" = NA,
            "genome" = "GRCh38.p13"
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
        object = metadata(object)[["md5"]],
        expected = "d9fd21961d3af5d25b39f724dc93ef3e"
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
        object = names(object),
        expected = as.character(mcols(object)[["txId"]])
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
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000456328.2"]),
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
            "ccdsId" = NA_character_,
            "description" = paste(
                "DEAD/H-box helicase 11 like 1 (pseudogene)",
                "[Source:HGNC Symbol;Acc:HGNC:37102]"
            ),
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972.5",
            "geneIdNoVersion" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "logicName" = "havana_homo_sapiens",
            "source" = "havana",
            "tag" = "basic",
            "txBiotype" = "processed_transcript",
            "txId" = "ENST00000456328.2",
            "txIdNoVersion" = "ENST00000456328",
            "txIdVersion" = "ENST00000456328.2",
            "txName" = "DDX11L1-202",
            "txSupportLevel" = "1",
            "type" = "lnc_RNA"
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
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["geneId"]])
    )
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
        object = vapply(
            X = as.data.frame(object["FBgn0031208"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "2L",
            "start" = "7529",
            "end" = "9484",
            "width" = "1956",
            "strand" = "+",
            "broadClass" = "other",
            "geneId" = "FBgn0031208",
            "geneName" = "CR11023",
            "source" = "FlyBase",
            "type" = "gene"
        )
    )
    expect_true(isSubset(
        x = c(
            "2L", "2R", "3L", "3R", "4", "X", "Y",
            "mitochondrion_genome", "rDNA"
        ),
        y = levels(seqnames(object))
    ))
    expect_true(all(is.na(seqlengths(object))))
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "r6.40"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "3563ac20aa9c6605a34227952c8e70bb"
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
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["txId"]])
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
    expect_identical(
        object = vapply(
            X = as.data.frame(object["FBtr0475186"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "2L",
            "start" = "7529",
            "end" = "9484",
            "width" = "1956",
            "strand" = "+",
            "broadClass" = "other",
            "geneId" = "FBgn0031208",
            "geneName" = "CR11023",
            "source" = "FlyBase",
            "txId" = "FBtr0475186",
            "txName" = "CR11023-RE",
            "type" = "pseudogene"
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
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["geneId"]])
    )
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
        object = vapply(
            X = as.data.frame(object["ENSG00000223972.5"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "11869",
            "end" = "14409",
            "width" = "2541",
            "strand" = "+",
            "broadClass" = "pseudo",
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972.5",
            "geneIdNoVersion" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "havanaGene" = "OTTHUMG00000000961.2",
            "hgncId" = "HGNC:37102",
            "level" = "2",
            "source" = "HAVANA",
            "tag" = NA_character_,
            "type" = "gene"
        )
    )
    expect_identical(
        object = levels(seqnames(object)),
        expected = paste0(
            "chr",
            c(
                seq(from = 1L, to = 22L),
                "X", "Y", "M"
            )
        )
    )
    expect_identical(
        object = as.data.frame(seqinfo(object))["chr1", , drop = TRUE],
        expected = list(
            "seqlengths" = 248956422L,
            "isCircular" = FALSE,
            "genome" = "hg38"
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
        object = metadata(object)[["md5"]],
        expected = "16fcae8ca8e488cd8056cf317d963407"
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
        object = names(object),
        expected = as.character(mcols(object)[["txId"]])
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
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000456328.2"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "11869",
            "end" = "14409",
            "width" = "2541",
            "strand" = "+",
            "broadClass" = "pseudo",
            "ccdsId" = NA_character_,
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972.5",
            "geneIdNoVersion" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "havanaGene" = "OTTHUMG00000000961.2",
            "havanaTranscript" = "OTTHUMT00000362751.1",
            "hgncId" = "HGNC:37102",
            "level" = "2",
            "ont" = NA_character_,
            "proteinId" = NA_character_,
            "source" = "HAVANA",
            "tag" = "basic",
            "txBiotype" = "processed_transcript",
            "txId" = "ENST00000456328.2",
            "txIdNoVersion" = "ENST00000456328",
            "txIdVersion" = "ENST00000456328.2",
            "txName" = "DDX11L1-202",
            "txSupportLevel" = "1",
            "type" = "transcript"
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
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["geneId"]])
    )
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
        object = vapply(
            X = as.data.frame(object["ENSG00000223972.5"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "11869",
            "end" = "14409",
            "width" = "2541",
            "strand" = "+",
            "broadClass" = "pseudo",
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972.5",
            "geneIdNoVersion" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "havanaGene" = "OTTHUMG00000000961.2",
            "hgncId" = "HGNC:37102",
            "level" = "2",
            "ont" = "character(0)",
            "source" = "HAVANA",
            "tag" = "character(0)",
            "type" = "gene"
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
        object = metadata(object)[["md5"]],
        expected = "9e1370a6002ee2eb27ad4a55f305de6c"
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
        object = names(object),
        expected = as.character(mcols(object)[["txId"]])
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
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000456328.2"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "11869",
            "end" = "14409",
            "width" = "2541",
            "strand" = "+",
            "broadClass" = "pseudo",
            "ccdsId" = NA_character_,
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972.5",
            "geneIdNoVersion" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "havanaGene" = "OTTHUMG00000000961.2",
            "havanaTranscript" = "OTTHUMT00000362751.1",
            "hgncId" = "HGNC:37102",
            "level" = "2",
            "ont" = "character(0)",
            "proteinId" = NA_character_,
            "source" = "HAVANA",
            "tag" = "basic",
            "txBiotype" = "processed_transcript",
            "txId" = "ENST00000456328.2",
            "txIdNoVersion" = "ENST00000456328",
            "txIdVersion" = "ENST00000456328.2",
            "txName" = "DDX11L1-202",
            "txSupportLevel" = "1",
            "type" = "transcript"
        )
    )
})



context("makeGRangesFromGFF : RefSeq")

skip_if_not(hasInternet())
file <- gffs[["refseq_grch38_gff3"]]

## FIXME Check the seqinfo handling here.
## → Getting `Seqinfo` from 5c0b178f46eb_GCF_000001405.39_GRCh38.p13_assembly_report.txt.
## ! Automatic `seqinfo` assignment failed.
## > seqnames(object[[1L]])
## factor-Rle of length 1 with 1 run
## Lengths:            1
## Values : NC_000012.12
## Levels(639): NC_000001.11 ... NW_021160031.1
## FIXME Check the levels here...

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "RefSeqGenes")
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
        object = vapply(
            X = as.data.frame(object[["A1BG"]][1L]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "NC_000019.10",
            "start" = "58345183",
            "end" = "58353492",
            "width" = "8310",
            "strand" = "-",
            "broadClass" = "coding",
            "description" = "alpha-1-B glycoprotein",
            "endRange" = "character(0)",
            "exception" = NA_character_,
            "experiment" = "character(0)",
            "gbkey" = "Gene",
            "geneBiotype" = "protein_coding",
            "geneId" = "A1BG",
            "geneName" = "A1BG",
            "geneSynonym" = "c(\"A1B\", \"ABG\", \"GAB\", \"HYST2477\")",
            "partial" = NA_character_,
            "pseudo" = NA_character_,
            "source" = "BestRefSeq",
            "startRange" = "character(0)",
            "translExcept" = "character(0)",
            "type" = "gene"
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
    expect_identical(
        object = vapply(
            X = as.data.frame(object[["NM_000014.6"]][1L]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "NC_000012.12",
            "start" = "9067708",
            "end" = "9115919",
            "width" = "48212",
            "strand" = "-",
            "broadClass" = "coding",
            "description" = "alpha-2-macroglobulin",
            "endRange" = "character(0)",
            "exception" = NA_character_,
            "experiment" = "character(0)",
            "gbkey" = "mRNA",
            "geneBiotype" = "protein_coding",
            "geneId" = "A2M",
            "geneName" = "A2M",
            "geneSynonym" = "character(0)",
            "inference" = NA_character_,
            "modelEvidence" = NA_character_,
            "partial" = NA_character_,
            "product" = "alpha-2-macroglobulin, transcript variant 1",
            "pseudo" = NA_character_,
            "source" = "BestRefSeq",
            "startRange" = "character(0)",
            "tag" = "MANE Select",
            "translExcept" = "character(0)",
            "txId" = "NM_000014.6",
            "txName" = "NM_000014.6",
            "type" = "mRNA"
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

    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "geneId" = "Rle",
            "geneName" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENSG00000223972"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "11869",
            "end" = "14409",
            "width" = "2541",
            "strand" = "+",
            "geneId" = "ENSG00000223972",
            "geneName" = "ENSG00000223972"
        )
    )
    expect_identical(
        object = as.data.frame(seqinfo(object))["chr1", , drop = TRUE],
        expected = list(
            "seqlengths" = 248956422L,
            "isCircular" = FALSE,
            "genome" = "hg38"
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
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000456328"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "11869",
            "end" = "14409",
            "width" = "2541",
            "strand" = "+",
            "geneId" = "ENSG00000223972",
            "geneName" = "ENSG00000223972",
            "txBiotype" = "transcript",
            "txChrom" = "chr1",
            "txEnd" = "14409",
            "txId" = "ENST00000456328",
            "txName" = "ENST00000456328",
            "txNumber" = "1",
            "txStart" = "11869",
            "txStrand" = "+"
        )
    )
})



context("makeGRangesFromGFF : WormBase")

skip_if_not(hasInternet())
file <- gffs[["wormbase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "WormBaseGenes")
    expect_identical(length(object), 46925L)
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["geneId"]])
    )
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
        object = vapply(
            X = as.data.frame(object["WBGene00022276"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "I",
            "start" = "11495",
            "end" = "16837",
            "width" = "5343",
            "strand" = "+",
            "broadClass" = "coding",
            "geneBiotype" = "protein_coding",
            "geneId" = "WBGene00022276",
            "geneName" = "nlp-40",
            "geneSource" = "WormBase",
            "geneVersion" = "1",
            "source" = "WormBase",
            "type" = "gene"
        )
    )
    expect_identical(
        object = levels(seqnames(object)),
        expected = c(
            "I",
            "II",
            "III",
            "IV",
            "V",
            "X",
            "MtDNA"
        )
    )
    expect_true(all(is.na(seqlengths(object))))
    expect_identical(
        object = metadata(object)[["file"]],
        expected = file
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "WS281"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "a63345f3410927caa1029a143b7e6acc"
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

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "WormBaseTranscripts")
    expect_identical(length(object), 59961L)
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["txId"]]))
    )
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
    expect_identical(
        object = vapply(
            X = as.data.frame(object["Y74C9A.2a.3"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "I",
            "start" = "11495",
            "end" = "16793",
            "width" = "5299",
            "strand" = "+",
            "broadClass" = "coding",
            "geneBiotype" = "protein_coding",
            "geneId" = "WBGene00022276",
            "geneName" = "nlp-40",
            "geneSource" = "WormBase",
            "geneVersion" = "1",
            "source" = "WormBase",
            "txBiotype" = "protein_coding",
            "txId" = "Y74C9A.2a.3",
            "txName" = "Y74C9A.2a.3",
            "txSource" = "WormBase",
            "type" = "transcript"
        )
    )
})



context("makeGRangesFromGFF : bcbio")

file <- file.path("cache", "ref-transcripts.gtf")

test_that("GTF genes", {
    object <- makeGRangesFromGFF(file = file, level = "genes")
    expect_s4_class(object, "EnsemblGenes")
    expect_length(object, n = 60L)
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["geneId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "geneSource" = "Rle",
            "source" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENSG00000223972"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "11869",
            "end" = "14409",
            "width" = "2541",
            "strand" = "+",
            "broadClass" = "pseudo",
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "geneSource" = "havana",
            "source" = "havana",
            "type" = "gene"
        )
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(file = file, level = "transcripts")
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, n = 167L)
    expect_identical(
        object = names(object),
        expected = as.character(mcols(object)[["txId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "ccdsId" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneIdVersion" = "Rle",
            "geneName" = "Rle",
            "geneSource" = "Rle",
            "source" = "Rle",
            "tag" = "Rle",
            "txBiotype" = "Rle",
            "txId" = "Rle",
            "txIdVersion" = "Rle",
            "txName" = "Rle",
            "txSource" = "Rle",
            "txSupportLevel" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000456328"]),
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "11869",
            "end" = "14409",
            "width" = "2541",
            "strand" = "+",
            "broadClass" = "pseudo",
            "ccdsId" = NA_character_,
            "geneBiotype" = "transcribed_unprocessed_pseudogene",
            "geneId" = "ENSG00000223972",
            "geneIdVersion" = "ENSG00000223972.5",
            "geneName" = "DDX11L1",
            "geneSource" = "havana",
            "source" = "havana",
            "tag" = "basic",
            "txBiotype" = "processed_transcript",
            "txId" = "ENST00000456328",
            "txIdVersion" = "ENST00000456328.2",
            "txName" = "DDX11L1-202",
            "txSource" = "havana",
            "txSupportLevel" = "1",
            "type" = "transcript"
        )
    )
})
