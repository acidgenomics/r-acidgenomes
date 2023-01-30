## FIXME Ensure that we check seqnames for Ensembl, GENCODE, RefSeq, UCSC
## Note that UCSC is currently failing.

test_that("Unsupported files", {
    for (file in gffs[c(
        "flybase_gff3",
        "wormbase_gff3"
    )]) {
        expect_error(
            object = makeGRangesFromGFF(file = file),
            regexp = "isSupportedGFF"
        )
    }
})

file <- gffs[["ensembl_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_length(object, 62703L)
    expect_named(
        object = object,
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
            "tag" = "CompressedCharacterList",
            "type" = "Rle"
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
            "description" = paste(
                "NFE2 like bZIP transcription factor 2",
                "[Source:HGNC Symbol;Acc:HGNC:7782]"
            ),
            "geneBiotype" = "protein_coding",
            "geneId" = "ENSG00000116044.17",
            "geneIdNoVersion" = "ENSG00000116044",
            "geneIdVersion" = "ENSG00000116044.17",
            "geneName" = "NFE2L2",
            "logicName" = "ensembl_havana_gene_homo_sapiens",
            "source" = "ensembl_havana",
            "tag" = "character(0)",
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
        object = metadata(object)[["url"]],
        expected = file
    )
    expect_true(
        isAFile(metadata(object)[["file"]])
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38.p13"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "0d9311600fa404a31f517b752002dcda"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 108L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "4cf5bfc008e3a9b8d22d7613d989ee72dd7424dca552ce14cc3d5a87e0ca8048" # nolint
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, 252300L)
    expect_named(
        object = object,
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
            "ccdsId" = "CCDS42782.1",
            "description" = paste(
                "NFE2 like bZIP transcription factor 2",
                "[Source:HGNC Symbol;Acc:HGNC:7782]"
            ),
            "geneBiotype" = "protein_coding",
            "geneId" = "ENSG00000116044.17",
            "geneIdNoVersion" = "ENSG00000116044",
            "geneIdVersion" = "ENSG00000116044.17",
            "geneName" = "NFE2L2",
            "logicName" = "ensembl_havana_gene_homo_sapiens",
            "source" = "ensembl_havana",
            "tag" = "c(\"basic\", \"Ensembl_canonical\", \"MANE_Select\")",
            "txBiotype" = "protein_coding",
            "txId" = "ENST00000397062.8",
            "txIdNoVersion" = "ENST00000397062",
            "txIdVersion" = "ENST00000397062.8",
            "txName" = "NFE2L2-201",
            "txSupportLevel" = "1 (assigned to previous version 7)",
            "type" = "mRNA"
        )
    )
})

file <- gffs[["ensembl_grch38_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_length(object, 61552L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["ENSG00000223972.5"]), # nolint
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
        object = metadata(object)[["url"]],
        expected = file
    )
    expect_true(
        isAFile(metadata(object)[["file"]])
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38.p13"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "9e4e55e43a9ecb8e4acfaa2266ac4e77"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 106L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "9f9cc8ea19f03659875a4a3b57220c379f0282b4c361fc7cb637e6f083967871" # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, 246511L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["ENST00000456328.2"]), # nolint
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

file <- gffs[["flybase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes"
    )
    expect_s4_class(object, "FlyBaseGenes")
    expect_length(object, 17872L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["FBgn0031208"]), # nolint
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
        object = metadata(object)[["url"]],
        expected = file
    )
    expect_true(
        isAFile(metadata(object)[["file"]])
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "r6.45"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "8910e1172661d508b43ca384fab0c618"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Drosophila melanogaster"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = "r6.45"
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "7e1d9c6b6ad1afe2458ab4b0d9202b6fddc0376081aff6a52d6ec07a8ef967e4" # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts"
    )
    expect_s4_class(object, "FlyBaseTranscripts")
    expect_length(object, 35664L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["FBtr0475186"]), # nolint
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

file <- gffs[["gencode_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_length(object, 61544L)
    expect_named(
        object = object,
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
            "tag" = "CompressedCharacterList",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENSG00000223972.5"]), # nolint
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
            "tag" = "character(0)",
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
            "genome" = "GRCh38"
        )
    )
    expect_identical(
        object = metadata(object)[["url"]],
        expected = file
    )
    expect_true(
        isAFile(metadata(object)[["file"]])
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "c3655ad23ed3844e7cf025cbb47d1372"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 40L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "aab8237aca2fce38fee709de6728c270fdf0c35e2938fddbb508285aae68349f" # nolint
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_length(object, 246624L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["ENST00000456328.2"]), # nolint
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

file <- gffs[["gencode_grch38_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_length(object, 61544L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["ENSG00000223972.5"]), # nolint
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
            "genome" = "GRCh38"
        )
    )
    expect_identical(
        object = metadata(object)[["url"]],
        expected = file
    )
    expect_true(
        isAFile(metadata(object)[["file"]])
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "14a867b82917c8c3006838c3a5053a3e"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 40L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "6254e1c52470d74c37f5a7969ee5b1b7debdd9d6bbd9bf722e65fc3d873f3104" # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_length(object, 246624L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["ENST00000456328.2"]), # nolint
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

file <- gffs[["refseq_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes"
    )
    expect_s4_class(object, "RefSeqGenes")
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "dbxref" = "CompressedCharacterList",
            "description" = "Rle",
            "entrezId" = "Rle",
            "exception" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "geneSynonym" = "CompressedCharacterList",
            "parentGeneId" = "Rle",
            "partial" = "Rle",
            "pseudo" = "Rle",
            "source" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object[["A1BG"]]),
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
            "dbxref" = "c(\"GeneID:1\", \"HGNC:HGNC:5\", \"MIM:138670\")",
            "description" = "alpha-1-B glycoprotein",
            "entrezId" = "1",
            "exception" = NA_character_,
            "geneBiotype" = "protein_coding",
            "geneId" = "A1BG",
            "geneName" = "A1BG",
            "geneSynonym" = "c(\"A1B\", \"ABG\", \"GAB\", \"HYST2477\")",
            "parentGeneId" = "A1BG",
            "partial" = NA_character_,
            "pseudo" = NA_character_,
            "source" = "BestRefSeq",
            "type" = "gene"
        )
    )
    expect_identical(
        object = lapply(
            X = as.data.frame(object[["AATF"]]),
            FUN = as.character
        ),
        expected = list(
            "seqnames" = c("NC_000017.11", "NT_187614.1"),
            "start" = c("36948954", "1185319"),
            "end" = c("37056871", "1293236"),
            "width" = rep("107918", 2L),
            "strand" = rep("+", 2L),
            "broadClass" = rep("coding", 2L),
            "dbxref" = c(
                "c(\"GeneID:26574\", \"HGNC:HGNC:19235\", \"MIM:608463\")",
                "c(\"GeneID:26574\", \"HGNC:HGNC:19235\", \"MIM:608463\")"
            ),
            "description" =
                rep("apoptosis antagonizing transcription factor", 2L),
            "entrezId" = rep("26574", 2L),
            "exception" = rep(NA_character_, 2L),
            "geneBiotype" = rep("protein_coding", 2L),
            "geneId" = rep("AATF", 2L),
            "geneName" = rep("AATF", 2L),
            "geneSynonym" =
                rep("c(\"BFR2\", \"CHE-1\", \"CHE1\", \"DED\")", 2L),
            "parentGeneId" = c("AATF", "AATF-2"),
            "partial" = rep(NA_character_, 2L),
            "pseudo" = rep(NA_character_, 2L),
            "source" = c("BestRefSeq%2CGnomon", "BestRefSeq"),
            "type" = rep("gene", 2L)
        )
    )
    expect_identical(
        object = as.data.frame(seqinfo(object))["NC_000001.11", , drop = TRUE],
        expected = list(
            "seqlengths" = 248956422L,
            "isCircular" = NA,
            "genome" = "GRCh38.p14"
        )
    )
    expect_identical(
        object = metadata(object)[["url"]],
        expected = file
    )
    expect_true(
        isAFile(metadata(object)[["file"]])
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38.p14"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts"
    )
    expect_s4_class(object, "RefSeqTranscripts")
    expect_true(all(grepl(
        pattern = "^[A-Z]{2}_[0-9]+\\.[0-9]+$",
        x = names(object)
    )))
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "dbxref" = "CompressedCharacterList",
            "description" = "Rle",
            "entrezId" = "Rle",
            "exception" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "geneSynonym" = "CompressedCharacterList",
            "inference" = "Rle",
            "modelEvidence" = "Rle",
            "parentGeneId" = "Rle",
            "partial" = "Rle",
            "product" = "Rle",
            "pseudo" = "Rle",
            "source" = "Rle",
            "tag" = "Rle",
            "txId" = "Rle",
            "txName" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object[["NM_000014.6"]]),
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
            "dbxref" = paste(
                "c(\"Ensembl:ENST00000318602.12\",",
                "\"GeneID:2\",",
                "\"Genbank:NM_000014.6\",",
                "\"HGNC:HGNC:7\",",
                "\"MIM:103950\")"
            ),
            "description" = "alpha-2-macroglobulin",
            "entrezId" = "2",
            "exception" = NA_character_,
            "geneBiotype" = "protein_coding",
            "geneId" = "A2M",
            "geneName" = "A2M",
            "geneSynonym" = "character(0)",
            "inference" = NA_character_,
            "modelEvidence" = NA_character_,
            "parentGeneId" = "A2M",
            "partial" = NA_character_,
            "product" = "alpha-2-macroglobulin, transcript variant 1",
            "pseudo" = NA_character_,
            "source" = "BestRefSeq",
            "tag" = "MANE Select",
            "txId" = "NM_000014.6",
            "txName" = "NM_000014.6",
            "type" = "mRNA"
        )
    )
})

file <- gffs[["refseq_grch38_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes"
    )
    expect_s4_class(object, "RefSeqGenes")
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "dbXref" = "Rle",
            "description" = "Rle",
            "entrezId" = "Rle",
            "exception" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "geneSynonym" = "Rle",
            "parentGeneId" = "Rle",
            "partial" = "Rle",
            "pseudo" = "Rle",
            "source" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object[["A1BG"]]),
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
            "dbXref" = "MIM:138670",
            "description" = "alpha-1-B glycoprotein",
            "entrezId" = "1",
            "exception" = NA_character_,
            "geneBiotype" = "protein_coding",
            "geneId" = "A1BG",
            "geneName" = "A1BG",
            "geneSynonym" = "HYST2477",
            "parentGeneId" = "A1BG",
            "partial" = NA_character_,
            "pseudo" = NA_character_,
            "source" = "BestRefSeq",
            "type" = "gene"
        )
    )
    expect_identical(
        object = lapply(
            X = as.data.frame(object[["AATF"]]),
            FUN = as.character
        ),
        expected = list(
            "seqnames" = c("NC_000017.11", "NT_187614.1"),
            "start" = c("36948954", "1185319"),
            "end" = c("37056871", "1293236"),
            "width" = rep("107918", 2L),
            "strand" = rep("+", 2L),
            "broadClass" = rep("coding", 2L),
            "dbXref" = rep("MIM:608463", 2L),
            "description" =
                rep("apoptosis antagonizing transcription factor", 2L),
            "entrezId" = rep("26574", 2L),
            "exception" = rep(NA_character_, 2L),
            "geneBiotype" = rep("protein_coding", 2L),
            "geneId" = rep("AATF", 2L),
            "geneName" = rep("AATF", 2L),
            "geneSynonym" = rep("DED", 2L),
            "parentGeneId" = c("AATF", "AATF_1"),
            "partial" = rep(NA_character_, 2L),
            "pseudo" = rep(NA_character_, 2L),
            "source" = c("BestRefSeq%2CGnomon", "BestRefSeq"),
            "type" = rep("gene", 2L)
        )
    )
    expect_identical(
        object = as.data.frame(seqinfo(object))["NC_000001.11", , drop = TRUE],
        expected = list(
            "seqlengths" = 248956422L,
            "isCircular" = NA,
            "genome" = "GRCh38.p14"
        )
    )
    expect_identical(
        object = metadata(object)[["url"]],
        expected = file
    )
    expect_true(
        isAFile(metadata(object)[["file"]])
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "GRCh38.p14"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts"
    )
    expect_s4_class(object, "RefSeqTranscripts")
    expect_true(all(grepl(
        pattern = "^[A-Z]{2}_[0-9]+\\.[0-9]+$",
        x = names(object)
    )))
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "Rle",
            "dbXref" = "Rle",
            "description" = "Rle",
            "entrezId" = "Rle",
            "exception" = "Rle",
            "geneBiotype" = "Rle",
            "geneId" = "Rle",
            "geneName" = "Rle",
            "inference" = "Rle",
            "modelEvidence" = "Rle",
            "parentGeneId" = "Rle",
            "partial" = "Rle",
            "product" = "Rle",
            "pseudo" = "Rle",
            "source" = "Rle",
            "tag" = "Rle",
            "txBiotype" = "Rle",
            "txId" = "Rle",
            "txName" = "Rle",
            "type" = "Rle"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object[["NM_000014.6"]]),
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
            "dbXref" = "GeneID:2",
            "description" = "alpha-2-macroglobulin",
            "entrezId" = "2",
            "exception" = NA_character_,
            "geneBiotype" = "protein_coding",
            "geneId" = "A2M",
            "geneName" = "A2M",
            "inference" = NA_character_,
            "modelEvidence" = NA_character_,
            "parentGeneId" = "A2M",
            "partial" = NA_character_,
            "product" = "alpha-2-macroglobulin, transcript variant 1",
            "pseudo" = NA_character_,
            "source" = "BestRefSeq",
            "tag" = "MANE Select",
            "txBiotype" = "mRNA",
            "txId" = "NM_000014.6",
            "txName" = "NM_000014.6",
            "type" = "transcript"
        )
    )
})

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
            X = as.data.frame(object["ENSG00000223972"]), # nolint
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
        object = metadata(object)[["url"]],
        expected = file
    )
    expect_true(
        isAFile(metadata(object)[["file"]])
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
            X = as.data.frame(object["ENST00000456328"]), # nolint
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

file <- gffs[["wormbase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes"
    )
    expect_s4_class(object, "WormBaseGenes")
    expect_length(object, 46925L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["WBGene00022276"]), # nolint
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
        object = metadata(object)[["url"]],
        expected = file
    )
    expect_true(
        isAFile(metadata(object)[["file"]])
    )
    expect_identical(
        object = metadata(object)[["genomeBuild"]],
        expected = "WS284"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "d1bcb4cc329c2eb23460fd53717dea5e"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Caenorhabditis elegans"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = "WS284"
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "c176ff7fcf1c24844c40df18602736ecc5aa6c3105d2d924227e2ba29a9273f0" # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts"
    )
    expect_s4_class(object, "WormBaseTranscripts")
    expect_length(object, 60128L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["txId"]])
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
            X = as.data.frame(object["Y74C9A.2a.3"]), # nolint
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

file <- file.path("cache", "ref-transcripts.gtf")

test_that("GTF genes", {
    object <- makeGRangesFromGFF(
        file = file,
        level = "genes"
    )
    expect_s4_class(object, "EnsemblGenes")
    expect_length(object, n = 60L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["ENSG00000223972"]), # nolint
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
    object <- makeGRangesFromGFF(
        file = file,
        level = "transcripts"
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, n = 167L)
    expect_named(
        object = object,
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
            X = as.data.frame(object["ENST00000456328"]), # nolint
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
