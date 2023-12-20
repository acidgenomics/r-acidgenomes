## FIXME Ensure we cover scaffold GTF files.

test_that("Unsupported files", {
    for (file in gffs[c("flybase_gff3", "wormbase_gff3")]) {
        expect_error(
            object = makeGRangesFromGff(file = file),
            regexp = "isSupportedGff"
        )
    }
})

file <- gffs[["ensembl_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGff(
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
            "broadClass" = "factor",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSynonyms" = "CompressedCharacterList",
            "logicName" = "factor",
            "ncbiGeneId" = "CompressedIntegerList",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "type" = "factor"
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
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
            "logicName" = "ensembl_havana_gene_homo_sapiens",
            "ncbiGeneId" = "4780",
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
    object <- makeGRangesFromGff(
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
            "broadClass" = "factor",
            "ccdsId" = "character",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSynonyms" = "CompressedCharacterList",
            "logicName" = "factor",
            "ncbiGeneId" = "CompressedIntegerList",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "txBiotype" = "factor",
            "txId" = "character",
            "txIdNoVersion" = "character",
            "txIdVersion" = "character",
            "txName" = "character",
            "txSupportLevel" = "factor",
            "type" = "factor"
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
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
            "logicName" = "ensembl_havana_gene_homo_sapiens",
            "ncbiGeneId" = "4780",
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
    object <- makeGRangesFromGff(
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
            "broadClass" = "factor",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSource" = "factor",
            "geneSynonyms" = "CompressedCharacterList",
            "ncbiGeneId" = "CompressedIntegerList",
            "source" = "factor",
            "type" = "factor"
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
            "geneSource" = "ensembl_havana",
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
            "ncbiGeneId" = "4780",
            "source" = "ensembl_havana",
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
        expected = "27864327a647f9fdae2cafff5f3ea45e"
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
        expected = "1ae26267acf522d8d9a94736e2c00b5ff1008d43623fc73644092016ae0c57e8" # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGff(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "EnsemblTranscripts")
    expect_length(object, 252301L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["txId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "ccdsId" = "character",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSource" = "factor",
            "geneSynonyms" = "CompressedCharacterList",
            "ncbiGeneId" = "CompressedIntegerList",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "txBiotype" = "factor",
            "txId" = "character",
            "txIdNoVersion" = "character",
            "txIdVersion" = "character",
            "txName" = "character",
            "txSource" = "factor",
            "txSupportLevel" = "factor",
            "type" = "factor"
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
            "ccdsId" = "CCDS42782",
            "description" = paste(
                "NFE2 like bZIP transcription factor 2",
                "[Source:HGNC Symbol;Acc:HGNC:7782]"
            ),
            "geneBiotype" = "protein_coding",
            "geneId" = "ENSG00000116044.17",
            "geneIdNoVersion" = "ENSG00000116044",
            "geneIdVersion" = "ENSG00000116044.17",
            "geneName" = "NFE2L2",
            "geneSource" = "ensembl_havana",
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
            "ncbiGeneId" = "4780",
            "source" = "ensembl_havana",
            "tag" = "MANE_Select",
            "txBiotype" = "protein_coding",
            "txId" = "ENST00000397062.8",
            "txIdNoVersion" = "ENST00000397062",
            "txIdVersion" = "ENST00000397062.8",
            "txName" = "NFE2L2-201",
            "txSource" = "ensembl_havana",
            "txSupportLevel" = "1 (assigned to previous version 7)",
            "type" = "transcript"
        )
    )
})

file <- gffs[["flybase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGff(
        file = file,
        level = "genes"
    )
    expect_s4_class(object, "FlybaseGenes")
    expect_length(object, 17896L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["geneId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "geneId" = "character",
            "geneName" = "character",
            "source" = "factor",
            "type" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["FBgn0262975"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "3R",
            "start" = "23185580",
            "end" = "23226711",
            "width" = "41132",
            "strand" = "-",
            "broadClass" = "other",
            "geneId" = "FBgn0262975",
            "geneName" = "cnc",
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
        expected = "r6.49"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "722ba353d2dd6d3036d747f1498a2e9f"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Drosophila melanogaster"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = "r6.49"
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "565b486cbfe78942127482018e4439a28a6b171d71624d8ffe2a3611e162a101" # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGff(
        file = file,
        level = "transcripts"
    )
    expect_s4_class(object, "FlybaseTranscripts")
    expect_length(object, 35723L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["txId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "geneId" = "character",
            "geneName" = "character",
            "source" = "factor",
            "txId" = "character",
            "txName" = "character",
            "type" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["FBtr0306744"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "3R",
            "start" = "23185580",
            "end" = "23191120",
            "width" = "5541",
            "strand" = "-",
            "broadClass" = "other",
            "geneId" = "FBgn0262975",
            "geneName" = "cnc",
            "source" = "FlyBase",
            "txId" = "FBtr0306744",
            "txName" = "cnc-RF",
            "type" = "mRNA"
        )
    )
})

file <- gffs[["gencode_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGff(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_length(object, 62696L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["geneId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "artifactualDuplication" = "character",
            "broadClass" = "factor",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSynonyms" = "CompressedCharacterList",
            "havanaGene" = "character",
            "hgncId" = "integer",
            "level" = "factor",
            "ncbiGeneId" = "CompressedIntegerList",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "type" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENSG00000116044.17"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr2",
            "start" = "177218667",
            "end" = "177392756",
            "width" = "174090",
            "strand" = "-",
            "artifactualDuplication" = NA_character_,
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
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
            "havanaGene" = "OTTHUMG00000133620.18",
            "hgncId" = "7782",
            "level" = "1",
            "ncbiGeneId" = "4780",
            "source" = "HAVANA",
            "tag" = "c(\"ncRNA_host\", \"overlapping_locus\")",
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
            "isCircular" = NA,
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
        expected = "b0199f8b1522f61093896e8d48750a0c"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 42L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "10d01b26b75d0142677d0bddf369af375c2e24f6b4c99ce04f408032c51ca432" # nolint
    )
})

test_that("GFF3 transcripts", {
    object <- makeGRangesFromGff(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_length(object, 252416L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["txId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "ccdsId" = "character",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSynonyms" = "CompressedCharacterList",
            "havanaGene" = "character",
            "havanaTranscript" = "character",
            "hgncId" = "integer",
            "level" = "factor",
            "ncbiGeneId" = "CompressedIntegerList",
            "proteinId" = "character",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "txBiotype" = "factor",
            "txId" = "character",
            "txIdNoVersion" = "character",
            "txIdVersion" = "character",
            "txName" = "character",
            "txSupportLevel" = "character",
            "type" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000397062.8"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr2",
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
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
            "havanaGene" = "OTTHUMG00000133620.18",
            "havanaTranscript" = "OTTHUMT00000257752.5",
            "hgncId" = "7782",
            "level" = "2",
            "ncbiGeneId" = "4780",
            "proteinId" = "ENSP00000380252.3",
            "source" = "HAVANA",
            "tag" = paste(
                "c(\"basic\",",
                "\"Ensembl_canonical\",",
                "\"MANE_Select\",",
                "\"appris_alternative_1\",",
                "\"CCDS\")"
            ),
            "txBiotype" = "protein_coding",
            "txId" = "ENST00000397062.8",
            "txIdNoVersion" = "ENST00000397062",
            "txIdVersion" = "ENST00000397062.8",
            "txName" = "NFE2L2-201",
            "txSupportLevel" = "1",
            "type" = "transcript"
        )
    )
})

file <- gffs[["gencode_grch38_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGff(
        file = file,
        level = "genes",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeGenes")
    expect_length(object, 62696L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["geneId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "artifactualDuplication" = "character",
            "broadClass" = "factor",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSynonyms" = "CompressedCharacterList",
            "havanaGene" = "character",
            "hgncId" = "integer",
            "level" = "factor",
            "ncbiGeneId" = "CompressedIntegerList",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "type" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENSG00000116044.17"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr2",
            "start" = "177218667",
            "end" = "177392756",
            "width" = "174090",
            "strand" = "-",
            "artifactualDuplication" = NA_character_,
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
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
            "havanaGene" = "OTTHUMG00000133620.18",
            "hgncId" = "7782",
            "level" = "1",
            "ncbiGeneId" = "4780",
            "source" = "HAVANA",
            "tag" = "c(\"overlapping_locus\")",
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
            "isCircular" = NA,
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
        expected = "b3a1bb9b1239d8075dd2b031817a36da"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = 42L
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "f59a5b38c6de1472430ccd0181927a01df0b4ad8d78f598e521badc9da80a76e" # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGff(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "GencodeTranscripts")
    expect_length(object, 252416L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["txId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "ccdsId" = "character",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdNoVersion" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSynonyms" = "CompressedCharacterList",
            "havanaGene" = "character",
            "havanaTranscript" = "character",
            "hgncId" = "integer",
            "level" = "factor",
            "ncbiGeneId" = "CompressedIntegerList",
            "proteinId" = "character",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "txBiotype" = "factor",
            "txId" = "character",
            "txIdNoVersion" = "character",
            "txIdVersion" = "character",
            "txName" = "character",
            "txSupportLevel" = "factor",
            "type" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000397062.8"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr2",
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
            "geneSynonyms" = "c(\"NRF-2\", \"NRF2\")",
            "havanaGene" = "OTTHUMG00000133620.18",
            "havanaTranscript" = "OTTHUMT00000257752.5",
            "hgncId" = "7782",
            "level" = "2",
            "ncbiGeneId" = "4780",
            "proteinId" = "ENSP00000380252.3",
            "source" = "HAVANA",
            "tag" = "c(\"CCDS\")",
            "txBiotype" = "protein_coding",
            "txId" = "ENST00000397062.8",
            "txIdNoVersion" = "ENST00000397062",
            "txIdVersion" = "ENST00000397062.8",
            "txName" = "NFE2L2-201",
            "txSupportLevel" = "1",
            "type" = "transcript"
        )
    )
})

file <- gffs[["refseq_grch38_gff3"]]

test_that("GFF3 genes", {
    object <- makeGRangesFromGff(
        file = file,
        level = "genes"
    )
    expect_s4_class(object, "RefseqGenes")
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "dbXref" = "CompressedCharacterList",
            "description" = "character",
            "exception" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneName" = "character",
            "geneSynonym" = "CompressedCharacterList",
            "ncbiGeneId" = "integer",
            "parentGeneId" = "character",
            "partial" = "logical",
            "pseudo" = "logical",
            "source" = "factor",
            "type" = "factor"
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
            "dbXref" = c(
                "c(\"GeneID:26574\", \"HGNC:HGNC:19235\", \"MIM:608463\")",
                "c(\"GeneID:26574\", \"HGNC:HGNC:19235\", \"MIM:608463\")"
            ),
            "description" =
                rep("apoptosis antagonizing transcription factor", 2L),
            "exception" = rep(NA_character_, 2L),
            "geneBiotype" = rep("protein_coding", 2L),
            "geneId" = rep("AATF", 2L),
            "geneName" = rep("AATF", 2L),
            "geneSynonym" =
                rep("c(\"BFR2\", \"CHE-1\", \"CHE1\", \"DED\")", 2L),
            "ncbiGeneId" = rep("26574", 2L),
            "parentGeneId" = c("AATF", "AATF-2"),
            "partial" = rep("FALSE", 2L),
            "pseudo" = rep("FALSE", 2L),
            "source" = rep("BestRefSeq/Gnomon", 2L),
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
    object <- makeGRangesFromGff(
        file = file,
        level = "transcripts"
    )
    expect_s4_class(object, "RefseqTranscripts")
    expect_true(all(grepl(
        pattern = "^[A-Z]{2}_[0-9]+\\.[0-9]+$",
        x = names(object)
    )))
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "dbXref" = "CompressedCharacterList",
            "description" = "character",
            "exception" = "factor",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneName" = "character",
            "geneSynonym" = "CompressedCharacterList",
            "inference" = "character",
            "modelEvidence" = "character",
            "ncbiGeneId" = "integer",
            "parentGeneId" = "character",
            "partial" = "logical",
            "product" = "character",
            "pseudo" = "logical",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "txId" = "character",
            "txName" = "character",
            "type" = "factor"
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
            "dbXref" = paste(
                "c(\"Ensembl:ENST00000318602.12\",",
                "\"GeneID:2\",",
                "\"GenBank:NM_000014.6\",",
                "\"HGNC:HGNC:7\",",
                "\"MIM:103950\")"
            ),
            "description" = "alpha-2-macroglobulin",
            "exception" = NA_character_,
            "geneBiotype" = "protein_coding",
            "geneId" = "A2M",
            "geneName" = "A2M",
            "geneSynonym" = "character(0)",
            "inference" = NA_character_,
            "modelEvidence" = NA_character_,
            "ncbiGeneId" = "2",
            "parentGeneId" = "A2M",
            "partial" = "FALSE",
            "product" = "alpha-2-macroglobulin, transcript variant 1",
            "pseudo" = "FALSE",
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
    object <- makeGRangesFromGff(
        file = file,
        level = "genes"
    )
    expect_s4_class(object, "RefseqGenes")
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "dbXref" = "CompressedCharacterList",
            "description" = "character",
            "exception" = "factor",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneName" = "character",
            "geneSynonym" = "character",
            "ncbiGeneId" = "integer",
            "parentGeneId" = "character",
            "partial" = "logical",
            "pseudo" = "logical",
            "source" = "factor",
            "type" = "factor"
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
            "exception" = rep(NA_character_, 2L),
            "geneBiotype" = rep("protein_coding", 2L),
            "geneId" = rep("AATF", 2L),
            "geneName" = rep("AATF", 2L),
            "geneSynonym" = rep("DED", 2L),
            "ncbiGeneId" = rep("26574", 2L),
            "parentGeneId" = c("AATF", "AATF_1"),
            "partial" = rep("FALSE", 2L),
            "pseudo" = rep("FALSE", 2L),
            "source" = c("BestRefSeq/Gnomon", "BestRefSeq/Gnomon"),
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
    object <- makeGRangesFromGff(
        file = file,
        level = "transcripts"
    )
    expect_s4_class(object, "RefseqTranscripts")
    expect_true(all(grepl(
        pattern = "^[A-Z]{2}_[0-9]+\\.[0-9]+$",
        x = names(object)
    )))
    expect_identical(
        object = lapply(mcols(object[[1L]]), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "dbXref" = "CompressedCharacterList",
            "description" = "character",
            "exception" = "factor",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneName" = "character",
            "inference" = "character",
            "modelEvidence" = "character",
            "ncbiGeneId" = "integer",
            "parentGeneId" = "character",
            "partial" = "logical",
            "product" = "character",
            "pseudo" = "logical",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "txBiotype" = "factor",
            "txId" = "character",
            "txName" = "character",
            "type" = "factor"
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
            "exception" = NA_character_,
            "geneBiotype" = "protein_coding",
            "geneId" = "A2M",
            "geneName" = "A2M",
            "inference" = NA_character_,
            "modelEvidence" = NA_character_,
            "ncbiGeneId" = "2",
            "parentGeneId" = "A2M",
            "partial" = FALSE,
            "product" = "alpha-2-macroglobulin, transcript variant 1",
            "pseudo" = FALSE,
            "source" = "BestRefSeq",
            "tag" = "MANE Select",
            "txBiotype" = "mRNA",
            "txId" = "NM_000014.6",
            "txName" = "NM_000014.6",
            "type" = "transcript"
        )
    )
})

file <- gffs[["ucsc_hg38_ncbirefseq_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGff(
        file = file,
        level = "genes",
        ignoreVersion = TRUE
    )
    expect_s4_class(object, "UcscGenes")
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "geneId" = "character",
            "geneName" = "character"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["A1BG"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr19",
            "start" = "58345183",
            "end" = "58353492",
            "width" = "8310",
            "strand" = "-",
            "geneId" = "A1BG",
            "geneName" = "A1BG"
        )
    )
    ## UCSC has changed some things related to hg38 that is currently causing
    ## issues with GenomeInfoDb (v1.34.7).
    ## https://github.com/Bioconductor/GenomeInfoDb/issues/83
    ## > expect_identical(
    ## >     object = as.data.frame(seqinfo(object))["chr1", , drop = TRUE],
    ## >     expected = list(
    ## >         "seqlengths" = 248956422L,
    ## >         "isCircular" = FALSE,
    ## >         "genome" = "hg38"
    ## >     )
    ## > )
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
    object <- makeGRangesFromGff(
        file = file,
        level = "transcripts",
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "UcscTranscripts")
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "geneId" = "character",
            "geneName" = "character",
            "txBiotype" = "factor",
            "txChrom" = "factor",
            "txEnd" = "integer",
            "txId" = "character",
            "txName" = "character",
            "txNumber" = "integer",
            "txStart" = "integer",
            "txStrand" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ABCA3P1_2"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr16",
            "start" = "21938865",
            "end" = "21940466",
            "width" = "1602",
            "strand" = "+",
            "geneId" = "ABCA3P1",
            "geneName" = "ABCA3P1",
            "txBiotype" = "transcript",
            "txChrom" = "chr16",
            "txEnd" = "21940466",
            "txId" = "ABCA3P1_2",
            "txName" = "ABCA3P1_2",
            "txNumber" = "144018",
            "txStart" = "21938865",
            "txStrand" = "+"
        )
    )
})

file <- gffs[["wormbase_gtf"]]

test_that("GTF genes", {
    object <- makeGRangesFromGff(
        file = file,
        level = "genes"
    )
    expect_s4_class(object, "WormbaseGenes")
    expect_length(object, 46928L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["geneId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneName" = "character",
            "geneSource" = "factor",
            "geneVersion" = "integer",
            "source" = "factor",
            "type" = "factor"
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
        expected = "WS287"
    )
    expect_identical(
        object = metadata(object)[["md5"]],
        expected = "5c45625e18ee2b5d4b9830da34cd01bf"
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Caenorhabditis elegans"
    )
    expect_identical(
        object = metadata(object)[["release"]],
        expected = "WS287"
    )
    expect_identical(
        object = metadata(object)[["sha256"]],
        expected = "6b8756bd1dfd5995e0c242d35bf3b3d3b04208a65cbf16bc7bfb9befb9a3684f" # nolint
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGff(
        file = file,
        level = "transcripts"
    )
    expect_s4_class(object, "WormbaseTranscripts")
    expect_length(object, 60136L)
    expect_named(
        object = object,
        expected = as.character(mcols(object)[["txId"]])
    )
    expect_identical(
        object = lapply(mcols(object), simpleClass),
        expected = list(
            "broadClass" = "factor",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneName" = "character",
            "geneSource" = "factor",
            "geneVersion" = "integer",
            "source" = "factor",
            "txBiotype" = "factor",
            "txId" = "character",
            "txName" = "character",
            "txSource" = "factor",
            "type" = "factor"
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

file <- file.path(cacheDir, "ref-transcripts.gtf")

test_that("GTF genes", {
    object <- makeGRangesFromGff(
        file = file,
        level = "genes",
        ignoreVersion = TRUE
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
            "broadClass" = "factor",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSource" = "factor",
            "geneSynonyms" = "CompressedCharacterList",
            "ncbiGeneId" = "CompressedIntegerList",
            "source" = "factor",
            "type" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENSG00000284332"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "30366",
            "end" = "30503",
            "width" = "138",
            "strand" = "+",
            "broadClass" = "small",
            "description" = paste(
                "microRNA 1302-2",
                "[Source:HGNC Symbol;Acc:HGNC:35294]"
            ),
            "geneBiotype" = "miRNA",
            "geneId" = "ENSG00000284332",
            "geneIdVersion" = "ENSG00000284332.1",
            "geneName" = "MIR1302-2",
            "geneSource" = "mirbase",
            "geneSynonyms" = "c(\"HSA-MIR-1302-2\", \"MIRN1302-2\")",
            "ncbiGeneId" = "100302278",
            "source" = "mirbase",
            "type" = "gene"
        )
    )
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = "Homo sapiens"
    )
})

test_that("GTF transcripts", {
    object <- makeGRangesFromGff(
        file = file,
        level = "transcripts",
        ignoreVersion = TRUE
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
            "broadClass" = "factor",
            "ccdsId" = "Rle",
            "description" = "character",
            "geneBiotype" = "factor",
            "geneId" = "character",
            "geneIdVersion" = "character",
            "geneName" = "character",
            "geneSource" = "factor",
            "geneSynonyms" = "CompressedCharacterList",
            "ncbiGeneId" = "CompressedIntegerList",
            "source" = "factor",
            "tag" = "CompressedCharacterList",
            "txBiotype" = "factor",
            "txId" = "character",
            "txIdVersion" = "character",
            "txName" = "character",
            "txSource" = "factor",
            "txSupportLevel" = "factor",
            "type" = "factor"
        )
    )
    expect_identical(
        object = vapply(
            X = as.data.frame(object["ENST00000607096"]), # nolint
            FUN = as.character,
            FUN.VALUE = character(1L)
        ),
        expected = c(
            "seqnames" = "chr1",
            "start" = "30366",
            "end" = "30503",
            "width" = "138",
            "strand" = "+",
            "broadClass" = "small",
            "ccdsId" = NA_character_,
            "description" = paste(
                "microRNA 1302-2",
                "[Source:HGNC Symbol;Acc:HGNC:35294]"
            ),
            "geneBiotype" = "miRNA",
            "geneId" = "ENSG00000284332",
            "geneIdVersion" = "ENSG00000284332.1",
            "geneName" = "MIR1302-2",
            "geneSource" = "mirbase",
            "geneSynonyms" = "c(\"HSA-MIR-1302-2\", \"MIRN1302-2\")",
            "ncbiGeneId" = "100302278",
            "source" = "mirbase",
            "tag" = "basic",
            "txBiotype" = "miRNA",
            "txId" = "ENST00000607096",
            "txIdVersion" = "ENST00000607096.1",
            "txName" = "MIR1302-2-201",
            "txSource" = "mirbase",
            "txSupportLevel" = NA_character_,
            "type" = "transcript"
        )
    )
})

## See related issues:
## - https://github.com/Bioconductor/GenomeInfoDb/issues/97
## - https://github.com/Bioconductor/GenomeInfoDb/issues/98

test_that("getChromInfoFromEnsembl Seqinfo failure", {
    file <- pasteUrl(
        "ftp.ensembl.org",
        "pub",
        "release-90",
        "gtf",
        "mus_musculus",
        "Mus_musculus.GRCm38.90.gtf.gz",
        protocol = "ftp"
    )
    object <- makeGRangesFromGff(file, level = "genes")
    expect_s4_class(object, "EnsemblGenes")
})
