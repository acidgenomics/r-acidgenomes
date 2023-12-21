test_that("Homo sapiens GRCh38", {
    testdir <- tempdir2()
    info <- downloadGencodeGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 44L,
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_identical(
        object = basename(outputDir),
        expected = "homo-sapiens-grch38-gencode-44"
    )
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
            "annotation.gff3.gz",
            "annotation.gtf.gz",
            "genes.rds",
            "genome.fa.gz",
            "metadata.rds",
            "transcriptome.fa.gz",
            "transcripts.rds",
            "tx2gene.csv.gz",
            "tx2gene.rds"
        )
    ))))
    genes <- import(file.path(outputDir, "genes.rds"))
    expect_s4_class(genes, "GencodeGenes")
    expect_length(genes, 62700L) # FIXME
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSG00000000003.16",
            "ENSG00000000005.6",
            "ENSG00000000419.14"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "GencodeTranscripts")
    expect_length(transcripts, 252835L) # FIXME
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c(
            "ENST00000000233.10",
            "ENST00000000412.8",
            "ENST00000000442.11"
        )
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    expect_identical(nrow(t2g), 252835L) # FIXME
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(nrow(t2g), 252835L) # FIXME
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    unlink2(testdir)
})

test_that("Homo sapiens GRCh37", {
    testdir <- tempdir2()
    info <- downloadGencodeGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh37",
        release = 44L,
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_identical(
        object = basename(outputDir),
        expected = "homo-sapiens-grch37-gencode-44"
    )
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
            "annotation.gff3.gz",
            "annotation.gtf.gz",
            "genes.rds",
            "genome.fa.gz",
            "metadata.rds",
            "transcriptome.fa.gz",
            "transcripts.rds",
            "tx2gene.csv.gz",
            "tx2gene.rds"
        )
    ))))
    genes <- import(file.path(outputDir, "genes.rds"))
    expect_s4_class(genes, "GencodeGenes")
    expect_length(genes, 64486L) # FIXME
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSG00000000003.16_9",
            "ENSG00000000005.6_6",
            "ENSG00000000419.14_12"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "GencodeTranscripts")
    expect_length(transcripts, 254413L) # FIXME
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c(
            "ENST00000000233.10_3",
            "ENST00000000412.3",
            "ENST00000000442.11_7"
        )
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    expect_identical(nrow(t2g), 254413L) # FIXME
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10_3",
            "geneId" = "ENSG00000004059.11_8"
        )
    )
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(nrow(t2g), 254413L) # FIXME
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10_3",
            "geneId" = "ENSG00000004059.11_8"
        )
    )
    unlink2(testdir)
})

test_that("Mus musculus GRCm39", {
    testdir <- tempdir2()
    info <- downloadGencodeGenome(
        organism = "Mus musculus",
        genomeBuild = "GRCm39",
        release = "M33",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_identical(
        object = basename(outputDir),
        expected = "mus-musculus-grcm39-gencode-m33"
    )
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
            "annotation.gff3.gz",
            "annotation.gtf.gz",
            "genes.rds",
            "genome.fa.gz",
            "metadata.rds",
            "transcriptome.fa.gz",
            "transcripts.rds",
            "tx2gene.csv.gz",
            "tx2gene.rds"
        )
    ))))
    genes <- import(file.path(outputDir, "genes.rds"))
    expect_s4_class(genes, "GencodeGenes")
    expect_length(genes, 56884L) # FIXME
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSMUSG00000000001.5",
            "ENSMUSG00000000003.16",
            "ENSMUSG00000000028.16"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "GencodeTranscripts")
    expect_length(transcripts, 149488L) # FIXME
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c(
            "ENSMUST00000000001.5",
            "ENSMUST00000000003.14",
            "ENSMUST00000000010.9"
        )
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    expect_identical(nrow(t2g), 149488L) # FIXME
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENSMUST00000000001.5",
            "geneId" = "ENSMUSG00000000001.5"
        )
    )
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(nrow(t2g), 149488L) # FIXME
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENSMUST00000000001.5",
            "geneId" = "ENSMUSG00000000001.5"
        )
    )
    unlink2(testdir)
})
