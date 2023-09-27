test_that("downloadEnsemblGenome", {
    testdir <- tempdir2()
    info <- downloadEnsemblGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 108L,
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
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
    expect_s4_class(genes, "EnsemblGenes")
    expect_length(genes, 62703L)
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSG00000000003.15",
            "ENSG00000000005.6",
            "ENSG00000000419.14"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "EnsemblTranscripts")
    expect_length(transcripts, 252301L)
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
    expect_identical(nrow(t2g), 272929L)
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
    expect_identical(nrow(t2g), 272929L)
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    unlink2(testdir)
})

test_that("downloadGencodeGenome", {
    testdir <- tempdir2()
    info <- downloadGencodeGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 42L,
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
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
    expect_length(genes, 62696L)
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSG00000000003.15",
            "ENSG00000000005.6",
            "ENSG00000000419.14"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "GencodeTranscripts")
    expect_length(transcripts, 252416L)
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
    expect_identical(nrow(t2g), 252416L)
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
    expect_identical(nrow(t2g), 252416L)
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    unlink2(testdir)
})

test_that("downloadRefseqGenome", {
    testdir <- tempdir2()
    info <- downloadRefseqGenome(
        organism = "Homo sapiens",
        taxonomicGroup = "vertebrate_mammalian",
        genomeBuild = "GCF_000001405.40_GRCh38.p14",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
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
    expect_s4_class(genes, "RefseqGenes")
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c("A1BG", "A1BG-AS1", "A1CF")
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "RefseqTranscripts")
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c("NM_000014.6", "NM_000015.3", "NM_000016.6")
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    t2g <- as.data.frame(t2g)
    aatfExpected <- data.frame(
        "txId" = c(
            "NM_001411094.1",
            "NM_012138.4",
            "XM_047435748.1",
            "XM_054329280.1"
        ),
        "geneId" = rep("AATF", 4L)
    )
    aatfCurrent <- t2g[t2g[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    aatfCurrent <- t2g[t2g[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    unlink2(testdir)
})

test_that("downloadUcscGenome", {
    testdir <- tempdir2()
    info <- downloadUcscGenome(
        organism = "Homo sapiens",
        genomeBuild = "hg38",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_true(all(file.exists(file.path(
        outputDir,
        c(
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
    expect_s4_class(genes, "UcscGenes")
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSG00000000003",
            "ENSG00000000005",
            "ENSG00000000419"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "UcscTranscripts")
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c(
            "ENST00000000233",
            "ENST00000000412",
            "ENST00000000442"
        )
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233",
            "geneId" = "ENSG00000004059"
        )
    )
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233",
            "geneId" = "ENSG00000004059"
        )
    )
    unlink2(testdir)
})
