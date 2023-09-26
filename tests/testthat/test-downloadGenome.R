test_that("downloadEnsemblGenome", {
    testdir <- tempdir2()
    info <- downloadEnsemblGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 108L,
        outputDir = testdir,
        cache = FALSE
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
            "txToGene.csv.gz",
            "txToGene.rds"
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
    txToGene <- import(file.path(outputDir, "txToGene.rds"))
    expect_s4_class(txToGene, "TxToGene")
    expect_identical(nrow(txToGene), 272929L)
    expect_identical(
        object = as.data.frame(txToGene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    txToGene <- import(
        con = file.path(outputDir, "txToGene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(nrow(txToGene), 272929L)
    expect_identical(
        object = as.data.frame(txToGene)[1L, , drop = TRUE],
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
        cache = FALSE
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
            "txToGene.csv.gz",
            "txToGene.rds"
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
    txToGene <- import(file.path(outputDir, "txToGene.rds"))
    expect_s4_class(txToGene, "TxToGene")
    expect_identical(nrow(txToGene), 252416L)
    expect_identical(
        object = as.data.frame(txToGene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    txToGene <- import(
        con = file.path(outputDir, "txToGene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(nrow(txToGene), 252416L)
    expect_identical(
        object = as.data.frame(txToGene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    unlink2(testdir)
})

test_that("downloadRefSeqGenome", {
    testdir <- tempdir2()
    info <- downloadRefSeqGenome(
        organism = "Homo sapiens",
        taxonomicGroup = "vertebrate_mammalian",
        genomeBuild = "GCF_000001405.40_GRCh38.p14",
        outputDir = testdir,
        cache = FALSE
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
            "txToGene.csv.gz",
            "txToGene.rds"
        )
    ))))
    genes <- import(file.path(outputDir, "genes.rds"))
    expect_s4_class(genes, "RefSeqGenes")
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c("A1BG", "A1BG-AS1", "A1CF")
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "RefSeqTranscripts")
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c("NM_000014.6", "NM_000015.3", "NM_000016.6")
    )
    txToGene <- import(file.path(outputDir, "txToGene.rds"))
    expect_s4_class(txToGene, "TxToGene")
    txToGene <- as.data.frame(txToGene)
    aatfExpected <- data.frame(
        "txId" = c(
            "NM_001411094.1",
            "NM_012138.4",
            "XM_047435748.1",
            "XM_054329280.1"
        ),
        "geneId" = rep("AATF", 4L)
    )
    aatfCurrent <- txToGene[txToGene[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    txToGene <- import(
        con = file.path(outputDir, "txToGene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    aatfCurrent <- txToGene[txToGene[, 2L] == "AATF", ]
    rownames(aatfCurrent) <- NULL
    expect_identical(aatfCurrent, aatfExpected)
    unlink2(testdir)
})

test_that("downloadUCSCGenome", {
    testdir <- tempdir2()
    info <- downloadUCSCGenome(
        organism = "Homo sapiens",
        genomeBuild = "hg38",
        outputDir = testdir,
        cache = FALSE
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
            "txToGene.csv.gz",
            "txToGene.rds"
        )
    ))))
    genes <- import(file.path(outputDir, "genes.rds"))
    expect_s4_class(genes, "UCSCGenes")
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSG00000000003",
            "ENSG00000000005",
            "ENSG00000000419"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "UCSCTranscripts")
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c(
            "ENST00000000233",
            "ENST00000000412",
            "ENST00000000442"
        )
    )
    txToGene <- import(file.path(outputDir, "txToGene.rds"))
    expect_s4_class(txToGene, "TxToGene")
    expect_identical(
        object = as.data.frame(txToGene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233",
            "geneId" = "ENSG00000004059"
        )
    )
    txToGene <- import(
        con = file.path(outputDir, "txToGene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(
        object = as.data.frame(txToGene)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233",
            "geneId" = "ENSG00000004059"
        )
    )
    unlink2(testdir)
})
