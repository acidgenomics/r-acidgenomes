## FIXME Cover fly and worm genome downloads here too.

test_that("Homo sapiens GRCh38", {
    testdir <- tempdir2()
    info <- downloadEnsemblGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 110L,
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_identical(
        object = basename(outputDir),
        expected = "homo-sapiens-grch38-ensembl-110"
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
    expect_s4_class(genes, "EnsemblGenes")
    expect_length(genes, 62754L) # FIXME
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSG00000000003.16",
            "ENSG00000000005.6",
            "ENSG00000000419.14"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "EnsemblTranscripts")
    expect_length(transcripts, 252894L) # FIXME
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
    expect_identical(nrow(t2g), 275741L) # FIXME
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
    expect_identical(nrow(t2g), 275741L) # FIXME
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
    info <- downloadEnsemblGenome(
        organism = "Homo sapiens",
        genomeBuild = "GRCh37",
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_identical(
        object = basename(outputDir),
        expected = "homo-sapiens-grch37-ensembl-87"
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
    expect_s4_class(genes, "EnsemblGenes")
    expect_length(genes, 57905L) # FIXME
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSG00000000003.10",
            "ENSG00000000005.5",
            "ENSG00000000419.8"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "EnsemblTranscripts")
    expect_length(transcripts, 196501L) # FIXME
    expect_identical(
        object = head(sort(names(transcripts)), n = 3L),
        expected = c(
            "ENST00000000233.5",
            "ENST00000000412.3",
            "ENST00000000442.6"
        )
    )
    t2g <- import(file.path(outputDir, "tx2gene.rds"))
    expect_s4_class(t2g, "TxToGene")
    expect_identical(nrow(t2g), 215170L) # FIXME
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.5",
            "geneId" = "ENSG00000004059.6"
        )
    )
    t2g <- import(
        con = file.path(outputDir, "tx2gene.csv.gz"),
        colnames = c("txId", "geneId")
    )
    expect_identical(nrow(t2g), 215170L) # FIXME
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.5",
            "geneId" = "ENSG00000004059.6"
        )
    )
    unlink2(testdir)
})

test_that("Mus musculus GRCm39", {
    testdir <- tempdir2()
    info <- downloadEnsemblGenome(
        organism = "Mus musculus",
        genomeBuild = "GRCm39",
        release = 110L,
        outputDir = testdir,
        cache = TRUE
    )
    outputDir <- info[["args"]][["outputDir"]]
    expect_true(dir.exists(outputDir))
    expect_identical(
        object = basename(outputDir),
        expected = "mus-musculus-grcm39-ensembl-110"
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
    expect_s4_class(genes, "EnsemblGenes")
    expect_length(genes, 56941L) # FIXME
    expect_identical(
        object = head(sort(names(genes)), n = 3L),
        expected = c(
            "ENSMUSG00000000001.5",
            "ENSMUSG00000000003.16",
            "ENSMUSG00000000028.16"
        )
    )
    transcripts <- import(file.path(outputDir, "transcripts.rds"))
    expect_s4_class(transcripts, "EnsemblTranscripts")
    expect_length(transcripts, 149547L) # FIXME
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
    expect_identical(nrow(t2g), 146263L) # FIXME
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
    expect_identical(nrow(t2g), 146263L) # FIXME
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENSMUST00000000001.5",
            "geneId" = "ENSMUSG00000000001.5"
        )
    )
    unlink2(testdir)
})
