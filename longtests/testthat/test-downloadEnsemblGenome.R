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
    expect_length(genes, 62754L)
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
    expect_length(transcripts, 252894L)
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
    expect_identical(nrow(t2g), 275741L)
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
    expect_identical(nrow(t2g), 275741L)
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
        )
    )
    unlink2(testdir)
})

# FIXME We're hitting an Ensembl FTP server file error related to Seqinfo:
#→ Making <GRanges> from GFF file (/private/var/folders/l1/8y8sjzmn15v49jgrqglghcfr0000gn/T/RtmpZXJGE7/RrELyzpqIc-170170865086322/homo-sapiens-grch37-ensembl-87/annotation/gtf/Homo_sapiens.GRCh37.87.gtf.gz).
#→ Importing /private/var/folders/l1/8y8sjzmn15v49jgrqglghcfr0000gn/T/RtmpZXJGE7/RrELyzpqIc-170170865086322/homo-sapiens-grch37-ensembl-87/annotation/gtf/Homo_sapiens.GRCh37.87.gtf.gz using rtracklayer::`import()`.
#Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#                  line 498 did not have 10 elements
#              Calls: makeGRangesFromGff ... simple_read_table -> do.call -> do.call -> <Anonymous> -> scan
#              Backtrace:
#                  ▆
#              1. └─AcidGenomes::makeGRangesFromGff(file)
#              2.   └─AcidGenomes:::.makeGRangesFromRtracklayer(...) at AcidGenomes/R/makeGRangesFromGff.R:269:13
#              3.     └─AcidGenomes:::.getSeqinfo(meta) at AcidGenomes/R/internal-rtracklayer.R:62:9
#              4.       └─AcidGenomes:::.getEnsemblSeqinfo(...) at AcidGenomes/R/internal-genomeinfodb.R:93:13
#              5.         ├─goalie::quietly(...) at AcidGenomes/R/internal-genomeinfodb.R:152:5
#              6.         │ └─base::withCallingHandlers(...)
#              7.         ├─BiocGenerics::do.call(what = getChromInfoFromEnsembl, args = args) at AcidGenomes/R/internal-genomeinfodb.R:153:9
#              8.         ├─base::do.call(what = getChromInfoFromEnsembl, args = args)
#              9.         └─GenomeInfoDb (local) `<fn>`(...)
#              10.           └─GenomeInfoDb:::.get_chrom_info_from_Ensembl_FTP(...)
#              11.             └─GenomeInfoDb:::fetch_seq_regions_from_Ensembl_FTP(...)
#              12.               └─GenomeInfoDb:::.fetch_synonyms_from_Ensembl_FTP(...)
#              13.                 └─GenomeInfoDb:::fetch_table_from_Ensembl_FTP(core_db_url, "external_db")
#              14.                   └─GenomeInfoDb:::fetch_table_from_url(...)
#              15.                     └─GenomeInfoDb:::simple_read_table(destfile, ...)
#              16.                       ├─BiocGenerics::do.call(read.table, args)
#              17.                       ├─base::do.call(read.table, args)
#              18.                       └─utils (local) `<fn>`(...)
#              19.                         └─base::scan(...)

## FIXME Check that this maps to release 87.
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
    expect_length(genes, 62754L)
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
    expect_length(transcripts, 252894L)
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
    expect_identical(nrow(t2g), 275741L)
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
    expect_identical(nrow(t2g), 275741L)
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENST00000000233.10",
            "geneId" = "ENSG00000004059.11"
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
    expect_length(genes, 56941L)
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
    expect_length(transcripts, 149547L)
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
    expect_identical(nrow(t2g), 146263L)
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
    expect_identical(nrow(t2g), 146263L)
    expect_identical(
        object = as.data.frame(t2g)[1L, , drop = TRUE],
        expected = list(
            "txId" = "ENSMUST00000000001.5",
            "geneId" = "ENSMUSG00000000001.5"
        )
    )
    unlink2(testdir)
})
