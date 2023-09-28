lst <- AcidDevTools::cacheTestFiles(
    pkg = .pkgName,
    files = c(
        "cellranger-genes.gtf.gz",
        "ensembl.gff3",
        "ensembl.gtf",
        "flybase.gtf",
        "gencode.gff3",
        "gencode.gtf",
        "refseq.gff3",
        "refseq.gtf",
        "ref-transcripts.gtf",
        "tx2gene.csv",
        "ucsc.gtf",
        "wormbase.gtf"
    )
)
cacheDir <- lst[["cacheDir"]]
rm(lst)
