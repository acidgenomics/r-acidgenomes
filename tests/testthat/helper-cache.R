if (!isTRUE(goalie::hasInternet())) {
    warning("No Internet connection detected.")
    return(invisible(NULL))
}
dir.create("cache", showWarnings = FALSE)
files <- c(
    "cellranger-genes.gtf.gz",
    "ensembl.gff3",
    "ensembl.gtf",
    "flybase.gtf",
    "gencode.gff3",
    "gencode.gtf",
    "refseq.gff3",
    "refseq.gtf",
    "tx2gene.csv",
    "ucsc.gtf",
    "wormbase.gtf"
)
mapply(
    FUN = function(remoteDir, file, envir) {
        destfile <- file.path("cache", file)
        if (!file.exists(destfile)) {
            utils::download.file(
                url = paste(remoteDir, file, sep = "/"),
                destfile = destfile
            )
        }
    },
    file = files,
    MoreArgs = list(
        "remoteDir" = AcidGenomesTestsURL,
        "envir" = environment()
    )
)
rm(files)
