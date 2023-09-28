if (!isTRUE(goalie::hasInternet())) {
    warning("No Internet connection detected.")
    return(invisible(NULL))
}
cacheDir <- file.path(
    tools::R_user_dir(package = .pkgName, which = "cache"),
    "testthat"
)
dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)
files <- c(
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
Map(
    f = function(remoteDir, file, envir) {
        destfile <- file.path(cacheDir, file)
        if (!file.exists(destfile)) {
            utils::download.file(
                url = paste(remoteDir, file, sep = "/"),
                destfile = destfile
            )
        }
    },
    file = files,
    MoreArgs = list(
        "remoteDir" = AcidGenomesTestsUrl,
        "envir" = environment()
    )
)
rm(files)
