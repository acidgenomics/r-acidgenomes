## nolint start
suppressPackageStartupMessages({
    library(devtools)
    library(AcidBase)
    library(S4Vectors)
    library(goalie)
    library(pipette)
    library(utils)
})
## nolint end
load_all(helpers = FALSE)
url <- pasteURL(
    "ftp.ncbi.nlm.nih.gov",
    "pub",
    "taxonomy",
    "taxdump.tar.gz",
    protocol = "ftp"
)
tarfile <- .cacheIt(url)
exdir <- tempdir2()
untar(tarfile = tarfile, exdir = exdir)
con <- file.path(exdir, "images.dmp")
df <- import(con = con, format = "tsv", colnames = FALSE)
unlink2(exdir)
df <- as(df, "DFrame")
df <- df[, c(15L, 3L)]
colnames(df) <- c("taxonomyId", "organism")
assert(
    hasNoDuplicates("taxonomyId"),
    hasNoDuplicates("organism")
)
df <- df[complete.cases(df), , drop = FALSE]
df[["taxonomyId"]] <- as.integer(df[["taxonomyId"]])
df[["organism"]] <- sub(
    pattern = "image:",
    replacement = "",
    x = df[["organism"]],
    fixed = TRUE
)
df <- df[order(df[["taxonomyId"]]), , drop = FALSE]
rownames(df) <- df[["taxonomyId"]]
saveRDS(df, "mapNcbiTaxId.rds")
export(df, "mapNcbiTaxId.csv")
