#' Convert UCSC build to Ensembl
#'
#' @export
#' @note Updated 2019-08-21.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `character`. Ensembl genome build as the value, UCSC build as the
#'   name. Stops on match failure.
#'
#' @seealso
#' - [UCSC hgGateway](https://genome.ucsc.edu/cgi-bin/hgGateway)
#' - [bcbio genome recipes](https://goo.gl/164J2P)
#'
#' @examples
#' from <- c("hg19", "hg38")
#' to <- convertUCSCBuildToEnsembl(from)
#' print(to)
convertUCSCBuildToEnsembl <- function(object) {
    assert(isCharacter(object))
    ucsc <- .ucscMatrix[, "ucsc"]
    ensembl <- .ucscMatrix[, "ensembl"]
    match <- match(x = object, table = ucsc)
    ## Stop on any match failure.
    if (any(is.na(match))) {
        stop(sprintf(
            "Failed to match UCSC to Ensembl: %s.",
            toString(object[which(is.na(match))], width = 100L)
        ))
    }
    out <- ensembl[match]
    names(out) <- ucsc[match]
    out
}



## nolint start
.ucscMatrix <- matrix(
    dimnames = list(NULL, c(
        "ucsc",       "ensembl",           "date",    "organism",                 "ensembldb"
    )),
    data = c(
        "araTha1",    "TAIR10",            "2011-02", "Arabidopsis thaliana",     FALSE,
        "canFam3",    "CanFam3.1",         "2011-09", "Canis familiaris",         TRUE,
        "ce10",       "WS220",             "2010-10", "Caenorhabditis elegans",   FALSE,
        "ce11",       "WBcel235",          "2013-02", "Caenorhabditis elegans",   TRUE,
        "danRer10",   "GRCz10",            "2014-09", "Danio rerio",              TRUE,
        "danRer11",   "GRCz11",            "2017-05", "Danio rerio",              TRUE,
        "dm3",        "BDGP5",             "2006-04", "Drosophila melanogaster",  FALSE,
        "dm6",        "BDGP6",             "2014-08", "Drosophila melanogaster",  TRUE,
        "galGal4",    "Galgal4",           "2011-11", "Gallus gallus",            FALSE,
        "galGal5",    "Gallus_gallus-5.0", "2015-12", "Gallus gallus",            TRUE,
        "galGal6",    "GRCg6a",            "2018-03", "Gallus gallus",            FALSE,
        "hg19",       "GRCh37",            "2009-02", "Homo sapiens",             TRUE,
        "hg38",       "GRCh38",            "2013-12", "Homo sapiens",             TRUE,
        "hg38-noalt", "GRCh38",            "2013-12", "Homo sapiens",             TRUE,
        "mm9",        "GRCm37",            "2007-07", "Mus musculus",             FALSE,
        "mm10",       "GRCm38",            "2011-12", "Mus musculus",             TRUE,
        "rn5",        "Rnor_5.0",          "2012-03", "Rattus norvegicus",        FALSE,
        "rn6",        "Rnor_6.0",          "2014-07", "Rattus norvegicus",        TRUE,
        "sacCer3",    "R64-1-1",           "2011-04", "Saccharomyces cerevisiae", TRUE,
        "susScr11",   "Sscrofa11.1",       "2017-02", "Sus scrofa",               TRUE,
        "xenTro2",    "JGI_4.1",           "2005-08", "Xenopus tropicalis",       FALSE,
        "xenTro3",    "JGI_4.2",           "2009-11", "Xenopus tropicalis",       TRUE
    ),
    ncol = 5L, byrow = TRUE
)
## nolint end
