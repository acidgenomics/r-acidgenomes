#' Detect organism
#'
#' Supports organism detection from Ensembl identifier or genome build.
#'
#' Only the first match is returned. We're using a while loop approach here so
#' we can skip transgenes or spike-ins. The function fails after a maximum of 50
#' unknowns, for speed.
#'
#' @section Supported organisms:
#'
#' - *Caenorhabditis elegans* (roundworm)
#' - *Danio rerio* (zebrafish)
#' - *Drosophila melanogaster* (fruitfly)
#' - *Gallus gallus* (chicken)
#' - *Homo sapiens* (human)
#' - *Mus musculus* (mouse)
#' - *Ovis aries* (sheep)
#' - *Rattus norvegicus* (rat)
#' - *Saccharomyces cerevisiae* (yeast)
#'
#' @export
#' @note [BiocGenerics::organism()] `character` method conflicts with annotate
#'   package, which gets loaded into the namespace when DESeq2 is attached.
#'   Instead, we're exporting the character method here as a separate function
#'   named `detectOrganism`.
#' @note Updated 2020-10-06.
#'
#' @param object `character`.
#'
#' @return `character(1)`.
#' Full latin organism name. Stops on match failure.
#'
#' @seealso `BiocGenerics::organism()`.
#'
#' @examples
#' ## Match by gene identifier.
#' detectOrganism("ENSG00000000003")
#'
#' ## Match by genome build.
#' detectOrganism("GRCh38")  # Ensembl
#' detectOrganism("hg38")    # UCSC
#'
#' ## Match by alternate organism name.
#' detectOrganism("H. sapiens")
#' detectOrganism("hsapiens")
#'
#' ## The function will skip transgenes/spike-ins until we find a match.
#' detectOrganism(c("EGFP", "TDTOMATO", "ENSG00000000003"))
#'
#' ## But it only returns the first match, if there are multiple genomes.
#' detectOrganism(c("ENSG00000000003", "ENSMUSG00000000001"))
detectOrganism <- function(object) {
    assert(isCharacter(object))
    data <- import(
        file = system.file(
            "extdata", "detectOrganism.rds",
            package = packageName()
        ),
        quiet = TRUE
    )
    ## Parse the character vector until we get a match.
    x <- NA_character_
    i <- 1L
    while (
        is.na(x) &&
        i <= min(length(object), 50L)
    ) {
        x <- .detectOrganism(string = object[[i]], data = data)
        i <- i + 1L
    }
    if (is.na(x)) {
        stop("Failed to detect organism.")
    } else {
        x
    }
}



## Updated 2020-10-06.
.detectOrganism <-  # nolint
    function(string, data) {
        ## Generate a logical matrix of grep matches.
        hits <- apply(
            X = data,
            MARGIN = 1L,
            FUN = function(row) {
                any(vapply(
                    X = row,
                    FUN = function(pattern) {
                        grepl(pattern = pattern, x = string, ignore.case = TRUE)
                    },
                    FUN.VALUE = integer(1L)
                ))
            }
        )
        ## Return organism name if there's a match, otherwise NA.
        ifelse(
            test = any(hits),
            yes = data[["organism"]][which(hits)],
            no = NA_character_
        )
    }
