## FIXME NEED TO INCLUDE IGNOREVERSION HERE...



#' Make a Gene2Symbol object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeGene2Symbol
#' @note Updated 2021-01-14.
#'
#' @inheritParams Gene2Symbol
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @return `Gene2Symbol`.
#'
#' @examples
#' ## makeGene2SymbolFromEnsembl ====
#' x <- makeGene2SymbolFromEnsembl(organism = "Homo sapiens")
#' print(x)
#'
#' ## makeTx2GeneFromEnsDb ====
#' if ("EnsDb.Hsapiens.v75" %in% rownames(installed.packages())) {
#'     x <- makeGene2SymbolFromEnsDb("EnsDb.Hsapiens.v75")
#'     print(x)
#' }
#'
#' ## makeGene2SymbolFromGFF ====
#' ## GTF
#' file <- file.path(AcidGenomesTestsURL, "example.gtf")
#' x <- makeGene2SymbolFromGFF(file)
#' print(x)
#'
#' ## GFF3
#' file <- file.path(AcidGenomesTestsURL, "example.gff3")
#' x <- makeGene2SymbolFromGFF(file)
#' print(x)
NULL



#' @describeIn makeGene2Symbol Make a `Gene2Symbol` object from Ensembl using
#'   an AnnotationHub lookup.
#' @export
## Updated 2021-01-14.
makeGene2SymbolFromEnsembl <-
    function(
        organism,
        ...
    ) {
        gr <- do.call(
            what = .makeGRangesFromEnsembl,
            args = matchArgsToDoCall(
                args = list(
                    level = "genes",
                    broadClass = FALSE,
                    synonyms = FALSE
                ),
                removeFormals = "format"
            )
        )
        Gene2Symbol(object = gr, format = match.arg(format))
    }

f <- formals(makeGRangesFromEnsembl)
f <- f[setdiff(names(f), c("level", "synonyms"))]
f[["format"]] <- formals(`Gene2Symbol,DataFrame`)[["format"]]
formals(makeGene2SymbolFromEnsembl) <- f



#' @describeIn makeGene2Symbol Make a `Gene2Symbol` object from an `EnsDb`
#'   object or annotation package.
#' @export
## Updated 2021-01-14.
makeGene2SymbolFromEnsDb <-
    function(object, format) {
        gr <- .makeGRangesFromEnsDb(
            object = object,
            broadClass = FALSE,
            synonyms = FALSE
        )
        Gene2Symbol(object = gr, format = match.arg(format))
    }

formals(makeGene2SymbolFromEnsDb)[["format"]] <-
    formals(makeGene2SymbolFromEnsembl)[["format"]]



## FIXME NEED TO CREATE A UNIT TEST TO ENSURE VERSION REMOVAL WORKS HERE.

#' @describeIn makeGene2Symbol Make a `Gene2Symbol` object from a GFF file.
#' @export
## Updated 2020-01-14.
makeGene2SymbolFromGFF <-
    function(
        file,
        ignoreVersion = TRUE,
        format
    ) {
        gr <- .makeGRangesFromGFF(
            file = file,
            level = "genes",
            ignoreVersion = ignoreVersion,
            broadClass = FALSE,
            synonyms = FALSE
        )
        Gene2Symbol(object = gr, format = match.arg(format))
    }

formals(makeGene2SymbolFromGFF)[["format"]] <-
    formals(makeGene2SymbolFromEnsembl)[["format"]]



#' @describeIn makeGene2Symbol GTF alias for `makeGene2SymbolFromGFF`.
#' @export
makeGene2SymbolFromGTF <- makeGene2SymbolFromGFF
