#' Make a Gene2Symbol object
#'
#' @section GFF/GTF file:
#'
#' Remote URLs and compressed files are supported.
#'
#' @name makeGene2Symbol
#' @note Updated 2020-10-05.
#'
#' @inheritParams Gene2Symbol
#' @inheritParams AcidRoxygen::params
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
#' file <- file.path(basejumpTestsURL, "example.gtf")
#' x <- makeGene2SymbolFromGFF(file)
#' print(x)
#'
#' ## GFF3
#' file <- file.path(basejumpTestsURL, "example.gff3")
#' x <- makeGene2SymbolFromGFF(file)
#' print(x)
NULL



#' @rdname makeGene2Symbol
#' @export
## Updated 2020-10-05.
makeGene2SymbolFromEnsembl <-
    function(organism) {
        format <- match.arg(format)
        gr <- do.call(
            what = makeGRangesFromEnsembl,
            args = matchArgsToDoCall(
                args = list(
                    level = "genes",
                    synonyms = FALSE
                ),
                removeFormals = "format"
            )
        )
        Gene2Symbol(gr, format = format)
    }

f <- formals(makeGRangesFromEnsembl)
f <- f[setdiff(names(f), "level")]
f[["format"]] <- formals(`Gene2Symbol,DataFrame`)[["format"]]
formals(makeGene2SymbolFromEnsembl) <- f



#' @rdname makeGene2Symbol
#' @export
## Updated 2020-10-05.
makeGene2SymbolFromEnsDb <-
    function(object, format) {
        format <- match.arg(format)
        gr <- makeGRangesFromEnsDb(
            object = object,
            synonyms = FALSE
        )
        Gene2Symbol(gr, format = format)
    }

formals(makeGene2SymbolFromEnsDb)[["format"]] <-
    formals(makeGene2SymbolFromEnsembl)[["format"]]



#' @rdname makeGene2Symbol
#' @export
## Updated 2020-10-05.
makeGene2SymbolFromGFF <-
    function(file, format) {
        format <- match.arg(format)
        gr <- makeGRangesFromGFF(
            file = file,
            level = "genes",
            synonyms = FALSE
        )
        Gene2Symbol(gr, format = format)
    }

formals(makeGene2SymbolFromGFF)[["format"]] <-
    formals(makeGene2SymbolFromEnsembl)[["format"]]



#' @rdname makeGene2Symbol
#' @usage NULL
#' @export
makeGene2SymbolFromGTF <- makeGene2SymbolFromGFF
