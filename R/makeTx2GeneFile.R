## Updated 2021-01-07.
#' @rdname makeTx2GeneFromFASTA
#' @param outputFile `character(1)`.
#'   Output file.
#' @export
makeTx2GeneFileFromFASTA <- function(
    file,
    outputFile = file.path(dirname(file), "tx2gene.csv.gz"),
    source
) {
    t2g <- makeTx2GeneFromFASTA(file = file, source = source)
    assert(is(t2g, "Tx2Gene"))
    export(object = t2g, file = outputFile, overwrite = TRUE)
}



## FIXME CAN WE GET REFSEQ TO A FLAT GRANGES?
## FIXME NOT SUPPORTED FOR REFSEQ.
## Error in (function (classes, fdef, mtable)  :
##               unable to find an inherited method for function 'Tx2Gene' for signature '"CompressedGRangesList"'
##           Calls: makeTx2GeneFromGFF -> Tx2Gene -> <Anonymous>

## Updated 2021-01-08.
#' @rdname makeTx2Gene
#' @export
makeTx2GeneFileFromGFF <- function(
    file,
    outputFile = file.path(dirname(file), "tx2gene.csv.gz")
) {
    t2g <- makeTx2GeneFromGFF(file)
    assert(is(t2g, "Tx2Gene"))
    export(object = t2g, file = outputFile, overwrite = TRUE)
}

#' @rdname makeTx2Gene
#' @export
makeTx2GeneFileFromGTF <- makeTx2GeneFileFromGFF
