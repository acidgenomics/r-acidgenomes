## Merge the gene-level annotations (`geneName`, `geneBiotype`) into a
## transcript-level GRanges object.
## Updated 2020-01-20.
.mergeGenesIntoTranscripts <- function(transcripts, genes) {
    ## Use `transcript_` prefix instead of `tx_` consistently.
    if (any(grepl("^tx_", colnames(mcols(transcripts))))) {
        colnames(mcols(transcripts)) <- gsub(
            pattern = "^tx_",
            replacement = "transcript_",
            x = colnames(mcols(transcripts))
        )
    }
    assert(
        is(transcripts, "GRanges"),
        is(genes, "GRanges"),
        ## Note that `hasValidNames()` will error on WormBase transcripts.
        hasNames(transcripts),
        hasNames(genes),
        isSubset("transcript_id", colnames(mcols(transcripts))),
        identical(names(transcripts), mcols(transcripts)[["transcript_id"]]),
        ## Don't proceed unless we have `gene_id` column to use for merge.
        isSubset("gene_id", colnames(mcols(transcripts))),
        isSubset("gene_id", colnames(mcols(genes)))
    )
    geneCols <- setdiff(
        x = colnames(mcols(genes)),
        y = colnames(mcols(transcripts))
    )
    ## Only attempt the merge if there's useful additional metadata to include.
    ## Note that base `merge()` can reorder rows, so be careful here.
    if (length(geneCols) > 0L) {
        cli_alert(sprintf(
            "Merging gene-level annotations: {.var %s}.",
            toString(camelCase(geneCols), width = 100L)
        ))
        geneCols <- c("gene_id", geneCols)
        ## x: transcripts; y: genes
        x <- mcols(transcripts)
        y <- mcols(genes)[, geneCols, drop = FALSE]
        assert(
            isSubset(
                x = unique(x[["gene_id"]]),
                y = unique(y[["gene_id"]])
            ),
            hasNoDuplicates(y[["gene_id"]])
        )
        merge <- leftJoin(x = x, y = y, by = "gene_id")
        assert(identical(x[["transcript_id"]], merge[["transcript_id"]]))
        mcols(transcripts) <- merge
    }
    transcripts
}
