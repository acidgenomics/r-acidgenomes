## WormBase
##
## WormBase identifier fix may be needed. WormBase GTF currently imports
## somewhat malformed, and the gene identifiers require additional sanitization
## to return correctly. Look for rows containing "Gene:" and "Transcript:" in ID
## columns.
##
## GTF:
## >  [1] "source"             "type"               "score"
## >  [4] "phase"              "gene_id"            "gene_source"
## >  [7] "gene_biotype"       "transcript_id"      "transcript_source"
## > [10] "transcript_biotype" "exon_number"        "exon_id"
## > [13] "protein_id"



.makeGenesFromWormBaseGTF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(x = "gene_id", y = colnames(mcols(object))),
        areDisjointSets(x = "gene_name", y = colnames(mcols(object)))
    )
    ## Sanitize the `gene_id` column, which can contain some malformed entries
    ## prefixed with "Gene:". We want to keep these entries.
    mcols(object)[["gene_id"]] <- gsub(
        pattern = "^Gene:",
        replacement = "",
        x = mcols(object)[["gene_id"]]
    )
    ## Now safe to only keep rows that match "WBGene" in `gene_id`.
    ## Note that WormBase GTF currently contains some malformed "Transcript:"
    ## entries that we want to drop with this step.
    keep <- grepl(
        pattern = "^WBGene[[:digit:]]{8}$",
        x = mcols(object)[["gene_id"]]
    )
    object <- object[keep, , drop = FALSE]
    ## Process using Ensembl conventions.
    object <- .makeGenesFromEnsemblGTF(object)
    object
}



.makeTranscriptsFromWormBaseGTF <- function(object) {
    assert(
        is(object, "GRanges"),
        isSubset(
            x = c("gene_id", "transcript_id"),
            y = colnames(mcols(object))
        ),
        areDisjointSets(
            x = c("gene_name", "transcript_name"),
            y = colnames(mcols(object))
        )
    )
    ## Sanitize the `gene_id` column, which can contain some malformed entries
    ## prefixed with "Gene:". We want to keep these entries.
    mcols(object)[["gene_id"]] <- gsub(
        pattern = "^Gene:",
        replacement = "",
        x = mcols(object)[["gene_id"]]
    )
    ## Process using Ensembl conventions.
    object <- .makeTranscriptsFromEnsemblGTF(object)
    object
}
