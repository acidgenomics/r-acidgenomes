#' @name params
#' @inherit AcidRoxygen::params return title
#' @keywords internal
#'
#' @param cache `logical(1)`.
#' Cache URLs locally, using BiocFileCache internally.
#'
#' @param extraMcols `logical(1)`.
#' Include extra metadata columns (e.g. `"broadClass"`).
#'
#' @param ignoreCase `logical(1)`.
#' Perform case-insensitive matching.
#'
#' @param ignoreVersion `logical(1)`.
#' Ignore identifier (e.g. transcript, gene) versions.
#' When applicable, the identifier containing version numbers will be stored
#' in `txIdVersion` and `geneIdVersion`, and the variants without versions
#' will be stored in `txId`, `txIdNoVersion`, `geneId`, and `geneIdNoVersion`.
NULL
