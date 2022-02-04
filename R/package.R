#' AcidGenomes
#'
#' Toolkit for downloading and processing genome annotations.
#'
#' @keywords internal
#'
#' @importClassesFrom GenomicRanges GenomicRanges GenomicRangesList
#' @importClassesFrom IRanges DataFrameList
#' @importClassesFrom S4Vectors DataFrame Vector
#'
#' @importMethodsFrom pipette as.DataFrame coerce
#'
#' @importFrom AcidBase download fileExt forceDetach initDir parentDir pasteURL
#'   printString realpath requireNamespaces showSlotInfo standardizeCall
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning dl h1
#'   toInlineString
#' @importFrom AcidGenerics as.DataFrame encode
#' @importFrom AcidPlyr leftJoin mutateAt rbindToDataFrame
#' @importFrom BiocGenerics %in% append as.data.frame do.call grep grepl
#'   is.unsorted lapply match order rbind setdiff sort unlist
#' @importFrom IRanges CharacterList DataFrameList IntegerList IRanges
#'   SplitDataFrameList gsub ranges trim
#' @importFrom GenomeInfoDb getChromInfoFromEnsembl Seqinfo genome genome<-
#'   seqinfo seqinfo<- seqlevels
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors DataFrame Rle complete.cases decode expand head mcols
#'   mcols<- metadata metadata<- na.omit split tail
#' @importFrom goalie allAreAtomic allAreFiles allAreMatchingFixed
#'   allAreMatchingRegex allAreNotMatchingRegex allAreURLs areDisjointSets
#'   areIntersectingSets areSetEqual assert bapply hasColnames hasDuplicates
#'   hasInternet hasLength hasNoDuplicates hasNames hasRownames hasRows
#'   hasValidNames isADir isAFile isAURL isAny isCharacter isFlag isInstalled
#'   isInt isMatchingFixed isMatchingRegex isNotMatchingFixed isOrganism
#'   isScalar isString isSubset isSystemCommand isWindows validate
#'   validateClasses
#' @importFrom httr GET content content_type stop_for_status
#' @importFrom methods as is isClass new setClass setGeneric setMethod
#'   setValidity signature validObject
#' @importFrom pipette as_tibble cacheURL export getJSON getURL getURLDirList
#'   import md5 removeNA sanitizeNA sha256
#' @importFrom stringr boundary str_extract str_match str_split_fixed
#' @importFrom syntactic camelCase kebabCase makeNames snakeCase upperCamelCase
#' @importFrom utils capture.output packageName packageVersion sessionInfo
"_PACKAGE"



#' @name params
#' @inherit AcidRoxygen::params return title
#' @keywords internal
#'
#' @param cache `logical(1)`.
#'   Cache URLs locally, using BiocFileCache internally.
#' @param ignoreVersion `logical(1)`.
#'   Ignore identifier (e.g. transcript, gene) versions.
#'   When applicable, the identifier containing version numbers will be stored
#'   in `txIdVersion` and `geneIdVersion`, and the variants without versions
#'   will be stored in `txId`, `txIdNoVersion`, `geneId`, and `geneIdNoVersion`.
#' @param synonyms `logical(1)`.
#'   Include gene synonyms.
#'   Queries the Ensembl web server, and is CPU intensive.
NULL
